#!/usr/bin/env python3
import pandas as pd
import re
from pathlib import Path
import sys
import tkinter as tk
from pandastable import Table

# Zentrale Konfiguration importieren
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, PLOTS_DIR


STEM = "Asgard__2026747_173"
OUT_DIR = TRNASCAN_OUTPUT_DIR

# KONFIGURATION: Verarbeitungsmodus
ALL = True  # True = alle Organismen verarbeiten, False = nur einen spezifischen Organismus
SINGLE_STEM = "Asgard__2026747_173"  # Nur relevant wenn ALL = False

#Debugging Mode
DEBUG = False  # FALSE um GUI- POPups zu vermeiden

if DEBUG:
    stem = "Asgard__2026747_173"
    for suffix in ["main", "detail", "struct", "mainstruct"]:
        f = OUT_DIR / f"{stem}_{suffix}.txt"
        print(f"[DEBUG] {f.name} existiert: {f.exists()}")

print(f"[INFO] Folgende Suche nach tRNAscan-SE Output in: {OUT_DIR}")
print(f"[INFO] Zukünftige Plots werden gespeichert in: {PLOTS_DIR}")

#Hilfsfunktion für Tabellenausgabe (pandas)
def show_df(df: pd.DataFrame, title: str = "Vorschau"):
    #Zeigt ein pandas-DataFrame in einer GUI mit pandastable.
    #Konvertiert mittlerweile auch Categorical-Spalten in Strings, um GUI-Probleme zu vermeiden -> gamechanger
    
    # Wenn Einzelorganismus-Modus, DEBUG-Verhalten beibehalten
    # Wenn Modus alle verarbeiten, nur Konsolen-Output (außer DEBUG explizit True)
    should_show_gui = DEBUG and (not ALL or len(df) < 100)  # GUI nur bei kleinen DataFrames im Batch-Modus
    
    if not should_show_gui:
        print(f"\n--- {title} ---")
        print(df.head())
        print(f"Shape: {df.shape}")
        return
        
    df_safe = df.copy()
    for col in df_safe.select_dtypes(include="category").columns:
        df_safe[col] = df_safe[col].astype("string")

    root = tk.Tk()
    root.title(title)
    frame = tk.Frame(root)
    frame.pack(fill="both", expand=True)
    pt = Table(frame, dataframe=df_safe, showtoolbar=True, showstatusbar=True)
    pt.show()
    root.mainloop()

def get_unique_stems(directory: Path) -> list[str]:
    #Findet alle eindeutigen Stems (Organismus-Identifier) im Verzeichnis.
    #Sucht nach *_main.txt Dateien und extrahiert den Stem-Teil.
    
    main_files = list(directory.glob("*_main.txt"))
    stems = []
    
    for file in main_files:
        # Entfernt "_main.txt" vom Dateinamen
        stem = file.stem.replace("_main", "")
        stems.append(stem)
    
    return sorted(stems)


###############################################################################
#MAIN Parser-Funktionen
###############################################################################
# Intron-Feld Parser
def _parse_intron_field(field: str | int | float) -> list[int]:
    if pd.isna(field) or str(field).strip() in {"", "0", "-"}:
        return []
    return [int(x) for x in str(field).split(",") if x.strip().isdigit()]
def _row_to_pairs(row) -> list[tuple[int, int]]:
    begins = _parse_intron_field(row["IntronBegin"])
    ends = _parse_intron_field(row["IntronEnd"])
    return list(zip(begins, ends))

# Main-Datei einlesen und parsen
def parse_main_file(stem: str, directory: Path) -> pd.DataFrame:
    file = directory / f"{stem}_main.txt"
    if not file.exists():
        raise FileNotFoundError(f"{file} nicht gefunden")

    if DEBUG:
        print(f"[DEBUG] Lese Datei: {file}")

        # Manuell gesetzter Header – robust gegen Dubletten
    cols_main = [
        "Sequence", "tRNA#", "Begin", "End", "Type", "Anticodon",
        "IntronBegin", "IntronEnd",
        "CMScore", "HMMScore", "Str2Score",
        "Isotype", "IsotypeScore", "Note"
    ]

    df = pd.read_csv(
        file,
        sep=r'\s+',
        skiprows=3,  # 1x Kommentar, 1x Header, 1x Bindestrich-Zeile
        names=cols_main,
        engine="python"
    )



    # Entferne die Dummy-Zeile mit Bindestrichen falls sie enthalten ist
    df = df[~df.iloc[:, 0].str.startswith("-")]

    # Umbenennung für Einheitlichkeit
    df = df.rename(columns={
        "Sequence": "Sequence",
        "tRNA": "tRNA#",
        "Begin": "Begin",
        "End": "End",
        "Type": "Type",
        "Anti": "Anticodon",
        "Intron": "IntronBegin",   # falls automatisch falsch gesetzt, ggf. korrigieren
        "Bounds": "IntronEnd",
        "Score.1": "HMMScore",
        "Score.2": "Str2Score",
        "Isotype": "Isotype",
        "Score.3": "IsotypeScore",
        "Note": "Note"
    })

    # Spalten explizit auswählen und benennen
    cols_required = [
        "Sequence", "tRNA#", "Begin", "End", "Type", "Anticodon",
        "IntronBegin", "IntronEnd", "CM", "HMMScore", "Str2Score",
        "Isotype", "IsotypeScore", "Note"
    ]
    df.columns = cols_required[:len(df.columns)]  # falls Note fehlt, kein Fehler

    # Superphylum aus Stem ableiten
    df["Superphylum"] = stem.split('_', 1)[0]

    # Numerische Spalten erzwingen
    for col in ["Begin", "End", "IsotypeScore", "CM", "HMMScore", "Str2Score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Intronverarbeitung
    df["IntronBegin"] = df["IntronBegin"].astype(str)
    df["IntronEnd"] = df["IntronEnd"].astype(str)
    df["IntronPairs"] = df.apply(_row_to_pairs, axis=1)
    df["IntronCount"] = df["IntronPairs"].apply(len)
    df["TotalIntronLen"] = df["IntronPairs"].apply(lambda ls: sum(e - s + 1 for s, e in ls))
    df["MatureLen"] = (df["End"] - df["Begin"]).abs() + 1 - df["TotalIntronLen"]
    df["Name"] = (stem)
    
    if DEBUG:
        show_df(df, f"Parsed: {stem}_main.txt")
    
    
    return df
    
################################################################################
# DETAIL Parser-Funktionen
################################################################################
def read_detail(path: Path) -> pd.DataFrame:
    #Liest *_detail.txt, egal wie viele Leerzeichen/Felder;
    #nutzt maxsplit, um alles Überzählige in die letzte Spalte ("Note") zu kippen.
    
    if not path.exists():
        print(f"Warning: {path} not found, returning empty DataFrame")
        return pd.DataFrame()
        
    with path.open() as f:
        # erste nicht-leere Zeile = Header
        for line in f:
            if line.strip() and not line.startswith("#"):
                header = re.split(r'\s+', line.strip())
                break

        ncols = len(header)
        rows = []
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            # maxsplit = ncols-1 → alles weitere landet gemeinsam in der letzten Spalte
            row = re.split(r'\s+', line.strip(), maxsplit=ncols-1)
            # auffüllen, falls Note-Spalte fehlt
            row += [""] * (ncols - len(row))
            rows.append(row)

    df = pd.DataFrame(rows, columns=header)
    return df

#Extraktion der Isotypen / Anticodon counts aus dem weirden Key-Value-Format
def parse_anticodon_counts_from_detail(detail_path: Path) -> pd.DataFrame:
    #Extrahiert aus dem Abschnitt 'Isotype / Anticodon Counts:' die Isotypen,
    #deren Anticodons und jeweilige Counts.
    #Gibt ein DataFrame mit Spalten: Type, Anticodon, Count 
    #Ausgangslage war das etwas weirde Key-Value-Paar-Format in der Detail.txt-Datei.
    
    if not detail_path.exists():
        return pd.DataFrame(columns=["Type", "Anticodon", "Count"])
        
    in_section = False
    rows = []

    with detail_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Isotype / Anticodon Counts:"):
                in_section = True
                continue

            if in_section:
                if re.match(r"^\w+", line):  # beginnt mit Isotyp
                    tokens = re.split(r'\s+', line)
                    if len(tokens) < 2:
                        continue

                    isotype = tokens[0]
                    i = 2  # skip isotype + count
                    while i < len(tokens):
                        token = tokens[i]
                        if token.endswith(":"):
                            anticodon = token[:-1]
                            count = 0
                            if i + 1 < len(tokens) and not tokens[i+1].endswith(":"):
                                try:
                                    count = int(tokens[i+1])
                                    i += 1
                                except ValueError:
                                    pass
                            rows.append({
                                "Type": isotype,
                                "Anticodon": anticodon,
                                "Count": count
                            })
                        i += 1

    df = pd.DataFrame(rows)

    if DEBUG:
        show_df(df, f"Parsed: {stem}_detail.txt")

    return df

    

################################################################################
# STRUCT Parser-Funktionen
################################################################################
def read_struct(path: Path) -> pd.DataFrame:
    #Parst *_struct.txt als Blockstruktur mit Sequenz, Struktur, Score usw.
    #Gibt ein DataFrame mit einer Zeile pro tRNA zurück.
    
    if not path.is_file():
        print(f"Warning: {path} not found, returning empty DataFrame")
        return pd.DataFrame()

    records = []
    current = {}

    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Start eines neuen Blocks: >Sequence.tRNAx (Begin-End)
            if re.match(r"^\S+\.trna\d+ \(\d+-\d+\)\s+Length:", line):
                # vorherigen Datensatz speichern (falls vorhanden)
                if current:
                    records.append(current)
                    current = {}

                m = re.match(r"^(\S+)\.trna(\d+) \((\d+)-(\d+)\)\s+Length:\s+(\d+)", line)
                if m:
                    current["Sequence"] = m.group(1)
                    current["tRNA#"] = int(m.group(2))
                    current["Begin"] = int(m.group(3))
                    current["End"] = int(m.group(4))
                    current["Length"] = int(m.group(5))
                continue

            # Type, Anticodon, Score
            if line.startswith("Type:"):
                m = re.match(r"Type:\s+(\S+)\s+Anticodon:\s+(\S+)\s+at\s+[\d-]+\s+\(\d+-\d+\)\s+Score:\s+([0-9.]+)", line)
                if m:
                    current["Type"] = m.group(1)
                    current["Anticodon"] = m.group(2)
                    current["Score"] = float(m.group(3))
                continue

            # Possible intron-Zeilen
            if line.startswith("Possible intron:"):
                intron = re.findall(r"Possible intron:\s+([\d-]+)", line)
                if "Introns" not in current:
                    current["Introns"] = []
                current["Introns"].extend(intron)
                continue

            # Scorezeile: HMM Sc, Str2Sc
            if "HMM Sc=" in line and "Sec struct Sc=" in line:
                m = re.search(r"HMM Sc=([0-9.]+)", line)
                if m:
                    current["HMMScore"] = float(m.group(1))
                m = re.search(r"Sec struct Sc=([0-9.]+)", line)
                if m:
                    current["Str2Score"] = float(m.group(1))
                continue

            # Strukturzeilen (Seq / Str / Pre)
            if line.startswith("Seq:"):
                current["Seq"] = line.replace("Seq:", "").strip()
                continue
            if line.startswith("Str:"):
                current["Str"] = line.replace("Str:", "").strip()
                continue
            if line.startswith("Pre:"):
                current["Pre"] = line.replace("Pre:", "").strip()
                continue

    # Letzten Block speichern
    if current:
        records.append(current)

    df = pd.DataFrame(records)

    # Dtypes setzen
    DTYPES_STRUCT = {
        "Sequence": "category",
        "tRNA#": "int16",
        "Begin": "int32",
        "End": "int32",
        "Type": "category",
        "Anticodon": "category",
        "Score": "float32",
        "HMMScore": "float32",
        "Str2Score": "float32",
        "Length": "int16",
        "Seq": "string",
        "Str": "string",
        "Pre": "string"
    }
    df = df.astype({k: v for k, v in DTYPES_STRUCT.items() if k in df.columns})
    
    if "Str" in df.columns:
        df["MatureLenStr"] = df["Str"].astype(str).apply(len)

    if DEBUG:
        show_df(df, f"Parsed: {stem}_struct.txt")

    return df
    
        #try:
    #    struct = read_struct(struct_file)
    #except Exception as e:
    #    print(f"  Fehler beim Lesen von {struct_file}: {e}")
    #    struct = pd.DataFrame()

    

################################################################################
# SOME DEF BEFORE MERGE
#################################################################################
def _parse_intron_field(field: str | int | float) -> list[int]:
    #Wandelt das Feld 'IntronBegin' oder 'IntronEnd' in eine Liste von
    #Integern um.  Beispiele:
    #    "37,58"  -> [37, 58]
    #    "0"      -> []
    #    0 (int)  -> []
    
    if pd.isna(field) or str(field).strip() in {"", "0", "-"}:
        return []
    return [int(x) for x in str(field).split(",") if x.strip().isdigit()]

def _row_to_pairs(row) -> list[tuple[int, int]]:
    #Erzeugt aus den Spalten IntronBegin / IntronEnd eine Liste von (start,end)-Tupeln.
    #Beispiel:
    #    IntronBegin = "37,58"
    #    IntronEnd   = "38,60"
    #    --> [(37,38), (58,60)]
    
    begins = _parse_intron_field(row["IntronBegin"])
    ends   = _parse_intron_field(row["IntronEnd"])
    return list(zip(begins, ends))

def calculate_avg_length_summary(df: pd.DataFrame) -> pd.DataFrame:
    #Berechnet den Durchschnitt der reifen tRNA-Längen (MatureLen)
    #über das gesamte DataFrame aka ohne Introns. 
    #Gibt ein DataFrame mit den Spalten: Superphylum, Type, AvgLength zurück.
    #updated! das läuft über die Spalte MatureLenStr (Zählung der Zeichen in Str aus Struct)
    # ne. läuft wieder über MatureLen
    if "Type" not in df.columns:
        df["Type"] = pd.NA  # Fallback

    avg = df["MatureLen"].mean()
    result = pd.DataFrame({
        "Superphylum": [df["Superphylum"].iloc[0] if "Superphylum" in df.columns else "Unknown"],
        "Type": [",".join(sorted(df["Type"].dropna().unique())) if "Type" in df.columns else "NA"],
        "AvgLength": [avg]
    })
    return result 




#################################################################################
# SINGLE / BATCH PROCESSING
# war wichtig zum Herantasten ans Parsen, ist noch drin für future adventures, aber spielt im weiteren verlauf keene rolle
#################################################################################
def process_single_organism(stem: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    #Verarbeitet einen einzelnen Organismus (stem) und gibt die kombinierten
    #Daten sowie die Durchschnittslängen-Zusammenfassung zurück.
    
    print(f"\nVerarbeite Organismus: {stem}")
    
    # Dateipfade definieren
    main_file = OUT_DIR / f"{stem}_main.txt"
    detail_file = OUT_DIR / f"{stem}_detail.txt"
    struct_file = OUT_DIR / f"{stem}_struct.txt"
    
    #Überprüfen, ob alle drei Dateien existieren
    missing_files = []
    if not main_file.exists():
        missing_files.append("main")
    if not detail_file.exists():
        missing_files.append("detail")
    if not struct_file.exists():
        missing_files.append("struct")
    
    if missing_files:
        print(f"  Warnung: Fehlende Dateien für {stem}: {', '.join(missing_files)}")
        return pd.DataFrame(), pd.DataFrame()
    
    try:
        main = parse_main_file(stem, OUT_DIR)
    except Exception as e:
        print(f"  Fehler beim Parsen der Main-Datei: {e}")
        return pd.DataFrame(), pd.DataFrame()

    try:
        struct = read_struct(struct_file)
    except Exception as e:
        print(f"  Fehler beim Lesen der Struct-Datei: {e}")
        struct = pd.DataFrame()

    # Anticodon counts parsen
    try:
        anticodon_counts = parse_anticodon_counts_from_detail(detail_file)
    except Exception as e:
        print(f"  Fehler beim Parsen der Anticodon-Counts: {e}")
        anticodon_counts = pd.DataFrame()
    
    # Intron-Analyse und Mature-Length-Berechnung
    try:
        # beide Spalten erst einmal als String behandeln, damit Split funktioniert
        main["IntronBegin"] = main["IntronBegin"].astype(str)
        main["IntronEnd"]   = main["IntronEnd"].astype(str)
        
        # Liste der (start,end)-Paare
        main["IntronPairs"] = main.apply(_row_to_pairs, axis=1)
        
        # Anzahl Introns
        main["IntronCount"] = main["IntronPairs"].apply(len)
        
        # Gesamtlänge aller Introns
        main["TotalIntronLen"] = main["IntronPairs"].apply(
            lambda ls: sum(e - s + 1 for s, e in ls)
        )
        
        # reife (mature) tRNA-Länge: |Begin-End| + 1  –  Summe Intronlängen
        main["MatureLen"] = (
            (main["End"].astype(int) - main["Begin"].astype(int)).abs() + 1
            - main["TotalIntronLen"]
        )
        
    except Exception as e:
        print(f"  Fehler bei der Intron-Analyse: {e}")
        return pd.DataFrame(), pd.DataFrame()


###############################################################################################################
#MERGE! MAIN-STRUCT
################################################################################################################
    print(f"Zeilen in main: {len(main)}")
    print(f"Zeilen in struct: {len(struct)}")
    
    #checkup, damit 
    main["tRNA#"] = main["tRNA#"].astype(int)
    struct["tRNA#"] = struct["tRNA#"].astype(int)

    if not struct.empty:
        try:
            mainstruct = pd.merge(
                main,
                struct[["Sequence", "tRNA#", "Str", "Seq"]],  #seq hinzugefügt
                on=["Sequence", "tRNA#"], 
                how="left"  # linke Seite = main, damit keine Zeilen verloren gehen
            )
            
            #Anheften der Spalte MatureLenStr, in welcher die (mature) nts gezählt werden aus der Anzahl an Zeichen in Str
            if "Str" in mainstruct.columns:
                mainstruct["MatureLenStr"] = mainstruct["Str"].astype(str).apply(len)

            #Anheften der Spalte CCA_Check; prüft, ob Str am Ende CCA traegt, gibt 1 oder 0 aus
            mainstruct["CCA_Check"] = mainstruct["Seq"].astype(str).apply(
                    lambda x: 1 if x.endswith("CCA") or x.endswith("CCA.") else 0
                )

            #Anheften der Spalte Species_ID, besteht nur aus der ID, ohne SUperphylum
            def extract_species_id(name):
                match = re.search(r'(\d.*\d)', name)
                return match.group(1) if match else None

            mainstruct["Species_ID"] = mainstruct["Name"].apply(extract_species_id)

            print(f"Zeilen in mainstruct: {len(mainstruct)}")
            assert len(mainstruct) == len(main), "Merge hat Zeilen vervielfacht!"


        except Exception as e:
            print(f"  Fehler beim Mergen von Main und Struct: {e}")
            mainstruct = main.copy()
    else:
        mainstruct = main.copy()
        mainstruct["Str"] = pd.NA  # Fallback-Spalte einfügen

    if "Str" in mainstruct.columns:
        mainstruct["StrParse"] = mainstruct["Str"].astype(str).apply(
            lambda s: s.replace(">", "(").replace("<", ")")
        )

    try:
        export_path = OUT_DIR / f"{stem}_mainstruct.txt"
        mainstruct.to_csv(export_path, sep="\t", index=False)
        print(f"  Exportiert: {export_path}")
    except Exception as e:
        print(f"  Fehler beim Export: {e}")
    
    # Durchschnittslängen-Zusammenfassung berechnen
    try:
        avg_length_summary = calculate_avg_length_summary(mainstruct)
        print(f"  Durchschnittliche Mature-Length: {avg_length_summary['AvgLength'].iloc[0]:.2f}")
    except Exception as e:
        print(f"  Fehler bei der Durchschnittsberechnung: {e}")
        avg_length_summary = pd.DataFrame()
    
    return mainstruct, avg_length_summary

###############################################################################
# HAUPTLOGIK: Entscheidet zwischen Einzel- und Batch-Verarbeitung 
# ->> aka the genesis of ALL_mainstruct.txt
#################################################################################

def main():
    #Hauptfunktion: Entscheidet basierend auf der ALL-Variable,
    #ob ein einzelner Organismus oder alle Organismen verarbeitet werden.
        
    if ALL:
        print("BATCH-MODUS: Verarbeite alle Organismen")
        process_all_organisms()
    else:
        print(f"EINZEL-MODUS: Verarbeite nur {SINGLE_STEM}")
        process_single_organism_interactive(SINGLE_STEM)

def process_all_organisms():
    #Verarbeitet alle Organismen im Batch-Modus (ursprüngliche main-Funktion)
    
    print("Suche nach verfügbaren Organismen...")
    stems = get_unique_stems(OUT_DIR)
    
    if not stems:
        print("Keine *_main.txt Dateien gefunden!")
        return
    
    print(f"Gefundene Organismen: {len(stems)}")
    for i, stem in enumerate(stems[:5]):  # Zeige die ersten 5
        print(f"  {i+1}. {stem}")
    if len(stems) > 5:
        print(f"  ... und {len(stems) - 5} weitere")
    
    # Listen für die Sammlung aller Ergebnisse
    all_mainstruct = []
    all_summaries = []
    
    # Schleife über alle Organismen
    for i, stem in enumerate(stems, 1):
        print(f"\n[{i}/{len(stems)}] Verarbeite: {stem}")
        
        try:
            mainstruct, summary = process_single_organism(stem)
            
            if not mainstruct.empty:
                all_mainstruct.append(mainstruct)
                print(f"{len(mainstruct)} tRNAs verarbeitet")
            
            if not summary.empty:
                all_summaries.append(summary)
                
        except Exception as e:
            print(f"Fehler bei {stem}: {e}")
            continue
    
    # Alle Ergebnisse kombinieren
    if all_mainstruct:
        print(f"\nKombiniere Ergebnisse von {len(all_mainstruct)} Organismen...")
        combined_mainstruct = pd.concat(all_mainstruct, ignore_index=True)
        
        # Export der kombinierten Daten
        combined_export_path = OUT_DIR / "ALL_combined_mainstruct.txt"
        combined_mainstruct.to_csv(combined_export_path, sep="\t", index=False)
        print(f"Alle kombinierten Daten exportiert: {combined_export_path}")
        print(f"Gesamt-tRNAs: {len(combined_mainstruct)}")
    
    if all_summaries:
        combined_summaries = pd.concat(all_summaries, ignore_index=True)
        
        # Export der Zusammenfassungen
        summary_export_path = OUT_DIR / "ALL_length_summaries.txt"
        combined_summaries.to_csv(summary_export_path, sep="\t", index=False)
        print(f"Alle Längen-Zusammenfassungen exportiert: {summary_export_path}")
        
        # Finale Statistiken
        print("\nAll them final statistics for overview purposes... because data <3")
        print(f"Durchschnittliche Mature-Length über alle Organismen: {combined_summaries['AvgLength'].mean():.2f}")
        print(f"Min: {combined_summaries['AvgLength'].min():.2f}")
        print(f"Max: {combined_summaries['AvgLength'].max():.2f}")
        print(f"Standardabweichung: {combined_summaries['AvgLength'].std():.2f}")
        if "IntronCount" in mainstruct.columns:
            avg_introns = mainstruct["IntronCount"].mean()
            print(f"Durchschnittliche Anzahl Introns pro tRNA: {avg_introns:.2f}")
        

        # Zeige Zusammenfassung (je nach DEBUG-Setting)
        show_df(combined_summaries, "Alle Längen-Zusammenfassungen")
    
    print("\nBatch-Verarbeitung yes check")

def process_single_organism_interactive(stem: str):
    #Verarbeitet einen einzelnen Organismus im interaktiven Modus 
    #(entspricht dem ursprünglichen Verhalten des Codes)
    
    print(f"Verarbeite Organismus: {stem}")
    
    # Verarbeitung durchführen
    mainstruct, summary = process_single_organism(stem)
    
    if mainstruct.empty:
        print(f"Keine Daten für {stem} gefunden oder Fehler beim Verarbeiten.")
        return
    
    print(f"{len(mainstruct)} tRNAs erfolgreich verarbeitet")
    
    # Im Einzelmodus: Detaillierte Anzeige der Zwischenergebnisse
    if DEBUG:
        print("\nDETAILLIERTE ERGEBNISSE")
        
        # Zeige die verschiedenen DataFrames wie im ursprünglichen Code
        # (Hier können Sie die original show_df Aufrufe aus Ihrem Code einfügen)
        show_df(mainstruct, f"Main + Struct kombiniert für {stem}")
        
        if not summary.empty:
            show_df(summary, f"Durchschnittslängen-Zusammenfassung für {stem}")
    
    # Zusammenfassung ausgeben
    if not summary.empty:
        avg_length = summary['AvgLength'].iloc[0]
        superphylum = summary['Superphylum'].iloc[0]
        print(f"\nUSAMMENFASSUNG für {stem}")
        print(f"Superphylum: {superphylum}")
        print(f"Anzahl tRNAs: {len(mainstruct)}")
        print(f"Durchschnittliche Mature-Length: {avg_length:.2f}")
        
        # Zusätzliche Statistiken für den Einzelmodus
        print(f"Median Mature-Length: {mainstruct['MatureLen'].median():.2f}")
        print(f"Min Mature-Length: {mainstruct['MatureLen'].min()}")
        print(f"Max Mature-Length: {mainstruct['MatureLen'].max()}")
        print(f"Anzahl verschiedener tRNA-Typen: {mainstruct['Type'].nunique()}")
        
        # Häufigste tRNA-Typen anzeigen
        type_counts = mainstruct['Type'].value_counts().head()
        print(f"\nHäufigste tRNA-Typen:")
        for trna_type, count in type_counts.items():
            print(f"  {trna_type}: {count}")

    print("\nBatch-Verarbeitung abgeschlossen")

    if DEBUG:
        show_df(mainstruct, "Mainstruct: Gesamte Ausgabedatei")

# Hauptausführung
if __name__ == "__main__":
    main()