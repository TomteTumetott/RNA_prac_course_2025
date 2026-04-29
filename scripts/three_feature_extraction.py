import pandas as pd
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np

# config einbinden
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR

# Struktur-Dekodierer
def extract_structural_features(dot_bracket_structure):
    #Extrahiert Features aus Dot-Bracket-Notation:
    #- creates short form on the base of symbol changes
    #- zählt die Anzahl an Stems ('(')
    #- gibt Positionen der Wechsel zurück
    
    if not dot_bracket_structure or not isinstance(dot_bracket_structure, str):
        return "", 0, []

    structure = dot_bracket_structure.strip()
    if not structure:
        return "", 0, []

    feature_string = ""
    feature_positions = []
    current_symbol = structure[0]
    feature_number = 1
    stem_count = 0

    # Initialsymbol
    feature_string += current_symbol
    if current_symbol == '(':
        stem_count += 1
    feature_positions.append((0, current_symbol, feature_number))

    for i in range(1, len(structure)):
        symbol = structure[i]
        if symbol != current_symbol:
            feature_number += 1
            current_symbol = symbol
            feature_string += symbol
            feature_positions.append((i, symbol, feature_number))
            if symbol == '(':
                stem_count += 1

    return feature_string, stem_count, feature_positions


def analyze_mainstruct_dotbracket(input_file: Path, output_file: Path = None):
    #Liest eine ALL_combined_mainstruct.txt Datei, analysiert Dot-Bracket-Strukturspalte (StrParse),
    #und hängt Spalten mit Feature-String und Stem-Zahl an.
    #
    if not input_file.exists():
        print(f"[FEHLER] Datei nicht gefunden: {input_file}")
        return

    df = pd.read_csv(input_file, sep="\t")
    if "StrParse" not in df.columns:
        print("[FEHLER] Erwartete Spalte 'StrParse' nicht gefunden.")
        return

    features = []
    stems = []

    for structure in df["StrParse"]:
        fstr, scount, _ = extract_structural_features(structure)
        features.append(fstr)
        stems.append(scount)

    df["StructureShort"] = features
    df["StemCount"] = stems

    print(f"[INFO] Struktur-Features berechnet für {len(df)} Einträge")

    #neue Spalte an die Datei anhängen
    if "Type" not in df.columns or "Anticodon" not in df.columns:
        print("[WARNUNG] Spalten 'Type' oder 'Anticodon' nicht gefunden.")
    
    features = []
    stems = []

    for structure in df["StrParse"]:
        fstr, scount, _ = extract_structural_features(structure)
        features.append(fstr)
        stems.append(scount)

    df["StructureShort"] = features
    df["StemCount"] = stems
    
    # Neue Spalte: Type und Anticodon kombiniert
    if "Type" in df.columns and "Anticodon" in df.columns:
        df["Type_Anticodon"] = df["Type"].astype(str) + "_" + df["Anticodon"].astype(str)
        print(f"[INFO] Spalte 'Type_Anticodon' hinzugefügt")

    #Anheften der Spalte Species_ID, besteht nur aus der ID, ohne SUperphylum
    def extract_species_id(name):
        match = re.search(r'(\d.*\d)', name)
        return match.group(1) if match else None

        mainstruct["ID"] = mainstruct["Name"].apply(extract_species_id)

    #speichern
    if output_file:
        df.to_csv(output_file, sep="\t", index=False)
        print(f"[INFO] Erweiterte Datei exportiert: {output_file}")

    return df



if __name__ == "__main__":
    input_path = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct.txt"
    output_path = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct_annotated.txt"
    noncanonical_path = TRNASCAN_OUTPUT_DIR / "ALL_noncanonical_mainstruct.txt"

    df_result = analyze_mainstruct_dotbracket(input_path, output_path)

# ausgabe der combined mainstruct Datei zur Übersicht
    if df_result is not None and output_path.exists():
        # Headerliste in Konsole ausgeben
        print("[INFO] Spaltenübersicht:")
        print(df_result.columns.to_list())

        # Plot der Stem‑Count‑Verteilung anzeigen
        values = df_result["StemCount"]

        # Bin-Kanten so wählen, dass jeder Integer zentriert ist
        min_val, max_val = values.min(), values.max()
        bin_edges = np.arange(min_val - 0.5, max_val + 1.5, 1)   #  …-0.5, 0.5, 1.5, …

        plt.figure(figsize=(8, 5))
        plt.hist(values, bins=bin_edges, edgecolor="black")       # align='mid' ist jetzt implizit
        plt.xticks(range(min_val, max_val + 1))                   # Ticks genau unter den Mittelpunkten
        plt.title("Distribution of Stem Counts")
        plt.xlabel("Stem Counts")
        plt.ylabel("Counts")
        plt.tight_layout()
        plt.show()

 
        if df_result is not None:
            # Kopfzeile zeigen
            print("[INFO] Spaltenübersicht:")
            print(df_result.columns.to_list())

            # Histogram (zentriert)
            values = df_result["StemCount"]
            min_val, max_val = values.min(), values.max()
            bin_edges = np.arange(min_val - 0.5, max_val + 1.5, 1)

            plt.figure(figsize=(8, 5))
            plt.hist(values, bins=bin_edges, edgecolor="black")
            plt.xticks(range(min_val, max_val + 1))
            plt.title("Verteilung der Stem-Anzahl")
            plt.xlabel("Stem-Anzahl")
            plt.ylabel("Häufigkeit")
            plt.tight_layout()
            plt.show()

            #   Nicht‑kanonische tRNAs 
            #   (≠ 4- oder 5-Stem‑Strukturen)
            # -----------------------------
            df_noncanonical = df_result[~df_result["StemCount"].isin([4, 5])]
            df_noncanonical.to_csv(noncanonical_path, sep="\t", index=False)

            print(
                f"[INFO] Non‑canonical Dataset gespeichert: {noncanonical_path} "
                f"({len(df_noncanonical)} Zeilen)"
            )
      
#     