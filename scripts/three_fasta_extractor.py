############hier entstehen Duplikate wie damit umgehen? 

#!/usr/bin/env python3
"""
tRNA FASTA Exporter für beste Kandidaten
Extrahiert Sequenzen und Sekundärstrukturen direkt aus ALL_best_candidates.txt
für mlocarna-Analyse

DATENQUELLE: ALL_best_candidates.txt (enthält bereits Sequenzen und Strukturen)
ZIELFORMAT: Extended FASTA mit Sekundärstrukturen (#FS Tag für mlocarna)
ORGANISATION: Nach Anticodon UND Struktur-Typ
"""

import pandas as pd
from pathlib import Path
import sys
import re
from collections import defaultdict
from typing import Dict, List, Tuple

# config einbinden
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, OUTPUT_DIR

# Konfiguration für Kandidaten-Ordner
CANDIDATES_DIR = TRNASCAN_OUTPUT_DIR / "candidates"

class tRNAFastaExporter:
    #Exportiert tRNA-Sequenzen und -Strukturen aus ALL_best_candidates.txt im mlocarna-Format
        
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.candidates_dir = output_dir / "candidates"
        self.candidates_dir.mkdir(exist_ok=True)
        
        print(f"[INFO] Kandidaten-Ordner: {self.candidates_dir}")
    
    def extract_sequence_and_structure(self, row: pd.Series) -> Tuple[str, str]:
        #Extrahiert Sequenz und Struktur aus einer Zeile von ALL_best_candidates.txt
        
        # Sequenz extrahieren - überprüfe verschiedene Spalten
        sequence = ""
        seq_source = ""
        
        # Mögliche Sequenz-Spalten in Prioritätsreihenfolge
        seq_columns = ['Seq', 'seq', 'sequence', 'Sequence']
        for col in seq_columns:
            if col in row and pd.notna(row[col]):
                raw_seq = str(row[col]).strip()
                if raw_seq and raw_seq != '' and 'N' not in raw_seq[:10]:  # Prüfe ersten 10 Zeichen auf N's
                    # Bereinige Sequenz
                    clean_seq = raw_seq.upper().replace('U', 'T').replace('-', '').replace(' ', '')
                    # Validiere dass es DNA/RNA Sequenz ist
                    valid_chars = set('ATCGN')
                    if all(c in valid_chars for c in clean_seq) and len(clean_seq) >= 50:
                        sequence = clean_seq
                        seq_source = col
                        break
        
        # Struktur extrahieren (bevorzuge StrParse, dann Str)
        structure = ""
        struct_source = ""
        
        struct_columns = ['StrParse', 'Str', 'str', 'structure']
        for col in struct_columns:
            if col in row and pd.notna(row[col]):
                raw_struct = str(row[col]).strip()
                if raw_struct and raw_struct != '':
                    # Bereinige Struktur
                    clean_struct = raw_struct.replace('>', ')').replace('<', '(').replace(' ', '')
                    # Validiere dass es Dot-Bracket Notation ist
                    valid_chars = set('().abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ[]{}')
                    if all(c in valid_chars for c in clean_struct):
                        structure = clean_struct
                        struct_source = col
                        break
        
        # Debugging Info für problematische Einträge
        if not sequence or len(sequence) < 50:
            print(f"  [WARNUNG] Keine gültige Sequenz gefunden für {row.get('Sequence', 'unknown')}.tRNA{row.get('tRNA#_x', '?')}")
            print(f"    Verfügbare Spalten: {[col for col in seq_columns if col in row]}")
            for col in seq_columns:
                if col in row and pd.notna(row[col]):
                    preview = str(row[col])[:50] + "..." if len(str(row[col])) > 50 else str(row[col])
                    print(f"    {col}: '{preview}'")
            
            # Letzter Versuch: Falls genomische Koordinaten verfügbar sind
            if not sequence:
                print(f"  [WARNUNG] Verwende Platzhalter-Sequenz für {row.get('Sequence', 'unknown')}.tRNA{row.get('tRNA#_x', '?')}")
                length = row.get('MatureLen', row.get('Length', 75))
                sequence = 'N' * int(length)
        else:
            print(f"  [OK] Sequenz aus Spalte '{seq_source}': {len(sequence)} nt")
        
        if not structure:
            print(f"  [WARNUNG] Keine gültige Struktur gefunden für {row.get('Sequence', 'unknown')}.tRNA{row.get('tRNA#_x', '?')}")
            structure = '.' * len(sequence)
        else:
            print(f"  [OK] Struktur aus Spalte '{struct_source}': {len(structure)} Zeichen")
        
        # Struktur-Länge an Sequenz-Länge anpassen
        if len(structure) != len(sequence):
            if len(structure) > len(sequence):
                structure = structure[:len(sequence)]
            else:
                structure = structure + '.' * (len(sequence) - len(structure))
        
        return sequence, structure
    
    def get_structure_type(self, structure: str) -> str:
        #Bestimmt den Struktur-Typ basierend auf der Sekundärstruktur
        #Falls bereits in den Daten verfügbar, nutze die vorhandene Information
        
        if not structure or pd.isna(structure):
            return "unknown"
        
        # Zähle Stems (Anzahl öffnender Klammern)
        stem_count = structure.count('(')
        
        # Klassifiziere nach Stem-Anzahl
        if stem_count == 4:
            return "4stem"
        elif stem_count == 5:
            return "5stem"
        else:
            return f"{stem_count}stem"
    
    def group_by_anticodon_and_structure(self, df: pd.DataFrame) -> Dict[str, Dict[str, pd.DataFrame]]:
        #Gruppiert tRNAs nach Anticodon UND Struktur-Typ
        
        print(f"[INFO] Gruppiere {len(df)} Kandidaten nach Anticodon und Struktur-Typ...")
        
        grouped = defaultdict(lambda: defaultdict(list))
        
        for idx, row in df.iterrows():
            anticodon = row.get('Anticodon', 'Unknown')
            
            # Bestimme Struktur-Typ
            if 'StemCount' in row and pd.notna(row['StemCount']):
                stem_count = int(row['StemCount'])
                if stem_count == 4:
                    structure_type = "4stem"
                elif stem_count == 5:
                    structure_type = "5stem"
                else:
                    structure_type = f"{stem_count}stem"
            else:
                # Fallback: aus Struktur ableiten
                _, structure = self.extract_sequence_and_structure(row)
                structure_type = self.get_structure_type(structure)
            
            grouped[anticodon][structure_type].append(row)
        
        # Konvertiere Listen zu DataFrames
        result = {}
        for anticodon, structure_data in grouped.items():
            result[anticodon] = {}
            for structure_type, rows in structure_data.items():
                result[anticodon][structure_type] = pd.DataFrame(rows)
        
        return result
    
    def try_extract_from_struct_files(self, df: pd.DataFrame) -> pd.DataFrame:
        #Fallback: Versuche Sequenzen aus originalen struct-Dateien zu extrahieren
        
        print(f"[FALLBACK] Extrahiere Sequenzen aus *_struct.txt Dateien...")
        
        # Finde verfügbare struct-Dateien
        struct_files = list(self.output_dir.glob("*_struct.txt"))
        print(f"[INFO] Gefunden: {len(struct_files)} struct-Dateien")
        
        if len(struct_files) == 0:
            print("[WARNUNG] Keine struct-Dateien gefunden!")
            return df
        
        # Lade struct-Daten
        struct_data_combined = []
        
        for struct_file in struct_files:
            try:
                struct_records = self.parse_struct_file(struct_file)
                if len(struct_records) > 0:
                    struct_data_combined.extend(struct_records)
                    print(f"  {struct_file.name}: {len(struct_records)} tRNAs geladen")
            except Exception as e:
                print(f"  [FEHLER] {struct_file.name}: {e}")
        
        if len(struct_data_combined) == 0:
            print("[WARNUNG] Keine Daten aus struct-Dateien extrahiert!")
            return df
        
        struct_df = pd.DataFrame(struct_data_combined)
        print(f"[INFO] Insgesamt {len(struct_df)} tRNAs aus struct-Dateien")
        
        # Merge mit besten Kandidaten
        enhanced_df = df.copy()
        sequences_added = 0
        structures_added = 0
        
        for idx, row in enhanced_df.iterrows():
            # Suche entsprechende tRNA in struct-Daten
            matches = struct_df[
                (struct_df['Sequence'] == row['Sequence']) &
                (struct_df['tRNA#'] == row.get('tRNA#_x', row.get('tRNA#_y', 1)))
            ]
            
            if len(matches) == 0:
                # Fallback: über Koordinaten
                matches = struct_df[
                    (struct_df['Sequence'] == row['Sequence']) &
                    (struct_df['Begin'] == row['Begin']) &
                    (struct_df['End'] == row['End'])
                ]
            
            if len(matches) > 0:
                match = matches.iloc[0]
                
                # Füge Sequenz hinzu falls nicht vorhanden
                if 'Seq' not in enhanced_df.columns or pd.isna(enhanced_df.at[idx, 'Seq']) or enhanced_df.at[idx, 'Seq'] == '':
                    if 'Seq' in match and pd.notna(match['Seq']):
                        enhanced_df.at[idx, 'Seq'] = match['Seq']
                        sequences_added += 1
                
                # Füge Struktur hinzu falls nicht vorhanden
                if 'Str' not in enhanced_df.columns or pd.isna(enhanced_df.at[idx, 'Str']) or enhanced_df.at[idx, 'Str'] == '':
                    if 'Str' in match and pd.notna(match['Str']):
                        enhanced_df.at[idx, 'Str'] = match['Str']
                        structures_added += 1
        
        print(f"[SUCCESS] {sequences_added} Sequenzen und {structures_added} Strukturen hinzugefügt")
        return enhanced_df
    
    def parse_struct_file(self, struct_file: Path) -> List[Dict]:
        #Parst tRNAscan-SE *_struct.txt Datei
        
        records = []
        current = {}

        with struct_file.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Start eines neuen Blocks
                if re.match(r"^\S+\.trna\d+ \(\d+-\d+\)\s+Length:", line):
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

                # Sequenz und Struktur
                if line.startswith("Seq:"):
                    current["Seq"] = line.replace("Seq:", "").strip()
                    continue
                if line.startswith("Str:"):
                    current["Str"] = line.replace("Str:", "").strip()
                    continue

        if current:
            records.append(current)

        return records
   
    def create_fasta_header(self, row: pd.Series) -> str:
        """
        Erstellt FASTA-Header für tRNA
        """
        superphylum = row["Superphylum"]
        species_id = row["Species_ID"]
        trna_type = row.get('Type_Anticodon')
        sequence_id = row.get('Sequence')
        trna_num = row.get('tRNA#_x', row.get('tRNA#_y', row.get('tRNA#', 1)))
        begin = row.get('Begin', 0)
        end = row.get('End', 0)
        score = row.get('IsotypeScore')
       #length = f"{row['value']}bp"
        return f">{superphylum}__{species_id}__{trna_type}__{sequence_id}.trna{trna_num}"
        return header



    def format_structure_for_mlocarna(self, structure: str) -> str:
        #Formatiert Struktur für mlocarna mit #FS Tag (Fixed Structure)
        
        if not structure or pd.isna(structure):
            return " #FS"
        
        # Bereinige Struktur
        clean_structure = str(structure).strip()
        
        # Normalisiere Klammern
        clean_structure = clean_structure.replace('>', ')').replace('<', '(')
        
        # Füge #FS Tag hinzu (Fixed Structure für mlocarna)
        return f"{clean_structure} #FS"
    
    def export_anticodon_structure_fasta(self, anticodon: str, structure_type: str, df: pd.DataFrame) -> Path:
        #Exportiert tRNAs eines Anticodons und Struktur-Typs in FASTA
        
        # Sicherer Dateiname
        safe_anticodon = re.sub(r'[^A-Za-z0-9]', '_', anticodon)
        safe_structure = re.sub(r'[^A-Za-z0-9]', '_', structure_type)
        
        fasta_file = self.candidates_dir / f"tRNA_{safe_anticodon}_{safe_structure}_candidates.fasta"
        
        sequences_written = 0
        
        with open(fasta_file, 'w') as f:
            # Header-Kommentar
            f.write(f"# tRNA candidates: Anticodon {anticodon}, Structure {structure_type}\n")
            f.write(f"# Total candidates: {len(df)}\n")
            f.write(f"# Source: ALL_best_candidates.txt (filtered best candidates)\n")
            f.write(f"# Format: Extended FASTA for mlocarna (Fixed Structure #FS)\n\n")
            
            for _, row in df.iterrows():
                try:
                    # Extrahiere Sequenz und Struktur
                    sequence, structure = self.extract_sequence_and_structure(row)
                    
                    if not sequence or len(sequence) < 30:  # Mindestlänge für tRNA
                        print(f"  [WARNUNG] Sequenz zu kurz oder leer für {row.get('Sequence', 'unknown')}.tRNA{row.get('tRNA#_x', '?')}")
                        continue
                    
                    # Erstelle FASTA-Eintrag
                    header = self.create_fasta_header(row)
                    formatted_structure = self.format_structure_for_mlocarna(structure)
                    
                    # Schreibe Extended FASTA Format
                    f.write(f"{header}\n")
                    f.write(f"{sequence}\n")
                    f.write(f"{formatted_structure}\n\n")
                    
                    sequences_written += 1
                    
                except Exception as e:
                    print(f"  [FEHLER] Beim Exportieren von {row.get('Sequence', 'unknown')}: {e}")
                    continue
        
        print(f"  {sequences_written} Sequenzen → {fasta_file.name}")
        return fasta_file
    
    def validate_data_completeness(self, df: pd.DataFrame) -> Dict[str, any]:
        #Validiert die Vollständigkeit der Daten in ALL_best_candidates.txt
        
        print(f"[VALIDATION] Überprüfe Datenvollständigkeit...")
        
        validation = {
            'total_rows': len(df),
            'has_sequence': False,
            'has_structure': False,
            'sequence_sources': [],
            'structure_sources': [],
            'valid_sequences': 0,
            'valid_structures': 0
        }
        
        # Überprüfe verfügbare Spalten
        available_columns = list(df.columns)
        print(f"[INFO] Verfügbare Spalten: {available_columns}")
        
        # Sequenz-Quellen
        seq_columns = [col for col in ['Seq', 'seq', 'sequence', 'Sequence'] if col in available_columns]
        if seq_columns:
            validation['sequence_sources'] = seq_columns
            
            # Prüfe Qualität der Sequenzen
            for col in seq_columns:
                non_null = df[col].notna().sum()
                non_empty = df[col].apply(lambda x: str(x).strip() != '' if pd.notna(x) else False).sum()
                
                # Prüfe auf echte DNA/RNA Sequenzen (nicht nur N's)
                valid_seqs = 0
                placeholder_seqs = 0
                
                for idx, val in df[col].items():
                    if pd.notna(val) and str(val).strip():
                        seq = str(val).upper().replace('U', 'T').replace('-', '').replace(' ', '')
                        if len(seq) >= 50:
                            # Prüfe ob es hauptsächlich N's sind (Platzhalter)
                            n_count = seq.count('N')
                            if n_count / len(seq) > 0.5:  # Mehr als 50% N's
                                placeholder_seqs += 1
                            else:
                                valid_seqs += 1
                
                print(f"  Sequenz-Spalte '{col}':")
                print(f"    Nicht-null: {non_null}, Nicht-leer: {non_empty}")
                print(f"    Gültige Sequenzen: {valid_seqs}, Platzhalter: {placeholder_seqs}")
                
                if valid_seqs > 0:
                    validation['has_sequence'] = True
                    validation['valid_sequences'] = max(validation['valid_sequences'], valid_seqs)
        
        # Struktur-Quellen
        struct_columns = [col for col in ['Str', 'StrParse', 'str', 'structure'] if col in available_columns]
        if struct_columns:
            validation['structure_sources'] = struct_columns
            
            # Prüfe Vollständigkeit
            for col in struct_columns:
                non_null = df[col].notna().sum()
                non_empty = df[col].apply(lambda x: str(x).strip() != '' if pd.notna(x) else False).sum()
                
                # Prüfe auf gültige Strukturen
                valid_structs = 0
                for idx, val in df[col].items():
                    if pd.notna(val) and str(val).strip():
                        struct = str(val).strip()
                        # Einfache Validierung: enthält Klammern oder Punkte
                        if any(c in struct for c in '().'):
                            valid_structs += 1
                
                print(f"  Struktur-Spalte '{col}':")
                print(f"    Nicht-null: {non_null}, Nicht-leer: {non_empty}")
                print(f"    Gültige Strukturen: {valid_structs}")
                
                if valid_structs > 0:
                    validation['has_structure'] = True
                    validation['valid_structures'] = max(validation['valid_structures'], valid_structs)
        
        # Zeige Sample-Daten zur Diagnose
        print(f"\n[SAMPLE DATA] Erste 3 Zeilen zur Diagnose:")
        for idx in range(min(3, len(df))):
            row = df.iloc[idx]
            print(f"  Zeile {idx + 1}: {row.get('Sequence', 'unknown')}.tRNA{row.get('tRNA#_x', '?')}")
            
            # Zeige verfügbare Sequenz-Daten
            for col in seq_columns:
                if col in row and pd.notna(row[col]):
                    seq_preview = str(row[col])[:30] + "..." if len(str(row[col])) > 30 else str(row[col])
                    print(f"    {col}: '{seq_preview}'")
            
            # Zeige verfügbare Struktur-Daten
            for col in struct_columns:
                if col in row and pd.notna(row[col]):
                    struct_preview = str(row[col])[:30] + "..." if len(str(row[col])) > 30 else str(row[col])
                    print(f"    {col}: '{struct_preview}'")
        
        # Warnungen
        if validation['valid_sequences'] == 0:
            print(f"\n[WARNUNG] Keine gültigen Sequenzen gefunden!")
            print(f"Mögliche Ursachen:")
            print(f"  - ALL_best_candidates.txt wurde ohne struct-Daten erstellt")
            print(f"  - Sequenzen sind in anderen Spalten gespeichert")
            print(f"  - Original struct-Dateien enthalten keine Sequenzen")
        
        if validation['valid_structures'] == 0:
            print(f"\n[WARNUNG] Keine gültigen Strukturen gefunden!")
        
        return validation
    
    def create_summary_report(self, grouped_data: Dict[str, Dict[str, pd.DataFrame]], validation: Dict) -> Path:
        #Erstellt Zusammenfassungsbericht
        
        report_file = self.candidates_dir / "extraction_summary.txt"
        
        total_candidates = sum(sum(len(struct_df) for struct_df in anticodon_data.values()) 
                             for anticodon_data in grouped_data.values())
        
        with open(report_file, 'w') as f:
            f.write("tRNA Best Candidates FASTA Export Summary\n")
            f.write("==========================================\n\n")
            f.write(f"Export Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Source: ALL_best_candidates.txt\n")
            f.write(f"Format: Extended FASTA with #FS tags for mlocarna\n\n")
            
            f.write("Data Validation:\n")
            f.write(f"  Total candidates processed: {validation['total_rows']}\n")
            f.write(f"  Sequences available: {validation['has_sequence']} ({validation['sequence_sources']})\n")
            f.write(f"  Structures available: {validation['has_structure']} ({validation['structure_sources']})\n\n")
            
            f.write(f"Export Results:\n")
            f.write(f"  Total candidates exported: {total_candidates}\n")
            f.write(f"  Unique anticodons: {len(grouped_data)}\n")
            f.write(f"  Total FASTA files: {sum(len(struct_data) for struct_data in grouped_data.values())}\n\n")
            
            f.write("Detailed Breakdown:\n")
            f.write("-" * 60 + "\n")
            
            for anticodon, structure_data in sorted(grouped_data.items()):
                anticodon_total = sum(len(df) for df in structure_data.values())
                f.write(f"\nAnticodon {anticodon}: {anticodon_total} candidates\n")
                
                for structure_type, df in sorted(structure_data.items()):
                    count = len(df)
                    superphyla = df['Superphylum'].unique() if 'Superphylum' in df.columns else ['Unknown']
                    
                    # Score-Statistiken
                    score_col = 'IsotypeScore' if 'IsotypeScore' in df.columns else 'Score'
                    if score_col in df.columns:
                        avg_score = df[score_col].mean()
                        score_info = f", avg_score={avg_score:.1f}"
                    else:
                        score_info = ""
                    
                    safe_anticodon = re.sub(r'[^A-Za-z0-9]', '_', anticodon)
                    safe_structure = re.sub(r'[^A-Za-z0-9]', '_', structure_type)
                    filename = f"tRNA_{safe_anticodon}_{safe_structure}_candidates.fasta"
                    
                    f.write(f"  {structure_type}: {count} candidates{score_info}\n")
                    f.write(f"    Superphyla: {', '.join(superphyla[:5])}\n")  # Limit für Übersichtlichkeit
                    f.write(f"    File: {filename}\n")
        
        return report_file
    
    def export_all_candidates(self, df: pd.DataFrame) -> Dict[str, Dict[str, Path]]:
        #Hauptfunktion: Exportiert alle besten Kandidaten aus ALL_best_candidates.txt
        
        print(f"[INFO] Starte FASTA-Export für {len(df)} beste Kandidaten...")
        
        # Validiere Daten
        validation = self.validate_data_completeness(df)
        
        if validation['valid_sequences'] == 0:
            print(f"\n[ALTERNATIVE] Versuche Sequenz-Extraktion aus struct-Dateien...")
            df = self.try_extract_from_struct_files(df)
            # Validiere erneut
            validation = self.validate_data_completeness(df)
        
        if validation['valid_sequences'] == 0:
            print("[WARNUNG] Immer noch keine gültigen Sequenzen verfügbar!")
            print("Es werden Platzhalter-Daten verwendet - dies ist nicht optimal für mlocarna!")
            response = input("Möchten Sie trotzdem fortfahren? (y/n): ")
            if response.lower() != 'y':
                print("Abbruch. Bitte überprüfen Sie die Datenquellen.")
                return {}
        
        # Gruppiere nach Anticodon und Struktur-Typ
        grouped_data = self.group_by_anticodon_and_structure(df)
        
        print(f"[INFO] Gruppierung: {len(grouped_data)} Anticodons")
        for anticodon, structure_data in grouped_data.items():
            structure_counts = {st: len(sdf) for st, sdf in structure_data.items()}
            print(f"  {anticodon}: {structure_counts}")
        
        # Exportiere FASTA-Dateien
        exported_files = {}
        
        for anticodon, structure_data in grouped_data.items():
            exported_files[anticodon] = {}
            
            for structure_type, struct_df in structure_data.items():
                print(f"[EXPORT] {anticodon} + {structure_type}: {len(struct_df)} Kandidaten")
                
                try:
                    fasta_file = self.export_anticodon_structure_fasta(anticodon, structure_type, struct_df)
                    exported_files[anticodon][structure_type] = fasta_file
                except Exception as e:
                    print(f"  [FEHLER] Export für {anticodon}/{structure_type} fehlgeschlagen: {e}")
        
        # Erstelle Zusammenfassungsbericht
        summary_file = self.create_summary_report(grouped_data, validation)
        print(f"[INFO] Zusammenfassungsbericht erstellt: {summary_file}")
        
        return exported_files


def update_config_file():
    #Aktualisiert die config.py Datei um den CANDIDATES_DIR
    
    config_file = Path(__file__).resolve().parent.parent / "config.py"
    
    if config_file.exists():
        with open(config_file, 'r') as f:
            content = f.read()
        
        if 'CANDIDATES_DIR' not in content:
            addition = f"\n# Kandidaten-Verzeichnis für beste tRNA-Kandidaten\nCANDIDATES_DIR = TRNASCAN_OUTPUT_DIR / \"candidates\"\n"
            
            with open(config_file, 'a') as f:
                f.write(addition)
            
            print(f"[INFO] Config-Datei aktualisiert: CANDIDATES_DIR hinzugefügt")
    else:
        print(f"[WARNUNG] Config-Datei nicht gefunden: {config_file}")


def validate_best_candidates_file(file_path: Path) -> bool:
    #Validiert die ALL_best_candidates.txt Eingabedatei
    
    if not file_path.exists():
        print(f"[FEHLER] ALL_best_candidates.txt nicht gefunden: {file_path}")
        print("Führen Sie zuerst die Best-Candidates-Filterung aus!")
        return False
    
    try:
        df = pd.read_csv(file_path, sep="\t", nrows=1)
        required_columns = ['Sequence', 'Anticodon', 'Begin', 'End']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"[FEHLER] Fehlende Spalten in ALL_best_candidates.txt: {missing_columns}")
            return False
        
        print(f"[INFO] ALL_best_candidates.txt erfolgreich validiert")
        return True
        
    except Exception as e:
        print(f"[FEHLER] Beim Validieren von ALL_best_candidates.txt: {e}")
        return False


def main():
    """
    Hauptfunktion: Exportiert beste tRNA-Kandidaten aus ALL_best_candidates.txt im mlocarna-Format
    
    WORKFLOW:
    1. Lädt ALL_best_candidates.txt (enthält bereits Sequenzen und Strukturen)
    2. Extrahiert Sequenzen aus 'Seq'-Spalte und Strukturen aus 'Str'/'StrParse'-Spalte
    3. Gruppiert nach Anticodon UND Struktur-Typ
    4. Exportiert mlocarna-kompatible FASTA-Dateien mit #FS Tags
    """
    print("=== tRNA FASTA Export für mlocarna ===\n")
    
    # EINGABEDATEI
    input_file = TRNASCAN_OUTPUT_DIR / "ALL_best_candidates.txt"
    
    print(f"[INFO] Datenquelle: {input_file}")
    print(f"[INFO] Sequenzen und Strukturen werden direkt aus dieser Datei extrahiert")
    
    if not validate_best_candidates_file(input_file):
        return
    
    # Aktualisiere Config
    update_config_file()
    
    # Lade gefilterte beste Kandidaten
    print(f"\n[VERARBEITUNG] Lade beste Kandidaten aus ALL_best_candidates.txt...")
    df = pd.read_csv(input_file, sep="\t")
    print(f"[INFO] {len(df)} gefilterte beste Kandidaten geladen")
    
    # Initialisiere Exporter
    exporter = tRNAFastaExporter(TRNASCAN_OUTPUT_DIR)
    
    print(f"\n[EXPORT] Beginne FASTA-Export...")
    
    # Hauptprozess: Export
    exported_files = exporter.export_all_candidates(df)
    
    # Zusammenfassung
    print(f"\n=== EXPORT ABGESCHLOSSEN ===")
    total_files = sum(len(structure_files) for structure_files in exported_files.values())
    print(f"Exportierte FASTA-Dateien: {total_files}")
    print(f"Ausgabe-Verzeichnis: {exporter.candidates_dir}")
    print(f"Format: Extended FASTA mit #FS Tags für mlocarna")
    print(f"")
    
    for anticodon, structure_files in exported_files.items():
        print(f"Anticodon {anticodon}:")
        for structure_type, fasta_file in structure_files.items():
            file_size = fasta_file.stat().st_size if fasta_file.exists() else 0
            print(f"  {structure_type}: {fasta_file.name} ({file_size} bytes)")
    
    print(f"\nFASTA-Dateien sind bereit für mlocarna-Analyse!")
    print(f"Sequenzen und Strukturen aus ALL_best_candidates.txt extrahiert")
    print(f"Gruppiert nach Anticodon UND Struktur-Typ")
    print(f"Fixed Structure Tags (#FS) für mlocarna hinzugefügt")


if __name__ == "__main__":
    main()


    #"#FS" wurde zu " #FS" geändert, um Kompatibilität mit mlocarna zu gewährleisten