#!/usr/bin/env python3

#Inspected mal die Struktur der sys.path.append(str(Path(__file__).resolve().parent.parent)) Datei, bevor wieder Quatsch passiert 

import sys
from pathlib import Path
import pandas as pd

# Import config
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import *

def inspect_file():
    #Zeigt detaillierte Informationen über die Dateistruktur
    filepath = TRNASCAN_OUTPUT_DIR / "ALL_noncanonical_mainstruct.txt"
    
    print(f"=== Inspektion von {filepath} ===\n")
    
    if not filepath.exists():
        print(f"FEHLER: Datei nicht gefunden!")
        return
    
    # Lese erste Zeilen
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    print(f"Datei enthält {len(lines)} Zeilen (inkl. Header)\n")
    
    # Zeige Header
    header = lines[0].strip().split('\t')
    print("HEADER:")
    print("-" * 80)
    for i, col in enumerate(header):
        print(f"  Spalte {i:2d}: {col}")
    print(f"\nGesamt: {len(header)} Spalten\n")
    
    # Zeige erste Datenzeilen
    print("ERSTE 3 DATENZEILEN:")
    print("-" * 80)
    for i in range(1, min(4, len(lines))):
        parts = lines[i].strip().split('\t')
        print(f"\nZeile {i}:")
        print(f"  Sequence ID: {parts[0] if len(parts) > 0 else 'N/A'}")
        print(f"  Type-Anticodon: {parts[4] if len(parts) > 4 else 'N/A'}-{parts[5] if len(parts) > 5 else 'N/A'}")
        print(f"  Position: {parts[2] if len(parts) > 2 else 'N/A'}-{parts[3] if len(parts) > 3 else 'N/A'}")
        print(f"  Superphylum: {parts[14] if len(parts) > 14 else 'N/A'}")
        print(f"  Stem Count: {parts[23] if len(parts) > 23 else 'N/A'}")
        print(f"  Structure: {parts[20][:50] if len(parts) > 20 else 'N/A'}...")
    
    # Lade als DataFrame für Statistiken
    try:
        df = pd.read_csv(filepath, sep='\t')
        
        print("\n\nSTATISTIKEN:")
        print("-" * 80)
        
        # Verteilung nach Superphylum
        print("\nVerteilung nach Superphylum:")
        for sp, count in df['Superphylum'].value_counts().items():
            print(f"  {sp}: {count}")
        
        # Verteilung nach Stem Count
        print("\nVerteilung nach Stem Count:")
        for sc, count in df['StemCount'].value_counts().sort_index().items():
            print(f"  {sc} Stems: {count}")
        
        # Top 10 Anticodons
        print("\nTop 10 Anticodons:")
        for ac, count in df['Anticodon'].value_counts().head(10).items():
            print(f"  {ac}: {count}")
        
        # Isotype Score Statistiken
        print("\nIsotype Score Statistiken:")
        print(f"  Min: {df['IsotypeScore'].min():.2f}")
        print(f"  Max: {df['IsotypeScore'].max():.2f}")
        print(f"  Mittelwert: {df['IsotypeScore'].mean():.2f}")
        print(f"  Median: {df['IsotypeScore'].median():.2f}")
        
        # Intron Statistiken
        intron_df = df[df['IntronCount'] > 0]
        print(f"\nIntron Statistiken:")
        print(f"  Kandidaten mit Introns: {len(intron_df)} ({len(intron_df)/len(df)*100:.1f}%)")
        if len(intron_df) > 0:
            print(f"  Durchschnittliche Intron-Anzahl: {intron_df['IntronCount'].mean():.2f}")
            print(f"  Max Intron-Anzahl: {intron_df['IntronCount'].max()}")
        
    except Exception as e:
        print(f"\nFehler beim Laden als DataFrame: {e}")
    
    print("\n" + "="*80)
    print("Inspektion doooone")

if __name__ == "__main__":
    inspect_file()