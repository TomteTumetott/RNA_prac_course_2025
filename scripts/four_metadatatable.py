#!/usr/bin/env python3
#Erstellt eine Metadaten-Tabelle für Isoleucine tRNAs mit korrigierter CCA-Enzym Zuordnung


import pandas as pd
from pathlib import Path
import sys

# Zentrale Konfiguration importieren
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, CCA_DIR, OUTPUT_DIR

def load_cca_enzyme_data():
    #Lädt alle CCA-Enzym Tabellen und kombiniert sie
    cca_data = []
    
    print("Lade CCA-Enzym Daten...")
    
    # Asgard
    asgard_file = CCA_DIR / "CCA_asgard_presence_absence.csv"
    if asgard_file.exists():
        asgard_df = pd.read_csv(asgard_file)
        if 'Accession' in asgard_df.columns:
            asgard_df.rename(columns={'Accession': 'Header'}, inplace=True)
        cca_data.append(asgard_df)
        print(f"  Asgard: {len(asgard_df)} Einträge geladen")
    
    # DPANN
    dpann_file = CCA_DIR / "CCA_DPANN_presence_absence.csv"
    if dpann_file.exists():
        dpann_df = pd.read_csv(dpann_file)
        cca_data.append(dpann_df)
        print(f"  DPANN: {len(dpann_df)} Einträge geladen")
    
    # Euryarchaeota
    eury_file = CCA_DIR / "CCA_Euryarchaeota_presence_absence.csv"
    if eury_file.exists():
        eury_df = pd.read_csv(eury_file)
        cca_data.append(eury_df)
        print(f"  Euryarchaeota: {len(eury_df)} Einträge geladen")
    
    # TACK
    tack_file = CCA_DIR / "CCA_TACK_presence_absence.csv"
    if tack_file.exists():
        tack_df = pd.read_csv(tack_file)
        cca_data.append(tack_df)
        print(f"  TACK: {len(tack_df)} Einträge geladen")
    
    # Kombiniere alle Daten
    if cca_data:
        all_cca = pd.concat(cca_data, ignore_index=True)
        all_cca['Header'] = all_cca['Header'].astype(str).str.replace('.fasta', '', regex=False)
        print(f"\nGesamt CCA-Enzym Einträge: {len(all_cca)}")
        
        # Stelle sicher dass CCA1 und CCA2 numerisch sind
        all_cca['CCA1'] = pd.to_numeric(all_cca['CCA1'], errors='coerce')
        all_cca['CCA2'] = pd.to_numeric(all_cca['CCA2'], errors='coerce')
        
        return all_cca
    else:
        print("Keine CCA-Enzym Daten gefunden!")
        return pd.DataFrame()

def create_ile_metadata_table():
    #Erstellt die Metadaten-Tabelle für Isoleucine tRNAs mit korrigierter CCA-Zuordnung
    
    # 1. Lade die kombinierte mainstruct Datei
    mainstruct_file = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct.txt"
    if not mainstruct_file.exists():
        print(f"Fehler: {mainstruct_file} nicht gefunden!")
        return None
    
    print("\nLade tRNA-Daten...")
    df = pd.read_csv(mainstruct_file, sep='\t')
    print(f"Gesamt tRNAs in Datei: {len(df)}")
    
    # 2. Filtere nur Isoleucine tRNAs
    ile_df = df[df['Type'] == 'Ile'].copy()
    print(f"Gefundene Isoleucine tRNAs: {len(ile_df)}")
    
    if len(ile_df) == 0:
        print("Keine Isoleucine tRNAs gefunden!")
        return None
    
    # 3. Extrahiere Species_ID
    ile_df['Species_ID'] = ile_df['Name'].str.split('__').str[1]
    
    # Debug: Zeige einige Species_IDs
    print("\nBeispiel Species_IDs aus Isoleucine tRNAs:")
    for sid in ile_df['Species_ID'].unique()[:5]:
        print(f"  {sid}")
    
    # 4. Lade CCA-Enzym Daten
    cca_enzymes = load_cca_enzyme_data()
    
    if cca_enzymes.empty:
        print("Keine CCA-Daten geladen!")
        return None
    
    # 5. Merge mit CCA-Daten
    print("\nFühre Merge durch...")
    merged_df = pd.merge(
        ile_df,
        cca_enzymes[['Header', 'CCA1', 'CCA2']],
        left_on='Species_ID',
        right_on='Header',
        how='left'
    )
    
    # 6. Überprüfe Merge-Ergebnis
    matches = merged_df['CCA1'].notna().sum()
    print(f"\nErfolgreiche CCA-Matches: {matches}/{len(merged_df)}")
    
    if matches > 0:
        print("\nBeispiele erfolgreicher Matches:")
        matched_examples = merged_df[merged_df['CCA1'].notna()].head(3)
        for idx, row in matched_examples.iterrows():
            print(f"  {row['Species_ID']}: CCA1={row['CCA1']}, CCA2={row['CCA2']}")
    
    # 7. Erstelle CCA-enzyme Spalte mit verbesserter Logik
    def get_cca_enzyme_type(row):
        #Bestimmt den CCA-Enzym Typ basierend auf CCA1 und CCA2
        # Wenn beide Werte fehlen -> NA
        if pd.isna(row['CCA1']) and pd.isna(row['CCA2']):
            return 'NA'
        
        # Konvertiere zu int, behandle NaN als 0
        try:
            cca1 = int(row['CCA1']) if pd.notna(row['CCA1']) else 0
            cca2 = int(row['CCA2']) if pd.notna(row['CCA2']) else 0
        except:
            return 'NA'
        
        # Bestimme Enzym-Typ
        if cca1 == 1 and cca2 == 1:
            return 'CCA1+CCA2'
        elif cca1 == 1:
            return 'CCA1'
        elif cca2 == 1:
            return 'CCA2'
        else:
            return 'None'
    
    merged_df['CCA-enzyme'] = merged_df.apply(get_cca_enzyme_type, axis=1)
    
    # Debug: Zeige CCA-enzyme Verteilung
    print("\nCCA-enzyme Verteilung:")
    print(merged_df['CCA-enzyme'].value_counts())
    
    # 8. Erstelle TipName
    merged_df['TipName'] = (
        merged_df['Superphylum'] + '__' +
        merged_df['Species_ID'] +  '__' +
        merged_df['Type'] + '_' + merged_df['Anticodon'] + '__' +
        merged_df['Sequence'] + '.trna' + merged_df['tRNA#'].astype(str)
    )
             
    # 10. Erstelle finale Tabelle
    final_table = pd.DataFrame({
        'CCA-enzyme': merged_df['CCA-enzyme'],
        'superphyla': merged_df['Superphylum'],
        'CCA-encoded': merged_df['CCA_Check'],
        'TipName': merged_df['TipName'],
        'tRNA-anticodon': merged_df['Type'] + '_' + merged_df['Anticodon'],
        'Species_ID': merged_df['Species_ID']
    })
    
    # 11. Sortiere
    final_table = final_table.sort_values(['superphyla', 'Species_ID'])
    
    # 12. Speichere
    output_file = OUTPUT_DIR / "TABLE_metadata_Isoleucine.txt"
    final_table.to_csv(output_file, sep='\t', index=False)
    print(f"\nTabelle gespeichert als: {output_file}")
    
    # 13. Zusammenfassung
    print("\n=== ZUSAMMENFASSUNG ===")
    print(f"Gesamt Isoleucine tRNAs: {len(final_table)}")
    
    print("\nVerteilung nach Superphylum:")
    print(final_table['superphyla'].value_counts())
    
    print("\nCCA-enzyme Verteilung:")
    enzyme_counts = final_table['CCA-enzyme'].value_counts()
    for enzyme_type, count in enzyme_counts.items():
        percentage = (count / len(final_table)) * 100
        print(f"  {enzyme_type}: {count} ({percentage:.1f}%)")
    
    print("\nCCA-encoded Verteilung:")
    cca_encoded = final_table['CCA-encoded'].value_counts()
    for val, count in cca_encoded.items():
        print(f"  {'Mit CCA' if val == 1 else 'Ohne CCA'}: {count}")
    
    # 14. Detaillierte Analyse der gematchten Daten
    matched_data = final_table[final_table['CCA-enzyme'] != 'NA']
    if len(matched_data) > 0:
        print(f"\nSpecies mit CCA-Enzym Daten: {len(matched_data.groupby('Species_ID'))}")
        print("\nCCA-Enzyme nach Superphylum:")
        crosstab = pd.crosstab(matched_data['superphyla'], matched_data['CCA-enzyme'])
        print(crosstab)
    
    return final_table

# Hauptausführung
if __name__ == "__main__":
    print("Erstelle Isoleucine tRNA Metadaten-Tabelle")
    
    result = create_ile_metadata_table()
    
    if result is not None:
        print("\nTabelle erfolgreich erstellt!")
        
        # Überprüfe ob CCA-enzyme korrekt befüllt wurde
        na_count = (result['CCA-enzyme'] == 'NA').sum()
        total = len(result)
        matched = total - na_count
        
        print(f"\CCA-Enzyme Status:")
        print(f"  Erfolgreich zugeordnet: {matched}")
        print(f"  Nicht gefunden (NA): {na_count}")
        print(f"  Erfolgsrate: {(matched/total)*100:.1f}%")