#!/usr/bin/env python3
#Berechnet Ensemble Defect (ED) Scores für nicht-kanonische tRNA-Kandidaten
#Verwendet RNA.fold_compound Algorithmus

import RNA
import pandas as pd
from pathlib import Path
import sys

# Import config
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import *


def parse_noncanonical_file(filepath):
    #Parst die Datei ALL_noncanonical_mainstruct.txt
    # Lade Datei
    df = pd.read_csv(filepath, sep='\t', dtype=str)
    
    # Numerische Spalten konvertieren
    int_cols = ['Begin', 'End', 'IntronCount', 'TotalIntronLen', 'MatureLen', 'StemCount']
    float_cols = ['CM', 'HMMScore', 'Str2Score', 'IsotypeScore']
    
    for col in int_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
    
    for col in float_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Erstelle dict für jeden Kandidaten
    candidates = df.to_dict(orient='records')
    
    print(f"Erfolgreich {len(candidates)} nicht-kanonische tRNA-Kandidaten geladen")
    return candidates


def calculate_ed_score(seq, reference_struct):
    #Berechnet den Ensemble Defect Score 
    #Verwendet ausschließlich den gegebenen Algorithmus
    
    if not seq or not reference_struct:
        return None
    
    # Konvertiere DNA zu RNA (T -> U)
    seq = seq.upper().replace('T', 'U')
    
    # Überprüfe ob Sequenz und Struktur gleich lang sind
    if len(seq) != len(reference_struct):
        print(f"Warnung: Sequenz ({len(seq)}) und Struktur ({len(reference_struct)}) haben unterschiedliche Längen")
        return None
    
    try:
        # create a fold compound for the sequence to predict various thermodynamic properties
        fc = RNA.fold_compound(seq)
        
        # compute MFE and corresponding structure
        ss, mfe = fc.mfe()
        
        # rescale Boltzmann factors to ensure numeric stability while computing partition function
        fc.exp_params_rescale(mfe)
        
        # compute partition function
        pp, dG = fc.pf()
        
        # compute ensemble defect with respect to reference structure
        ed = fc.ensemble_defect(reference_struct)
        
        return ed
    except Exception as e:
        print(f"Fehler bei ED-Berechnung: {e}")
        return None


def main():
    
    # Eingabedatei
    noncanonical_file = TRNASCAN_OUTPUT_DIR / "ALL_noncanonical_mainstruct.txt"
    
    if not noncanonical_file.exists():
        print(f"Fehler: {noncanonical_file} nicht gefunden!")
        sys.exit(1)
    
    # Parse nicht-kanonische Kandidaten
    print("Lese nicht-kanonische tRNA-Kandidaten...")
    candidates = parse_noncanonical_file(noncanonical_file)
    print(f"Gefunden: {len(candidates)} Kandidaten")
    
    # Berechne ED Scores
    results = []
    good_candidates = []
    errors = []
    
    print("\nBerechne ED Scores...")
    for i, candidate in enumerate(candidates):
        if i % 100 == 0:
            print(f"  Verarbeitet: {i}/{len(candidates)}")
        
        # Hole Sequenz und Referenzstruktur
        seq = candidate.get('Seq')
        reference_struct = candidate.get('StrParse')
        
        if not seq or not reference_struct:
            errors.append({
                'sequence_id': candidate.get('Sequence'),
                'error': 'Fehlende Sequenz oder Struktur'
            })
            continue
        
        # Berechne ED Score
        ed_score = calculate_ed_score(seq, reference_struct)
        
        if ed_score is not None:
            # Füge ED Score zum Kandidaten hinzu
            candidate['ed_score'] = ed_score
            candidate['sequence_rna'] = seq.upper().replace('T', 'U')
            results.append(candidate)
            
            # Gute Kandidaten (ED ≤ 0.2)
            if ed_score <= 0.2:
                good_candidates.append(candidate)
                print(f"  Guter Kandidat gefunden: {candidate.get('Sequence')} "
                      f"({candidate.get('Type')}-{candidate.get('Anticodon')}) "
                      f"ED={ed_score:.3f}")
        else:
            errors.append({
                'sequence_id': candidate.get('Sequence'),
                'seq_length': len(seq) if seq else 0,
                'struct_length': len(reference_struct) if reference_struct else 0,
                'error': 'ED-Berechnung fehlgeschlagen'
            })
    
    # Speichere Ergebnisse
    output_dir = OUTPUT_DIR / "ed_scores"
    output_dir.mkdir(exist_ok=True)
    
    # Alle Ergebnisse
    output_file = output_dir / "ALL_noncanonical_ED_scores.tsv"
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, sep='\t', index=False)
    print(f"\nErgebnisse gespeichert in: {output_file}")
    
    # Gute Kandidaten
    if good_candidates:
        good_output = output_dir / "noncanonical_good_candidates_ED.tsv"
        good_df = pd.DataFrame(good_candidates)
        good_df.to_csv(good_output, sep='\t', index=False)
        print(f"Gute Kandidaten (ED ≤ 0.2) gespeichert in: {good_output}")
        
        # FASTA für gute Kandidaten
        fasta_output = output_dir / "good_candidates_ed.fasta"
        with open(fasta_output, 'w') as f:
            for candidate in good_candidates:
                header = (f">{candidate.get('Sequence')}_{candidate.get('tRNA#_x')}_"
                         f"{candidate.get('Type')}-{candidate.get('Anticodon')}_"
                         f"ED={candidate['ed_score']:.3f}")
                f.write(f"{header}\n{candidate['sequence_rna']}\n")
        print(f"FASTA-Datei gespeichert in: {fasta_output}")
    
    # Fehler protokollieren
    if errors:
        error_output = output_dir / "ed_calculation_errors.tsv"
        error_df = pd.DataFrame(errors)
        error_df.to_csv(error_output, sep='\t', index=False)
        print(f"\nFehler protokolliert in: {error_output}")
    
    # Zusammenfassung
    print(f"\n=== Zusammenfassung ===")
    print(f"Gesamt verarbeitet: {len(candidates)}")
    print(f"Erfolgreich berechnet: {len(results)}")
    print(f"Fehler: {len(errors)}")
    print(f"Gute Kandidaten (ED ≤ 0.2): {len(good_candidates)}")
    
    if results:
        ed_scores = [r['ed_score'] for r in results]
        print(f"\nED Score Statistiken:")
        print(f"  Minimum: {min(ed_scores):.3f}")
        print(f"  Maximum: {max(ed_scores):.3f}")
        print(f"  Durchschnitt: {sum(ed_scores)/len(ed_scores):.3f}")
        print(f"  Median: {sorted(ed_scores)[len(ed_scores)//2]:.3f}")
        
        # Verteilung
        print("\nED Score Verteilung:")
        ranges = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 1.0), (1.0, float('inf'))]
        for low, high in ranges:
            count = sum(1 for ed in ed_scores if low <= ed < high)
            if count > 0:
                print(f"  {low:.1f} - {high:.1f}: {count} ({count/len(ed_scores)*100:.1f}%)")


if __name__ == "__main__":
    main()