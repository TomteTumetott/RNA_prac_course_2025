#!/usr/bin/env python3
#Re-fold recovered tRNA candidates with RNAfold using constraints


import pandas as pd
import RNA
from pathlib import Path
import sys

# Konfigurationspfad einbinden
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import ED_DIR, OUTPUT_DIR

# Eingabedatei (gute Kandidaten mit ED < 0.2)
input_file = ED_DIR / "noncanonical_good_candidates_ED.tsv"

if not input_file.exists():
    print(f"[ERROR] Eingabedatei nicht gefunden: {input_file}")
    sys.exit(1)

# Datei laden
df = pd.read_csv(input_file, sep="\t")
if df.empty:
    print("[ERROR] Eingabedatei ist leer.")
    sys.exit(1)

def find_constraint_positions(sequence):
    #Constraints: Position 9 (A/G) + erste zwei Us im D-Arm (Position 10–24)
    constraints = []
    if len(sequence) > 8 and sequence[8] in ['A', 'G']:
        constraints.append(8)
    u_count = 0
    for i in range(10, min(25, len(sequence))):
        if sequence[i] == 'U' and u_count < 2:
            constraints.append(i)
            u_count += 1
        if u_count == 2:
            break
    return constraints

def refold_with_constraints(seq_rna):
    # Fold-Compound ganz normal initialisieren
    fc = RNA.fold_compound(seq_rna)

    # Initiales MFE zum Reskalieren
    ss, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()

    # Constraint-String vorbereiten
    constraint_string = ['.'] * len(seq_rna)
    for idx in find_constraint_positions(seq_rna):
        constraint_string[idx] = 'x'
    constraint_str = ''.join(constraint_string)

    # Constraints setzen (als Direktiven-String)
    fc.constraints_add(constraint_str, RNA.CONSTRAINT_DB_DEFAULT)

    # Neue Struktur mit Constraints
    new_structure, new_mfe = fc.mfe()

    # Zähle grob die Stems anhand Klammer-Blöcke
    simplified = []
    last = ''
    for char in new_structure:
        if char in '().' and char != last:
            simplified.append(char)
        last = char
    stem_count = simplified.count('(')
    is_canonical = stem_count in [4, 5]

    return new_structure, new_mfe, stem_count, is_canonical



# Refolding wird jetzt durchgeführt nach der Vorrede
results = []
for idx, row in df.iterrows():
    seq = row['sequence_rna']
    struct, mfe, stems, is_canon = refold_with_constraints(seq)
    results.append({
        'Sequence': row['Sequence'],
        'Type': row['Type'],
        'Anticodon': row['Anticodon'],
        'RefoldedStructure': struct,
        'RefoldedMFE': mfe,
        'StemCount': stems,
        'Canonical': is_canon
    })

# safe that
refolded_df = pd.DataFrame(results)

output_all = ED_DIR / "noncanonical_refolded_all.tsv"
refolded_df.to_csv(output_all, sep="\t", index=False)
print(f"[INFO] Alle Refolding-Ergebnisse gespeichert: {output_all}")

# output? file -> doch kanonische dabei ? 
canonical_df = refolded_df[refolded_df["Canonical"] == True]
output_canon = ED_DIR / "noncanonical_refolded_canonical.tsv"
canonical_df.to_csv(output_canon, sep="\t", index=False)
print(f"[INFO] Kanonische Strukturen gespeichert: {output_canon}")
