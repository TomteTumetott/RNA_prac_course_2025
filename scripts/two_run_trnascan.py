#!/usr/bin/env python3
import subprocess
from pathlib import Path
import sys

# Füge das Hauptverzeichnis zum Python-Pfad hinzu
sys.path.append(str(Path(__file__).resolve().parent.parent))

from config import GENOMES_DIR, TRNASCAN_OUTPUT_DIR
#####################################################################################
# tRNA_Scan mit allen special einstellungen, 
# struct detail main output (txt)
#####################################################################################
for fasta_file in GENOMES_DIR.glob("*.fasta"):

    stem        = fasta_file.stem
    outfile     = TRNASCAN_OUTPUT_DIR / f"{stem}_main.txt"
    detailfile  = TRNASCAN_OUTPUT_DIR / f"{stem}_detail.txt"
    structfile  = TRNASCAN_OUTPUT_DIR / f"{stem}_struct.txt"
    fastafile   = TRNASCAN_OUTPUT_DIR / f"{stem}_trnas.fa"

    subprocess.run([
        "tRNAscan-SE",
        "-A",                              # archaeale Modelle
        "-H",                              # Score-Aufschlüsselung
        "--pseudogene",                    # Pseudogen-Check
        "--isotype",   
        "--detail",                        # Isotyp-spezifische Scores
        "-o", str(outfile),                # Haupttabelle
        "-f", str(structfile),             # Sekundärstrukturen
        "-m", str(detailfile),             # Detail/Stats
        "-a", str(fastafile),              # FASTA der tRNAs
        str(fasta_file)                    # Eingabe-FASTA
    ], check=True)

print("Alle passenden Dateien verarbeitet.")