#!/usr/bin/env python3
"""
Zentrale Konfigurationsdatei für das tRNA-Analyse Projekt
"""
from pathlib import Path
import sys


# Basis-Verzeichnis (wo config.py liegt)
BASE_DIR = Path(__file__).resolve().parent

# Datenverzeichnisse
DATA_DIR = BASE_DIR / "data"
GENOMES_DIR = DATA_DIR / "genomes"
CCA_DIR = DATA_DIR / "CCA_enzymes"

# Output-Verzeichnisse
OUTPUT_DIR = BASE_DIR / "output"
TRNASCAN_OUTPUT_DIR = OUTPUT_DIR / "trnascan_output" #output/trnascan_output
PLOTS_DIR = OUTPUT_DIR / "plots"
ED_DIR = OUTPUT_DIR / "ed_scores"
CANDIDATES_DIR = TRNASCAN_OUTPUT_DIR / "candidates"
TREE_DIR = PLOTS_DIR / "trees"

# mlocarna Verzeichnisse
MLOCARNA_OUTPUT_DIR = OUTPUT_DIR / "mlocarna_results"
MLOCARNA_ALIGNMENTS_DIR = MLOCARNA_OUTPUT_DIR / "alignments"
MLOCARNA_CONSENSUS_DIR = MLOCARNA_OUTPUT_DIR / "consensus_structures"
MLOCARNA_LOGS_DIR = MLOCARNA_OUTPUT_DIR / "logs"

# Skript-Verzeichnis
SCRIPTS_DIR = BASE_DIR / "scripts"

# Stelle sicher, dass alle Verzeichnisse existieren
for dir_path in [GENOMES_DIR, TRNASCAN_OUTPUT_DIR, PLOTS_DIR, CANDIDATES_DIR,
                 MLOCARNA_OUTPUT_DIR, MLOCARNA_ALIGNMENTS_DIR, 
                 MLOCARNA_CONSENSUS_DIR, MLOCARNA_LOGS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)
 

#verzeichnisse vorhanden? 
TEST_DIR = BASE_DIR / "test_output"
TEST_DIR.mkdir(parents=True, exist_ok=True)
sys.path.append(str(Path(__file__).resolve().parent.parent))
print("Verzeichnis existiert parently:", TEST_DIR.exists())