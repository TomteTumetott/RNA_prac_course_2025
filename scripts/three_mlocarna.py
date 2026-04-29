#!/usr/bin/env python3
#Führt mlocarna auf allen tRNA-Kandidaten-FASTA-Dateien aus
#und speichert die Alignments im Stockholm-Format sowie die Konsensus-Strukturen.


import subprocess
import sys
from pathlib import Path
import re
from datetime import datetime

# Füge das Hauptverzeichnis zum Python-Pfad hinzu
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import (CANDIDATES_DIR, MLOCARNA_OUTPUT_DIR, 
                   MLOCARNA_ALIGNMENTS_DIR, MLOCARNA_CONSENSUS_DIR, 
                   MLOCARNA_LOGS_DIR)

# Konfiguration
DEBUG = True
DRY_RUN = False  # Wenn True, werden Befehle nur angezeigt, nicht ausgeführt

def find_fasta_files(directory: Path) -> list[Path]:
    #Findet alle FASTA-Dateien im angegebenen Verzeichnis.
    
    fasta_files = list(directory.glob("*.fasta")) + list(directory.glob("*.fa"))
    return sorted(fasta_files)

def extract_consensus_structure(log_content: str) -> str:
    #Extrahiert die Konsensus-Struktur aus der mlocarna/alifold Ausgabe.
    #Sucht nach dem alifold Output im Log.
    
    consensus_lines = []
    in_alifold_section = False
    
    for line in log_content.split('\n'):
        # Suche nach alifold Ausgabe
        if 'alifold' in line.lower() or 'consensus' in line.lower():
            in_alifold_section = True
        
        # Sammle Konsensus-Struktur (Zeilen mit Klammern und Punkten)
        if in_alifold_section and re.match(r'^[.\(\)]+\s*\(.*\)\s*$', line.strip()):
            consensus_lines.append(line.strip())
        
        # Oder suche nach typischen RNA-Struktur-Patterns
        if re.match(r'^[ACGTU]+\s+[.\(\)]+\s*[-\d.]+', line):
            consensus_lines.append(line.strip())
    
    return '\n'.join(consensus_lines)

def run_mlocarna(fasta_file: Path) -> dict:
    #Führt mlocarna auf einer FASTA-Datei aus.
    #Returns:
    #    dict mit Statusinfos und Pfaden zu Output-Dateien
    
    # Basis-Namen für Output-Dateien
    base_name = fasta_file.stem
    
    # Output-Pfade
    output_dir = MLOCARNA_OUTPUT_DIR / base_name
    output_dir.mkdir(exist_ok=True)
    
    alignment_file = MLOCARNA_ALIGNMENTS_DIR / f"{base_name}.stockholm"
    consensus_file = MLOCARNA_CONSENSUS_DIR / f"{base_name}_consensus.txt"
    log_file = MLOCARNA_LOGS_DIR / f"{base_name}.log"
    
    # mlocarna Befehl konstruieren
    cmd = [
        "mlocarna",
        "--stockholm",  # Stockholm Format Output
        str(fasta_file),
        "--tgtdir", str(output_dir)  # Zielverzeichnis für alle Outputs
    ]
    
    if DEBUG:
        print(f"\nVerarbeite: {fasta_file.name}")
        print(f"Befehl: {' '.join(cmd)}")
    
    if DRY_RUN:
        return {
            "status": "dry_run",
            "fasta": str(fasta_file),
            "alignment": str(alignment_file),
            "consensus": str(consensus_file),
            "log": str(log_file)
        }
    
    try:
        # Führe mlocarna aus und erfasse stdout/stderr
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Speichere das komplette Log
        with open(log_file, 'w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Started: {datetime.now()}\n")
            f.write("\n=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)
        
        # Extrahiere Konsensus-Struktur aus der Ausgabe
        consensus_structure = extract_consensus_structure(result.stdout + result.stderr)
        
        if consensus_structure:
            with open(consensus_file, 'w') as f:
                f.write(f"# Consensus structure for {base_name}\n")
                f.write(f"# Generated: {datetime.now()}\n")
                f.write(f"# Source: {fasta_file.name}\n\n")
                f.write(consensus_structure)
        
        # Kopiere Stockholm-Alignment
        stockholm_source = output_dir / "results" / "result.aln"
        if stockholm_source.exists():
            with open(stockholm_source, 'r') as src:
                content = src.read()
            with open(alignment_file, 'w') as dst:
                dst.write(content)
        
        return {
            "status": "success",
            "fasta": str(fasta_file),
            "alignment": str(alignment_file),
            "consensus": str(consensus_file),
            "log": str(log_file),
            "output_dir": str(output_dir)
        }
        
    except subprocess.CalledProcessError as e:
        error_msg = f"mlocarna fehlgeschlagen: {e}"
        if DEBUG:
            print(f"  FEHLER: {error_msg}")
            print(f"  STDOUT: {e.stdout}")
            print(f"  STDERR: {e.stderr}")
        
        # Speichere Fehler-Log
        with open(log_file, 'w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Started: {datetime.now()}\n")
            f.write(f"ERROR: {error_msg}\n")
            f.write("\n=== STDOUT ===\n")
            f.write(e.stdout if e.stdout else "")
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr if e.stderr else "")
        
        return {
            "status": "error",
            "fasta": str(fasta_file),
            "error": error_msg,
            "log": str(log_file)
        }

def create_summary_report(results: list[dict]) -> Path:
    #Erstellt einen zusammenfassenden Bericht über alle mlocarna-Läufe.
    
    report_file = MLOCARNA_OUTPUT_DIR / f"mlocarna_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    successful = [r for r in results if r["status"] == "success"]
    failed = [r for r in results if r["status"] == "error"]
    
    with open(report_file, 'w') as f:
        f.write(f"mlocarna Zusammenfassung\n")
        f.write(f"========================\n")
        f.write(f"Generiert: {datetime.now()}\n\n")
        
        f.write(f"Gesamt: {len(results)} FASTA-Dateien verarbeitet\n")
        f.write(f"Erfolgreich: {len(successful)}\n")
        f.write(f"Fehlgeschlagen: {len(failed)}\n\n")
        
        if successful:
            f.write("Erfolgreich verarbeitete Dateien:\n")
            f.write("-" * 50 + "\n")
            for r in successful:
                f.write(f"\n{Path(r['fasta']).name}\n")
                f.write(f"  Alignment: {r['alignment']}\n")
                f.write(f"  Konsensus: {r['consensus']}\n")
                f.write(f"  Output-Dir: {r['output_dir']}\n")
        
        if failed:
            f.write("\n\nFehlgeschlagene Dateien:\n")
            f.write("-" * 50 + "\n")
            for r in failed:
                f.write(f"\n{Path(r['fasta']).name}\n")
                f.write(f"  Fehler: {r['error']}\n")
                f.write(f"  Log: {r['log']}\n")
    
    return report_file

def main():
    #Hauptfunktion: Führt mlocarna auf allen FASTA-Dateien aus.
    
    print(f"mlocarna Runner")
    print(f"===============")
    print(f"Suche FASTA-Dateien in: {CANDIDATES_DIR}")
    
    # Finde alle FASTA-Dateien
    fasta_files = find_fasta_files(CANDIDATES_DIR)
    
    if not fasta_files:
        print("Keine FASTA-Dateien gefunden!")
        return
    
    print(f"Gefunden: {len(fasta_files)} FASTA-Dateien")
    
    ## Bestätigung bei vielen Dateien
    #if len(fasta_files) > 10 and not DRY_RUN:
    #    response = input(f"\nMöchten Sie wirklich {len(fasta_files)} mlocarna-Jobs starten? (y/n): ")
    #    if response.lower() != 'yy':
    #        print("Abgebrochen.")
    #        return
    
    # Verarbeite alle FASTA-Dateien
    results = []
    for i, fasta_file in enumerate(fasta_files, 1):
        print(f"\n[{i}/{len(fasta_files)}] {fasta_file.name}")
        result = run_mlocarna(fasta_file)
        results.append(result)
        
        if result["status"] == "success":
            print(f"Erfolgreich")
        else:
            print(f"Fehler: {result.get('error', 'Unbekannt')}")
    
    # Erstelle Zusammenfassung
    report_file = create_summary_report(results)
    print(f"\n\nZusammenfassung gespeichert in: {report_file}")
    
    # Statistiken
    successful = sum(1 for r in results if r["status"] == "success")
    print(f"\nErgebnis: {successful}/{len(results)} erfolgreich verarbeitet")

if __name__ == "__main__":
    main()