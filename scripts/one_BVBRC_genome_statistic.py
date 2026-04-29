#!/usr/bin/env python3
"""
Auswertung BV-BRC-Genome:
1. Qualitätsdaten laden
2. fehlende Genome-Größe & Genanzahl notfalls aus FASTA errechnen
3. tRNA-Tabelle anhängen
4. Plots erzeugen
"""

import os, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys

# Add parent directory to path to find config.py
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

# Now import config
import config

# Use config variables
DATA_DIR = config.DATA_DIR
GENOMES_DIR = config.GENOMES_DIR
PLOTS_DIR = config.PLOTS_DIR
OUTPUT_DIR = config.OUTPUT_DIR
BASE_DIR = config.BASE_DIR

print(f"[INFO] Arbeitsverzeichnis: {Path.cwd()}")
print(f"[INFO] Script directory: {Path(__file__).parent}")
print(f"[INFO] Parent directory: {parent_dir}")
print(f"[INFO] Basisverzeichnis: {BASE_DIR}")
print(f"[INFO] Datenverzeichnis: {DATA_DIR}")
print(f"[INFO] Ausgabeverzeichnis: {PLOTS_DIR}")


# Daten-Pfade definieren
qual_path = DATA_DIR / "genome_quality_report.tsv"  # CheckM2-Qualitätsdaten 
trna_path = DATA_DIR / "bvbrc_trna.tsv"            # matchings von IDs und tRNA-Counts

print(f"[INFO] Arbeitsverzeichnis: {Path.cwd()}")
print(f"[INFO] Datenverzeichnis: {DATA_DIR}")
print(f"[INFO] Ausgabeverzeichnis: {PLOTS_DIR}")

# Daten laden
qual = pd.read_csv(qual_path, sep="\t")
trna = pd.read_csv(trna_path, sep="\t")  # Auflistung beider Dateien als Pandas DataFrames

# Spalten harmonisieren (Benennung)
qual = qual.rename(columns={
    "Name": "genome_id",
    "Genome_Size": "genome_size_bp",         # Größe in bp
    "Total_Coding_Sequences": "gene_count"
})


# ID-Format vereinheitlichen: Unterstrich → Punkt (passt zu trna.tsv)
qual["genome_id"] = qual["genome_id"].str.replace("_", ".", regex=False)

# numerische Spalten erzwingen
for c in ["Completeness", "Contamination", "genome_size_bp", "gene_count"]:
    qual[c] = pd.to_numeric(qual[c], errors="coerce")

trna["genome_id"]   = trna["genome_id"].astype(str).str.strip()
trna["trna_count"]  = pd.to_numeric(trna["trna_count"], errors="coerce")

# fehlende Werte aus FASTA-Dateien ergänzen, öffne Fasta-Datei und zähle Anzahl der Gene durch Zählung der Header, wenn keine counts vorhanden sind. 
def fasta_stats(path: Path):
    length = genes = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                genes += 1
            else:
                length += len(line.strip())
    return length, genes


#print qual als tabelle 
fig, ax = plt.subplots(figsize=(12, 0.4 * len(qual)))  # Höhe dynamisch je nach Zeilenzahl
ax.axis('off')  # Achsen verbergen
tbl = ax.table(cellText=qual.values, colLabels=qual.columns, loc='center', cellLoc='left')
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
tbl.scale(1, 1.5)  # Spaltenbreite, Zeilenhöhe

plt.tight_layout()
plt.show()

# Speichere qual als CSV-Datei 
#qual.to_csv("qualtable_full.csv", index=False)

print("[INFO] Qualitätsdaten gespeichert als qualtable_full.csv")
print("[INFO] Fasta-Dateien werden ausgewertet...")

#Merge mit tRNA
df = qual.merge(trna, on="genome_id", how="left")

#  Qualitätsklassen nach CheckM2 nach Bowers et al. 2017
conditions = [
    (df["Completeness"] >= 90) & (df["Contamination"] <= 5),
    (df["Completeness"] >= 50) & (df["Contamination"] <= 10),
    (df["Completeness"] < 50) | (df["Contamination"] > 10) 
    ]
df["quality"] = np.select(conditions, ["high", "medium", "low"], default="low")

# Plot 1 QC-Scat
plt.figure(figsize=(8,6))  # Größerer Plot für bessere Sichtbarkeit
plt.scatter(
    df["Completeness"],
    df["Contamination"],
    c=df["quality"].map({"high":"green","medium":"orange","low":"red"}),
    alpha=0.7,  # Etwas transparent für überlappende Punkte
    s=60        # Größere Punkte
)
plt.axvline(90, ls="--", color="gray", label="90% Completeness threshold")
plt.axhline(5, ls="--", color="gray", label="5% Contamination threshold")
plt.xlabel("Completeness (%)")
plt.ylabel("Contamination (%)")
plt.title(f"CheckM2 Genome Quality Assessment according to Bowers et al. 2017 for BV-BRC genomes (n={len(df)})")

# Legende hinzufügen
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='green', label='High quality'),
    Patch(facecolor='orange', label='Medium quality'), 
    Patch(facecolor='red', label='Low quality')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.savefig("plot_completeness_vs_contamination.png", dpi=350, bbox_inches='tight')
plt.show()  # Zeigt den Plot an
print("[INFO] QC-Scatterplot gespeichert und angezeigt.")

#Korrelations-Hilfsfunktion
def correlation_plot(x, y, xlabel, ylabel, fname):
    ok = ~x.isna() & ~y.isna()
    if ok.sum() == 0:     # keine Daten
        return
    plt.figure(figsize=(5,4))
    plt.scatter(x[ok], y[ok], alpha=0.8)
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs {xlabel} (n={ok.sum()})")
    plt.tight_layout(); plt.savefig(fname, dpi=350); plt.close()

# Plot 2 & 3: tRNA-Korrelation ohne Logarithmierung
correlation_plot(
    df["genome_size_bp"]/1e6,
    df["trna_count"],
    "Genome Size (Mbp)",
    "tRNA Count",
    "genome_size_vs_trna_count.png"
)

correlation_plot(
    df["gene_count"],
    df["trna_count"],
    "Total Gene Count",
    "tRNA Count",
    "gene_count_vs_trna_count.png"
)

#Logarithmierung erscheint sinnvoll, da die Daten oft exponentiell verteilt sind. 
def loglog_correlation_plot(x, y, xlabel, ylabel, fname):
    """
    Scatter-Plot auf log10-Skalen + Power-Law-Trendlinie.
    y  ≈  C · x^a     (im log10-Raum  log10(y) = a·log10(x) + b)
    """
    # nur positive und nicht-NaN-Werte zulassen
    ok = (x > 0) & (y > 0) & (~x.isna()) & (~y.isna())
    if ok.sum() < 2:          # zu wenig Punkte für Fit
        return

    x_ok, y_ok = x[ok], y[ok]

    # --- Scatterplot --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(x_ok, y_ok, alpha=0.7)

    # --- Achsen log10 -------------------------------------------------------
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} vs {xlabel} (n={ok.sum()})")

    # --- Lineares Fit im log10-Raum ----------------------------------------
    logx = np.log10(x_ok)
    logy = np.log10(y_ok)
    a, b = np.polyfit(logx, logy, deg=1)        # Steigung & Achsenabschnitt
    C = 10**b                                   # C = 10^b

    # Power-Law-Fit hinzufügen  y_fit = C · x^a
    x_line = np.logspace(np.log10(x_ok.min()), np.log10(x_ok.max()), 200)
    y_line = C * x_line**a
    ax.plot(x_line, y_line, color="red", lw=2,
            label=f"y ≈ {C:.2g} · x^{a:.2f}")
    ax.legend()

    plt.tight_layout()
    plt.savefig(fname, dpi=350)
    plt.close()

#  (A) tRNA-Zahl vs Genome-Größe  (beide Achsen log10)
loglog_correlation_plot(
    df["genome_size_bp"] / 1e6,      # Mbp
    df["trna_count"],
    "Genome Size (Mbp, log₁₀)",
    "tRNA Count (log₁₀)",
    "linearcurve_logax_genome_size_vs_trna_count.png"
)

#  (B) tRNA-Zahl vs Gesamtgene  (beide Achsen log10)
loglog_correlation_plot(
    df["gene_count"],
    df["trna_count"],
    "Total Gene Count (log₁₀)",
    "tRNA Count (log₁₀)",
    "linearcurve_logax_gene_count_vs_trna_count.png"
)


#eventuell ist eine Log-Fit-Kurve sinnvoller als POWER-Law-Fit.

#In vielen Bakterien­genomen wächst die tRNA-Zahl langsamer als linear mit der Genomgröße.
#Ein logarith­misches Modell y=alog⁡(x)+by=alog(x)+b bildet einen solchen Sättigungs­trend
#natürlich ab.



def log_correlation_plot(x, y, xlabel, ylabel, fname):
    
    """
    Semi-log-Scatterplot + logarithmischer Fit
       Modell:  y = a · log10(x) + b
       Achsen:  x   log10-skaliert
                y   linear
    """
    
    ok = (x > 0) & (y > 0) & (~x.isna()) & (~y.isna())
    if ok.sum() < 2:
        return                        # zu wenig Daten für einen Fit

    x_ok, y_ok = x[ok], y[ok]


    # Scatter ­
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(x_ok, y_ok, alpha=0.7)

    #Achse auf log10-Skala stellen
    #ax.set_xscale("log")
    #ax.set_yscale("log")


    # Log-Fit 
    logx   = np.log10(x_ok)
    a, b   = np.polyfit(logx, y_ok, deg=1)      # y = a·log10(x) + b
    x_line = np.logspace(np.log10(x_ok.min()),
                         np.log10(x_ok.max()), 300)
    y_line = a * np.log10(x_line) + b
    ax.plot(x_line, y_line, color="red", lw=2,
            label=f"y ≈ {a:.2f}·log₁₀(x) + {b:.2f}")
    ax.legend()

    # Achsen & Titel 
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} vs {xlabel} (n={ok.sum()})")

    plt.tight_layout()
    plt.show()
    try:
        output_path = PLOTS_DIR / fname
        plt.savefig(output_path, dpi=350, bbox_inches='tight')
        print(f"[INFO] Plot gespeichert: {output_path}")
    except Exception as e:
        print(f"[ERROR] Konnte Plot {fname} nicht speichern: {e}")
    
    plt.close() 
    

# (A) tRNA-Zahl vs. Genom-Größe
log_correlation_plot(
    df["genome_size_bp"] / 1e6,     # Mbp
    df["trna_count"],
    "log 10 (Genome Size) (Mbp)",
    "log10 (BV-BRC predicted number of tRNA)",
    "log_genome_size_vs_trna_count.png"
)
print("[INFO] Plot tRNA-Anzahl vs Genomgrösse saved!")

# (B) tRNA-Zahl vs. Gesamt­genanzahl
log_correlation_plot(
    df["gene_count"],
    df["trna_count"],
    "Total Gene Count",
    "BV-BRC predicted number of tRNA",
    "log_gene_count_vs_trna_count.png"
)

print("[INFO] Plot tRNA-Anzahl vs Gesamtgenzahl saved!")
print("[INFO] Alle Plots gespeichert.")


#1 tRNA-Count schein verloren gegangen zu sein, wo ist er ? 
# ► IDs in beiden Ausgangstabellen
ids_qual = set(qual["genome_id"])
ids_trna = set(trna["genome_id"])

# ► Schnittmenge, Differenzen
common   = ids_qual & ids_trna          # sollte 40 sein
only_q   = ids_qual - ids_trna          # in qual, aber NICHT in trna
only_t   = ids_trna - ids_qual          # in trna, aber NICHT in qual

print("Gemeinsame Genome_IDs in qual und trna :", len(common))
print("Nur in qual    :", only_q)
print("Nur in trna    :", only_t)
#fixed
#edit: 2599821_20 in qual umbenennen in 2599821_2, sonst missmatch!!!

print(f"Plots will be saved in: {PLOTS_DIR}")
try:
    output_path = PLOTS_DIR / "genome_size_vs_trna_count.png"
    plt.savefig(output_path, dpi=350)
    print(f"Plot saved successfully: {output_path}")
except Exception as e:
    print(f"Error saving plot: {e}")