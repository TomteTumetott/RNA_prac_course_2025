import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Importiere zentrale Konfiguration
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, PLOTS_DIR

# Eingabedatei
input_file = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct.txt"

# Einlesen
df = pd.read_csv(input_file, sep="\t")

# IntronPairs sicherstellen (eval nur wenn String)
df["IntronPairs"] = df["IntronPairs"].apply(lambda x: eval(x) if isinstance(x, str) else x)

# Funktion: relative Intron-Startposition (tRNA-intern, 0-basiert)
def get_intron_positions(row):
    start, end = row["Begin"], row["End"]
    strand = 1 if end > start else -1
    mature_range = list(range(start, end + strand, strand))
    mature_map = {gpos: i for i, gpos in enumerate(mature_range)}
    intron_positions = []
    for gstart, _ in row["IntronPairs"]:
        if gstart in mature_map:
            rel_pos = mature_map[gstart]
            if 0 <= rel_pos < 100:
                intron_positions.append(rel_pos)
    return intron_positions

# Berechne relative Positionen
df["RelativeIntronPositions"] = df.apply(get_intron_positions, axis=1)

# Explodiere nach Positionen
exploded = df.explode("RelativeIntronPositions")
exploded = exploded.dropna(subset=["RelativeIntronPositions"])
exploded["RelativeIntronPositions"] = exploded["RelativeIntronPositions"].astype(int)

# Pivot-Tabelle: Zeilen = Superphylum (Phylum), Spalten = Position
heatmap_data = exploded.pivot_table(
    index="Superphylum",
    columns="RelativeIntronPositions",
    aggfunc="size",
    fill_value=0
)

# Nur Positionen 0–99 (für einheitliche Darstellung)
heatmap_data = heatmap_data.reindex(columns=range(100), fill_value=0)

# Überprüfe den Wertebereich für sinnvolle Ticks
max_value = heatmap_data.values.max()
print(f"[DEBUG] Maximum value in heatmap: {max_value}")

# Bestimme sinnvolle Tick-Schritte basierend auf dem Wertebereich
if max_value <= 10:
    tick_step = 1
elif max_value <= 50:
    tick_step = 5
elif max_value <= 100:
    tick_step = 10
elif max_value <= 500:
    tick_step = 50
else:
    tick_step = 100

# Erstelle Ticks
colorbar_ticks = np.arange(0, max_value + 1, tick_step)
if colorbar_ticks[-1] < max_value:
    colorbar_ticks = np.append(colorbar_ticks, max_value)

print(f"[DEBUG] Colorbar ticks: {colorbar_ticks}")

# Plot erstellen
plt.figure(figsize=(18, 6))
ax = sns.heatmap(
    heatmap_data,
    cmap="plasma",  # dunkles = viele Introns, hell = wenige
    cbar_kws={
        "label": "Number of tRNAs with Introns",
        "ticks": colorbar_ticks
    },
    linewidths=0.2,
    linecolor="black"
)

# x-Achsenbeschriftung zentrieren
xtick_locs = np.arange(heatmap_data.shape[1]) + 0.5
ax.set_xticks(xtick_locs)
ax.set_xticklabels(heatmap_data.columns)

# Achsentitel
ax.set_xlabel("Position along mature tRNA", fontsize=12)
ax.set_ylabel("Phylum", fontsize=12)
ax.set_title("Positional Heatmap of Intron Locations (absolute tRNA coordinates)", fontsize=16)

# Layout optimieren
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=9)
plt.tight_layout()

# Speichern und anzeigen
plt.savefig(PLOTS_DIR / "heatmap_relative_intron_positions.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"[INFO] Heatmap saved to: {PLOTS_DIR / 'heatmap_relative_intron_positions.png'}")



###############################################################################################
# MeanIsotypeScore 
###############################################################################################
# Datei einlesen
input_file = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct.txt"
df = pd.read_csv(input_file, sep="\t")

# Kombination sicherstellen
if "Type_Anticodon" not in df.columns:
    df["Type_Anticodon"] = df["Type"].astype(str) + "_" + df["Anticodon"].astype(str)

# Gruppieren
summary = (
    df.groupby(["Superphylum", "Type_Anticodon"])["IsotypeScore"]
    .agg(["mean", "std"])
    .reset_index()
)

# Vorbereitung
superphyla = summary["Superphylum"].unique()
fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.subplots_adjust(top=0.92, hspace=0.25)  # Etwas Platz nach oben schaffen

axes = axes.flatten()
bar_color = "#4c72b0"

# Ausreißer vorbereiten
for i, sp in enumerate(superphyla):
    ax = axes[i]
    sub_summary = summary[summary["Superphylum"] == sp]
    sub_data = df[df["Superphylum"] == sp]

    bars = ax.bar(
        sub_summary["Type_Anticodon"],
        sub_summary["mean"],
        color = bar_color,
        edgecolor="black"
    )

    # Kleine SD-Balken am oberen Ende
    for bar, std in zip(bars, sub_summary["std"]):
        height = bar.get_height()
        center = bar.get_x() + bar.get_width() / 2
        ax.plot([center, center], [height, height + std], color="orange", linewidth=1.2)
        ax.plot([center - 0.1, center + 0.1], [height + std]*2, color="orange", linewidth=1)

    # Streupunkte (Ausreißer)
    #sns.stripplot(
    #    data=sub_data,
    #    x="Type_Anticodon",
    #    y="IsotypeScore",
    #    ax=ax,
    #    color="black",
    #    size=3,
    #    jitter=True,
    #    alpha=0.5
    #)

    ax.set_title(f"Mean Isotype Score per Anticodon with tRNA Type and Std Dev: {sp}")
    ax.set_xlabel("Anticodon (tRNA-Type)")
    ax.set_ylabel("Mean Isotype Score")
    ax.tick_params(axis='x', rotation=90, labelsize=6)

plt.suptitle("Isotype Scores per Type-Anticodon", fontsize=16, y=0.98)


plt.savefig(PLOTS_DIR / "barplot_isotype_scores_by_superphylum.png", dpi=300, bbox_inches="tight")
plt.show()



##########################################################
#average mature anticodon size
##########################################################
# Datei einlesen
input_file = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct.txt"
df = pd.read_csv(input_file, sep="\t")

# Sicherstellen, dass die notwendigen Spalten vorhanden sind
if "Type_Anticodon" not in df.columns:
    df["Type_Anticodon"] = df["Type"].astype(str) + "_" + df["Anticodon"].astype(str)

# Berechnung der tRNA-Länge aus End und Begin
df["tRNA_Length"] = (df["End"].astype(int) - df["Begin"].astype(int)).abs() + 1

# Gruppierung für Mittelwert und Standardabweichung
summary = (
    df.groupby(["Superphylum", "Type_Anticodon"])["tRNA_Length"]
    .agg(["mean", "std"])
    .reset_index()
)

# Ploteinstellungen
superphyla = summary["Superphylum"].unique()
fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.subplots_adjust(top=0.9, hspace= 0.25)  # Etwas Platz nach oben schaffen
fig.suptitle("tRNA-Anticodon Lengths by Superphylum", fontsize=16, y=0.98)


axes = axes.flatten()
bar_color = "#419abd"

# Plot für jedes Superphylum
for i, sp in enumerate(superphyla):
    ax = axes[i]
    sub_summary = summary[summary["Superphylum"] == sp]
    sub_data = df[df["Superphylum"] == sp]

    bars = ax.bar(
        sub_summary["Type_Anticodon"],
        sub_summary["mean"],
        color=bar_color,
        edgecolor="black"
    )

    # Fehlerbalken für Standardabweichung
    for bar, std in zip(bars, sub_summary["std"]):
        height = bar.get_height()
        center = bar.get_x() + bar.get_width() / 2
        ax.plot([center, center], [height, height + std], color="orange", linewidth=1.2)
        ax.plot([center - 0.1, center + 0.1], [height + std]*2, color="orange", linewidth=1)

    # Einzelpunkte (Ausreißer)
    sns.stripplot(
        data=sub_data,
        x="Type_Anticodon",
        y="tRNA_Length",
        ax=ax,
        color="black",
        size=3,
        jitter=True,
        alpha=0.5
    )

    ax.set_title(f"Average mature tRNA-anticodon size: {sp}")
    ax.set_xlabel("tRNA-Anticodon")
    ax.set_ylabel("Average mature tRNA-anticodon size (nts)")
    ax.tick_params(axis='x', rotation=90, labelsize=6)
    
plt.savefig(PLOTS_DIR / "barplot_trna_lengths_by_superphylum.png", dpi=300, bbox_inches="tight")
plt.show()




#####################################################################
# Average IntronCount per Superphyla
#####################################################################

if "Type_Anticodon" not in df.columns:
    df["Type_Anticodon"] = df["Type"].astype(str) + "_" + df["Anticodon"].astype(str)

# IntronCount sicherstellen
df["IntronCount"] = pd.to_numeric(df["IntronCount"], errors="coerce")

# Gruppierung: Mittelwert und Standardabweichung der IntronCount pro Type_Anticodon und Superphylum
summary = (
    df.groupby(["Superphylum", "Type_Anticodon"])["IntronCount"]
    .agg(["mean", "std"])
    .reset_index()
)

# Ploteinstellungen
superphyla = summary["Superphylum"].unique()
fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.subplots_adjust(top=0.9, hspace=0.25)
fig.suptitle("Average Intron Number per tRNA-Anticodon by Superphylum", fontsize=16, y=0.98)

axes = axes.flatten()
bar_color = "#41bdb7"

# Plot für jedes Superphylum
for i, sp in enumerate(superphyla):
    ax = axes[i]
    sub_summary = summary[summary["Superphylum"] == sp]
    sub_data = df[df["Superphylum"] == sp]

    bars = ax.bar(
        sub_summary["Type_Anticodon"],
        sub_summary["mean"],
        color=bar_color,
        edgecolor="black"
    )

    # Fehlerbalken für Standardabweichung
    for bar, std in zip(bars, sub_summary["std"]):
        height = bar.get_height()
        center = bar.get_x() + bar.get_width() / 2
        ax.plot([center, center], [height, height + std], color="orange", linewidth=1.2)
        ax.plot([center - 0.1, center + 0.1], [height + std]*2, color="orange", linewidth=1)

    # Einzelpunkte (Ausreißer)
    #sns.stripplot(
    #    data=sub_data,
    #    x="Type_Anticodon",
    #    y="IntronCount",
    #    ax=ax,
    #    color="black",
    #    size=3,
    #    jitter=True,
    #    alpha=0.5
    #)

    ax.set_title(f"Avg. Intron Count per tRNA-Anticodon: {sp}")
    ax.set_xlabel("tRNA-Anticodon")
    ax.set_ylabel("Average Intron Count")
    ax.tick_params(axis='x', rotation=90, labelsize=6)

plt.savefig(PLOTS_DIR / "barplot_introncount_by_superphylum.png", dpi=300, bbox_inches="tight")
plt.show()

