#!/usr/bin/env python3
# Archaeal tRNA Structure Heatmap Generator als sns seaborn heatmap

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# config.py einbinden
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, PLOTS_DIR

# Einstellbare Parameter am Anfang des Skripts (nach den Imports)
# Diese können angepasst werden, um die Plotgröße zu optimieren
INCH_PER_GENOME = 0.15  # Höhe pro Genom in inches
INCH_PER_ANTICODON = 0.15  # Breite pro Anticodon in inches
MAX_FIG_HEIGHT = 60  # Maximale Höhe in inches
MAX_FIG_WIDTH = 40   # Maximale Breite in inches
MIN_FIG_HEIGHT = 8   # Minimale Höhe in inches
MIN_FIG_WIDTH = 12   # Minimale Breite in inches

# Eingabedatei
input_file = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct_annotated.txt"

# Einlesen
df = pd.read_csv(input_file, sep="\t")
print("Spaltenübersicht:", df.columns.tolist())

# Debug: Zeige die tatsächliche Struktur der Daten
print("\nBeispiel der ersten 3 Zeilen:")
print(df.head(3).to_string())

# Überprüfe welche Spalten für Name und Taxonomie vorhanden sind
if 'Name' in df.columns:
    print(f"\nBeispiel Namen: {df['Name'].unique()[:5].tolist()}")
elif 'name' in df.columns:
    # Falls 'name' statt 'Name' existiert, umbenennen
    df.rename(columns={'name': 'Name'}, inplace=True)
    print(f"\nBeispiel Namen (aus 'name' Spalte): {df['Name'].unique()[:5].tolist()}")
    
if 'Superphylum' in df.columns:
    print(f"\nSuperphylum Werte: {df['Superphylum'].unique()}")
    print(f"Anzahl pro Superphylum: \n{df['Superphylum'].value_counts()}")

# Überprüfe ob Type_Anticodon Spalte existiert
if 'Type_Anticodon' not in df.columns:
    print("\n[WARNUNG] 'Type_Anticodon' Spalte nicht gefunden.")
    if 'Type' in df.columns and 'Anticodon' in df.columns:
        # Erstelle Type_Anticodon aus Type und Anticodon
        df['Type_Anticodon'] = df['Type'].astype(str) + '_' + df['Anticodon'].astype(str)
        print("Erstelle 'Type_Anticodon' aus 'Type' und 'Anticodon' Spalten")
    else:
        print("[FEHLER] Kann Type_Anticodon nicht erstellen - fehlende Spalten")

# Strukturklassifikation
def classify(stem_count):
    if stem_count == 4:
        return "4-stems"
    elif stem_count == 5:
        return "5-stems" 
    elif stem_count in [0, 1, 2, 3, 6, 7]:
        return "out-ss"
    return "absence"

df["StructureClass"] = df["StemCount"].apply(classify)

# Taxonomische Gruppierung definieren
# Option 1: Falls eine Spalte 'Superphylum' oder 'TaxGroup' existiert
if 'Superphylum' in df.columns:
    df['TaxonomicGroup'] = df['Superphylum']
    print("\nVerwende 'Superphylum' Spalte für taxonomische Gruppierung")
    print(f"Gefundene Gruppen: {df['Superphylum'].unique()}")
    
    # Zeige die Verteilung
    print("\nAnzahl Einträge pro Gruppe:")
    print(df['Superphylum'].value_counts().to_string())
elif 'TaxGroup' in df.columns:
    df['TaxonomicGroup'] = df['TaxGroup']
    print("Verwende 'TaxGroup' Spalte für taxonomische Gruppierung")
    
    df['TaxonomicGroup'] = df['Name'].apply(assign_taxonomic_group)
    
    # Manuelle Korrekturen für spezifische Genome (falls nötig)
    manual_corrections = {
        # Beispiel: 'Candidatus_Heimdallarchaeota_archaeon_S012_metabat1_scaf2bin.003': 'Asgard',
        # Füge hier problematische Genome hinzu
    }
    
    for genome, group in manual_corrections.items():
        df.loc[df['Name'] == genome, 'TaxonomicGroup'] = group
    
    # Zeige Gruppenzuordnung
    print("\nGruppenzuordnung:")
    print(df['TaxonomicGroup'].value_counts())

# Pivot: Zeilen = Name (Genome), Spalten = Type_Anticodon
if 'Name' not in df.columns:
    print("[FEHLER] 'Name' Spalte nicht gefunden für Pivot!")
    # Versuche alternative Spalten
    if 'Sequence' in df.columns:
        print("Verwende 'Sequence' als Index")
        index_col = 'Sequence'
    else:
        print("Verwende erste Spalte als Index")
        index_col = df.columns[0]
else:
    index_col = 'Name'

heatmap_df = df.pivot_table(
    index=index_col, 
    columns="Type_Anticodon", 
    values="StructureClass", 
    aggfunc="first", 
    fill_value="absence"
)

# Taxonomische Gruppen für Sortierung
tax_groups_df = df[[index_col, 'TaxonomicGroup']].drop_duplicates()
tax_groups_df = tax_groups_df.set_index(index_col)

# Sortiere Genome nach taxonomischer Gruppe und dann alphabetisch
genome_order = []
group_positions = {}
current_pos = 0

# Stelle sicher, dass alle Gruppen berücksichtigt werden
all_groups = ['Asgard', 'DPANN', 'Euryarchaeota', 'TACK']
for group in all_groups:
    group_genomes = tax_groups_df[tax_groups_df['TaxonomicGroup'] == group].index.tolist()
    group_genomes.sort()  # Alphabetisch sortieren
    
    if group_genomes:
        group_positions[group] = (current_pos, current_pos + len(group_genomes) - 1)
        current_pos += len(group_genomes)
        genome_order.extend(group_genomes)
        print(f"{group}: {len(group_genomes)} Genome")
    else:
        print(f"[WARNUNG] Keine Genome in Gruppe {group} gefunden!")

# Füge 'Other' Gruppe hinzu, falls vorhanden
other_genomes = tax_groups_df[tax_groups_df['TaxonomicGroup'] == 'Other'].index.tolist()
if other_genomes:
    other_genomes.sort()
    group_positions['Other'] = (current_pos, current_pos + len(other_genomes) - 1)
    genome_order.extend(other_genomes)
    print(f"Other: {len(other_genomes)} Genome")

# Sortiere Heatmap nach Genomordnung
heatmap_df = heatmap_df.reindex(genome_order)

# Sortiere Spalten alphabetisch
heatmap_df = heatmap_df.reindex(sorted(heatmap_df.columns), axis=1)

# Mapping auf Farbklassen
structure_order = ["4-stems", "5-stems", "absence", "out-ss"]
color_map = {
    "4-stems": "#fdbf6f",   # orange
    "5-stems": "#a6cee3",   # hellblau
    "out-ss": "#1b9e77",    # dunkelgrün
    "absence": "#ffffff"    # weiß 
}

# Erstelle eine Liste der Farben in der richtigen Reihenfolge
colors = [color_map[k] for k in structure_order]

# Werte in numerisch übersetzen
class_to_int = {cls: idx for idx, cls in enumerate(structure_order)}
numeric_df = heatmap_df.replace(class_to_int)

# Dynamische Größenberechnung basierend auf Anzahl der Genome und Anticodons
n_genomes = len(heatmap_df)
n_anticodons = len(heatmap_df.columns)

# Berechne optimale Figurgrößen
fig_height = max(MIN_FIG_HEIGHT, min(MAX_FIG_HEIGHT, n_genomes * INCH_PER_GENOME))
fig_width = max(MIN_FIG_WIDTH, min(MAX_FIG_WIDTH, n_anticodons * INCH_PER_ANTICODON + 4))

print(f"Erstelle Heatmap für {n_genomes} Genome und {n_anticodons} Anticodons")
print(f"Figurgröße: {fig_width:.1f} x {fig_height:.1f} inches")

# Plot mit dynamischer Größe
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Heatmap erstellen mit angepassten Einstellungen
sns.heatmap(
    numeric_df,
    cmap=colors,
    linewidths=0.1 if n_genomes > 100 else 0.3,  # Dünnere Linien bei vielen Genomen
    linecolor='black',
    cbar=False,
    vmin=0,
    vmax=len(structure_order)-1,
    ax=ax,
    rasterized=n_genomes > 200  # Rasterisierung bei sehr großen Datensätzen
)

# Titel und Achsenbeschriftungen mit angepassten Schriftgrößen
title_size = min(16, max(12, fig_width * 0.5))
label_size = min(12, max(10, fig_width * 0.4))
tick_size = min(8, max(6, fig_height * 0.15))

#x-Achsenbeschriftung zentrieren auf die Mitte der Felder 
xtick_positions = np.arange(len(numeric_df.columns)) + 0.5
ax.set_xticks(xtick_positions)
ax.set_xticklabels(numeric_df.columns, rotation=90, ha='center', fontsize=tick_size)

ax.set_title("tRNA results archaeal genome: tRNA Presence/Absence categories", 
             fontsize=title_size, pad=20)
ax.set_xlabel("tRNA anticodons", fontsize=label_size)
ax.set_ylabel("Archaeal genomes", fontsize=label_size)

# X-Achse Rotation mit angepasster Schrift
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=tick_size)
#ax.set_yticklabels(ax.get_yticklabels(), fontsize=tick_size)

# Horizontale Trennlinien zwischen Gruppen
for group, (start, end) in group_positions.items():
    if end < len(genome_order) - 1:  # Nicht nach der letzten Gruppe
        ax.axhline(y=end + 1, color='black', linewidth=2)

# Gruppenlabels am rechten Rand mit angepasster Schriftgröße
group_label_size = min(10, max(8, fig_height * 0.2))
for group, (start, end) in group_positions.items():
    middle = (start + end) / 2
    ax.text(len(heatmap_df.columns) + 0.5, middle, group, 
            verticalalignment='center', fontsize=group_label_size, weight='bold')

# Legende erstellen
legend_elements = [
    Line2D([0], [0], marker='s', color='w', label='4-stems',
           markerfacecolor=color_map["4-stems"], markersize=15),
    Line2D([0], [0], marker='s', color='w', label='5-stems',
           markerfacecolor=color_map["5-stems"], markersize=15),
    Line2D([0], [0], marker='s', color='w', label='absence',
           markerfacecolor=color_map["absence"], markersize=15, markeredgecolor='black'),
    Line2D([0], [0], marker='s', color='w', label='out-ss',
           markerfacecolor=color_map["out-ss"], markersize=15)
]

# Legende positionieren
legend = ax.legend(handles=legend_elements, 
                   title='tRNA secondary structures',
                   title_fontsize=12,
                   fontsize=11,
                   loc='center left',
                   bbox_to_anchor=(1.05, 0.5),
                   frameon=True)

# Falls zu viele Genome, zeige nur jedes n-te Label
if n_genomes > 50:
    # Behalte alle Labels, aber mache sie kleiner
    y_labels = ax.get_yticklabels()
    # Optional: Zeige nur jedes 2. oder 3. Label
    # for i, label in enumerate(y_labels):
    #     if i % 2 != 0:  # Zeige nur jedes 2. Label
    #         label.set_visible(False)
    
# Falls zu viele Anticodons, rotiere stärker
#if n_anticodons > 40:
#    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right', fontsize=tick_size-1)

# Anpassung für bessere Sichtbarkeit
plt.subplots_adjust(left=0.15, right=0.90, top=0.95, bottom=0.10)

# Layout optimieren
try:
    plt.tight_layout()
except:
    pass  # Falls tight_layout fehlschlägt, verwende die manuellen Einstellungen

# Speichern mit höherer DPI für bessere Qualität
output_path = PLOTS_DIR / "archaeal_trna_heatmap.png"

# Bei sehr großen Plots, reduziere DPI um Speicher zu sparen
save_dpi = 300 if n_genomes < 100 else 200

plt.savefig(output_path, dpi=save_dpi, bbox_inches='tight')
print(f"[INFO] Heatmap gespeichert unter: {output_path} (DPI: {save_dpi})")

# Statistiken ausgeben
print(f"\nAnzahl Genome: {len(heatmap_df)}")
print(f"Anzahl Anticodons: {len(heatmap_df.columns)}")
print(f"\nGenome pro Gruppe:")
for group, (start, end) in group_positions.items():
    print(f"  {group}: {end - start + 1}")
print(f"\nStrukturverteilung:")
print(df["StructureClass"].value_counts())

# Debug: Zeige Genome ohne ASGARD-Zuordnung, die möglicherweise ASGARD sein sollten
print("\n[DEBUG] Genome in 'Other' Gruppe (möglicherweise falsch zugeordnet):")
other_df = tax_groups_df[tax_groups_df['TaxonomicGroup'] == 'Other']
if not other_df.empty:
    for genome in other_df.index[:10]:  # Erste 10 zeigen
        print(f"  - {genome}")

plt.show()


