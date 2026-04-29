#!/usr/bin/env python3
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TREE_DIR, OUTPUT_DIR

# Funktion zum Laden des Newick-Strings
def load_newick_from_tree_file(tree_file: Path) -> str:
    #Lädt den Newick-String aus einer .tree-Datei (erste nicht-leere Zeile).
    
    with tree_file.open() as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    raise ValueError("Kein Newick-String in der Datei gefunden.")

# Newick-Dateipfad (PLEASE ADAPT MANUALLY)
TREE_FILE = OUTPUT_DIR / 'mlocarna_results' / 'tRNA_GAT_4stem_candidates' / 'results' / 'result.tree'
#TREE_FILE = OUTPUT_DIR / 'mlocarna_results' / 'tRNA_TAT_4stem_candidates' / 'results' / 'result.tree'


# Metadaten-Dateipfad
METADATA_FILE = OUTPUT_DIR / 'TABLE_metadata_Isoleucine.txt'

# ANTICODON-TYPE für Plot-Benennung (ADAPT manually!)
ANTICODON_TYPE = "Ile_GAT"  
#ANTICODON_TYPE = "Ile_TAT"  

def load_metadata() -> pd.DataFrame:
    #Lädt die Metadatentabelle aus der Datei
    
    if not METADATA_FILE.exists():
        raise FileNotFoundError(f"Metadaten-Datei nicht gefunden: {METADATA_FILE}")
    
    print(f"[INFO] Lade Metadaten aus: {METADATA_FILE}")
    metadata_df = pd.read_csv(METADATA_FILE, sep='\t')
    print(f"[INFO] {len(metadata_df)} Einträge geladen")
    print(f"[INFO] Verfügbare Anticodons in Metadaten: {sorted(metadata_df['tRNA-anticodon'].unique())}")
    print(f"[INFO] Verwende manuell gesetzten Anticodon-Type: {ANTICODON_TYPE}")
    
    # Debug: Zeige erste paar TipNames
    print(f"[DEBUG] Beispiel TipNames in Metadaten:")
    for i, tip in enumerate(metadata_df['TipName'].head(3)):
        print(f"  {i+1}: {tip}")
    
    return metadata_df

def debug_tree_names(tree, metadata_df):
    #Debug-Funktion: Vergleicht Baum-Namen mit Metadaten-Namen
    
    leaf_names = [leaf.name for leaf in tree.get_leaves()]
    metadata_tipnames = set(metadata_df['TipName'].tolist())
    
    print(f"\n[DEBUG] Baum-Analyse:")
    print(f"  Anzahl Blätter im Baum: {len(leaf_names)}")
    print(f"  Anzahl TipNames in Metadaten: {len(metadata_tipnames)}")
    
    print(f"\n[DEBUG] Erste paar Blatt-Namen im Baum:")
    for i, name in enumerate(leaf_names[:3]):
        print(f"  {i+1}: '{name}'")
    
    # Prüfe Übereinstimmungen
    matches = set(leaf_names) & metadata_tipnames
    tree_only = set(leaf_names) - metadata_tipnames
    meta_only = metadata_tipnames - set(leaf_names)
    
    print(f"\n[DEBUG] Namens-Übereinstimmung:")
    print(f"  Übereinstimmende Namen: {len(matches)}")
    print(f"  Nur im Baum: {len(tree_only)}")
    print(f"  Nur in Metadaten: {len(meta_only)}")
    
    if tree_only:
        print(f"\n[DEBUG] Namen nur im Baum (erste 3):")
        for name in list(tree_only)[:3]:
            print(f"  '{name}'")
    
    if meta_only:
        print(f"\n[DEBUG] Namen nur in Metadaten (erste 3):")
        for name in list(meta_only)[:3]:
            print(f"  '{name}'")
    
    # Erstelle Name-Mapping basierend auf Übereinstimmungen
    if len(matches) > 0:
        # Direkte Matches gefunden - erstelle 1:1 Mapping
        name_mapping = {name: name for name in matches}
        print(f"\n[INFO] Direkte Übereinstimmungen gefunden: {len(name_mapping)}")
        return True, name_mapping
    else:
        # Keine direkten Matches - versuche flexible Zuordnung
        print(f"\n[INFO] Keine direkten Übereinstimmungen - versuche flexible Zuordnung...")
        name_mapping = try_name_matching(leaf_names, metadata_tipnames)
        return len(name_mapping) > 0, name_mapping
    """
    #Versucht flexible Namens-Zuordnung zwischen Baum und Metadaten
    
    print(f"\n[INFO] Versuche flexible Namens-Zuordnung...")
    
    # Erstelle Mapping-Dictionary
    name_mapping = {}
    metadata_tipnames_list = list(metadata_tipnames)
    
    # Strategie 1: Exakte Übereinstimmung
    for tree_name in tree_names:
        if tree_name in metadata_tipnames:
            name_mapping[tree_name] = tree_name
    
    # Strategie 2: Case-insensitive Matching (speziell .trna vs .tRNA)
    for tree_name in tree_names:
        if tree_name not in name_mapping:
            # Ersetze .trna durch .tRNA (kleines t zu großem T)
            if '.trna' in tree_name:
                trna_name = tree_name.replace('.trna', '.tRNA')
                if trna_name in metadata_tipnames:
                    name_mapping[tree_name] = trna_name
                    print(f"[DEBUG] Case-Matching: '{tree_name}' -> '{trna_name}'")
            
            # Umgekehrt: .tRNA zu .trna
            elif '.tRNA' in tree_name:
                trna_name = tree_name.replace('.tRNA', '.trna')
                if trna_name in metadata_tipnames:
                    name_mapping[tree_name] = trna_name
                    print(f"[DEBUG] Case-Matching: '{tree_name}' -> '{trna_name}'")
    
    # Strategie 3: Entferne Quotes falls vorhanden
    for tree_name in tree_names:
        if tree_name not in name_mapping:
            clean_name = tree_name.strip('\'"')
            if clean_name in metadata_tipnames:
                name_mapping[tree_name] = clean_name
                print(f"[DEBUG] Quote-Matching: '{tree_name}' -> '{clean_name}'")
    
    # Strategie 4: Ersetze Punkte durch Unterstriche
    for tree_name in tree_names:
        if tree_name not in name_mapping:
            underscore_name = tree_name.replace('.', '_')
            if underscore_name in metadata_tipnames:
                name_mapping[tree_name] = underscore_name
                print(f"[DEBUG] Punkt-Matching: '{tree_name}' -> '{underscore_name}'")
    
    # Strategie 5: Entferne alles nach dem letzten Punkt (falls es eine Erweiterung ist)
    for tree_name in tree_names:
        if tree_name not in name_mapping:
            if '.' in tree_name:
                base_name = tree_name.rsplit('.', 1)[0]
                if base_name in metadata_tipnames:
                    name_mapping[tree_name] = base_name
                    print(f"[DEBUG] Extension-Matching: '{tree_name}' -> '{base_name}'")
    
    print(f"[INFO] Erfolgreiche Zuordnungen: {len(name_mapping)} von {len(tree_names)}")
    
    if len(name_mapping) > 0:
        print(f"[DEBUG] Erste 3 erfolgreiche Zuordnungen:")
        for i, (tree_name, meta_name) in enumerate(list(name_mapping.items())[:3]):
            print(f"  {i+1}: '{tree_name}' -> '{meta_name}'")
    
    return name_mapping
    """
def create_tree_only():
    #Erstellt nur den phylogenetischen Baum mit gefärbten Ästen (ohne Heatmap)
    
    # Parse Metadaten
    metadata_df = load_metadata()
    
    # Verwende manuellen Anticodon-Type für Titel und Dateinamen
    anticodon_title = ANTICODON_TYPE
    
    # Erstelle Metadaten-Dictionary
    metadata_dict = {}
    for _, row in metadata_df.iterrows():
        metadata_dict[row['TipName']] = {
            'cca_enzyme': row['CCA-enzyme'],
            'superphyla': row['superphyla'],
            'cca_encoded': row['CCA-encoded'],
            'species_id': row['Species_ID']
        }
    
    # Farbschema für Superphyla
    colors = {
        'Asgard': '#e41a1c',      # Rot
        'DPANN': "#459dcc",        # Blau
        'Euryarchaeota': "#3ebc3a", # Grün
        'TACK': "#9533a3"          # Lila
    }
    
    # Lade Baum
    NEWICK_TREE = load_newick_from_tree_file(TREE_FILE)
    NEWICK_TREE = NEWICK_TREE.replace("\\.", ".")
    tree = Tree(NEWICK_TREE, format=1)
    
    # TreeStyle für rechteckige Darstellung
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "r"  # Rechteckig
    ts.scale = 20
    ts.branch_vertical_margin = 5
    ts.show_scale = False
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 1
    ts.guiding_lines_color = "#aaaaaa"
    
    # Titel für den gesamten Baum
    title_face = TextFace(f"tRNA phylogenetic tree {anticodon_title}", fsize=14, bold=True)
    ts.title.add_face(title_face, column=0)
    
    # Funktion um Superphylum eines Blattes zu finden
    def get_leaf_superphylum(node):
        #Findet das Superphylum für einen Knoten basierend auf seinen Blättern
        if node.is_leaf():
            if node.name in metadata_dict:
                return metadata_dict[node.name]['superphyla']
        else:
            # Sammle alle Superphyla der Blätter
            superphyla = set()
            for leaf in node.get_leaves():
                if leaf.name in metadata_dict:
                    superphyla.add(metadata_dict[leaf.name]['superphyla'])
            # Wenn alle Blätter das gleiche Superphylum haben, gib es zurück
            if len(superphyla) == 1:
                return list(superphyla)[0]
        return None
    
    # Style für jeden Knoten
    for node in tree.traverse():
        nstyle = NodeStyle()
        
        # Bestimme Farbe basierend auf Superphylum
        superphylum = get_leaf_superphylum(node)
        
        if node.is_leaf() and node.name in metadata_dict:
            meta = metadata_dict[node.name]
            
            # Setze Ast-Farbe
            if superphylum:
                nstyle["vt_line_color"] = colors[superphylum]
                nstyle["hz_line_color"] = colors[superphylum]
                nstyle["vt_line_width"] = 3
                nstyle["hz_line_width"] = 3
            
            nstyle["size"] = 0  # Kein Knoten-Symbol für Blätter
            
            # Label mit Species Info
            label = f"{meta['superphyla']}_{meta['species_id']}"
            label_color = colors.get(meta['superphyla'], "black")
            name_face = TextFace(label, fsize=10, fgcolor=label_color, bold=True)
            node.add_face(name_face, column=0, position="aligned")
            
        else:
            # Interne Knoten
            nstyle["size"] = 6
            nstyle["fgcolor"] = "black"
            
            # Färbe Äste wenn möglich
            if superphylum:
                nstyle["vt_line_color"] = colors[superphylum]
                nstyle["hz_line_color"] = colors[superphylum]
                nstyle["vt_line_width"] = 3
                nstyle["hz_line_width"] = 3
        
        node.set_style(nstyle)
    
    # Speichere Baum
    output_path_tree = TREE_DIR / f"trna_tree_only_{anticodon_title}.png"
    tree.render(str(output_path_tree), tree_style=ts, dpi=400, w=800, h=1400)
    print(f"[INFO] Baum gespeichert als: {output_path_tree}")
    
    # Extrahiere die korrekte Blatt-Reihenfolge aus dem gerenderten Baum
    leaf_order = [leaf.name for leaf in tree.get_leaves()]
    
    # Rückgabe der Blatt-Reihenfolge für Heatmap
    return leaf_order

def create_heatmap_only(leaf_order):
    #Erstellt nur die Heatmap mit matplotlib/seaborn
    
    # Parse Metadaten
    metadata_df = load_metadata()
    
    # Verwende manuellen Anticodon-Type für Titel und Dateinamen
    anticodon_title = ANTICODON_TYPE
    
    # Erstelle Metadaten-Dictionary
    metadata_dict = {}
    for _, row in metadata_df.iterrows():
        # Bestimme CCA1, CCA2, CCA12 Status
        cca_status = {'CCA1': 0, 'CCA2': 0, 'CCA12': 0}
        if row['CCA-enzyme'] == 'CCA1':
            cca_status['CCA1'] = 1
        elif row['CCA-enzyme'] == 'CCA2':
            cca_status['CCA2'] = 1
        elif row['CCA-enzyme'] == 'CCA1+CCA2':
            cca_status['CCA12'] = 1
        
        metadata_dict[row['TipName']] = {
            'cca_enzyme': row['CCA-enzyme'],
            'superphyla': row['superphyla'],
            'cca_encoded': row['CCA-encoded'],
            'species_id': row['Species_ID'],
            'cca_status': cca_status
        }
    
    # Erstelle Heatmap-Daten in der Reihenfolge des Baums
    heatmap_data = []
    taxa_labels = []
    
    for tip_name in leaf_order:
        if tip_name in metadata_dict:
            meta = metadata_dict[tip_name]
            # Label für Y-Achse
            label = f"{meta['superphyla']}_{meta['species_id']}"
            taxa_labels.append(label)
            
            # Heatmap-Zeile: CCA1, CCA2, CCA12, CCA-encoded
            row = [
                meta['cca_status']['CCA1'],
                meta['cca_status']['CCA2'], 
                meta['cca_status']['CCA12'],
                meta['cca_encoded']
            ]
            heatmap_data.append(row)
    
    # Konvertiere zu numpy array
    heatmap_array = np.array(heatmap_data)
    
    # Spalten-Labels
    column_labels = ['CCA1', 'CCA2', 'CCA12', 'CCA-encoded']
    
    # Custom Colormap für die verschiedenen Spalten
    fig, ax = plt.subplots(figsize=(3, 12))
    
    # Farben für die einzelnen Zellen
    colors_map = {
        0: '#f0f0f0',  # Weiß/Grau für absent
        1: '#ccbf2c',  # Gelb für CCA1, CCA2, CCA12
        2: '#4CAF50'   # Grün für CCA-encoded
    }
    
    # Erstelle farbige Heatmap
    for i in range(len(taxa_labels)):
        for j in range(len(column_labels)):
            value = heatmap_array[i, j]
            
            # Bestimme Farbe
            if j == 3:  # CCA-encoded Spalte
                color = '#4CAF50' if value == 1 else '#f0f0f0'
            else:  # CCA1, CCA2, CCA12 Spalten
                color = '#ccbf2c' if value == 1 else '#f0f0f0'
            
            # Zeichne Rechteck - i direkt verwenden (nicht umkehren)
            rect = patches.Rectangle((j, len(taxa_labels)-i-1), 1, 1, 
                                   linewidth=0.5, edgecolor='black', 
                                   facecolor=color)
            ax.add_patch(rect)
    
    # Achsen konfigurieren
    ax.set_xlim(0, len(column_labels))
    ax.set_ylim(0, len(taxa_labels))
    
    # Labels setzen
    ax.set_xticks(np.arange(len(column_labels)) + 0.5)
    ax.set_xticklabels(column_labels, fontsize=10, fontweight='bold')
    ax.set_yticks(np.arange(len(taxa_labels)) + 0.5)
    # Taxa-Labels in der gleichen Reihenfolge wie der Baum (von oben nach unten)
    ax.set_yticklabels(list(reversed(taxa_labels)), fontsize=8)
    
    # Titel
    ax.set_title(f'CCA Heatmap {anticodon_title}', fontsize=12, fontweight='bold', pad=20)
    
    # X-Achse oben
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    # Entferne Achsenstriche
    ax.tick_params(length=0)
    
    # Grid entfernen
    ax.grid(False)
    
    # Layout anpassen
    plt.tight_layout()
    
    # Speichern
    output_path_heatmap = TREE_DIR / f"trna_heatmap_only_{anticodon_title}.png"
    plt.savefig(output_path_heatmap, dpi=400, bbox_inches='tight')
    print(f"[INFO] Heatmap gespeichert als: {output_path_heatmap}")
    
    plt.show()

def create_combined_visualization():
    #Erstellt sowohl Baum als auch Heatmap und kombiniert sie in einem Bild
    
    print("Erstelle phylogenetischen Baum...")
    leaf_order = create_tree_only()
    
    print("Erstelle Heatmap...")
    create_heatmap_only(leaf_order)
    
    print("Kombiniere Baum und Heatmap...")
    combine_tree_and_heatmap()
    
    print("Fertig! Kombinierte Visualisierung erstellt.")

def combine_tree_and_heatmap():
    #Kombiniert den gespeicherten Baum und die Heatmap in einem Bild
    
    import matplotlib.image as mpimg
    
    # Verwende manuellen Anticodon-Type für Titel und Dateinamen
    anticodon_title = ANTICODON_TYPE
    
    # Lade die gespeicherten Bilder
    tree_path = TREE_DIR / f"trna_tree_only_{anticodon_title}.png"
    heatmap_path = TREE_DIR / f"trna_heatmap_only_{anticodon_title}.png"
    
    if not tree_path.exists() or not heatmap_path.exists():
        print("Fehler: Tree oder Heatmap Datei nicht gefunden!")
        return
    
    # Lade Bilder
    tree_img = mpimg.imread(tree_path)
    heatmap_img = mpimg.imread(heatmap_path)
    
    # Parse Metadaten für die Anzahl der Taxa
    metadata_df = pd.read_csv(StringIO(METADATA), sep='\t')
    n_taxa = len(metadata_df)
    
    # Erstelle kombinierte Figure
    fig = plt.figure(figsize=(16, 12))
    
    # Berechne Subplot-Verhältnisse: Tree breiter als Heatmap
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.02)
    
    # Tree subplot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.imshow(tree_img)
    ax1.axis('off')
    ax1.set_title('Phylogenetic Tree', fontsize=14, fontweight='bold', pad=20)
    
    # Heatmap subplot  
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(heatmap_img)
    ax2.axis('off')
    ax2.set_title('CCA Heatmap', fontsize=14, fontweight='bold', pad=20)
    
    # Haupttitel
    fig.suptitle(f'tRNA Phylogenetic Tree with CCA Enzyme Heatmap ({anticodon_title})', 
                 fontsize=16, fontweight='bold', y=0.95)
    
    # Layout anpassen
    plt.tight_layout()
    
    # Speichere kombiniertes Bild
    output_path_combined = TREE_DIR / f"trna_tree_heatmap_combined_{anticodon_title}.png"
    plt.savefig(output_path_combined, dpi=400, bbox_inches='tight')
    print(f"[INFO] Kombiniertes Bild gespeichert als: {output_path_combined}")
    
    plt.show()
    plt.close()

def create_perfectly_aligned_visualization():
   
    print("Erstelle Visualisierung...")
    
    # Parse Metadaten
    metadata_df = load_metadata()
    
    # Verwende manuellen Anticodon-Type für Titel und Dateinamen
    anticodon_title = ANTICODON_TYPE
    print(f"[INFO] Verwende manuellen Anticodon-Type: {anticodon_title}")
    
    # Lade Baum
    NEWICK_TREE = load_newick_from_tree_file(TREE_FILE)
    NEWICK_TREE = NEWICK_TREE.replace("\\.", ".")
    tree = Tree(NEWICK_TREE, format=1)
    
    # DEBUG: Vergleiche Namen
    has_matches, name_mapping = debug_tree_names(tree, metadata_df)
    
    if not has_matches:
        print("\n[FEHLER] Keine übereinstimmenden Namen gefunden!")
        print("Das führt zu einem leeren/farblosen Baum.")
        print("Bitte überprüfen Sie:")
        print("1. Dass die Tree-Datei zu den Metadaten passt")
        print("2. Die Namensformate in beiden Dateien")
        return
    
    # Erstelle Metadaten-Dictionary mit Name-Mapping
    metadata_dict = {}
    for _, row in metadata_df.iterrows():
        # Bestimme CCA1, CCA2, CCA12 Status
        cca_status = {'CCA1': 0, 'CCA2': 0, 'CCA12': 0}
        if row['CCA-enzyme'] == 'CCA1':
            cca_status['CCA1'] = 1
        elif row['CCA-enzyme'] == 'CCA2':
            cca_status['CCA2'] = 1
        elif row['CCA-enzyme'] == 'CCA1+CCA2':
            cca_status['CCA12'] = 1
        
        # Erstelle Einträge für alle Tree-Namen, die auf diesen TipName mappen
        tip_name = row['TipName']
        for tree_name, mapped_name in name_mapping.items():
            if mapped_name == tip_name:
                metadata_dict[tree_name] = {
                    'cca_enzyme': row['CCA-enzyme'],
                    'superphyla': row['superphyla'],
                    'cca_encoded': row['CCA-encoded'],
                    'species_id': row['Species_ID'],
                    'cca_status': cca_status
                }
    
    # Farbschema für Superphyla
    colors = {
        'Asgard': '#e41a1c',      # Rot
        'DPANN': "#459dcc",        # Blau
        'Euryarchaeota': "#3ebc3a", # Grün
        'TACK': "#9533a3"          # Lila
    }
    
    # Extrahiere leaf order
    leaf_order = [leaf.name for leaf in tree.get_leaves()]
    print(f"[INFO] Baum hat {len(leaf_order)} Blätter")
    
    print(f"\n[INFO] Metadaten-Dictionary erstellt mit {len(metadata_dict)} Einträgen")
    
      
    # Extrahiere leaf order
    leaf_order = [leaf.name for leaf in tree.get_leaves()]
    print(f"[INFO] Baum hat {len(leaf_order)} Blätter")
    
    # Prüfe wie viele Matches wir haben
    matched_leaves = [name for name in leaf_order if name in metadata_dict]
    print(f"[INFO] {len(matched_leaves)} von {len(leaf_order)} Blättern haben Metadaten")
    
    # Erstelle Figure mit zwei Subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 12), 
                                   gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.05})
    
    # Baum
    # Für Ausrichtung: Erstelle den Tree erst normal und speichere als PNG
    
    # TreeStyle für rechteckige Darstellung
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "r"  # Rechteckig
    ts.scale = 20
    ts.branch_vertical_margin = 5
    ts.show_scale = False
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 1
    ts.guiding_lines_color = "#aaaaaa"
    
    # Funktion um Superphylum eines Blattes zu finden
    def get_leaf_superphylum(node):
        if node.is_leaf():
            if node.name in metadata_dict:
                return metadata_dict[node.name]['superphyla']
        else:
            superphyla = set()
            for leaf in node.get_leaves():
                if leaf.name in metadata_dict:
                    superphyla.add(metadata_dict[leaf.name]['superphyla'])
            if len(superphyla) == 1:
                return list(superphyla)[0]
        return None
    
    # Style für jeden Knoten
    for node in tree.traverse():
        nstyle = NodeStyle()
        superphylum = get_leaf_superphylum(node)
        
        if node.is_leaf() and node.name in metadata_dict:
            meta = metadata_dict[node.name]
            
            if superphylum:
                nstyle["vt_line_color"] = colors[superphylum]
                nstyle["hz_line_color"] = colors[superphylum]
                nstyle["vt_line_width"] = 3
                nstyle["hz_line_width"] = 3
            
            nstyle["size"] = 0
            
            label = f"{meta['superphyla']}_{meta['species_id']}"
            label_color = colors.get(meta['superphyla'], "black")
            name_face = TextFace(label, fsize=10, fgcolor=label_color, bold=True)
            node.add_face(name_face, column=0, position="aligned")
            
        else:
            nstyle["size"] = 6
            nstyle["fgcolor"] = "black"
            
            if superphylum:
                nstyle["vt_line_color"] = colors[superphylum]
                nstyle["hz_line_color"] = colors[superphylum]
                nstyle["vt_line_width"] = 3
                nstyle["hz_line_width"] = 3
        
        node.set_style(nstyle)
    
    # Speichere Tree temporär und lade in subplot
    temp_tree_path = TREE_DIR / f"temp_tree_{anticodon_title}.png"
    tree.render(str(temp_tree_path), tree_style=ts, dpi=400, w=1200, h=1400)
    
    # Lade Tree-Bild in ersten subplot
    import matplotlib.image as mpimg
    tree_img = mpimg.imread(temp_tree_path)
    ax1.imshow(tree_img)
    ax1.axis('off')
    ax1.set_title(f'Phylogenetic Tree {anticodon_title}', fontsize=14, fontweight='bold')
    
    
    # HEATMAP TEIL
  
    
    # Erstelle Heatmap-Daten in der Reihenfolge des Baums
    heatmap_data = []
    taxa_labels = []
    
    for tip_name in leaf_order:
        if tip_name in metadata_dict:
            meta = metadata_dict[tip_name]
            label = f"{meta['superphyla']}_{meta['species_id']}"
            taxa_labels.append(label)
            
            row = [
                meta['cca_status']['CCA1'],
                meta['cca_status']['CCA2'], 
                meta['cca_status']['CCA12'],
                meta['cca_encoded']
            ]
            heatmap_data.append(row)
    
    heatmap_array = np.array(heatmap_data)
    column_labels = ['CCA1', 'CCA2', 'CCA12', 'CCA-encoded']
    
    # Erstelle farbige Heatmap im zweiten subplot
    for i in range(len(taxa_labels)):
        for j in range(len(column_labels)):
            value = heatmap_array[i, j]
            
            if j == 3:  # CCA-encoded Spalte
                color = '#4CAF50' if value == 1 else '#f0f0f0'
            else:  # CCA1, CCA2, CCA12 Spalten
                color = '#ccbf2c' if value == 1 else '#f0f0f0'
            
            rect = patches.Rectangle((j, len(taxa_labels)-i-1), 1, 1, 
                                   linewidth=0.5, edgecolor='black', 
                                   facecolor=color)
            ax2.add_patch(rect)
    
    # Konfiguriere Heatmap-Achsen
    ax2.set_xlim(0, len(column_labels))
    ax2.set_ylim(0, len(taxa_labels))
    ax2.set_xticks(np.arange(len(column_labels)) + 0.5)
    ax2.set_xticklabels(column_labels, fontsize=10, fontweight='bold')
    ax2.set_yticks([])  # Keine Y-Labels, da Tree schon Labels hat
    ax2.set_title(f'CCA Heatmap {anticodon_title}', fontsize=14, fontweight='bold')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.tick_params(length=0)
    ax2.grid(False)
    
    # Haupttitel
    fig.suptitle(f'tRNA Phylogenetic Tree with CCA Enzyme Heatmap ({anticodon_title})', 
                 fontsize=16, fontweight='bold')
    
    # Layout anpassen
    plt.tight_layout()
    
    # Speichere kombiniertes Bild
    output_path_perfect = TREE_DIR / f"trna_tree_heatmap_perfect_aligned_{anticodon_title}.png"
    plt.savefig(output_path_perfect, dpi=400, bbox_inches='tight')
    print(f"[INFO] Bild gespeichert als: {output_path_perfect}")
    
    # Entferne temporäre Datei
    if temp_tree_path.exists():
        temp_tree_path.unlink()
    
    plt.show()
    plt.close()

if __name__ == "__main__":
    print("Erstelle kombinierte Baum- und Heatmap-Visualisierung...")
    print(f"[INFO] Anticodon-Type für Plot-Namen: {ANTICODON_TYPE}")
    print("       (Zum Ändern: ANTICODON_TYPE Variable am Anfang des Skripts anpassen)")
    
    # Prüfe ob Metadaten-Datei existiert
    if not METADATA_FILE.exists():
        print(f"FEHLER: Metadaten-Datei nicht gefunden: {METADATA_FILE}")
        print("Bitte stellen Sie sicher, dass die Datei 'TABLE_metadata_Isoleucine.txt' im OUTPUT_DIR liegt.")
        sys.exit(1)
    
    try:
        create_perfectly_aligned_visualization()
    except Exception as e:
        print(f"FEHLER bei der Visualisierung: {e}")
        print("Überprüfen Sie die Metadaten-Datei und Tree-Datei.")
    
 