#!/usr/bin/env python3
#tRNA Best Candidates Filter
#Anwendung der spezifizierten Filterkriterien für beste tRNA-Kandidaten mit Superphylum-Aufschlüsselung

import pandas as pd
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
from typing import Dict, List, Tuple, Optional

# config einbinden
sys.path.append(str(Path(__file__).resolve().parent.parent))
from config import TRNASCAN_OUTPUT_DIR, PLOTS_DIR, OUTPUT_DIR

class tRNALoopAnalyzer:
    #Analysiert tRNA-Sekundärstrukturen und extrahiert Loop-Informationen
    
    
    def __init__(self):
        # Konfigurierbare Thresholds
        self.thresholds = {
            'min_stems': 4,
            'max_stems': 5,
            'max_size': 110,
            'min_anticodon_loop': 7,
            'max_anticodon_loop': 14,
            'max_d_loop': 20,
            'max_t_loop': 10,
            'max_v_loop': 6,
            'min_isotype_score': 60.0
        }
    
    def analyze_loops_from_structure(self, dot_bracket: str) -> Dict[str, int]:
        #Extrahiert Loop-Größen aus Dot-Bracket-Notation
        #Vereinfachte Heuristik zur Loop-Identifikation
        
        if not dot_bracket or pd.isna(dot_bracket):
            return {
                'anticodon_loop': 0,
                'd_loop': 0,
                't_loop': 0,
                'v_loop': 0,
                'total_loops': 0
            }
        
        structure = str(dot_bracket).replace('>', ')').replace('<', '(')
        
        # Finde alle Loop-Regionen (aufeinanderfolgende Punkte zwischen Stems)
        loops = []
        current_loop = ""
        in_loop = False
        
        for char in structure:
            if char == '.':
                if not in_loop:
                    in_loop = True
                    current_loop = "."
                else:
                    current_loop += "."
            else:  # '(' oder ')'
                if in_loop:
                    loops.append(len(current_loop))
                    current_loop = ""
                    in_loop = False
        
        # Letzter Loop falls String mit Loop endet
        if in_loop and current_loop:
            loops.append(len(current_loop))
        
        # Sortiere Loops nach Größe (größter zuerst)
        loops.sort(reverse=True)
        
        # Heuristik zur Loop-Zuordnung basierend auf typischen tRNA-Mustern:
        # - Größter Loop ist meist Anticodon-Loop (7-14 nt)
        # - Zweitgrößter ist oft D-Loop
        # - Weitere Loops sind T-Loop und V-Loop
        
        result = {
            'anticodon_loop': loops[0] if len(loops) > 0 else 0,
            'd_loop': loops[1] if len(loops) > 1 else 0,
            't_loop': loops[2] if len(loops) > 2 else 0,
            'v_loop': loops[3] if len(loops) > 3 else 0,
            'total_loops': len(loops)
        }
        
        # Validierung: Anticodon-Loop sollte typischerweise der größte sein
        # Falls nicht, versuche bessere Zuordnung
        if len(loops) >= 2 and result['anticodon_loop'] < 7:
            # Suche nach Loop im typischen Anticodon-Bereich
            anticodon_candidates = [l for l in loops if 7 <= l <= 14]
            if anticodon_candidates:
                result['anticodon_loop'] = max(anticodon_candidates)
                # Entferne diesen aus der Liste für andere Zuordnungen
                remaining_loops = [l for l in loops if l != result['anticodon_loop']]
                result['d_loop'] = remaining_loops[0] if len(remaining_loops) > 0 else 0
                result['t_loop'] = remaining_loops[1] if len(remaining_loops) > 1 else 0
                result['v_loop'] = remaining_loops[2] if len(remaining_loops) > 2 else 0
        
        return result
    
    def check_cca_terminus(self, sequence: str) -> bool:
        #Überprüft das Vorhandensein des 3' CCA-Terminus
        
        if not sequence or pd.isna(sequence):
            return False
        
        seq = str(sequence).upper().replace('U', 'T')
        return seq.endswith('CCA')
    
    def apply_best_candidate_filters(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, int], Dict[str, Dict[str, int]]]:
        #Wendet alle Filterkriterien für die besten tRNA-Kandidaten an 
        
        print("[INFO] Starte Filterung für beste tRNA-Kandidaten...")
        
        initial_count = len(df)
        filter_stats = {'initial': initial_count}
        
        # Initialisiere Superphylum-spezifische Statistiken
        superphyla = df['Superphylum'].unique() if 'Superphylum' in df.columns else ['Unknown']
        superphylum_stats = {}
        for sp in superphyla:
            superphylum_stats[sp] = {}
        
        # Kopie für Bearbeitung
        filtered_df = df.copy()
        
        print(f"[INFO] Ausgangsdaten: {initial_count} tRNAs")
        print(f"[INFO] Gefundene Superphyla: {list(superphyla)}")
        
        # Speichere initiale Verteilung
        if 'Superphylum' in filtered_df.columns:
            initial_dist = filtered_df['Superphylum'].value_counts()
            for sp in superphyla:
                superphylum_stats[sp]['initial'] = initial_dist.get(sp, 0)
                print(f"  {sp}: {superphylum_stats[sp]['initial']} tRNAs")
        
        # Filter 1: 4 oder 5 Stems
        print("\n[FILTER 1] Stem-Anzahl: 4 oder 5 Stems")
        filtered_df = filtered_df[filtered_df['StemCount'].isin([4, 5])]
        filter_stats['stem_filter'] = len(filtered_df)
        
        if 'Superphylum' in filtered_df.columns:
            current_dist = filtered_df['Superphylum'].value_counts()
            for sp in superphyla:
                superphylum_stats[sp]['stem_filter'] = current_dist.get(sp, 0)
                print(f"  {sp}: {superphylum_stats[sp]['stem_filter']} tRNAs")
        
        print(f"  Gesamt nach Stem-Filter: {len(filtered_df)} tRNAs ({len(filtered_df)/initial_count*100:.1f}%)")
        
        # Filter 2: Länge
        print("\n[FILTER 2] Größe unter ['max_size']")
        if 'MatureLen' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["MatureLen"] < self.thresholds['max_size']]
            filter_stats['size_filter'] = len(filtered_df)
            
            if 'Superphylum' in filtered_df.columns:
                current_dist = filtered_df['Superphylum'].value_counts()
                for sp in superphyla:
                    superphylum_stats[sp]['size_filter'] = current_dist.get(sp, 0)
                    print(f"  {sp}: {superphylum_stats[sp]['size_filter']} tRNAs")
                    
            print(f"  Gesamt nach Größen-Filter: {len(filtered_df)} tRNAs ({len(filtered_df)/initial_count*100:.1f}%)")
        else:
            print("  [WARNUNG] 'MatureLen' nicht vorhanden, Filter übersprungen")
            filter_stats['size_filter'] = len(filtered_df)
            if 'Superphylum' in filtered_df.columns:
                for sp in superphyla:
                    superphylum_stats[sp]['size_filter'] = superphylum_stats[sp]['stem_filter']
                
        # Filter 3: Keine Pseudogene
        print("\n[FILTER 3] Eliminiere Pseudogene")
        if 'Type' in filtered_df.columns:
            before_pseudo = len(filtered_df)
            filtered_df = filtered_df[filtered_df['Type'] != 'Pseudo']
            filter_stats['pseudo_filter'] = len(filtered_df)
            eliminated = before_pseudo - len(filtered_df)
            
            if 'Superphylum' in filtered_df.columns:
                current_dist = filtered_df['Superphylum'].value_counts()
                for sp in superphyla:
                    superphylum_stats[sp]['pseudo_filter'] = current_dist.get(sp, 0)
                    print(f"  {sp}: {superphylum_stats[sp]['pseudo_filter']} tRNAs")
            
            print(f"  {eliminated} Pseudogene eliminiert")
            print(f"  Gesamt nach Pseudogen-Filter: {len(filtered_df)} tRNAs ({len(filtered_df)/initial_count*100:.1f}%)")
        else:
            print("  [WARNUNG] Keine Type-Spalte gefunden, Filter übersprungen")
            filter_stats['pseudo_filter'] = len(filtered_df)
            if 'Superphylum' in filtered_df.columns:
                for sp in superphyla:
                    superphylum_stats[sp]['pseudo_filter'] = superphylum_stats[sp]['size_filter']
        
        # Filter 4: Isotype Score
        print("\n[FILTER 4] Isotype Score >= ['min_isotype_score']")
        if 'IsotypeScore' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['IsotypeScore'] >= self.thresholds['min_isotype_score']]
            filter_stats['isotype_filter'] = len(filtered_df)
            
            if 'Superphylum' in filtered_df.columns:
                current_dist = filtered_df['Superphylum'].value_counts()
                for sp in superphyla:
                    superphylum_stats[sp]['isotype_filter'] = current_dist.get(sp, 0)
                    print(f"  {sp}: {superphylum_stats[sp]['isotype_filter']} tRNAs")
                    
            print(f"  Gesamt nach Isotype-Score-Filter: {len(filtered_df)} tRNAs ({len(filtered_df)/initial_count*100:.1f}%)")
        else:
            print("  [WARNUNG] Keine IsotypeScore-Spalte gefunden, Filter übersprungen")
            filter_stats['isotype_filter'] = len(filtered_df)
            if 'Superphylum' in filtered_df.columns:
                for sp in superphyla:
                    superphylum_stats[sp]['isotype_filter'] = superphylum_stats[sp]['pseudo_filter']
        
        # Filter 5-7: Loop-Analysen
        print("\n[FILTER 5-7] Loop-Größen-Analyse...")
        if 'StrParse' in filtered_df.columns or 'Str' in filtered_df.columns:
            struct_col = 'StrParse' if 'StrParse' in filtered_df.columns else 'Str'
            
            # Analysiere Loops für alle verbleibenden tRNAs
            loop_analyses = []
            for idx, row in filtered_df.iterrows():
                loop_info = self.analyze_loops_from_structure(row[struct_col])
                loop_analyses.append(loop_info)
            
            # Füge Loop-Informationen als neue Spalten hinzu
            for key in ['anticodon_loop', 'd_loop', 't_loop', 'v_loop', 'total_loops']:
                filtered_df[key] = [analysis[key] for analysis in loop_analyses]
            
            before_loop = len(filtered_df)
            
            # Anticodon-Loop Filter
            filtered_df = filtered_df[
                (filtered_df['anticodon_loop'] >= self.thresholds['min_anticodon_loop']) &
                (filtered_df['anticodon_loop'] <= self.thresholds['max_anticodon_loop'])
            ]
            
            # D-Loop Filter
            filtered_df = filtered_df[filtered_df['d_loop'] <= self.thresholds['max_d_loop']]
            
            # T-Loop Filter  
            filtered_df = filtered_df[filtered_df['t_loop'] <= self.thresholds['max_t_loop']]
            
            # V-Loop Filter
            filtered_df = filtered_df[filtered_df['v_loop'] <= self.thresholds['max_v_loop']]
            
            filter_stats['loop_filter'] = len(filtered_df)
            eliminated = before_loop - len(filtered_df)
            
            if 'Superphylum' in filtered_df.columns:
                current_dist = filtered_df['Superphylum'].value_counts()
                for sp in superphyla:
                    superphylum_stats[sp]['loop_filter'] = current_dist.get(sp, 0)
                    print(f"  {sp}: {superphylum_stats[sp]['loop_filter']} tRNAs")
            
            print(f"  {eliminated} tRNAs durch Loop-Filter eliminiert")
            print(f"  Gesamt nach Loop-Filtern: {len(filtered_df)} tRNAs ({len(filtered_df)/initial_count*100:.1f}%)")
        else:
            print("  [WARNUNG] Keine Struktur-Spalte gefunden, Loop-Filter übersprungen")
            filter_stats['loop_filter'] = len(filtered_df)
            if 'Superphylum' in filtered_df.columns:
                for sp in superphyla:
                    superphylum_stats[sp]['loop_filter'] = superphylum_stats[sp]['isotype_filter']
        
        # Finale Statistiken
        filter_stats['final'] = len(filtered_df)
        if 'Superphylum' in filtered_df.columns:
            final_dist = filtered_df['Superphylum'].value_counts()
            for sp in superphyla:
                superphylum_stats[sp]['final'] = final_dist.get(sp, 0)
        
        print(f"\n[ZUSAMMENFASSUNG] Filterung abgeschlossen:")
        print(f"  Ausgangsdaten: {filter_stats['initial']}")
        print(f"  Nach allen Filtern: {filter_stats['final']} ({filter_stats['final']/filter_stats['initial']*100:.1f}%)")
        
        if 'Superphylum' in filtered_df.columns:
            print(f"\n[FINALE VERTEILUNG]:")
            for sp in superphyla:
                print(f"  {sp}: {superphylum_stats[sp]['final']} tRNAs")
        
        return filtered_df, filter_stats, superphylum_stats

def create_filter_summary_plot(filter_stats: Dict[str, int], output_dir: Path):
    #Erstellt ein ganz okayisches Diagramm der Filterungsschritte
    
    steps = ['initial', 'stem_filter', 'size_filter', 'pseudo_filter', 'isotype_filter', 'loop_filter', 'final']
    step_labels = ['Start', 'only 4-5 Stems', 'Size <110 nts', 'No Pseudo', 'IsotypeScore ≥60', 'Loop-Filter', 'Final']
    counts = [filter_stats.get(step, 0) for step in steps]
    
    plt.figure(figsize=(12, 6))
    plt.plot(step_labels, counts, marker='o', linewidth=2, markersize=8)
    plt.title('tRNA-Filter Process: Best Candidates')
    plt.xlabel('Filter Step')
    plt.ylabel('Number of tRNAs')
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)
    
    # Annotiere Werte
    for i, count in enumerate(counts):
        plt.annotate(f'{count}', (i, count), textcoords="offset points", xytext=(0,10), ha='center')
    
    plt.tight_layout()
    
    plot_path = PLOTS_DIR / "best_candidates_filtering_process.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"[INFO] Filterungsdiagramm gespeichert: {plot_path}")


def create_filter_summary_plot_by_superphylum(filter_stats: Dict[str, int], superphylum_stats: Dict[str, Dict[str, int]], output_dir: Path):
    #Erstellt ein Diagramm der Filterungsschritte aufgeteilt nach Superphylum
    
    steps = ['initial', 'stem_filter', 'size_filter', 'pseudo_filter', 'isotype_filter', 'loop_filter', 'final']
    step_labels = ['Start', '4-5 Stems', 'Size <110', 'No Pseudo', 'IsotypeScore ≥60', 'Loop-Filter', 'Final']
    
    # Definiere Farben für Superphyla
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    superphyla = list(superphylum_stats.keys())
    
    plt.figure(figsize=(14, 8))
    
    # Sammle alle Counts für y-Achsen-Skalierung
    all_counts = []
    
    # Plotte für jedes Superphylum eine separate Linie
    for i, sp in enumerate(superphyla):
        counts = [superphylum_stats[sp].get(step, 0) for step in steps]
        all_counts.extend(counts)
        plt.plot(step_labels, counts, marker='o', linewidth=3, markersize=9, 
                color=colors[i % len(colors)], label=sp)
        
        # Annotiere Endwerte
        final_count = counts[-1]
        if final_count > 0:
            plt.annotate(f'{final_count}', (len(step_labels)-1, final_count), 
                        textcoords="offset points", xytext=(8,0), ha='left',
                        color=colors[i % len(colors)], fontweight='bold', fontsize=11)
    
    plt.title('tRNA Filter Process: Best Candidates by Superphylum', fontsize=16, fontweight='bold')
    plt.xlabel('Filter Step', fontsize=13)
    plt.ylabel('Number of tRNAs', fontsize=13)
    plt.xticks(rotation=45, fontsize=11)
    plt.yticks(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    
    # Setze y-Achse mit etwas Spielraum basierend auf den tatsächlichen Superphylum-Daten
    max_count = max(all_counts) if all_counts else 100
    plt.ylim(0, max_count * 1.1)
    plt.tight_layout()

    plot_path = PLOTS_DIR / "best_candidates_filtering_by_superphylum.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    print(f"[INFO] Superphylum-spezifisches Filterungsdiagramm gespeichert: {plot_path}")
    plt.show()

def create_superphylum_summary_table(superphylum_stats: Dict[str, Dict[str, int]], output_dir: Path):
    #Erstellt eine Tabelle mit der Zusammenfassung der Filterung pro Superphylum
    
    steps = ['initial', 'stem_filter', 'size_filter', 'pseudo_filter', 'isotype_filter', 'loop_filter', 'final']
    step_labels = ['Start', '4-5 Stems', 'Size <110', 'No Pseudo', 'IsotypeScore ≥60', 'Loop-Filter', 'Final']
    
    # Erstelle DataFrame für die Tabelle
    summary_data = []
    for sp in superphylum_stats.keys():
        row = {'Superphylum': sp}
        for step, label in zip(steps, step_labels):
            row[label] = superphylum_stats[sp].get(step, 0)
        
        # Berechne Prozentsätze
        initial = superphylum_stats[sp].get('initial', 1)
        final = superphylum_stats[sp].get('final', 0)
        row['Retention_%'] = f"{final/initial*100:.1f}%" if initial > 0 else "0%"
        
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    
    # Speichere als CSV
    table_path = OUTPUT_DIR / "superphylum_filtering_summary.csv"
    summary_df.to_csv(table_path, index=False)
    print(f"[INFO] Superphylum-Zusammenfassung gespeichert: {table_path}")
    
    # Zeige Tabelle an
    print(f"\n[SUPERPHYLUM SUMMARY TABLE]")
    print(summary_df.to_string(index=False))
    
    return summary_df


def analyze_best_candidates_quality(df: pd.DataFrame, output_dir: Path):
    #Analysiert die Qualität der finalen besten Kandidaten
    
    print("\n[QUALITÄTSANALYSE] Beste tRNA-Kandidaten")
    
    # Grundstatistiken
    print(f"Finale Anzahl bester Kandidaten: {len(df)}")
    
    if len(df) == 0:
        print("Keine Kandidaten nach Filterung übrig!")
        return
    
    # Verteilung nach Superphylum
    if 'Superphylum' in df.columns:
        print("\nVerteilung nach Superphylum:")
        superphylum_counts = df['Superphylum'].value_counts()
        print(superphylum_counts)
    
    # Verteilung nach tRNA-Typ
    if 'Type' in df.columns:
        print(f"\nVerteilung nach tRNA-Typ:")
        type_counts = df['Type'].value_counts()
        print(type_counts.head(10))
    
    # Score-Statistiken
    score_columns = ['HMMScore', 'Str2Score', 'IsotypeScore']
    print(f"\nScore-Statistiken:")
    for col in score_columns:
        if col in df.columns:
            print(f"  {col}: Mittel={df[col].mean():.1f}, Median={df[col].median():.1f}, Min={df[col].min():.1f}, Max={df[col].max():.1f}")
    
    # Loop-Statistiken
    loop_columns = ['anticodon_loop', 'd_loop', 't_loop', 'v_loop']
    print(f"\nLoop-Größen-Statistiken:")
    for col in loop_columns:
        if col in df.columns:
            print(f"  {col}: Mittel={df[col].mean():.1f}, Median={df[col].median():.1f}, Bereich={df[col].min()}-{df[col].max()}")
    
    # Erstelle Qualitätsplot
    if 'Superphylum' in df.columns:
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Plot 1: Score-Verteilungen nach Superphylum
        if 'IsotypeScore' in df.columns:
            df.boxplot(column='IsotypeScore', by='Superphylum', ax=axes[0,0])
            axes[0,0].set_title('IsotypeScore by Superphylum')
            axes[0,0].set_xlabel('Superphylum')
        
        # Plot 2: Loop-Größen-Verteilungen
        if 'anticodon_loop' in df.columns:
            df.boxplot(column='anticodon_loop', by='Superphylum', ax=axes[0,1])
            axes[0,1].set_title('Anticodon Loop Size by Superphylum')
            axes[0,1].set_xlabel('Superphylum')
        
        # Plot 3: Stem Count Verteilung
        if 'StemCount' in df.columns:
            df.groupby('Superphylum')['StemCount'].value_counts().unstack().plot(kind='bar', ax=axes[1,0])
            axes[1,0].set_title('Stem Count Distribution by Superphylum')
            axes[1,0].set_xlabel('Superphylum')
            axes[1,0].legend(title='Stem Count')
        
        # Plot 4: tRNA Typ Verteilung
        if 'Type' in df.columns:
            type_counts_by_sp = df.groupby(['Superphylum', 'Type']).size().unstack(fill_value=0)
            type_counts_by_sp.plot(kind='bar', ax=axes[1,1], stacked=True)
            axes[1,1].set_title('tRNA Type Distribution by Superphylum')
            axes[1,1].set_xlabel('Superphylum')
            axes[1,1].legend(title='tRNA Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        plot_path = PLOTS_DIR / "best_candidates_quality_analysis.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"[INFO] Qualitätsdiagramm gespeichert: {plot_path}")
    

def main():
    #Hauptfunktion: Lädt annotierte Daten und wendet Best-Candidate-Filter an
    
    # Dateipfade
    input_path = TRNASCAN_OUTPUT_DIR / "ALL_combined_mainstruct_annotated.txt"
    best_candidates_path = TRNASCAN_OUTPUT_DIR / "ALL_best_candidates.txt"
    rejected_candidates_path = TRNASCAN_OUTPUT_DIR / "ALL_rejected_by_filters.txt"
    
    print(f"[INFO] Lade annotierte Daten von: {input_path}")
    
    if not input_path.exists():
        print(f"[FEHLER] Eingabedatei nicht gefunden: {input_path}")
        print("Führen Sie zuerst die Strukturanalyse aus!")
        return
    
    # Lade Daten
    df = pd.read_csv(input_path, sep="\t")
    print(f"[INFO] {len(df)} tRNAs geladen")
    print(f"[INFO] Verfügbare Spalten: {df.columns.tolist()}")
    
    # Initialisiere Analyzer
    analyzer = tRNALoopAnalyzer()
    
    # Wende Filter an (mit Superphylum-Tracking)
    best_candidates, filter_stats, superphylum_stats = analyzer.apply_best_candidate_filters(df)
    
    # Speichere beste Kandidaten
    best_candidates.to_csv(best_candidates_path, sep="\t", index=False)
    print(f"\n[EXPORT] Beste Kandidaten gespeichert: {best_candidates_path}")
    
    # Speichere abgelehnte Kandidaten (für weitere Analyse)
    if len(best_candidates) < len(df):
        # Finde IDs der besten Kandidaten
        if 'Sequence' in df.columns and 'tRNA#_x' in df.columns:
            best_ids = set(zip(best_candidates['Sequence'], best_candidates['tRNA#_x']))
            rejected = df[~df.apply(lambda row: (row['Sequence'], row['tRNA#_x']) in best_ids, axis=1)]
        else:
            # Fallback: verwende Index
            rejected = df.drop(best_candidates.index)
        
        rejected.to_csv(rejected_candidates_path, sep="\t", index=False)
        print(f"[EXPORT] Abgelehnte Kandidaten gespeichert: {rejected_candidates_path}")
    
    # Erstelle Visualisierungen
    create_filter_summary_plot(filter_stats, TRNASCAN_OUTPUT_DIR)
    create_filter_summary_plot_by_superphylum(filter_stats, superphylum_stats, TRNASCAN_OUTPUT_DIR)
    create_superphylum_summary_table(superphylum_stats, TRNASCAN_OUTPUT_DIR)
    
    # Qualitätsanalyse
    analyze_best_candidates_quality(best_candidates, TRNASCAN_OUTPUT_DIR)
    
    print(f"\n[ABGESCHLOSSEN] Best Candidates Filtering")
    print(f"  Eingabe: {len(df)} tRNAs")
    print(f"  Beste Kandidaten: {len(best_candidates)} tRNAs ({len(best_candidates)/len(df)*100:.1f}%)")
    print(f"  Abgelehnt: {len(df) - len(best_candidates)} tRNAs")


if __name__ == "__main__":
    main()