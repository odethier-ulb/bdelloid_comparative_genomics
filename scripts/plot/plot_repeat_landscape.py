import json
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
from matplotlib.patches import Patch


"""
Plot repeat landscape from RepeatMasker divsum file, 
skip Unknown, Satellite, Segmental, Structural_RNA, Simple_repeat, rRNA and tRNA.
"""

#----- From RepeatMasker createRepeatLandscape.pl -----

NAME_MAP = { 
    "DNA/Chompy": "DNA",
    "DNA/CMC-Chapaev": "DNA/CMC",
    "DNA/CMC-Chapaev-3": "DNA/CMC",
    "DNA/CMC-EnSpm": "DNA/CMC",
    "DNA/CMC-Transib": "DNA/CMC",
    "DNA/En-Spm": "DNA/CMC",
    "DNA/Crypton-3": "DNA/Crypton",
    "DNA/PIF-Harbinger": "DNA/Harbinger",
    "DNA/PIF-ISL2EU": "DNA/Harbinger",
    "DNA/Tourist": "DNA/Harbinger",
    "DNA/AcHobo": "DNA/hAT",
    "DNA/Charlie": "DNA/hAT",
    "DNA/Chompy1": "DNA/hAT",
    "DNA/MER1_type": "DNA/hAT",
    "DNA/Tip100": "DNA/hAT",
    "DNA/hAT-Ac": "DNA/hAT",
    "DNA/hAT-Blackjack": "DNA/hAT",
    "DNA/hAT-Charlie": "DNA/hAT",
    "DNA/hAT-Tag1": "DNA/hAT",
    "DNA/hAT-Tip100": "DNA/hAT",
    "DNA/hAT-hATw": "DNA/hAT",
    "DNA/hAT-hobo": "DNA/hAT",
    "DNA/hAT_Tol2": "DNA/hAT",
    "DNA/Kolobok-IS4EU": "DNA/Kolobok",
    "DNA/Kolobok-T2": "DNA/Kolobok",
    "DNA/T2": "DNA/Kolobok",
    "DNA/MULE-MuDR": "DNA/MULE",
    "DNA/MULE-NOF": "DNA/MULE",
    "DNA/MuDR": "DNA/MULE",
    "DNA/piggyBac": "DNA/PiggyBac",
    "DNA/MER2_type": "DNA/TcMar",
    "DNA/Mariner": "DNA/TcMar",
    "DNA/Pogo": "DNA/TcMar",
    "DNA/Stowaway": "DNA/TcMar",
    "DNA/Tc1": "DNA/TcMar",
    "DNA/Tc2": "DNA/TcMar",
    "DNA/Tc4": "DNA/TcMar",
    "DNA/TcMar-Fot1": "DNA/TcMar",
    "DNA/TcMar-ISRm11": "DNA/TcMar",
    "DNA/TcMar-Mariner": "DNA/TcMar",
    "DNA/TcMar-Pogo": "DNA/TcMar",
    "DNA/TcMar-Tc1": "DNA/TcMar",
    "DNA/TcMar-Tc2": "DNA/TcMar",
    "DNA/TcMar-Tigger": "DNA/TcMar",
    "DNA/Tigger": "DNA/TcMar",    
    "DNA/Helitron": "RC/Helitron",
    "LTR/DIRS1": "LTR/DIRS",
    "LTR/ERV-Foamy": "LTR/ERVL",
    "LTR/ERV-Lenti": "LTR/ERV",
    "LTR/ERVL-MaLR": "LTR/ERVL",
    "LTR/Gypsy-Troyka": "LTR/Gypsy",
    "LTR/MaLR": "LTR/ERVL",
    "LINE/CR1-Zenon": "LINE/CR1",
    "LINE/I": "LINE/Jockey-I",
    "LINE/Jockey": "LINE/Jockey-I",
    "LINE/L1-Tx1": "LINE/L1",
    "LINE/R2-Hero": "LINE/R2",
    "LINE/RTE-BovB": "LINE/RTE",
    "LINE/RTE-RTE": "LINE/RTE",
    "LINE/RTE-X": "LINE/RTE",
    "LINE/telomeric": "LINE/Jockey-I",
    "LINE/Penlope": "PLE",
    "PLE/Hydra": "PLE",
    "PLE/Naiad": "PLE",
    "PLE/Athena": "PLE",
    "PLE/Poseidon": "PLE",
    "PLE/Nematis": "PLE",
    "PLE/Neptune": "PLE",
    "PLE/Coprina": "PLE",
    "SINE/B2": "SINE/tRNA",
    "SINE/B4": "SINE/tRNA-Alu",
    "SINE/BovA": "SINE/tRNA-RTE",
    "SINE/C": "SINE/tRNA",
    "SINE/Core": "SINE",
    "SINE/ID": "SINE/tRNA",
    "SINE/Lys": "SINE/tRNA",
    "SINE/MERMAID": "SINE/tRNA-V",
    "SINE/RTE-BovB": "SINE/RTE",
    "SINE/tRNA-Glu": "SINE/tRNA",
    "SINE/tRNA-Lys": "SINE/tRNA",
    "SINE/V": "SINE/tRNA-V",
    "Unknown/Y-chromosome": "Unknown",

    "DNA/CMC-3": "DNA/CMC",
    "DNA/CMC-Mirage": "DNA/CMC",
    "DNA/Dada": "DNA/Dada",
    "DNA/Kolobok-Hydra": "DNA/Kolobok",
    "DNA/MULE-F": "DNA/MULE",
    "DNA/TcMar-Stowaway": "DNA/TcMar",
    "DNA/TcMar-Tc4": "DNA/TcMar",
    "DNA/TcMar-m44": "DNA/TcMar",
    "DNA/hAT-Pegasus": "DNA/hAT",
    "DNA/hAT-Tol2": "DNA/hAT",
    "DNA/hAT-hAT1": "DNA/hAT",
    "DNA/hAT-hAT5": "DNA/hAT",
    "DNA/hAT-hAT6": "DNA/hAT",
    "DNA/hAT-hAT19": "DNA/hAT",
    "LINE/CRE": "LINE/CRE",
    "LINE/Jockey-I-I": "LINE/Jockey-I",
    "LINE/I-Jockey": "LINE/Jockey-I",
    "LINE/R2-NeSL": "LINE/R2",
    "LTR/Copia(Xen1)": "LTR/Copia",
    "LTR/ERV4": "LTR/ERV",

    # Mistake should have been LTR/Gypsy to begin with
    "LTR/Ginger": "LTR/Gypsy",
    "Retroposon": "Other",
    "Other/Composite": "Other",
    "SINE/5S-Deu-L2": "SINE/5S",
    "SINE/5S-Sauria-RTE": "SINE/5S",
    "SINE/L2": "SINE",
    "SINE/Mermaid": "SINE/tRNA",
    "SINE/R2": "SINE",
    "SINE/U": "SINE/U",
    "SINE/Core-RTE": "SINE/RTE",
    "SINE/tRNA-C": "SINE/tRNA",
    "SINE/tRNA-Core-RTE": "SINE/tRNA",
    "SINE/tRNA-Core": "SINE/tRNA",
    "SINE/tRNA-Deu-CR1": "SINE/Deu",
    "SINE/tRNA-Deu-L2": "SINE/Deu",
    "SINE/tRNA-Deu": "SINE/Deu",
    "SINE/tRNA-Jockey": "SINE/tRNA",
    "SINE/tRNA-L2": "SINE/tRNA",
    "SINE/tRNA-Rex": "SINE/tRNA",
    "SINE/tRNA-Sauria-L2": "SINE/tRNA",
    "SINE/tRNA-Sauria-RTE": "SINE/tRNA",
    "SINE/tRNA-Sauria": "SINE/tRNA",
    "SINE/Sauria": "SINE/tRNA",
    "SINE/tRNA-V-Core-L2": "SINE/tRNA",

    # No longer in the database
    "SINE/tRNAore-RTE": "SINE/tRNA",
    "SINE/tRNAore": "SINE/tRNA",
    "tRNA": "Structural_RNA",
    "scRNA": "Structural_RNA",
    "RNA": "Structural_RNA",
    "snRNA": "Structural_RNA",
    "rRNA": "Structural_RNA",
    "Satellite/centromeric": "Satellite",
    "Satellite/telomeric": "Satellite",
    "Satellite/acromeric": "Satellite",
    
    # Custom from this script
    "LINE/Penelope": "PLE", 
    "LINE/Penelope-Athena": "PLE",
    "DNA/RC": "RC/Helitron"
}

GRAPH_LABELS = {
    "Unknown": "#999999",
    "Other": "#4D4D4D",
    "DNA/Academ": "#FF0000",
    "DNA/CMC": "#FF200B",
    "DNA/Crypton": "#FF3115",
    "DNA/Ginger": "#FF3D1E",
    "DNA/Harbinger": "#FF4825",
    "DNA/hAT": "#FF512D",
    "DNA/Kolobok": "#FF5A34",
    "DNA/Maverick": "#FF623B",
    "DNA": "#FF6A42",
    "DNA/Merlin": "#FF7149",
    "DNA/MULE": "#FF7850",
    "DNA/P": "#FF7F57",
    "DNA/PiggyBac": "#FF865E",
    "DNA/Sola": "#FF8D65",
    "DNA/TcMar": "#FF936C",
    "DNA/Transib": "#FF9972",
    "DNA/Zator": "#FF9F79",
    "DNA/Dada": "#FFCFBC",
    "RC": "#ff2987",
    "RC/Helitron": "#ff2987",
    "PLE": "#ffce6e",
    "Retroposon/SVA": "#FF4D4D",
    "LTR/DIRS": "#006400",
    "LTR/Ngaro": "#197214",
    "LTR/Pao": "#2A8024",
    "LTR/Copia": "#3A8F33",
    "LTR/Gypsy": "#489E42",
    "LTR/ERVL": "#57AE51",
    "LTR": "#65BD61",
    "LTR/ERV1": "#73CD70",
    "LTR/ERV": "#81DD80",
    "LTR/ERVK": "#90ED90",
    "LINE/L1": "#00008B",
    "LINE": "#4594c7",
    "LINE/RTE": "#38299A",
    "LINE/CR1": "#483AA2",
    "LINE/Rex-Babar": "#554BAA",
    "LINE/L2": "#625CB1",
    "LINE/Proto2": "#6E6DB9",
    "LINE/LOA": "#797EC0",
    "LINE/R1": "#848FC8",
    "LINE/Jockey-I": "#8FA1CF",
    "LINE/Dong-R4": "#99B3D7",
    "LINE/R2": "#A3C5DE",
    "LINE/CRE": "#C1D9FF",
    "SINE": "#9F1FF0",
    "SINE/5S": "#A637F1",
    "SINE/7SL": "#AD49F2",
    "SINE/Alu": "#B358F3",
    "SINE/tRNA": "#B966F4",
    "SINE/tRNA-Alu": "#BF74F4",
    "SINE/tRNA-RTE": "#C481F5",
    "SINE/RTE": "#C98EF6",
    "SINE/Deu": "#CE9BF7",
    "SINE/tRNA-V": "#D3A7F7",
    "SINE/MIR": "#D7B4F8",
    "SINE/U": "#DFCDF9",
    "SINE/tRNA-7SL": "#E2D9F9",
    "SINE/tRNA-CR1": "#E5E5F9",
}

#----------------------------------------------


BDELLOIDS = ['a_ricciae', 
             'a_sp_wild', 
             'a_vaga', 
             'h_sp_wild', 
             'm_quadricornifera', 
             'r_rotatoria', 
             'p_roseola']

#----------------------------------------------


def parse_landscape_divsum(file_path: str) -> pd.DataFrame:
    """
    Parses the landscape divsum file to extract the coverage for each repeat class and divergence.
    
    Args:
        file_path (str): Path to the landscape divsum file.
        
    Returns:
        pd.DataFrame: DataFrame containing the coverage data.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('Coverage for each repeat class and divergence'):
            return pd.read_csv(StringIO(''.join(lines[i + 1:])), sep=' ')
    raise ValueError("No coverage section found in the file.")
    

def get_name_map(name_map: dict, graph_label: dict, df_columns: list) -> dict:
    """Updates the name map to include any columns in the dataframe that are not already mapped.
    Args:
        name_map (dict): Existing name map.
        graph_label (dict): Graph labels.
        df_columns (list): Columns in a dataframe.
        
    Returns:
        dict: Updated name map.
    """
    new_map = name_map.copy()
    for col in df_columns:
        if col == 'Div':
            continue
        if not col in name_map and col in graph_label:
            new_map[col] = col
        if not col in name_map and not col in graph_label and col.split('/')[0] in graph_label:
            new_map[col] = col.split('/')[0]
    return new_map
    
   
def plot_repeat_landscape(repeat_landscape_file: str, output_file: str, species_details_file: str) -> None:
    """
    Plots the repeat landscape from the given file.
    
    Args:
        repeat_landscape_file (str): Path to the repeat landscape file.
        output_file (str): Path to save the output plot.
        species_details_file (str): Path to the species details JSON file.
    """
    with open(species_details_file, 'r') as f:
        species_details = json.load(f)
    assembly_size = species_details.get('assembly_size', 1)
    scientific_name = species_details.get('scientific_name')
    
    df = parse_landscape_divsum(repeat_landscape_file)

    # Update name map to include current unmapped columns
    updated_name_map = get_name_map(NAME_MAP, GRAPH_LABELS, df.columns)
    
    for col in df.columns:
        if col != 'Div' and col not in updated_name_map:
            print(f"Warning: Column '{col}' not found in NAME_MAP or GRAPH_LABELS. It will be ignored.")
    
    # Filter dataframe to include only relevant columns
    df = df[[col for col in df.columns if col == 'Div' or col in updated_name_map]]
    
    # Sum columns based on updated name map
    df.rename(columns=updated_name_map, inplace=True)
    df = df.groupby(df.columns, axis=1).sum()
    
    # change the column order to match GRAPH_LABELS in reverse order
    ordered_cols = ['Div'] + [rc for rc in reversed(GRAPH_LABELS.keys()) if rc in df.columns and rc != 'Div']
    df = df[ordered_cols]
    
    # Prepare data for plotting
    df_melted = df.melt(id_vars=['Div'], var_name='Repeat_Class', value_name='Coverage')
    df_melted['Coverage'] = (df_melted['Coverage'] / assembly_size) * 100
    df_melted = df_melted[~df_melted['Repeat_Class'].isin([
        'Satellite', 
        'Segmental', 
        'Structural_RNA', 
        'Simple_repeat', 
        'rRNA', 
        'tRNA', 
        'Unknown',
        'Unnamed'])]
    
    # Plotting        
    plt.figure(figsize=(10, 6))
    sns.set_palette(sns.color_palette([GRAPH_LABELS.get(rc, '#000000') for rc in df_melted['Repeat_Class'].unique()]))
    sns.histplot(
        data=df_melted,
        x='Div',
        weights='Coverage',     
        multiple='stack',    
        hue='Repeat_Class',
        linewidth=0,      
        edgecolor=None,    
        shrink=1,
        discrete=True
    ) 
    legend_elements = [
        Patch(facecolor=GRAPH_LABELS['DNA'], label='DNA transposons'),
        Patch(facecolor=GRAPH_LABELS['RC/Helitron'], label='Rolling circle'),
        Patch(facecolor=GRAPH_LABELS['PLE'], label='PLE'),
        Patch(facecolor=GRAPH_LABELS['LTR'], label='LTR'),
        Patch(facecolor=GRAPH_LABELS['LINE'], label='LINE'),
        Patch(facecolor=GRAPH_LABELS['SINE'], label='SINE'),
    ]
    plt.legend(handles=legend_elements, title='Repeat type', loc='upper right', 
               frameon=True, facecolor='#f2f2f2', edgecolor='none', fontsize=12, title_fontsize=13)               
    plt.xlim(-2, 50)
    if species_details.get('short_name') in BDELLOIDS:
        plt.ylim(0, 1)
    plt.xlabel('Divergence (Kimura distance)', fontsize=14)
    plt.ylabel('Genome span (%)', fontsize=14)
    plt.title(scientific_name, fontsize=16, fontstyle='italic')
    plt.gca().set_axisbelow(True)
    plt.grid(axis='y', color='lightgrey', linestyle='--', linewidth=0.7)
    plt.grid(axis='x', color='lightgrey', linestyle='--', linewidth=0.7)
    plt.tight_layout()
    plt.savefig(output_file, format='svg')  
    

##----- MAIN -----

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot repeat landscape from RepeatMasker divsum file.")
    parser.add_argument("--rm_divsum", type=str, help="Path to the repeat landscape divsum file.")
    parser.add_argument("--output", type=str, help="Path to save the output plot (SVG format).")
    parser.add_argument("--sp_details", type=str, help="Path to the species details JSON file.")
    args = parser.parse_args()
    plot_repeat_landscape(args.rm_divsum, args.output, args.sp_details)
    
     
# # TEST
# sp = 'm_quadricornifera'
# landscape_divsum = f'/home/odethier/Desktop/repeat_landscape/{sp}/{sp}.assembly.fa.divsum'
# output_plot = f'/home/odethier/Desktop/repeat_landscape/{sp}/{sp}_repeat_landscape_plot.svg'
# sp_details = f'/home/odethier/Desktop/repeat_landscape/{sp}/{sp}.details.json'
# plot_repeat_landscape(landscape_divsum, output_plot, sp_details)