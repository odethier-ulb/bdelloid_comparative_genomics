import re
import sys
import json
import pandas as pd
from pycirclize import Circos

"""
Create a circos plot of genomic features.

Each genomic feature (TEs, HGT candidates, Genes) is plotted across windows of a predefined size. 
The plotted value represents the feature density within that window, calculated as the percentage of bases 
covered by the respective feature.
"""

WINDOW_SIZE = 100000 # plot features in windows of this size
TE_REPEATS = ['DNA', 'RC', 'PLE', 'PLE-Athena', 'LTR', 'LINE', 'SINE'] # known TE classes
CHR_NAME_PATTERN = r'(?:(?<=\D)|(?<=^))(\d+)$' # pattern to extract chromosome number from chromosome name
COLORS = ['#9ac8e0', '#a2d7a3', '#fdb07a']
YTICKS_100 = [0, 20, 40, 60, 80, 100]
YTICKS_25 = [0, 5, 10, 15, 20, 25]


def __get_chr_to_plot(sp_fai, sp_details):
    """Return an ordered list of chromosome names to plot with their size and coresponding plot name.

    Args:
        sp_fai (str): Path to the species FAI file from Samtools faidx.
        sp_details (str): Path to the species details file.
    """
    # Parse the FAI file to get chromosome names and sizes
    with open(sp_fai, 'r') as f:
        chrom_sizes = {line.split('\t')[0]: int(line.split('\t')[1]) for line in f}
    # Parse the species details file to get the chromosomes to plot 
    with open(sp_details, 'r') as f:
        details = json.load(f)
    # Combine the two to get the chromosome names, sizes, and plot names
    chr_to_plot = []
    for chr in details['circos_order']:
        if chr in chrom_sizes:
            match = re.search(CHR_NAME_PATTERN, chr)
            plot_name = f'{details["accronym"]}{match.group(1) if match else chr}'         
            chr_to_plot.append((chr, chrom_sizes[chr], plot_name))
    return chr_to_plot


def __get_re_density(te_bed, chr, chr_size):
    """Get the repetitive elements density computed by window size along a chromosome.

    Args:
        te_bed (str): Path to the RepeatMasker TE bed.
        chr (str): Chromosome name.
        chr_size (str): Chromosome size.
    """
    # Init
    chr_mask = [False] * chr_size # mask of bases already covered to prevent multiple counting 
    x = range(0, chr_size + WINDOW_SIZE, WINDOW_SIZE)
    y_te = [0] * len(x) # known transposable elements
    y_re = [0] * len(x) # all repetitive elements
    df = pd.read_csv(te_bed, sep='\t', header=None, 
                     names=['chr', 'start', 'end', 'skip1', 'skip2', 'skip3', 'class', 'skip4', 'skip5', 'skip6'])
    
    # Retrieve the density
    for _, row in df[df['chr'] == chr].sort_values('start').iterrows():
        for i in range(row['start'], row['end']):
            if chr_mask[i]:
                continue
            idx = i // WINDOW_SIZE
            if idx < len(y_re):
                if row['class'] in TE_REPEATS:
                    y_te[idx] += 1
                y_re[idx] += 1
            chr_mask[i] = True
            
    # Convert to percentage
    y_te = [i / WINDOW_SIZE * 100 for i in y_te]
    y_re = [i / WINDOW_SIZE * 100 for i in y_re]
    
    # Arrange values for plotting
    x_final, y_te_final, y_re_final = [], [], []
    for i in range(len(x)):
        if x[i] >= chr_size:
            break 
        x_final.extend([x[i], min(chr_size, x[i] + WINDOW_SIZE - 1)])
        y_te_final.extend([y_te[i], y_te[i]])
        y_re_final.extend([y_re[i], y_re[i]])
    return x_final, y_te_final, y_re_final  
            
        
def __get_gene_density(mcscanx_gff, chr, chr_size, htgs, hgt_only=False):
    """Get the gene density computed by window size along a chromosome.

    Args:
        mcscanx_gff (str): Path to the MCScanX GFF file.
        chr (str): Chromosome name.
        chr_size (int): Chromosome size
        htgs (set): Set of HTGs.
        hgt_only (bool, optional): Plot only HTGs.
    """
    # Init
    x = range(0, chr_size + WINDOW_SIZE, WINDOW_SIZE)
    y = [0] * len(x)
    df = pd.read_csv(mcscanx_gff, sep='\t', header=None, names=['chr', 'gene', 'start', 'end'])
     
    # Retrieve the density
    for _, row in df[df['chr'] == chr].sort_values('start').iterrows():
        # -- check if we consider the gene --#
        gene_name = row['gene'].split(f'{chr}-')[-1]
        if (hgt_only and gene_name not in htgs) or (not hgt_only and gene_name in htgs):
            continue 
        # -- fill the concerned window(s) -- #
        for i in range(row['start'], row['end'] + 1):
            y[min(i // WINDOW_SIZE, len(y) - 1)] += 1
                      
    # Convert to percentage
    y = [i / WINDOW_SIZE * 100 for i in y]
    
    # Arrange values for plotting
    x_final, y_final = [], []
    for i in range(len(x)):
        if x[i] >= chr_size:
            break 
        x_final.extend([x[i], min(chr_size, x[i] + WINDOW_SIZE - 1)])
        y_final.extend([y[i], y[i]])
    return x_final, y_final


def __get_htgs(hgt_results):
    """Return a set of HTGs 

    Args:
        hgt_results (str): Path to the HGT results file.
    """
    df = pd.read_csv(hgt_results, sep='\t', header=0)
    return set(df[df['hgt_class'] >= 3]['prot_id'])


def __get_links(mcscanx_gff, mcscanx_col):
    """Retrieve syntenic links from MCScanX results.
    
    Args:
        mcscanx_gff (str): MScanX input file.
        mcscanx_collinearity (str): MScanX collinearity output file.
    """
    # Parse input file to get gene positions
    df_gff, gene_pos = pd.read_csv(mcscanx_gff, sep='\t', header=None, names=['chr', 'gene', 'start', 'end']), {}
    for _, row in df_gff.iterrows():
        gene_pos[row['gene']] = (row['chr'], row['start'], row['end'])
    # Parse collinearity file to get syntenic 'blocks'
    col_blocks = []
    with open(mcscanx_col, 'r') as f:
        for line in f:
            if line.startswith('## Alignment'):
                col_blocks.append([])
            elif not line.startswith('#'):
                col_blocks[-1].append(line)
    # Transform blocks into regions
    regions = []
    for block in col_blocks:
        blk_start = block[0].split('\t')
        start_offset = 0 if len(blk_start) == 4 else 1 
        blk_end = block[-1].split('\t')
        end_offset = 0 if len(blk_start) == 4 else 1 
        # region 1
        chr_1 = gene_pos[blk_start[1 + start_offset]][0]
        start_1 = gene_pos[blk_start[1 + start_offset]][1]
        end_1 = gene_pos[blk_end[1 + end_offset]][2]
        # region 2
        chr_2 = gene_pos[blk_start[2 + start_offset]][0]
        start_2 = gene_pos[blk_start[2 + start_offset]][1]
        end_2 = gene_pos[blk_end[2 + end_offset]][2]
        # add regions
        regions.append([(chr_1, start_1, end_1), (chr_2, start_2, end_2)]) 
    return regions


def plot_circos(sp_fai, sp_details, hgt_results, te_bed, mcscanx_gff, mcscanx_col, out_file):
    """Plot a circos plot of genomic features.

    Args:
        sp_fai (str): Path to the species FAI file from Samtools faidx.
        sp_details (str): Path to the species details file.
        hgt_results (str): Path to the HGT results file.
        te_bed (str): Path to the TE BED file.
        mcscanx_gff (str): Path to the MCScanX GFF file.
        mcscanx_col (str): Path to the MCScanX collinearity result file.
        out_file (str): Path to the output file.
    """
    # Retrieve data
    chr_to_plot = __get_chr_to_plot(sp_fai, sp_details)
    htgs = __get_htgs(hgt_results)
    
    # Init Circos object
    circos = Circos({chr[-1]: chr[1] for chr in chr_to_plot}, space=3, start=0, end=359)
    
    # Draw features for each chromosome
    for chr in chr_to_plot:
        sector = circos.get_sector(chr[-1])
        
        # -- chr name -- #
        name_track = sector.add_track((91, 100))
        name_track.text(sector.name, r=91, size=9, ha="center", va="center")
        
        # --  TEs -- #
        te_track = sector.add_track((81, 89), r_pad_ratio=0.1)
        te_track.axis()
        te_track.grid()
        te_track.yticks(YTICKS_100, labels=YTICKS_100, vmin=0, vmax=100, side='right', 
                        tick_length=.3, label_size=3, label_margin=.2)
        x, y_te, y_re = __get_re_density(te_bed, chr[0], chr[1])
        te_track.line(x, y_re, vmax=YTICKS_100[-1], color="gray", alpha=0.8, lw=0)
        te_track.fill_between(x, y_re, vmax=YTICKS_100[-1], color="gray", alpha=0.6)
        te_track.line(x, y_te, vmax=YTICKS_100[-1], color="black", alpha=0.8, lw=0.5)
        te_track.fill_between(x, y_te, vmax=YTICKS_100[-1], color='#4eaf4b', alpha=0.8)
        
        # -- Genes -- #
        gene_track = sector.add_track((71, 79), r_pad_ratio=0.1)
        gene_track.axis()
        gene_track.grid()
        gene_track.yticks(YTICKS_100, labels=YTICKS_100, vmin=0, vmax=YTICKS_100[-1], side='right', 
                        tick_length=.3, label_size=3, label_margin=.2)
        x, y = __get_gene_density(mcscanx_gff, chr[-1], chr[1], htgs, hgt_only=False)
        gene_track.line(x, y, vmax=YTICKS_100[-1], color="black", alpha=0.8, lw=0.5)
        gene_track.fill_between(x, y, vmax=YTICKS_100[-1], color='#387fb7', alpha=0.8)
        
        # -- HGTs -- #
        hgt_track = sector.add_track((61, 69), r_pad_ratio=0.1)
        hgt_track.axis()
        hgt_track.grid()
        hgt_track.yticks(YTICKS_25, labels=YTICKS_25, vmin=0, vmax=YTICKS_25[-1], side='right', 
                        tick_length=.3, label_size=3, label_margin=.2)
        x, y = __get_gene_density(mcscanx_gff, chr[-1], chr[1], htgs, hgt_only=True)
        y = [min(i, YTICKS_25[-1]) for i in y] # keep values in range
        hgt_track.line(x, y, vmax=YTICKS_25[-1], color="black", alpha=0.8, lw=0.5)
        hgt_track.fill_between(x, y, vmax=YTICKS_25[-1], color='#ec1920', alpha=0.8)
        
        # -- Coordinates -- #
        track2 = sector.add_track((54.8, 55))
        track2.axis(fc="black")
        interval = 2e6
        track2.xticks_by_interval(interval, outer=True, label_size=4,
                                 label_formatter=lambda v: f"{v / 1e6:.0f}M",
                                 text_kws=dict(alpha=0.7), label_orientation="vertical")
        
        # -- Track label -- #
        if sector.name == circos.sectors[0].name:
            circos.text("TEs", r=te_track.r_center, deg=-1.5, color="black", size=6, rotation=-90, rotation_mode='anchor')
            circos.text("Genes", r=gene_track.r_center, deg=-1.5, color="black", size=6, rotation=-90, rotation_mode='anchor')
            circos.text("HGTs", r=hgt_track.r_center, deg=-1.5, color="black", size=6, rotation=-90, rotation_mode='anchor')
        
    # Draw collinearity blocks
    for regions in __get_links(mcscanx_gff, mcscanx_col):
        color = COLORS[int(regions[0][0][2:]) % 3]
        circos.link(regions[0], regions[1], color=color, r1=53, r2=53, alpha=0.7, ec="black", lw=0.5)
    
    # Save
    _ = circos.plotfig()
    circos.savefig(out_file, pad_inches=0)



#----------------------- RUN -----------------------#

# Input files
sp_fai = 'a_ricciae.assembly.fa.fai'
sp_details =  'a_ricciae.details.json'
hgt_results = 'a_ricciae_hgt_results_nr.tsv'
te_bed = 'a_ricciae_rm.bed'
mcscanx_gff = 'a_ricciae_mcscanx.gff'
mcscanx_col = 'a_ricciae_mcscanx.collinearity'
out_file = 'a_ricciae_circos.svg'

plot_circos(sp_fai, sp_details, hgt_results, te_bed, mcscanx_gff, mcscanx_col, out_file)
