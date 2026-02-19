import os.path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

"""
Plot the coverage along a sequence from a bedgraph file
with a smoothed line and data spread (min-max).
"""

WINDOW_SIZE = 100_000


def plot_coverage(bedgraph_file, scaffolds, out_dir):
    print('Loading begraph...')
    bedgraph_df = pd.read_csv(bedgraph_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])
    
    for scaffold_name in scaffolds:
        print(f'Processing scaffold {scaffold_name}...')
    
        df = bedgraph_df[bedgraph_df['chrom'] == scaffold_name].copy() # Use .copy() to avoid SettingWithCopyWarning
        
        # Calculate rolling statistics (Mean and Min-Max)
        df["smoothed"] = df["coverage"].rolling(window=WINDOW_SIZE, center=True).mean()
        df["min"] = df["coverage"].rolling(window=WINDOW_SIZE, center=True).min()
        df["max"] = df["coverage"].rolling(window=WINDOW_SIZE, center=True).max()
        
        # Calculate center position for plotting
        df["x_pos"] = (df['start'] + df['end']) / 2
        
        # Plotting
        plt.figure(figsize=(10, 5))
        
        # Plot smoothed coverage line
        sns.lineplot(x='x_pos', 
                    y='smoothed', 
                    data=df, 
                    color='blue', 
                    lw=1.5, 
                    alpha=0.6, 
                    label='Coverage (window size: 100kb)', 
                    rasterized=True
        )
        
        # Fill min-max area
        plt.fill_between(
            df['x_pos'], 
            df['min'], 
            df['max'], 
            color='blue', 
            alpha=0.2,          
            edgecolor='none',   
            label='Coverage spread (min-max)',
            rasterized=True 
        )
    
        # Plot Average Line
        plt.axhline(y=df['coverage'].mean(), color='gray', linestyle='--', lw=1, label='Average coverage')
        
        # Formatting
        plt.title(f'{scaffold_name} coverage', fontsize=16)
        plt.xlabel('Position (bp)', fontsize=14)
        plt.ylabel('Coverage', fontsize=14)
        plt.ylim(0, 500)
        
        plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(5_000_000))
        plt.gca().set_axisbelow(True)

        formatter = ticker.FuncFormatter(lambda x, _: f'{x/1e6:.0f}M' if x != 0 else '0')
        plt.gca().xaxis.set_major_formatter(formatter)
        
        plt.grid(axis='y', color='lightgrey', linestyle='--', linewidth=0.7)
        plt.grid(axis='x', color='lightgrey', linestyle='--', linewidth=0.7)
        
        plt.legend(frameon=True, facecolor='#f2f2f2', edgecolor='none', loc='upper right') 
        plt.tight_layout()
        
        plt.savefig(os.path.join(out_dir, f'{scaffold_name}_coverage.svg'), dpi=300)
        plt.close()


# Get bedgraph:
# minimap2 -t 12 -ax lr:hq assembly.fa ont.fq.gz --secondary=no | samtools sort -@12 -O BAM -o ont_aligned.bam -
# bedtools genomecov -bg -ibam ont_aligned.bam > ont_coverage.bedgraph


species = 'h_sp_wild'
bedgraph_file = f'/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/metrics/species/{species}/coverage/{species}.ont_coverage.bedgraph'
scaffolds = ['Chrom_1', 'Chrom_2', 'Chrom_3', 'Chrom_4', 'Chrom_5', 'Chrom_6']
out_dir = f'/mnt/sdb1/Olivier/comparative_genomics/pipeline/data/metrics/species/{species}/coverage'
plot_coverage(bedgraph_file, scaffolds, out_dir)
