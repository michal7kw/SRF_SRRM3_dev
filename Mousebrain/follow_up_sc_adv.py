import loompy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict
from scipy.stats import zscore
import gc
import os

# set the working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain')

def classify_detailed_cell_type(row: pd.Series) -> str:
    """
    Classify cells into detailed types matching Figure 6f.
    
    Args:
        row: Pandas Series containing cell metadata
        
    Returns:
        str: Classified cell type
    """
    cluster = str(row['ClusterName'])
    class_name = str(row['Class'])
    
    if cluster.startswith('Rgl'):
        if any(x in str(row['Subclass']).lower() for x in ['early', 'neural progenitor']):
            return 'Early radial glia'
        return 'Late radial glia'
    elif cluster.startswith('Gbl'):
        return 'Glioblasts'
    elif cluster.startswith(('OPC', 'PreOPC')):
        return 'OPCs'
    elif cluster.startswith('Nbl') or class_name == 'Neuroblast':
        return 'Neuroblasts'
    elif cluster.startswith('Neur') or class_name == 'Neuron':
        return 'Neurons'
    return 'Other'

def get_gene_data(loom_path: str, genes: List[str]) -> pd.DataFrame:
    """Extract data for multiple genes with UMAP coordinates."""
    all_data = []
    
    with loompy.connect(loom_path) as ds:
        for gene in genes:
            gene_idx = np.where(ds.ra.Gene == gene)[0]
            if len(gene_idx) == 0:
                print(f"Warning: {gene} not found in dataset")
                continue
                
            gene_expr = ds[gene_idx[0], :]
            
            # Create data for this gene
            gene_data = pd.DataFrame({
                'Gene': gene,
                'Expression': gene_expr,
                'Expression_zscore': zscore(gene_expr, nan_policy='omit'),
                'PseudoAge': ds.ca['PseudoAge'],
                'Class': ds.ca['Class'],
                'ClusterName': ds.ca['ClusterName'],
                'Tissue': ds.ca['Tissue'],
                'Region': ds.ca['Region'],
                'Subclass': ds.ca['Subclass'],
                'UMAP1': ds.ca['UMAP'][:, 0],
                'UMAP2': ds.ca['UMAP'][:, 1]
            })
            
            all_data.append(gene_data)
    
    return pd.concat(all_data, ignore_index=True)

def plot_umap_clusters(df: pd.DataFrame, output_dir: str):
    """Create UMAP visualization with cluster labels."""
    plt.figure(figsize=(12, 10))
    
    # Plot UMAP points colored by cell type
    scatter = plt.scatter(df['UMAP1'], df['UMAP2'], 
                         c=pd.factorize(df['DetailedType'])[0],
                         cmap='tab20', s=1, alpha=0.6)
    
    # Add legend
    handles = [plt.scatter([], [], c=c, label=l) 
              for c, l in zip(plt.cm.tab20(np.linspace(0, 1, len(df['DetailedType'].unique()))),
                            df['DetailedType'].unique())]
    plt.legend(handles=handles, title='Cell Type', 
              bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.title('UMAP Visualization of Cell Types')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umap_clusters.pdf', bbox_inches='tight')
    plt.close()

def plot_gene_expression_umap(df: pd.DataFrame, output_dir: str):
    """Create UMAP visualization colored by gene expression."""
    for gene in df['Gene'].unique():
        gene_df = df[df['Gene'] == gene]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
        # Raw expression
        scatter1 = ax1.scatter(gene_df['UMAP1'], gene_df['UMAP2'],
                             c=gene_df['Expression'], 
                             cmap='viridis', s=1, alpha=0.6)
        plt.colorbar(scatter1, ax=ax1, label='Expression')
        ax1.set_title(f'{gene} Raw Expression')
        ax1.set_xlabel('UMAP1')
        ax1.set_ylabel('UMAP2')
        
        # Z-scored expression
        scatter2 = ax2.scatter(gene_df['UMAP1'], gene_df['UMAP2'],
                             c=gene_df['Expression_zscore'],
                             cmap='viridis', s=1, alpha=0.6)
        plt.colorbar(scatter2, ax=ax2, label='Z-score')
        ax2.set_title(f'{gene} Z-scored Expression')
        ax2.set_xlabel('UMAP1')
        ax2.set_ylabel('UMAP2')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_umap_expression.pdf', 
                   bbox_inches='tight')
        plt.close()

def plot_expression_timeline(df: pd.DataFrame, output_dir: str):
    """Create timeline plots similar to Figure 9."""
    for gene in df['Gene'].unique():
        gene_df = df[df['Gene'] == gene]
        
        # Create figure with grid specification
        fig = plt.figure(figsize=(15, 10))
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
        
        # Timeline plot (top)
        ax1 = fig.add_subplot(gs[0, :])
        for cell_type in gene_df['DetailedType'].unique():
            cell_data = gene_df[gene_df['DetailedType'] == cell_type]
            
            # Calculate binned average using pd.cut
            bins = np.arange(cell_data['PseudoAge'].min(), 
                           cell_data['PseudoAge'].max() + 0.2, 0.2)
            cell_data['PseudoAge_bin'] = pd.cut(cell_data['PseudoAge'], bins)
            rolling_mean = cell_data.groupby('PseudoAge_bin')['Expression'].mean()
            
            # Plot using bin centers for x-axis
            bin_centers = [(bin.left + bin.right)/2 for bin in rolling_mean.index]
            ax1.plot(bin_centers, rolling_mean.values, label=cell_type, alpha=0.8)
            
        ax1.set_title(f'{gene} Expression Timeline')
        ax1.set_xlabel('Pseudotime')
        ax1.set_ylabel('Expression')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Distribution plot (bottom left)
        ax2 = fig.add_subplot(gs[1, 0])
        sns.violinplot(data=gene_df, x='DetailedType', y='Expression',
                      ax=ax2, cut=0)
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
        ax2.set_title('Expression Distribution')
        
        # Heatmap (bottom right)
        ax3 = fig.add_subplot(gs[1, 1])
        pivot_df = gene_df.pivot_table(
            values='Expression',
            index='DetailedType',
            columns=pd.qcut(gene_df['PseudoAge'], 10, labels=False),
            aggfunc='mean'
        )
        sns.heatmap(pivot_df, cmap='viridis', ax=ax3)
        ax3.set_title('Expression Heatmap')
        ax3.set_xlabel('Pseudotime Bins')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_timeline.pdf', 
                   bbox_inches='tight', dpi=300)
        plt.close()

# Define genes to analyze (including markers from Fig 9)
genes = [
    'Srrm3', 
    'Srrm4',
    # Astrocyte markers
    'Slc6a11', 'Aqp4', 'Aldoc', 'Hes5', 'Hes1',
    # OPC markers
    'Pdgfra', 'Sox10', 'Olig2', 'Olig1',
    # Neuronal markers
    'Neurod2', 'Tubb3', 'Neurog2', 'Eomes',
    # Radial glia markers
    'Pax6', 'Sox2'
]

# Main execution
loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.loom"
plots_dir = 'plots_trajectories_refined'

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Get data and process
df = get_gene_data(loom_path, genes)
df['DetailedType'] = df.apply(classify_detailed_cell_type, axis=1)
df = df[df['DetailedType'] != 'Other']
df = df[df['Tissue'] == 'ForebrainDorsal']

# Create visualizations
# plot_umap_clusters(df, plots_dir)
# plot_gene_expression_umap(df, plots_dir)
plot_expression_timeline(df, plots_dir)

# Save processed data
df.to_csv(f'{plots_dir}/gene_expression_data.csv', index=False)

# Generate summary statistics
stats = df.groupby(['Gene', 'DetailedType']).agg({
    'Expression': ['mean', 'std', 'count', 
                  lambda x: (x > 0).mean() * 100],
    'Expression_zscore': ['mean', 'std']
}).round(3)
stats.to_csv(f'{plots_dir}/expression_statistics.csv')

print("Analysis complete. Results saved to:", plots_dir)
gc.collect()