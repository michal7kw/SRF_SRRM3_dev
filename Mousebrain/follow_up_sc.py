import loompy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict
import gc
import os


# Set working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain')

def get_gene_data(loom_path: str, genes: List[str]) -> pd.DataFrame:
    """
    Extract data for multiple genes from loom file.
    """
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
                'PseudoAge': ds.ca['PseudoAge'],
                'Class': ds.ca['Class'],
                'ClusterName': ds.ca['ClusterName'],
                'Tissue': ds.ca['Tissue'],
                'Region': ds.ca['Region'],
                'Subclass': ds.ca['Subclass']
            })
            
            all_data.append(gene_data)
    
    return pd.concat(all_data, ignore_index=True)

def classify_detailed_cell_type(row: pd.Series) -> str:
    """
    Classify cells into detailed types matching Figure 6f.
    """
    cluster = str(row['ClusterName'])
    if cluster.startswith('Rgl'):
        if 'early' in str(row['Subclass']).lower():
            return 'Early radial glia'
        return 'Late radial glia'
    elif cluster.startswith('Gbl'):
        return 'Glioblasts'
    elif cluster.startswith('OPC') or cluster.startswith('PreOPC'):
        return 'OPCs'
    elif cluster.startswith('Nbl'):
        return 'Neuroblasts'
    elif cluster.startswith('Neur'):
        return 'Neurons'
    return 'Other'

def plot_expression_dot_plot(df: pd.DataFrame, output_dir: str):
    """
    Create dot plot similar to Figure 6e.
    """
    # Calculate statistics for each gene-celltype combination
    stats = df.groupby(['Gene', 'DetailedType']).agg({
        'Expression': ['mean', lambda x: (x > 0).mean() * 100]
    }).reset_index()
    stats.columns = ['Gene', 'DetailedType', 'mean_expr', 'pct_expressing']
    
    # Order cell types and genes
    cell_type_order = ['Early radial glia', 'Late radial glia', 'Glioblasts', 
                       'OPCs', 'Neuroblasts', 'Neurons']
    
    # Create dot plot
    plt.figure(figsize=(12, 8))
    
    # Create scatter plot where size represents percentage expressing
    # and color represents mean expression
    for idx, gene in enumerate(stats['Gene'].unique()):
        gene_data = stats[stats['Gene'] == gene]
        
        plt.scatter(gene_data['DetailedType'], [idx] * len(gene_data), 
                   s=gene_data['pct_expressing'] * 5,  # Scale size
                   c=gene_data['mean_expr'], cmap='viridis',
                   alpha=0.7)
    
    plt.yticks(range(len(stats['Gene'].unique())), stats['Gene'].unique())
    plt.xticks(rotation=45, ha='right')
    plt.title('Gene Expression Patterns Across Cell Types')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/gene_expression_dotplot.pdf', bbox_inches='tight')
    plt.close()

def plot_expression_trajectories(df: pd.DataFrame, output_dir: str):
    """
    Create trajectory plots showing gene expression across pseudotime.
    """
    cell_type_order = ['Early radial glia', 'Late radial glia', 'Glioblasts', 
                       'OPCs', 'Neuroblasts', 'Neurons']
    
    # Plot for each gene
    for gene in df['Gene'].unique():
        gene_df = df[df['Gene'] == gene]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), height_ratios=[2, 1])
        
        # Trajectory plot
        sns.lineplot(data=gene_df, x='PseudoAge', y='Expression', 
                    hue='DetailedType', ax=ax1)
        ax1.set_title(f'{gene} Expression Trajectory')
        ax1.set_xlabel('Pseudotime')
        ax1.set_ylabel('Expression')
        
        # Box plot
        sns.boxplot(data=gene_df, x='DetailedType', y='Expression', 
                   order=cell_type_order, ax=ax2)
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_trajectory.pdf', bbox_inches='tight')
        plt.close()

# Define genes to analyze
genes = [
    'Srrm3', 
    'Srrm4',
    # Marker genes from Figure 6e
    'Pax6',    # Radial glia marker
    'Sox2',    # Neural progenitor marker
    'Neurog2', # Neuroblast marker
    'Eomes',   # Intermediate progenitor marker
    'Neurod2', # Neuronal differentiation marker
    'Tubb3',   # Neuron marker
]

# Main execution
loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.loom"
plots_dir = 'plots_trajectories_sc'

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Get data
df = get_gene_data(loom_path, genes)

# Add detailed cell type classification
df['DetailedType'] = df.apply(classify_detailed_cell_type, axis=1)

# Filter for relevant cell types
df = df[df['DetailedType'] != 'Other']

# Filter for forebrain dorsal tissue if needed
df = df[df['Tissue'] == 'ForebrainDorsal']

# Create plots
plot_expression_dot_plot(df, plots_dir)
plot_expression_trajectories(df, plots_dir)

# Save processed data
df.to_csv(f'{plots_dir}/gene_expression_data.csv', index=False)

# Generate summary statistics
stats = df.groupby(['Gene', 'DetailedType']).agg({
    'Expression': ['mean', 'std', 'count', lambda x: (x > 0).mean() * 100]
}).round(3)
stats.columns = ['mean', 'std', 'n_cells', 'percent_expressing']
stats.to_csv(f'{plots_dir}/expression_statistics.csv')

print("Analysis complete. Results saved to:", plots_dir)
gc.collect()