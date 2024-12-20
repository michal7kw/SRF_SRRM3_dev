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

def get_gene_expression(ds, gene_name: str) -> np.ndarray:
    """Extract expression data for a given gene."""
    gene_idx = np.where(ds.ra['Gene'] == gene_name)[0]
    if len(gene_idx) == 0:
        raise ValueError(f"{gene_name} not found in dataset")
    return ds[gene_idx[0], :]

def extract_age(age_col: str) -> float:
    """Extract age value from column name."""
    return float(age_col.replace('Age_e', '').replace('_0', '').replace('_', '.'))

def classify_cell_type(cluster_name: str) -> str:
    """Map cluster names to main cell types matching Figure 6f."""
    if cluster_name.startswith('Rgl'):
        return 'Radial glia'
    elif cluster_name.startswith('Nbl'):
        return 'Neuroblasts'
    elif cluster_name.startswith('Neur'):
        return 'Neurons'
    elif cluster_name.startswith('OPC') or cluster_name.startswith('PreOPC'):
        return 'OPCs'
    elif cluster_name.startswith('Gbl'):
        return 'Glioblasts'
    else:
        return 'Other'

def extract_cell_data(ds, genes: List[str]) -> pd.DataFrame:
    """
    Extract expression data for multiple genes with detailed cell metadata.
    """
    data = []
    
    # Get age columns
    age_cols = [col for col in ds.ca.keys() if col.startswith('Age_e')]
    
    for gene in genes:
        try:
            gene_expr = get_gene_expression(ds, gene)
            
            # For each cell
            for cell_idx in range(ds.shape[1]):
                # Get age
                cell_age = None
                for age_col in age_cols:
                    if ds.ca[age_col][cell_idx] == 1:
                        cell_age = extract_age(age_col)
                        break
                
                if cell_age is not None:
                    cluster_name = ds.ca['ClusterName'][cell_idx]
                    main_type = classify_cell_type(str(cluster_name))
                    
                    if main_type != 'Other':  # Only include relevant cell types
                        data.append({
                            'Gene': gene,
                            'Expression': gene_expr[cell_idx],
                            'Age': cell_age,
                            'Class': ds.ca['Class'][cell_idx],
                            'MainType': main_type,
                            'ClusterName': cluster_name,
                            'Subclass': ds.ca['Subclass'][cell_idx]
                        })
                    
        except ValueError as e:
            print(f"Warning: {e}")
            continue
            
    return pd.DataFrame(data)

def plot_expression_analysis(df: pd.DataFrame, output_dir: str):
    """Create comprehensive plots showing gene expression patterns."""
    
    # Set consistent order for cell types matching developmental progression
    cell_type_order = ['Radial glia', 'Glioblasts', 'OPCs', 'Neuroblasts', 'Neurons']
    
    # 1. Expression trajectory across ages for each cell type
    cell_types = df['MainType'].unique()
    
    for cell_type in cell_types:
        plt.figure(figsize=(15, 8))
        cell_type_df = df[df['MainType'] == cell_type]
        
        for gene in cell_type_df['Gene'].unique():
            gene_df = cell_type_df[cell_type_df['Gene'] == gene]
            sns.lineplot(data=gene_df, x='Age', y='Expression', 
                        label=gene, marker='o')
    
        plt.title(f'Gene Expression Trajectories in {cell_type}')
        plt.xlabel('Embryonic Day')
        plt.ylabel('Expression Level')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Genes')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/expression_trajectories_{cell_type}.pdf', 
                   bbox_inches='tight')
        plt.close()
    
    # 2. Heatmap for each gene
    for gene in df['Gene'].unique():
        gene_df = df[df['Gene'] == gene]
        
        # Create pivot table for age vs cell type
        pivot_df = gene_df.pivot_table(
            values='Expression',
            index='MainType',
            columns='Age',
            aggfunc='mean'
        ).reindex(cell_type_order)
        
        plt.figure(figsize=(15, 6))
        sns.heatmap(pivot_df, cmap='viridis', annot=True, fmt='.2f',
                    cbar_kws={'label': 'Mean Expression'})
        plt.title(f'{gene} Expression by Cell Type and Age')
        plt.xlabel('Embryonic Day')
        plt.ylabel('Cell Type')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_heatmap.pdf', bbox_inches='tight')
        plt.close()
    
    # 3. Violin plots showing distribution
    plt.figure(figsize=(15, 8))
    for gene in df['Gene'].unique():
        gene_df = df[df['Gene'] == gene]
        
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=gene_df, x='MainType', y='Expression', order=cell_type_order)
        plt.title(f'{gene} Expression Distribution by Cell Type')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{gene}_distribution.pdf', bbox_inches='tight')
        plt.close()

# Define genes to analyze
genes = [
    'Srrm3', 
    'Srrm4',
    # Key markers matching Figure 6e
    'Pax6',       # Radial glia marker
    'Sox2',       # Neural progenitor marker
    'Neurog2',    # Neuroblast marker
    'Eomes',      # Intermediate progenitor marker
    'Neurod2',    # Neuronal differentiation marker
    'Tubb3',      # Neuron marker
    'Olig2',      # OPC marker
    'Sox10'       # Glial lineage marker
]

# Main execution
loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.agg.loom"
plots_dir = 'plots_trajectories'

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Process data and create visualizations
with loompy.connect(loom_path) as ds:
    df = extract_cell_data(ds, genes)

# Create plots
plot_expression_analysis(df, plots_dir)

# Generate detailed statistics
stats_df = df.groupby(['Gene', 'MainType', 'Age'])['Expression'].agg([
    'mean', 'std', 'count',
    ('percent_expressing', lambda x: (x > 0).mean() * 100)
]).round(3)

# Save results
df.to_csv(f'{plots_dir}/full_expression_data.csv', index=False)
stats_df.to_csv(f'{plots_dir}/expression_statistics.csv')

print("\nAnalysis complete. Files saved to:", plots_dir)
gc.collect()