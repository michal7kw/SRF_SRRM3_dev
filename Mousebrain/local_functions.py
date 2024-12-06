import loompy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from typing import Dict

# %%
def get_data(loom_path: str, gene_name: str) -> Dict:
    """Extract data from the loom file
    
    Args:
        loom_path: Path to the .loom file containing single-cell data
        gene_name: Name of the gene to extract
    Returns:
        Dictionary containing expression data and relevant metadata
    """
    with loompy.connect(loom_path) as ds:
            
        # Find gene index
        gene_indices = np.where(ds.ra.Gene == gene_name)[0]
        if len(gene_indices) == 0:
            raise ValueError(f"{gene_name} gene not found in dataset")
        gene_idx = gene_indices[0]
        
        # Extract expression data
        gene_expr = ds[gene_idx, :].astype(np.float32)
        
        # Extract relevant metadata
        metadata = {
            'Age': np.array(ds.ca.Age),
            'Tissue': np.array(ds.ca.Tissue), 
            'Sex': np.array(ds.ca.Sex),
            'Subclass': np.array(ds.ca.Subclass),
            'Class': np.array(ds.ca.Class),
            'Region': np.array(ds.ca.Region),
            'PseudoAge': np.array(ds.ca.PseudoAge),
            'Expression': gene_expr,
            'CellID': np.array(ds.ca.CellID),
            'ClusterName': np.array(ds.ca.ClusterName),
            'TotalUMI': np.array(ds.ca.TotalUMI),
            'NGenes': np.array(ds.ca.NGenes)
        }
        
        return metadata

# %%
def plot_temporal_cell_type_heatmap(data: Dict, time_column: str):
    """Plot heatmap of expression across cell types and time"""
    df = pd.DataFrame({
        time_column: data[time_column],
        'Cell Type': data['Class'],
        'Expression': data['Expression']
    })
    
    # Calculate mean expression for each cell type at each time point
    pivot_table = df.pivot_table(
        values='Expression',
        index='Cell Type',
        columns=time_column,
        aggfunc='mean'
    )
    
    # Get min and max values from the data
    vmin = pivot_table.min().min()
    vmax = pivot_table.max().max()
    
    plt.figure(figsize=(15, 10))
    sns.heatmap(pivot_table, cmap='viridis', center=vmin+(vmax-vmin)/2, vmin=vmin, vmax=vmax, cbar_kws={'label': 'Expression'})
    plt.title(f'Expression Across Cell Types and {time_column}')
    plt.xlabel(time_column)
    plt.ylabel('Cell Type')
    plt.tight_layout()
    # plt.savefig(f'plots/expression_temporal_cell_type_heatmap_{time_column}.pdf')
    # plt.close()

# %%
def plot_cell_type_expression(data: Dict):
    """Plot expression across cell types"""
    plt.figure(figsize=(15, 6))
    
    # Create box plot
    df = pd.DataFrame({
        'Cell Type': data['Class'],
        'Expression': data['Expression']
    })
    
    # Calculate mean expression per cell type for sorting
    cell_type_means = df.groupby('Cell Type')['Expression'].mean().sort_values(ascending=False)
    
    sns.boxplot(data=df, x='Cell Type', y='Expression', order=cell_type_means.index)
    plt.title('Expression Across Cell Types')
    plt.xlabel('Cell Type')
    plt.ylabel('Expression Level')
    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.savefig('plots_sc/expression_cell_type_expression.pdf')
    # plt.close()

# %%
def plot_temporal_expression(data: Dict, time_column: str):
    """Plot expression across developmental time"""
    plt.figure(figsize=(12, 6))
    
    # Create violin plot
    df = pd.DataFrame({
        'Age': data[time_column],
        'Expression': data['Expression']
    })
    
    # Sort the age values to ensure proper ordering on x-axis
    age_order = sorted(df['Age'].unique())
    
    sns.boxplot(data=df, x='Age', y='Expression', order=age_order)
    plt.title('Expression Across Developmental Time')
    plt.xlabel('Pseudo Age')
    plt.ylabel('Expression Level')
    plt.xticks(rotation=45)
    plt.tight_layout()
    # plt.savefig('plots_sc/expression_temporal_expression.pdf')
    # plt.close()

# %%
def examine_loom_file(loom_path: str) -> None:
    """
    Examine and print basic information about a loom file.
    
    Args:
        loom_path: Path to the loom file
    """
    # Connect to the loom file (this doesn't load it entirely into memory)
    with loompy.connect(loom_path) as ds:
        # Basic information about the dataset
        print("\nDataset shape:", ds.shape)
        print(f"Number of cells: {ds.shape[1]:,}")
        print(f"Number of genes: {ds.shape[0]:,}")
        
        # Examine column attributes (cell metadata)
        print("\nColumn attributes (cell metadata):")
        for attr in ds.ca.keys():
            print(f"- {attr}: {ds.ca[attr].dtype}")
        
        # Examine row attributes (gene metadata)
        print("\nRow attributes (gene metadata):")
        for attr in ds.ra.keys():
            print(f"- {attr}: {ds.ra[attr].dtype}")
        
        # Get a small sample of the expression matrix (first 5 genes, first 5 cells)
        print("\nSample of expression matrix (5x5):")
        sample_matrix = ds[:5, :5]
        print(sample_matrix)
        
        # Get some basic statistics
        print("\nBasic statistics:")
        print(f"Mean expression: {np.mean(sample_matrix):.4f}")
        print(f"Median expression: {np.median(sample_matrix):.4f}")
        print(f"Sparsity: {(sample_matrix == 0).sum() / sample_matrix.size:.2%}")

    # Force garbage collection
    gc.collect()

def plot_regional_expression(data: Dict):
    """Plot expression across brain regions"""
    plt.figure(figsize=(15, 6))
    
    # Create box plot
    df = pd.DataFrame({
        'Region': data['Region'],
        'Expression': data['Expression']
    })
    
    # Calculate mean expression per region for sorting
    region_means = df.groupby('Region')['Expression'].mean().sort_values(ascending=False)
    
    sns.boxplot(data=df, x='Region', y='Expression', order=region_means.index)
    plt.title('Expression Across Brain Regions')
    plt.xlabel('Brain Region')
    plt.ylabel('Expression Level')
    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.savefig('plots_sc/expression_regional_expression.pdf')
    # plt.close()