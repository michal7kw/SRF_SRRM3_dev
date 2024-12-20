# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python (snakemake)
#     language: python
#     name: snakemake
# ---

# # Environment setup

# +
# Import required libraries
import loompy
import numpy as np
import pandas as pd
import gc  # for garbage collection
import os
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List
import warnings
warnings.filterwarnings('ignore')

# Set working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain')

# Set plotting style
plt.style.use('default')  # Use default matplotlib style
sns.set_theme()  # Apply seaborn styling
sns.set_palette("husl")

# Create plots directory if it doesn't exist
plots_dir = 'plots_clusters'
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)
# -

loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.agg.loom"

# Verify file exists
print(f"File exists: {os.path.exists(loom_path)}")
print(f"File size: {os.path.getsize(loom_path) / (1024**3):.2f} GB")

# # Functions definitions

# Function to get age from column name
def extract_age(age_col):
    return float(age_col.replace('Age_e', '').replace('_0', '').replace('_', '.'))


def extract_expression_data(loom_path: str, gene: str = 'Srrm3', tissues: List[str] = ['Tissue_ForebrainDorsal']) -> List[Dict]:
    """
    Extract expression data for a given gene from loom file, focusing on specific tissue regions.
    
    Args:
        loom_path: Path to the loom file
        gene: Gene name to analyze (default: 'Srrm3')
        tissues: List of tissue column names to filter on (default: ['Tissue_ForebrainDorsal'])
        
    Returns:
        List of dictionaries containing age and expression data for each cell
    """
    # Connect to the loom file
    with loompy.connect(loom_path) as ds:
        # Find gene index
        gene_idx = np.where(ds.ra['Gene'] == gene)[0]
        if len(gene_idx) == 0:
            raise ValueError(f"{gene} gene not found in the dataset")
        
        # Get gene expression data
        gene_expression = ds[gene_idx[0], :]
        
        # Get age columns
        age_cols = [col for col in ds.ca.keys() if col.startswith('Age_e')]
        
        # Create a dictionary to store age and expression data
        age_expression_data = []
        
        # For each cell
        for cell_idx in range(ds.shape[1]):
            # Check if cell belongs to any of the specified tissues
            if any(ds.ca[tissue][cell_idx] == 1 for tissue in tissues):
                # Find which age this cell belongs to
                cell_age = None
                for age_col in age_cols:
                    if ds.ca[age_col][cell_idx] == 1:
                        cell_age = extract_age(age_col)
                        break
                
                if cell_age is not None:
                    age_expression_data.append({
                        'Age': cell_age,
                        'Expression': gene_expression[cell_idx],
                        'Cluster': ds.ca['ClusterName'][cell_idx],
                        'Class': ds.ca['Class'][cell_idx],
                        'Subclass': ds.ca['Subclass'][cell_idx],
                        'Tissue': [tissue for tissue in tissues if ds.ca[tissue][cell_idx] == 1][0]
                    })
                    
    return age_expression_data


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


# # Data exploration

examine_loom_file(loom_path)

# **Available tissues:**
#
# - Tissue_All: int64
# - Tissue_Forebrain: int64
# - Tissue_ForebrainDorsal: int64
# - Tissue_ForebrainVentral: int64
# - Tissue_ForebrainVentroLateral: int64
# - Tissue_ForebrainVentroThalamic: int64
# - Tissue_Head: int64
# - Tissue_Hindbrain: int64
# - Tissue_Midbrain: int64
# - Tissue_MidbrainDorsal: int64
# - Tissue_MidbrainVentral: int64

# ## Expression format

ds = loompy.connect(loom_path)
matrix = ds[:, :]
print("Matrix shape:", matrix.shape)
print("Value range:", np.min(matrix), "-", np.max(matrix))
print("Zero proportion:", np.sum(matrix == 0) / matrix.size)
print("Mean expression:", np.mean(matrix))
print("Mean expression (excluding zeros):", np.mean(matrix[matrix > 0]))


plt.figure(figsize=(8,4))
sns.histplot(matrix[matrix > 0].flatten(), bins=50, log=True)
plt.xlabel('Expression values')
plt.ylabel('Count')

# +
# Get Srrm3 index
srrm3_idx = np.where(ds.ra.Gene == 'Srrm3')[0][0]
srrm3_expr = matrix[srrm3_idx, :]

plt.figure(figsize=(8,4))
sns.histplot(srrm3_expr[srrm3_expr > 0], bins=30)
plt.xlabel('Srrm3 expression')
plt.ylabel('Count')

# Calculate mean expression excluding zeros
srrm3_mean_nonzero = np.mean(srrm3_expr[srrm3_expr > 0])
print(f"Mean Srrm3 expression (excluding zeros): {srrm3_mean_nonzero:.3f}")
# -

# Calculate proportion of non-zero expression
srrm3_expressed = np.sum(srrm3_expr > 0) / len(srrm3_expr)
print(f"Srrm3 is expressed in {srrm3_expressed:.2%} of samples")

# +
age_expression_data = extract_expression_data(loom_path, tissues=['Tissue_ForebrainDorsal'])

# Convert to DataFrame
df = pd.DataFrame(age_expression_data)
print(df.shape)
# -

df.head()

print(list(df['Class'].unique()))
print(list(df['Subclass'].unique()))
print(list(df['Cluster'].unique()))
print(list(df['Tissue'].unique()))


print(list(df[df['Class'] == 'Gastrulation']["Subclass"].unique()))
print(list(df[df['Class'] == 'Neuron']["Subclass"].unique()))

neuronal_subclasses = list(df[df['Class'] == 'Neuron']["Subclass"].unique())
# Count occurrences of each neuronal subclass
neuronal_subclass_counts = df[df['Class'] == 'Neuron']['Subclass'].value_counts()
print("\nNeuronal subclass counts:")
print(neuronal_subclass_counts)


# Filter for neuronal lineages (you might need to adjust these filters based on your data)
df['Is_Neuronal'] = df['Class'] == 'Neuron'
df_neurons = df[df['Is_Neuronal']]

df_neurons.head()

# +
# Aggregate expression data by neuronal subclass
pivot_df = df_neurons.groupby('Subclass')['Expression'].mean().reset_index()

# Sort by expression value and set 'Subclass' as the index
pivot_df = pivot_df.sort_values('Expression', ascending=False)
pivot_df.set_index('Subclass', inplace=True)

# Create the plot
plt.figure(figsize=(10, 8))
sns.heatmap(pivot_df, 
            cmap='viridis',  
            annot=True,      # Show values in cells
            fmt='.2f',       # Format to 2 decimal places
            cbar_kws={'label': 'Mean Expression'})

# Customize the plot
plt.title('Mean Srrm3 Expression by Neuronal Subclass', fontsize=14, pad=20)
plt.xlabel('Mean Expression Level', fontsize=12)
plt.ylabel('Neuronal Subclass', fontsize=12)

# Adjust layout
plt.tight_layout()

# Save figures
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_heatmap_no_age.pdf'), 
            bbox_inches='tight', 
            dpi=300)
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_heatmap_no_age.png'), 
            bbox_inches='tight', 
            dpi=300)

plt.show()

# +
# Calculate mean expression by age and class
mean_by_age_class = df_neurons.groupby(['Age', 'Class'])['Expression'].agg(['mean', 'std']).reset_index()

# Create line plot for different neuronal classes
plt.figure(figsize=(12, 6))
for class_name in mean_by_age_class['Class'].unique():
    class_data = mean_by_age_class[mean_by_age_class['Class'] == class_name]
    plt.plot(class_data['Age'], class_data['mean'], marker='o', label=class_name)
    plt.fill_between(class_data['Age'], 
                    class_data['mean'] - class_data['std'],
                    class_data['mean'] + class_data['std'],
                    alpha=0.2)

plt.title('Mean Srrm3 Expression in Cortical Neuronal Lineages')
plt.xlabel('Embryonic Day')
plt.ylabel('Mean Expression Level')
plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_lineplot.pdf'), 
            bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_lineplot.png'), 
            bbox_inches='tight', dpi=300)
plt.show()

# +
# Create heatmap of expression by age and neuronal subclass
pivot_df = df_neurons.pivot_table(
    values='Expression',
    index='Subclass',
    columns='Age',
    aggfunc='mean'
)

plt.figure(figsize=(15, 10))
sns.heatmap(pivot_df, cmap='viridis', annot=True, fmt='.2f', 
            cbar_kws={'label': 'Mean Expression'})
plt.title('Srrm3 Expression by Neuronal Subclass and Age in Cortical Regions')
plt.xlabel('Embryonic Day')
plt.ylabel('Neuronal Subclass')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_heatmap.pdf'), 
            bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_heatmap.png'), 
            bbox_inches='tight', dpi=300)
plt.show()

# +
# Print summary statistics
print("\nSummary Statistics for Srrm3 Expression in Cortical Neuronal Lineages:")
print("\nBy Age and Class:")
summary_stats = df_neurons.groupby(['Age', 'Class'])['Expression'].describe()
print(summary_stats)

# Calculate the percentage of expressing cells in each class
print("\nPercentage of cells expressing Srrm3 (Expression > 0):")
expression_percentage = df_neurons.groupby('Class').agg({
    'Expression': lambda x: (x > 0).mean() * 100
}).round(2)
print(expression_percentage)

# +
plt.figure(figsize=(15, 6))
sns.boxplot(data=df_neurons, x='Subclass', y='Expression')
plt.xticks(rotation=45, ha='right')
plt.title('Srrm3 Expression Distribution by Neuronal Subclass in Cortical Regions')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_boxplot.pdf'), 
            bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_boxplot.png'), 
            bbox_inches='tight', dpi=300)
plt.show()

# Save the summary statistics to a text file
with open(os.path.join(plots_dir, 'srrm3_cortical_statistics.txt'), 'w') as f:
    f.write("Summary Statistics for Srrm3 Expression in Cortical Neuronal Lineages:\n\n")
    f.write("By Age and Class:\n")
    f.write(summary_stats.to_string())
    f.write("\n\nPercentage of cells expressing Srrm3 (Expression > 0):\n")
    f.write(expression_percentage.to_string())

# Force garbage collection
gc.collect()
