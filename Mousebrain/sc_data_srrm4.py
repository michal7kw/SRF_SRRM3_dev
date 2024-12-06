# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python (snakemake)
#     language: python
#     name: snakemake
# ---

# # Environment

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

# Add local_functions directory to Python path
from local_functions import *

# Set plotting style
plt.style.use('default')  # Use default matplotlib style
sns.set_theme()  # Apply seaborn styling
sns.set_palette("husl")
# -

loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.loom"

# Verify file exists
print(f"File exists: {os.path.exists(loom_path)}")
print(f"File size: {os.path.getsize(loom_path) / (1024**3):.2f} GB")

# Create plots directory if it doesn't exist
os.makedirs('plots_sc', exist_ok=True)

# # Define functions

# **Column attributes**
# - Age: object
# - BTSNE: float64
# - CellCycle: float64
# - CellID: object
# - Cell_Conc: int64
# - Chemistry: object
# - ChipID: object
# - Class: object
# - ClusterName: object
# - Clusters: int64
# - Date_Captured: object
# - DonorID: object
# - DoubletFinderPCA: float64
# - HPF: float64
# - HPF_LogPP: float64
# - HPF_theta: float64
# - IsCycling: uint8
# - Label: object
# - Location_E9_E11: object
# - NCellsCluster: int64
# - NGenes: float64
# - Num_Pooled_Animals: int64
# - PCA: float64
# - PCR_Cycles: int64
# - Plug_Date: object
# - Project: object
# - PseudoAge: float64
# - PseudoTissue: object
# - Region: object
# - SampleID: object
# - SampleName: object
# - Sample_Index: object
# - Sex: object
# - Species: object
# - Split: int64
# - Strain: object
# - Subclass: object
# - TSNE: float64
# - Target_Num_Cells: float64
# - Tissue: object
# - TotalUMI: float64
# - Transcriptome: object
# - UMAP: float32
# - UMAP3D: float32
# - cDNA_Lib_Ok: object
# - ngperul_cDNA: object


# # Data exploration

examine_loom_file(loom_path)

data = get_data(loom_path, 'Srrm4')

data_df = pd.DataFrame(data)

def create_pseudo_cells(df, group_by_column):
    """
    Create pseudo-cells by aggregating data grouped by a specific column.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe containing single-cell data
    group_by_column : str
        Column name to group the cells by
        
    Returns:
    --------
    pandas.DataFrame
        Aggregated pseudo-cell data with mean/median values for numeric columns
    """
    # Numeric columns to aggregate
    numeric_cols = ['PseudoAge', 'Expression', 'TotalUMI', 'NGenes']
    
    # Initialize aggregation dictionary with numeric columns
    agg_dict = {
        'PseudoAge': 'mean',  # Average pseudo age
        'Expression': 'mean',  # Average expression
        'TotalUMI': 'mean',  # Average UMI count
        'NGenes': 'mean',  # Average number of genes
        'CellID': 'count'  # Count of cells in each group
    }
    
    # Add categorical columns to aggregation dict, excluding the grouping column
    categorical_cols = ['Tissue', 'Sex', 'Subclass', 'Region', 'Age', 'Class']
    for col in categorical_cols:
        if col in df.columns and col != group_by_column:
            agg_dict[col] = lambda x: x.mode().iloc[0]
    
    # Group by the specified column and aggregate
    pseudo_cells = df.groupby(group_by_column).agg(agg_dict).reset_index()
    
    # Rename the CellID column to CellCount
    pseudo_cells = pseudo_cells.rename(columns={'CellID': 'CellCount'})
    
    # Round numeric columns to 2 decimal places
    for col in numeric_cols:
        if col in pseudo_cells.columns:
            pseudo_cells[col] = pseudo_cells[col].round(2)
            
    return pseudo_cells

data_df.head()

# Print unique values for each column
print("Age:", list(data_df.Age.unique()))
print("Tissue:", list(data_df.Tissue.unique()))
print("Sex:", list(data_df.Sex.unique()))
print("Subclass:", list(data_df.Subclass.unique()))
print("Class:", list(data_df.Class.unique()))
print("Region:", list(data_df.Region.unique()))
print("PseudoAge:", list(data_df.PseudoAge.unique()))
print("ClusterName:", list(data_df.ClusterName.unique()))

plt.figure(figsize=(10, 6))
plt.hist(data_df.PseudoAge, bins=30, edgecolor='black')
plt.xlabel('PseudoAge')
plt.ylabel('Frequency')
plt.title('Distribution of PseudoAge')
plt.show()

# +
# Create 10 bins for PseudoAge
min_age = data_df.PseudoAge.min()
max_age = data_df.PseudoAge.max()
bin_size = (max_age - min_age) / 10

# Create bins
bins = [min_age + i * bin_size for i in range(11)]

# Create new column with binned PseudoAge
data_df['PseudoAge_Binned'] = pd.cut(data_df.PseudoAge, bins=bins, labels=list(range(1,11)))
# -

plt.figure(figsize=(10, 6))
plt.hist(data_df.PseudoAge_Binned.astype(str), bins=10, edgecolor='black')
plt.xlabel('PseudoAge')
plt.ylabel('Frequency')
plt.title('Distribution of PseudoAge')
plt.show()

# # Data subsets

type(data_df.Subclass)

data_df_forebrain_dorsal = data_df[data_df.Tissue == 'ForebrainDorsal']
print(data_df_forebrain_dorsal.shape)

# +
print(list(data_df_forebrain_dorsal.Subclass.unique()))

print(list(data_df_forebrain_dorsal.Class.unique()))
# -


data_df_forebrain_dorsal[data_df_forebrain_dorsal.Class.isin(['Gastrulation', 'Neural tube', 'Neural crest', 'Neuroblast', 'Neuron'])]['Class'].value_counts()

data_df_forebrain_dorsal_neurons = data_df_forebrain_dorsal[data_df_forebrain_dorsal.Class.isin(['Neuron'])]

# # Plots

# ## Cell type expression

plot_cell_type_expression(data_df_forebrain_dorsal[data_df_forebrain_dorsal.Class.isin(['Neuroblast', 'Neuron'])])  

pseudo_cells_by_class = create_pseudo_cells(data_df_forebrain_dorsal[data_df_forebrain_dorsal.Class.isin(['Neuroblast', 'Neuron'])], 'Class')
pseudo_cells_by_class.head()

# ## Temporal expression

print("Age value counts:")
print(data_df_forebrain_dorsal_neurons.Age.value_counts())
print("\nPseudoAge_Binned value counts:")
print(data_df_forebrain_dorsal_neurons.PseudoAge_Binned.value_counts())


plot_temporal_expression(data_df_forebrain_dorsal_neurons, "Age")

pseudo_cells_by_age = create_pseudo_cells(data_df_forebrain_dorsal_neurons, 'Age')
pseudo_cells_by_age.head()

plot_temporal_cell_type_heatmap(pseudo_cells_by_age,"Age")

plot_regional_expression(data_df)

pseudo_cells_by_region = create_pseudo_cells(data_df, 'Region')
pseudo_cells_by_region.head()

plot_regional_expression(pseudo_cells_by_region[pseudo_cells_by_region.Class.isin(['Neuron'])])

plot_temporal_cell_type_heatmap(data_df, "Age")
plot_temporal_cell_type_heatmap(data_df, "PseudoAge")
plot_temporal_cell_type_heatmap(data_df, "PseudoAge_Binned")

# # Statistics

# Basic statistics
print("\nSrrm4 Expression Statistics:")
print(f"Mean expression: {np.mean(data['Expression']):.4f}")
print(f"Median expression: {np.median(data['Expression']):.4f}")
print(f"Percentage of cells expressing Srrm4: {(data['Expression'] > 0).mean() * 100:.2f}%")

# +
# Top expressing cell types
df = pd.DataFrame({
    'Cell Type': data['Class'],
    'Expression': data['Expression']
})

print("\nTop 5 Cell Types by Mean Srrm4 Expression:")
print(df.groupby('Cell Type')['Expression']
        .agg(['mean', 'count'])
        .sort_values('mean', ascending=False)
        .head())


