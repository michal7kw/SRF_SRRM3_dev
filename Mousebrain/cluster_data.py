# %%
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

# %%
loom_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/dev_all.agg.loom"

# %%
# Verify file exists
print(f"File exists: {os.path.exists(loom_path)}")
print(f"File size: {os.path.getsize(loom_path) / (1024**3):.2f} GB")

# %%
# Function to get age from column name
def extract_age(age_col):
    return float(age_col.replace('Age_e', '').replace('_0', '').replace('_', '.'))

# Connect to the loom file
with loompy.connect(loom_path) as ds:
    # Find Srrm3 gene index
    srrm3_idx = np.where(ds.ra['Gene'] == 'Srrm3')[0]
    if len(srrm3_idx) == 0:
        raise ValueError("Srrm3 gene not found in the dataset")
    
    # Get Srrm3 expression data
    srrm3_expression = ds[srrm3_idx[0], :]
    
    # Get age columns
    age_cols = [col for col in ds.ca.keys() if col.startswith('Age_e')]
    
    # Create a dictionary to store age and expression data
    age_expression_data = []
    
    # For each cell
    for cell_idx in range(ds.shape[1]):
        # Only include cells from cortical regions (Forebrain Dorsal)
        if ds.ca['Tissue_ForebrainDorsal'][cell_idx] == 1:
            # Find which age this cell belongs to
            cell_age = None
            for age_col in age_cols:
                if ds.ca[age_col][cell_idx] == 1:
                    cell_age = extract_age(age_col)
                    break
            
            if cell_age is not None:
                age_expression_data.append({
                    'Age': cell_age,
                    'Expression': srrm3_expression[cell_idx],
                    'Cluster': ds.ca['ClusterName'][cell_idx],
                    'Class': ds.ca['Class'][cell_idx],
                    'Subclass': ds.ca['Subclass'][cell_idx]
                })

# Convert to DataFrame
df = pd.DataFrame(age_expression_data)

# Filter for neuronal lineages (you might need to adjust these filters based on your data)
neuronal_keywords = ['neuron', 'Neural', 'Neuroblast', 'RG', 'IP', 'Projection']
df['Is_Neuronal'] = df['Class'].str.contains('|'.join(neuronal_keywords), case=False)
df_neurons = df[df['Is_Neuronal']]

# %%
# Create violin plot of expression by age for neuronal lineages
plt.figure(figsize=(12, 6))
sns.violinplot(data=df_neurons, x='Age', y='Expression', hue='Class')
plt.title('Srrm3 Expression in Cortical Neuronal Lineages Across Development')
plt.xlabel('Embryonic Day')
plt.ylabel('Expression Level')
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_violin.pdf'), 
            bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(plots_dir, 'srrm3_cortical_violin.png'), 
            bbox_inches='tight', dpi=300)
plt.show()

# %%
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

# %%
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

# %%
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

# %%
# Optional: Box plot showing expression distribution by neuronal subclass
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

# %%