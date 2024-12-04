import numpy as np
import umap
import os

# Get script directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Load data using absolute path
input_file = os.path.join(script_dir, "scanpy_input.tab")
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input file not found: {input_file}")
    
pca = np.loadtxt(input_file)

# Run UMAP
reducer = umap.UMAP(random_state=42)
embedding = reducer.fit_transform(pca)

# Save results with absolute path
output_file = os.path.join(script_dir, "umap.tab")
np.savetxt(output_file, embedding, delimiter="\t")

