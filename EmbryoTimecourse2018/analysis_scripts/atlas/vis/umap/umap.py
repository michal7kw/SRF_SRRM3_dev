import numpy as np
import umap

# Load data
pca = np.loadtxt("scanpy_input.tab")

# Run UMAP
reducer = umap.UMAP(random_state=42)
embedding = reducer.fit_transform(pca)

# Save results
np.savetxt("umap.tab", embedding, delimiter="\t")

