import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# Read h5ad data
adata = sc.read("sceN1.h5ad")

# Ensure celltype is category dtype
adata.obs['celltype'] = adata.obs['celltype'].astype('category')

# Print basic information
print(adata.obs['celltype'].head())
print(adata.obs.dtypes)

#######################################
# DiffusionMap calculation
#######################################
# Compute neighbors
sc.pp.neighbors(adata, n_pcs=20)

# Compute Diffusion Map
sc.tl.diffmap(adata)
print(adata.obsm['X_diffmap'])

# 2D visualization
sc.pl.diffmap(adata)
sc.pl.diffmap(adata, color='celltype', use_raw=False)

# 3D visualization
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(
    adata.obsm['X_diffmap'][:, 0],
    adata.obsm['X_diffmap'][:, 1],
    adata.obsm['X_diffmap'][:, 2]
)
ax.set_xlabel('DC1')
ax.set_ylabel('DC2')
ax.set_zlabel('DC3')
plt.show()

#######################################
# PAGA trajectory inference
#######################################
sc.tl.paga(adata, groups='celltype')
sc.pl.paga(adata)

#######################################
# Subset selected cell types and set root
#######################################
subset = adata[adata.obs['celltype'].isin(['monocyte', 'macrophage_like'])].copy()
root_cell = subset.obs.query('celltype == "monocyte"').index[0]
subset.uns['iroot'] = subset.obs_names.get_loc(root_cell)

#######################################
# Compute Diffusion Pseudotime (DPT)
#######################################
sc.tl.dpt(subset)

print(subset.obs['dpt_pseudotime'])
print(subset.obs[['celltype', 'dpt_pseudotime']].head())
print(subset.obs[['celltype', 'dpt_pseudotime']].tail())

# UMAP visualization of DPT
sc.pl.umap(subset, color='dpt_pseudotime')

#######################################
# Density plot for selected cell types
#######################################
selected_cells = subset.obs[subset.obs['celltype'].isin(['monocyte', 'macrophage_like'])]
g = sns.FacetGrid(
    selected_cells, row="celltype", aspect=4, height=2,
    sharex=True, sharey=False, row_order=['monocyte', 'macrophage_like']
)
g.map(sns.kdeplot, "dpt_pseudotime", fill=True, color="skyblue", alpha=0.7)
g.set_titles("{row_name}")
g.set_axis_labels("DPT", "Density")
g.fig.subplots_adjust(hspace=0.5)
plt.show()

#######################################
# Export results to CSV (for R)
#######################################
selected_cells[['dpt_pseudotime', 'celltype']].to_csv("selected_cells.csv", index=False)
