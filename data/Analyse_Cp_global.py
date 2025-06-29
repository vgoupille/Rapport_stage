# %%
import sys
import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir("/Users/valentingoupille/Documents/Rapport_stage/")

# %%
# Define path to AnnData objects
path_to_anndata = "data/anndata_objects_with_gene_metadata/"

# Check available files
print("Available AnnData files:")
print(os.listdir(path_to_anndata))

# %%
# Load the original AnnData object
adata_original = ad.read_h5ad(path_to_anndata + "final_adata_filtered.h5ad")
print(f"Original data shape: {adata_original.shape}")

# %%
# Create a working copy
adata = adata_original.copy()

# %%
# View gene metadata
print("Gene metadata columns:")
print(adata.var.columns.tolist())

# %%
# Create a test copy for analysis
test = adata.copy()

# %%
# Store raw data
test.raw = test.copy()
print(f"Raw data shape: {test.raw.X.shape}")

# %%
# View gene filtering results
list_pass = [
    "gene_pass_union_OD_5pct",
    "gene_pass_intersection_OD_5pct",
    "gene_pass_union_OD_10pct",
    "gene_pass_intersection_OD_10pct",
    "gene_pass_union_OD_15pct",
    "gene_pass_intersection_OD_15pct",
]

# Display first few rows
print("First rows of filtering results:")
print(test.var[list_pass].head())

# Display value counts for each threshold
print("\nValue counts for each filtering threshold:")
for col in list_pass:
    print(f"\n{col}:")
    print(test.var[col].value_counts())

# %%
# Filter genes based on gene_pass_union_OD_5pct
adata_union = test[:, test.var["gene_pass_union_OD_5pct"] == True].copy()
print(f"Number of genes after filtering: {adata_union.shape[1]}")

# %%
# Verify filtering results
print("Verification of filtering:")
print(adata_union.var["gene_pass_union_OD_5pct"].value_counts())


# %% filter to only keep Culture Medium M9F
adata_union_m9f = adata_union[adata_union.obs["CultureMedium"] == "M9F"]
print(f"Number of genes after filtering: {adata_union_m9f.shape[1]}")

# %%
# filter to only keep Culture Medium M9
adata_union_m9 = adata_union[adata_union.obs["CultureMedium"] == "M9"]
print(f"Number of genes after filtering: {adata_union_m9.shape[1]}")


# %%
# Create a test copy for analysis
adata_test = adata_union_m9f.copy()

# %%
# Normalize data
sc.pp.normalize_total(adata_test, target_sum=2000, inplace=True)
sc.pp.log1p(adata_test)
print("Data normalized and log-transformed")
sc.pp.scale(adata_test, max_value=5)


# # %%
# # Identify highly variable genes
# sc.pp.highly_variable_genes(adata_test, n_top_genes=200)
# print(
#     f"📊 Nombre de gènes hautement variables: {adata_test.var['highly_variable'].sum():,}"
# )

# # Plot highly variable genes
# sc.pl.highly_variable_genes(adata_test)

# # %%
# # Keep only highly variable genes
# adata_test = adata_test[:, adata_test.var["highly_variable"]]
# print(f"Final data shape after HVG filtering: {adata_test.shape}")

# %%
# Perform PCA
sc.tl.pca(adata_test)
print("PCA computed")

# %%
# Plot PCA variance ratio
sc.pl.pca_variance_ratio(adata_test, n_pcs=50, log=True)

# %%
# Plot PCA colored by different variables
sc.pl.pca(adata_test, color="CultureMedium")
sc.pl.pca(adata_test, color=["OD", "CultureMedium", "RepBio", "total_counts"])

# %%
# Compute neighborhood graph
sc.pp.neighbors(
    adata_test,
    n_neighbors=10,  # Number of neighbors (5-50)
    n_pcs=10,  # Number of PCs to use (10-100)
    use_rep="X_pca",  # Representation to use
    metric="euclidean",  # Distance metric
)
print("Neighborhood graph computed")

# %%
# Compute UMAP
sc.tl.umap(adata_test)
print("UMAP computed")

# %%
# Perform Leiden clustering
sc.tl.leiden(
    adata_test,
    resolution=0.5,  # Resolution (0.1-2.0)
    n_iterations=2,  # Iterations (1-10)
    random_state=42,  # Reproducibility
)
print("Leiden clustering performed")

# %%
# Plot UMAP with different colorings
sc.pl.umap(
    adata_test, color=["leiden", "CultureMedium", "RepBio", "OD", "total_counts"]
)

# %%
