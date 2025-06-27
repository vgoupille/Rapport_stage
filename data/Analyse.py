# %%
import sys
import os
import pprint as pp
import numpy as np

os.chdir("/Users/valentingoupille/Documents/Rapport_stage/")


# %%

import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData


# %% mettre le chemin dans une variable
path_to_anndata = "data/anndata_objects_with_gene_metadata/"
# %% see the dir with the anndata objects
os.listdir(path_to_anndata)
# %%
adata_original = ad.read_h5ad(path_to_anndata + "final_adata_filtered.h5ad")


# %%
adata = adata_original.copy()


# %%
adata.var
# %%
test = adata.copy()

# %%
test.raw = test.copy()

# %%
test.raw.X.shape

# %%
test.raw.X


# # %%
# # 1. Normalisation par cellule
# sc.pp.normalize_total(test, target_sum=1e4, inplace=True)
# print(f"   ✅ Normalisation terminée")
# # view the normalized data
# test.X

# # Vérification après normalisation
# print(f"   Valeurs min/max: {test.X.min():.2f} / {test.X.max():.2f}")
# print(f"   Valeurs non-nulles: {(test.X > 0).sum()}")
# test.norm = test.copy()
# test.norm.X.shape
# test.norm.X

# # %%


# # %%
# test.obs


# # %%

# # 2. Transformation log permet de moins avoir des valeurs trop grandes
# sc.pp.log1p(test)
# print(f"   ✅ Transformation log terminée")

# print(f"   Valeurs min/max: {test.X.min():.2f} / {test.X.max():.2f}")
# print(f"   Valeurs non-nulles: {(test.X > 0).sum()}")
# test.log = test.copy()
# test.log.X.shape
# test.log.X


# # %%
# # 3. Scaling (optionnel, mais recommandé pour PCA)
# sc.pp.scale(test, max_value=10)
# print(f"   ✅ Scaling terminé")
# print(f"   Valeurs min/max: {test.X.min():.2f} / {test.X.max():.2f}")
# print(f"   Valeurs non-nulles: {(test.X > 0).sum()}")
# test.scale = test.copy()
# test.scale.X.shape
# test.scale.X


# # %%
# test
# test.var
# # %%
# # %%
# # # Sélection des gènes hautement variables
# # Problème avec seurat_v3 - utiliser une méthode plus robuste
# print(f"🔍 Tentative de sélection des gènes hautement variables...")

# # Option 1: Méthode simple et robuste (recommandée pour les petits datasets)
# print(f"🔄 Utilisation de la méthode simple...")
# sc.pp.highly_variable_genes(test, flavor="seurat", n_top_genes=100)
# print(f"✅ Méthode seurat réussie")

# # Option 2: Si vous voulez essayer seurat_v3 (peut échouer)
# # try:
# #     sc.pp.highly_variable_genes(test, flavor="seurat_v3")
# #     print(f"✅ Méthode seurat_v3 réussie")
# # except Exception as e:
# #     print(f"❌ Seurat v3 échoué: {e}")
# #     # Fallback vers seurat classique
# #     sc.pp.highly_variable_genes(test, flavor="seurat", n_top_genes=2000)
# #     print(f"✅ Méthode seurat classique réussie")

# # %%
# # # Vérifier les gènes hautement variables
# print(f"📊 Nombre de gènes hautement variables: {test.var['highly_variable'].sum():,}")
# # %%
# # %%
# # # Plot des gènes hautement variables
# sc.pl.highly_variable_genes(test)

# # %%
# # %%
# # PCA
# sc.tl.pca(test)

# # %%
# test.var


# # %%
# # Plot PCA
# sc.pl.pca(test, color="CultureMedium")
# sc.pl.pca(test, color="ODt")
# sc.pl.pca(test, color="Group")
# sc.pl.pca(test, color="OD")
# sc.pl.pca(test, color="total_counts")

# # %%
# test.obs

# # %%
# test.var

# # # %%
# # sc.pp.calculate_qc_metrics(
# #     test, inplace=True
# # )

# test.obs
# # %%
# test.var

# # %%
# test.obs


# # %%
# # Plot variance expliquée par les composantes principales
# sc.pl.pca_variance_ratio(test, n_pcs=50, log=True)

# # %%
# # Calcul des voisins
# sc.pp.neighbors(
#     test,
#     n_neighbors=10,
# )


# # %%
# # UMAP
# sc.tl.umap(test, min_dist=0.5)

# # %%
# # Plot UMAP coloré par groupe
# sc.pl.umap(test, color="Group", size=15)
# sc.pl.umap(test, color="OD", size=15)
# sc.pl.umap(test, color="total_counts", size=15)

# # %%
# # Plot UMAP coloré par groupe

# sc.pl.umap(test, color="CultureMedium", size=15)
# sc.pl.umap(test, color="RepBio", size=15)
# sc.pl.umap(test, color="OD", size=15)
# sc.pl.umap(test, color="total_counts", size=15)

# # %%
# # Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
# sc.tl.leiden(test, flavor="igraph", n_iterations=2)

# # %%
# sc.pl.umap(test, color="leiden", size=10)

# # %%
# sc.pl.umap(test, color="leiden", size=10)


#  %%
adata.obs
# %%
adata.var
# %%
adata.var["gene_pass_union_OD_5pct"].value_counts()
# %% Filter genes based on gene_pass_union_OD_5pct
adata_union = adata[:, adata.var["gene_pass_union_OD_5pct"] == True].copy()
print(f"Number of genes after filtering: {adata_union.shape[1]}")

# %%
adata_intersect = adata[:, adata.var["gene_pass_intersection_OD_5pct"] == True].copy()
print(f"Number of genes after filtering: {adata_intersect.shape[1]}")

# %%
adata_intersect.var["gene_pass_intersection_OD_5pct"].value_counts()


# %%
sc.pp.normalize_total(adata_intersect, target_sum=1e4, inplace=True)
sc.pp.log1p(adata_intersect)
sc.pp.scale(adata_intersect)
sc.tl.pca(adata_intersect)

# %%
adata_intersect.obsm["X_pca"]
# %%
sc.pl.pca_variance_ratio(adata_intersect, n_pcs=50, log=True)

# %%
sc.pl.pca(adata_intersect, color="Group")
sc.pl.pca(adata_intersect, color="CultureMedium")
sc.pl.pca(adata_intersect, color="total_counts")
# %%

# %%
# Try different numbers of neighbors to find optimal clustering
for n in [5, 10, 15, 20, 30]:
    print(f"\nTesting with {n} neighbors:")
    sc.pp.neighbors(adata_intersect, n_neighbors=n)
    # Calculate UMAP and Leiden clustering to evaluate results
    sc.tl.umap(adata_intersect)
    sc.tl.leiden(adata_intersect)
    print(f"Number of clusters found: {len(adata_intersect.obs['leiden'].unique())}")
# %%
sc.pp.neighbors(
    adata_intersect,
    n_neighbors=5,  # Nombre de voisins (5-50)
    n_pcs=10,  # Nombre de PCs à utiliser (10-100)
    use_rep="X_pca",  # Représentation à utiliser
    metric="euclidean",  # Métrique de distance
)
# %%
sc.tl.umap(
    adata_intersect,
    min_dist=0.5,  # Distance minimale (0.1-1.0)
    spread=1.0,  # Spread (0.5-2.0)
    n_components=10,  # Dimensions (2 ou 3)
    random_state=42,  # Reproductibilité
)
# %%
sc.tl.leiden(
    adata_intersect,
    resolution=0.9,  # Résolution (0.1-2.0)
    n_iterations=2,  # Itérations (1-10)
    random_state=42,  # Reproductibilité
)
# %%
sc.pl.umap(adata_intersect, color="leiden")
# %%
sc.pl.umap(adata_intersect, color="Group")
# %%
sc.pl.umap(adata_intersect, color="CultureMedium")
# %%
adata_intersect.obs.columns
# %%
sc.pl.umap(adata_intersect, color="OD_measured")
# %%
sc.pl.umap(adata_intersect, color="total_counts")
# %%
sc.pl.umap(adata_intersect, color="n_genes_by_counts")


# %%
# Test de différents paramètres
results = []

# Test différentes combinaisons
for n_neighbors in [5, 10, 15, 20]:
    for n_pcs in [20, 30, 50]:
        for resolution in [0.5, 1.0, 1.5]:

            # Recalculer neighbors
            sc.pp.neighbors(adata_intersect, n_neighbors=n_neighbors, n_pcs=n_pcs)

            # Recalculer UMAP
            sc.tl.umap(adata_intersect, min_dist=0.5)

            # Recalculer clustering
            sc.tl.leiden(adata_intersect, resolution=resolution)

            # Évaluer les résultats
            n_clusters = len(adata_intersect.obs["leiden"].unique())

            results.append(
                {
                    "n_neighbors": n_neighbors,
                    "n_pcs": n_pcs,
                    "resolution": resolution,
                    "n_clusters": n_clusters,
                }
            )

            print(
                f"Neighbors: {n_neighbors}, PCs: {n_pcs}, Resolution: {resolution} -> {n_clusters} clusters"
            )

# %%
adata_intersect.var

# Remplacer les noms de gènes dans l'objet AnnData
adata_intersect.var["original_gene_ids"] = (
    adata_intersect.var_names
)  # Optionnel : garder les anciens noms
adata_intersect.var_names = adata_intersect.var["final_gene_name"]

print("✅ Les noms de gènes ont été remplacés avec succès par 'final_gene_name'.")
# %%
# Plot genes expression levels across leiden clusters
sc.pl.dotplot(
    adata_intersect,
    var_names=adata_intersect.var_names,  # Use all genes
    groupby="leiden",
    standard_scale="var",  # Scale variance across genes
    dendrogram=True,  # Add hierarchical clustering
)
# %%
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata_intersect, groupby="leiden", method="wilcoxon")
# %%
sc.pl.rank_genes_groups(adata_intersect, n_genes=10)
# %%
sc.pl.rank_genes_groups_dotplot(adata_intersect)
# %%
sc.pl.rank_genes_groups_heatmap(adata_intersect)
# %%
sc.pl.rank_genes_groups_matrixplot(adata_intersect)
# %%
sc.get.rank_genes_groups_df(adata_intersect, group="2").head(5)
# %%
