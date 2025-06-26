# %%
import sys
import os
import pprint as pp
import numpy as np

os.chdir("/Users/valentingoupille/Documents/Rapport_stage/")
# %%


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
adata_original = ad.read_h5ad(path_to_anndata + "GeneFull_raw_unique_mult.h5ad")

# %%
adata = adata_original.copy()

# %% Load annotation file and update Chromosome column
print("Loading annotation file...")
annotation_df = pd.read_csv("data/extracted_features_no_duplicates.tsv", sep="\t")
print(
    f"Annotation file loaded with {len(annotation_df)} rows and columns: {annotation_df.columns.tolist()}"
)

# Check if final_gene_name column exists in annotation
if "final_gene_name" in annotation_df.columns:
    print("'final_gene_name' column found in annotation file")

    # Create a mapping from final_gene_name to seq_type
    gene_to_seqtype = dict(
        zip(annotation_df["final_gene_name"], annotation_df["seq_type"])
    )

    # Check how many genes in adata.var_names are in the annotation
    genes_in_annotation = [gene for gene in adata.var_names if gene in gene_to_seqtype]
    print(
        f"Number of genes in adata that are in annotation: {len(genes_in_annotation)} out of {len(adata.var_names)}"
    )

    # Update the Chromosome column in adata.var
    if "Chromosome" in adata.var.columns:
        print("Updating existing Chromosome column...")
        # Create a new chromosome mapping
        new_chromosome_values = []
        for gene in adata.var_names:
            if gene in gene_to_seqtype:
                new_chromosome_values.append(gene_to_seqtype[gene])
            else:
                new_chromosome_values.append("unknown")

        adata.var["Chromosome"] = new_chromosome_values
        print("Chromosome column updated successfully!")

        # Show unique values in the updated Chromosome column
        unique_chromosomes = adata.var["Chromosome"].unique()
        print(f"Unique values in updated Chromosome column: {unique_chromosomes}")

        # Count genes per chromosome
        chromosome_counts = adata.var["Chromosome"].value_counts()
        print("Number of genes per chromosome:")
        for chrom, count in chromosome_counts.items():
            print(f"  {chrom}: {count} genes")
    else:
        print("Chromosome column not found in adata.var, creating new column...")
        # Create new chromosome mapping
        new_chromosome_values = []
        for gene in adata.var_names:
            if gene in gene_to_seqtype:
                new_chromosome_values.append(gene_to_seqtype[gene])
            else:
                new_chromosome_values.append("unknown")

        adata.var["Chromosome"] = new_chromosome_values
        print("New Chromosome column created successfully!")

        # Show unique values in the new Chromosome column
        unique_chromosomes = adata.var["Chromosome"].unique()
        print(f"Unique values in new Chromosome column: {unique_chromosomes}")

        # Count genes per chromosome
        chromosome_counts = adata.var["Chromosome"].value_counts()
        print("Number of genes per chromosome:")
        for chrom, count in chromosome_counts.items():
            print(f"  {chrom}: {count} genes")
else:
    print("'final_gene_name' column not found in annotation file!")
    print("Available columns in annotation file:")
    pp.pprint(annotation_df.columns.tolist())


# Create new Group variable by processing OD values
print("🔧 Creating new Group variable from OD column...")


def create_group_from_od(od_value):
    """
    Create Group variable from OD value.
    Extracts the medium and OD parts, removing the middle letter (A, B, C).

    Example: "OD M9_C_OD3" -> "M9_OD3"
    Example: "M9_C_OD3" -> "M9_OD3"
    """
    try:
        # Remove "OD " prefix if present
        if od_value.startswith("OD "):
            od_value = od_value[3:]  # Remove "OD "

        # Split by underscore to get parts
        parts = od_value.split("_")

        if len(parts) >= 3:
            # Get the medium part (first part) and OD part (last part)
            medium_part = parts[0]  # e.g., "M9"
            od_part = parts[-1]  # e.g., "OD3"

            # Combine them, removing the middle letter (A, B, C)
            group = f"{medium_part}_{od_part}"
            return group
        else:
            # If format is different, return original value
            return od_value
    except:
        # If any error occurs, return original value
        return od_value


# Apply the function to create Group variable
adata.obs["Group"] = adata.obs["OD"].apply(create_group_from_od)

# Display the mapping and results
print(f"📋 OD to Group mapping:")
od_to_group = {}
for od_val in adata.obs["OD"].unique():
    group_val = create_group_from_od(od_val)
    od_to_group[od_val] = group_val
    print(f"   {od_val} -> {group_val}")

print(f"\n🧬 New Group column created!")
print(f"🔢 Number of unique Group values: {len(adata.obs['Group'].unique())}")
print(f"📋 Unique Group values: {sorted(adata.obs['Group'].unique())}")

# Show distribution of cells per Group
print(f"\n📊 Distribution of cells per Group:")
group_counts = adata.obs["Group"].value_counts()
for group, count in group_counts.items():
    print(f"   {group}: {count:,} cells")


# %%
for col in ["CultureMedium", "OD", "RepBio", "RepTech", "Group"]:
    print(f"🧬 Column: {col}")
    print(f"🔢 Number of unique values: {len(adata.obs[col].unique())}")
    print(f"📋 Unique values: {adata.obs[col].unique()}\n")


# %%
# Calcul de statistiques simples
adata.obs["n_counts"] = adata.X.sum(axis=1).A1  # total counts
adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1  # nombre de gènes détectés


# %%
adata.obs

# %%
adata.obs["n_counts"].describe()
# %%
adata.obs["n_genes"].describe()


# # %% filter cells by n_counts
# sc.pp.filter_cells(adata, min_counts=100)


# sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)


# %% filter 100UMI per cell
# Filter cells with n_counts < 100
print(f"📊 Before filtering: {adata.n_obs} cells")

# Create mask for cells to keep (n_counts >= 100)
cells_to_keep = adata.obs["n_counts"] >= 100

# Apply filtering
adata_filtered = adata[cells_to_keep, :].copy()

print(f"📊 After filtering (n_counts >= 100): {adata_filtered.n_obs} cells")
print(f"🗑️ Removed {adata.n_obs - adata_filtered.n_obs} cells")

# Update adata to the filtered version
adata = adata_filtered

# Show summary statistics after filtering
print(f"\n📈 Summary statistics after filtering:")
print(f"   Mean n_counts: {adata.obs['n_counts'].mean():.2f}")
print(f"   Median n_counts: {adata.obs['n_counts'].median():.2f}")
print(f"   Min n_counts: {adata.obs['n_counts'].min():.2f}")
print(f"   Max n_counts: {adata.obs['n_counts'].max():.2f}")


# # %%
# sc.pl.violin(adata, ["n_counts"], groupby="CultureMedium", multi_panel=True)
# sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)
# # %%
# sc.pl.violin(adata, ["n_counts"], groupby="OD", multi_panel=True)
# # %%
# sc.pl.violin(adata, ["n_counts"], groupby="RepBio", multi_panel=True)
# # %%
# sc.pl.violin(adata, ["n_counts"], groupby="RepTech", multi_panel=True)


# %%
# Get unique culture mediums and create color map
culture_mediums = adata.obs["CultureMedium"].unique()
# Create color map with red for M9 and green for M9F
color_map = {}
for medium in culture_mediums:
    if "M9" in medium and "F" not in medium:  # M9
        color_map[medium] = "#FF6B6B"  # Red
    elif "M9F" in medium:  # M9F
        color_map[medium] = "#2ECC71"  # Green
    else:  # Other mediums
        color_map[medium] = "#95A5A6"  # Gray for other

# Print statistics by culture medium
print(f"\n🍽️ Statistics by Culture Medium:")
for medium in culture_mediums:
    mask = adata.obs["CultureMedium"] == medium
    medium_data = adata.obs.loc[mask]
    print(f"\n   {medium}:")
    print(f"     Cells: {len(medium_data):,}")
    print(f"     Mean n_counts: {medium_data['n_counts'].mean():.0f}")
    print(f"     Mean n_genes: {medium_data['n_genes'].mean():.0f}")

# %%
# Create separate scatter plots for each culture medium with density plots below
# Reorder culture mediums: M9F first (left), M9 second (right)
culture_mediums_ordered = []
for medium in culture_mediums:
    if "M9F" in medium:
        culture_mediums_ordered.append(medium)
for medium in culture_mediums:
    if "M9" in medium and "F" not in medium:
        culture_mediums_ordered.append(medium)
for medium in culture_mediums:
    if "M9" not in medium:
        culture_mediums_ordered.append(medium)

n_mediums = len(culture_mediums_ordered)
n_cols = min(3, n_mediums)  # Maximum 3 columns
n_rows = 2 * n_mediums  # 2 rows per medium (scatter + density)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3 * n_rows))
if n_rows == 1:
    axes = axes.reshape(1, -1)
if n_cols == 1:
    axes = axes.reshape(-1, 1)

# Flatten axes for easier iteration
axes_flat = axes.flatten()

# Calculate global min/max for consistent scaling
global_x_min = adata.obs["n_counts"].min()
global_x_max = adata.obs["n_counts"].max()
global_y_min = adata.obs["n_genes"].min()
global_y_max = adata.obs["n_genes"].max()

# Add 5% padding to global limits
x_padding = (global_x_max - global_x_min) * 0.05
y_padding = (global_y_max - global_y_min) * 0.05

global_x_lim = (global_x_min - x_padding, global_x_max + x_padding)
global_y_lim = (global_y_min - y_padding, global_y_max + y_padding)

for i, medium in enumerate(culture_mediums_ordered):
    # Calculate correct indices for 2-row layout
    col_idx = i % n_cols
    scatter_row = (i // n_cols) * 2
    density_row = (i // n_cols) * 2 + 1

    # Get data for this medium
    mask = adata.obs["CultureMedium"] == medium
    medium_data = adata.obs.loc[mask]

    # Create scatter plot (top row)
    scatter_idx = scatter_row * n_cols + col_idx
    ax_scatter = axes_flat[scatter_idx]
    ax_scatter.scatter(
        medium_data["n_counts"],
        medium_data["n_genes"],
        alpha=0.7,
        s=8,  # Reduced point size from 15 to 8
        c=color_map[medium],
        edgecolors="none",
    )

    # Calculate correlation for this medium
    if len(medium_data) > 1:
        corr = np.corrcoef(medium_data["n_counts"], medium_data["n_genes"])[0, 1]
        # Position correlation text in upper right corner
        ax_scatter.text(
            0.98,
            0.98,
            f"Correlation: {corr:.3f}",
            transform=ax_scatter.transAxes,
            fontsize=9,
            ha="right",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
            ),
        )

    # Add barcode count in upper left corner
    barcode_text = f"Barcodes: {len(medium_data):,}"
    ax_scatter.text(
        0.02,
        0.98,
        barcode_text,
        transform=ax_scatter.transAxes,
        fontsize=8,
        ha="left",
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
        ),
    )

    # Customize scatter subplot
    ax_scatter.set_xlabel("Total Counts", fontsize=10, fontweight="bold")
    ax_scatter.set_ylabel("Number of Genes", fontsize=10, fontweight="bold")
    ax_scatter.set_title(
        f"{medium} - Scatter Plot", fontsize=12, fontweight="bold", pad=10
    )
    ax_scatter.grid(True, alpha=0.3, linestyle="--")

    # Set consistent axis limits for all scatter subplots
    ax_scatter.set_xlim(global_x_lim)
    ax_scatter.set_ylim(global_y_lim)

    # Create density plot (bottom row)
    density_idx = density_row * n_cols + col_idx
    ax_density = axes_flat[density_idx]

    # Create density plot for n_counts with global x limits
    ax_density.hist(
        medium_data["n_counts"],
        bins=30,
        alpha=0.7,
        color=color_map[medium],
        edgecolor="black",
        linewidth=0.5,
        density=True,
        range=global_x_lim,  # Use global x limits for consistent scaling
    )

    # Add KDE curve with global x range
    from scipy.stats import gaussian_kde

    if len(medium_data) > 1:
        kde = gaussian_kde(medium_data["n_counts"])
        x_range = np.linspace(
            global_x_lim[0], global_x_lim[1], 100
        )  # Use global x range
        ax_density.plot(x_range, kde(x_range), color="black", linewidth=2, label="KDE")

    # Customize density subplot
    ax_density.set_xlabel("Total Counts", fontsize=10, fontweight="bold")
    ax_density.set_ylabel("Density", fontsize=10, fontweight="bold")
    ax_density.set_title(
        f"{medium} - Distribution", fontsize=12, fontweight="bold", pad=10
    )
    ax_density.grid(True, alpha=0.3, linestyle="--")
    ax_density.legend(fontsize=8)

    # Set consistent x-axis limits for all density plots
    ax_density.set_xlim(global_x_lim)

# Hide empty subplots
for i in range(n_mediums * 2, len(axes_flat)):
    axes_flat[i].set_visible(False)

plt.suptitle(
    "Cell Quality by Culture Medium with Density Distributions (Filtered: n_counts ≥ 100)",
    fontsize=16,
    fontweight="bold",
    y=0.98,
)
plt.tight_layout()
plt.show()

# Print detailed statistics for each medium
print(f"\n📊 Detailed Statistics by Culture Medium (After filtering: n_counts ≥ 100):")
for medium in culture_mediums_ordered:
    mask = adata.obs["CultureMedium"] == medium
    medium_data = adata.obs.loc[mask]

    print(f"\n🍽️ {medium}:")
    print(f"   📊 Total cells: {len(medium_data):,}")
    print(
        f"   📈 n_counts - Mean: {medium_data['n_counts'].mean():.0f}, Median: {medium_data['n_counts'].median():.0f}"
    )
    print(
        f"   📈 n_counts - Min: {medium_data['n_counts'].min():.0f}, Max: {medium_data['n_counts'].max():.0f}"
    )
    print(
        f"   🧬 n_genes - Mean: {medium_data['n_genes'].mean():.0f}, Median: {medium_data['n_genes'].median():.0f}"
    )
    print(
        f"   🧬 n_genes - Min: {medium_data['n_genes'].min():.0f}, Max: {medium_data['n_genes'].max():.0f}"
    )

    if len(medium_data) > 1:
        corr = np.corrcoef(medium_data["n_counts"], medium_data["n_genes"])[0, 1]
        print(f"   🔗 Correlation n_counts vs n_genes: {corr:.3f}")

# %%
# %%
# Create violin plots comparing M9 and M9F
# Filter data for M9 and M9F only
m9_m9f_mask = adata.obs["CultureMedium"].isin(["M9", "M9F"])
adata_m9_m9f = adata[m9_m9f_mask].copy()

if len(adata_m9_m9f) > 0:
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Violin plot for n_counts (top)
    sns.violinplot(
        data=adata_m9_m9f.obs,
        x="CultureMedium",
        y="n_counts",
        ax=ax1,
        palette=["#FF6B6B", "#2ECC71"],
    )
    ax1.set_title(
        "Distribution of Total Counts per Barcode (cell) by Culture Medium",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    ax1.set_xlabel("Culture Medium", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Total Counts per Barcode", fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.3, linestyle="--")

    # Add statistics on the plot
    for i, medium in enumerate(["M9", "M9F"]):
        if medium in adata_m9_m9f.obs["CultureMedium"].values:
            medium_data = adata_m9_m9f.obs[adata_m9_m9f.obs["CultureMedium"] == medium]
            median_val = medium_data["n_counts"].median()
            q1_val = medium_data["n_counts"].quantile(0.25)
            q3_val = medium_data["n_counts"].quantile(0.75)

            # Position M9 on the left, M9F on the right
            if medium == "M9":
                ax1.text(
                    0.02,
                    0.98,
                    f"{medium}:\nMedian: {median_val:.0f}\nQ1: {q1_val:.0f}\nQ3: {q3_val:.0f}",
                    transform=ax1.transAxes,
                    fontsize=10,
                    ha="left",
                    va="top",
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        facecolor="white",
                        alpha=0.9,
                        edgecolor="gray",
                    ),
                )
            else:  # M9F
                ax1.text(
                    0.98,
                    0.98,
                    f"{medium}:\nMedian: {median_val:.0f}\nQ1: {q1_val:.0f}\nQ3: {q3_val:.0f}",
                    transform=ax1.transAxes,
                    fontsize=10,
                    ha="right",
                    va="top",
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        facecolor="white",
                        alpha=0.9,
                        edgecolor="gray",
                    ),
                )

    # Violin plot for n_genes (bottom)
    sns.violinplot(
        data=adata_m9_m9f.obs,
        x="CultureMedium",
        y="n_genes",
        ax=ax2,
        palette=["#FF6B6B", "#2ECC71"],
    )
    ax2.set_title(
        "Distribution of Genes Detected per Barcode (cell) by Culture Medium",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    ax2.set_xlabel("Culture Medium", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Genes per Barcode", fontsize=12, fontweight="bold")
    ax2.grid(True, alpha=0.3, linestyle="--")

    # Add statistics on the plot
    for i, medium in enumerate(["M9", "M9F"]):
        if medium in adata_m9_m9f.obs["CultureMedium"].values:
            medium_data = adata_m9_m9f.obs[adata_m9_m9f.obs["CultureMedium"] == medium]
            median_val = medium_data["n_genes"].median()
            q1_val = medium_data["n_genes"].quantile(0.25)
            q3_val = medium_data["n_genes"].quantile(0.75)

            # Position M9 on the left, M9F on the right
            if medium == "M9":
                ax2.text(
                    0.02,
                    0.98,
                    f"{medium}:\nMedian: {median_val:.0f}\nQ1: {q1_val:.0f}\nQ3: {q3_val:.0f}",
                    transform=ax2.transAxes,
                    fontsize=10,
                    ha="left",
                    va="top",
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        facecolor="white",
                        alpha=0.9,
                        edgecolor="gray",
                    ),
                )
            else:  # M9F
                ax2.text(
                    0.98,
                    0.98,
                    f"{medium}:\nMedian: {median_val:.0f}\nQ1: {q1_val:.0f}\nQ3: {q3_val:.0f}",
                    transform=ax2.transAxes,
                    fontsize=10,
                    ha="right",
                    va="top",
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        facecolor="white",
                        alpha=0.9,
                        edgecolor="gray",
                    ),
                )

    plt.suptitle(
        "Comparison of Cell Quality Metrics: M9 vs M9F (Filtered: n_counts ≥ 100)",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout()
    plt.show()

    # Print comparison statistics
    print(f"\n🎻 Violin Plot Comparison: M9 vs M9F")
    print(f"📊 Total cells in comparison: {len(adata_m9_m9f):,}")

    for medium in ["M9", "M9F"]:
        if medium in adata_m9_m9f.obs["CultureMedium"].values:
            medium_data = adata_m9_m9f.obs[adata_m9_m9f.obs["CultureMedium"] == medium]
            print(f"\n🍽️ {medium}:")
            print(f"   📊 Cells: {len(medium_data):,}")
            print(
                f"   📈 n_counts - Mean: {medium_data['n_counts'].mean():.0f}, Median: {medium_data['n_counts'].median():.0f}"
            )
            print(f"   📈 n_counts - Std: {medium_data['n_counts'].std():.0f}")
            print(
                f"   🧬 n_genes - Mean: {medium_data['n_genes'].mean():.0f}, Median: {medium_data['n_genes'].median():.0f}"
            )
            print(f"   🧬 n_genes - Std: {medium_data['n_genes'].std():.0f}")
else:
    print("⚠️ No M9 or M9F data found for violin plot comparison")


# %%
# Create beautiful plot of number of cells per culture medium
plt.figure(figsize=(14, 8))

# Count cells per culture medium BEFORE filtering
cell_counts_before = adata_original.obs["CultureMedium"].value_counts()

# Count cells per culture medium AFTER filtering
cell_counts_after = adata.obs["CultureMedium"].value_counts()

# Reorder culture mediums: M9F first (left), M9 second (right)
culture_mediums_ordered = []
for medium in cell_counts_before.index:
    if "M9F" in medium:
        culture_mediums_ordered.append(medium)
for medium in cell_counts_before.index:
    if "M9" in medium and "F" not in medium:
        culture_mediums_ordered.append(medium)
for medium in cell_counts_before.index:
    if "M9" not in medium:
        culture_mediums_ordered.append(medium)

# Reorder the data according to the new order
cell_counts_before_ordered = cell_counts_before.reindex(culture_mediums_ordered)
cell_counts_after_ordered = cell_counts_after.reindex(culture_mediums_ordered)

# Create color palette
colors = []
for medium in culture_mediums_ordered:
    if "M9" in medium and "F" not in medium:  # M9
        colors.append("#FF6B6B")  # Red
    elif "M9F" in medium:  # M9F
        colors.append("#2ECC71")  # Green
    else:  # Other mediums
        colors.append("#95A5A6")  # Gray for other

# Set up the bar positions
x = np.arange(len(culture_mediums_ordered))
width = 0.35

# Create bar plot
bars_before = plt.bar(
    x - width / 2,
    cell_counts_before_ordered.values,
    width,
    color=colors,
    alpha=0.7,
    edgecolor="black",
    linewidth=1,
)
bars_after = plt.bar(
    x + width / 2,
    cell_counts_after_ordered.values,
    width,
    color=colors,
    alpha=0.9,
    edgecolor="black",
    linewidth=1,
)

# Customize the plot
plt.title(
    "Number of Barcodes per Culture Medium: Not Filtered vs Filtered (n_counts > 100)",
    fontsize=16,
    fontweight="bold",
    pad=20,
)
plt.xlabel("Culture Medium", fontsize=14, fontweight="bold")
plt.ylabel("Number of Barcodes", fontsize=14, fontweight="bold")

# Set x-axis labels
plt.xticks(
    x, culture_mediums_ordered, rotation=0, ha="center", fontsize=12, fontweight="bold"
)

# Add "Before" and "After" labels under each bar group
for i, medium in enumerate(culture_mediums_ordered):
    plt.text(
        x[i] - width / 2,
        -max(cell_counts_before_ordered.values) * 0.05,
        "Not filtered",
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )
    plt.text(
        x[i] + width / 2,
        -max(cell_counts_before_ordered.values) * 0.05,
        "Filtered\n(n_counts > 100)",
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )

# Add value labels on top of each bar
for bar, count in zip(bars_before, cell_counts_before_ordered.values):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max(cell_counts_before_ordered.values) * 0.01,
        f"{count:,}",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

for bar, count in zip(bars_after, cell_counts_after_ordered.values):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max(cell_counts_before_ordered.values) * 0.01,
        f"{count:,}",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

# Add percentage labels inside bars
total_before = cell_counts_before_ordered.sum()
total_after = cell_counts_after_ordered.sum()

for bar, count in zip(bars_before, cell_counts_before_ordered.values):
    percentage = (count / total_before) * 100
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() / 2,
        f"{percentage:.1f}%",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        color="black",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
        ),
    )

for bar, count in zip(bars_after, cell_counts_after_ordered.values):
    percentage = (count / total_after) * 100
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() / 2,
        f"{percentage:.1f}%",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        color="black",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
        ),
    )

# Add grid for better readability
plt.grid(True, alpha=0.3, linestyle="--", axis="y")

# Add total counts in the upper right corner
total_text = f"Total Not filtered: {total_before:,}\nTotal Filtered: {total_after:,}\nRemoved: {total_before - total_after:,}"
plt.text(
    0.98,
    0.98,
    total_text,
    transform=plt.gca().transAxes,
    fontsize=11,
    fontweight="bold",
    ha="right",
    va="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"),
)

# Show the plot
plt.show()

# Print detailed statistics
print(f"\n📊 Cell Count Summary: Before vs After Filtering (n_counts ≥ 100):")
print(f"📈 Total barcodes before filtering: {total_before:,}")
print(f"📈 Total barcodes after filtering: {total_after:,}")
print(
    f"🗑️ Barcodes removed: {total_before - total_after:,} ({(total_before - total_after)/total_before*100:.1f}%)"
)

print(f"\n🍽️ Breakdown by Culture Medium:")
for medium in culture_mediums_ordered:
    before_count = cell_counts_before_ordered[medium]
    after_count = cell_counts_after_ordered.get(medium, 0)
    removed_count = before_count - after_count
    removed_percentage = (removed_count / before_count) * 100 if before_count > 0 else 0

    print(f"\n   {medium}:")
    print(f"     Before: {before_count:,} barcodes")
    print(f"     After: {after_count:,} barcodes")
    print(f"     Removed: {removed_count:,} barcodes ({removed_percentage:.1f}%)")

# %%
adata

adata.var["Gene Type"].value_counts()


# %% filter


def select_cells_by_group(
    adata: AnnData,
    target_total_cells: int,
    group_by: str = "Condition",
    sort_by: str = "n_counts",  # souvent stocké dans adata.obs
    random: bool = False,
    seed: int = None,
    descending: bool = True,
):
    meta = adata.obs.copy()
    meta["cell_barcode"] = meta.index

    if group_by not in meta.columns:
        raise ValueError(f"Group column '{group_by}' not found in adata.obs.")
    if sort_by not in meta.columns:
        raise ValueError(f"Sort column '{sort_by}' not found in adata.obs.")

    n_groups = meta[group_by].nunique()
    target_per_group = int(target_total_cells // n_groups)

    print(f"🎯 Target total cells: {target_total_cells}")
    print(f"📁 Group column: '{group_by}' — {n_groups} groups")
    print(f"🧮 Sort by: '{sort_by}' | Descending: {descending}")
    print(f"🔢 Cells per group: {target_per_group}")
    if random and seed is not None:
        print(f"🎲 Random sampling with seed {seed}")

    selected_cells = []

    for group, df_group in meta.groupby(group_by):
        n_available = df_group.shape[0]
        print(f"🔍 Group '{group}': {n_available} cells available")

        if n_available <= target_per_group:
            print(f"✅ Keeping all {n_available} cells")
            selected = df_group["cell_barcode"].tolist()
        else:
            if random:
                selected = df_group.sample(n=target_per_group, random_state=seed)[
                    "cell_barcode"
                ].tolist()
            else:
                df_sorted = df_group.sort_values(by=sort_by, ascending=not descending)
                selected = df_sorted.head(target_per_group)["cell_barcode"].tolist()
        selected_cells.extend(selected)

    adata_filtered = adata[selected_cells].copy()

    # Résumé par groupe
    summary = adata_filtered.obs.groupby(group_by).size().reset_index(name="n_cells")

    return {"adata": adata_filtered, "summary": summary}


# %% remove the 6 top ncounts cells
print("🗑️ Removing the 6 cells with highest n_counts...")

# Get the 3 highest n_counts values
top_6_ncounts = adata.obs["n_counts"].nlargest(6)
print(f"📊 Top 6 n_counts values to remove:")
for i, (idx, value) in enumerate(top_6_ncounts.items(), 1):
    print(f"   {i}. Cell {idx}: {value:,} counts")

# Create mask to exclude these 3 cells
cells_to_remove = top_6_ncounts.index
cells_to_keep = ~adata.obs.index.isin(cells_to_remove)

# Apply filtering
adata = adata[cells_to_keep].copy()

print(f"📊 After removing top 6 cells: {adata.n_obs:,} cells")
print(f"🗑️ Removed {len(cells_to_remove)} cells with highest n_counts")

# %% remove the highest M9F_OD2 cell
print("🗑️ Removing the highest M9F_OD2 cell...")
# Get the highest n_counts value for M9F_OD2
highest_m9f_od2 = adata.obs[adata.obs["Group"] == "M9F_OD2"]["n_counts"].nlargest(1)
print(f"📊 Highest n_counts for M9F_OD2: {highest_m9f_od2.iloc[0]:,} counts")
# Create mask to exclude this cell
cells_to_remove = highest_m9f_od2.index

adata = adata[~adata.obs.index.isin(cells_to_remove)].copy()


# %%
sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)


# %% remove the highest M9_OD2 cell
print("🗑️ Removing the highest M9_OD2 cell...")
# Get the highest n_counts value for M9_OD2
highest_m9_od2 = adata.obs[adata.obs["Group"] == "M9_OD2"]["n_counts"].nlargest(1)
print(f"📊 Highest n_counts for M9_OD2: {highest_m9_od2.iloc[0]:,} counts")
# Create mask to exclude this cell
cells_to_remove = highest_m9_od2.index

adata = adata[~adata.obs.index.isin(cells_to_remove)].copy()

# %%
sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)

# %% remove the highest M9_OD1 cell
print("🗑️ Removing the highest M9_OD1 cell...")
# Get the highest n_counts value for M9_OD1
highest_m9_od1 = adata.obs[adata.obs["Group"] == "M9_OD1"]["n_counts"].nlargest(1)
print(f"📊 Highest n_counts for M9_OD1: {highest_m9_od1.iloc[0]:,} counts")
# Create mask to exclude this cell
cells_to_remove = highest_m9_od1.index

adata = adata[~adata.obs.index.isin(cells_to_remove)].copy()

# %%
sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)

# %%
adata.obs["n_counts"].describe()
# %% n cell by group
adata.obs["Group"].value_counts()
# %%


# %%
# Ajoute les totaux d'ARN si nécessaire
adata.obs["n_counts"] = (
    adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
)

# Sous-échantillonne 1000 cellules par condition
result = select_cells_by_group(
    adata=adata,
    target_total_cells=3000,
    group_by="RepTech",
    sort_by="n_counts",
    random=False,
    seed=42,
    descending=True,
)

adata_filtered = result["adata"]
print(result["summary"])

adata = adata_filtered.copy()

# %%
sc.pl.violin(adata, ["n_counts"], groupby="Group", multi_panel=True)

# %%
sc.pl.violin(adata, ["n_counts"], groupby="CultureMedium", multi_panel=True)

# %%
sc.pl.violin(adata, ["n_counts"], groupby="OD", multi_panel=True, rotation=45)

# %%
sc.pl.violin(
    adata,
    keys=["n_counts"],
    groupby="OD",
    rotation=45,
    multi_panel=True,
    stripplot=True,  # optionnel
    jitter=True,  # optionnel
)

# %%
import pandas as pd

# Extraire les colonnes nécessaires
df = adata.obs[["n_counts", "RepBio", "ODt"]].copy()

# S'assurer que les colonnes sont au bon type
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)
# %%
import seaborn as sns
import matplotlib.pyplot as plt

# Créer le catplot
g = sns.catplot(
    data=df,
    x="RepBio",
    y="n_counts",
    col="ODt",  # Facet par Odt
    kind="violin",
    sharey=False,  # Facultatif : échelles indépendantes pour chaque panel
    col_wrap=3,  # Si beaucoup de niveaux de Odt
    height=4,
    aspect=1,
)

# Rotation des labels x
for ax in g.axes.flatten():
    ax.tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.show()

# %% # Palette dynamique basée sur le préfixe (M9 vs M9F)
unique_groups = df["RepBio"].unique()
palette = {grp: ("red" if grp.startswith("M9_") else "green") for grp in unique_groups}

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Extraire les colonnes nécessaires
df = adata.obs[["n_counts", "RepBio", "ODt"]].copy()
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)

# Compter le nombre de cellules par groupe
cell_counts = df.groupby(["ODt", "RepBio"]).size().reset_index(name="n_cells")

# Palette dynamique
unique_groups = df["RepBio"].unique()
palette = {grp: ("red" if grp.startswith("M9_") else "green") for grp in unique_groups}

# Créer le catplot
g = sns.catplot(
    data=df,
    x="RepBio",
    y="n_counts",
    col="ODt",
    hue="RepBio",
    kind="violin",
    palette=palette,
    sharey=False,
    col_wrap=3,
    height=4,
    aspect=1,
)

# Annoter les plots avec le nombre de cellules
for ax in g.axes.flatten():
    odt = ax.get_title().split(" = ")[-1]
    for tick, label in enumerate(ax.get_xticklabels()):
        repbio = label.get_text()
        count_row = cell_counts[
            (cell_counts["ODt"] == odt) & (cell_counts["RepBio"] == repbio)
        ]
        if not count_row.empty:
            n = int(count_row["n_cells"].values[0])
            ax.text(
                tick,
                ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                f"n={n}",
                ha="center",
                va="top",
                fontsize=9,
                color=palette[repbio],
            )

# Rotation des étiquettes x
for ax in g.axes.flatten():
    ax.tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.show()


# %% see min max n counts by group
adata.obs["n_counts"].describe()
# %%  see min max n counts by group OD
adata.obs.groupby("OD")["n_counts"].describe()

# %% create a new variable Reptechg qui est le suffixe de RepTech ..._1
adata.obs["RepTech_suffix"] = adata.obs["RepTech"].str.split("_").str[-1]
# %%
sc.pl.violin(adata, ["n_counts"], groupby="RepTech_suffix", multi_panel=True)


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Créer la colonne suffixe
adata.obs["RepTech_suffix"] = adata.obs["RepTech"].astype(str).str.split("_").str[-1]

# Extraire les colonnes nécessaires
df = adata.obs[["n_counts", "RepBio", "ODt", "RepTech_suffix"]].copy()
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)
df["RepTech_suffix"] = df["RepTech_suffix"].astype(str)

# Compter les cellules par groupe
cell_counts = df.groupby(["ODt", "RepBio"]).size().reset_index(name="n_cells")

# Palette dynamique pour RepBio
unique_groups = df["RepBio"].unique()
palette = {grp: ("red" if grp.startswith("M9_") else "green") for grp in unique_groups}

# Marqueurs différents pour chaque RepTech_suffix
unique_suffixes = sorted(df["RepTech_suffix"].unique())
reptech_markers = {
    val: marker
    for val, marker in zip(
        unique_suffixes,
        ["o", "s", "D", "^", "v", "P", "X", "*"],  # jusqu'à 8 formes distinctes
    )
}

# Base plot : violin plot
g = sns.catplot(
    data=df,
    x="RepBio",
    y="n_counts",
    col="ODt",
    hue="RepBio",
    kind="violin",
    palette=palette,
    sharey=True,
    col_wrap=3,
    height=4,
    aspect=1,
    inner=None,
)

# Ajouter les points avec jitter et distinction par RepTech_suffix, avec offset manuel
for ax in g.axes.flatten():
    odt = ax.get_title().split(" = ")[-1]
    data_subset = df[df["ODt"] == odt]

    # Obtenir positions x des catégories RepBio
    x_ticks = ax.get_xticks()
    x_labels = [tick.get_text() for tick in ax.get_xticklabels()]
    repbio_to_x = dict(zip(x_labels, x_ticks))

    for j, suffix in enumerate(unique_suffixes):
        offset = (j - len(unique_suffixes) / 2) * 0.1
        subdata = data_subset[data_subset["RepTech_suffix"] == suffix]

        for repbio in subdata["RepBio"].unique():
            sub_subdata = subdata[subdata["RepBio"] == repbio]
            x = repbio_to_x.get(repbio)
            if x is None:
                continue
            x_offset = x + offset

            ax.scatter(
                np.full(len(sub_subdata), x_offset),
                sub_subdata["n_counts"],
                marker=reptech_markers[suffix],
                color="black",
                edgecolor="white",
                linewidth=0.4,
                s=50,
                alpha=0.8,
            )

# Supprimer légendes automatiques multiples
for ax in g.axes.flatten():
    if ax.legend_:
        ax.legend_.remove()

# Légende manuelle pour RepTech_suffix
legend_elements = [
    Line2D(
        [0],
        [0],
        marker=reptech_markers[suffix],
        color="black",
        linestyle="",
        label=f"RepTech_{suffix}",
        markersize=7,
        markerfacecolor="black",
    )
    for suffix in unique_suffixes
]
g.axes[0].legend(handles=legend_elements, title="RepTech_suffix", loc="upper right")

# Ajouter n=... en haut (y = 5000), en noir
for ax in g.axes.flatten():
    odt = ax.get_title().split(" = ")[-1]
    for tick, label in enumerate(ax.get_xticklabels()):
        repbio = label.get_text()
        count_row = cell_counts[
            (cell_counts["ODt"] == odt) & (cell_counts["RepBio"] == repbio)
        ]
        if not count_row.empty:
            n = int(count_row["n_cells"].values[0])
            ax.text(
                tick,
                5000,
                f"n={n}",
                ha="center",
                va="bottom",
                fontsize=9,
                color="black",
                weight="bold",
            )

# Rotation des étiquettes
for ax in g.axes.flatten():
    ax.tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.show()

# %% create a new variable RepBio_suffix which is the suffix of RepBio ..._1
adata.obs["RepBio_suffix"] = adata.obs["RepBio"].str.split("_").str[-1]
# %%

# Create dataframe from AnnData object
df = pd.DataFrame(adata.obs)
df["n_counts"] = adata.obs["n_counts"]

# Plot using seaborn catplot
sns.catplot(
    data=df, x="RepBio", y="n_counts", hue="RepTech_suffix", col="ODt", aspect=0.5
)


# %% correct


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Preprocessing
df = adata.obs.copy()
df["n_counts"] = adata.obs["n_counts"]
df["RepTech_suffix"] = df["RepTech"].astype(str).str.split("_").str[-1]
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)
df["RepTech_suffix"] = df["RepTech_suffix"].astype(str)

# Violin plot sans légende
g = sns.catplot(
    data=df,
    x="RepBio",
    y="n_counts",
    hue=None,  # Ne pas générer de légende ici
    col="ODt",
    kind="violin",
    dodge=True,
    col_wrap=3,
    height=4,
    aspect=1,
    sharey=True,
    inner=None,
    linewidth=1,
)

# Stripplot avec légende récupérée depuis un seul axe
handles, labels = None, None
for ax in g.axes.flatten():
    odt = ax.get_title().split(" = ")[-1]
    data_subset = df[df["ODt"] == odt]

    strip = sns.stripplot(
        data=data_subset,
        x="RepBio",
        y="n_counts",
        hue="RepTech_suffix",
        dodge=True,
        jitter=0.25,
        size=4,
        ax=ax,
        alpha=0.8,
        edgecolor="white",
        linewidth=0.4,
        palette="dark",
    )

    # Récupère les handles une seule fois
    if handles is None and ax.get_legend() is not None:
        handles, labels = ax.get_legend_handles_labels()

    ax.legend_.remove()

# Affichage de la légende des points seulement
if handles is not None:
    g.fig.legend(
        handles=handles,
        labels=labels,
        title="RepTech_suffix",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )

# Rotation des étiquettes
for ax in g.axes.flatten():
    ax.tick_params(axis="x", rotation=45)

plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.show()


# %%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

# Preprocessing
df = adata.obs.copy()
df["n_counts"] = adata.obs["n_counts"]
df["RepTech_suffix"] = df["RepTech"].astype(str).str.split("_").str[-1]
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)
df["RepTech_suffix"] = df["RepTech_suffix"].astype(str)

# Paramètres
col_wrap = 3
palette = sns.color_palette("dark")
repbios = sorted(df["RepBio"].unique())
odts = sorted(df["ODt"].unique())
repbio_to_pos = {repbio: i for i, repbio in enumerate(repbios)}

# Setup grid
n_cols = col_wrap
n_rows = int(np.ceil(len(odts) / col_wrap))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows), sharey=True)
axes = axes.flatten()

handles, labels = None, None

# Boucle par facette (ODt)
for i, odt in enumerate(odts):
    ax = axes[i]
    data_odt = df[df["ODt"] == odt]

    for j, repbio in enumerate(repbios):
        data_subset = data_odt[data_odt["RepBio"] == repbio]
        values = data_subset["n_counts"].dropna()

        if len(values) < 2:
            continue  # KDE échoue avec <2 points

        # Densité limitée au min/max
        kde = gaussian_kde(values, bw_method=0.3)
        x_vals = np.linspace(values.min(), values.max(), 100)
        density = kde(x_vals)

        # Normaliser pour que tous les violons aient la même largeur max
        width = 0.4
        density = density / density.max() * width

        # Tracer le "violon"
        ax.fill_betweenx(
            x_vals, j - density, j + density, color="lightblue", alpha=0.7, linewidth=1
        )

    # Stripplot pour les points
    sns.stripplot(
        data=data_odt,
        x="RepBio",
        y="n_counts",
        hue="RepTech_suffix",
        dodge=True,
        jitter=0.25,
        size=4,
        ax=ax,
        alpha=0.9,
        edgecolor="white",
        linewidth=0.4,
        palette=palette,
        order=repbios,
    )

    # Récupère la légende une fois
    if handles is None and ax.get_legend() is not None:
        handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend():
        ax.legend_.remove()

    ax.set_title(f"ODt = {odt}")
    ax.set_xticks(range(len(repbios)))
    ax.set_xticklabels(repbios, rotation=45)

# Enlève les axes inutilisés si nb ODt < n_rows × n_cols
for j in range(len(odts), len(axes)):
    fig.delaxes(axes[j])

# Légende commune
if handles is not None:
    fig.legend(
        handles,
        labels,
        title="RepTech_suffix",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )

plt.tight_layout()
plt.show()


# %%
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

# Preprocessing
df = adata.obs.copy()
df["n_counts"] = adata.obs["n_counts"]
df["RepTech_suffix"] = df["RepTech"].astype(str).str.split("_").str[-1]
df["ODt"] = df["ODt"].astype(str)
df["RepBio"] = df["RepBio"].astype(str)
df["RepTech_suffix"] = df["RepTech_suffix"].astype(str)

# Paramètres
col_wrap = 3
palette = sns.color_palette("dark")
repbios = sorted(df["RepBio"].unique())
odts = sorted(df["ODt"].unique())
repbio_to_pos = {repbio: i for i, repbio in enumerate(repbios)}

# Setup grid
n_cols = col_wrap
n_rows = int(np.ceil(len(odts) / col_wrap))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows), sharey=True)
axes = axes.flatten()

handles, labels = None, None

for i, odt in enumerate(odts):
    ax = axes[i]
    data_odt = df[df["ODt"] == odt]

    for j, repbio in enumerate(repbios):
        data_subset = data_odt[data_odt["RepBio"] == repbio]
        values = data_subset["n_counts"].dropna()

        if len(values) < 2:
            continue

        # KDE clippée
        kde = gaussian_kde(values, bw_method=0.3)
        x_vals = np.linspace(values.min(), values.max(), 100)
        density = kde(x_vals)
        width = 0.4
        density = density / density.max() * width

        # Tracé du violon
        ax.fill_betweenx(
            x_vals, j - density, j + density, color="lightblue", alpha=0.7, linewidth=1
        )

        # ➕ Affiche le nombre de points
        ax.text(
            j,
            -100,
            f"(n = {len(values)})",
            ha="center",
            va="bottom",
            fontsize=6,
            fontweight="bold",
            color="black",
        )

    # Stripplot
    sns.stripplot(
        data=data_odt,
        x="RepBio",
        y="n_counts",
        hue="RepTech_suffix",
        dodge=True,
        jitter=0.25,
        size=4,
        ax=ax,
        alpha=0.9,
        edgecolor="white",
        linewidth=0.4,
        palette=palette,
        order=repbios,
    )

    if handles is None and ax.get_legend() is not None:
        handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend():
        ax.legend_.remove()

    ax.set_title(f"ODt = {odt}")
    ax.set_xticks(range(len(repbios)))
    ax.set_xticklabels(repbios, rotation=0)  # Changed rotation from 45 to 0

# Supprimer axes vides
for j in range(len(odts), len(axes)):
    fig.delaxes(axes[j])

# Légende commune
if handles is not None:
    fig.legend(
        handles,
        labels,
        title="RepTech_suffix",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )

plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.show()


# %%
# Create separate scatter plots for each culture medium with density plots below
# Reorder culture mediums: M9F first (left), M9 second (right)
culture_mediums_ordered = []
for medium in culture_mediums:
    if "M9F" in medium:
        culture_mediums_ordered.append(medium)
for medium in culture_mediums:
    if "M9" in medium and "F" not in medium:
        culture_mediums_ordered.append(medium)
for medium in culture_mediums:
    if "M9" not in medium:
        culture_mediums_ordered.append(medium)

n_mediums = len(culture_mediums_ordered)
n_cols = min(3, n_mediums)  # Maximum 3 columns
n_rows = 2 * n_mediums  # 2 rows per medium (scatter + density)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3 * n_rows))
if n_rows == 1:
    axes = axes.reshape(1, -1)
if n_cols == 1:
    axes = axes.reshape(-1, 1)

# Flatten axes for easier iteration
axes_flat = axes.flatten()

# Calculate global min/max for consistent scaling
global_x_min = adata.obs["n_counts"].min()
global_x_max = adata.obs["n_counts"].max()
global_y_min = adata.obs["n_genes"].min()
global_y_max = adata.obs["n_genes"].max()

# Add 5% padding to global limits
x_padding = (global_x_max - global_x_min) * 0.05
y_padding = (global_y_max - global_y_min) * 0.05

global_x_lim = (global_x_min - x_padding, global_x_max + x_padding)
global_y_lim = (global_y_min - y_padding, global_y_max + y_padding)

for i, medium in enumerate(culture_mediums_ordered):
    # Calculate correct indices for 2-row layout
    col_idx = i % n_cols
    scatter_row = (i // n_cols) * 2
    density_row = (i // n_cols) * 2 + 1

    # Get data for this medium
    mask = adata.obs["CultureMedium"] == medium
    medium_data = adata.obs.loc[mask]

    # Create scatter plot (top row)
    scatter_idx = scatter_row * n_cols + col_idx
    ax_scatter = axes_flat[scatter_idx]
    ax_scatter.scatter(
        medium_data["n_counts"],
        medium_data["n_genes"],
        alpha=0.7,
        s=8,  # Reduced point size from 15 to 8
        c=color_map[medium],
        edgecolors="none",
    )

    # Calculate correlation for this medium
    if len(medium_data) > 1:
        corr = np.corrcoef(medium_data["n_counts"], medium_data["n_genes"])[0, 1]
        # Position correlation text in upper right corner
        ax_scatter.text(
            0.98,
            0.98,
            f"Correlation: {corr:.3f}",
            transform=ax_scatter.transAxes,
            fontsize=9,
            ha="right",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
            ),
        )

    # Add barcode count in upper left corner
    barcode_text = f"Barcodes: {len(medium_data):,}"
    ax_scatter.text(
        0.02,
        0.98,
        barcode_text,
        transform=ax_scatter.transAxes,
        fontsize=8,
        ha="left",
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray"
        ),
    )

    # Customize scatter subplot
    ax_scatter.set_xlabel("Total Counts", fontsize=10, fontweight="bold")
    ax_scatter.set_ylabel("Number of Genes", fontsize=10, fontweight="bold")
    ax_scatter.set_title(
        f"{medium} - Scatter Plot", fontsize=12, fontweight="bold", pad=10
    )
    ax_scatter.grid(True, alpha=0.3, linestyle="--")

    # Set consistent axis limits for all scatter subplots
    ax_scatter.set_xlim(global_x_lim)
    ax_scatter.set_ylim(global_y_lim)

    # Create density plot (bottom row)
    density_idx = density_row * n_cols + col_idx
    ax_density = axes_flat[density_idx]

    # Create density plot for n_counts with global x limits
    ax_density.hist(
        medium_data["n_counts"],
        bins=30,
        alpha=0.7,
        color=color_map[medium],
        edgecolor="black",
        linewidth=0.5,
        density=True,
        range=global_x_lim,  # Use global x limits for consistent scaling
    )

    # Add KDE curve with global x range
    from scipy.stats import gaussian_kde

    if len(medium_data) > 1:
        kde = gaussian_kde(medium_data["n_counts"])
        x_range = np.linspace(
            global_x_lim[0], global_x_lim[1], 100
        )  # Use global x range
        ax_density.plot(x_range, kde(x_range), color="black", linewidth=2, label="KDE")

    # Customize density subplot
    ax_density.set_xlabel("Total Counts", fontsize=10, fontweight="bold")
    ax_density.set_ylabel("Density", fontsize=10, fontweight="bold")
    ax_density.set_title(
        f"{medium} - Distribution", fontsize=12, fontweight="bold", pad=10
    )
    ax_density.grid(True, alpha=0.3, linestyle="--")
    ax_density.legend(fontsize=8)

    # Set consistent x-axis limits for all density plots
    ax_density.set_xlim(global_x_lim)

# Hide empty subplots
for i in range(n_mediums * 2, len(axes_flat)):
    axes_flat[i].set_visible(False)

plt.suptitle(
    "Cell Quality by Culture Medium with Density Distributions (Filtered: n_counts ≥ 100)",
    fontsize=16,
    fontweight="bold",
    y=0.98,
)
plt.tight_layout()
plt.show()

# Print detailed statistics for each medium
print(f"\n📊 Detailed Statistics by Culture Medium (After filtering: n_counts ≥ 100):")
for medium in culture_mediums_ordered:
    mask = adata.obs["CultureMedium"] == medium
    medium_data = adata.obs.loc[mask]

    print(f"\n🍽️ {medium}:")
    print(f"   📊 Total cells: {len(medium_data):,}")
    print(
        f"   📈 n_counts - Mean: {medium_data['n_counts'].mean():.0f}, Median: {medium_data['n_counts'].median():.0f}"
    )
    print(
        f"   📈 n_counts - Min: {medium_data['n_counts'].min():.0f}, Max: {medium_data['n_counts'].max():.0f}"
    )
    print(
        f"   🧬 n_genes - Mean: {medium_data['n_genes'].mean():.0f}, Median: {medium_data['n_genes'].median():.0f}"
    )
    print(
        f"   🧬 n_genes - Min: {medium_data['n_genes'].min():.0f}, Max: {medium_data['n_genes'].max():.0f}"
    )

    if len(medium_data) > 1:
        corr = np.corrcoef(medium_data["n_counts"], medium_data["n_genes"])[0, 1]
        print(f"   🔗 Correlation n_counts vs n_genes: {corr:.3f}")

# %%
