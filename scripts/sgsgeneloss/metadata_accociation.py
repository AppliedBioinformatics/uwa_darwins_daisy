import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Optional
from tabulate import tabulate
from umap import UMAP
from pathlib import Path


# Load data/metadata.
pav_df = pd.read_csv("../../data/sgsgeneloss/pav_matrix.csv", index_col=0, header=0)
meta_df = pd.read_excel("../../metadata/raw_sample_metadata.xlsx", index_col=0, header=0)

# Do some formatting.
pav_df.columns = pav_df.columns.str.replace('_merged_all$', '', regex=True)
pav_df_t = pav_df.T
#meta_df = meta_df.iloc[:-1].drop(columns="comment")
meta_df = meta_df.loc[meta_df.index.intersection(pav_df_t.index)]

print(pav_df_t.shape)
print(meta_df.shape)

# Some plots to show association of samples between islands:
def plt_island_umap(pav_df_t: pd.DataFrame, meta_df: pd.DataFrame, outfile: Optional[str] = None):
    # Fit UMAP on binary PAV data using Jaccard distance
    umap_model = UMAP(n_components=2, random_state=42, metric="jaccard")
    umap_coords = umap_model.fit_transform(pav_df_t)

    # Create UMAP dataframe for plotting
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"], index=meta_df.index)
    umap_df["Island"] = meta_df["Island"]

    # Assign consistent colors to islands
    islands = umap_df["Island"].unique()
    colors = plt.cm.tab10.colors  # or use another colormap like plt.cm.Set2.colors
    color_map = {island: colors[i % len(colors)] for i, island in enumerate(islands)}

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for island in islands:
        subset = umap_df[umap_df["Island"] == island]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            label=island,
            color=color_map[island],
            edgecolor='black',
            s=60,
            alpha=0.8
        )

    # Formatting
    ax.set_title("UMAP of PAV Matrix", fontsize=14)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    ax.legend(title="Island", loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.tight_layout()

    if outfile:
        outfile_path = Path(outfile)
        plt.savefig(outfile_path)

    return plt


def plt_island_umap_3_dims(pav_df_t: pd.DataFrame, meta_df: pd.DataFrame, outfile: Optional[str] = None):
    umap_model = UMAP(n_components=3, random_state=42, metric="hamming", n_neighbors=7)
    umap_coords = umap_model.fit_transform(pav_df_t)

    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2", "UMAP3"], index=meta_df.index)
    umap_df["Island"] = meta_df["Island"]

    islands = umap_df["Island"].unique()
    colors = plt.cm.tab10.colors
    color_map = {island: colors[i % len(colors)] for i, island in enumerate(islands)}

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for island in islands:
        subset = umap_df[umap_df["Island"] == island]

        z = subset["UMAP3"]
        min_size, max_size = 20, 100
        sizes = min_size + (z - z.min()) / (z.max() - z.min()) * (max_size - min_size)

        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            subset["UMAP3"],
            color=color_map[island],
            edgecolor='black',
            s=sizes,
            alpha=0.8
        )

    ax.set_xlim(umap_df["UMAP1"].min(), umap_df["UMAP1"].max())
    ax.set_ylim(umap_df["UMAP2"].min(), umap_df["UMAP2"].max())
    ax.set_zlim(umap_df["UMAP3"].min(), umap_df["UMAP3"].max())
    ax.margins(0)

    ax.set_title("3D UMAP of PAV Matrix", fontsize=14)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_zlabel("UMAP3")

    plt.tight_layout()

    if outfile:
        outfile_path = Path(outfile)
        plt.savefig(outfile_path, bbox_inches='tight')

    return plt

def plt_geo_map(outfile: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-91.75, -89, -1.5, 1], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.5)
    ax.gridlines(draw_labels=True)

    # OPTIONAL: Color by 'Island' if it exists
    if "Island" in meta_df.columns:
        islands = meta_df["Island"].unique()
        colors = plt.cm.tab10.colors  # Up to 10 colors
        color_map = {island: colors[i % len(colors)] for i, island in enumerate(islands)}

        for island in islands:
            subset = meta_df[meta_df["Island"] == island]
            ax.scatter(
                subset["Longitude"], subset["Latitude"],
                transform=ccrs.PlateCarree(),
                label=island,
                color=color_map[island],
                edgecolor="black",
                s=60,
                alpha=0.8
            )
    else:
        # Just plot all samples as red if no group
        ax.scatter(
            meta_df["Longitude"], meta_df["Latitude"],
            transform=ccrs.PlateCarree(),
            color="red",
            edgecolor="black",
            s=60,
            label="Samples"
        )

    # Move legend outside the plot
    ax.legend(
        title="Island",
        loc='center left',
        bbox_to_anchor=(1.25, 0.5),
        borderaxespad=0
    )

    # Title and layout
    ax.set_title("Galápagos Islands with Sample Locations", fontsize=14)
    plt.tight_layout()

    if outfile:
        outfile_path = Path(outfile)
        plt.savefig(outfile_path)

def filter_pavs_by_prevalance(pav_df: pd.DataFrame, min: float, max: float):
    df = pav_df.copy()
    prevalence = df.mean(axis=1)
    filtered_genes = prevalence[(prevalence >= min) & (prevalence <= max)].index
    subset = df.loc[filtered_genes]
    print(f"Prevalence subset contains {len(subset)} samples")
    print(tabulate(subset, headers="keys", tablefmt="psql"))
    return subset

def plt_non_core_heatmap(pav_df: pd.DataFrame, meta_df: pd.DataFrame, outfile: Optional[str] = None) -> None:
    df = pav_df.copy()
    mdf = meta_df.copy()
    nc_gene_mask = df.columns[(pav_df != 1).any(axis=0)]
    subset = df[nc_gene_mask]

    # Sort Columns.
    sorted_samples = mdf.sort_values("Island").index.tolist()
    subset = subset[sorted_samples]

    # Build Columns colors
    islands = mdf["Island"].unique()
    colors = plt.cm.tab10.colors
    color_map = {island: colors[i % len(colors)] for i, island in enumerate(islands)}
    col_colors = mdf["Island"].map(color_map)

    sns.clustermap(subset,
                col_cluster=False,
                row_cluster=False,
                cmap="Greys",
                figsize=(20, 20),
                cbar_kws={"label": "Presence (1) / Absence (0)"},
                col_colors=col_colors,
                )
    plt.title("PAV Heatmap of Non-Core genes.")
    plt.xlabel("Features")
    plt.ylabel("Samples")
    plt.show()

def plt_umap_with_color(
    pav_df_t: pd.DataFrame,
    meta_df: pd.DataFrame,
    color_col: str,
    outfile: Optional[str] = None
    ):
    """
    Create a UMAP plot of binary PAV data, coloring points based on a specified column in meta_df.

    Parameters
    ----------
    pav_df_t : pd.DataFrame
        Binary PAV data (samples as rows, features as columns).
    meta_df : pd.DataFrame
        Metadata DataFrame with an index matching pav_df_t.
    color_col : str
        Column name in meta_df to use for coloring the points.
    outfile : str, optional
        If provided, saves the plot to this path.
    """
    # Fit UMAP on binary PAV data using Jaccard distance
    umap_model = UMAP(n_components=2, random_state=42, metric="hamming", n_neighbors=7)
    umap_coords = umap_model.fit_transform(pav_df_t)

    # Create UMAP dataframe for plotting
    umap_df = pd.DataFrame(umap_coords, columns=["UMAP1", "UMAP2"], index=meta_df.index)
    umap_df[color_col] = meta_df[color_col]

    # Assign consistent colors to unique groups
    groups = umap_df[color_col].unique()
    colors = plt.cm.tab10.colors  # or another colormap like plt.cm.Set2.colors
    color_map = {group: colors[i % len(colors)] for i, group in enumerate(groups)}

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for group in groups:
        subset = umap_df[umap_df[color_col] == group]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            label=group,
            color=color_map[group],
            edgecolor='none',
            s=60,
            alpha=0.8
        )

    # Formatting
    ax.set_title("UMAP of PAV Matrix", fontsize=14)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    ax.legend(title=color_col, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.tight_layout()

    if outfile:
        outfile_path = Path(outfile)
        plt.savefig(outfile_path)

    return plt

if __name__ == "__main__":
    plt_island_umap_3_dims(pav_df_t=pav_df_t, meta_df=meta_df, outfile="../../plots/sgsgeneloss/umap_3d.png" )
    plt_umap_with_color(pav_df_t, meta_df, color_col="Island", outfile="../../plots/sgsgeneloss/umap_2d_ff.png")

    #plt_non_core_heatmap(pav_df, meta_df)