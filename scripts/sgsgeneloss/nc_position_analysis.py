import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Function to extract chromosome and co-ordinates of each gene.
def build_position_table(excov_file:Path) -> pd.DataFrame:
    """Generates a pandas dataframe with rows as genes, and columns for start, end and length in BP."""
    df = pd.read_csv(excov_file, usecols=["ID", "chromosome", "start_position", "end_postion"])
    df = df.set_index("ID")
    df["length"] = df["end_postion"] - df["start_position"]
    return df

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def plot_non_core_genes(filtered, chrs_df, figsize=(10, 12), point_size=16):
    """
    Plot non-core gene positions along chromosomes (publication-ready).
    X-axis is scaled in base pairs (shown as Mbp).
    """

    # Sort chromosomes numerically by suffix
    sorted_chroms = sorted(chrs_df['chromosome'], key=lambda x: int(x.split("_")[1]))

    fig, ax = plt.subplots(figsize=figsize)

    for i, chrom in enumerate(sorted_chroms):
        chrom_len = chrs_df.loc[chrs_df['chromosome'] == chrom, 'end'].values[0]
        subdf = filtered[filtered['chromosome'] == chrom]

        # Chromosome backbone
        ax.hlines(
            y=i,
            xmin=0,
            xmax=chrom_len,
            color="black",
            linewidth=1.5,
            alpha=0.6,
            zorder=1
        )

        # Non-core gene positions
        ax.scatter(
            subdf['start_position'],
            [i] * len(subdf),
            s=point_size,
            alpha=0.7,
            color="#d62728",
            edgecolor="black",
            linewidth=0.2,
            zorder=2
        )

    # Y-axis formatting
    ax.set_yticks(range(len(sorted_chroms)))
    ax.set_yticklabels(sorted_chroms, fontsize=10)

    # X-axis formatting in Mbp
    ax.set_xlabel("Genomic position (Mbp)", fontsize=16)
    ax.set_ylabel("Chromosome", fontsize=16)
    ax.set_title("Distribution of Non-core Genes Along Chromosomes", fontsize=16, weight="bold")

    # Convert base pairs to Mbp on ticks
    def bp_to_mbp(x, pos):
        return f"{x/1e6:.0f}"

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(bp_to_mbp))
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    # Light grid for readability
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    return fig

def plot_combined_percentage_histogram(filtered, chrs_df, bin_size=5):
    """
    Plot a combined histogram of non-core gene density across all chromosomes,
    scaled to percentage chromosome length (publication-ready).
    """

    # Map chromosome lengths
    chr_len_dict = dict(zip(chrs_df['chromosome'], chrs_df['len']))

    # Compute percentage position for each gene
    filtered = filtered.copy()
    filtered['percent_pos'] = (
        filtered['start_position'] / filtered['chromosome'].map(chr_len_dict) * 100
    )

    # Define bins in % of chromosome length
    bins = np.arange(0, 100 + bin_size, bin_size)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(
        filtered['percent_pos'],
        bins=bins,
        color="#d62728",
        alpha=0.75,
        edgecolor="black",
        linewidth=0.6,
    )

    # Formatting
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_xlabel("Position along chromosome (% of length)", fontsize=12)
    ax.set_ylabel("Number of non-core genes", fontsize=12)
    ax.set_title(
        f"Normalised Non-core Gene Density Across Chromosomes\n(bin = {bin_size}%)",
        fontsize=14,
        weight="bold",
    )

    ax.tick_params(axis="both", labelsize=10)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    return fig

# Get PAV matrix and make core/non-core summary.
DATA_FOLDER = Path("../../data/sgsgeneloss/")
pav_df = pd.read_csv(DATA_FOLDER/ "pav_matrix.csv", index_col=0, header=0)
classification = pav_df.all(axis=1).map({True: "core", False: "non-core"})

# Build summary dataframe
summary_df = pd.DataFrame({
    "classification": classification
})

# Get contig (chromosome), start position and end position and relative chromosome end in one dataframe..
position_df = build_position_table( DATA_FOLDER/ "SGSGL_results/AHA6_30_merged_all.excov")
chrs_df = pd.read_csv(DATA_FOLDER/"SGSGL_results/chrs.csv", header=0)
chrs_df = chrs_df.rename(columns={"chr": "chromosome"})
chrs_df = chrs_df[chrs_df['chromosome'].str.startswith("chr")]
summary_df = summary_df.merge(position_df, left_index=True, right_index=True)
summary_df = summary_df.merge(chrs_df, on="chromosome", how="left")

# Filter for non-core genes and chromosomes only (not scaffolds).
filtered = summary_df[
    (summary_df['classification'] == 'non-core') &
    (summary_df['chromosome'].str.startswith("chr"))
]

# Sort chromosomes numerically
sorted_chroms = sorted(chrs_df['chromosome'], key=lambda x: int(x.split("_")[1]))

# Plot
plot_non_core_genes(filtered, chrs_df).show()
plot_combined_percentage_histogram(filtered, chrs_df, bin_size=0.5).show()
