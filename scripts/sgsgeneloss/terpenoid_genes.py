import pandas as pd
from pathlib import Path
from popcolors import ordered_pop_colors
import seaborn as sns
import matplotlib.pyplot as plt
import math
import matplotlib.colors as mcolors
import numpy as np
from tabulate import tabulate
from scripts.functional_annotation.combine_figures import grid_cols

# Imports.
DATA_FOLDER = Path("../../data/")
NC_DATASET = DATA_FOLDER / "functional_annotation/noncore_go_merged_diamond_results_uniprot.csv"
PAV_MATRIX = DATA_FOLDER / "sgsgeneloss_/pav_matrix.csv"
LOC_DATA = DATA_FOLDER / "sgsgeneloss_/gene_length_table.csv"
METADATA = Path("../../metadata/raw_sample_metadata.xlsx")


# Pandas-ify - indexes are the same.
nc_df = pd.read_csv(NC_DATASET, header=0, index_col=0)
pav_df = pd.read_csv(PAV_MATRIX, header=0, index_col=0)
loc_df = pd.read_csv(LOC_DATA, header=0, index_col=0)
meta_df = pd.read_excel(METADATA, header=0, index_col=0)

# Get indexes of GO term for "terpenoid biosynthesis" - GO:0016114
go_terms_to_include = [
    "GO:0016114",
    "GO:0008299",
    "GO:0006721",
    "GO:0006720",
    "GO:0120251",
    "GO:0120252",
    "GO:0006629",
    "GO:0016102",
    "GO:0009699",
    "GO:0044550",
    "GO:0048544",
    "GO:0016101",
    "GO:0019748",
    "GO:0008610",
    "GO:0046246",
    "GO:0008037",
    "GO:0042214",
    "GO:0009698",
    "GO:0009805",
    "GO:0031408",
    "GO:0009804",
    "GO:0006952",
    "GO:0031407",
    "GO:0005576",
    "GO:0048046",
    "GO:0004497",
    "GO:0016491",
    "GO:0046906",
    "GO:0020037",
    "GO:0005506",
    "GO:0016705",
    "GO:0003824",
    "GO:0010333",
    "GO:0016838",
    "GO:0080043",
    "GO:0030246",
    "GO:0080044",
    "GO:0106310",
    "GO:0046527",
    "GO:0003674",
    "GO:0004674",
    "GO:0035251",
    "GO:0016712",
    "GO:0004672",
    "GO:0016709"
]
pattern = "|".join(go_terms_to_include)


# Filter.
filtered_df = nc_df[nc_df["go_terms"].str.contains(pattern, na=False)]
#filtered_df = nc_df
i_filtered = filtered_df.index
filtered_pav = pav_df.loc[i_filtered]

# To plot ALL genes;
#filtered_pav = pav_df

# Sort by co-ordinate using start position of each feature.
sorted_pav_df = (filtered_pav.merge(loc_df[["start_position"]], left_index=True, right_index=True)
                        .sort_values("start_position")
                       )

# To plot ALL genes:
#functional_pav_df = sorted_pav_df


# Re annotate index using functional annotation.
functional_pav_df = (sorted_pav_df.merge(nc_df[["subject_id"]], left_index=True, right_index=True))
functional_pav_df = functional_pav_df.set_index("subject_id", drop=True)
functional_pav_df.index = functional_pav_df.index.str.split("|", n=2).str[-1]

# Make index start position.
#functional_pav_df = functional_pav_df.set_index("start_position", drop=True)


# Plot heatmap.
# Sort columns by color/group (or by population name)
sorted_cols = [col for col in ordered_pop_colors.keys() if col in functional_pav_df.columns]
functional_pav_df = functional_pav_df[sorted_cols]

# Color bars for legend.
cmap = mcolors.ListedColormap(["white", "darkgreen"])  # absent=white, present=darkgreen
bounds = [0, 0.5, 1]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# --- Figure setup ---
plt.figure(figsize=(16, 18))
ax = sns.heatmap(
    functional_pav_df,
    cmap=cmap,
    norm=norm,
    cbar_kws={
        "label": "Presence / Absence",
        "ticks": [0.25, 0.75],
        "orientation": "vertical",
        "shrink": 0.7
    },
    linewidths=1,
    linecolor='white',
    xticklabels=True,
    yticklabels=True,
)

# --- Y-axis ticks (rounded & human-readable) ---
#num_ticks = 10
#start = functional_pav_df.index.min()
#end = functional_pav_df.index.max()
#tick_positions = np.linspace(0, functional_pav_df.shape[0]-1, num_ticks)

#def human_format(num):
#    """Format large numbers with k/M/B suffixes."""
#    for unit in ['', 'k', 'M', 'B']:
#        if abs(num) < 1000:
#            return f"{int(num)}{unit}"
#        num /= 1000.0
#    return f"{num:.1f}B"

#tick_labels_raw = np.linspace(math.floor(start/1000)*1000, math.ceil(end/1000)*1000, num_ticks, dtype=int)
#tick_labels = [human_format(x) for x in tick_labels_raw]

#ax.set_yticks(tick_positions)
#ax.set_yticklabels(tick_labels, fontsize=20)
#ax.set_ylabel("Position in Genome", fontsize=20)

# --- X-axis customization: color labels by population ---
for tick_label in ax.get_xticklabels():
    sample = tick_label.get_text()
    if sample in ordered_pop_colors:
        tick_label.set_color(ordered_pop_colors[sample])
    tick_label.set_rotation(90)
    tick_label.set_fontsize(18)

#ax.set_xlabel("Sample (grouped by island)", fontsize=20)
plt.title("PAV Heatmap (Grouped by Population)", fontsize=24, pad=20)


# --- Colorbar tweaks ---
cbar = ax.collections[0].colorbar
cbar.set_ticklabels(["Absent", "Present"])
cbar.ax.tick_params(labelsize=14)
cbar.set_label("Presence / Absence", fontsize=20)

# --- Tight layout for publication ---
ax.set_xticklabels(ax.get_xticklabels(), fontsize=20)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=20)
plt.tight_layout()
plt.show()

# Filter by GO then look at metadata associations.
def plt_asccociation_habitat(pav_df: pd.DataFrame, meta_df: pd.DataFrame) -> None:
    df = pav_df.copy()
    meta_df = meta_df.copy()

    # Calculate number of present genes per sample
    df["n_present_genes"] = df.sum(axis=1)

    # Merge with metadata to get habitat
    merged_df = df.merge(meta_df[["Climate"]], left_index=True, right_index=True)

    # Print mean number of present genes per habitat
    mean_per_habitat = merged_df.groupby("Climate")["n_present_genes"].mean()
    print("Average number of present genes per habitat:")
    print(mean_per_habitat)

    # Plot boxplot
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=merged_df, x="Climate", y="n_present_genes")
    sns.stripplot(data=merged_df, x="Climate", y="n_present_genes",
                  color="black", alpha=0.5)
    plt.ylabel("Number of Non-Core Genes Present")
    plt.title("Gene Counts by Climate")
    plt.show()


sorted_pav_df.drop("start_position", inplace=True, axis=1)
plt_asccociation_habitat(sorted_pav_df.T, meta_df)