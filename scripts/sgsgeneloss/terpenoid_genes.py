import pandas as pd
from pathlib import Path
from popcolors import ordered_pop_colors
import seaborn as sns
import matplotlib.pyplot as plt
import math
import matplotlib.colors as mcolors
import numpy as np

from scripts.functional_annotation.combine_figures import grid_cols

# Imports.
DATA_FOLDER = Path("../../data/")
NC_DATASET = DATA_FOLDER / "functional_annotation/noncore_go_merged_diamond_results_uniprot.csv"
PAV_MATRIX = DATA_FOLDER / "sgsgeneloss/pav_matrix.csv"
LOC_DATA = DATA_FOLDER / "sgsgeneloss/gene_length_table.csv"
METADATA = Path("../../metadata/raw_sample_metadata.xlsx")


# Pandas-ify - indexes are the same.
nc_df = pd.read_csv(NC_DATASET, header=0, index_col=0)
pav_df = pd.read_csv(PAV_MATRIX, header=0, index_col=0)
loc_df = pd.read_csv(LOC_DATA, header=0, index_col=0)
meta_df = pd.read_excel(METADATA, header=0, index_col=0)

# Get indexes of GO term for "terpenoid biosynthesis" - GO:0016114
go_terms_to_include = [
    "GO:0016114",  # Terpenoid biosynthetic process
    "GO:0008299",  # Isoprenoid biosynthetic process
    "GO:0006721",  # Terpenoid metabolic process
    "GO:0006720",  # Isoprenoid metabolic process
    "GO:0120251",  # Hydrocarbon biosynthetic process
    #"GO:0120252",  # Hydrocarbon metabolic process
    #"GO:0006629",  # Lipid metabolic process
    #"GO:0016102",  # Diterpenoid biosynthetic process
    #"GO:0009699",  # Phenylpropanoid biosynthetic process
    #"GO:0044550",  # Secondary metabolite biosynthetic process
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