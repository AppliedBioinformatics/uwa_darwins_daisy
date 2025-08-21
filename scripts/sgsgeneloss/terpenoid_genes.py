import pandas as pd
from pathlib import Path
from popcolors import ordered_pop_colors
import seaborn as sns
import matplotlib.pyplot as plt

# Imports.
DATA_FOLDER = Path("../../data/")
NC_DATASET = DATA_FOLDER / "functional_annotation/noncore_go_merged_diamond_results_uniprot.csv"
PAV_MATRIX = DATA_FOLDER / "sgsgeneloss/pav_matrix.csv"
METADATA = Path("../../metadata/raw_sample_metadata.xlsx")


# Pandas-ify - indexes are the same.
nc_df = pd.read_csv(NC_DATASET, header=0, index_col=0)
pav_df = pd.read_csv(PAV_MATRIX, header=0, index_col=0)
meta_df = pd.read_excel(METADATA, header=0, index_col=0)

# Get indexes of GO term for "terpenoid biosynthesis" - GO:0016114
go_terms_to_include = [
    "GO:0016114",  # Terpenoid biosynthetic process
    "GO:0008299",  # Isoprenoid biosynthetic process
    "GO:0006721",  # Terpenoid metabolic process
    "GO:0006720",  # Isoprenoid metabolic process
    "GO:0120251",  # Hydrocarbon biosynthetic process
    "GO:0120252",  # Hydrocarbon metabolic process
    "GO:0006629",  # Lipid metabolic process
    "GO:0016102",  # Diterpenoid biosynthetic process
    "GO:0009699",  # Phenylpropanoid biosynthetic process
    "GO:0044550",  # Secondary metabolite biosynthetic process
]
pattern = "|".join(go_terms_to_include)
#filtered_df = nc_df[nc_df["go_terms"].str.contains(pattern, na=False)]
filtered_df = nc_df
i_filtered = filtered_df.index

filtered_pav = pav_df.loc[i_filtered]
filtered_pav.index = filtered_pav.index.map(nc_df["subject_id"])

# Plot heatmap.
# Sort columns by color/group (or by population name)
sorted_cols = [col for col in ordered_pop_colors.keys() if col in filtered_pav.columns]
filtered_pav = filtered_pav[sorted_cols]

# Plot simple heatmap (no clustering)
plt.figure(figsize=(12, 6))
ax = sns.heatmap(
    filtered_pav,
    cmap="viridis",
    cbar_kws={"label": "Presence/Absence"},
    xticklabels=True,
    yticklabels=False
)

# Add colored bars above columns
for tick_label in ax.get_xticklabels():
    sample = tick_label.get_text()
    if sample in ordered_pop_colors:
        tick_label.set_color(ordered_pop_colors[sample])

plt.ylabel("Genes / Features")
plt.xlabel("Samples grouped by Island")
plt.title("PAV Heatmap (Grouped by Population)")
plt.tight_layout()
plt.show()