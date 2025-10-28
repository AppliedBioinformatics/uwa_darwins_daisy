from pathlib import Path
from sklearn.metrics import pairwise_distances
from popcolors import order, pop_colors
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Helper funcs.
def match_similarity(u, v):
    return np.sum(u == v) / len(u)

# Load PAV.
DATA_FOLDER = Path("../../data/")
PAV_MATRIX = DATA_FOLDER / "sgsgeneloss_/pav_matrix.csv"

# To pandas + filter + transpose
pav_df = pd.read_csv(PAV_MATRIX, header=0, index_col=0)
nc_pav_df = pav_df[~(pav_df == 1).all(axis=1)]
nc_pav_t = nc_pav_df.T

# Compute similarity matrix
similarity_matrix = pd.DataFrame(
    pairwise_distances(nc_pav_t, metric=lambda u, v: 1 - match_similarity(u, v)),
    index=nc_pav_t.index,
    columns=nc_pav_t.index
)

similarity_matrix = 1 - similarity_matrix
similarity_matrix = similarity_matrix.loc[order, order]


# Visualise.
plt.figure(figsize=(12, 9))
ax = sns.heatmap(similarity_matrix, cmap='viridis', vmin=0.8, vmax=1)

# Color tick labels.
for tick_label in ax.get_xticklabels():
    sample = tick_label.get_text()
    tick_label.set_color(pop_colors[sample])

for tick_label in ax.get_yticklabels():
    sample = tick_label.get_text()
    tick_label.set_color(pop_colors[sample])

# Overlay sample of max similarity:
for i, sample in enumerate(similarity_matrix.index):
    row = similarity_matrix.loc[sample].drop(sample)   # exclude self
    partner = row.idxmax()
    j = similarity_matrix.columns.get_loc(partner)

    rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='red', linewidth=2)
    ax.add_patch(rect)

plt.title("Pairwise sample similarity")
plt.show()