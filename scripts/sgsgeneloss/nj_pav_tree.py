import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj
from pathlib import Path
from tabulate import tabulate

# Read data.
pav_file = Path("../../data/sgsgeneloss/pav_matrix.csv")
pav_df = pd.read_csv(pav_file, header=0, index_col=0)
pav_df.columns = pav_df.columns.str.replace('_merged_all', '', regex=False)
pav_df_T = pav_df.T

# Compute hamming distance.
dist_array = pdist(pav_df_T, 'hamming')
dist_matrix = squareform(dist_array)

# Build distance matrix with sample names.
sample_names = pav_df_T.index.tolist()
dm = DistanceMatrix(dist_matrix, sample_names)

# Build NJ-tree.
tree = nj(dm)
tree = tree.root_at_midpoint()
tree.write("../../data/phylo_pav_tree.nwk")
print(tabulate(pav_df.head(), headers="keys"))