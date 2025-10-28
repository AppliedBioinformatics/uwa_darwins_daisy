import pandas as pd
from pathlib import Path
from geopy.distance import geodesic
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import ast

from scripts.sgsgeneloss.popcolors import pop_colors, island_colors

# Configuration & Globals
GO_TERMS = [
    # "GO:0016114",  # Terpenoid biosynthetic process
    # "GO:0008299",  # Isoprenoid biosynthetic process
    # "GO:0006721",  # Terpenoid metabolic process
    # "GO:0006720",  # Isoprenoid metabolic process
    # "GO:0120251",  # Hydrocarbon biosynthetic process
    # "GO:0120252",  # Hydrocarbon metabolic process
    # "GO:0006629",  # Lipid metabolic process
    # "GO:0016102",  # Diterpenoid biosynthetic process
    # "GO:0009699",  # Phenylpropanoid biosynthetic process
    # "GO:0044550",  # Secondary metabolite biosynthetic process
    #'"GO:0016114",
    #"GO:0008299",
    #"GO:0006721",
    #"GO:0006720",
    #"GO:0120251",
    # "GO:0120252",
    # "GO:0006629",
    # "GO:0016102",
    # "GO:0009699",
    # "GO:0044550",
    # "GO:0048544",
    # "GO:0016101",
    # "GO:0019748",
    # "GO:0008610",
    # "GO:0046246",
    # "GO:0008037",
    # "GO:0042214",
    # "GO:0009698",
    # "GO:0009805",
    # "GO:0031408",
    # "GO:0009804",
    # "GO:0006952",
    # "GO:0031407",
    # "GO:0005576",
    # "GO:0048046",
    # "GO:0004497",
    # "GO:0016491",
    # "GO:0046906",
    # "GO:0020037",
    # "GO:0005506",
    # "GO:0016705",
    # "GO:0003824",
    # "GO:0010333",
    # "GO:0016838",
    # "GO:0080043",
    # "GO:0030246",
    # "GO:0080044",
    # "GO:0106310",
    # "GO:0046527",
    # "GO:0003674",
    # "GO:0004674",
    # "GO:0035251",
    # "GO:0016712",
    # "GO:0004672",
    # "GO:0016709"
]

PAV_FILE = Path("../../data/sgsgeneloss/pav_matrix.csv")
DMND_FILE = Path("../../data/functional_annotation/noncore_go_merged_diamond_results_uniprot.csv")
META_FILE = Path("../../metadata/raw_sample_metadata.xlsx")
ORIGIN = (-0.252, -90.718)  # Santiago

def load_data():
    """Load PAV matrix, diamond results, and metadata."""
    pav_df = pd.read_csv(PAV_FILE, index_col=0)
    dmnd_df = pd.read_csv(DMND_FILE, index_col=0)
    dmnd_df.index = dmnd_df.index.str.replace(r'-mRNA-1$', '', regex=True)
    meta_df = pd.read_excel(META_FILE, index_col=0)
    return pav_df, dmnd_df, meta_df

def filter_dmnd_by_go_terms(dmnd_df, go_terms):
    """Filter diamond dataframe to only contain genes with specified GO terms."""
    dmnd_df = dmnd_df.copy()
    dmnd_df["go_terms"] = dmnd_df["go_terms"].apply(ast.literal_eval)

    # Uncomment below to filter by GO terms if needed:
    dmnd_df = dmnd_df[dmnd_df["go_terms"].apply(lambda x: any(term in go_terms for term in x))]
    return dmnd_df

def filter_pav_by_dmnd(pav_df, dmnd_df):
    """Filter PAV dataframe to only contain rows present in filtered diamond dataframe."""
    pav_df_filtered = pav_df.loc[dmnd_df.index]
    assert len(pav_df_filtered) == len(dmnd_df)
    return pav_df_filtered

def merge_with_metadata(pav_df_filtered, meta_df):
    """Summarize present genes per sample and merge with metadata."""
    present_genes_per_column = (pav_df_filtered == 1).sum(axis=0)
    present_summary = present_genes_per_column.reset_index()
    present_summary.columns = ["Sample", "Total_Present_Genes"]
    present_summary['Sample'] = present_summary['Sample'].str.replace('_merged_all', '', regex=False)
    merged_df = pd.merge(
        meta_df,
        present_summary,
        left_on='sampleID',
        right_on='Sample'
    )
    return merged_df

def compute_distance(row, origin):
    """Compute geodesic distance from origin to sample point."""
    sample_point = (row['Latitude'], row['Longitude'])
    return geodesic(origin, sample_point).kilometers

def add_distance_column(merged_df, origin):
    merged_df = merged_df.copy()
    merged_df["distance_km"] = merged_df.apply(lambda row: compute_distance(row, origin), axis=1)
    return merged_df

def plot_gene_count_vs_distance(merged_df):
    plt.figure(figsize=(10, 5))

    sns.scatterplot(
        x="distance_km",
        y="Total_Present_Genes",
        hue="Island",
        data=merged_df,
        palette=island_colors,
        s=80,
        edgecolor="black"
    )
    sns.regplot(
        x="distance_km",
        y="Total_Present_Genes",
        data=merged_df,
        scatter=False,
        color="black",
        line_kws={"lw": 2, "ls": "--"}
    )
    plt.title("GO:GO:0016114 - Terpenoid biosynthetic process")
    plt.xlabel("Distance from Origin - Santiago. (km)")
    plt.ylabel("Total Present Genes.")
    plt.legend(title="Island", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("../../plots/go_enrichment/top1_bp_ncore_presence_v_distance_frm_santiago.png")
    plt.show()

def run_regression_analysis(merged_df):
    X = sm.add_constant(merged_df["distance_km"])
    y = merged_df["Total_Present_Genes"]
    model = sm.OLS(y, X).fit()
    print("\n--- OLS Regression Results ---")
    print(model.summary())
    poisson_model = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    print("\n--- Poisson Regression Results ---")
    print(poisson_model.summary())

def main():
    pav_df, dmnd_df, meta_df = load_data()
    dmnd_df_filtered = filter_dmnd_by_go_terms(dmnd_df, GO_TERMS)
    pav_df_filtered = filter_pav_by_dmnd(pav_df, dmnd_df_filtered)
    merged_df = merge_with_metadata(pav_df_filtered, meta_df)
    merged_df = add_distance_column(merged_df, ORIGIN)
    plot_gene_count_vs_distance(merged_df)
    run_regression_analysis(merged_df)

if __name__ == "__main__":
    main()