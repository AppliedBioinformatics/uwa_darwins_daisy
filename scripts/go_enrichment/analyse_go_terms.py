import pandas as pd
from collections import Counter
from pathlib import Path
import ast
from itertools import chain
import matplotlib.pyplot as plt
import seaborn as sns

# GLOBALS
DATA_DIR = Path("../../data/functional_annotation/")
GO_DIAMOND_DF_PATH = Path(DATA_DIR / "go_merged_diamond_results_uniprot.tsv")
GO_NAMES_MAP_DF_PATH = Path(DATA_DIR / "go_terms_human_readable.tsv")
GO_HR_PATH = Path(DATA_DIR / "goterms_human_readable.tsv")
PAV_DF_PATH = Path("../../data/sgsgeneloss/pav_matrix.csv")

# Read and format Go-ID : Go-Names data.
go_map_df =pd.read_csv(GO_HR_PATH, sep="\t")
go_id_to_name = dict(zip(go_map_df["GO_ID"], go_map_df["Name"]))
go_id_to_ns = dict(zip(go_map_df["GO_ID"], go_map_df["Namespace"]))

# Read and format diamond results with GO terms.
god_df = pd.read_csv(GO_DIAMOND_DF_PATH, sep="\t", header=0)
god_df["go_terms"] = god_df["go_terms"].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x
    )

# Generate a count table for each Go term:
def get_top_n_go_terms_df(df: pd.DataFrame) -> pd.DataFrame:

    # Flatten and count freq.
    all_go_terms = list(chain.from_iterable(df["go_terms"]))
    go_counts = Counter(all_go_terms)
    go_freq_df = pd.DataFrame(go_counts.items(), columns=["term", "count"])
    go_freq_df["go_name"] = go_freq_df["term"].map(go_id_to_name).fillna("Unknown")
    go_freq_df["Namespace"] = go_freq_df["term"].map(go_id_to_ns).fillna("Unknown")
    count_df = go_freq_df.sort_values("count", ascending=False)

    return count_df


# Visualisation.
def truncate_label(label, max_length=30):
    return label if len(label) <= max_length else label[:max_length] + "..."

def plt_top_n_go_terms(df: pd.DataFrame, n: int) -> plt.Figure:
    _df = df.copy()
    _df['short_name'] = _df['go_name'].apply(lambda x: truncate_label(x))

    plt.figure(figsize=(10, 6))
    sns.barplot(data=_df.head(n), x="count", y="short_name", palette="viridis")
    plt.title(f"Top {n} Most Frequent GO Terms - Core genes.")
    plt.xlabel("Frequency")
    plt.ylabel("GO Term (Human-Readable)")
    plt.tight_layout()
    plt.savefig("../../plots/functional_annotation/diamond_hits_core_most_freq_go_terms.png", dpi=300, bbox_inches="tight")

def plt_top_n_go_terms_by_namespace(df: pd.DataFrame) -> plt.Figure:

    _df = df.copy()
    _df['short_name'] = _df['go_name'].apply(lambda x: truncate_label(x))

    top_n_per_ns = (
        _df.groupby("Namespace", group_keys=False)
        .apply(lambda g: g.sort_values("count", ascending=False).head(10))
    )

    namespace_label_map = {
        "biological_process": "Biological Process",
        "molecular_function": "Molecular Function",
        "cellular_component": "Cellular Component",
        "Unknown": "Unknown"
    }

    top_n_per_ns["Namespace"] = top_n_per_ns["Namespace"].map(namespace_label_map)
    plt.figure(figsize=(12, 8))
    sns.barplot(data=top_n_per_ns, x="count", y="short_name", hue="Namespace", dodge=False)
    plt.title("Top 10 GO Terms per Namespace - Core genes.")
    plt.xlabel("Frequency")
    plt.ylabel("GO Term (Human-Readable)")
    plt.legend(title="GO Namespace")
    plt.tight_layout()
    plt.savefig("../../plots/functional_annotation/diamond_hits_core_most_freq_go_terms_by_namespace.png", dpi=300, bbox_inches="tight")

if __name__=="__main__":

    # 1) Get a list of core genes and a list of non-core genes from the pav dataframe.
    pav_df = pd.read_csv(PAV_DF_PATH, index_col=0)
    core_ids = pav_df.index[pav_df.sum(axis=1) == pav_df.shape[1]]
    ncore_ids = pav_df.index[pav_df.sum(axis=1) < pav_df.shape[1]]

    # 2) Filter the go_diamond dataframe using those ID's to create a dataframe of core and
    # non-core genes.
    god_df["query_id"] = god_df["query_id"].str.replace("-mRNA-1$", "", regex=True)
    core_god_df = god_df[god_df["query_id"].isin(core_ids)]
    ncore_god_df = god_df[god_df["query_id"].isin(ncore_ids)]

    # 3) Pass those back into the visualisation plots above to view the top GO scores for present and abcent genes.
    core_god_df.to_csv(DATA_DIR / "core_go_merged_diamond_results_uniprot.csv", index=False)
    ncore_god_df.to_csv(DATA_DIR / "noncore_go_merged_diamond_results_uniprot.csv", index=False)

