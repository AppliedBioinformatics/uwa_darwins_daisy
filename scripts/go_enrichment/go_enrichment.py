import pandas as pd
import ast
from pathlib import Path
from goatools.obo_parser import GODag
import numpy as np
from goatools.go_enrichment import GOEnrichmentStudy
import seaborn as sns
import matplotlib.pyplot as plt

# GLOBALS
DATA_FOLDER = Path("../../data/functional_annotation")
BG_DATASET = DATA_FOLDER / "go_merged_diamond_results_uniprot.tsv"
NC_DATASET = DATA_FOLDER / "noncore_go_merged_diamond_results_uniprot.csv"

# Load background dataset and non-core dataset as pandas dataframes.
bg_df = pd.read_csv(BG_DATASET, sep="\t", header=0)
nc_df = pd.read_csv(NC_DATASET, header=0)

unique_subject_ids = nc_df["subject_id"].nunique()
print(f"Number of unique subject_id values: {unique_subject_ids}")


def get_gene2go(df: pd.DataFrame) -> dict:
    """Convert a DataFrame with 'subject_id' and 'go_terms' columns into a gene2go dict."""
    gene2go = {}

    for _, row in df.iterrows():
        # Use the UniProt accession from the subject_id (e.g. sp|Q8RWD5|SCD2_ARATH → Q8RWD5)
        parts = row["subject_id"].split("|")
        gene_id = parts[1]

        # Safely parse go_terms if stored as a string
        go_list = row["go_terms"]
        if isinstance(go_list, str):
            try:
                go_list = ast.literal_eval(go_list)
            except (SyntaxError, ValueError):
                go_list = []

        if isinstance(go_list, list) and go_list:
            gene2go[gene_id] = set(go_list)

    return gene2go

def get_enrichment_results_df(go_results: list) -> pd.DataFrame:
    df_results = pd.DataFrame([
        {
            "GO": r.GO,
            "name": r.name,
            "namespace": r.NS,
            "study_count": r.study_count,
            "study_n": r.study_n,
            "pop_count": r.pop_count,
            "pop_n": r.pop_n,
            "p_uncorrected": r.p_uncorrected,
            "p_fdr_bh": r.p_fdr_bh
        }
        for r in results
    ])

    return df_results

def truncate_label(label: str, max_length: int = 40) -> str:
    return label if len(label) <= max_length else label[:max_length - 3] + "..."

def plt_top_n_go_enriched_terms_by_namespace(df: pd.DataFrame, n: int = 15) -> plt.Figure:
    """
    Plot top N enriched GO terms per namespace using -log10(FDR) from GOATOOLS results.
    Expects columns: ['name', 'namespace', 'p_fdr_bh']
    """

    _df = df.copy()
    _df["-log10(FDR)"] = -np.log10(_df["p_fdr_bh"])
    _df["short_name"] = _df["name"].apply(truncate_label)

    # Get top N terms per namespace
    top_n_per_ns = (
        _df.groupby("namespace", group_keys=False)
        .apply(lambda g: g.sort_values("p_fdr_bh").head(n))
    )

    # Rename namespaces for plot
    namespace_label_map = {
        "biological_process": "Biological Process",
        "molecular_function": "Molecular Function",
        "cellular_component": "Cellular Component",
        "BP": "Biological Process",
        "MF": "Molecular Function",
        "CC": "Cellular Component",
        "Unknown": "Unknown"
    }
    top_n_per_ns["Namespace"] = top_n_per_ns["namespace"].map(namespace_label_map)

    # Plot
    plt.figure(figsize=(12, 10))
    sns.barplot(
        data=top_n_per_ns,
        x="-log10(FDR)",
        y="short_name",
        hue="Namespace",
        dodge=False
    )

    plt.title(f"Top {n} Enriched GO Terms per Namespace (Absent Genes)")
    plt.xlabel("-log10(FDR)")
    plt.ylabel("GO Term")
    plt.legend(title="GO Namespace")
    plt.tight_layout()
    return plt.gcf()


if __name__ == "__main__":
    bg_gene2go = get_gene2go(bg_df)
    nc_gene2go = get_gene2go(nc_df)

    # For GOatools.
    pop_genes = list(bg_gene2go.keys())
    study_genes = list(nc_gene2go.keys())
    obodag = GODag(DATA_FOLDER / "go-basic.obo")

    # Run analysis:
    goea = GOEnrichmentStudy(
        pop = pop_genes,  # background genes
        assoc=bg_gene2go,  # gene2go mapping
        obo_dag=obodag,  # GO DAG
        methods=["fdr_bh"],  # adjust p-values with Benjamini-Hochberg
        alpha=0.01,  # significance threshold
        propagate_counts=True
    )

    results = goea.run_study(study_genes)
    results_df = get_enrichment_results_df(results)
    sig_df = results_df[results_df["p_fdr_bh"] < 0.05].copy()
    sig_df.to_csv(DATA_FOLDER / "goatools_go_enrichment_results.csv", index=False)

    # Format for visualisation
    sig_df["-log10(FDR)"] = -np.log10(sig_df["p_fdr_bh"])
    sig_df["short_name"] = sig_df["name"].apply(lambda x: x if len(x) < 40 else x[:37] + "...")

