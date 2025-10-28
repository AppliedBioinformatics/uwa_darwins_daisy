import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from tabulate import tabulate
from go_enrichment_.goatools_enrichment_study import *

if __name__ == "__main__":
    # Globals.
    DATA_FOLDER = Path("../../data/functional_annotation")
    BG_DATASET = DATA_FOLDER / "go_merged_diamond_results_uniprot.tsv"
    NC_DATASET = DATA_FOLDER / "noncore_go_merged_diamond_results_uniprot.csv"
    OBODAG = DATA_FOLDER / "go-basic.obo"

    # Pandas-ify.
    bg_df = pd.read_csv(BG_DATASET, sep="\t", header=0)
    nc_df = pd.read_csv(NC_DATASET, header=0)

    # Run study.
    results_df = run_goatools_enrichment_study(bg_set=bg_df, study_set=nc_df, go_obo=OBODAG)

    # Get enriched significant results.
    over_df = results_df.query("p_fdr_bh < 0.01 and enrichment == 'e'")

    # Build wordcloud.
    wc_fig = plt_wordcloud(over_df, max_words=45, subset="e")
    wc_fig.show()
    wc_fig.savefig("../../plots/wordcloud.png")

    # Plot top n by namespace figure.
    fig = plt_top_n_go_enriched_terms_by_namespace(over_df)
    fig.show()

    #over_df.to_csv(DATA_FOLDER / "goatools_go_enrichment_results_overrep.csv", index=False)
    #under_df = sig_df[sig_df['enrichment'] == 'p'].copy()
    #under_df.to_csv(DATA_FOLDER / "goatools_go_enrichment_results_underrep.csv", index=False)

    # Format for visualisation
    #under_df["-log10(FDR)"] = -np.log10(sig_df["p_fdr_bh"])
    #under_df["short_name"] = under_df["name"].apply(lambda x: x if len(x) < 40 else x[:37] + "...")
    #under_df = under_df.sort_values("-log10(FDR)", ascending=False)