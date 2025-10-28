# Builds a GO network using networkx to visualise enriched GO terms.
from pathlib import Path
from go_enrichment_.go_network import build_go_network, plot_go_network
import pandas as pd

if __name__ == "__main__":
    DATA_FOLDER = Path("../../data/functional_annotation")
    go_df = pd.read_csv(DATA_FOLDER/"goatools_go_enrichment_results_overrep.csv", header=0)

    # Build plot.
    network_x_graph = build_go_network(go_df=go_df)
    plt = plot_go_network(graph=network_x_graph, go_df=go_df, plot_labels=True, pvalue=0.01)
    plt.show()


