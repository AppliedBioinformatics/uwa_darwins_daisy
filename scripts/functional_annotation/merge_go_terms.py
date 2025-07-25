import obonet
import pandas as pd
from pathlib import Path
from annotation_analysis import DIAMOND_TSV_HEADERS
import csv

# Data dir
data_dir = Path("../../data/functional_annotation")

# Load diamond results.
diamond_df = pd.read_csv(data_dir / "diamond_results_protein.tsv", sep="\t", header=None)
diamond_df.columns = DIAMOND_TSV_HEADERS
diamond_df["sseqid"] = diamond_df["subject_id"].str.extract(r'\|([A-Z0-9]+)\|')

# Load GO file.
go_df = pd.read_csv(data_dir / "go_mapping.tsv", sep="\t", header=None, names= ["sseqid", "go_raw"])
go_df["go_terms"] = go_df["go_raw"].fillna("").apply(lambda x: [term.strip() for term in x.split(";")] if x else [])

# Merge Diamond results with the GO term table.
merged_df = diamond_df.merge(go_df[["sseqid", "go_terms"]], on="sseqid", how="left")
merged_df["go_terms"] = merged_df["go_terms"].apply(lambda x: x if isinstance(x, list) else [])

merged_df.to_csv(data_dir / "go_merged_diamond_results_uniprot.tsv", sep="\t", index=False)

# --- Generate GoID: Human readable mapping file.
graph = obonet.read_obo(data_dir / "go-basic.obo")
with open(data_dir / "go_terms_human_readable.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["GO_ID", "Name", "Namespace"])
    for go_id, data in graph.nodes(data=True):
        if go_id.startswith("GO:"):
            w.writerow([go_id, data.get("name", ""), data.get("namespace", "")])