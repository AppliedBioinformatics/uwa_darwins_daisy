import pandas as pd
from pathlib import Path

# This script merges taxonomy data with blastn results based on taxonomyID.
masurca_folder = Path("../../data/masurca")
blastn_file = masurca_folder / "final_blastn_results.tsv"
taxonomy_file = masurca_folder / "taxdata.tsv"

# Format blastn df.
blastn_headers = ["squeryid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle", "staxid", "NA1", "NA2"]
blastn_df = pd.read_csv(blastn_file, sep="\t", names=blastn_headers).iloc[:, :-2]

# Format taxdata df.
taxdata_headers = ["staxid", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]
taxdata_df = pd.read_csv(taxonomy_file, sep="\t", names=taxdata_headers)

# Inner Join on staxid.
blastn_df['staxid'] = blastn_df['staxid'].astype(str)
taxdata_df['staxid'] = taxdata_df['staxid'].astype(str)
merged_df = blastn_df.merge(taxdata_df, on="staxid", how="inner")

# Check all rows still present.
print(f"Total rows in merged dataframe: {len(merged_df)}")
strepto_rows = merged_df[merged_df['phylum'].str.lower() == 'streptophyta']
print(f"Number of rows with phylum 'Streptophyta': {len(strepto_rows)}")
