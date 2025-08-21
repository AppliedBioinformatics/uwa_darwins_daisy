import pandas as pd
from tabulate import tabulate
from great_tables import GT


# Function to extract the uniprot accession from "subject_id".
def extract_accession(s):
    # If format is like "sp|F8S1H3|C7BL1_HELAN", extract middle part
    if '|' in s:
        parts = s.split('|')
        if len(parts) > 1:
            return parts[1]
    # Otherwise return the string as-is
    return s

# Import Non-core Go Merged gene list.
df = pd.read_csv("../../data/functional_annotation/noncore_go_merged_diamond_results_uniprot.csv")

# GO TERMS TO FILTER.
go_terms = ["GO:0016114"]
filtered_df = df[df["go_terms"].str.contains("GO:0016114", na=False)]

# Generate summary table.
summary = filtered_df.groupby("subject_id").agg(
    n_copies_in_gff=("subject_id", "count"),
    avg_bit_score=("bit_score", "mean"),
    reference_species=('sscinames', lambda x: ', '.join(sorted(set(x)))),
    description=('stitle', "first"),
).reset_index().sort_values("n_copies_in_gff", ascending=False)

# Ensure row index is integer-based (starting from 0)
summary = summary.reset_index(drop=True)
summary["subject_id"] = summary["subject_id"].str.split('|').str[-1]

# Add uniprot links.
summary["description"] = summary["description"].str.extract(r'^\S+\s(.+?)\sOS=')[0]
#summary['uniprot_link'] = summary['uniprot_accession'].apply(
#    lambda x: f'https://www.uniprot.org/uniprot/{x}'
#)

print(tabulate(summary, headers="keys", tablefmt="psql"))

gt = (
    GT(summary)
    .tab_header(
        title="Summary of GO:0016114 flagged genes.",
        subtitle="GO:0016114 - Terpenoid Biosynthetic Process."
    )
    .cols_label(
        subject_id="Uniprot Accession",
        n_copies_in_gff="Number of copies",
        avg_bit_score="Average bit score",
        reference_species="Reference species",
        description="Description",
    )
    .fmt_number(columns="avg_bit_score", decimals=0)
    .data_color(
        columns="avg_bit_score",
        palette=["white", "lightblue"],
        domain=[summary["avg_bit_score"].min(), summary["avg_bit_score"].max()],
        na_color="lightgray",
        autocolor_text=True,
        reverse=False,
    )
    .data_color(
        columns="n_copies_in_gff",
        palette=["white", "lightblue"],
        domain=[summary["n_copies_in_gff"].min(), summary["n_copies_in_gff"].max()],
        autocolor_text=True,
        reverse=False
    )

    .opt_table_font(font="Arial")
)

gt.show()

