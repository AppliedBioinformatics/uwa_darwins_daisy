import pandas as pd
from scripts.go_enrichment.non_core_presence_v_distance import filter_pav_by_dmnd, filter_dmnd_by_go_terms
from tabulate import tabulate
from great_tables import GT, style, loc
from popcolors import pop_colors

# GO Terms to filter by:
GO_TERMS = [
    "GO:0016114",  # Terpenoid biosynthetic process
    "GO:0008299",  # Isoprenoid biosynthetic process
    "GO:0006721",  # Terpenoid metabolic process
    "GO:0006720",  # Isoprenoid metabolic process
    "GO:0120251",  # Hydrocarbon biosynthetic process
    "GO:0120252",  # Hydrocarbon metabolic process
    "GO:0006629",  # Lipid metabolic process
    "GO:0016102",  # Diterpenoid biosynthetic process
    "GO:0009699",  # Phenylpropanoid biosynthetic process
    "GO:0044550",  # Secondary metabolite biosynthetic process
]

# Read and transpose.
pav_df = pd.read_csv("../../data/sgsgeneloss/pav_matrix.csv", header=0, index_col=0)
dmnd_df = pd.read_csv("../../data/functional_annotation/noncore_go_merged_diamond_results_uniprot.csv", index_col=0)
nc_pav_df = pav_df[(pav_df == 0).any(axis=1)]
nc_pav_df_t = nc_pav_df.T

# Build summary table.
nc_summary_df = nc_pav_df_t.sum(axis=1)
nc_summary_df = nc_summary_df.to_frame(name="All non core (n = 1988).")


# Filter by Top 10 BP GO term.
dmnd_filtered_df = filter_dmnd_by_go_terms(dmnd_df, GO_TERMS)
pav_filtered_df = filter_pav_by_dmnd(pav_df, dmnd_filtered_df)

filtered_summary_df = pav_filtered_df.T.sum(axis=1)
filtered_summary_df = filtered_summary_df.to_frame(name="Top 10 enriched GO:BP's (n = 91).")

combined_summary_df = nc_summary_df.merge(filtered_summary_df, left_index=True, right_index=True)
combined_summary_df = combined_summary_df.reset_index().rename(columns={"index": "Sample"})

# Great table-ify
gf_table = (
            GT(combined_summary_df)
            .tab_header(
                title="Gene Presence of Non-core features.",
                subtitle="Number of genes present in all non core and subset of top 10 enriched GO biological processes."
            )
            .tab_spanner(
                label="Gene Counts,",
                columns=["All non core (n = 1988).", "Top 10 enriched GO:BP's (n = 91)."]
            )
            .data_color(
                columns="All non core (n = 1988).",
                palette=["white", "lightblue"],
                domain=[combined_summary_df["All non core (n = 1988)."].min(), combined_summary_df["All non core (n = 1988)."].max()],
                na_color="lightgray",
                autocolor_text=True,
                reverse=False,
            )
            .data_color(
                columns="Top 10 enriched GO:BP's (n = 91).",
                palette=["white", "lightblue"],
                domain=[combined_summary_df["Top 10 enriched GO:BP's (n = 91)."].min(),
                        combined_summary_df["Top 10 enriched GO:BP's (n = 91)."].max()],
                na_color="lightgray",
                autocolor_text=True,
                reverse=False,
            )
)

gf_table.show()