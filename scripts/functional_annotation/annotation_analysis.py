import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

# GLOBALS
TOTAL_GFF_FEATURES = 43093 # (Number of features described in maker gff file).
SIG_THRESHOLD = 1e-5
DIAMOND_PROT_RESULTS_FILE = Path("../../data/functional_annotation/diamond_results_protein.tsv")
OUTDIR_PATH = Path("../../plots/functional_annotation/")
DIAMOND_TSV_HEADERS = ["query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens", "q_start",
                       "q_end", "s_start", "s_end", "evalue", "bit_score", "staxids", "stitle", "sscinames",
                       "sskingdoms", "skingdoms", "sphylums"]

# Function for saving figures.
def save_plotly_figure(fig: go.Figure, name:str, extension: str = "png") -> None:
    file_path = OUTDIR_PATH / f"{name}"
    fig.write_image(f"{file_path}.{extension}", width=1200, height=1000, scale=3)
    print(f"Figure saved to {file_path}.{extension}")

# Functions for plots.
def plt_identity_scores_hist(ddf = pd.DataFrame) -> px.histogram:
    df = ddf.copy()
    fig = px.histogram(df, x="identity", nbins=50,
                       title="Distribution of Identity Scores (DIAMOND Hits)",
                       labels={"identity": "Percent Identity"},
                       template="plotly_white")
    fig.update_layout(bargap=0.1)
    return df

def plt_bit_score_hist(ddf = pd.DataFrame) -> px.histogram:
    df = ddf.copy()
    plt = px.histogram(df,
                       x="bit_score",
                       nbins=100,
                       title="Distribution of Bit Score (DIAMOND Hits)",
                       labels={"bit_score": "Bit Score"},
                       template="plotly_white"
                   )
    plt.update_layout(bargap=0.1)
    return plt

def plt_alignment_length_hist(ddf = pd.DataFrame) -> px.histogram:
    df = ddf.copy()
    plt = px.histogram(df,
                       x="alignment_length",
                       nbins=100,
                       title="Distribution of Alignment Length (DIAMOND Hits)",
                       labels={"alignment_length": "Alignment Length"},
                       template="plotly_white"
                        )
    plt.update_layout(bargap=0.1)
    return plt

def plt_identity_v_bitscore_sct(ddf = pd.DataFrame) -> px.scatter:
    df = ddf.copy()
    plt = px.scatter(df,
                     x="bit_score",
                     log_x=True,
                     y="identity",
                     title="Diamond Hits Bit-score by Identity.",
                     labels={"bit_score": "Bit Score", "identity": "Percent Identity"},
                     opacity=0.6,
                     template="plotly_white")
    plt.update_traces(marker=dict(size=4))

def plt_alignment_length_v_identity_sct(ddf = pd.DataFrame) -> px.scatter:
    df = ddf.copy()
    plt = px.scatter(ddf,
                     x="alignment_length",
                     y="identity",
                     title="Diamond Hits Alignment Length by Identity.",
                     labels={"alignment_length": "Alignment Length", "identity": "Percent Identity"},
                     opacity=0.6,
                     template="plotly_white"
                     )
    plt.update_traces(marker=dict(size=4))
    return plt

def plt_counts_per_kingdom(ddf = pd.DataFrame) -> px.bar:
    df = ddf.copy()

    #Reformat and Tidy.
    df["skingdoms"] = df["skingdoms"].replace("0", "Uncategorised")
    kingdom_counts = df["skingdoms"].value_counts().reset_index()
    kingdom_counts.columns = ["Kingdom", "Count"]

    plt = px.bar(
        kingdom_counts,
        x="Kingdom",
        y="Count",
        log_y=True,
        title="Frequency of Top BLAST Hits by Kingdom.",
        text="Count",
        color="Kingdom",
        color_discrete_sequence=px.colors.qualitative.Vivid
    )
    plt.update_traces(textposition="outside")
    plt.update_layout(
        xaxis_title="Taxonomic Kingdom",
        yaxis_title="Number of Hits",
        uniformtext_minsize=8,
        uniformtext_mode='hide'
    )
    return plt

def plt_counts_per_species(ddf = pd.DataFrame, n: int = 10) -> px.bar:
    df = ddf.copy()

    # Format, bin and tidy.
    df["sscinames"] = df["sscinames"].replace("0", "Uncategorised").fillna("Uncategorised")
    species_counts = df["sscinames"].value_counts()
    top_species = species_counts.nlargest(n - 1)
    other_count = species_counts.iloc[n - 1:].sum()
    plot_df = pd.concat([top_species, pd.Series({"Other": other_count})]).reset_index()
    plot_df.columns = ["Species", "Count"]

    plt = px.bar(
        plot_df,
        x="Species",
        y="Count",
        log_y=True,
        title="Top 9 Species by DIAMOND BLASTP Hits (Others Grouped)",
        text="Count",
        color="Species",
        color_discrete_sequence=px.colors.qualitative.Set3
    )

    plt.update_traces(textposition="outside")
    plt.update_layout(
        xaxis_title="Species",
        yaxis_title="Number of Hits",
        uniformtext_minsize=8,
        uniformtext_mode='hide'
    )

    return plt

if __name__=="__main__":
    # Read the file into a pandas dataframe.
    ddf = pd.read_csv(DIAMOND_PROT_RESULTS_FILE, sep="\t", names=DIAMOND_TSV_HEADERS)

    # Summary metrics.
    n_total_hits = len(ddf)
    n_unique_queries = ddf["query_id"].nunique()

    # Significant hits based on revalue.
    sig_hits = ddf[ddf["evalue"] <= SIG_THRESHOLD]
    n_sig = len(sig_hits)

    # Build plots.
    fig = plt_counts_per_species(ddf=ddf, n=15)
    save_plotly_figure(fig, "diamond_hits_counts_per_species")
