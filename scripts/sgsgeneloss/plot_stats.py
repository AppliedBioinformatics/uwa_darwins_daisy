from pathlib import Path
from pickle import FALSE
from typing import List, Optional
from sklearn.decomposition import PCA
import pandas as pd
import plotly.io as pio
import plotly.graph_objs as go
import plotly.express as px
import random
from tabulate import tabulate

# Format.
def find_all_files(folder_path: Path, file_type: str) -> List[Path]:
    """Searches a directory for all files with the specified filetype and returns a list of Path objects for
     those files."""

    root = folder_path
    files = list(root.rglob(f"*{file_type}"))
    return files

def parse_stat_file(file:Path) -> dict:
    """Parses the stats.txt output from SGSGeneLoss into a Python dictionary for a single <sample>stats.txt file."""
    result = {}

    # Parse bam file name.
    lines = file.read_text().splitlines()
    bam_name = lines[1].split(":", 1)[1].strip()
    result["sample_id"] = bam_name.split("_sorted.bam")[0]

    # Parse other metrics.
    for line in lines[8:18]:
        parts = line.split(": ", 1)
        key = str(parts[0].strip()).replace(" ", "_")
        value = parts[1].strip()
        result[key] = round(float(value), 3)

    assert len(result) == 11
    return result

def build_stats_df(stat_file_paths: List[Path]) -> pd.DataFrame:
    """Generates a pandas dataframe containing sample.stats for SGSGeneLoss outputs."""
    records  = [parse_stat_file(f) for f in stat_file_paths]
    return pd.DataFrame(records)

def extract_present_subsample(pav_df: pd.DataFrame, n:int = None) -> List[str]:
    """Identifies a list of "core" genes using the PAV matrix, will randomly subset the data based on user input."""
    df = pav_df.copy()
    core_gene_mask = df.eq(1).all(axis=1)
    core_genes = df.index[core_gene_mask].tolist()
    print(len(core_genes))

    if n is not None:
        return random.sample(core_genes, n)

    return core_genes

def build_gene_length_table(excov_file:Path) -> pd.DataFrame:
    """Generates a pandas dataframe with rows as genes, and columns for start, end and length in BP."""
    df = pd.read_csv(excov_file, usecols=["ID", "start_position", "end_postion"])
    df = df.set_index("ID")
    df["length"] = df["end_postion"] - df["start_position"]
    return df


# Plots.
def plt_pct_total_genes_lost_hist(raw_df: pd.DataFrame, nbins: Optional[int] = 50) -> go.Figure:
    df = raw_df.copy()
    df["pct_genes_lost"] = (df["Total_number_of_genes_lost"] / df["Total_number_of_genes"]) * 100

    fig = px.histogram(df,
                       x="pct_genes_lost",
                       nbins=nbins,
                       title=f"Percentage of Genes Lost ({df['Total_number_of_genes'].max()}).",
                       labels={"percent_genes_lost": "Percentage of Genes Lost"},
                       )

    fig.update_layout(
        xaxis_title="% Genes Lost",
        yaxis_title=f"Number of Samples ({len(df)}).",
        bargap=0.1
    )

    return fig

def plt_total_genes_lost_hist(raw_df: pd.DataFrame, nbins: Optional[int] = 100) -> go.Figure:
    df = raw_df.copy()
    fig = px.histogram(df,
                       x="Total_number_of_genes_lost",
                       nbins=nbins,
                       title=f"Total number of Genes Lost ({df['Total_number_of_genes'].max()}).",
                       labels={"Total_number_of_genes_lost": "Total number of Genes Lost"},
                       )

    fig.update_layout(
        xaxis_title="Total Genes Lost",
        yaxis_title=f"Number of Samples ({len(df)}).",
        bargap=0.1
    )

    return fig

def plt_genes_lost_avg_gene_length_sct(raw_df: pd.DataFrame ) -> go.Figure:
    df = raw_df.copy()

    fig = px.scatter(
        df,
        x="Total_number_of_genes_lost",
        y="Average_gene_length,_genes_lost",
        hover_name="sample_id",
        title="Lost Genes: Total vs Average Length.",
        labels={"Total_number_of_genes_lost": "Total number of Genes Lost",
                "Average_gene_length,_genes_lost": "Average Length of lost genes."},
    )

    return fig

def plt_lost_gene_sizes_box(raw_df: pd.DataFrame) -> go.Figure:
    """Generates a box plot visualising lost gene sizes."""
    df = raw_df.copy()
    fig = px.box(
        df,
        y="Average_gene_length,_genes_lost",
        points="all",
        hover_name="sample_id",
        title="Distribution of average length of lost genes per Sample.",
        labels={"Average_gene_length,_genes_lost": "Average Length of lost genes."},
    )

    return fig

def plt_lost_vs_present_avg_gene_length_box(raw_df: pd.DataFrame) -> go.Figure:
    col_lost = 'Average_gene_length,_genes_lost'
    col_not_lost = 'Average_gene_length,_genes_not_lost'

    df = raw_df.copy()
    df_long = df.melt(id_vars="sample_id",
                      value_vars=[col_lost, col_not_lost],
                      var_name="Gene_Status",
                      value_name="average_gene_length")

    df_long["Gene_Status"] = df_long["Gene_Status"].map({
        col_lost: "Lost Genes",
        col_not_lost: "Retained Genes",
    })

    fig = px.box(
        df_long,
        x='Gene_Status',
        y='average_gene_length',
        points='all',
        color='Gene_Status',
        hover_name='sample_id',
        title='Comparison of Average Gene Lengths: Lost vs Retained'
    )

    fig.update_layout(
        yaxis_title='Average Gene Length',
        xaxis_title='Gene Status',
        boxmode='group'
    )

    return fig

def plt_total_read_count_v_number_present_genes(raw_df: pd.DataFrame) -> go.Figure:
    bam_stats_df = pd.read_csv("../../data/bam_stats/bam_stats_df.csv")
    bam_stats_df = bam_stats_df.rename(columns={'file': 'sample_id'})
    print(tabulate(bam_stats_df, headers='keys', tablefmt='psql'))
    merged_df = pd.merge(raw_df, bam_stats_df, on='sample_id')
    merged_df["number_present_genes"] = merged_df["Total_number_of_genes"] - merged_df["Total_number_of_genes_lost"]

    fig = px.scatter(
        merged_df,
        x="reads_mapped",
        y="number_present_genes",
        log_y=True,
        hover_name="sample_id",
        title="Number of genes genes marked present V Total read count.",
        labels={"reads_mapped": "Total Number of reads mapped.",
                "number_present_genes": "Number of genes marked present.",}
    )

    return fig

def plt_pav_matrix_pca(pav_df: pd.DataFrame) -> go.Figure:
    X = pav_df.T
    X = X.loc[:, (X != X.iloc[0]).any()]
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)
    pca_df = pd.DataFrame({
        'Sample': X.index,
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1]
    })
    fig = px.scatter(
        pca_df,
        x='PC1',
        y='PC2',
        title='PCA of PAV Matrix',
        hover_name="Sample",
        labels={
            'PC1': f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%)",
            'PC2': f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%)"
        }
    )
    fig.update_traces(textposition='top center')
    return fig

def plt_presence_v_coverage(pav_df: pd.DataFrame, cov_df: pd.DataFrame) -> go.Figure:
    assert (pav_df.index == cov_df.index).all()

    # Create summary dataframe.
    summary_df = pd.DataFrame({
        "sample": pav_df.columns,
        "num_present_genes": pav_df.sum(),
        "average_read_depth": cov_df.mean()
    })

    fig = px.scatter(
        summary_df,
        x="num_present_genes",
        y="average_read_depth",
        hover_name="sample",
        title="Genes Present vs. Average read depth",
        labels={"average_read_depth": "Average Read Depth",
                "num_present_genes": "Number of Present Genes (Log10)"
                },
    )

    return fig.update_xaxes(type="log")

def plt_core_gene_length_box(len_df: pd.DataFrame, core_gene_list: List[str]) -> go.Figure:
    """Creates a box plot to visualise the gene sizes of genes present in all samples (core genes)."""
    df = len_df.copy()
    df['is_core'] = df.index.isin(core_gene_list)

    fig = px.box(
        df,
        x="is_core",
        y="length",
        log_y=True,
        color="is_core",
        title="Core vs Non-core Gene Length.",
        labels={"is_core": "Feature set", "length": "Feature Length (BP)"},
        points="all",

    )

    fig.update_traces(boxmean=True)
    fig.update_xaxes(
        tickvals=[True, False],
        ticktext=[
            f"Core Genes (n={df['is_core'].sum()})",
            f"Non-Core Genes (n={(~df['is_core']).sum()})"
        ]
    )

    return fig

def build_report(fig_list: List[go.Figure], outfile: str) -> None:
    """Combines all plotly figures passed as the argument into a single html report. Keeps interactivity for each
    plot, so users can scroll between charts rather than having to open new pages."""


    # Convert to html.
    html_chunks = [pio.to_html(fig, full_html=False, include_plotlyjs=False) for fig in fig_list]
    full_html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{
                font-family: sans-serif;
                margin: 2em;
            }}
            .plot {{
                margin-bottom: 60px;
                page-break-after: always;
            }}
        </style>
    </head>
    <body>
        <h1>02/07/2025 Darwin Daisies SGSGeneLoss (Default Params) :  </h1>
        {"".join(f'<div class="plot">{chunk}</div>' for chunk in html_chunks)}
    </body>
    </html>
    """

    # Save.
    outfile_path = Path(outfile)
    outfile_path.write_text(full_html)
    print(f"Report written to: {outfile_path.resolve()}")


if __name__ == "__main__":

    # Settings.
    pio.templates.default = "plotly_dark"

    # Build .stats dataframe.
    stats_files =find_all_files(Path("../../sgsgeneloss_results/"), file_type=".txt")
    raw_df = build_stats_df(stats_files)
    raw_df["sample_id"] = raw_df["sample_id"].str.removesuffix("_merged.bam")

    # Build pav dataframe and remove outlier samples.
    pav_df = pd.read_csv("../../data/sgsgeneloss/pav_matrix.csv", index_col=0)

    # Extract core genes list (present in 100% of individuals).
    core_gene_list = extract_present_subsample(pav_df=pav_df)

    # Build gene length dataframe.
    len_df = build_gene_length_table(excov_file=Path("../../sgsgeneloss_results/HALM12_19_merged_all.excov"))

    # Build coverage dataframe.
    cov_df = pd.read_csv("../../data/sgsgeneloss/cov_matrix.csv", index_col=0)

    # Call functions and generate report.
    figs = [plt_total_genes_lost_hist(raw_df),
            plt_pct_total_genes_lost_hist(raw_df),
            plt_total_read_count_v_number_present_genes(raw_df),
            plt_genes_lost_avg_gene_length_sct(raw_df),
            plt_lost_vs_present_avg_gene_length_box(raw_df),
            plt_lost_gene_sizes_box(raw_df),
            plt_pav_matrix_pca(pav_df),
            plt_presence_v_coverage(pav_df, cov_df),
            plt_core_gene_length_box(len_df, core_gene_list)]

    build_report(figs, "../../reports/sgsgeneloss_report.html")