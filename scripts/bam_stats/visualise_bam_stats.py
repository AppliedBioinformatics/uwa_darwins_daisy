from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from typing import List

# Funcs.
def parse_samtools_stats(filepath: str) -> dict:
    """
    Parse key metrics from a samtools stats file into a dictionary.
    Uses pathlib for file handling.
    """

    filepath = Path(filepath)

    fields_of_interest = {
        'raw total sequences': ('sequences', lambda x: int(x.split()[0])),
        'filtered sequences': ('filtered_sequences', lambda x: int(x.split()[0])),
        'reads mapped': ('reads_mapped', lambda x: int(x.split()[0])),
        'reads properly paired': ('reads_properly_paired', lambda x: int(x.split()[0])),
        'error rate': ('error_rate', float),
        'average length': ('average_length', float),
        'insert size average': ('insert_size_average', float),
        'insert size standard deviation': ('insert_size_stddev', float)
    }

    metrics = {'file': filepath.stem.replace('_sorted', '')}

    for _, (field, _) in fields_of_interest.items():
        metrics[field] = None

        with filepath.open() as f:
            for line in f:
                if line.startswith('SN'):
                    try:
                        _, rest = line.split('\t', 1)
                        key, value = rest.split(':', 1)
                        key = key.strip()
                        value = value.strip()
                        if key in fields_of_interest:
                            dict_key, caster = fields_of_interest[key]
                            metrics[dict_key] = caster(value)
                    except ValueError:
                        # If line is malformed, skip gracefully
                        continue

        return metrics

def build_bam_stats_dataframe(stats_folder_path:str) -> pd.DataFrame:
    """
    Parses all .stats files in a folder and returns a pandas dataframe with a number of metrics for each sample
    as a pandas dataframe. This dataframe can then be visualised using plotly express.
    """
    stats_folder = Path(stats_folder_path)
    files = list(stats_folder.glob('*.stats'))
    print(f"Found {len(files)} .stats files")
    results = [parse_samtools_stats(f) for f in files]
    return pd.DataFrame(results)

def build_pct_mapped_bar_plot(df: pd.DataFrame) -> px.bar:
    """Build a plot showing the percentage of mapped reads per sample in a pandas dataframe."""

    fig = px.bar(
        df,
        x='file',
        y='percent_reads_mapped',
        title='Percentage of Reads Mapped per Sample',
        labels={'file': 'Sample', 'percent_mapped': 'Mapped (%)'},
        template='plotly_dark'
    )

    fig.update_layout(yaxis_range=[90, 100])
    return fig

def build_pct_mapped_box_plot(df: pd.DataFrame) -> px.box:
    """Build a plot showing the distribution of mapped reads per sample in a pandas dataframe."""

    fig = px.box(df, y='percent_reads_mapped', title='Distribution of % Reads Mapped', template='plotly_dark')
    return fig

def build_pct_mapped_hist_plot(df: pd.DataFrame) -> px.histogram:
    """Build a plot showing the distribution of mapped reads per sample in a pandas dataframe."""
    fig = px.histogram(
        df,
        x='percent_reads_mapped',
        nbins=20,
        title='Distribution of % Reads Mapped',
        labels={'percent_reads_mapped': 'Mapped (%)', "count": "Number of samples."},
        template='plotly_dark'
    )

    return fig

def build_total_mapped_reads_plot(df: pd.DataFrame) -> px.bar:
    fig = go.Figure(data=[
        go.Bar(name="Total Number of Reads", x=df["file"], y=df["sequences"]),
        go.Bar(name="Reads Mapped", x=df["file"], y=df["reads_mapped"])],
    )

    fig.update_layout(title="Total Reads/Total Reads mapped.",
                      xaxis_title="Sample",
                      yaxis_title="Read Count",
                      barmode="group",
                      template = 'plotly_dark'
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
        <h1>Darwin's Daisy bowtie2 alignment summary statistics:  </h1>
        {"".join(f'<div class="plot">{chunk}</div>' for chunk in html_chunks)}
    </body>
    </html>
    """

    # Save.
    outfile_path = Path(outfile)
    outfile_path.write_text(full_html)
    print(f"Report written to: {outfile_path.resolve()}")

if __name__ == "__main__":

    df = build_bam_stats_dataframe(stats_folder_path="../../data/bam_stats/")
    df["percent_reads_mapped"] = df["reads_mapped"] / df["sequences"] * 100

    # Build figs.
    figs = [
        build_pct_mapped_bar_plot(df),
        build_pct_mapped_box_plot(df),
        build_pct_mapped_hist_plot(df),
        build_total_mapped_reads_plot(df),
    ]

    build_report(figs, "../../reports/alignment_summary_report.html")






