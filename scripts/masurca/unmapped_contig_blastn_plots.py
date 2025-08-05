# This script is used to generate some summary plots to assess levels of contamination of the unmapped contigs
# blastn query.

import pandas as pd
from tabulate import tabulate
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Filter funcs.
def save_filtered_ids(blast_df, outfile_path: Path, column="phylum", value="Streptophyta") -> None:
    filtered_df = blast_df[blast_df[column].str.lower() == value.lower()]
    filtered_df['squeryid'].to_csv(outfile_path, index=False, header=False)
    print(f"Saved {filtered_df['squeryid']} unique squeryid's where {column} == '{value}' to {outfile_path}")


# Plot funcs.
def plt_superkingdom_distribution(blast_df):
    # Collapse superkingdom categories
    blast_df = blast_df.copy()  # avoid modifying original
    blast_df['superkingdom'] = blast_df['superkingdom'].apply(
        lambda x: x if x in ['Eukaryota', 'Bacteria'] else 'Other'
    )

    # Create plot
    fig, ax = plt.subplots(figsize=(6, 6))
    blast_df['superkingdom'].value_counts().plot(
        kind='pie', autopct='%1.1f%%', colors=sns.color_palette('Set2'), ax=ax
    )
    ax.set_title("Distribution of Hits (Bacteria, Eukaryota, Other)")
    ax.set_ylabel("")  # remove ylabel
    return fig

def plt_top_phyla(blast_df):
    top_phyla = blast_df['phylum'].value_counts().nlargest(10)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(x=top_phyla.values, y=top_phyla.index, palette='viridis', ax=ax)
    ax.set_title("Top 10 Phyla by Hit Count")
    ax.set_xscale('log')
    ax.set_xlabel("Number of Hits")
    ax.set_ylabel("Phylum")
    return fig


def plt_top_species_streptophyta(blast_df, n=10):
    """
    Filters Streptophyta hits, groups by species, pools outside top n as 'Other', and plots counts.
    """
    # Filter for Streptophyta
    strepto_df = blast_df[blast_df['phylum'].str.lower() == 'streptophyta']

    # Count species
    species_counts = strepto_df['species'].value_counts()

    # Pool species outside the top n
    if len(species_counts) > n:
        top_species = species_counts.nlargest(n)
        other_count = species_counts.iloc[n:].sum()
        species_counts = pd.concat([top_species, pd.Series({'Other': other_count})])

    # Plot
    fig, ax = plt.subplots(figsize=(16, 6))
    sns.barplot(x=species_counts.values, y=species_counts.index, palette='viridis', ax=ax)
    ax.set_title(f"Top {n} Streptophyta Species (Others Pooled)")
    ax.set_xlabel("Number of Hits")
    ax.set_xscale('log')
    ax.set_ylabel("Species")
    return fig



if __name__=="__main__":
    blast_df = pd.read_csv("../../data/masurca/blastn_with_taxonomy.csv")
    print(tabulate(blast_df.head(), headers="keys"))
    plant_fasta_path = "../../data/masurca/plant_fasta_ids.txt"
    save_filtered_ids(blast_df, plant_fasta_path)
    plt_superkingdom_distribution(blast_df).savefig("../../plots/masurca/blastn_superkingdom_distribution.png")
    plt_top_phyla(blast_df).savefig("../../plots/masurca/blastn_top_phyla.png")
    plt_top_species_streptophyta(blast_df, n=20).savefig("../../plots/masurca/blastn_top_species_streptophyta.png")



