# This script was used to extract high-confidence gene models from the maker-round 1 output.gff.

from pathlib import Path
from gffutils import create_db, FeatureDB


# Funcs.
def parse_gff_file(gff_file: Path) -> FeatureDB:

    if not gff_path.exists():
        raise FileNotFoundError(f"GFF file not found: {gff_path}")

    # Setup.
    db_dir = gff_path.parent
    db_dir = db_dir or gff_path.parent
    db_path = db_dir / f"{gff_path.stem}.db"

    if not db_path.exists():
        print(f"Creating database for {gff_path.name}...")
        create_db(
            data=str(gff_path),
            dbfn=str(db_path),
            force=True,
            keep_order=True,
            merge_strategy="create_unique",
            sort_attribute_values=True
        )
    else:
        print(f"Using existing database: {db_path.name}")

    return FeatureDB(str(db_path), keep_order=True)

def extract_high_conf_hits(db: FeatureDB, max_aed: float = 0.25, min_len: int = 50) -> list:
    high_conf_genes = []

    for gene in db.features_of_type("mRNA"):
        aed = float(gene.attributes.get('_AED', [1])[0])
        print(f"AED: {aed}")
        gene_len = gene.end - gene.start + 1

        if aed <= max_aed and gene_len >= min_len:
            high_conf_genes.append(gene)

        else:
            print("Gene Not considered High conf hit.")

    print(f"Found {len(high_conf_genes)} high-confidence genes (AED ≤ {max_aed}, length ≥ {min_len} bp)")
    return high_conf_genes

def write_high_conf_gff(db: FeatureDB, high_conf_genes: list, out_gff: Path) -> None:
    with open(out_gff, "w") as f:
        for gene in high_conf_genes:
            # Write the gene itself
            f.write(str(gene) + "\n")

            # Get all child features (mRNA, CDS, exon, UTR, etc.)
            for child in db.children(gene, level=None, order_by="start"):
                f.write(str(child) + "\n")

    print(f"Wrote {len(high_conf_genes)} high-confidence genes to {out_gff}")


if __name__=="__main__":
    gff_path = Path("../../data/maker_round_2/maker_round1_all.gff")
    db = parse_gff_file(gff_path)
    high_conf_genes = extract_high_conf_hits(db)
    write_high_conf_gff(db, high_conf_genes, Path("./high_conf_genes.gff"))


