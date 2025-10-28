from typing import List
from pathlib import Path
import pandas as pd
from functools import reduce
from tabulate import tabulate

def parse_folder_for_samples(folder_path: Path) -> List[Path]:
    """Parses a parent folder and returns a list of Path objects for each child folder found within it."""
    return [p for p in folder_path.iterdir() if p.is_dir()]

def merge_excovs_in_folder(sample_folder_path: Path, outdir_path: Path) -> None:
    """Parses the folder path supplied as the argument for .excov files and merges them into a single output file."""
    sample_name = sample_folder_path.name
    out_file = outdir_path / f"{sample_name}_all.excov"

    header_written = False
    with out_file.open('w') as out:
        for excov_file in sorted(sample_folder_path.glob("*.excov")):
            with excov_file.open() as f:
                for i, line in enumerate(f):
                    if i == 0 and header_written:
                        continue  # skip header if already written
                    out.write(line)
            header_written = True

    print(f"Merged: {sample_name} -> {out_file}")

def merge_excovs_for_sample_set(master_folder: str, outdir_folder: str) -> None:
    """Will create a merged .excov file for all samples found within the highest level folder."""

    # Turn strings into paths.
    master_folder_path = Path(master_folder)
    outdir_folder_path = Path(outdir_folder)
    outdir_folder_path.mkdir(parents=True, exist_ok=True)


    # .excov files for each sample and merge, send results to outdir.
    sample_file_paths = parse_folder_for_samples(master_folder_path)
    for sample in sample_file_paths:
        merge_excovs_in_folder(sample, outdir_folder_path)

    print(f"Created merged .excov file for {len(sample_file_paths)} samples within {master_folder_path.name}.")

def create_pav_matrix(infolder:str) -> pd.DataFrame:
    """Build a PAV matrix from sample_merged.excov files"""

    # Make str into Path objects.
    infolder_path = Path(infolder)

    # Get file list to include in matrix.
    if not infolder_path:
        raise FileNotFoundError(f"Input folder does not exist: {infolder}")

    excov_files = sorted(infolder_path.glob("*.excov"))
    if not excov_files:
        raise FileNotFoundError(f"No .excov files found in {infolder}")

    # Only interested in geneID and is_lost columns.
    pav_dfs = []
    for file_path in excov_files:
        sample_name = file_path.stem
        df = pd.read_csv(file_path, engine="pyarrow")

        df = df[["ID", "is_lost"]].copy()
        df["is_lost"] = df["is_lost"].map({"PRESENT": 1, "LOST": 0})
        df = df.rename(columns={"is_lost": sample_name})
        pav_dfs.append(df)

    # Merge on 'ID' (gene IDs)
    pav_matrix = reduce(lambda left, right: pd.merge(left, right, on="ID", how="outer"), pav_dfs)
    pav_matrix = pav_matrix.fillna(0).set_index("ID").astype(int)

    return pav_matrix

def create_coverage_matrix(infolder:str) -> pd.DataFrame:
    """Builds a coverage matrix for all samples in the input folder."""
    infolder_path = Path(infolder)
    excov_files = sorted(infolder_path.glob("*.excov"))

    coverage_dfs = []
    for file_path in excov_files:
        sample_name = file_path.stem
        df = pd.read_csv(file_path, engine="pyarrow")

        # Select geneID and coverage.
        df = df[["ID", "ave_cove_depth_gene"]].copy()
        df = df.rename(columns={"ave_cove_depth_gene": sample_name})
        coverage_dfs.append(df)

    coverage_matrix = reduce(lambda left, right: pd.merge(left, right, on="ID", how="outer"), coverage_dfs)
    coverage_matrix = coverage_matrix.fillna(0).set_index("ID").astype(float)

    return coverage_matrix

if __name__ == "__main__":
    merge_excovs_for_sample_set(master_folder=r"../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_yacon", outdir_folder="../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_yacon_merged")
    pav_df = create_pav_matrix(infolder="../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_yacon_merged")
    pav_df.to_csv("../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_pav_matrix.csv", index=True)
    # Create coverage matrix.
    cov_df = create_coverage_matrix(infolder="../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_yacon_merged")
    cov_df.to_csv("../../data/sgsgeneloss_/250408_mainland_species/run_2/mainland_cov_matrix.csv", index=True)


