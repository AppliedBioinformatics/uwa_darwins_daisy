#!/usr/bin/env python3
import argparse
import logging
from typing import NoReturn
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def clean_fasta(input_file: str, output_file: str) -> None:
    """
    Cleans a protein FASTA file by removing invalid characters (e.g., '.') from sequences.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to write the cleaned FASTA file.
    """
    try:
        with open(output_file, "w") as out_f:
            for record in SeqIO.parse(input_file, "fasta"):
                clean_seq = str(record.seq).replace(".", "")
                clean_record = SeqRecord(Seq(clean_seq), id=record.id, description=record.description)
                SeqIO.write(clean_record, out_f, "fasta")
        logging.info(f"Cleaned file written to: {output_file}")
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}")
        raise
    except Exception as e:
        logging.error(f"An error occurred while cleaning the FASTA file: {e}")
        raise


def main() -> NoReturn:
    """
    Parses command-line arguments and runs the FASTA cleaning process.
    """
    parser = argparse.ArgumentParser(
        description="Clean protein FASTA file by removing invalid characters (e.g., '.') from sequences."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output cleaned FASTA file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (debug) logging output"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    try:
        clean_fasta(args.input, args.output)
    except Exception:
        logging.error("Cleaning failed.")
        exit(1)
    exit(0)


if __name__ == "__main__":
    main()
