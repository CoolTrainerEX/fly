"""Fly main module."""

import json
import sys
from io import TextIOWrapper
from pathlib import Path
from zipfile import ZipFile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def gc_content(record: SeqRecord) -> float:
    """Gets the GC Content of the record.

    Args:
        record (SeqRecord): Record to analyze

    Returns:
        float: GC Content (%)
    """
    return (
        (record.lower().count("g") + record.lower().count("c"))
        / len(record.seq or "")
        * 100
    )


def main():
    """Fly main function."""
    data_path = Path("ncbi_dataset/data")
    records: list[SeqRecord] = []

    with (
        ZipFile("ncbi_dataset.zip") as f,
        f.open((data_path / "dataset_catalog.json").as_posix()) as catalog,
    ):
        file_paths = [
            Path(file["filePath"])
            for assembly in json.load(catalog)["assemblies"]
            for file in assembly["files"]
            if file["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA"
        ]

        for file_path in file_paths:
            with TextIOWrapper(f.open((data_path / file_path).as_posix())) as file:
                records += SeqIO.parse(file, "fasta")

    while True:
        input_id = input("ID: ")

        try:
            print(
                "GC Content:",
                f"{gc_content(next(record for record in records if record.id == input_id)):.2f}",
            )
            break
        except StopIteration:
            print("ID not found", file=sys.stderr)


if __name__ == "__main__":
    main()
