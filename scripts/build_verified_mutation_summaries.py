"""Build processed summaries from verified mutation evidence only."""

from collections import defaultdict
import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
OUTPUT_DIR = ROOT / "data" / "processed"

VERIFIED_EVIDENCE_PATH = OUTPUT_DIR / "verified_mutation_evidence.tsv"
SPECIES_SUMMARY_PATH = OUTPUT_DIR / "mutation_summary_by_species.tsv"
HERBICIDE_SUMMARY_PATH = OUTPUT_DIR / "mutation_summary_by_herbicide.tsv"


def read_rows(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing source table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def write_rows(path, fieldnames, rows):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def split_terms(value):
    return [term.strip() for term in value.split(";") if term.strip()]


def sorted_join(values):
    return "; ".join(sorted(value for value in values if value))


def build_species_summary(rows):
    grouped = defaultdict(lambda: {"mutations": set(), "papers": set(), "row_count": 0})
    for row in rows:
        species = row.get("species", "").strip()
        grouped[species]["row_count"] += 1
        grouped[species]["mutations"].add(row.get("mutation", "").strip())
        grouped[species]["papers"].add(row.get("paper_id", "").strip())

    summary_rows = []
    for species in sorted(grouped):
        group = grouped[species]
        summary_rows.append(
            {
                "species": species,
                "verified_evidence_rows": str(group["row_count"]),
                "unique_mutations": str(len(group["mutations"])),
                "mutations": sorted_join(group["mutations"]),
                "paper_ids": sorted_join(group["papers"]),
            }
        )
    return summary_rows


def build_herbicide_summary(rows):
    grouped = defaultdict(
        lambda: {"mutations": set(), "species": set(), "papers": set(), "row_count": 0}
    )
    for row in rows:
        herbicides = split_terms(row.get("specific_herbicides_tested", ""))
        for herbicide in herbicides:
            grouped[herbicide]["row_count"] += 1
            grouped[herbicide]["mutations"].add(row.get("mutation", "").strip())
            grouped[herbicide]["species"].add(row.get("species", "").strip())
            grouped[herbicide]["papers"].add(row.get("paper_id", "").strip())

    summary_rows = []
    for herbicide in sorted(grouped):
        group = grouped[herbicide]
        summary_rows.append(
            {
                "herbicide": herbicide,
                "verified_evidence_rows": str(group["row_count"]),
                "unique_mutations": str(len(group["mutations"])),
                "mutations": sorted_join(group["mutations"]),
                "species": sorted_join(group["species"]),
                "paper_ids": sorted_join(group["papers"]),
            }
        )
    return summary_rows


def main():
    fieldnames, rows = read_rows(SOURCE_PATH)
    verified_rows = [
        row for row in rows if row.get("verification_status", "").strip() == "verified"
    ]

    write_rows(VERIFIED_EVIDENCE_PATH, fieldnames, verified_rows)
    write_rows(
        SPECIES_SUMMARY_PATH,
        [
            "species",
            "verified_evidence_rows",
            "unique_mutations",
            "mutations",
            "paper_ids",
        ],
        build_species_summary(verified_rows),
    )
    write_rows(
        HERBICIDE_SUMMARY_PATH,
        [
            "herbicide",
            "verified_evidence_rows",
            "unique_mutations",
            "mutations",
            "species",
            "paper_ids",
        ],
        build_herbicide_summary(verified_rows),
    )

    print(f"verified rows written: {len(verified_rows)}")
    print(f"generated: {VERIFIED_EVIDENCE_PATH}")
    print(f"generated: {SPECIES_SUMMARY_PATH}")
    print(f"generated: {HERBICIDE_SUMMARY_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
