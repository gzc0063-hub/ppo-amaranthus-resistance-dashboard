"""Validate the curated mutation-level evidence table.

The table is manually curated; this script checks schema, required metadata,
controlled values, duplicate keys, and conservative warning conditions.
"""

from collections import Counter
import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"

REQUIRED_COLUMNS = [
    "mutation_id",
    "species",
    "gene_or_isoform",
    "protein_accession",
    "nucleotide_accession",
    "mutation",
    "wildtype_residue",
    "mutant_residue",
    "residue_position",
    "mutation_type",
    "herbicide_class",
    "specific_herbicides_tested",
    "resistance_phenotype",
    "evidence_type",
    "citation_key",
    "paper_id",
    "doi",
    "pmid",
    "figure_or_table",
    "verification_status",
    "notes",
]

ALLOWED_VERIFICATION_STATUSES = {"verified", "needs_manual_check", "rejected"}

VERIFIED_REQUIRED_FIELDS = [
    "species",
    "gene_or_isoform",
    "mutation",
    "mutation_type",
    "herbicide_class",
    "specific_herbicides_tested",
    "resistance_phenotype",
    "evidence_type",
    "notes",
]


def read_table(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing mutation evidence table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def split_terms(value):
    return [term.strip() for term in value.split(";") if term.strip()]


def main():
    fieldnames, rows = read_table(TABLE_PATH)
    errors = []
    warnings = []

    if fieldnames != REQUIRED_COLUMNS:
        errors.append("required columns do not match mutation_evidence_table.tsv schema")

    status_counts = Counter()
    duplicate_keys = Counter()
    mutations = set()
    species_values = set()
    herbicides = set()

    for row_number, row in enumerate(rows, start=2):
        status = row.get("verification_status", "").strip()
        status_counts[status] += 1

        if status not in ALLOWED_VERIFICATION_STATUSES:
            errors.append(
                f"row {row_number} has invalid verification_status: {status}"
            )

        mutation = row.get("mutation", "").strip()
        species = row.get("species", "").strip()
        herbicide_text = row.get("specific_herbicides_tested", "").strip()
        if mutation:
            mutations.add(mutation)
        if species:
            species_values.add(species)
        for herbicide in split_terms(herbicide_text):
            herbicides.add(herbicide)

        key = (
            row.get("paper_id", "").strip(),
            species,
            row.get("gene_or_isoform", "").strip(),
            mutation,
            herbicide_text,
        )
        duplicate_keys[key] += 1

        if status == "verified":
            missing = [
                field for field in VERIFIED_REQUIRED_FIELDS if not row.get(field, "").strip()
            ]
            if not row.get("citation_key", "").strip() and not row.get("paper_id", "").strip():
                missing.append("citation_key or paper_id")
            if not row.get("doi", "").strip() and not row.get("pmid", "").strip():
                missing.append("doi or pmid")
            if missing:
                errors.append(
                    f"row {row_number} verified row missing required fields: "
                    + ", ".join(missing)
                )
            if ";" in species:
                warnings.append(f"row {row_number} verified species contains semicolon")
            if ";" in row.get("gene_or_isoform", ""):
                warnings.append(
                    f"row {row_number} verified gene_or_isoform contains semicolon"
                )

        if status == "rejected" and not row.get("notes", "").strip():
            errors.append(f"row {row_number} rejected row is missing notes")

    duplicate_count = sum(1 for count in duplicate_keys.values() if count > 1)
    for key, count in duplicate_keys.items():
        if count > 1:
            warnings.append(
                "duplicate evidence key appears "
                f"{count} times: paper_id={key[0]}, species={key[1]}, "
                f"gene_or_isoform={key[2]}, mutation={key[3]}, "
                f"specific_herbicides_tested={key[4]}"
            )

    print("mutation evidence validation summary")
    print(f"total rows: {len(rows)}")
    print(f"verified rows: {status_counts['verified']}")
    print(f"needs_manual_check rows: {status_counts['needs_manual_check']}")
    print(f"rejected rows: {status_counts['rejected']}")
    print(f"unique mutations: {len(mutations)}")
    print(f"unique species: {len(species_values)}")
    print(f"unique herbicides: {len(herbicides)}")
    print(f"duplicate count: {duplicate_count}")
    print(f"warning count: {len(warnings)}")

    if warnings:
        print("warnings:")
        for warning in warnings:
            print(f"  - {warning}")

    if errors:
        print("validation errors:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
