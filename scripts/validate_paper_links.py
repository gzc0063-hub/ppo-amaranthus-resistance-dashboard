"""Validate curated paper-link metadata without downloading paper content."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PAPER_LINKS_PATH = ROOT / "references" / "paper_links.tsv"

REQUIRED_COLUMNS = [
    "paper_id",
    "paper_title",
    "authors",
    "year",
    "journal",
    "doi",
    "pmid",
    "url",
    "source_database",
    "species_mentioned",
    "ppo_gene_or_isoform",
    "mutation_terms_mentioned",
    "herbicide_terms_mentioned",
    "evidence_type",
    "verification_status",
    "notes",
]

ALLOWED_VERIFICATION_STATUS = {"verified", "needs_manual_check", "rejected"}


def clean(value: str | None) -> str:
    return (value or "").strip()


def row_has_content(row: dict[str, str | None]) -> bool:
    return any(clean(value) for value in row.values())


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str | None]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [row for row in reader if row_has_content(row)]


def validate_rows(fieldnames: list[str], rows: list[dict[str, str | None]]) -> tuple[list[str], list[str], Counter]:
    errors: list[str] = []
    warnings: list[str] = []
    counts: Counter = Counter()

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing_columns:
        errors.append(f"Missing required columns: {', '.join(missing_columns)}")
        return errors, warnings, counts

    seen_paper_ids: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        paper_id = clean(row.get("paper_id"))
        url = clean(row.get("url"))
        status = clean(row.get("verification_status"))
        doi = clean(row.get("doi"))
        pmid = clean(row.get("pmid"))
        species = clean(row.get("species_mentioned"))
        mutation_terms = clean(row.get("mutation_terms_mentioned"))
        evidence_type = clean(row.get("evidence_type"))

        if not paper_id:
            errors.append(f"Row {row_number}: paper_id is required.")
        elif paper_id in seen_paper_ids:
            errors.append(f"Row {row_number}: duplicate paper_id '{paper_id}'.")
        else:
            seen_paper_ids.add(paper_id)

        if not url:
            errors.append(f"Row {row_number}: url is required.")

        if status not in ALLOWED_VERIFICATION_STATUS:
            errors.append(
                "Row "
                f"{row_number}: verification_status must be one of "
                f"{', '.join(sorted(ALLOWED_VERIFICATION_STATUS))}."
            )
        else:
            counts[status] += 1

        if not doi and not pmid:
            counts["missing_doi_and_pmid"] += 1
            warnings.append(f"Row {row_number}: doi and pmid are both empty.")

        if not species:
            counts["missing_species"] += 1
            warnings.append(f"Row {row_number}: species_mentioned is empty.")

        if not mutation_terms:
            counts["missing_mutation_terms"] += 1
            warnings.append(f"Row {row_number}: mutation_terms_mentioned is empty.")

        if not evidence_type:
            warnings.append(f"Row {row_number}: evidence_type is empty.")

    counts["total_rows"] = len(rows)
    return errors, warnings, counts


def print_summary(errors: list[str], warnings: list[str], counts: Counter) -> None:
    print("Paper link validation summary")
    print("=============================")
    print(f"total rows: {counts['total_rows']}")
    print(f"number verified: {counts['verified']}")
    print(f"number needs_manual_check: {counts['needs_manual_check']}")
    print(f"number rejected: {counts['rejected']}")
    print(f"number with missing DOI and PMID: {counts['missing_doi_and_pmid']}")
    print(f"number with missing species: {counts['missing_species']}")
    print(f"number with missing mutation terms: {counts['missing_mutation_terms']}")

    if warnings:
        print("\nWarnings")
        print("--------")
        for warning in warnings:
            print(f"WARNING: {warning}")

    if errors:
        print("\nErrors")
        print("------")
        for error in errors:
            print(f"ERROR: {error}")


def main() -> int:
    if not PAPER_LINKS_PATH.exists():
        print(f"ERROR: Missing input file: {PAPER_LINKS_PATH}", file=sys.stderr)
        return 1

    fieldnames, rows = read_rows(PAPER_LINKS_PATH)
    errors, warnings, counts = validate_rows(fieldnames, rows)
    print_summary(errors, warnings, counts)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
