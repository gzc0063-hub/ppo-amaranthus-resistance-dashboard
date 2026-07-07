"""Validate manual mutation review decisions.

This script checks the review gate between candidate rows and the final
mutation evidence table. It does not promote rows.
"""

from collections import Counter
import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"

REVIEW_PATH = REFERENCES / "mutation_review_decisions.tsv"
CANDIDATE_PATH = REFERENCES / "mutation_candidate_table.tsv"

REQUIRED_COLUMNS = [
    "candidate_id",
    "paper_id",
    "mutation",
    "species",
    "gene_or_isoform",
    "herbicides_reported_in_paper",
    "evidence_methods",
    "experimental_confirmation",
    "current_status",
    "review_decision",
    "review_notes",
    "promote_to_final_table",
]

ALLOWED_REVIEW_DECISIONS = {"promote", "reject", "needs_more_review"}
ALLOWED_PROMOTE_VALUES = {"yes", "no"}


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def main():
    errors = []
    fieldnames, rows = read_tsv(REVIEW_PATH)
    if fieldnames != REQUIRED_COLUMNS:
        errors.append("mutation_review_decisions.tsv columns do not match required schema")

    candidate_fieldnames, candidate_rows = read_tsv(CANDIDATE_PATH)
    if "candidate_id" not in candidate_fieldnames:
        errors.append("mutation_candidate_table.tsv is missing candidate_id")
        candidate_ids = set()
    else:
        candidate_ids = {
            row["candidate_id"].strip()
            for row in candidate_rows
            if row.get("candidate_id", "").strip()
        }

    seen_candidate_ids = set()
    review_counts = Counter()
    promote_counts = Counter()

    for row_number, row in enumerate(rows, start=2):
        candidate_id = row.get("candidate_id", "").strip()
        review_decision = row.get("review_decision", "").strip()
        promote_to_final_table = row.get("promote_to_final_table", "").strip()

        if not candidate_id:
            errors.append(f"row {row_number} is missing candidate_id")
        elif candidate_id not in candidate_ids:
            errors.append(
                f"row {row_number} candidate_id is not in mutation_candidate_table.tsv: "
                f"{candidate_id}"
            )
        elif candidate_id in seen_candidate_ids:
            errors.append(f"row {row_number} duplicate candidate_id: {candidate_id}")
        seen_candidate_ids.add(candidate_id)

        if review_decision not in ALLOWED_REVIEW_DECISIONS:
            errors.append(
                f"row {row_number} has invalid review_decision: {review_decision}"
            )
        if promote_to_final_table not in ALLOWED_PROMOTE_VALUES:
            errors.append(
                f"row {row_number} has invalid promote_to_final_table: "
                f"{promote_to_final_table}"
            )

        if review_decision and promote_to_final_table:
            review_counts[review_decision] += 1
            promote_counts[promote_to_final_table] += 1

    missing_review_rows = sorted(candidate_ids - seen_candidate_ids)
    if missing_review_rows:
        errors.append(
            "mutation_review_decisions.tsv is missing candidate_id values: "
            + ", ".join(missing_review_rows)
        )

    print(f"review rows: {len(rows)}")
    print("review_decision counts:")
    for value in sorted(ALLOWED_REVIEW_DECISIONS):
        print(f"  {value}: {review_counts[value]}")
    print("promote_to_final_table counts:")
    for value in sorted(ALLOWED_PROMOTE_VALUES):
        print(f"  {value}: {promote_counts[value]}")

    if errors:
        print("validation errors:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
