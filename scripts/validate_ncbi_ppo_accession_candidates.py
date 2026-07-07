#!/usr/bin/env python3
"""Validate the NCBI PPO accession candidate manual-review TSV."""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TSV_PATH = ROOT / "sequences" / "ncbi_ppo_accession_candidates.tsv"

REQUIRED_COLUMNS = [
    "candidate_id",
    "species",
    "query",
    "database",
    "accession_id",
    "record_title",
    "sequence_length_aa",
    "linked_nucleotide_accession",
    "gene_symbol",
    "organism",
    "pubmed_id",
    "doi",
    "url",
    "candidate_reason",
    "manual_review_status",
    "notes",
]

ALLOWED_REVIEW_STATUSES = {"not_started", "needs_manual_check", "verified", "rejected"}
PPO2_PATTERN = re.compile(r"\b(PPO2|PPX2|PPX2L)\b", re.IGNORECASE)
WARNING_PATTERN = re.compile(r"\b(PPO1|PPX1)\b|chloroplastic", re.IGNORECASE)


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"Candidate TSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def row_text(row: dict[str, str]) -> str:
    return " ".join(
        [
            row.get("record_title", ""),
            row.get("gene_symbol", ""),
            row.get("candidate_reason", ""),
            row.get("notes", ""),
        ]
    )


def validate(fieldnames: list[str], rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    if fieldnames != REQUIRED_COLUMNS:
        errors.append("Required columns are missing or out of order.")

    invalid_statuses = sorted(
        {
            row.get("manual_review_status", "")
            for row in rows
            if row.get("manual_review_status", "") not in ALLOWED_REVIEW_STATUSES
        }
    )
    if invalid_statuses:
        errors.append(f"Invalid manual_review_status values: {', '.join(invalid_statuses)}")

    duplicate_count = count_duplicate_accessions(rows)
    if duplicate_count:
        errors.append(f"Duplicate accession IDs found: {duplicate_count}")
    return errors


def count_duplicate_accessions(rows: list[dict[str, str]]) -> int:
    counts = Counter(row.get("accession_id", "") for row in rows if row.get("accession_id", ""))
    return sum(count - 1 for count in counts.values() if count > 1)


def print_summary(rows: list[dict[str, str]]) -> None:
    species_counts = Counter(row.get("species", "") for row in rows)
    ppo2_count = sum(1 for row in rows if PPO2_PATTERN.search(row_text(row)))
    warning_rows = [row for row in rows if WARNING_PATTERN.search(row_text(row))]

    print(f"total candidate rows: {len(rows)}")
    print("candidate rows by species:")
    for species in sorted(species_counts):
        print(f"  {species}: {species_counts[species]}")
    print(f"duplicate accession count: {count_duplicate_accessions(rows)}")
    print(f"candidates containing PPO2/PPX2/PPX2L: {ppo2_count}")
    print(f"candidates containing PPO1/PPX1/chloroplastic warnings: {len(warning_rows)}")
    for row in warning_rows:
        print(f"  warning: {row.get('accession_id', '')} {row.get('record_title', '')}")


def main() -> int:
    try:
        fieldnames, rows = load_rows(TSV_PATH)
        errors = validate(fieldnames, rows)
    except ValueError as exc:
        print(f"validation failed: {exc}")
        return 1

    print_summary(rows)
    if errors:
        for error in errors:
            print(f"validation error: {error}")
        print("validation result: failed")
        return 1

    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
