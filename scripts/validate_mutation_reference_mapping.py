#!/usr/bin/env python3
"""Validate verified mutation-to-reference PPO2/PPX2 mapping table."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"

REQUIRED_COLUMNS = [
    "mapping_id",
    "species",
    "gene_or_isoform",
    "mutation",
    "literature_position",
    "literature_wildtype_residue",
    "literature_mutant_residue",
    "reference_accession",
    "reference_sequence_length",
    "reference_residue_at_literature_position",
    "position_match_status",
    "alternative_position_checked",
    "alternative_reference_residue",
    "final_reference_position",
    "mapping_status",
    "notes",
]

ALLOWED_POSITION_MATCH_STATUS = {
    "matches_reference",
    "does_not_match_reference",
    "needs_manual_check",
}
ALLOWED_MAPPING_STATUS = {"confirmed", "needs_manual_check", "rejected"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def validate() -> tuple[list[str], list[dict[str, str]]]:
    errors: list[str] = []
    try:
        fieldnames, rows = read_tsv(MAPPING_PATH)
    except ValueError as exc:
        return [str(exc)], []

    if fieldnames != REQUIRED_COLUMNS:
        errors.append("mutation reference mapping columns are missing or out of order")

    if not EVIDENCE_PATH.is_file():
        errors.append(f"required input missing: {EVIDENCE_PATH}")

    for row in rows:
        if row["position_match_status"] not in ALLOWED_POSITION_MATCH_STATUS:
            errors.append(f"invalid position_match_status: {row.get('mapping_id', '')}")
        if row["mapping_status"] not in ALLOWED_MAPPING_STATUS:
            errors.append(f"invalid mapping_status: {row.get('mapping_id', '')}")

    return errors, rows


def main() -> int:
    errors, rows = validate()
    status_counts = Counter(row.get("mapping_status", "") for row in rows)
    mismatch_by_species = Counter(
        row.get("species", "")
        for row in rows
        if row.get("position_match_status") == "does_not_match_reference"
    )
    manual_review = [
        row.get("mutation", "")
        for row in rows
        if row.get("mapping_status") == "needs_manual_check"
    ]

    print(f"total mappings: {len(rows)}")
    print(f"confirmed mappings: {status_counts.get('confirmed', 0)}")
    print(f"needs_manual_check mappings: {status_counts.get('needs_manual_check', 0)}")
    print(f"rejected mappings: {status_counts.get('rejected', 0)}")
    print("mismatches by species:")
    if mismatch_by_species:
        for species, count in sorted(mismatch_by_species.items()):
            print(f"  {species}: {count}")
    else:
        print("  none")
    print("mutations needing manual review:")
    if manual_review:
        for mutation in manual_review:
            print(f"  {mutation}")
    else:
        print("  none")

    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"validation error: {error}")
        return 1

    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
