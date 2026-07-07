#!/usr/bin/env python3
"""Validate PPO2/PPX2 shortlist and sequence accession tables."""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_PATH = ROOT / "sequences" / "ncbi_ppo_accession_candidates.tsv"
SHORTLIST_PATH = ROOT / "sequences" / "ppo_accession_shortlist.tsv"
ACCESSION_PATH = ROOT / "sequences" / "ppo_sequence_accessions.tsv"
FINAL_REFERENCE_PATH = ROOT / "sequences" / "final_reference_sequences.tsv"
RESIDUE_MAPPING_TEMPLATE_PATH = ROOT / "sequences" / "residue_mapping_template.tsv"

SHORTLIST_COLUMNS = [
    "shortlist_id",
    "species",
    "accession_id",
    "record_title",
    "gene_symbol",
    "sequence_length_aa",
    "linked_nucleotide_accession",
    "pubmed_id",
    "url",
    "recommended_use",
    "reason",
    "review_status",
    "notes",
]

ACCESSION_COLUMNS = [
    "sequence_id",
    "species",
    "accession_id",
    "record_title",
    "gene_symbol",
    "sequence_length_aa",
    "linked_nucleotide_accession",
    "pubmed_id",
    "url",
    "verification_status",
    "used_for_modeling",
    "used_for_residue_mapping",
    "is_reference_sequence",
    "notes",
]

FINAL_REFERENCE_COLUMNS = [
    "species",
    "gene_or_isoform",
    "reference_accession",
    "linked_nucleotide_accession",
    "sequence_length_aa",
    "recommended_use",
    "verification_status",
    "notes",
]

RESIDUE_MAPPING_TEMPLATE_COLUMNS = [
    "species",
    "gene_or_isoform",
    "reference_accession",
    "reference_sequence_id",
    "numbering_confirmed",
    "mapping_status",
    "notes",
]

ALLOWED_RECOMMENDED_USE = {
    "residue_mapping_candidate",
    "modeling_candidate",
    "background_only",
    "reject",
}
ALLOWED_REVIEW_STATUS = {"needs_manual_check", "verified", "rejected"}
ALLOWED_VERIFICATION_STATUS = {"needs_manual_check", "verified", "rejected"}
FINAL_REFERENCE_ACCESSIONS = {"ATE88443.1", "ABD52326.1", "QDQ68833.1"}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def duplicate_count(rows: list[dict[str, str]], key: str) -> int:
    counts = Counter(row.get(key, "") for row in rows if row.get(key, ""))
    return sum(count - 1 for count in counts.values() if count > 1)


def validate() -> tuple[list[str], dict[str, object]]:
    candidate_fields, candidate_rows = read_tsv(CANDIDATE_PATH)
    shortlist_fields, shortlist_rows = read_tsv(SHORTLIST_PATH)
    accession_fields, accession_rows = read_tsv(ACCESSION_PATH)
    final_reference_fields, final_reference_rows = read_tsv(FINAL_REFERENCE_PATH)
    template_fields, template_rows = read_tsv(RESIDUE_MAPPING_TEMPLATE_PATH)
    del candidate_fields

    errors: list[str] = []
    if shortlist_fields != SHORTLIST_COLUMNS:
        errors.append("shortlist columns are missing or out of order")
    if accession_fields != ACCESSION_COLUMNS:
        errors.append("sequence accession columns are missing or out of order")
    if final_reference_fields != FINAL_REFERENCE_COLUMNS:
        errors.append("final reference columns are missing or out of order")
    if template_fields != RESIDUE_MAPPING_TEMPLATE_COLUMNS:
        errors.append("residue mapping template columns are missing or out of order")

    species_counts = Counter(row["species"] for row in shortlist_rows)
    if any(count > 3 for count in species_counts.values()):
        errors.append("shortlist contains more than 3 rows for at least one species")

    if duplicate_count(shortlist_rows, "accession_id"):
        errors.append("shortlist contains duplicate accession IDs")
    if duplicate_count(accession_rows, "accession_id"):
        errors.append("sequence accession table contains duplicate accession IDs")

    for row in shortlist_rows:
        if row["recommended_use"] not in ALLOWED_RECOMMENDED_USE:
            errors.append(f"invalid recommended_use: {row['accession_id']}")
        if row["review_status"] not in ALLOWED_REVIEW_STATUS:
            errors.append(f"invalid review_status: {row['accession_id']}")
        if row["review_status"] == "verified":
            errors.append(f"shortlist row is incorrectly marked verified: {row['accession_id']}")

    kept_accessions = {
        row["accession_id"]
        for row in shortlist_rows
        if row["recommended_use"] != "reject" and row["review_status"] != "rejected"
    }
    accession_accessions = {row["accession_id"] for row in accession_rows}
    if kept_accessions != accession_accessions:
        errors.append("sequence accession table does not match non-rejected shortlist accessions")

    for row in accession_rows:
        if row["verification_status"] not in ALLOWED_VERIFICATION_STATUS:
            errors.append(f"invalid verification_status: {row['accession_id']}")
        if row["used_for_modeling"] not in {"no", "yes"}:
            errors.append(f"invalid used_for_modeling: {row['accession_id']}")
        if row["used_for_residue_mapping"] != "yes":
            errors.append(f"used_for_residue_mapping must be yes: {row['accession_id']}")
        if row["is_reference_sequence"] not in {"needs_manual_check", "yes", "no"}:
            errors.append(f"invalid is_reference_sequence: {row['accession_id']}")

    verified_accessions = {row["accession_id"] for row in accession_rows if row["verification_status"] == "verified"}
    if verified_accessions != FINAL_REFERENCE_ACCESSIONS:
        errors.append("verified accession set does not match final working references")
    for row in accession_rows:
        if row["accession_id"] in FINAL_REFERENCE_ACCESSIONS:
            if row["used_for_modeling"] != "yes" or row["is_reference_sequence"] != "yes":
                errors.append(f"final reference flags are incomplete: {row['accession_id']}")
        elif row["used_for_modeling"] == "yes" or row["is_reference_sequence"] == "yes":
            errors.append(f"non-reference accession is marked as modeling/reference: {row['accession_id']}")

    final_reference_accessions = {row["reference_accession"] for row in final_reference_rows}
    if final_reference_accessions != FINAL_REFERENCE_ACCESSIONS:
        errors.append("final reference table does not contain exactly the selected references")
    if any(row["verification_status"] != "verified" for row in final_reference_rows):
        errors.append("final reference table contains non-verified rows")

    template_accessions = {row["reference_accession"] for row in template_rows}
    if template_accessions != FINAL_REFERENCE_ACCESSIONS:
        errors.append("residue mapping template does not contain exactly the selected references")
    if any(row["numbering_confirmed"] != "needs_manual_check" for row in template_rows):
        errors.append("residue mapping template numbering_confirmed must remain needs_manual_check")
    if any(row["mapping_status"] != "needs_manual_check" for row in template_rows):
        errors.append("residue mapping template mapping_status must remain needs_manual_check")

    rejected_count = len(candidate_rows) - len(kept_accessions)
    summary = {
        "species_counts": species_counts,
        "kept_accessions": sorted(kept_accessions),
        "rejected_count": rejected_count,
        "verified_accessions": sorted(FINAL_REFERENCE_ACCESSIONS),
    }
    return errors, summary


def main() -> int:
    try:
        errors, summary = validate()
    except ValueError as exc:
        print(f"validation result: failed")
        print(f"validation error: {exc}")
        return 1

    print("shortlisted rows by species:")
    for species, count in sorted(summary["species_counts"].items()):
        print(f"  {species}: {count}")
    print("accession IDs kept:")
    for accession_id in summary["kept_accessions"]:
        print(f"  {accession_id}")
    print(f"rejected count: {summary['rejected_count']}")
    print("verified reference accessions:")
    for accession_id in summary["verified_accessions"]:
        print(f"  {accession_id}")

    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"validation error: {error}")
        return 1

    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
