#!/usr/bin/env python3
"""Validate WT and confirmed-mutant PPO2/PPX2 FASTA outputs."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_FASTA_PATH = ROOT / "sequences" / "fasta" / "shortlisted_ppo2_candidates.fasta"
WT_FASTA_PATH = ROOT / "sequences" / "fasta" / "final_reference_ppo2_wildtype.fasta"
MUTANT_FASTA_PATH = ROOT / "sequences" / "fasta" / "confirmed_mutant_ppo2_sequences.fasta"
MANIFEST_PATH = ROOT / "sequences" / "mutant_sequence_manifest.tsv"
MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"

MANIFEST_COLUMNS = [
    "mutant_sequence_id",
    "species",
    "reference_accession",
    "mutation",
    "reference_position",
    "wildtype_residue",
    "mutant_residue",
    "sequence_length_aa",
    "source_paper_ids",
    "mapping_status",
    "fasta_file",
    "notes",
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"missing TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def parse_fasta(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise ValueError(f"missing FASTA: {path}")
    records: dict[str, str] = {}
    header = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records[header.split()[0]] = "".join(chunks)
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header:
        records[header.split()[0]] = "".join(chunks)
    return records


def validate() -> tuple[list[str], dict[str, object]]:
    errors: list[str] = []
    source_sequences = parse_fasta(SOURCE_FASTA_PATH)
    wt_sequences = parse_fasta(WT_FASTA_PATH)
    mutant_sequences = parse_fasta(MUTANT_FASTA_PATH)
    manifest_fields, manifest_rows = read_tsv(MANIFEST_PATH)
    _, mapping_rows = read_tsv(MAPPING_PATH)

    if manifest_fields != MANIFEST_COLUMNS:
        errors.append("manifest columns are missing or out of order")

    confirmed_rows = [row for row in mapping_rows if row["mapping_status"] == "confirmed"]
    confirmed_mutation_counts = Counter(row["mutation"] for row in confirmed_rows)
    manifest_mutation_counts = Counter(row["mutation"] for row in manifest_rows)

    if len(wt_sequences) != 3:
        errors.append("WT FASTA must contain 3 reference sequences")
    if len(mutant_sequences) != len(confirmed_rows):
        errors.append("mutant FASTA count does not match confirmed mapping count")
    if len(manifest_rows) != len(confirmed_rows):
        errors.append("manifest row count does not match confirmed mapping count")
    if confirmed_mutation_counts != manifest_mutation_counts:
        errors.append("manifest mutations do not match confirmed mapping mutations")

    for row in manifest_rows:
        if row["mapping_status"] != "confirmed":
            errors.append(f"non-confirmed mapping used: {row['mutant_sequence_id']}")
        if row["mutation"] in {"R98G", "R98M"}:
            errors.append(f"manual-review R98 mutation was included: {row['mutation']}")
        if row["mutant_sequence_id"] not in mutant_sequences:
            errors.append(f"missing mutant FASTA record: {row['mutant_sequence_id']}")
            continue

        reference_sequence = source_sequences[row["reference_accession"]]
        mutant_sequence = mutant_sequences[row["mutant_sequence_id"]]
        position = int(row["reference_position"])
        wildtype = row["wildtype_residue"]
        mutant = row["mutant_residue"]

        if reference_sequence[position - 1] != wildtype:
            errors.append(f"reference residue mismatch for {row['mutant_sequence_id']}")

        if row["mutation"].startswith("Delta") or not mutant:
            if len(mutant_sequence) != len(reference_sequence) - 1:
                errors.append(f"deletion length is incorrect for {row['mutant_sequence_id']}")
            expected = reference_sequence[: position - 1] + reference_sequence[position:]
            if mutant_sequence != expected:
                errors.append(f"deletion sequence is incorrect for {row['mutant_sequence_id']}")
        else:
            if len(mutant_sequence) != len(reference_sequence):
                errors.append(f"substitution length changed for {row['mutant_sequence_id']}")
            if mutant_sequence[position - 1] != mutant:
                errors.append(f"mutant residue is incorrect for {row['mutant_sequence_id']}")
            expected = reference_sequence[: position - 1] + mutant + reference_sequence[position:]
            if mutant_sequence != expected:
                errors.append(f"substitution sequence is incorrect for {row['mutant_sequence_id']}")

    summary = {
        "wt_count": len(wt_sequences),
        "mutant_count": len(mutant_sequences),
        "included_mutations": sorted(manifest_mutation_counts),
        "excluded_mutations": sorted(
            {
                row["mutation"]
                for row in mapping_rows
                if row["mapping_status"] != "confirmed"
            }
        ),
    }
    return errors, summary


def main() -> int:
    try:
        errors, summary = validate()
    except ValueError as exc:
        print("validation result: failed")
        print(f"validation error: {exc}")
        return 1

    print(f"WT sequences created: {summary['wt_count']}")
    print(f"mutant sequences created: {summary['mutant_count']}")
    print("mutations included:")
    for mutation in summary["included_mutations"]:
        print(f"  {mutation}")
    print("mutations excluded:")
    if summary["excluded_mutations"]:
        for mutation in summary["excluded_mutations"]:
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
