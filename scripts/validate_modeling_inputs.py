#!/usr/bin/env python3
"""Validate prepared ColabFold PPO2/PPX2 modeling inputs."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "modeling_inputs" / "colabfold" / "modeling_input_manifest.tsv"

REQUIRED_COLUMNS = [
    "model_input_id",
    "species",
    "reference_accession",
    "mutation",
    "sequence_type",
    "sequence_length_aa",
    "input_fasta",
    "modeling_status",
    "notes",
]
ALLOWED_SEQUENCE_TYPES = {"wildtype", "mutant"}
ALLOWED_MODELING_STATUSES = {"not_started", "submitted", "completed", "failed", "needs_manual_check"}


def parse_fasta_length(path: Path) -> int:
    sequence = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            sequence.append(line.strip())
    return len("".join(sequence))


def validate() -> tuple[list[str], dict[str, int]]:
    errors: list[str] = []
    if not MANIFEST_PATH.is_file():
        return [f"manifest missing: {MANIFEST_PATH}"], {"wildtype": 0, "mutant": 0, "total": 0}

    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    if fieldnames != REQUIRED_COLUMNS:
        errors.append("manifest columns are missing or out of order")

    for row in rows:
        if row["sequence_type"] not in ALLOWED_SEQUENCE_TYPES:
            errors.append(f"invalid sequence_type: {row['model_input_id']}")
        if row["modeling_status"] not in ALLOWED_MODELING_STATUSES:
            errors.append(f"invalid modeling_status: {row['model_input_id']}")
        if row["modeling_status"] != "not_started":
            errors.append(f"modeling_status must initially be not_started: {row['model_input_id']}")
        if row["mutation"] in {"R98G", "R98M"}:
            errors.append(f"manual-review R98 mutation included: {row['mutation']}")

        fasta_path = ROOT / row["input_fasta"]
        if not fasta_path.is_file():
            errors.append(f"input FASTA missing: {row['input_fasta']}")
            continue
        if parse_fasta_length(fasta_path) != int(row["sequence_length_aa"]):
            errors.append(f"sequence length mismatch: {row['model_input_id']}")

    summary = {
        "wildtype": sum(1 for row in rows if row.get("sequence_type") == "wildtype"),
        "mutant": sum(1 for row in rows if row.get("sequence_type") == "mutant"),
        "total": len(rows),
    }
    return errors, summary


def main() -> int:
    errors, summary = validate()
    print(f"wildtype FASTA files created: {summary['wildtype']}")
    print(f"mutant FASTA files created: {summary['mutant']}")
    print(f"total modeling inputs: {summary['total']}")
    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"validation error: {error}")
        return 1
    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
