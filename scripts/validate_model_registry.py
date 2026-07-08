"""Validate structures/model_registry.tsv for pilot ColabFold models."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"

REQUIRED_COLUMNS = [
    "model_id",
    "species",
    "reference_accession",
    "mutation",
    "sequence_type",
    "modeling_tool",
    "modeling_mode",
    "input_fasta",
    "external_output_zip_path",
    "best_model_file_name",
    "mean_plddt",
    "pae_file_present",
    "model_status",
    "notes",
]
ALLOWED_SEQUENCE_TYPES = {"wildtype", "mutant"}
ALLOWED_MODEL_STATUSES = {"completed", "failed", "needs_manual_check", "not_started"}
ALLOWED_PAE_VALUES = {"yes", "no", ""}


def read_rows() -> tuple[list[str] | None, list[dict[str, str]]]:
    if not REGISTRY_PATH.exists():
        return None, []
    with REGISTRY_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def validate() -> tuple[list[dict[str, str]], list[str]]:
    fieldnames, rows = read_rows()
    errors = []
    if fieldnames != REQUIRED_COLUMNS:
        errors.append(f"model_registry.tsv columns do not match required schema: {fieldnames}")
        return rows, errors

    seen_ids = set()
    for line_number, row in enumerate(rows, start=2):
        model_id = row.get("model_id", "").strip()
        if not model_id:
            errors.append(f"line {line_number}: missing model_id")
        elif model_id in seen_ids:
            errors.append(f"line {line_number}: duplicate model_id {model_id}")
        seen_ids.add(model_id)

        sequence_type = row.get("sequence_type", "").strip()
        if sequence_type not in ALLOWED_SEQUENCE_TYPES:
            errors.append(f"line {line_number}: invalid sequence_type {sequence_type}")

        status = row.get("model_status", "").strip()
        if status not in ALLOWED_MODEL_STATUSES:
            errors.append(f"line {line_number}: invalid model_status {status}")

        pae_value = row.get("pae_file_present", "").strip()
        if pae_value not in ALLOWED_PAE_VALUES:
            errors.append(f"line {line_number}: invalid pae_file_present {pae_value}")

        zip_value = row.get("external_output_zip_path", "").strip()
        zip_path = Path(zip_value) if zip_value else None
        if status == "completed" and (zip_path is None or not zip_path.exists()):
            errors.append(
                f"line {line_number}: completed row external ZIP path does not exist: "
                f"{zip_value}"
            )

        if status == "completed" and not row.get("best_model_file_name", "").strip():
            errors.append(f"line {line_number}: completed row missing best_model_file_name")
        if status == "completed" and not row.get("mean_plddt", "").strip():
            errors.append(f"line {line_number}: completed row missing mean_plddt")

    return rows, errors


def main() -> int:
    rows, errors = validate()
    status_counts = Counter(row.get("model_status", "") for row in rows)
    print(f"total model rows: {len(rows)}")
    print(f"completed model rows: {status_counts.get('completed', 0)}")
    print(f"needs_manual_check model rows: {status_counts.get('needs_manual_check', 0)}")
    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"- {error}")
        return 1
    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
