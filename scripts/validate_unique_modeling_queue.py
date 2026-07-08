"""Validate the unique ColabFold modeling queue."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
QUEUE_PATH = ROOT / "modeling_inputs" / "colabfold" / "unique_modeling_queue.tsv"

REQUIRED_COLUMNS = [
    "queue_id",
    "species",
    "reference_accession",
    "mutation",
    "sequence_type",
    "representative_input_fasta",
    "source_mutant_sequence_ids",
    "source_paper_ids",
    "already_modeled",
    "modeling_priority",
    "modeling_status",
    "notes",
]
ALLOWED_ALREADY_MODELED = {"yes", "no"}
ALLOWED_MODELING_PRIORITY = {"high", "medium", "low"}
ALLOWED_MODELING_STATUS = {"not_started", "completed", "needs_manual_check"}
ALLOWED_SEQUENCE_TYPE = {"wildtype", "mutant"}
EXCLUDED_MUTATIONS = {"R98G", "R98M"}


def read_rows() -> tuple[list[str] | None, list[dict[str, str]]]:
    if not QUEUE_PATH.exists():
        return None, []
    with QUEUE_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def path_exists(path_text: str) -> bool:
    path = Path(path_text)
    if path.is_absolute():
        return path.exists()
    return (ROOT / path).exists()


def validate() -> tuple[list[dict[str, str]], list[str]]:
    fieldnames, rows = read_rows()
    errors = []
    if fieldnames != REQUIRED_COLUMNS:
        errors.append(f"unique_modeling_queue.tsv columns do not match required schema: {fieldnames}")
        return rows, errors

    seen_keys = set()
    seen_queue_ids = set()
    apal_v361a_completed = False
    for line_number, row in enumerate(rows, start=2):
        queue_id = row.get("queue_id", "").strip()
        if not queue_id:
            errors.append(f"line {line_number}: missing queue_id")
        elif queue_id in seen_queue_ids:
            errors.append(f"line {line_number}: duplicate queue_id {queue_id}")
        seen_queue_ids.add(queue_id)

        key = (
            row.get("species", "").strip(),
            row.get("reference_accession", "").strip(),
            row.get("mutation", "").strip(),
        )
        if key in seen_keys:
            errors.append(
                f"line {line_number}: duplicate species/reference_accession/mutation {key}"
            )
        seen_keys.add(key)

        mutation = row.get("mutation", "").strip()
        if mutation in EXCLUDED_MUTATIONS:
            errors.append(f"line {line_number}: excluded mutation present: {mutation}")

        sequence_type = row.get("sequence_type", "").strip()
        if sequence_type not in ALLOWED_SEQUENCE_TYPE:
            errors.append(f"line {line_number}: invalid sequence_type {sequence_type}")

        already_modeled = row.get("already_modeled", "").strip()
        if already_modeled not in ALLOWED_ALREADY_MODELED:
            errors.append(f"line {line_number}: invalid already_modeled {already_modeled}")

        priority = row.get("modeling_priority", "").strip()
        if priority not in ALLOWED_MODELING_PRIORITY:
            errors.append(f"line {line_number}: invalid modeling_priority {priority}")

        status = row.get("modeling_status", "").strip()
        if status not in ALLOWED_MODELING_STATUS:
            errors.append(f"line {line_number}: invalid modeling_status {status}")

        if already_modeled == "yes" and status != "completed":
            errors.append(f"line {line_number}: already_modeled=yes must have completed status")
        if status == "completed" and already_modeled != "yes":
            errors.append(f"line {line_number}: completed status must have already_modeled=yes")

        fasta_path = row.get("representative_input_fasta", "").strip()
        if not fasta_path or not path_exists(fasta_path):
            errors.append(
                f"line {line_number}: representative_input_fasta missing or not found: {fasta_path}"
            )

        if key == ("Amaranthus palmeri", "ATE88443.1", "V361A"):
            apal_v361a_completed = status == "completed" and already_modeled == "yes"

    if not apal_v361a_completed:
        errors.append("APAL V361A is not marked completed/already modeled")

    return rows, errors


def suggested_job_name(row: dict[str, str]) -> str:
    species = row["species"].replace(" ", "_")
    mutation = row["mutation"].replace(" ", "_")
    return f"{species}_{row['reference_accession']}_{mutation}"


def main() -> int:
    rows, errors = validate()
    status_counts = Counter(row.get("modeling_status", "") for row in rows)
    print(f"total unique queue rows: {len(rows)}")
    print(f"already completed rows: {status_counts.get('completed', 0)}")
    print(f"remaining models to run: {status_counts.get('not_started', 0)}")
    remaining = [suggested_job_name(row) for row in rows if row.get("modeling_status") == "not_started"]
    print("remaining model names/jobnames suggested: " + ("; ".join(remaining) if remaining else "none"))
    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"- {error}")
        return 1
    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())