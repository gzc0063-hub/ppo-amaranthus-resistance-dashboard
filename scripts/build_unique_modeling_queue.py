"""Build a unique ColabFold modeling queue from validated local inputs."""
from __future__ import annotations

import csv
from collections import OrderedDict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "modeling_inputs" / "colabfold" / "modeling_input_manifest.tsv"
MUTANT_MANIFEST_PATH = ROOT / "sequences" / "mutant_sequence_manifest.tsv"
MODEL_REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
QUEUE_PATH = ROOT / "modeling_inputs" / "colabfold" / "unique_modeling_queue.tsv"

FIELDNAMES = [
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

EXCLUDED_MUTATIONS = {"R98G", "R98M"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key_for(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("species", "").strip(),
        row.get("reference_accession", "").strip(),
        row.get("mutation", "").strip(),
    )


def sort_join(values: list[str]) -> str:
    seen = []
    for value in values:
        cleaned = value.strip()
        if cleaned and cleaned not in seen:
            seen.append(cleaned)
    return ";".join(seen)


def modeled_keys(registry_rows: list[dict[str, str]]) -> set[tuple[str, str, str]]:
    return {
        key_for(row)
        for row in registry_rows
        if row.get("model_status", "").strip() == "completed"
    }


def build_queue_rows() -> list[dict[str, str]]:
    manifest_rows = read_tsv(MANIFEST_PATH)
    mutant_rows = read_tsv(MUTANT_MANIFEST_PATH)
    registry_rows = read_tsv(MODEL_REGISTRY_PATH)
    completed_keys = modeled_keys(registry_rows)

    mutant_details: dict[tuple[str, str, str], dict[str, list[str]]] = {}
    for row in mutant_rows:
        mutation = row.get("mutation", "").strip()
        if mutation in EXCLUDED_MUTATIONS:
            continue
        key = key_for(row)
        details = mutant_details.setdefault(key, {"ids": [], "papers": []})
        details["ids"].append(row.get("mutant_sequence_id", ""))
        details["papers"].extend(row.get("source_paper_ids", "").split(";"))

    queue_by_key: OrderedDict[tuple[str, str, str], dict[str, str]] = OrderedDict()
    for row in manifest_rows:
        mutation = row.get("mutation", "").strip()
        if mutation in EXCLUDED_MUTATIONS:
            continue
        key = key_for(row)
        if key in queue_by_key:
            continue

        sequence_type = row.get("sequence_type", "").strip()
        is_completed = key in completed_keys
        details = mutant_details.get(key, {"ids": [], "papers": []})
        queue_by_key[key] = {
            "queue_id": "",
            "species": key[0],
            "reference_accession": key[1],
            "mutation": key[2],
            "sequence_type": sequence_type,
            "representative_input_fasta": row.get("input_fasta", "").strip(),
            "source_mutant_sequence_ids": sort_join(details["ids"]),
            "source_paper_ids": sort_join(details["papers"]),
            "already_modeled": "yes" if is_completed else "no",
            "modeling_priority": "low" if is_completed else "high",
            "modeling_status": "completed" if is_completed else "not_started",
            "notes": (
                "Already completed in pilot ColabFold model registry."
                if is_completed
                else "Unique confirmed mutation input pending ColabFold modeling."
            ),
        }

    rows = list(queue_by_key.values())
    for index, row in enumerate(rows, start=1):
        row["queue_id"] = f"QUEUE_{index:03d}"
    return rows


def write_queue(rows: list[dict[str, str]]) -> None:
    QUEUE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with QUEUE_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    rows = build_queue_rows()
    write_queue(rows)
    completed = sum(1 for row in rows if row["modeling_status"] == "completed")
    remaining = sum(1 for row in rows if row["modeling_status"] == "not_started")
    print(f"total unique queue rows: {len(rows)}")
    print(f"already completed rows: {completed}")
    print(f"remaining models to run: {remaining}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())