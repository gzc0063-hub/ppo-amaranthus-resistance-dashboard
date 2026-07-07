#!/usr/bin/env python3
"""Prepare split ColabFold FASTA inputs from validated PPO2/PPX2 FASTAs."""

from __future__ import annotations

import csv
import re
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WT_FASTA_PATH = ROOT / "sequences" / "fasta" / "final_reference_ppo2_wildtype.fasta"
MUTANT_FASTA_PATH = ROOT / "sequences" / "fasta" / "confirmed_mutant_ppo2_sequences.fasta"
MANIFEST_PATH = ROOT / "sequences" / "mutant_sequence_manifest.tsv"
OUTPUT_ROOT = ROOT / "modeling_inputs" / "colabfold"
WT_OUTPUT_DIR = OUTPUT_ROOT / "wildtype"
MUTANT_OUTPUT_DIR = OUTPUT_ROOT / "mutants"
MODELING_MANIFEST_PATH = OUTPUT_ROOT / "modeling_input_manifest.tsv"

MANIFEST_COLUMNS = [
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


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header:
        records.append((header, "".join(chunks)))
    return records


def parse_header_fields(header: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for token in header.split():
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def read_manifest_by_sequence_id() -> dict[str, dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        return {
            row["mutant_sequence_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def safe_token(value: str) -> str:
    value = value.replace("Delta ", "Delta_")
    value = value.replace("/", "_")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_")


def write_single_fasta(path: Path, header: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")
        handle.write("\n".join(textwrap.wrap(sequence, width=70)))
        handle.write("\n")


def main() -> int:
    WT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MUTANT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    model_index = 1

    for header, sequence in parse_fasta(WT_FASTA_PATH):
        fields = parse_header_fields(header)
        species = fields["species"]
        accession = fields["reference_accession"]
        file_name = f"{safe_token(species)}_{safe_token(accession)}_WT.fasta"
        output_path = WT_OUTPUT_DIR / file_name
        write_single_fasta(output_path, header, sequence)
        rows.append(
            {
                "model_input_id": f"MODEL_INPUT_{model_index:03d}",
                "species": species.replace("_", " "),
                "reference_accession": accession,
                "mutation": "WT",
                "sequence_type": "wildtype",
                "sequence_length_aa": str(len(sequence)),
                "input_fasta": str(output_path.relative_to(ROOT)),
                "modeling_status": "not_started",
                "notes": "Wildtype/reference PPO2/PPX2 sequence prepared for ColabFold input.",
            }
        )
        model_index += 1

    manifest_by_id = read_manifest_by_sequence_id()
    for header, sequence in parse_fasta(MUTANT_FASTA_PATH):
        sequence_id = header.split()[0]
        manifest = manifest_by_id[sequence_id]
        species = manifest["species"]
        accession = manifest["reference_accession"]
        mutation = manifest["mutation"]
        if mutation in {"R98G", "R98M"}:
            continue
        file_name = f"{safe_token(species)}_{safe_token(accession)}_{safe_token(mutation)}_{sequence_id}.fasta"
        output_path = MUTANT_OUTPUT_DIR / file_name
        write_single_fasta(output_path, header, sequence)
        rows.append(
            {
                "model_input_id": f"MODEL_INPUT_{model_index:03d}",
                "species": species,
                "reference_accession": accession,
                "mutation": mutation,
                "sequence_type": "mutant",
                "sequence_length_aa": str(len(sequence)),
                "input_fasta": str(output_path.relative_to(ROOT)),
                "modeling_status": "not_started",
                "notes": f"Confirmed mapping mutant sequence {sequence_id} prepared for ColabFold input.",
            }
        )
        model_index += 1

    with MODELING_MANIFEST_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    wt_count = sum(1 for row in rows if row["sequence_type"] == "wildtype")
    mutant_count = sum(1 for row in rows if row["sequence_type"] == "mutant")
    print(f"wildtype FASTA files created: {wt_count}")
    print(f"mutant FASTA files created: {mutant_count}")
    print(f"total modeling inputs: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
