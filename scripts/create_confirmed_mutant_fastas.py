#!/usr/bin/env python3
"""Create WT and confirmed-mutant PPO2/PPX2 FASTA files."""

from __future__ import annotations

import csv
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FINAL_REFERENCE_PATH = ROOT / "sequences" / "final_reference_sequences.tsv"
SOURCE_FASTA_PATH = ROOT / "sequences" / "fasta" / "shortlisted_ppo2_candidates.fasta"
MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"
EVIDENCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
WT_FASTA_PATH = ROOT / "sequences" / "fasta" / "final_reference_ppo2_wildtype.fasta"
MUTANT_FASTA_PATH = ROOT / "sequences" / "fasta" / "confirmed_mutant_ppo2_sequences.fasta"
MANIFEST_PATH = ROOT / "sequences" / "mutant_sequence_manifest.tsv"

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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records[header] = "".join(chunks)
            header = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if header:
        records[header] = "".join(chunks)
    return records


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    handle.write("\n".join(textwrap.wrap(sequence, width=70)))
    handle.write("\n")


def apply_mutation(reference_sequence: str, mapping: dict[str, str]) -> tuple[str, str]:
    position = int(mapping["final_reference_position"])
    wildtype = mapping["literature_wildtype_residue"]
    mutant = mapping["literature_mutant_residue"]
    observed = reference_sequence[position - 1]
    if observed != wildtype:
        raise ValueError(
            f"{mapping['mapping_id']} expected {wildtype} at {position}, found {observed}"
        )

    if mapping["mutation"].startswith("Delta") or not mutant:
        return reference_sequence[: position - 1] + reference_sequence[position:], "deletion"
    return reference_sequence[: position - 1] + mutant + reference_sequence[position:], "substitution"


def paper_ids_by_verified_mapping() -> dict[str, str]:
    if not EVIDENCE_PATH.is_file():
        return {}
    verified_rows = [
        row for row in read_tsv(EVIDENCE_PATH) if row.get("verification_status") == "verified"
    ]
    return {
        f"MUT_MAP_{index:03d}": row.get("paper_id", "")
        for index, row in enumerate(verified_rows, start=1)
    }


def safe_value(value: str) -> str:
    return value.replace(" ", "_").replace(";", ",")


def main() -> int:
    references = read_tsv(FINAL_REFERENCE_PATH)
    source_sequences = parse_fasta(SOURCE_FASTA_PATH)
    mappings = [
        row for row in read_tsv(MAPPING_PATH) if row.get("mapping_status") == "confirmed"
    ]
    paper_ids = paper_ids_by_verified_mapping()

    WT_FASTA_PATH.parent.mkdir(parents=True, exist_ok=True)
    with WT_FASTA_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        for row in references:
            accession = row["reference_accession"]
            header = (
                f"species={safe_value(row['species'])} "
                f"reference_accession={accession} "
                f"gene_or_isoform={safe_value(row['gene_or_isoform'])} "
                "status=wildtype_reference"
            )
            write_fasta_record(handle, header, source_sequences[accession])

    manifest_rows: list[dict[str, str]] = []
    with MUTANT_FASTA_PATH.open("w", encoding="utf-8", newline="\n") as fasta_handle:
        for index, mapping in enumerate(mappings, start=1):
            mutant_sequence_id = f"MUTSEQ_{index:03d}"
            reference_sequence = source_sequences[mapping["reference_accession"]]
            mutant_sequence, mutation_kind = apply_mutation(reference_sequence, mapping)
            paper_id = paper_ids.get(mapping["mapping_id"], "")
            header = (
                f"{mutant_sequence_id} "
                f"species={safe_value(mapping['species'])} "
                f"reference_accession={mapping['reference_accession']} "
                f"mutation={mapping['mutation']} "
                f"source_position={mapping['final_reference_position']} "
                f"paper_id={paper_id or 'not_available'} "
                "status=confirmed_mapping"
            )
            write_fasta_record(fasta_handle, header, mutant_sequence)
            manifest_rows.append(
                {
                    "mutant_sequence_id": mutant_sequence_id,
                    "species": mapping["species"],
                    "reference_accession": mapping["reference_accession"],
                    "mutation": mapping["mutation"],
                    "reference_position": mapping["final_reference_position"],
                    "wildtype_residue": mapping["literature_wildtype_residue"],
                    "mutant_residue": mapping["literature_mutant_residue"],
                    "sequence_length_aa": str(len(mutant_sequence)),
                    "source_paper_ids": paper_id,
                    "mapping_status": mapping["mapping_status"],
                    "fasta_file": str(MUTANT_FASTA_PATH.relative_to(ROOT)),
                    "notes": f"Created from {mapping['mapping_id']} confirmed {mutation_kind} mapping.",
                }
            )

    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"WT sequences created: {len(references)}")
    print(f"mutant sequences created: {len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
