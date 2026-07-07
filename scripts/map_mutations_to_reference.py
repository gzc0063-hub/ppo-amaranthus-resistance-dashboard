#!/usr/bin/env python3
"""Map verified PPO2/PPX2 mutation evidence rows onto final references."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_PATH = ROOT / "references" / "mutation_evidence_table.tsv"
FINAL_REFERENCE_PATH = ROOT / "sequences" / "final_reference_sequences.tsv"
FASTA_PATH = ROOT / "sequences" / "fasta" / "shortlisted_ppo2_candidates.fasta"
MAPPING_PATH = ROOT / "sequences" / "mutation_reference_mapping.tsv"

FIELDNAMES = [
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


def residue_at(sequence: str, position_text: str) -> str:
    if not position_text.isdigit():
        return ""
    position = int(position_text)
    if position < 1 or position > len(sequence):
        return ""
    return sequence[position - 1]


def mutation_parts(row: dict[str, str]) -> tuple[str, str, str]:
    wildtype = row["wildtype_residue"]
    mutant = row["mutant_residue"]
    position = row["residue_position"]
    if wildtype and position:
        return wildtype, position, mutant

    match = re.match(r"([A-Z]+)(\d+)([A-Z]+)?", row["mutation"])
    if match:
        return match.group(1), match.group(2), match.group(3) or ""
    return wildtype, position, mutant


def map_row(
    index: int,
    evidence_row: dict[str, str],
    reference_by_species: dict[str, dict[str, str]],
    sequences: dict[str, str],
) -> dict[str, str]:
    species = evidence_row["species"]
    reference = reference_by_species[species]
    accession = reference["reference_accession"]
    sequence = sequences[accession]
    wildtype, position, mutant = mutation_parts(evidence_row)
    reference_residue = residue_at(sequence, position)

    alternative_position = ""
    alternative_residue = ""
    final_position = position
    position_status = "needs_manual_check"
    mapping_status = "needs_manual_check"
    notes: list[str] = []

    is_r98 = evidence_row["mutation"] in {"R98G", "R98M"} or (wildtype == "R" and position == "98")
    if is_r98:
        alternative_position = "128"
        alternative_residue = residue_at(sequence, alternative_position)
        if reference_residue == wildtype:
            position_status = "matches_reference"
            mapping_status = "confirmed"
            final_position = position
            notes.append("R98 literature position directly matches reference residue.")
        elif alternative_residue == wildtype:
            position_status = "needs_manual_check"
            mapping_status = "needs_manual_check"
            final_position = alternative_position
            notes.append(
                "Literature R98 does not match full-length reference position 98, but position 128 is R; likely signal-peptide numbering offset, requires manual confirmation."
            )
        else:
            position_status = "does_not_match_reference"
            mapping_status = "needs_manual_check"
            notes.append("Neither position 98 nor alternative position 128 matches literature wildtype residue.")
    elif evidence_row["mutation_type"] == "deletion" and wildtype == "G" and position == "210":
        if reference_residue == "G":
            position_status = "matches_reference"
            mapping_status = "confirmed"
            notes.append("Reference has glycine at literature Delta G210 position.")
        else:
            position_status = "does_not_match_reference"
            mapping_status = "needs_manual_check"
            notes.append("Reference does not have glycine at literature Delta G210 position.")
    else:
        if reference_residue == wildtype:
            position_status = "matches_reference"
            mapping_status = "confirmed"
            notes.append("Reference residue matches literature wildtype residue at reported position.")
        else:
            position_status = "does_not_match_reference"
            mapping_status = "needs_manual_check"
            notes.append("Reference residue does not match literature wildtype residue at reported position.")

    return {
        "mapping_id": f"MUT_MAP_{index:03d}",
        "species": species,
        "gene_or_isoform": reference["gene_or_isoform"],
        "mutation": evidence_row["mutation"],
        "literature_position": position,
        "literature_wildtype_residue": wildtype,
        "literature_mutant_residue": mutant,
        "reference_accession": accession,
        "reference_sequence_length": str(len(sequence)),
        "reference_residue_at_literature_position": reference_residue,
        "position_match_status": position_status,
        "alternative_position_checked": alternative_position,
        "alternative_reference_residue": alternative_residue,
        "final_reference_position": final_position,
        "mapping_status": mapping_status,
        "notes": " ".join(notes),
    }


def main() -> int:
    evidence_rows = [
        row for row in read_tsv(EVIDENCE_PATH) if row.get("verification_status") == "verified"
    ]
    references = read_tsv(FINAL_REFERENCE_PATH)
    reference_by_species = {row["species"]: row for row in references}
    sequences = parse_fasta(FASTA_PATH)

    mappings = [
        map_row(index, row, reference_by_species, sequences)
        for index, row in enumerate(evidence_rows, start=1)
    ]

    with MAPPING_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(mappings)

    print(f"wrote mappings: {len(mappings)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
