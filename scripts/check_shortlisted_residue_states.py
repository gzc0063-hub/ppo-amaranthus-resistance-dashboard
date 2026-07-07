#!/usr/bin/env python3
"""Fetch shortlisted PPO2/PPX2 proteins and check requested residue states."""

from __future__ import annotations

import csv
import sys
import textwrap
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SHORTLIST_PATH = ROOT / "sequences" / "ppo_accession_shortlist.tsv"
FASTA_PATH = ROOT / "sequences" / "fasta" / "shortlisted_ppo2_candidates.fasta"
RESIDUE_PATH = ROOT / "sequences" / "shortlisted_residue_state_check.tsv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ACCESSIONS = [
    "ATE88443.1",
    "QIA20308.1",
    "QIA20315.1",
    "QDQ68833.1",
    "QDQ68834.1",
    "QXF78614.1",
    "ABD52326.1",
    "ABD52328.1",
    "ABD52329.1",
]

RESIDUE_COLUMNS = [
    "species",
    "accession_id",
    "sequence_length_aa",
    "residue_98",
    "residue_128",
    "residue_210_window",
    "has_glycine_210",
    "residue_361",
    "residue_399",
    "apparent_mutation_state",
    "recommended_role",
    "notes",
]


def read_shortlist() -> dict[str, dict[str, str]]:
    with SHORTLIST_PATH.open(newline="", encoding="utf-8") as handle:
        rows = {row["accession_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    missing = [accession for accession in ACCESSIONS if accession not in rows]
    if missing:
        raise ValueError(f"shortlist missing accessions: {', '.join(missing)}")
    return rows


def fetch_fasta(accessions: list[str]) -> str:
    params = urllib.parse.urlencode(
        {
            "db": "protein",
            "id": ",".join(accessions),
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    request = urllib.request.Request(
        f"{EUTILS}/efetch.fcgi?{params}",
        headers={"User-Agent": "codex-shortlisted-ppo2-residue-check/1.0"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def parse_fasta(fasta_text: str) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = ""
    chunks: list[str] = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if header:
                accession = header.split()[0].lstrip(">").split("|")[-1]
                records[accession] = (header, "".join(chunks))
            header = line
            chunks = []
        else:
            chunks.append(line.strip())
    if header:
        accession = header.split()[0].lstrip(">").split("|")[-1]
        records[accession] = (header, "".join(chunks))
    return records


def write_ordered_fasta(records: dict[str, tuple[str, str]]) -> None:
    FASTA_PATH.parent.mkdir(parents=True, exist_ok=True)
    with FASTA_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in ACCESSIONS:
            header, sequence = records[accession]
            handle.write(f"{header}\n")
            handle.write("\n".join(textwrap.wrap(sequence, width=70)))
            handle.write("\n")


def residue(sequence: str, position: int) -> str:
    if position < 1 or position > len(sequence):
        return ""
    return sequence[position - 1]


def residue_window(sequence: str, center: int = 210, flank: int = 5) -> str:
    start = max(1, center - flank)
    end = min(len(sequence), center + flank)
    return f"{start}-{end}:{sequence[start - 1:end]}"


def has_expected_glycine_210(sequence: str) -> str:
    window = residue_window(sequence).split(":", 1)[1]
    return "yes" if "CGGDP" in window else "no"


def classify_rows(rows: list[dict[str, str]]) -> None:
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)

    for species_rows in by_species.values():
        for row in species_rows:
            row["recommended_role"] = "residue_mapping_reference"

    palmeri = {row["accession_id"]: row for row in by_species["Amaranthus palmeri"]}
    if "ATE88443.1" in palmeri:
        palmeri["ATE88443.1"]["recommended_role"] = "wildtype_model_candidate"
    for accession in ("QIA20308.1", "QIA20315.1"):
        if accession in palmeri:
            palmeri[accession]["recommended_role"] = "needs_manual_check"
            palmeri[accession]["notes"] = (
                "Full-length A. palmeri PPX2L; compare publication allele context before modeling."
            )

    retroflexus = {row["accession_id"]: row for row in by_species["Amaranthus retroflexus"]}
    if "QDQ68833.1" in retroflexus:
        retroflexus["QDQ68833.1"]["recommended_role"] = "wildtype_model_candidate"

    tuberculatus = {row["accession_id"]: row for row in by_species["Amaranthus tuberculatus"]}
    if "ABD52326.1" in tuberculatus:
        tuberculatus["ABD52326.1"]["recommended_role"] = "wildtype_model_candidate"
    for row in by_species["Amaranthus tuberculatus"]:
        if row["has_glycine_210"] == "no":
            row["recommended_role"] = "mutant_allele_reference"
            row["notes"] = "A. tuberculatus candidate lacks the expected glycine in the 210-region motif; inspect alignment for DeltaG210 context."

    for row in rows:
        if int(row["sequence_length_aa"]) not in range(530, 541):
            row["recommended_role"] = "needs_manual_check"
            row["notes"] = "Sequence length differs from expected full-length PPO2/PPX2 range; numbering must be checked."


def build_residue_rows(
    shortlist_rows: dict[str, dict[str, str]],
    records: dict[str, tuple[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for accession in ACCESSIONS:
        sequence = records[accession][1]
        row = shortlist_rows[accession]
        has_glycine_210 = has_expected_glycine_210(sequence)
        state_notes: list[str] = []
        if len(sequence) != int(row["sequence_length_aa"]):
            state_notes.append("Fetched FASTA length differs from shortlist length.")

        apparent_state = "no requested-site mutation detected by direct numbering"
        if row["species"] == "Amaranthus tuberculatus" and has_glycine_210 == "no":
            apparent_state = "possible DeltaG210-related allele by 210-region motif comparison"
        if row["species"] == "Amaranthus palmeri" and residue(sequence, 399) == "A":
            apparent_state = "possible G399A allele by direct numbering"

        rows.append(
            {
                "species": row["species"],
                "accession_id": accession,
                "sequence_length_aa": str(len(sequence)),
                "residue_98": residue(sequence, 98),
                "residue_128": residue(sequence, 128),
                "residue_210_window": residue_window(sequence),
                "has_glycine_210": has_glycine_210,
                "residue_361": residue(sequence, 361),
                "residue_399": residue(sequence, 399),
                "apparent_mutation_state": apparent_state,
                "recommended_role": "needs_manual_check",
                "notes": " ".join(state_notes),
            }
        )
    classify_rows(rows)
    return rows


def write_residue_table(rows: list[dict[str, str]]) -> None:
    with RESIDUE_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESIDUE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    try:
        shortlist_rows = read_shortlist()
        fasta_text = fetch_fasta(ACCESSIONS)
        records = parse_fasta(fasta_text)
        missing = [accession for accession in ACCESSIONS if accession not in records]
        if missing:
            raise ValueError(f"NCBI FASTA response missing accessions: {', '.join(missing)}")
        write_ordered_fasta(records)
        write_residue_table(build_residue_rows(shortlist_rows, records))
    except Exception as exc:
        print(f"residue-state check failed: {exc}", file=sys.stderr)
        return 1

    print(f"FASTA saved: {FASTA_PATH}")
    print(f"Residue check table created: {RESIDUE_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
