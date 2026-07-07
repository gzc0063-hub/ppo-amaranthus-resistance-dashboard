#!/usr/bin/env python3
"""Collect candidate Amaranthus PPO2/PPX2 protein accessions from NCBI.

The output is intentionally a manual-review queue. This script does not verify
accessions, build models, or alter mutation evidence tables.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = ROOT / "sequences" / "ncbi_ppo_accession_candidates.tsv"
DATABASE = "protein"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RETMAX = 20

FIELDNAMES = [
    "candidate_id",
    "species",
    "query",
    "database",
    "accession_id",
    "record_title",
    "sequence_length_aa",
    "linked_nucleotide_accession",
    "gene_symbol",
    "organism",
    "pubmed_id",
    "doi",
    "url",
    "candidate_reason",
    "manual_review_status",
    "notes",
]

QUERIES = [
    ("Amaranthus palmeri", "Amaranthus palmeri mitochondrial protoporphyrinogen oxidase"),
    ("Amaranthus palmeri", "Amaranthus palmeri PPX2 protoporphyrinogen oxidase"),
    ("Amaranthus tuberculatus", "Amaranthus tuberculatus mitochondrial protoporphyrinogen oxidase"),
    ("Amaranthus tuberculatus", "Amaranthus tuberculatus PPX2 protoporphyrinogen oxidase"),
    ("Amaranthus retroflexus", "Amaranthus retroflexus mitochondrial protoporphyrinogen oxidase"),
    ("Amaranthus retroflexus", "Amaranthus retroflexus PPX2 protoporphyrinogen oxidase"),
]

INCLUDE_TERMS = [
    "mitochondrial protoporphyrinogen",
    "protoporphyrinogen oxidase",
    "ppo2",
    "ppx2",
    "ppx2l",
]

WRONG_ISOFORM_PATTERNS = [
    re.compile(r"\bPPO1\b", re.IGNORECASE),
    re.compile(r"\bPPX1\b", re.IGNORECASE),
    re.compile(r"chloroplastic", re.IGNORECASE),
]

KNOWN_CANDIDATE = {
    "candidate_id": "",
    "species": "Amaranthus palmeri",
    "query": "manual known candidate",
    "database": DATABASE,
    "accession_id": "QIA20315.1",
    "record_title": "mitochondrial protoporphyrinogen IX oxidase [Amaranthus palmeri]",
    "sequence_length_aa": "535",
    "linked_nucleotide_accession": "MK408978.1",
    "gene_symbol": "PPX2L",
    "organism": "Amaranthus palmeri",
    "pubmed_id": "31156659",
    "doi": "",
    "url": "https://www.ncbi.nlm.nih.gov/protein/QIA20315.1",
    "candidate_reason": "Known candidate requested for manual review.",
    "manual_review_status": "needs_manual_check",
    "notes": (
        "Candidate PPX2/PPO2 sequence from resistant A. palmeri population; "
        "useful for residue mapping but must be checked before modeling."
    ),
}


def ncbi_get(endpoint: str, params: dict[str, str | int]) -> str:
    query = urllib.parse.urlencode(params)
    url = f"{EUTILS}/{endpoint}?{query}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "codex-ncbi-ppo-accession-candidate-search/1.0"},
    )
    with urllib.request.urlopen(request, timeout=20) as response:
        return response.read().decode("utf-8", errors="replace")


def esearch_ids(query: str) -> list[str]:
    payload = ncbi_get(
        "esearch.fcgi",
        {
            "db": DATABASE,
            "term": query,
            "retmode": "json",
            "retmax": RETMAX,
        },
    )
    data = json.loads(payload)
    return data.get("esearchresult", {}).get("idlist", [])


def efetch_genbank(ids: list[str]) -> str:
    return ncbi_get(
        "efetch.fcgi",
        {
            "db": DATABASE,
            "id": ",".join(ids),
            "rettype": "gb",
            "retmode": "text",
        },
    )


def fetch_pubmed_dois(pubmed_ids: set[str]) -> dict[str, str]:
    if not pubmed_ids:
        return {}
    try:
        payload = ncbi_get(
            "esummary.fcgi",
            {
                "db": "pubmed",
                "id": ",".join(sorted(pubmed_ids)),
                "retmode": "json",
            },
        )
        data = json.loads(payload)
    except (OSError, ValueError, urllib.error.URLError) as exc:
        print(f"Warning: could not fetch PubMed DOI metadata: {exc}", file=sys.stderr)
        return {}

    dois: dict[str, str] = {}
    result = data.get("result", {})
    for pubmed_id in pubmed_ids:
        summary = result.get(pubmed_id, {})
        for article_id in summary.get("articleids", []):
            if article_id.get("idtype") == "doi" and article_id.get("value"):
                dois[pubmed_id] = article_id["value"]
                break
    return dois


def normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip()


def parse_definition(lines: list[str]) -> str:
    chunks: list[str] = []
    collecting = False
    for line in lines:
        if line.startswith("DEFINITION"):
            collecting = True
            chunks.append(line[12:].strip())
            continue
        if collecting:
            if line.startswith("ACCESSION"):
                break
            chunks.append(line.strip())
    return normalize_text(" ".join(chunks))


def parse_locus_length(lines: list[str]) -> str:
    for line in lines:
        if line.startswith("LOCUS"):
            match = re.search(r"\b(\d+)\s+aa\b", line)
            if match:
                return match.group(1)
    return ""


def parse_first_line_value(lines: list[str], prefix: str) -> str:
    for line in lines:
        if line.startswith(prefix):
            return normalize_text(line[len(prefix) :])
    return ""


def parse_qualifier(record_text: str, key: str) -> str:
    match = re.search(rf"/{re.escape(key)}=\"([^\"]+)\"", record_text, re.DOTALL)
    if not match:
        return ""
    return normalize_text(match.group(1))


def parse_genbank_records(text: str, species: str, query: str) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for raw_record in text.split("\n//"):
        record_text = raw_record.strip()
        if not record_text:
            continue
        lines = record_text.splitlines()
        accession_id = ""
        version = parse_first_line_value(lines, "VERSION")
        if version:
            accession_id = version.split()[0]
        if not accession_id:
            accession = parse_first_line_value(lines, "ACCESSION")
            accession_id = accession.split()[0] if accession else ""

        title = parse_definition(lines)
        product = parse_qualifier(record_text, "product")
        gene_symbol = parse_qualifier(record_text, "gene")
        organism = parse_first_line_value(lines, "SOURCE")
        pubmed_match = re.search(r"^\s*PUBMED\s+(\d+)", record_text, re.MULTILINE)
        dbsource_match = re.search(r"DBSOURCE\s+accession\s+([A-Z]{1,4}_?\d+(?:\.\d+)?)", record_text)
        searchable = f"{title} {product}"

        candidates.append(
            {
                "candidate_id": "",
                "species": species,
                "query": query,
                "database": DATABASE,
                "accession_id": accession_id,
                "record_title": title,
                "sequence_length_aa": parse_locus_length(lines),
                "linked_nucleotide_accession": dbsource_match.group(1) if dbsource_match else "",
                "gene_symbol": gene_symbol,
                "organism": organism,
                "pubmed_id": pubmed_match.group(1) if pubmed_match else "",
                "doi": "",
                "url": f"https://www.ncbi.nlm.nih.gov/protein/{accession_id}" if accession_id else "",
                "candidate_reason": "Matched NCBI Protein search term for PPO2/PPX2 manual review.",
                "manual_review_status": "needs_manual_check",
                "notes": "",
                "_searchable": searchable,
            }
        )
    return candidates


def is_relevant(candidate: dict[str, str]) -> bool:
    searchable = candidate.get("_searchable", "").lower()
    return any(term in searchable for term in INCLUDE_TERMS)


def is_wrong_isoform(candidate: dict[str, str]) -> bool:
    searchable = candidate.get("_searchable", "")
    return any(pattern.search(searchable) for pattern in WRONG_ISOFORM_PATTERNS)


def collect_candidates() -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    warnings: list[str] = []
    pubmed_ids: set[str] = set()

    for species, query in QUERIES:
        try:
            ids = esearch_ids(query)
            time.sleep(0.34)
            if not ids:
                warnings.append(f"No NCBI Protein results for query: {query}")
                continue
            records = parse_genbank_records(efetch_genbank(ids), species, query)
            time.sleep(0.34)
        except (OSError, ValueError, urllib.error.URLError) as exc:
            warnings.append(f"NCBI request failed for query '{query}': {exc}")
            continue

        relevant = [row for row in records if is_relevant(row)]
        preferred = [row for row in relevant if not is_wrong_isoform(row)]
        if preferred:
            selected = preferred
        else:
            selected = relevant
            for row in selected:
                row["notes"] = "possible wrong isoform"

        if not selected:
            warnings.append(f"No PPO2/PPX2-like candidates retained for query: {query}")
        rows.extend(selected)
        pubmed_ids.update(row["pubmed_id"] for row in selected if row["pubmed_id"])

    rows.append(KNOWN_CANDIDATE.copy())
    pubmed_ids.add(KNOWN_CANDIDATE["pubmed_id"])

    doi_by_pubmed = fetch_pubmed_dois(pubmed_ids)
    for row in rows:
        if row["pubmed_id"] and not row["doi"]:
            row["doi"] = doi_by_pubmed.get(row["pubmed_id"], "")

    return deduplicate_rows(rows), warnings


def deduplicate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    deduped: dict[str, dict[str, str]] = {}
    for row in rows:
        accession_id = row.get("accession_id", "")
        if not accession_id:
            continue
        clean_row = {field: row.get(field, "") for field in FIELDNAMES}
        if accession_id not in deduped:
            deduped[accession_id] = clean_row
            continue
        existing = deduped[accession_id]
        if clean_row["query"] and clean_row["query"] not in existing["query"].split(" | "):
            existing["query"] = " | ".join(filter(None, [existing["query"], clean_row["query"]]))
        if clean_row["notes"] and clean_row["notes"] not in existing["notes"]:
            existing["notes"] = "; ".join(filter(None, [existing["notes"], clean_row["notes"]]))
    return list(deduped.values())


def assign_candidate_ids(rows: list[dict[str, str]]) -> None:
    rows.sort(key=lambda row: (row["species"], row["accession_id"]))
    for index, row in enumerate(rows, start=1):
        row["candidate_id"] = f"NCBI_PPO_{index:03d}"


def write_candidates(rows: list[dict[str, str]]) -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    rows, warnings = collect_candidates()
    assign_candidate_ids(rows)
    write_candidates(rows)

    print(f"Wrote {len(rows)} candidate rows to {OUTPUT_PATH}")
    for warning in warnings:
        print(f"Warning: {warning}")
    if warnings:
        print("NCBI accession collection completed with warnings; all rows require manual review.")
    else:
        print("NCBI accession collection completed; all rows require manual review.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
