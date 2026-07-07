"""Validate PDF-derived paper-level evidence curation."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"
PDF_EVIDENCE_PATH = REFERENCES / "pdf_evidence_extraction.tsv"
PAPER_LINKS_PATH = REFERENCES / "paper_links.tsv"

REQUIRED_COLUMNS = [
    "paper_id",
    "pdf_file_name",
    "paper_title",
    "authors",
    "year",
    "doi",
    "pmid",
    "species_studied",
    "weed_common_name",
    "ppo_gene_or_isoform",
    "mutations_reported",
    "herbicides_tested",
    "resistance_context",
    "evidence_methods",
    "experimental_confirmation",
    "key_result_summary",
    "figure_or_table_reference",
    "page_note",
    "confidence_level",
    "needs_manual_review",
    "notes",
]

HIGH_PRIORITY_PAPER_IDS = {
    "PAPER_003",
    "PAPER_006",
    "PAPER_008",
    "PAPER_016",
    "PAPER_017",
    "PAPER_020",
    "PAPER_024",
    "PAPER_030",
    "PAPER_035",
}

ALLOWED_CONFIDENCE = {"high", "medium", "low", "needs_manual_check"}
ALLOWED_MANUAL_REVIEW = {"yes", "no", "maybe"}


def clean(value: str | None) -> str:
    return (value or "").strip()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str | None]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def main() -> int:
    errors: list[str] = []
    counts: Counter = Counter()

    if not PDF_EVIDENCE_PATH.exists():
        print(f"ERROR: Missing PDF evidence table: {PDF_EVIDENCE_PATH}", file=sys.stderr)
        return 1
    if not PAPER_LINKS_PATH.exists():
        print(f"ERROR: Missing paper links table: {PAPER_LINKS_PATH}", file=sys.stderr)
        return 1

    paper_link_fields, paper_link_rows = read_tsv(PAPER_LINKS_PATH)
    evidence_fields, evidence_rows = read_tsv(PDF_EVIDENCE_PATH)

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in evidence_fields]
    if missing_columns:
        errors.append(f"Missing required columns: {', '.join(missing_columns)}")

    if "paper_id" not in paper_link_fields:
        errors.append("paper_links.tsv is missing paper_id column.")

    paper_link_ids = {
        clean(row.get("paper_id")) for row in paper_link_rows if clean(row.get("paper_id"))
    }

    for row_number, row in enumerate(evidence_rows, start=2):
        paper_id = clean(row.get("paper_id"))
        confidence = clean(row.get("confidence_level"))
        manual_review = clean(row.get("needs_manual_review"))

        if not paper_id:
            errors.append(f"Row {row_number}: paper_id is required.")
        elif paper_id not in paper_link_ids:
            errors.append(f"Row {row_number}: paper_id '{paper_id}' is not in paper_links.tsv.")
        elif paper_id not in HIGH_PRIORITY_PAPER_IDS:
            errors.append(f"Row {row_number}: paper_id '{paper_id}' is not one of the 9 high-priority papers.")

        if confidence not in ALLOWED_CONFIDENCE:
            errors.append(f"Row {row_number}: invalid confidence_level '{confidence}'.")
        else:
            counts[f"confidence:{confidence}"] += 1

        if manual_review not in ALLOWED_MANUAL_REVIEW:
            errors.append(f"Row {row_number}: invalid needs_manual_review '{manual_review}'.")
        else:
            counts[f"needs_manual_review:{manual_review}"] += 1

        if clean(row.get("mutations_reported")):
            counts["mutations_reported_not_blank"] += 1
        if clean(row.get("ppo_gene_or_isoform")):
            counts["ppo_gene_or_isoform_not_blank"] += 1
        if clean(row.get("experimental_confirmation")) == "yes":
            counts["experimental_confirmation_yes"] += 1

    counts["total_rows"] = len(evidence_rows)

    print("PDF evidence extraction validation summary")
    print("==========================================")
    print(f"total rows: {counts['total_rows']}")
    print("confidence counts:")
    for value in ["high", "medium", "low", "needs_manual_check"]:
        print(f"  {value}: {counts[f'confidence:{value}']}")
    print("needs_manual_review counts:")
    for value in ["yes", "no", "maybe"]:
        print(f"  {value}: {counts[f'needs_manual_review:{value}']}")
    print(f"papers with mutations_reported not blank: {counts['mutations_reported_not_blank']}")
    print(f"papers with ppo_gene_or_isoform not blank: {counts['ppo_gene_or_isoform_not_blank']}")
    print(f"papers with experimental_confirmation = yes: {counts['experimental_confirmation_yes']}")

    if errors:
        print()
        print("Errors")
        print("------")
        for error in errors:
            print(f"ERROR: {error}")

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
