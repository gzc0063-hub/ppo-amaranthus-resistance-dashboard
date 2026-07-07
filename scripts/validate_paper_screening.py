"""Validate paper screening metadata without downloading paper content."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"
PAPER_LINKS_PATH = REFERENCES / "paper_links.tsv"
PAPER_SCREENING_PATH = REFERENCES / "paper_screening.tsv"

REQUIRED_COLUMNS = [
    "paper_id",
    "paper_title",
    "authors",
    "year",
    "doi",
    "url",
    "screening_status",
    "priority_level",
    "reason_for_status",
    "species_focus",
    "mentions_amaranthus",
    "mentions_waterhemp",
    "mentions_palmer_amaranth",
    "mentions_ppo",
    "mentions_target_site_resistance",
    "mentions_mutation",
    "mutation_terms_seen",
    "likely_use_for_project",
    "needs_pdf_review",
    "notes",
]

ALLOWED_SCREENING_STATUS = {
    "usable_for_mutation_table",
    "background_only",
    "not_relevant",
    "needs_manual_check",
}
ALLOWED_PRIORITY_LEVEL = {"high", "medium", "low", "needs_manual_check"}
ALLOWED_LIKELY_USE = {
    "mutation_evidence",
    "background_intro",
    "methods_reference",
    "tool_reference",
    "not_useful",
    "needs_manual_check",
}
ALLOWED_NEEDS_PDF_REVIEW = {"yes", "no", "maybe"}
ALLOWED_YES_NO = {"yes", "no", "needs_manual_check"}

YES_NO_COLUMNS = [
    "mentions_amaranthus",
    "mentions_waterhemp",
    "mentions_palmer_amaranth",
    "mentions_ppo",
    "mentions_target_site_resistance",
    "mentions_mutation",
]


def clean(value: str | None) -> str:
    return (value or "").strip()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str | None]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def validate() -> tuple[list[str], Counter, list[dict[str, str | None]]]:
    errors: list[str] = []
    counts: Counter = Counter()

    if not PAPER_LINKS_PATH.exists():
        return [f"Missing paper links file: {PAPER_LINKS_PATH}"], counts, []
    if not PAPER_SCREENING_PATH.exists():
        return [f"Missing paper screening file: {PAPER_SCREENING_PATH}"], counts, []

    _, paper_link_rows = read_tsv(PAPER_LINKS_PATH)
    fieldnames, screening_rows = read_tsv(PAPER_SCREENING_PATH)

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing_columns:
        errors.append(f"Missing required columns: {', '.join(missing_columns)}")
        return errors, counts, screening_rows

    paper_ids = {clean(row.get("paper_id")) for row in paper_link_rows if clean(row.get("paper_id"))}
    for row_number, row in enumerate(screening_rows, start=2):
        paper_id = clean(row.get("paper_id"))
        screening_status = clean(row.get("screening_status"))
        priority_level = clean(row.get("priority_level"))
        likely_use = clean(row.get("likely_use_for_project"))
        needs_pdf_review = clean(row.get("needs_pdf_review"))

        if not paper_id:
            errors.append(f"Row {row_number}: paper_id is required.")
        elif paper_id not in paper_ids:
            errors.append(f"Row {row_number}: paper_id '{paper_id}' is not in paper_links.tsv.")

        if screening_status not in ALLOWED_SCREENING_STATUS:
            errors.append(f"Row {row_number}: invalid screening_status '{screening_status}'.")
        else:
            counts[f"screening_status:{screening_status}"] += 1

        if priority_level not in ALLOWED_PRIORITY_LEVEL:
            errors.append(f"Row {row_number}: invalid priority_level '{priority_level}'.")
        else:
            counts[f"priority_level:{priority_level}"] += 1

        if likely_use not in ALLOWED_LIKELY_USE:
            errors.append(f"Row {row_number}: invalid likely_use_for_project '{likely_use}'.")
        else:
            counts[f"likely_use_for_project:{likely_use}"] += 1

        if needs_pdf_review not in ALLOWED_NEEDS_PDF_REVIEW:
            errors.append(f"Row {row_number}: invalid needs_pdf_review '{needs_pdf_review}'.")
        else:
            counts[f"needs_pdf_review:{needs_pdf_review}"] += 1

        for column in YES_NO_COLUMNS:
            value = clean(row.get(column))
            if value not in ALLOWED_YES_NO:
                errors.append(f"Row {row_number}: invalid {column} '{value}'.")

    counts["total_rows"] = len(screening_rows)
    return errors, counts, screening_rows


def print_counts(title: str, counts: Counter, prefix: str, allowed_values: list[str]) -> None:
    print(title)
    print("-" * len(title))
    for value in allowed_values:
        print(f"{value}: {counts[f'{prefix}:{value}']}")


def main() -> int:
    errors, counts, rows = validate()

    print("Paper screening validation summary")
    print("==================================")
    print(f"total rows: {counts['total_rows']}")
    print()
    print_counts(
        "Counts by screening_status",
        counts,
        "screening_status",
        ["usable_for_mutation_table", "background_only", "not_relevant", "needs_manual_check"],
    )
    print()
    print_counts(
        "Counts by priority_level",
        counts,
        "priority_level",
        ["high", "medium", "low", "needs_manual_check"],
    )
    print()
    print_counts(
        "Counts by likely_use_for_project",
        counts,
        "likely_use_for_project",
        [
            "mutation_evidence",
            "background_intro",
            "methods_reference",
            "tool_reference",
            "not_useful",
            "needs_manual_check",
        ],
    )
    print()
    print_counts(
        "Counts by needs_pdf_review",
        counts,
        "needs_pdf_review",
        ["yes", "no", "maybe"],
    )

    high_priority = [
        row for row in rows
        if clean(row.get("priority_level")) == "high"
        and clean(row.get("needs_pdf_review")) == "yes"
    ]
    print()
    print("High-priority papers needing PDF review")
    print("---------------------------------------")
    if high_priority:
        for row in high_priority:
            print(f"{clean(row.get('paper_id'))}: {clean(row.get('paper_title'))}")
    else:
        print("(none)")

    if errors:
        print()
        print("Errors")
        print("------")
        for error in errors:
            print(f"ERROR: {error}")

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
