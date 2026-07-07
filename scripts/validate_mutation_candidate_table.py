"""Validate conservative mutation candidate rows.

This validator checks the candidate table only. It does not promote any
candidate to verified evidence and does not read or download PDFs.
"""

from collections import Counter, defaultdict
import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"

TABLE_PATH = REFERENCES / "mutation_candidate_table.tsv"
PDF_EVIDENCE_PATH = REFERENCES / "pdf_evidence_extraction.tsv"

REQUIRED_COLUMNS = [
    "candidate_id",
    "paper_id",
    "citation_key",
    "paper_title",
    "species",
    "weed_common_name",
    "gene_or_isoform",
    "mutation",
    "wildtype_residue",
    "mutant_residue",
    "residue_position",
    "mutation_type",
    "herbicides_reported_in_paper",
    "evidence_methods",
    "experimental_confirmation",
    "figure_or_table_reference",
    "page_note",
    "confidence_level",
    "verification_status",
    "manual_review_needed",
    "notes",
]

ALLOWED_MUTATION_TYPES = {
    "substitution",
    "deletion",
    "insertion",
    "needs_manual_check",
}
ALLOWED_VERIFICATION_STATUSES = {"verified", "needs_manual_check", "rejected"}
ALLOWED_MANUAL_REVIEW_VALUES = {"yes", "no"}


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def main():
    errors = []

    fieldnames, rows = read_tsv(TABLE_PATH)
    if fieldnames != REQUIRED_COLUMNS:
        errors.append(
            "mutation_candidate_table.tsv columns do not match required schema"
        )

    pdf_fieldnames, pdf_rows = read_tsv(PDF_EVIDENCE_PATH)
    if "paper_id" not in pdf_fieldnames:
        errors.append("pdf_evidence_extraction.tsv is missing paper_id")
        pdf_paper_ids = set()
    else:
        pdf_paper_ids = {row["paper_id"].strip() for row in pdf_rows if row["paper_id"].strip()}

    mutation_counts = Counter()
    paper_counts = Counter()
    papers_by_mutation = defaultdict(set)
    rows_missing_species = []
    rows_missing_gene_or_isoform = []

    for row_number, row in enumerate(rows, start=2):
        paper_id = row.get("paper_id", "").strip()
        mutation = row.get("mutation", "").strip()
        verification_status = row.get("verification_status", "").strip()
        notes = row.get("notes", "").strip()
        mutation_type = row.get("mutation_type", "").strip()
        manual_review_needed = row.get("manual_review_needed", "").strip()
        species = row.get("species", "").strip()
        gene_or_isoform = row.get("gene_or_isoform", "").strip()

        for required_field in ["paper_id", "mutation", "verification_status", "notes"]:
            if not row.get(required_field, "").strip():
                errors.append(f"row {row_number} is missing {required_field}")

        if paper_id and paper_id not in pdf_paper_ids:
            errors.append(
                f"row {row_number} has paper_id not found in pdf evidence: {paper_id}"
            )
        if mutation_type not in ALLOWED_MUTATION_TYPES:
            errors.append(
                f"row {row_number} has invalid mutation_type: {mutation_type}"
            )
        if verification_status not in ALLOWED_VERIFICATION_STATUSES:
            errors.append(
                f"row {row_number} has invalid verification_status: "
                f"{verification_status}"
            )
        if manual_review_needed not in ALLOWED_MANUAL_REVIEW_VALUES:
            errors.append(
                f"row {row_number} has invalid manual_review_needed: "
                f"{manual_review_needed}"
            )
        if verification_status == "verified":
            errors.append(f"row {row_number} is marked verified; candidates cannot be verified")

        if mutation:
            mutation_counts[mutation] += 1
            if paper_id:
                papers_by_mutation[mutation].add(paper_id)
        if paper_id:
            paper_counts[paper_id] += 1
        if not species or species == "needs_manual_check":
            rows_missing_species.append(row_number)
        if not gene_or_isoform or gene_or_isoform == "needs_manual_check":
            rows_missing_gene_or_isoform.append(row_number)

    print(f"total candidate rows: {len(rows)}")
    print(f"unique mutations: {len(mutation_counts)}")
    print(f"unique papers: {len(paper_counts)}")
    print("mutation counts:")
    for mutation, count in sorted(mutation_counts.items()):
        papers = ", ".join(sorted(papers_by_mutation[mutation]))
        print(f"  {mutation}: {count} ({papers})")
    print("paper counts:")
    for paper_id, count in sorted(paper_counts.items()):
        print(f"  {paper_id}: {count}")
    print(f"rows with missing species: {len(rows_missing_species)}")
    if rows_missing_species:
        print("  " + ", ".join(str(row) for row in rows_missing_species))
    print(f"rows with missing gene_or_isoform: {len(rows_missing_gene_or_isoform)}")
    if rows_missing_gene_or_isoform:
        print("  " + ", ".join(str(row) for row in rows_missing_gene_or_isoform))

    if errors:
        print("validation errors:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
