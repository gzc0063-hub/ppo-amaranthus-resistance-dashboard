"""Append manually promoted mutation candidates to the final evidence table.

Rows are promoted only when mutation_review_decisions.tsv explicitly contains:
review_decision = promote and promote_to_final_table = yes.
"""

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"

REVIEW_PATH = REFERENCES / "mutation_review_decisions.tsv"
CANDIDATE_PATH = REFERENCES / "mutation_candidate_table.tsv"
LITERATURE_PATH = REFERENCES / "literature_matrix.tsv"
PDF_EVIDENCE_PATH = REFERENCES / "pdf_evidence_extraction.tsv"
FINAL_PATH = REFERENCES / "mutation_evidence_table.tsv"

REQUIRED_FINAL_COLUMNS = [
    "mutation_id",
    "species",
    "gene_or_isoform",
    "protein_accession",
    "nucleotide_accession",
    "mutation",
    "wildtype_residue",
    "mutant_residue",
    "residue_position",
    "mutation_type",
    "herbicide_class",
    "specific_herbicides_tested",
    "resistance_phenotype",
    "evidence_type",
    "citation_key",
    "paper_id",
    "doi",
    "pmid",
    "figure_or_table",
    "verification_status",
    "notes",
]

REQUIRED_PROMOTION_FIELDS = [
    "species",
    "gene_or_isoform",
    "mutation",
    "citation_key",
    "paper_id",
    "evidence_type",
    "doi",
    "notes",
]


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def row_key(row):
    return (
        row.get("paper_id", "").strip(),
        row.get("mutation", "").strip(),
        row.get("species", "").strip(),
        row.get("gene_or_isoform", "").strip(),
    )


def next_mutation_id(existing_rows):
    max_id = 0
    prefix = "MUT_"
    for row in existing_rows:
        mutation_id = row.get("mutation_id", "").strip()
        if mutation_id.startswith(prefix):
            suffix = mutation_id[len(prefix) :]
            if suffix.isdigit():
                max_id = max(max_id, int(suffix))
    return max_id + 1


def indexed_by(rows, key):
    return {row.get(key, "").strip(): row for row in rows if row.get(key, "").strip()}


def main():
    final_fieldnames, final_rows = read_tsv(FINAL_PATH)
    if final_fieldnames != REQUIRED_FINAL_COLUMNS:
        print("mutation_evidence_table.tsv columns do not match required schema")
        return 1

    _, review_rows = read_tsv(REVIEW_PATH)
    _, candidate_rows = read_tsv(CANDIDATE_PATH)
    _, literature_rows = read_tsv(LITERATURE_PATH)
    _, pdf_rows = read_tsv(PDF_EVIDENCE_PATH)

    candidates_by_id = indexed_by(candidate_rows, "candidate_id")
    literature_by_paper = indexed_by(literature_rows, "paper_id")
    pdf_by_paper = indexed_by(pdf_rows, "paper_id")
    existing_keys = {row_key(row) for row in final_rows}

    promoted = 0
    skipped = 0
    next_id = next_mutation_id(final_rows)

    for review in review_rows:
        if review.get("review_decision", "").strip() != "promote":
            continue
        if review.get("promote_to_final_table", "").strip() != "yes":
            continue

        candidate_id = review.get("candidate_id", "").strip()
        candidate = candidates_by_id.get(candidate_id)
        if not candidate:
            print(f"warning: skipping {candidate_id}; candidate_id not found")
            skipped += 1
            continue

        paper_id = candidate.get("paper_id", "").strip()
        literature = literature_by_paper.get(paper_id, {})
        pdf = pdf_by_paper.get(paper_id, {})

        final_row = {
            "mutation_id": f"MUT_{next_id:03d}",
            "species": candidate.get("species", "").strip(),
            "gene_or_isoform": candidate.get("gene_or_isoform", "").strip(),
            "protein_accession": "",
            "nucleotide_accession": "",
            "mutation": candidate.get("mutation", "").strip(),
            "wildtype_residue": candidate.get("wildtype_residue", "").strip(),
            "mutant_residue": candidate.get("mutant_residue", "").strip(),
            "residue_position": candidate.get("residue_position", "").strip(),
            "mutation_type": candidate.get("mutation_type", "").strip(),
            "herbicide_class": "PPO-inhibiting herbicide",
            "specific_herbicides_tested": candidate.get(
                "herbicides_reported_in_paper", ""
            ).strip(),
            "resistance_phenotype": pdf.get("resistance_context", "").strip(),
            "evidence_type": literature.get("evidence_type", "").strip(),
            "citation_key": candidate.get("citation_key", "").strip()
            or literature.get("citation_key", "").strip(),
            "paper_id": paper_id,
            "doi": literature.get("doi", "").strip() or pdf.get("doi", "").strip(),
            "pmid": literature.get("pmid", "").strip() or pdf.get("pmid", "").strip(),
            "figure_or_table": candidate.get("figure_or_table_reference", "").strip(),
            "verification_status": "verified",
            "notes": review.get("review_notes", "").strip()
            or "Promoted after manual review decision.",
        }

        missing = [
            field for field in REQUIRED_PROMOTION_FIELDS if not final_row.get(field, "")
        ]
        if missing:
            print(
                f"warning: skipping {candidate_id}; missing required final-table "
                f"fields: {', '.join(missing)}"
            )
            skipped += 1
            continue

        key = row_key(final_row)
        if key in existing_keys:
            print(f"warning: skipping {candidate_id}; duplicate final-table key")
            skipped += 1
            continue

        final_rows.append(final_row)
        existing_keys.add(key)
        next_id += 1
        promoted += 1

    if promoted:
        write_tsv(FINAL_PATH, final_fieldnames, final_rows)

    print(f"rows promoted: {promoted}")
    print(f"rows skipped: {skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
