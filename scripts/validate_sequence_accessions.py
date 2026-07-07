"""Validate sequence accession and residue-mapping curation tables.

This script does not fetch sequences, add accessions, or build models.
"""

import csv
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]

TARGETS_PATH = ROOT / "data" / "processed" / "required_sequence_targets.tsv"
NORMALIZED_TARGETS_PATH = (
    ROOT / "data" / "processed" / "required_sequence_targets_normalized.tsv"
)
ACCESSIONS_PATH = ROOT / "sequences" / "ppo_sequence_accessions.tsv"
MAPPING_PATH = ROOT / "sequences" / "residue_mapping_template.tsv"
NORMALIZATION_PATH = ROOT / "sequences" / "gene_isoform_normalization.tsv"
NORMALIZATION_SCRIPT_PATH = ROOT / "scripts" / "apply_gene_isoform_normalization.py"

TARGET_COLUMNS = [
    "target_id",
    "species",
    "gene_or_isoform",
    "mutations_needed",
    "herbicides_linked",
    "priority",
    "sequence_needed",
    "notes",
]

NORMALIZED_TARGET_COLUMNS = [
    "target_id",
    "species",
    "original_gene_or_isoform",
    "gene_or_isoform",
    "mutations_needed",
    "herbicides_linked",
    "priority",
    "sequence_needed",
    "notes",
]

NORMALIZATION_COLUMNS = [
    "raw_gene_or_isoform",
    "standard_gene_or_isoform",
    "gene_symbol",
    "protein_name",
    "target_organelle",
    "normalization_status",
    "notes",
]

ACCESSION_COLUMNS = [
    "sequence_id",
    "species",
    "common_name",
    "gene_or_isoform",
    "accession_type",
    "accession_id",
    "database",
    "sequence_length_aa",
    "sequence_length_nt",
    "is_reference_sequence",
    "used_for_modeling",
    "used_for_residue_mapping",
    "source_paper_id",
    "source_citation_key",
    "doi",
    "url",
    "verification_status",
    "notes",
]

MAPPING_COLUMNS = [
    "mapping_id",
    "species",
    "gene_or_isoform",
    "mutation",
    "wildtype_residue",
    "residue_position",
    "mutant_residue",
    "reference_sequence_id",
    "reference_accession",
    "alignment_file",
    "numbering_confirmed",
    "mapping_status",
    "notes",
]

ALLOWED_ACCESSION_TYPES = {"protein", "nucleotide", "both", "needs_manual_check"}
ALLOWED_DATABASES = {
    "NCBI",
    "UniProt",
    "Phytozome",
    "EnsemblPlants",
    "other",
    "needs_manual_check",
}
ALLOWED_YES_NO = {"yes", "no", "needs_manual_check"}
ALLOWED_VERIFICATION_STATUSES = {"verified", "needs_manual_check", "rejected"}
ALLOWED_NUMBERING_CONFIRMED = {"yes", "no", "needs_manual_check"}
ALLOWED_MAPPING_STATUSES = {"confirmed", "needs_manual_check", "rejected"}
ALLOWED_PRIORITIES = {"high", "medium"}
ALLOWED_NORMALIZATION_STATUSES = {"accepted", "needs_manual_check", "rejected"}


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def check_columns(errors, table_name, fieldnames, expected):
    if fieldnames != expected:
        errors.append(f"{table_name} columns do not match required schema")


def main():
    errors = []
    warnings = []

    target_fieldnames, target_rows = read_tsv(TARGETS_PATH)
    normalization_fieldnames, normalization_rows = read_tsv(NORMALIZATION_PATH)
    accession_fieldnames, accession_rows = read_tsv(ACCESSIONS_PATH)
    mapping_fieldnames, mapping_rows = read_tsv(MAPPING_PATH)

    check_columns(
        errors, "required_sequence_targets.tsv", target_fieldnames, TARGET_COLUMNS
    )
    check_columns(
        errors,
        "gene_isoform_normalization.tsv",
        normalization_fieldnames,
        NORMALIZATION_COLUMNS,
    )
    check_columns(
        errors, "ppo_sequence_accessions.tsv", accession_fieldnames, ACCESSION_COLUMNS
    )
    check_columns(
        errors, "residue_mapping_template.tsv", mapping_fieldnames, MAPPING_COLUMNS
    )

    if not NORMALIZATION_SCRIPT_PATH.exists():
        errors.append(f"missing normalization script: {NORMALIZATION_SCRIPT_PATH}")
    else:
        result = subprocess.run(
            [sys.executable, str(NORMALIZATION_SCRIPT_PATH)],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            errors.append(
                "apply_gene_isoform_normalization.py failed: "
                + (result.stderr.strip() or result.stdout.strip())
            )

    normalized_fieldnames, normalized_rows = read_tsv(NORMALIZED_TARGETS_PATH)
    check_columns(
        errors,
        "required_sequence_targets_normalized.tsv",
        normalized_fieldnames,
        NORMALIZED_TARGET_COLUMNS,
    )

    target_keys = {
        (row.get("species", "").strip(), row.get("gene_or_isoform", "").strip())
        for row in target_rows
    }
    species_values = {
        row.get("species", "").strip()
        for row in target_rows
        if row.get("species", "").strip()
    }
    gene_values = {
        row.get("gene_or_isoform", "").strip()
        for row in target_rows
        if row.get("gene_or_isoform", "").strip()
    }
    normalized_gene_values = {
        row.get("gene_or_isoform", "").strip()
        for row in normalized_rows
        if row.get("gene_or_isoform", "").strip()
    }

    for row_number, row in enumerate(normalization_rows, start=2):
        status = row.get("normalization_status", "").strip()
        if status not in ALLOWED_NORMALIZATION_STATUSES:
            errors.append(
                f"gene_isoform_normalization.tsv:{row_number} "
                f"invalid normalization_status"
            )
        if not row.get("standard_gene_or_isoform", "").strip():
            errors.append(
                f"gene_isoform_normalization.tsv:{row_number} "
                "missing standard_gene_or_isoform"
            )
        if status == "needs_manual_check":
            warnings.append(
                "gene_isoform_normalization.tsv:"
                f"{row_number} normalization requires manual check for "
                f"{row.get('raw_gene_or_isoform', '').strip()}"
            )

    for row_number, row in enumerate(normalized_rows, start=2):
        if not row.get("gene_or_isoform", "").strip():
            errors.append(
                "required_sequence_targets_normalized.tsv:"
                f"{row_number} has empty normalized gene_or_isoform"
            )

    for row_number, row in enumerate(target_rows, start=2):
        priority = row.get("priority", "").strip()
        sequence_needed = row.get("sequence_needed", "").strip()
        if priority not in ALLOWED_PRIORITIES:
            errors.append(f"required_sequence_targets.tsv:{row_number} invalid priority")
        if sequence_needed not in ALLOWED_YES_NO:
            errors.append(
                f"required_sequence_targets.tsv:{row_number} invalid sequence_needed"
            )

    for row_number, row in enumerate(accession_rows, start=2):
        accession_type = row.get("accession_type", "").strip()
        database = row.get("database", "").strip()
        verification_status = row.get("verification_status", "").strip()
        if accession_type not in ALLOWED_ACCESSION_TYPES:
            errors.append(
                f"ppo_sequence_accessions.tsv:{row_number} invalid accession_type"
            )
        if database not in ALLOWED_DATABASES:
            errors.append(f"ppo_sequence_accessions.tsv:{row_number} invalid database")
        if verification_status not in ALLOWED_VERIFICATION_STATUSES:
            errors.append(
                f"ppo_sequence_accessions.tsv:{row_number} invalid verification_status"
            )
        for field in [
            "is_reference_sequence",
            "used_for_modeling",
            "used_for_residue_mapping",
        ]:
            value = row.get(field, "").strip()
            if value not in ALLOWED_YES_NO:
                errors.append(f"ppo_sequence_accessions.tsv:{row_number} invalid {field}")
        if verification_status == "verified":
            missing = [
                field
                for field in [
                    "species",
                    "gene_or_isoform",
                    "accession_id",
                    "database",
                    "accession_type",
                    "notes",
                ]
                if not row.get(field, "").strip()
            ]
            if missing:
                errors.append(
                    f"ppo_sequence_accessions.tsv:{row_number} verified accession "
                    f"missing: {', '.join(missing)}"
                )

    mappings_needing_confirmation = 0
    for row_number, row in enumerate(mapping_rows, start=2):
        key = (row.get("species", "").strip(), row.get("gene_or_isoform", "").strip())
        numbering_confirmed = row.get("numbering_confirmed", "").strip()
        mapping_status = row.get("mapping_status", "").strip()
        if key not in target_keys:
            errors.append(
                "residue_mapping_template.tsv:"
                f"{row_number} species + gene_or_isoform not in required targets"
            )
        if numbering_confirmed not in ALLOWED_NUMBERING_CONFIRMED:
            errors.append(
                f"residue_mapping_template.tsv:{row_number} invalid numbering_confirmed"
            )
        if mapping_status not in ALLOWED_MAPPING_STATUSES:
            errors.append(
                f"residue_mapping_template.tsv:{row_number} invalid mapping_status"
            )
        if (
            numbering_confirmed == "needs_manual_check"
            or mapping_status == "needs_manual_check"
        ):
            mappings_needing_confirmation += 1

    print("sequence accession validation summary")
    print(f"required sequence targets: {len(target_rows)}")
    print(f"unique species: {len(species_values)}")
    print(f"unique gene_or_isoform values: {len(gene_values)}")
    print(f"normalization rows: {len(normalization_rows)}")
    print(f"normalized sequence targets: {len(normalized_rows)}")
    print(
        "remaining unique normalized gene_or_isoform values: "
        + "; ".join(sorted(normalized_gene_values))
    )
    print(f"accession rows: {len(accession_rows)}")
    print(f"residue mapping rows: {len(mapping_rows)}")
    print(f"residue mappings needing confirmation: {mappings_needing_confirmation}")
    print(f"warning count: {len(warnings)}")

    if warnings:
        print("warnings:")
        for warning in warnings:
            print(f"  - {warning}")

    if errors:
        print("validation errors:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("validation passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
