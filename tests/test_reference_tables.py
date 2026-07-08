from pathlib import Path
import csv
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"
SCRIPTS = ROOT / "scripts"
PROCESSED = ROOT / "data" / "processed"
SEQUENCES = ROOT / "sequences"

MODEL_REGISTRY_COLUMNS = [
    "model_id",
    "species",
    "reference_accession",
    "mutation",
    "sequence_type",
    "modeling_tool",
    "modeling_mode",
    "input_fasta",
    "external_output_zip_path",
    "best_model_file_name",
    "mean_plddt",
    "pae_file_present",
    "model_status",
    "notes",
]
ALLOWED_VERIFICATION_VALUES = {"verified", "needs_manual_check", "rejected"}

REQUIRED_COLUMNS = {
    "paper_links.tsv": [
        "paper_id",
        "paper_title",
        "authors",
        "year",
        "journal",
        "doi",
        "pmid",
        "url",
        "source_database",
        "species_mentioned",
        "ppo_gene_or_isoform",
        "mutation_terms_mentioned",
        "herbicide_terms_mentioned",
        "evidence_type",
        "verification_status",
        "notes",
    ],
    "literature_matrix.tsv": [
        "citation_key",
        "paper_id",
        "authors",
        "year",
        "title",
        "journal",
        "doi",
        "pmid",
        "species",
        "gene_or_isoform",
        "mutation_reported",
        "herbicide_or_class",
        "evidence_type",
        "method_used",
        "main_finding",
        "verified_status",
        "manual_check_notes",
    ],
    "mutation_evidence_table.tsv": [
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
    ],
    "paper_screening.tsv": [
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
    ],
    "pdf_evidence_extraction.tsv": [
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
    ],
    "mutation_candidate_table.tsv": [
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
    ],
    "mutation_review_decisions.tsv": [
        "candidate_id",
        "paper_id",
        "mutation",
        "species",
        "gene_or_isoform",
        "herbicides_reported_in_paper",
        "evidence_methods",
        "experimental_confirmation",
        "current_status",
        "review_decision",
        "review_notes",
        "promote_to_final_table",
    ],
    "required_sequence_targets.tsv": [
        "target_id",
        "species",
        "gene_or_isoform",
        "mutations_needed",
        "herbicides_linked",
        "priority",
        "sequence_needed",
        "notes",
    ],
    "required_sequence_targets_normalized.tsv": [
        "target_id",
        "species",
        "original_gene_or_isoform",
        "gene_or_isoform",
        "mutations_needed",
        "herbicides_linked",
        "priority",
        "sequence_needed",
        "notes",
    ],
    "gene_isoform_normalization.tsv": [
        "raw_gene_or_isoform",
        "standard_gene_or_isoform",
        "gene_symbol",
        "protein_name",
        "target_organelle",
        "normalization_status",
        "notes",
    ],
    "ppo_sequence_accessions.tsv": [
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
    ],
    "residue_mapping_template.tsv": [
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
    ],
}


class ReferenceTableTests(unittest.TestCase):
    def read_tsv(self, name):
        path = REFERENCES / name
        self.assertTrue(path.exists(), f"Missing reference table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return reader.fieldnames, list(reader)

    def test_paper_links_has_required_columns(self):
        fieldnames, _ = self.read_tsv("paper_links.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["paper_links.tsv"])

    def test_paper_links_template_has_same_columns_as_paper_links(self):
        fieldnames, _ = self.read_tsv("paper_links_template.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["paper_links.tsv"])

    def test_validate_paper_links_script_exists(self):
        path = SCRIPTS / "validate_paper_links.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_paper_screening_has_required_columns(self):
        fieldnames, _ = self.read_tsv("paper_screening.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["paper_screening.tsv"])

    def test_validate_paper_screening_script_exists(self):
        path = SCRIPTS / "validate_paper_screening.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_pdf_evidence_extraction_has_required_columns(self):
        fieldnames, _ = self.read_tsv("pdf_evidence_extraction.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["pdf_evidence_extraction.tsv"])

    def test_validate_pdf_evidence_extraction_script_exists(self):
        path = SCRIPTS / "validate_pdf_evidence_extraction.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_mutation_candidate_table_has_required_columns(self):
        fieldnames, _ = self.read_tsv("mutation_candidate_table.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["mutation_candidate_table.tsv"])

    def test_validate_mutation_candidate_table_script_exists(self):
        path = SCRIPTS / "validate_mutation_candidate_table.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_mutation_candidate_rows_are_not_marked_verified(self):
        _, rows = self.read_tsv("mutation_candidate_table.tsv")
        for row_number, row in enumerate(rows, start=2):
            self.assertNotEqual(
                row.get("verification_status", "").strip(),
                "verified",
                "mutation_candidate_table.tsv:"
                f"{row_number} must stay needs_manual_check until manual review "
                "confirms species, PPO isoform, mutation, resistance evidence, "
                "and citation metadata",
            )

    def test_mutation_review_decisions_has_required_columns(self):
        fieldnames, _ = self.read_tsv("mutation_review_decisions.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["mutation_review_decisions.tsv"])

    def test_validate_mutation_review_decisions_script_exists(self):
        path = SCRIPTS / "validate_mutation_review_decisions.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_promote_reviewed_mutations_script_exists(self):
        path = SCRIPTS / "promote_reviewed_mutations.py"
        self.assertTrue(path.exists(), f"Missing promotion script: {path}")

    def test_literature_matrix_has_required_columns(self):
        fieldnames, _ = self.read_tsv("literature_matrix.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["literature_matrix.tsv"])

    def test_mutation_evidence_table_has_required_columns(self):
        fieldnames, _ = self.read_tsv("mutation_evidence_table.tsv")
        self.assertEqual(fieldnames, REQUIRED_COLUMNS["mutation_evidence_table.tsv"])

    def test_validate_mutation_evidence_table_script_exists(self):
        path = SCRIPTS / "validate_mutation_evidence_table.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_build_verified_mutation_summaries_script_exists(self):
        path = SCRIPTS / "build_verified_mutation_summaries.py"
        self.assertTrue(path.exists(), f"Missing summary builder script: {path}")

    def test_verified_mutation_summary_files_can_be_generated(self):
        script = SCRIPTS / "build_verified_mutation_summaries.py"
        self.assertTrue(script.exists(), f"Missing summary builder script: {script}")
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            result.returncode,
            0,
            f"Summary builder failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )
        for name in [
            "verified_mutation_evidence.tsv",
            "mutation_summary_by_species.tsv",
            "mutation_summary_by_herbicide.tsv",
        ]:
            path = PROCESSED / name
            self.assertTrue(path.exists(), f"Missing generated summary file: {path}")
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertTrue(reader.fieldnames, f"Generated file has no header: {path}")

    def test_required_sequence_targets_has_required_columns(self):
        path = PROCESSED / "required_sequence_targets.tsv"
        self.assertTrue(path.exists(), f"Missing processed table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, REQUIRED_COLUMNS["required_sequence_targets.tsv"])

    def test_ppo_sequence_accessions_has_required_columns(self):
        path = SEQUENCES / "ppo_sequence_accessions.tsv"
        self.assertTrue(path.exists(), f"Missing sequence table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, REQUIRED_COLUMNS["ppo_sequence_accessions.tsv"])

    def test_residue_mapping_template_has_required_columns(self):
        path = SEQUENCES / "residue_mapping_template.tsv"
        self.assertTrue(path.exists(), f"Missing residue mapping table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, REQUIRED_COLUMNS["residue_mapping_template.tsv"])

    def test_validate_sequence_accessions_script_exists(self):
        path = SCRIPTS / "validate_sequence_accessions.py"
        self.assertTrue(path.exists(), f"Missing validation script: {path}")

    def test_gene_isoform_normalization_has_required_columns(self):
        path = SEQUENCES / "gene_isoform_normalization.tsv"
        self.assertTrue(path.exists(), f"Missing normalization table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(
                reader.fieldnames, REQUIRED_COLUMNS["gene_isoform_normalization.tsv"]
            )

    def test_apply_gene_isoform_normalization_script_exists(self):
        path = SCRIPTS / "apply_gene_isoform_normalization.py"
        self.assertTrue(path.exists(), f"Missing normalization script: {path}")

    def test_required_sequence_targets_normalized_has_required_columns(self):
        path = PROCESSED / "required_sequence_targets_normalized.tsv"
        self.assertTrue(path.exists(), f"Missing normalized target table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(
                reader.fieldnames,
                REQUIRED_COLUMNS["required_sequence_targets_normalized.tsv"],
            )

    def test_verification_status_values_are_controlled(self):
        for table_name in [
            "paper_links.tsv",
            "literature_matrix.tsv",
            "mutation_evidence_table.tsv",
        ]:
            _, rows = self.read_tsv(table_name)
            status_column = (
                "verified_status"
                if table_name == "literature_matrix.tsv"
                else "verification_status"
            )
            for row_number, row in enumerate(rows, start=2):
                value = row.get(status_column, "").strip()
                self.assertTrue(
                    not value or value in ALLOWED_VERIFICATION_VALUES,
                    f"{table_name}:{row_number} has invalid {status_column}: {value}",
                )

    def test_verified_rows_include_required_citation_fields(self):
        required_by_table = {
            "paper_links.tsv": ["paper_id", "paper_title", "authors", "year"],
            "literature_matrix.tsv": ["citation_key", "paper_id", "authors", "year", "title"],
            "mutation_evidence_table.tsv": ["citation_key", "paper_id"],
        }
        for table_name, required_fields in required_by_table.items():
            _, rows = self.read_tsv(table_name)
            status_column = (
                "verified_status"
                if table_name == "literature_matrix.tsv"
                else "verification_status"
            )
            for row_number, row in enumerate(rows, start=2):
                if row.get(status_column, "").strip() == "verified":
                    missing = [
                        field for field in required_fields if not row.get(field, "").strip()
                    ]
                    self.assertFalse(
                        missing,
                        f"{table_name}:{row_number} verified row is missing citation "
                        f"fields: {', '.join(missing)}",
                    )
                    has_locator = any(
                        row.get(field, "").strip()
                        for field in ("doi", "pmid", "url")
                        if field in row
                    )
                    self.assertTrue(
                        has_locator,
                        f"{table_name}:{row_number} verified row needs doi, pmid, "
                        "or url metadata",
                    )

    def test_verified_mutation_rows_include_required_scientific_fields(self):
        _, rows = self.read_tsv("mutation_evidence_table.tsv")
        required_fields = [
            "species",
            "gene_or_isoform",
            "mutation",
            "citation_key",
            "paper_id",
            "evidence_type",
            "notes",
        ]
        for row_number, row in enumerate(rows, start=2):
            if row.get("verification_status", "").strip() == "verified":
                missing = [field for field in required_fields if not row.get(field, "").strip()]
                self.assertFalse(
                    missing,
                    "mutation_evidence_table.tsv:"
                    f"{row_number} verified mutation row is missing: {', '.join(missing)}",
                )

    def test_model_registry_has_required_columns(self):
        path = ROOT / "structures" / "model_registry.tsv"
        self.assertTrue(path.exists(), f"Missing model registry table: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, MODEL_REGISTRY_COLUMNS)

    def test_colabfold_registry_scripts_exist(self):
        for script_name in [
            "inspect_colabfold_outputs.py",
            "validate_model_registry.py",
        ]:
            path = SCRIPTS / script_name
            self.assertTrue(path.exists(), f"Missing ColabFold registry script: {path}")

if __name__ == "__main__":
    unittest.main()

