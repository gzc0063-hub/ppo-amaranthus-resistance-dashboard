"""Create and inspect the ColabFold PPO2/PPX2 model registry.

The script reads ColabFold result ZIP files in-place and writes only metadata to
structures/model_registry.tsv. ZIP contents are not extracted into the repo.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import statistics
import zipfile


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
MANIFEST_PATH = ROOT / "modeling_inputs" / "colabfold" / "modeling_input_manifest.tsv"
QUEUE_PATH = ROOT / "modeling_inputs" / "colabfold" / "unique_modeling_queue.tsv"
QUALITY_SUMMARY_PATH = ROOT / "data" / "processed" / "model_quality_summary.tsv"
EXTERNAL_OUTPUT_DIR = Path(
    r"C:\Users\gzc0063\OneDrive - Auburn University\Documents\PPO_ColabFold_outputs_raw"
)

FIELDNAMES = [
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
    "dashboard_pdb_path",
]

QUALITY_SUMMARY_FIELDNAMES = [
    "species",
    "mutation",
    "sequence_type",
    "reference_accession",
    "mean_plddt",
    "pae_file_present",
    "model_status",
    "best_model_file_name",
]

EXPECTED_MODELS = [
    {
        "model_id": "APAL_WT",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "WT",
        "sequence_type": "wildtype",
        "zip_pattern": "APAL_WT_ATE88443_f2ebe_0.result.zip",
    },
    {
        "model_id": "ATUB_WT",
        "species": "Amaranthus tuberculatus",
        "reference_accession": "ABD52326.1",
        "mutation": "WT",
        "sequence_type": "wildtype",
        "zip_pattern": "ATUB_WT_ABD52326_8c7b7.result.zip",
    },
    {
        "model_id": "ARET_WT",
        "species": "Amaranthus retroflexus",
        "reference_accession": "QDQ68833.1",
        "mutation": "WT",
        "sequence_type": "wildtype",
        "zip_pattern": "ARET_WT_QDQ68833_b13ed.result.zip",
    },
    {
        "model_id": "APAL_V361A",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "V361A",
        "sequence_type": "mutant",
        "zip_pattern": "APAL_V361A_ATE88443_303f1.result.zip",
    },
    {
        "model_id": "APAL_DeltaG210",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "Delta G210",
        "sequence_type": "mutant",
        "zip_pattern": "APAL_DeltaG210_89355.result.zip",
    },
    {
        "model_id": "APAL_R128G",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "R128G",
        "sequence_type": "mutant",
        "zip_pattern": "APAL_R128G_96ab9.result.zip",
    },
    {
        "model_id": "APAL_R128M",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "R128M",
        "sequence_type": "mutant",
        "zip_pattern": "APAL_R128M_4a006.result.zip",
    },
    {
        "model_id": "ATUB_DeltaG210",
        "species": "Amaranthus tuberculatus",
        "reference_accession": "ABD52326.1",
        "mutation": "Delta G210",
        "sequence_type": "mutant",
        "zip_pattern": "ATUB_DeltaG210_77dd0.result.zip",
    },
    {
        "model_id": "ATUB_R128G",
        "species": "Amaranthus tuberculatus",
        "reference_accession": "ABD52326.1",
        "mutation": "R128G",
        "sequence_type": "mutant",
        "zip_pattern": "ATUB_R128G_77cf5.result.zip",
    },
    {
        "model_id": "APAL_G399A",
        "species": "Amaranthus palmeri",
        "reference_accession": "ATE88443.1",
        "mutation": "G399A",
        "sequence_type": "mutant",
        "zip_pattern": "APAL_G399A*.result.zip",
    },
    {
        "model_id": "ATUB_R128I",
        "species": "Amaranthus tuberculatus",
        "reference_accession": "ABD52326.1",
        "mutation": "R128I",
        "sequence_type": "mutant",
        "zip_pattern": "ATUB_R128I_abb70.result.zip",
    },
    {
        "model_id": "ARET_R128G",
        "species": "Amaranthus retroflexus",
        "reference_accession": "QDQ68833.1",
        "mutation": "R128G",
        "sequence_type": "mutant",
        "zip_pattern": "ARET_R128G_546f0.result.zip",
    },
]


def read_modeling_inputs() -> list[dict[str, str]]:
    rows = []
    for path in (QUEUE_PATH, MANIFEST_PATH):
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def find_input_fasta(model: dict[str, str], manifest_rows: list[dict[str, str]]) -> str:
    for row in manifest_rows:
        if (
            row.get("species") == model["species"]
            and row.get("reference_accession") == model["reference_accession"]
            and row.get("mutation") == model["mutation"]
            and row.get("sequence_type") == model["sequence_type"]
        ):
            return row.get("representative_input_fasta") or row.get("input_fasta", "")
    return ""


def mean_plddt_from_scores(zip_handle: zipfile.ZipFile, scores_name: str) -> str:
    data = json.loads(zip_handle.read(scores_name).decode("utf-8"))
    if isinstance(data.get("mean_plddt"), (int, float)):
        return f"{float(data['mean_plddt']):.2f}"
    plddt = data.get("plddt")
    if isinstance(plddt, list) and plddt:
        values = [float(value) for value in plddt]
        return f"{statistics.fmean(values):.2f}"
    return ""


def resolve_zip_path(zip_pattern: str) -> Path | None:
    matches = sorted(EXTERNAL_OUTPUT_DIR.glob(zip_pattern))
    if not matches:
        return None
    return matches[0]


def inspect_zip(zip_path: Path | None) -> dict[str, str]:
    result = {
        "best_model_file_name": "",
        "score_json_file_name": "",
        "mean_plddt": "",
        "pae_file_present": "no",
        "model_status": "needs_manual_check",
        "notes": "ColabFold ZIP not found at external path.",
    }
    if zip_path is None or not zip_path.exists():
        return result

    try:
        with zipfile.ZipFile(zip_path) as zip_handle:
            names = [name for name in zip_handle.namelist() if not name.endswith("/")]
            rank_001_pdbs = sorted(
                name for name in names if name.endswith(".pdb") and "rank_001" in name
            )
            rank_001_scores = sorted(
                name
                for name in names
                if name.endswith(".json") and "scores_rank_001" in name
            )
            pae_present = any(
                ("predicted_aligned_error" in name.lower() or "pae" in name.lower())
                for name in names
            )

            result["best_model_file_name"] = rank_001_pdbs[0] if rank_001_pdbs else ""
            result["score_json_file_name"] = (
                rank_001_scores[0] if rank_001_scores else ""
            )
            result["pae_file_present"] = "yes" if pae_present else "no"

            notes = []
            if not rank_001_pdbs:
                notes.append("rank_001 PDB not found")
            if not rank_001_scores:
                notes.append("rank_001 score JSON not found")

            if rank_001_scores:
                try:
                    result["mean_plddt"] = mean_plddt_from_scores(
                        zip_handle, rank_001_scores[0]
                    )
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    notes.append(f"mean pLDDT could not be parsed: {exc}")

            if rank_001_pdbs and rank_001_scores:
                if result["mean_plddt"]:
                    result["model_status"] = "completed"
                    notes.append(f"rank_001 score JSON: {result['score_json_file_name']}")
                    notes.append("ColabFold model metadata inspected from external ZIP.")
                else:
                    result["model_status"] = "needs_manual_check"
                    notes.append("Confidence JSON found but mean pLDDT was not parseable.")
            else:
                result["model_status"] = "needs_manual_check"

            result["notes"] = "; ".join(notes)
    except zipfile.BadZipFile:
        result["model_status"] = "failed"
        result["notes"] = "External output ZIP is not readable as a ZIP archive."

    return result


def build_registry_rows() -> list[dict[str, str]]:
    manifest_rows = read_modeling_inputs()
    rows = []
    for model in EXPECTED_MODELS:
        zip_path = resolve_zip_path(model["zip_pattern"])
        inspected = inspect_zip(zip_path)
        row = {field: "" for field in FIELDNAMES}
        row.update(
            {
                "model_id": model["model_id"],
                "species": model["species"],
                "reference_accession": model["reference_accession"],
                "mutation": model["mutation"],
                "sequence_type": model["sequence_type"],
                "modeling_tool": "ColabFold",
                "modeling_mode": "AlphaFold2_mmseqs2",
                "input_fasta": find_input_fasta(model, manifest_rows),
                "external_output_zip_path": str(zip_path) if zip_path else "",
            }
        )
        row.update({field: inspected[field] for field in FIELDNAMES if field in inspected})
        rows.append(row)
    return rows


def write_registry(rows: list[dict[str, str]]) -> None:
    REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with REGISTRY_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_quality_summary(rows: list[dict[str, str]]) -> None:
    QUALITY_SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with QUALITY_SUMMARY_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=QUALITY_SUMMARY_FIELDNAMES, delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "species": row["species"],
                    "mutation": row["mutation"],
                    "sequence_type": row["sequence_type"],
                    "reference_accession": row["reference_accession"],
                    "mean_plddt": row["mean_plddt"],
                    "pae_file_present": row["pae_file_present"],
                    "model_status": row["model_status"],
                    "best_model_file_name": row["best_model_file_name"],
                }
            )


def main() -> int:
    rows = build_registry_rows()
    write_registry(rows)
    write_quality_summary(rows)
    detected = sum(
        1
        for row in rows
        if row["external_output_zip_path"] and Path(row["external_output_zip_path"]).exists()
    )
    completed = sum(1 for row in rows if row["model_status"] == "completed")
    warnings = [row for row in rows if row["model_status"] != "completed"]
    print(f"ZIPs detected: {detected}")
    print(f"models completed: {completed}")
    if warnings:
        print("warnings:")
        for row in warnings:
            print(f"- {row['model_id']}: {row['notes']}")
    else:
        print("warnings: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
