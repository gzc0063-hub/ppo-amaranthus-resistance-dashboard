"""Extract dashboard PDB assets from completed ColabFold ZIP outputs.

Only the registry-selected rank_001 PDB file is extracted from each ZIP. ZIPs
remain outside the repository.
"""
from __future__ import annotations

import csv
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
DASHBOARD_MODELS_DIR = ROOT / "structures" / "dashboard_models"


def read_registry() -> tuple[list[str], list[dict[str, str]]]:
    with REGISTRY_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_registry(fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with REGISTRY_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def extract_best_pdb(row: dict[str, str]) -> str:
    if row.get("model_status", "").strip() != "completed":
        return ""

    zip_path = Path(row.get("external_output_zip_path", "").strip())
    best_model = row.get("best_model_file_name", "").strip()
    if not zip_path.exists() or not best_model:
        return ""

    output_path = DASHBOARD_MODELS_DIR / Path(best_model).name
    with zipfile.ZipFile(zip_path) as zip_handle:
        with zip_handle.open(best_model) as source, output_path.open("wb") as target:
            target.write(source.read())

    return output_path.relative_to(ROOT).as_posix()


def main() -> int:
    fieldnames, rows = read_registry()
    if "dashboard_pdb_path" not in fieldnames:
        fieldnames.append("dashboard_pdb_path")

    DASHBOARD_MODELS_DIR.mkdir(parents=True, exist_ok=True)

    extracted = []
    missing = []
    for row in rows:
        row.setdefault("dashboard_pdb_path", "")
        path = extract_best_pdb(row)
        row["dashboard_pdb_path"] = path
        if path:
            extracted.append(path)
        elif row.get("model_status", "").strip() == "completed":
            missing.append(row.get("model_id", ""))

    write_registry(fieldnames, rows)
    print(f"PDB files extracted: {len(extracted)}")
    for path in extracted:
        print(f"- {path}")
    if missing:
        print("completed models missing PDB:")
        for model_id in missing:
            print(f"- {model_id}")
    else:
        print("completed models missing PDB: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
