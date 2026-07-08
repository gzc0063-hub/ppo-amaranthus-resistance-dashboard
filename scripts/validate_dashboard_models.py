"""Validate dashboard-ready 3D model assets and Streamlit viewer wiring."""
from __future__ import annotations

import ast
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_MODELS_DIR = ROOT / "structures" / "dashboard_models"
REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
APP_PATH = ROOT / "app.py"


def read_registry() -> list[dict[str, str]]:
    with REGISTRY_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_app_viewer(errors: list[str]) -> None:
    try:
        tree = ast.parse(APP_PATH.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        errors.append(f"app.py syntax error: {exc}")
        return

    imports_py3dmol = any(
        isinstance(node, ast.Import)
        and any(alias.name == "py3Dmol" for alias in node.names)
        for node in ast.walk(tree)
    )
    uses_components_html = any(
        isinstance(node, ast.Attribute)
        and node.attr == "html"
        and isinstance(node.value, ast.Name)
        and node.value.id == "components"
        for node in ast.walk(tree)
    )
    uses_py3dmol_view = any(
        isinstance(node, ast.Attribute)
        and node.attr == "view"
        and isinstance(node.value, ast.Name)
        and node.value.id == "py3Dmol"
        for node in ast.walk(tree)
    )

    if not imports_py3dmol:
        errors.append("app.py does not import py3Dmol")
    if not uses_py3dmol_view:
        errors.append("app.py does not call py3Dmol.view")
    if not uses_components_html:
        errors.append("app.py does not embed viewer HTML with Streamlit components")


def validate_dashboard_models() -> tuple[list[str], list[str], list[str]]:
    errors = []
    available = []
    missing = []

    if not DASHBOARD_MODELS_DIR.exists():
        errors.append("structures/dashboard_models/ does not exist")
        return available, missing, errors

    rows = read_registry()
    for line_number, row in enumerate(rows, start=2):
        model_id = row.get("model_id", "").strip()
        status = row.get("model_status", "").strip()
        dashboard_path = row.get("dashboard_pdb_path", "").strip()

        if model_id == "APAL_G399A" and status == "needs_manual_check" and not dashboard_path:
            missing.append(model_id)
            continue

        if status != "completed":
            continue

        if not dashboard_path:
            errors.append(f"line {line_number}: completed model missing dashboard_pdb_path")
            missing.append(model_id)
            continue

        pdb_path = ROOT / dashboard_path
        if not pdb_path.exists():
            errors.append(f"line {line_number}: dashboard PDB does not exist: {dashboard_path}")
            missing.append(model_id)
            continue
        if pdb_path.suffix.lower() != ".pdb":
            errors.append(f"line {line_number}: dashboard model is not a PDB: {dashboard_path}")
        if pdb_path.stat().st_size == 0:
            errors.append(f"line {line_number}: dashboard PDB is empty: {dashboard_path}")
            missing.append(model_id)
            continue

        available.append(model_id)

    validate_app_viewer(errors)
    return available, missing, errors


def main() -> int:
    available, missing, errors = validate_dashboard_models()
    print(f"models available for 3D dashboard viewing: {len(available)}")
    for model_id in available:
        print(f"- {model_id}")
    if missing:
        print("models missing PDB:")
        for model_id in missing:
            print(f"- {model_id}")
    else:
        print("models missing PDB: none")

    if errors:
        print("validation result: failed")
        for error in errors:
            print(f"- {error}")
        return 1

    print("validation result: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
