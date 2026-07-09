"""Extract ColabFold dashboard figure assets from completed model ZIP outputs.

Only the small ColabFold PNG summary figures are extracted. ZIP files remain
outside the repository and are not copied.
"""
from __future__ import annotations

import csv
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "structures" / "model_registry.tsv"
DASHBOARD_FIGURES_DIR = ROOT / "structures" / "dashboard_figures"
FIGURE_SUFFIXES = {
    "coverage": "_coverage.png",
    "plddt": "_plddt.png",
    "pae": "_pae.png",
}


def read_registry() -> list[dict[str, str]]:
    with REGISTRY_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_figures_for_model(row: dict[str, str]) -> list[str]:
    if row.get("model_status", "").strip() != "completed":
        return []

    zip_path = Path(row.get("external_output_zip_path", "").strip())
    if not zip_path.exists():
        return []

    model_id = row.get("model_id", "").strip()
    extracted = []
    with zipfile.ZipFile(zip_path) as zip_handle:
        names = [name for name in zip_handle.namelist() if not name.endswith("/")]
        for figure_key, suffix in FIGURE_SUFFIXES.items():
            matches = sorted(name for name in names if name.endswith(suffix))
            if not matches:
                continue
            output_path = DASHBOARD_FIGURES_DIR / f"{model_id}_{figure_key}.png"
            with zip_handle.open(matches[0]) as source, output_path.open("wb") as target:
                target.write(source.read())
            extracted.append(output_path.relative_to(ROOT).as_posix())
    return extracted


def main() -> int:
    DASHBOARD_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_registry()
    extracted = []
    missing = []
    for row in rows:
        model_id = row.get("model_id", "")
        model_figures = extract_figures_for_model(row)
        extracted.extend(model_figures)
        if row.get("model_status", "").strip() == "completed" and len(model_figures) < 3:
            missing.append(model_id)

    print(f"dashboard figures extracted: {len(extracted)}")
    for path in extracted:
        print(f"- {path}")
    if missing:
        print("completed models with missing figure types:")
        for model_id in missing:
            print(f"- {model_id}")
    else:
        print("completed models with missing figure types: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
