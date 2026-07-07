"""Apply controlled PPO gene/isoform label normalization.

This script does not fetch accessions or sequences. It creates a normalized
copy of required_sequence_targets.tsv while preserving the original label.
"""

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
NORMALIZATION_PATH = ROOT / "sequences" / "gene_isoform_normalization.tsv"
INPUT_PATH = ROOT / "data" / "processed" / "required_sequence_targets.tsv"
OUTPUT_PATH = ROOT / "data" / "processed" / "required_sequence_targets_normalized.tsv"

OUTPUT_COLUMNS = [
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


def read_rows(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    normalization_rows = read_rows(NORMALIZATION_PATH)
    target_rows = read_rows(INPUT_PATH)
    normalization = {
        row.get("raw_gene_or_isoform", "").strip(): row.get(
            "standard_gene_or_isoform", ""
        ).strip()
        for row in normalization_rows
        if row.get("raw_gene_or_isoform", "").strip()
    }

    output_rows = []
    for row in target_rows:
        original_gene = row.get("gene_or_isoform", "").strip()
        normalized_gene = normalization.get(original_gene, original_gene)
        output_rows.append(
            {
                "target_id": row.get("target_id", "").strip(),
                "species": row.get("species", "").strip(),
                "original_gene_or_isoform": original_gene,
                "gene_or_isoform": normalized_gene,
                "mutations_needed": row.get("mutations_needed", "").strip(),
                "herbicides_linked": row.get("herbicides_linked", "").strip(),
                "priority": row.get("priority", "").strip(),
                "sequence_needed": row.get("sequence_needed", "").strip(),
                "notes": row.get("notes", "").strip(),
            }
        )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    unique_normalized = sorted(
        {row["gene_or_isoform"] for row in output_rows if row["gene_or_isoform"]}
    )
    print(f"normalized sequence targets created: {len(output_rows)}")
    print(
        "remaining unique normalized gene_or_isoform values: "
        + "; ".join(unique_normalized)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
