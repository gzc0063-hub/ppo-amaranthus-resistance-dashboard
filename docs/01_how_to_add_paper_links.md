# How to Add Paper Links

Use `references/paper_links.tsv` to record paper links and citation metadata for later manual evidence extraction. Keep the file tab-separated and preserve the header exactly.

Do not download papers into this repository. Papers should stay as links, citation metadata, and manually curated notes only.

## Add Rows

1. Open `references/paper_links.tsv`.
2. Add one row per paper link under the header.
3. Assign a stable `paper_id`, such as `PAPER_001`, `PAPER_002`, and so on.
4. Paste the paper URL into the `url` column.
5. Fill citation metadata only when you know it from a reliable source.
6. Leave unknown fields blank rather than guessing.
7. Set `verification_status` to `needs_manual_check` for newly added rows.

Use `references/paper_links_template.tsv` as a formatting example only. Its example rows are placeholders and are not evidence.

## Validate Rows

Run the validator from the repository root:

```powershell
python scripts/validate_paper_links.py
```

If `python` is not available on your PATH, run it with any local Python 3 interpreter:

```powershell
path\to\python.exe scripts/validate_paper_links.py
```

The validator checks required columns, non-empty and unique `paper_id` values, non-empty URLs, controlled `verification_status` values, and common missing metadata warnings.

## Verification Status

Keep `verification_status` as `needs_manual_check` until the paper has been manually reviewed.

Use `verified` only after species, PPO gene or isoform, mutation, herbicide or resistance evidence, and citation metadata are confirmed from the source. Do not infer mutations, species, resistance phenotypes, accession numbers, DOI values, or PubMed IDs from titles or incomplete metadata.

Use `rejected` for links that have been reviewed and determined not to support this project's evidence needs.
