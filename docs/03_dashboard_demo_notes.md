# Dashboard Demo Notes

## How to run the dashboard

Install the minimal dashboard dependencies:

```powershell
pip install -r requirements.txt
```

Run the Streamlit app from the repository root:

```powershell
streamlit run app.py
```

## What is currently implemented

The MVP dashboard reads the curated project TSV files and presents four tabs:

- Overview metrics for verified mutation evidence rows, species, mutations, completed ColabFold models, and models needing manual review.
- Mutation Evidence filters and tables from `references/mutation_evidence_table.tsv`.
- Protein Models summaries from `structures/model_registry.tsv` and `data/processed/model_quality_summary.tsv`.
- Future Work / Next Phase items with all docking-related work marked as planned.

## What is intentionally not implemented yet

This version does not run docking, does not report docking scores, and does not interpret binding affinity or field resistance from predicted structures. It also does not change biological evidence calls or copy model ZIP/PDB files into the repository.

`APAL_G399A` remains `needs_manual_check` because the corresponding ColabFold ZIP was not found in the external output folder.

## Suggested 5-minute presentation script

Start by explaining that the dashboard is a prototype for organizing PPO-inhibitor resistance evidence in Amaranthus species. Emphasize that the current version integrates manually curated literature evidence, sequence accession curation, residue mapping, and ColabFold model metadata.

In the Overview tab, point out the evidence and model counts. Use this as the project status snapshot: literature evidence is curated and model metadata is available, while docking remains a future phase.

In the Mutation Evidence tab, demonstrate filtering by species, mutation, herbicide, and verification status. Explain that verified rows are separated from candidates that still require manual review.

In the Protein Models tab, show the completed ColabFold models, mean pLDDT values, PAE availability, and the `APAL_G399A` manual-check row. State clearly that predicted structures are not being treated as resistance proof.

Close with the Future Work tab. Walk through ligand preparation, receptor setup, docking validation, WT/mutant docking comparisons, binding affinity summaries, cross-resistance interpretation, and final dashboard integration as planned next steps.
