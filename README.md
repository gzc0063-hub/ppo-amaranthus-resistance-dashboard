# PPO-Inhibitor Resistance Dashboard for *Amaranthus* spp.

[Launch the dashboard](https://ppo-amaranthus-resistance-dashboard.streamlit.app/)

This repository supports a research dashboard for PPO-inhibitor target-site resistance in *Amaranthus* species. The dashboard brings together curated literature evidence, PPO2/PPX2 reference sequence mapping, mutation-to-reference residue mapping, ColabFold-predicted protein structure metadata, model confidence summaries, and interactive 3D structure viewing.

## Current Prototype Scope

Implemented:

- verified mutation evidence curation
- species and mutation summary views
- PPO2/PPX2 reference sequence and mutation mapping tables
- ColabFold model registry and model quality summaries
- interactive single-model 3D viewing
- separate side-by-side 3D model comparison
- aligned overlay comparison with mutation-site focus
- ColabFold summary figures for pLDDT, PAE, and coverage

Not implemented yet:

- herbicide ligand preparation
- receptor preparation for docking
- binding-pocket/grid definition
- docking validation
- docking against WT and mutant proteins
- binding-affinity summaries
- cross-resistance prediction interpretation

Docking results are not included in this version. Predicted protein structures are used for structural exploration only and should not be interpreted as proof of herbicide binding, field resistance, or cross-resistance.

## Public 3D Model Availability

The full ColabFold ZIP outputs and raw modeling folders are not stored in this repository. For the public dashboard, only the 12 small extracted rank_001 dashboard PDB files are published separately as GitHub Release assets:

- release tag: `dashboard-pdb-v1`
- release page: <https://github.com/gzc0063-hub/ppo-amaranthus-resistance-dashboard/releases/tag/dashboard-pdb-v1>

The Streamlit app first looks for local dashboard PDB files in `structures/dashboard_models/`. If those files are not present, it loads the matching public PDB asset from the GitHub Release.

## Run Locally

From the repository root:

```powershell
pip install -r requirements.txt
streamlit run app.py
```

If the `python` or `streamlit` commands are not on PATH, use the Python executable from your local environment:

```powershell
python -m streamlit run app.py
```

## Scientific Guardrails

- Do not invent mutation names, residue numbers, accessions, herbicides, docking scores, or biological conclusions.
- Verified resistance evidence and computational predictions must remain clearly separated.
- Do not treat docking scores or predicted structures as proof of field resistance.
- Do not mix species, PPO isoforms, or residue numbering systems unless the mapping is verified.
- Do not store paper PDFs, full ColabFold ZIP files, or large generated modeling folders in this repository.
