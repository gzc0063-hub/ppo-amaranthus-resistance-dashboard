# Dashboard Demo Notes

## How to run the dashboard

Install the dashboard dependencies:

```powershell
pip install -r requirements.txt
```

Run the Streamlit app from the repository root:

```powershell
streamlit run app.py
```

## New UI/UX highlights

- Polished hero section with clear prototype and in-progress status badges.
- Sidebar quick-jump panel for Overview, Mutation Evidence, Protein Models, Current Status, and Future Work / Next Phase.
- Presentation-ready metric cards for evidence rows, species, mutations, completed models, and manual-check models.
- Narrative overview cards explaining the objective, why Amaranthus resistance matters, and the biological question.
- Mutation Evidence tab with filters, readable tables, and simple evidence-count charts.
- Protein Models tab with improved model filtering, selected-model metadata cards, and interactive py3Dmol 3D viewing.
- Dedicated Current Status tab separating completed work, in-progress work, and not-yet-implemented work.
- Future Work tab grouped into ligand preparation, receptor/docking setup, docking analysis, and decision-support integration.

## What to show during class

Start on the Overview tab. Point to the badges first: this is a prototype, it is in progress, and docking is not yet implemented. Then use the metric cards to summarize the current project state.

Move to Mutation Evidence and demonstrate a filter by species or mutation. Explain that this section is literature-based and that biological evidence calls are not being changed by the dashboard.

Move to Protein Models and select a completed model. Show the metadata card, then rotate or zoom the 3D structure viewer. Explain that the structures are ColabFold predictions used for exploration, not proof of resistance or binding affinity.

Open Current Status to make the project boundary clear: evidence curation, mapping, model generation, and visualization are complete; ligand preparation and docking are planned next.

End with Future Work / Next Phase. Use it as the transition from the current prototype to the planned docking and decision-support workflow.

## How to explain the project status

Suggested wording:

This dashboard is already functional as a research prototype. It integrates curated PPO target-site resistance evidence, mutation-to-reference mapping, and predicted PPO2/PPX2 protein structures for Amaranthus species. The next phase is docking, but docking results and binding-affinity interpretation are intentionally not included yet.

## What is intentionally not implemented yet

This version does not run docking, does not report docking scores, and does not interpret binding affinity or field resistance from predicted structures. It also does not change biological evidence calls or copy ColabFold ZIP files into the repository.

`APAL_G399A` remains `needs_manual_check` because the corresponding completed model ZIP/PDB is not available for dashboard visualization.

## Suggested 5-minute presentation script

Start by saying the project is a prototype dashboard for PPO-inhibitor resistance in Amaranthus species. Emphasize that it is evidence-first and separates confirmed literature evidence from future computational predictions.

In Overview, explain the project objective, why waterhemp and Palmer amaranth matter, and what the current metrics show.

In Mutation Evidence, filter the table and show that the evidence is traceable to curated literature rows.

In Protein Models, choose one model and show the interactive 3D structure. Say clearly that these are predicted structures from ColabFold and that docking is not implemented in this version.

In Current Status, explain what is complete, what is in progress, and what is not yet implemented.

Close in Future Work by explaining that the next scientific step is to add ligand preparation, receptor setup, validated docking, and cautious interpretation of possible cross-resistance risk.
