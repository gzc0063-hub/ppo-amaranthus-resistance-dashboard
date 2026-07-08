# Future Work Notes

## A. Ligand preparation

1. Herbicide ligand table with PubChem IDs  
Create a curated ligand table that links herbicide names to verified PubChem identifiers and source notes.

2. Ligand 3D structure preparation  
Prepare standardized ligand 3D structures for later docking, with provenance and preparation settings recorded.

## B. Receptor/docking setup

3. Receptor preparation from PPO models  
Prepare selected ColabFold PPO2/PPX2 receptor structures for docking while preserving model provenance.

4. Binding pocket/grid definition  
Define and document the docking search space around the relevant PPO inhibitor binding region.

5. FAD/cofactor decision  
Decide whether and how to include FAD or other cofactors in receptor preparation, with the rationale documented.

## C. Docking analysis

6. Docking validation  
Validate the docking workflow before interpreting any WT or mutant docking outputs.

7. Docking against WT and mutant proteins  
Run the validated workflow against selected WT and mutant PPO models only after setup decisions are finalized.

8. Binding affinity summary  
Summarize docking outputs after they exist, without treating docking scores as proof of field resistance.

9. Interpretation of cross-resistance predictions  
Interpret docking patterns cautiously as computational predictions, separate from experimentally confirmed resistance evidence.

## D. Dashboard integration

10. Dashboard integration of docking results  
Add validated docking result tables and visual summaries to the dashboard after docking analysis is complete.

## Suggested class presentation sentence

The next phase is to move from curated resistance evidence and predicted PPO structures into a carefully validated docking workflow, while keeping computational predictions clearly separate from experimentally confirmed resistance evidence.
