# ColabFold Modeling Protocol

Use the FASTA files in `modeling_inputs/colabfold/` as inputs for ColabFold. Run one sequence at a time, or use batch mode if the active ColabFold interface supports it.

Upload the selected FASTA file, run the prediction with the chosen ColabFold settings, and save output PDB/mmCIF files plus confidence files outside the repo first. Record model output paths and confidence metrics later in a separate results table.

Do not interpret docking, binding, or herbicide-response behavior from these models until sequence validation, model confidence, and downstream validation are recorded.
