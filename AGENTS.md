# AGENTS.md

Instructions for future Codex tasks in this repository:

- Follow the scientific rules in the project request and repository documentation.
- Never invent mutation names, residue numbers, accession numbers, DOI numbers, PubMed IDs, herbicide names, docking scores, or biological conclusions.
- Use `needs_manual_check` when uncertain.
- Do not create mutation rows unless citation metadata is available.
- Do not promote all mutation candidates automatically. A mutation may be marked
  `verified` only after manual review confirms clear species, PPO gene/isoform,
  mutation identity, resistance evidence, and citation metadata for that
  specific mutation-paper combination.
- Clearly separate experimentally confirmed resistance evidence from computational predictions.
- Do not treat docking scores as proof of field resistance.
- Do not mix waterhemp and Palmer amaranth residue numbering unless there is a verified alignment or mapping table.
- Do not mix PPO1 and PPO2 unless the gene or isoform is clearly verified.
- Prefer small, testable scripts.
- Always preserve tab-separated table schemas.
- Do not add large generated files to git.
- Do not download papers unless explicitly instructed later.
- Do not download or store PDFs in the repository.
