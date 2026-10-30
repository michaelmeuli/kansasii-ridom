# kansasii-ridom

Reads the Ridom SeqSphere+ GUI exports (`data/imm/ridom-*`), the source of truth for cgMLST/MST/clusters. Only what is easy to get wrong:

- Never edit a Ridom export and never write into `data/imm/`; derived tables go to `output/kansasii-ridom/results/`. New exports are recorded in the export log of `kansasii-docs/docs/RIDOM.md`.
- Do not rebuild anything from `kansasii-mst` / `kansasii-cgmlst` (obsolete). Ridom's scheme (~1,400 `MKAN_RS*` loci) is not the chewBBACA scheme (3,646), so distances and clusters are not comparable.
- Ridom exports are `;`-separated with a BOM-less header like `"Sample ID";"Species"`; `00` sniffs the separator. Ridom ids are `M081`, controls by `LABEL`, references `GCF_*` (`ridom_id` in `common.py`).
- Cohort = `screening_map_results.csv` `excluded` (keep, 174) or the pool (175, adds M055). The `cohort` column of `config/expected_exports.tsv` must match what was sent to Ridom: `01_export_ridom.py` packs the pool (182 with the references), not `keep`.
- `01_export_ridom.py` is the Ridom input and must stay identical to the old kansasii-cgmlst export except for M055: same names, metadata columns and unmodified assemblies. Do not strip plasmids (Ridom does). It runs via `submit_export.sbatch`, not `submit_all.sbatch`.
- A new export format gets its `id_column` filled in `config/expected_exports.tsv` before a script reads it; `00` reports UNCHECKED until then.
- SNP distances have no Ridom equivalent (RIDOM.md, "Not replaceable"); do not reintroduce them here.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy` from the repo root; tests: `pytest` in `env_ridom`. Run real jobs with `sbatch submit_all.sbatch`. Commit directly on `master`; push only when asked.
