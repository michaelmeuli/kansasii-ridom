# kansasii-ridom

Reads the Ridom SeqSphere+ GUI exports (`data/imm/ridom-*`), the source of truth for cgMLST/MST/clusters. Only what is easy to get wrong:

- Never edit a Ridom export and never write into `data/imm/`; derived tables go to `output/kansasii-ridom/results/`. New exports are recorded in the export log of `kansasii-docs/docs/RIDOM.md`.
- Do not rebuild anything from `kansasii-mst` / `kansasii-cgmlst` (not used for the report). Ridom's scheme (~1,400 `MKAN_RS*` loci) is not the chewBBACA scheme (3,647), so distances and clusters are not comparable.
- Ridom exports are `;`-separated with a BOM-less header like `"Sample ID";"Species"`; `00` sniffs the separator. Ridom ids are `M081`, controls by `LABEL`, references `GCF_*` (`ridom_id` in `common.py`).
- M055 is sent to Ridom (pool) but Ridom leaves it out of the distances and clusters, so the derived tables (species table, paper numbers) count only the `keep` cohort; M055 stays in `qc.tsv` as `pool_only`.
- Cohort = `screening_map_results.csv` `excluded` (keep, 174) or the pool (175, adds M055). The `cohort` column of `config/expected_exports.tsv` must match what was sent to Ridom: `01_export_ridom.py` packs the pool (182 with the references), not `keep`.
- `01_export_ridom.py` is the Ridom input: keep its names (`M001` ..., controls by `LABEL`, references by accession), metadata columns and unmodified assemblies stable, or the Ridom project no longer matches. Do not strip plasmids (Ridom does): the references are the unmodified `<acc>_query.fasta`, so GCF_000157895.3 (ATCC 12478) has 2 contigs, 6,577,228 bp. It runs via `submit_export.sbatch`, not `submit_all.sbatch`.
- A new export format gets its `id_column` filled in `config/expected_exports.tsv` before a script reads it; `00` reports UNCHECKED until then. Do not ask for a GUI export that a script can derive: the list holds only the exports of the sources.
- SNP distances have no Ridom equivalent (RIDOM.md, "Not available from Ridom"); do not reintroduce them here.
- Type check: `~/data/conda/envs/env_typecheck/bin/mypy` from the repo root; tests: `pytest` in `env_ridom`. Run real jobs with `sbatch submit_all.sbatch`. Commit directly on `master`; push only when asked.
