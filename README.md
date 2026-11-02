# kansasii-ridom

Reads the exports of **Ridom SeqSphere+** (the source of truth for cgMLST, MST and clusters of the *M. kansasii* complex, decision
2026-10-30) and turns them into the tables and figures of the final report. It replaces `kansasii-mst` and `kansasii-cgmlst`, which are obsolete.
Why, which exports are needed and what Ridom cannot replace: `kansasii-docs/docs/RIDOM.md`.

Ridom is run by hand in the GUI. Its exports go **unedited** to `data/imm/ridom-<content>-mkan329[-<subset>].<ext>`; this repo only reads them.
The input to Ridom is `output/kansasii-ridom/ridom_assemblies.{tar.gz,zip}` (`scripts/01_export_ridom.py`): the same export as the earlier
`kansasii-cgmlst/scripts/10_export_ridom.py` plus M055 (182 genomes), except that the ATCC 12478 reference GCF_000157895.3 keeps its plasmid contig (plasmids are removed by Ridom, not here).

```bash
sbatch submit_all.sbatch      # export check (00) + all derived tables (02-08) -> output/kansasii-ridom/results/
sbatch submit_export.sbatch   # archives for Ridom (scripts/01, ~20 min) -> output/kansasii-ridom/
```

| Script | Does |
|---|---|
| `00_check_exports.py` | compares `data/imm/ridom-*` with `config/expected_exports.tsv` (the 12 exports of RIDOM.md) and the isolates each must hold -> `export_status.tsv` (MISSING / INCOMPLETE / EXTRA / OK / UNCHECKED); `--strict` exits 1 unless all are fine |
| `02_distance_matrices.py` | `ridom-distance-complex-mkan329.meg` (SeqSphere+ MEGA, pairwise ignore) -> `distance_<subset>.tsv` for each subset of `DISTANCE_SUBSETS` (`SUBSETS` + one per minor species; submatrices of the complex one) + `distance_summary.tsv` |
| `03_scheme.py` | `ridom-scheme-mkan329.tasktemplate` (gzip XML) -> core targets, core / accessory / discarded classes, Target Definer filters, funnel, genomes, settings, `scheme_summary.tsv`, `scheme_parameters.txt`; checks the comparison-table loci = core targets |
| `04_qc.py` | comparison table -> `qc.tsv` (good targets %, genome size, missing targets; consistency with Ridom's own columns), `qc_summary.tsv`, `locus_stats.tsv`, `target_failures.tsv` |
| `05_clusters.py` | `ridom-clusters.csv` + distances -> `clusters.tsv` (all isolates incl. singletons, nearest neighbour), `cluster_summary.tsv`, `cluster_threshold.tsv` (Ridom's clusters = single linkage at <= 5 alleles; exit 1 if no threshold matches) |
| `06_mst_trees.py` | distances -> `mst_nodes.tsv`, `mst_edges.tsv` (node = isolates at distance 0), `tree_upgma_<subset>.nwk` per subset |
| `07_control_distances.py` | controls (NR 128/132/133) vs. same-species references and nearest isolates -> `control_distances.tsv` |
| `08_paper_numbers.py` | `paper_numbers.md` (methods / results figures with their source file) and `paper_species_table.tsv` |
| `01_export_ridom.py` | `ridom_assemblies.{tar.gz,zip}` for SeqSphere+ (`ridom/assemblies/<name>.fasta`, `ridom/metadata.tsv`): the pool of 175 complex isolates incl. M055 (names `M001` ..., controls by `LABEL`) + the 7 GTDB references from `GTDB_REPRESENTATIVES` (`<acc>_query.fasta`, unchanged). Written to `.part` and renamed when complete |

Only the MST figures (SVG / PDF per subset) are still a pure GUI export; everything else the paper needs is derived from the four exports in hand (`kansasii-docs/docs/RIDOM.md`). Failure reasons, new-allele flags and target coordinates are not in any export.

## Layout

```
scripts/            00_check_exports.py, 01_export_ridom.py, 02-08 (derived tables), common.py, graph.py (clusters, MST, UPGMA), scheme.py (task template) (cohort, subsets, ridom_id, REFERENCES); kansasii_paths.py = shared constants (vendored)
config/             expected_exports.tsv: the five GUI exports and the cohort each must cover
tests/              pytest
submit_all.sbatch   runs scripts/00 (log in output/kansasii-ridom/logs/); submit_export.sbatch runs scripts/01
```

**Environment:** `environment.yml` (`env_ridom`: pandas, pytest). Type check: `~/data/conda/envs/env_typecheck/bin/mypy`.

**Cohort:** the isolates sent to Ridom are the complex species without `excluded` in `screening_map_results.csv` (174: the 3 complex control
strains NR 128/132/133 are kept, named by `LABEL`) plus the 7 public GTDB references (`GCF_*`) = 181. Ridom ids: `Mkan329-081` is `M081`
(`ridom_id()` in `common.py`, same rule as in `immensekansasii/scripts/06`). The archive for Ridom covers the whole `pool` (175 + 7 = 182), so it includes M055, which is
excluded for low chewBBACA coverage only; Ridom gives its own good-target value for it. The 181 genomes of the earlier export are unchanged.

Paths come from `scripts/kansasii_paths.py` (same names in every repo) and are flags of the scripts (`--imm-dir`, `--screening-map-results`,
`--out`). Edit `kansasii-docs/templates/kansasii_paths.py`, not the copy (`bash repos/check_sync.sh --fix`).

## Provenance

The sbatch files write `MANIFEST.txt` (`ridom_assemblies.MANIFEST.txt` next to the archives) (commit, clean/dirty flag, script, conda env hash, inputs) into the output directory when it
finishes; see `kansasii-docs/docs/PROVENANCE.md`. The GUI exports have no manifest: record each in the export log of `RIDOM.md`.
Project index: `kansasii-docs/README.md`.
