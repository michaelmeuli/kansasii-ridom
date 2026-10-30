# kansasii-ridom

Reads the exports of **Ridom SeqSphere+** (the source of truth for cgMLST, MST and clusters of the *M. kansasii* complex, decision
2026-10-30) and turns them into the tables and figures of the final report. It replaces `kansasii-mst` and `kansasii-cgmlst`, which are obsolete.
Why, which exports are needed and what Ridom cannot replace: `kansasii-docs/docs/RIDOM.md`.

Ridom is run by hand in the GUI. Its exports go **unedited** to `data/imm/ridom-<content>-mkan329[-<subset>].<ext>`; this repo only reads them.
The input to Ridom is still `output/kansasii-cgmlst/ridom_assemblies.{tar.gz,zip}` (from `kansasii-cgmlst/scripts/10_export_ridom.py`).

```bash
sbatch submit_all.sbatch   # scripts/00 .. in order (outputs: output/kansasii-ridom/results/)
```

| Script | Does |
|---|---|
| `00_check_exports.py` | compares `data/imm/ridom-*` with `config/expected_exports.tsv` (the 12 exports of RIDOM.md) and the isolates each must hold -> `export_status.tsv` (MISSING / INCOMPLETE / EXTRA / OK / UNCHECKED); `--strict` exits 1 unless all are fine |

Next steps (not written yet, they need the exports): distance matrices and MST node tables per subset, cluster table joined to
`screening_map_results.csv`, report figures. Add each as the next numbered script when its export format has been seen, and fill
`id_column` for it in `config/expected_exports.tsv`.

## Layout

```
scripts/            00_check_exports.py, common.py (cohort, subsets, ridom_id); kansasii_paths.py = shared constants (vendored)
config/             expected_exports.tsv: the exports still needed and the cohort each must cover
tests/              pytest
submit_all.sbatch   runs scripts/00.. in order (log in output/kansasii-ridom/logs/)
```

**Environment:** `environment.yml` (`env_ridom`: pandas, pytest). Type check: `~/data/conda/envs/env_typecheck/bin/mypy`.

**Cohort:** the isolates sent to Ridom are the complex species without `excluded` in `screening_map_results.csv` (174: the 3 complex control
strains NR 128/132/133 are kept, named by `LABEL`) plus the 7 public GTDB references (`GCF_*`) = 181. Ridom ids: `Mkan329-081` is `M081`
(`ridom_id()` in `common.py`, same rule as in `immensekansasii/scripts/06`). Some exports must also cover M055 (the `pool`, 175), which is
excluded for low chewBBACA coverage only; Ridom needs to give its own good-target value.

Paths come from `scripts/kansasii_paths.py` (same names in every repo) and are flags of the scripts (`--imm-dir`, `--screening-map-results`,
`--out`). Edit `kansasii-docs/templates/kansasii_paths.py`, not the copy (`bash repos/check_sync.sh --fix`).

## Provenance

The sbatch file writes `MANIFEST.txt` (commit, clean/dirty flag, script, conda env hash, inputs) into the output directory when it
finishes; see `kansasii-docs/docs/PROVENANCE.md`. The GUI exports have no manifest: record each in the export log of `RIDOM.md`.
Project index: `kansasii-docs/README.md`.
