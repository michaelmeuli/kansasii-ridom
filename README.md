# kansasii-ridom

Reads the exports of **Ridom SeqSphere+** (the source of truth for cgMLST, MST and clusters of the *M. kansasii* complex, decision
2026-10-30) and turns them into the tables and figures of the final report. It replaces `kansasii-mst` and `kansasii-cgmlst`, which are obsolete.
Why, which exports are needed and what Ridom cannot replace: `kansasii-docs/docs/RIDOM.md`.

Ridom is run by hand in the GUI. Its exports go **unedited** to `data/imm/ridom-<content>-mkan329[-<subset>].<ext>`; this repo only reads them.
The input to Ridom is `output/kansasii-ridom/ridom_assemblies.{tar.gz,zip}` (`scripts/01_export_ridom.py`): the same export as the earlier
`kansasii-cgmlst/scripts/10_export_ridom.py` plus M055 (182 genomes); plasmids are removed by Ridom, not here.

```bash
sbatch submit_all.sbatch      # export check (scripts/00) -> output/kansasii-ridom/results/
sbatch submit_export.sbatch   # archives for Ridom (scripts/01, ~20 min) -> output/kansasii-ridom/
```

| Script | Does |
|---|---|
| `00_check_exports.py` | compares `data/imm/ridom-*` with `config/expected_exports.tsv` (the 12 exports of RIDOM.md) and the isolates each must hold -> `export_status.tsv` (MISSING / INCOMPLETE / EXTRA / OK / UNCHECKED); `--strict` exits 1 unless all are fine |
| `01_export_ridom.py` | `ridom_assemblies.{tar.gz,zip}` for SeqSphere+ (`ridom/assemblies/<name>.fasta`, `ridom/metadata.tsv`): the pool of 175 complex isolates incl. M055 (names `M001` ..., controls by `LABEL`) + the 7 GTDB references from `GTDB_REPRESENTATIVES` (`<acc>_query.fasta`, unchanged). Written to `.part` and renamed when complete |

Next steps (not written yet, they need the exports): distance matrices and MST node tables per subset, cluster table joined to
`screening_map_results.csv`, report figures. Add each as the next numbered script when its export format has been seen, and fill
`id_column` for it in `config/expected_exports.tsv`.

## Layout

```
scripts/            00_check_exports.py, 01_export_ridom.py, common.py (cohort, subsets, ridom_id, REFERENCES); kansasii_paths.py = shared constants (vendored)
config/             expected_exports.tsv: the exports still needed and the cohort each must cover
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
