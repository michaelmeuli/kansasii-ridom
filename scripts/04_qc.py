"""Per-isolate QC and per-locus statistics from the Ridom comparison table (exports 7 and 8 of RIDOM.md).

Reads data/imm/ridom-comparison-table-mkan329.csv (pool incl. M055 + the 7 references). Ridom writes no allele for a target
it could not call and does not export the reason, so a failure here is "no allele". Writes to output/kansasii-ridom/results/:
  qc.tsv               per isolate: kind, species, good targets %, genome size, missing / called targets (recomputed and Ridom's)
  qc_summary.tsv       per species and for all isolates: n, good targets and genome size (median, min, max)
  locus_stats.tsv      per target: isolates without allele, % missing, number of distinct alleles
  target_failures.tsv  isolate x target without allele (long form)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import kind_of, load_cohort, read_comparison
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS


def qc_table(meta: pd.DataFrame, alleles: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    n_loci = alleles.shape[1]
    missing = (alleles == "").sum(axis=1)
    sp = dict(zip(cohort.ridom_id, cohort.species))
    keep = dict(zip(cohort.ridom_id, cohort.keep))
    out = pd.DataFrame({
        "sample": meta.index,
        "kind": [kind_of(s, cohort) for s in meta.index],
        "species": [sp.get(s, meta.species_ridom[s].replace("M. ", "")) for s in meta.index],
        "cohort": ["reference" if s.startswith("GCF_") else ("keep" if keep.get(s) else "pool_only") for s in meta.index],
        "good_targets_pct": meta.good_pct.to_numpy(),
        "genome_mbases": meta.genome_mbases.to_numpy(),
        "targets_total": n_loci,
        "targets_missing": missing.to_numpy(),
        "targets_called": (n_loci - missing).to_numpy(),
        "missing_ridom": meta.missing_ridom.to_numpy(),
    })
    out["consistent"] = (out.targets_missing == out.missing_ridom) & (
        (100 * out.targets_called / n_loci).round(1) == out.good_targets_pct)
    return out


def summary(qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [("all", qc), *[(str(k), g) for k, g in qc.groupby("species")]]
    for name, g in groups:
        rows.append({"species": name, "n": len(g), "n_references": int((g.kind == "reference").sum()),
                     "good_median": g.good_targets_pct.median(), "good_min": g.good_targets_pct.min(), "good_max": g.good_targets_pct.max(),
                     "size_median_mb": g.genome_mbases.median(), "size_min_mb": g.genome_mbases.min(), "size_max_mb": g.genome_mbases.max()})
    return pd.DataFrame(rows)


def locus_stats(alleles: pd.DataFrame) -> pd.DataFrame:
    called = alleles != ""
    return pd.DataFrame({
        "locus_tag": alleles.columns,
        "isolates": len(alleles),
        "n_missing": (~called).sum().to_numpy(),
        "missing_pct": (100 * (~called).mean()).round(2).to_numpy(),
        "n_alleles": [alleles[c][called[c]].nunique() for c in alleles.columns],
    })


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--comparison-table", type=Path, help="default: <imm-dir>/ridom-comparison-table-mkan329.csv")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()

    meta, alleles = read_comparison(a.comparison_table or a.imm_dir / "ridom-comparison-table-mkan329.csv")
    qc = qc_table(meta, alleles, load_cohort(a.screening_map_results))
    a.out.mkdir(parents=True, exist_ok=True)
    qc.to_csv(a.out / "qc.tsv", sep="\t", index=False)
    s = summary(qc)
    s.to_csv(a.out / "qc_summary.tsv", sep="\t", index=False)
    locus_stats(alleles).to_csv(a.out / "locus_stats.tsv", sep="\t", index=False)
    rows, cols = np.nonzero((alleles == "").to_numpy())
    pd.DataFrame({"sample": alleles.index[rows], "locus_tag": alleles.columns[cols]}).to_csv(
        a.out / "target_failures.tsv", sep="\t", index=False)
    print(s.to_string(index=False))
    bad = qc[~qc.consistent]
    if len(bad):
        print(f"ERROR: Ridom's missing count / good-target % differ from the table for: {', '.join(bad['sample'])}", file=sys.stderr)
    return 1 if len(bad) else 0


if __name__ == "__main__":
    sys.exit(main())
