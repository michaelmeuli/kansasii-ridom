"""Allele-distance matrices per subset from the complex-wide Ridom export (export 2 of RIDOM.md).

Reads data/imm/ridom-distance-complex-mkan329.meg (SeqSphere+ MEGA, lower triangle). Missing values were handled by
Ridom when it computed the matrix (pairwise ignore, see the export log of RIDOM.md); under pairwise ignore the distance of
a pair does not depend on the other isolates, so a subset matrix is the submatrix of the complex one.

Writes to output/kansasii-ridom/results/: distance_<subset>.tsv (square, Sample ID in the first column, subsets of
DISTANCE_SUBSETS in common.py: the complex, the MST subsets incl. the combined `minor`, and each minor species; references are in the subset of their species) and distance_summary.tsv (n, pairs, min / median /
max distance and number of identical pairs per subset).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import DISTANCE_SUBSETS, load_cohort, read_meg, species_of
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS


def subset_matrix(dist: pd.DataFrame, species: dict[str, str], wanted: list[str] | None) -> pd.DataFrame:
    """Rows and columns of `dist` whose isolate is of one of the `wanted` species (None = all)."""
    if wanted is None:
        return dist
    keep = [i for i in dist.index if species.get(i) in wanted]
    return dist.loc[keep, keep]


def summarize(name: str, d: pd.DataFrame) -> dict[str, object]:
    pairs = d.to_numpy()[np.triu_indices(len(d), 1)]
    return {"subset": name, "n": len(d), "pairs": len(pairs),
            "min": int(pairs.min()) if len(pairs) else "", "median": float(np.median(pairs)) if len(pairs) else "",
            "max": int(pairs.max()) if len(pairs) else "", "identical_pairs": int((pairs == 0).sum())}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--meg", type=Path, help="default: <imm-dir>/ridom-distance-complex-mkan329.meg")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()

    dist = read_meg(a.meg or a.imm_dir / "ridom-distance-complex-mkan329.meg")
    cohort = load_cohort(a.screening_map_results)
    species = species_of(cohort)
    unknown = [i for i in dist.index if i not in species]
    if unknown:
        print(f"ERROR: {len(unknown)} matrix labels not in the cohort or references: {', '.join(unknown[:8])}", file=sys.stderr)
        return 1
    a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, wanted in DISTANCE_SUBSETS.items():
        d = subset_matrix(dist, species, wanted)
        d.rename_axis("Sample ID").to_csv(a.out / f"distance_{name}.tsv", sep="\t")
        rows.append(summarize(name, d))
    summary = pd.DataFrame(rows)
    summary.to_csv(a.out / "distance_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
