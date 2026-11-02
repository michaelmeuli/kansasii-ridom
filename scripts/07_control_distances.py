"""Distances of the control strains to the public references and to their nearest isolates (export 9 of RIDOM.md).

Controls are the cohort strains with a LABEL (the type / reference strains NR 128, 132, 133); references are the 7 GTDB genomes (GCF_*).
Reads the distance matrix of export 2. Writes output/kansasii-ridom/results/control_distances.tsv: per control the references of the
same species (comparator_kind `reference`) and its NEAREST nearest others (`nearest`), with the allele distance.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from common import REFERENCES, kind_of, load_cohort, read_meg, species_of
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS

NEAREST = 5


def control_rows(dist: pd.DataFrame, cohort: pd.DataFrame, nearest: int = NEAREST) -> pd.DataFrame:
    species = species_of(cohort)
    rows = []
    for c in (i for i in dist.index if kind_of(i, cohort) == "control"):
        refs = [r for r in dist.index if r in REFERENCES and REFERENCES[r] == species[c]]
        for r in refs:
            rows.append({"control": c, "species": species[c], "comparator": r, "comparator_kind": "reference", "distance": int(dist.loc[c, r])})
        others = dist.loc[c].drop(c).drop(refs).sort_values(kind="stable").head(nearest)
        rows += [{"control": c, "species": species[c], "comparator": o, "comparator_kind": "nearest", "distance": int(v)}
                 for o, v in others.items()]
    return pd.DataFrame(rows, columns=["control", "species", "comparator", "comparator_kind", "distance"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()
    out = control_rows(read_meg(a.imm_dir / "ridom-distance-complex-mkan329.meg"), load_cohort(a.screening_map_results))
    a.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out / "control_distances.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    return 0 if len(out) else 1


if __name__ == "__main__":
    sys.exit(main())
