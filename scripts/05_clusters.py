"""Cluster table of all isolates incl. singletons, and the distance threshold behind Ridom's MST clusters (export 5 of RIDOM.md).

Reads data/imm/ridom-clusters.csv (`-Group-`, empty = not in a cluster) and the distance matrix of export 2. Ridom's clusters are
single-linkage groups at an allele-distance threshold; the script finds the threshold at which they are reproduced exactly
(`cluster_threshold.tsv`) and uses it to explain every isolate: cluster, size, nearest neighbour. Writes to output/kansasii-ridom/results/:
  cluster_threshold.tsv  per threshold: clusters, clustered isolates, identical to Ridom's grouping (true / false)
  clusters.tsv           every isolate: kind, species, Ridom cluster ("singleton" if none), cluster size, nearest neighbour and distance
  cluster_summary.tsv    per cluster: size, species, diameter (largest pairwise distance), largest nearest-neighbour step, distance to outside
Exit 1 if no threshold reproduces Ridom's clusters.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import kind_of, load_cohort, read_clusters, read_meg, species_of
from graph import components
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS

MAX_THRESHOLD = 60


def groups_at(d: np.ndarray, ids: list[str], t: float) -> set[frozenset[str]]:
    lab = components(d, t)
    g: dict[int, set[str]] = {}
    for i, k in enumerate(lab):
        g.setdefault(k, set()).add(ids[i])
    return {frozenset(v) for v in g.values() if len(v) > 1}


def ridom_groups(clusters: pd.Series[str], ids: list[str]) -> set[frozenset[str]]:
    g: dict[str, set[str]] = {}
    for s in ids:
        c = clusters.get(s, "")
        if c:
            g.setdefault(c, set()).add(s)
    return {frozenset(v) for v in g.values()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()

    dist = read_meg(a.imm_dir / "ridom-distance-matrix.meg")
    clusters = read_clusters(a.imm_dir / "ridom-clusters.csv")
    cohort = load_cohort(a.screening_map_results)
    species = species_of(cohort)
    ids = list(dist.index)
    d = dist.to_numpy(dtype=float)
    missing = [i for i in ids if i not in clusters.index]
    if missing:
        print(f"ERROR: {len(missing)} isolates of the matrix are not in ridom-clusters.csv: {', '.join(missing[:8])}", file=sys.stderr)
        return 1
    target = ridom_groups(clusters, ids)

    rows = []
    for t in range(MAX_THRESHOLD + 1):
        g = groups_at(d, ids, t)
        rows.append({"threshold": t, "clusters": len(g), "clustered_isolates": sum(map(len, g)), "same_as_ridom": g == target})
    th = pd.DataFrame(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    th.to_csv(a.out / "cluster_threshold.tsv", sep="\t", index=False)
    hit = th[th.same_as_ridom]
    if hit.empty:
        print("ERROR: no single-linkage threshold reproduces Ridom's clusters (see cluster_threshold.tsv)", file=sys.stderr)
        return 1
    print(f"Ridom clusters = single linkage at threshold {int(hit.threshold.min())}..{int(hit.threshold.max())} alleles")

    off = d + np.diag(np.full(len(d), np.inf))
    nn = off.argmin(axis=1)
    label = {i: clusters[i] or "singleton" for i in ids}
    size = {c: sum(1 for i in ids if label[i] == c) for c in set(label.values())}
    out = pd.DataFrame({
        "sample": ids, "kind": [kind_of(i, cohort) for i in ids], "species": [species[i] for i in ids],
        "cluster": [label[i] for i in ids],
        "cluster_size": [size[label[i]] if label[i] != "singleton" else 1 for i in ids],
        "nearest_neighbour": [ids[j] for j in nn], "nn_distance": [int(d[k, j]) for k, j in enumerate(nn)],
    })
    out.to_csv(a.out / "clusters.tsv", sep="\t", index=False)

    srows = []
    pos = {i: k for k, i in enumerate(ids)}
    for c in sorted((c for c in size if c != "singleton"), key=lambda x: int(x.rsplit(" ", 1)[1])):
        m = [pos[i] for i in ids if label[i] == c]
        out_idx = [k for k in range(len(ids)) if k not in m]
        sub = d[np.ix_(m, m)]
        step = max(float(off[np.ix_(m, m)][k].min()) for k in range(len(m)))
        srows.append({"cluster": c, "n": len(m), "species": ";".join(sorted({species[ids[k]] for k in m})),
                      "members": ";".join(ids[k] for k in m), "diameter": int(sub.max()), "max_nn_step": int(step),
                      "min_distance_outside": int(d[np.ix_(m, out_idx)].min()) if out_idx else ""})
    pd.DataFrame(srows).to_csv(a.out / "cluster_summary.tsv", sep="\t", index=False)
    print(f"{len(srows)} clusters, {sum(int(str(r['n'])) for r in srows)} of {len(ids)} isolates clustered, {len(ids) - sum(int(str(r['n'])) for r in srows)} singletons")
    return 0


if __name__ == "__main__":
    sys.exit(main())
