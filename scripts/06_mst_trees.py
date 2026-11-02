"""MST nodes / edges and UPGMA trees per subset, computed from the Ridom allele distances (exports 3 and 10 of RIDOM.md).

SeqSphere+ draws the MST in the GUI; its distances are those of the .meg export, so the graph data is recomputed here. A node is a set of
isolates at distance 0 (identical on the targets they share). The MST has the same total length as Ridom's; among equally short trees the
choice of edges can differ (ties go to the lowest index). Writes to output/kansasii-ridom/results/:
  mst_nodes.tsv   subset, node, members, n, species, ridom cluster
  mst_edges.tsv   subset, node_a, node_b, distance
  tree_upgma_<subset>.nwk   UPGMA tree (branch length = half the distance), subsets of SUBSETS with at least 2 isolates
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import SUBSETS, load_cohort, read_clusters, read_meg, species_of
from graph import components, mst_edges, upgma_newick
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS


def node_graph(d: np.ndarray, ids: list[str]) -> tuple[list[list[str]], list[tuple[int, int, float]]]:
    """Nodes (isolates at distance 0, lowest index first) and the MST between them (distance = smallest member distance)."""
    lab = components(d, 0)
    order = sorted(set(lab))
    nodes = [[ids[i] for i in range(len(ids)) if lab[i] == k] for k in order]
    n = len(order)
    dn = np.zeros((n, n))
    for x in range(n):
        for y in range(n):
            if x != y:
                dn[x, y] = d[np.ix_([i for i in range(len(ids)) if lab[i] == order[x]],
                                    [i for i in range(len(ids)) if lab[i] == order[y]])].min()
    return nodes, mst_edges(dn)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()

    dist = read_meg(a.imm_dir / "ridom-distance-matrix.meg")
    species = species_of(load_cohort(a.screening_map_results))
    clusters = read_clusters(a.imm_dir / "ridom-clusters.csv")
    a.out.mkdir(parents=True, exist_ok=True)
    nrows, erows = [], []
    for name, wanted in SUBSETS.items():
        ids = [i for i in dist.index if wanted is None or species.get(i) in wanted]
        d = dist.loc[ids, ids].to_numpy(dtype=float)
        nodes, edges = node_graph(d, ids)
        for k, m in enumerate(nodes):
            nrows.append({"subset": name, "node": f"{name}-N{k + 1}", "members": ";".join(m), "n": len(m),
                          "species": ";".join(sorted({species[i] for i in m})),
                          "cluster": ";".join(sorted({clusters.get(i, "") for i in m} - {""}))})
        erows += [{"subset": name, "node_a": f"{name}-N{x + 1}", "node_b": f"{name}-N{y + 1}", "distance": int(w)} for x, y, w in edges]
        if len(ids) >= 2:
            (a.out / f"tree_upgma_{name}.nwk").write_text(upgma_newick(d, ids) + "\n")
        print(f"{name}: {len(ids)} isolates, {len(nodes)} nodes, MST length {int(sum(w for _, _, w in edges))}")
    pd.DataFrame(nrows).to_csv(a.out / "mst_nodes.tsv", sep="\t", index=False)
    pd.DataFrame(erows).to_csv(a.out / "mst_edges.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
