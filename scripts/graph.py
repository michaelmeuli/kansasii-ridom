"""Small graph algorithms on allele-distance matrices (numpy only): single-linkage clusters, MST, UPGMA tree."""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

Matrix = npt.NDArray[np.float64]


def components(d: Matrix, threshold: float) -> list[int]:
    """Single-linkage clusters: isolates joined by a chain of pairs with distance <= threshold. Label = first member's index."""
    n = len(d)
    label = list(range(n))
    for i in range(n):
        if label[i] != i:
            continue
        stack = [i]
        while stack:
            u = stack.pop()
            for v in np.flatnonzero(d[u] <= threshold):
                if label[int(v)] == int(v) and int(v) != i:
                    label[int(v)] = i
                    stack.append(int(v))
    return label


def mst_edges(d: Matrix) -> list[tuple[int, int, float]]:
    """Minimum spanning tree (Prim, start at 0). Ties go to the lowest index, so the result is deterministic.

    Equal-distance alternatives exist in allele data; SeqSphere+ may draw a different but equally short tree.
    """
    n = len(d)
    if n < 2:
        return []
    in_tree = np.zeros(n, dtype=bool)
    in_tree[0] = True
    best = d[0].copy()
    parent = np.zeros(n, dtype=int)
    edges: list[tuple[int, int, float]] = []
    for _ in range(n - 1):
        cand = np.where(in_tree, np.inf, best)
        v = int(np.argmin(cand))
        edges.append((int(parent[v]), v, float(cand[v])))
        in_tree[v] = True
        closer = (d[v] < best) & ~in_tree
        best[closer] = d[v][closer]
        parent[closer] = v
    return edges


def upgma_newick(d: Matrix, labels: list[str]) -> str:
    """UPGMA (average linkage) tree as Newick; branch length = half the distance between the merged clusters."""
    n = len(labels)
    if n == 0:
        return ";"
    if n == 1:
        return f"{labels[0]};"
    work = d.astype(float).copy()
    np.fill_diagonal(work, np.inf)
    nodes: list[str | None] = list(labels)
    height = [0.0] * n
    size = [1] * n
    active = list(range(n))
    while len(active) > 1:
        sub = work[np.ix_(active, active)]
        k = int(np.argmin(sub))
        a, b = active[k // len(active)], active[k % len(active)]
        if a > b:
            a, b = b, a
        h = work[a, b] / 2
        nodes[a] = f"({nodes[a]}:{h - height[a]:g},{nodes[b]}:{h - height[b]:g})"
        for c in active:
            if c not in (a, b):
                work[a, c] = work[c, a] = (work[a, c] * size[a] + work[b, c] * size[b]) / (size[a] + size[b])
        size[a] += size[b]
        height[a] = h
        active.remove(b)
        nodes[b] = None
    return f"{nodes[active[0]]};"
