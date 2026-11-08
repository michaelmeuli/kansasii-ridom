"""Graph algorithms, task-template parser and the derived tables on small hand-made data."""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from common import read_clusters, read_comparison  # noqa: E402
from graph import components, mst_edges, upgma_newick  # noqa: E402
from scheme import parse_filters, parse_funnel, parse_genomes, read_tasktemplate  # noqa: E402


def _load(name: str, file: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    assert spec is not None and spec.loader is not None
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


qc_mod = _load("qc", "04_qc.py")
ctl_mod = _load("ctl", "07_control_distances.py")
mst_mod = _load("mst", "06_mst_trees.py")
paper = _load("paper", "08_paper_numbers.py")

D = np.array([[0, 1, 5, 9], [1, 0, 6, 9], [5, 6, 0, 2], [9, 9, 2, 0]], dtype=float)


def test_components() -> None:
    assert components(D, 1) == [0, 0, 2, 3]
    assert components(D, 2) == [0, 0, 2, 2]
    assert len(set(components(D, 5))) == 1  # 1-0-2 chain at 5 joins 3 via 2-3 = 2
    assert components(D, 0) == [0, 1, 2, 3]


def test_mst_and_nodes() -> None:
    e = mst_edges(D)
    assert sum(w for _, _, w in e) == 1 + 5 + 2 and len(e) == 3
    assert mst_edges(D[:1, :1]) == []
    d0 = np.array([[0, 0, 3], [0, 0, 4], [3, 4, 0]], dtype=float)
    nodes, edges = mst_mod.node_graph(d0, ["a", "b", "c"])
    assert nodes == [["a", "b"], ["c"]] and edges == [(0, 1, 3.0)]


def test_upgma() -> None:
    assert upgma_newick(D, ["a", "b", "c", "d"]).count("(") == 3
    t = upgma_newick(np.array([[0, 4], [4, 0]], dtype=float), ["x", "y"])
    assert t == "(x:2,y:2);"
    assert upgma_newick(np.zeros((1, 1)), ["x"]) == "x;"


REPORT = """cgMLST Target Definer Results

Seed Genome Sources:
 * GCF_1.1 (01-JAN-2026), 1000 bases, 5 genes with CDS (Mycobacterium kansasii ATCC 12478, complete sequence)

Penetration Query Genomes (1):
 * GCF_2.1, 900 bases, 3 contigs (Mycobacterium persicum X)

Exclude Sequences (1):
 * NC_1.1 (01-JAN-2026), 100 bases, 2 genes with CDS (plasmid)

Seed Genome Filters:
 * Minimum Length Filter (requires >=50 bases)

DETAILS

Progress:
 * 5 targets found in seed genome GCF_1.1
 * 3 targets after blasting against GCF_2.1

TARGET DETAILS

Seed Genome Start Codon Filter: filtered out 2 targets (discarded)
MKAN_RS00001, MKAN_RS00002

Query Genome BLAST Search: filtered out 1 targets (moved to Accessory)
MKAN_RS00003
Messages for GCF_1.1:
 * Skipped MKAN_RS00009 because location is not contiguous.
"""


def _template(tmp: Path) -> Path:
    def target(tag: str, seq: str) -> str:
        return (f'<PCRBankEntry type="template"><Name>{tag}</Name><Comment><Properties><Key>GenBank locus_tag</Key><Value>{tag}</Value>'
                f'<Key>GenBank product</Key><Value>prod {tag}</Value></Properties></Comment>'
                f'<Ridom><ReferenceSequence><Name>Ref {tag}</Name><Bases>{seq}</Bases></ReferenceSequence><CodonTable>11</CodonTable>'
                f'<RefAreas><RefArea><Forward>1</Forward><ContainsStart>1</ContainsStart><ContainsStop>1</ContainsStop></RefArea></RefAreas>'
                f'</Ridom></PCRBankEntry>')
    xml = (f'<?xml version="1.0"?><TaskTemplate><PCRBankEntry type="group" created="2026-01-01" creator="rt"><Name>scheme</Name>'
           f'<NgsParams><scan_refseq_identity>90.0</scan_refseq_identity></NgsParams>{target("MKAN_RS00004", "ATGAAATGA")}'
           f'{target("MKAN_RS00005", "ATGTAA")}</PCRBankEntry><AdditionalValues><map><item key="Bundle" value="b 2/1"/></map></AdditionalValues>'
           f'<Comment>{REPORT}</Comment></TaskTemplate>')
    p = tmp / "s.tasktemplate"
    with gzip.open(p, "wb") as fh:
        fh.write(xml.encode())
    return p


def test_tasktemplate(tmp_path: Path) -> None:
    s = read_tasktemplate(_template(tmp_path))
    assert s.name == "scheme" and s.created == "2026-01-01" and s.additional["Bundle"] == "b 2/1"
    assert [t["locus_tag"] for t in s.targets] == ["MKAN_RS00004", "MKAN_RS00005"]
    assert s.targets[0]["length_bp"] == "9" and s.targets[0]["strand"] == "+" and s.targets[0]["product"] == "prod MKAN_RS00004"
    f = parse_filters(s.report)
    assert [(x["fate"], x["n_reported"], x["loci"]) for x in f] == [
        ("discarded", 2, ["MKAN_RS00001", "MKAN_RS00002"]), ("moved to Accessory", 1, ["MKAN_RS00003"])]  # Messages block not a locus
    assert parse_funnel(s.report) == [("found in seed genome GCF_1.1", 5), ("after blasting against GCF_2.1", 3)]
    g = parse_genomes(s.report)
    assert [(x["role"], x["accession"], x["bases"]) for x in g] == [
        ("seed", "GCF_1.1 (01-JAN-2026)", "1000"), ("query", "GCF_2.1", "900"), ("exclude", "NC_1.1 (01-JAN-2026)", "100")]


def test_tasktemplate_without_report(tmp_path: Path) -> None:
    p = tmp_path / "x.tasktemplate"
    with gzip.open(p, "wb") as fh:
        fh.write(b'<TaskTemplate><PCRBankEntry type="group"><Name>s</Name></PCRBankEntry></TaskTemplate>')
    with pytest.raises(ValueError, match="Target Definer"):
        read_tasktemplate(p)


def _table(tmp: Path) -> Path:
    p = tmp / "c.csv"
    p.write_text('"#Missing values in Distance Columns";"Perc. Good Targets";"Approximated Genome Size (Mbases)";"Sample ID";"Species";"L1";"L2";"L3";"L4";\n'
                 '0;100.0;6.4;GCF_000157895.3;M. kansasii;1;1;2;3;\n'
                 '1;75.0;6.2;M001;M. kansasii;1;;2;3;\n')
    return p


def test_comparison_and_qc(tmp_path: Path) -> None:
    meta, al = read_comparison(_table(tmp_path))
    assert list(al.columns) == ["L1", "L2", "L3", "L4"] and al.loc["M001", "L2"] == ""
    assert meta.loc["M001", "missing_ridom"] == 1
    cohort = pd.DataFrame({"sample": ["Mkan329-001"], "ridom_id": ["M001"], "label": [""], "species": ["kansasii"], "keep": [True]})
    qc = qc_mod.qc_table(meta, al, cohort)
    assert qc.consistent.all() and list(qc.kind) == ["reference", "isolate"] and list(qc.cohort) == ["reference", "keep"]
    ls = qc_mod.locus_stats(al)
    assert list(ls.n_missing) == [0, 1, 0, 0] and list(ls.n_alleles) == [1, 1, 1, 1]
    meta.loc["M001", "good_pct"] = 80.0
    assert not qc_mod.qc_table(meta, al, cohort).consistent.all()
    meta.loc["M001", "good_pct"] = 75.04  # within Ridom's one-decimal rounding
    assert qc_mod.qc_table(meta, al, cohort).consistent.all()


def test_clusters_file(tmp_path: Path) -> None:
    p = tmp_path / "cl.csv"
    p.write_text("Perc. Good Targets,Approximated Genome Size (Mbases),Sample ID,Species,-Group-\n96,6,M001,x,MST Cluster 1\n96,6,M002,x,\n")
    c = read_clusters(p)
    assert c["M001"] == "MST Cluster 1" and c["M002"] == ""


def test_control_rows() -> None:
    ids = ["GCF_000157895.3", "M_kansasii_ATCC_12478", "M001", "M002"]
    dist = pd.DataFrame(np.array([[0, 2, 7, 8], [2, 0, 4, 5], [7, 4, 0, 1], [8, 5, 1, 0]]), index=ids, columns=ids)
    cohort = pd.DataFrame({"sample": ["a", "b", "c"], "ridom_id": ids[1:], "label": ["ATCC 12478", "", ""],
                           "species": ["kansasii"] * 3, "keep": [True] * 3})
    r = ctl_mod.control_rows(dist, cohort, nearest=1)
    assert r[["comparator", "comparator_kind", "distance"]].values.tolist() == [
        ["GCF_000157895.3", "reference", 2], ["M001", "nearest", 4]]


def test_species_table() -> None:
    qc = pd.DataFrame({"sample": ["M001", "M002", "M055", "GCF_000157895.3"], "kind": ["isolate", "isolate", "isolate", "reference"],
                       "species": ["kansasii"] * 4, "cohort": ["keep", "keep", "pool_only", "reference"],
                       "good_targets_pct": [99.0, 98.0, 60.0, 100.0], "genome_mbases": [6.4, 6.3, 6.0, 6.4]})
    ds = pd.DataFrame({"subset": ["kansasii"], "median": [9.0], "min": [0], "max": [20], "identical_pairs": [3],
                       "study_median": [3.0], "study_min": [1], "study_max": [5], "study_identical_pairs": [0]})
    cl = pd.DataFrame({"sample": ["M001", "M002", "GCF_000157895.3"], "kind": ["isolate", "isolate", "reference"], "species": ["kansasii"] * 3,
                       "cluster": ["MST Cluster 1", "singleton", "MST Cluster 1"]})
    t = paper.species_table(qc, ds, cl)
    k = t[t.species == "kansasii"].iloc[0]
    assert k.isolates == 2 and k.references == 1 and k.clusters == 1 and k.singletons == 1  # M055 not counted
    assert k.good_targets_pct == "98.5 (98.0-99.0)"
    assert (k.distance_median, k.distance_min, k.distance_max, k.identical_pairs) == (3.0, 1, 5, 0)  # without the reference and the control strain
