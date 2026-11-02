"""MEGA reader and subset slicing on a hand-made 4-taxon matrix."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from common import meg_label, read_meg, species_of  # noqa: E402

_spec = importlib.util.spec_from_file_location("distance_matrices", SCRIPTS / "02_distance_matrices.py")
assert _spec is not None and _spec.loader is not None
dm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dm)
check_spec = importlib.util.spec_from_file_location("check_exports", SCRIPTS / "00_check_exports.py")
assert check_spec is not None and check_spec.loader is not None
check = importlib.util.module_from_spec(check_spec)
check_spec.loader.exec_module(check)

MEG = """#mega
!TITLE  test;
!Format DataType=distance DataFormat=LowerLeft NTaxa=4;

#M001
#M002
#GCF_900566075_1
#M003

 1
 10 11
 2 12 13
"""


def _write(tmp: Path, text: str = MEG) -> Path:
    p = tmp / "d.meg"
    p.write_text(text)
    return p


def test_meg_label() -> None:
    assert meg_label("GCF_000157895_3") == "GCF_000157895.3"
    assert meg_label("M_kansasii_ATCC_12478") == "M_kansasii_ATCC_12478"


def test_read_meg(tmp_path: Path) -> None:
    d = read_meg(_write(tmp_path))
    assert list(d.index) == ["M001", "M002", "GCF_900566075.1", "M003"]
    assert d.loc["M002", "M001"] == d.loc["M001", "M002"] == 1
    assert d.loc["GCF_900566075.1", "M002"] == 11 and d.loc["M003", "GCF_900566075.1"] == 13
    assert (d.to_numpy().diagonal() == 0).all() and (d.to_numpy() == d.to_numpy().T).all()


def test_read_meg_rejects_bad_files(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="values"):
        read_meg(_write(tmp_path, MEG.replace(" 2 12 13", " 2 12")))
    with pytest.raises(ValueError, match="NTaxa"):
        read_meg(_write(tmp_path, MEG.replace("NTaxa=4", "NTaxa=5")))
    with pytest.raises(ValueError, match="not a MEGA"):
        read_meg(_write(tmp_path, MEG.replace("#mega", "#nexus", 1)))


def test_subset_and_summary(tmp_path: Path) -> None:
    import pandas as pd
    d = read_meg(_write(tmp_path))
    cohort = pd.DataFrame({"ridom_id": ["M001", "M002", "M003"], "species": ["kansasii", "kansasii", "persicum"]})
    sp = species_of(cohort)
    k = dm.subset_matrix(d, sp, ["kansasii"])
    assert list(k.index) == ["M001", "M002"] and k.loc["M001", "M002"] == 1  # same value as in the complex matrix
    assert list(dm.subset_matrix(d, sp, ["pseudokansasii"]).index) == ["GCF_900566075.1"]  # reference by its species
    s = dm.summarize("kansasii", k)
    assert (s["n"], s["pairs"], s["min"], s["max"], s["identical_pairs"]) == (2, 1, 1, 1, 0)
    assert dm.summarize("pseudokansasii", dm.subset_matrix(d, sp, ["pseudokansasii"]))["pairs"] == 0


def test_check_file_meg(tmp_path: Path) -> None:
    p = _write(tmp_path)
    assert check.check_file(p, "labels", {"M001", "M002", "M003"}, True)[0] == "OK"
    assert check.check_file(p, "labels", {"M001", "M002", "M003", "M004"}, True)[0] == "INCOMPLETE"
    assert check.check_file(_write(tmp_path, "garbage"), "labels", set(), False)[0] == "INCOMPLETE"
