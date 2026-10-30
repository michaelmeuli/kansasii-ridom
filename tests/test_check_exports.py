"""Export check on small hand-made tables."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from common import load_cohort, ridom_id  # noqa: E402

_spec = importlib.util.spec_from_file_location("check_exports", SCRIPTS / "00_check_exports.py")
assert _spec is not None and _spec.loader is not None
check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check)


def test_ridom_id() -> None:
    assert ridom_id("Mkan329-081", "") == "M081"
    assert ridom_id("Mkan329-128", "M. kansasii ATCC 12478") == "M_kansasii_ATCC_12478"


def _results(tmp: Path) -> Path:
    p = tmp / "r.csv"
    pd.DataFrame({"PROBENNUMMER": ["Mkan329-001", "Mkan329-002", "Mkan329-055", "Mkan329-010", "Mkan329-126"],
                  "LABEL": ["", "", "", "", "M. tuberculosis H37Rv"],
                  "species": ["kansasii", "persicum", "kansasii", "kansasii", "tuberculosis"],
                  "contamination_flag": ["", "", "", "contaminated", ""],
                  "excluded": ["false", "false", "true", "true", "true"]}).to_csv(p, index=False)
    return p


def test_cohort_keep_and_pool(tmp_path: Path) -> None:
    c = load_cohort(_results(tmp_path))
    assert set(c.ridom_id) == {"M001", "M002", "M055"}  # contaminated and non-complex dropped
    assert set(c[c.keep].ridom_id) == {"M001", "M002"}
    assert check.expected_ids(c, "pool", None) == {"M001", "M002", "M055"}
    assert check.expected_ids(c, "keep", ["persicum"]) == {"M002"}


def test_check_file(tmp_path: Path) -> None:
    f = tmp_path / "t.csv"
    f.write_text("Sample ID;Species\nM001;x\nGCF_1;y\nM009;z\n")
    assert check.check_file(tmp_path / "none.csv", "Sample ID", set(), False)[0] == "MISSING"
    assert check.check_file(f, "", set(), False)[0] == "UNCHECKED"
    assert check.check_file(f, "Sample ID", {"M001", "M002"}, True)[0] == "INCOMPLETE"
    assert check.check_file(f, "Sample ID", {"M001"}, True)[0] == "EXTRA"  # M009 is not in the cohort
    assert check.check_file(f, "Sample ID", {"M001", "M009"}, True)[0] == "OK"
