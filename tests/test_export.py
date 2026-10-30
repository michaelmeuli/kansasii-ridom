"""Ridom archive: pool naming (M055 in), references, archive round trip on tiny FASTA files."""
from __future__ import annotations

import importlib.util
import sys
import tarfile
import zipfile
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from common import REFERENCES, load_cohort  # noqa: E402

_spec = importlib.util.spec_from_file_location("export_ridom", SCRIPTS / "01_export_ridom.py")
assert _spec is not None and _spec.loader is not None
export = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(export)


def _fixture(tmp: Path) -> tuple[Path, Path, Path]:
    res = tmp / "r.csv"
    pd.DataFrame({"PROBENNUMMER": ["Mkan329-001", "Mkan329-055", "Mkan329-010", "Mkan329-126", "Mkan329-128"],
                  "LABEL": ["", "", "", "M. tuberculosis H37Rv", "M. kansasii ATCC 12478"],
                  "species": ["kansasii", "kansasii", "kansasii", "tuberculosis", "kansasii"],
                  "contamination_flag": ["", "", "contaminated", "", ""],
                  "excluded": ["false", "true", "true", "true", "false"]}).to_csv(res, index=False)
    asm = tmp / "asm"
    for s in ("Mkan329-001", "Mkan329-055", "Mkan329-128"):
        (asm / s / "1_unicycler").mkdir(parents=True)
        (asm / s / "1_unicycler" / f"{s}.fasta").write_text(">c\nACGT\n")
        (asm / s / "5_typing/mlst").mkdir(parents=True)
        (asm / s / "5_typing/mlst" / f"mlst_output_{s}.tsv").write_text("FILE\tSCHEME\tST\nx\ty\t12\n")
    refs = tmp / "refs"
    refs.mkdir()
    for acc in REFERENCES:
        (refs / f"{acc}_query.fasta").write_text(">c\nACGT\n")
    return res, asm, refs


def test_pool_includes_m055(tmp_path: Path) -> None:
    res, asm, refs = _fixture(tmp_path)
    chk = pd.Series({"Mkan329-001": "0.1", "Mkan329-055": "0.2", "Mkan329-128": "0.3"})
    df = export.build_table(load_cohort(res), chk, asm, refs)
    cohort = df[~df.is_reference]
    assert sorted(cohort.name) == ["M001", "M055", "M_kansasii_ATCC_12478"]  # contaminated and non-complex out
    assert int(df.is_reference.sum()) == 7 and int(df.is_control.sum()) == 1
    assert df.set_index("name").at["M055", "mlst_ST"] == "12"


def test_archives_round_trip(tmp_path: Path) -> None:
    res, asm, refs = _fixture(tmp_path)
    chk = pd.Series({"Mkan329-001": "0.1", "Mkan329-055": "0.2", "Mkan329-128": "0.3"})
    df = export.build_table(load_cohort(res), chk, asm, refs)
    tar, zp = tmp_path / "a.tar.gz", tmp_path / "a.zip"
    export.write_archives(df, tar, zp)
    with tarfile.open(tar) as tf:
        names = sorted(m.name for m in tf.getmembers() if m.name.endswith(".fasta"))
    with zipfile.ZipFile(zp) as zf:
        assert sorted(i.filename for i in zf.infolist() if i.filename.endswith(".fasta")) == names
    assert "ridom/assemblies/M055.fasta" in names and len(names) == 10
    assert not tar.with_name(tar.name + ".part").exists()
