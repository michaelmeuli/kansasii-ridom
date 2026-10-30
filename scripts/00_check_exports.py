"""Compare the Ridom exports in data/imm with the list of exports still needed and the isolates they must hold.

Writes output/kansasii-ridom/results/export_status.tsv: one row per expected file (status MISSING, OK, INCOMPLETE,
EXTRA, UNCHECKED). A file with an id column is checked against the cohort (screening_map_results.csv `excluded`), the
7 GTDB references are expected as GCF_* ids. Ridom's own species names (M. kansasii) are not used for the check.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from common import REFERENCE_PREFIX, SUBSETS, load_cohort
from kansasii_paths import IMM_DIR, OUT_RIDOM, SCREENING_MAP_RESULTS

EXPECTED = Path(__file__).resolve().parents[1] / "config/expected_exports.tsv"
SHOWN = 8  # ids listed per status row


def expected_ids(cohort: pd.DataFrame, which: str, species: list[str] | None) -> set[str]:
    c = cohort if which == "pool" else cohort[cohort.keep]
    if species is not None:
        c = c[c.species.isin(species)]
    return set(c.ridom_id)


def check_file(path: Path, id_column: str, want: set[str], refs: bool) -> tuple[str, str]:
    """Status and detail of one export."""
    if not path.exists():
        return "MISSING", ""
    if not id_column:
        return "UNCHECKED", "format not seen yet, only presence checked"
    header = path.read_text(encoding="utf-8-sig").splitlines()[0]
    sep = ";" if header.count(";") > header.count(",") else ","
    ids = set(pd.read_csv(path, sep=sep, encoding="utf-8-sig", usecols=[id_column], dtype=str)[id_column].dropna())
    have_refs = {i for i in ids if i.startswith(REFERENCE_PREFIX)}
    missing = sorted(want - ids)
    extra = sorted(ids - want - (have_refs if refs else set()))
    if refs and not have_refs:
        missing.append(f"{REFERENCE_PREFIX}* references")
    n = f"{len(ids)} ids"
    if missing:
        return "INCOMPLETE", f"{n}; {len(missing)} missing: {', '.join(missing[:SHOWN])}{' ...' if len(missing) > SHOWN else ''}"
    if extra:
        return "EXTRA", f"{n}; {len(extra)} not in cohort: {', '.join(extra[:SHOWN])}"
    return "OK", n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    ap.add_argument("--strict", action="store_true", help="exit 1 unless every expected export is OK or UNCHECKED")
    a = ap.parse_args()

    cohort = load_cohort(a.screening_map_results)
    spec = pd.read_csv(EXPECTED, sep="\t", comment="#", dtype=str, keep_default_na=False)
    rows = []
    for r in spec.to_dict("records"):
        subsets = list(SUBSETS) if r["subsets"] == "all" else [""]
        refs = r["refs"] == "yes"
        for s in subsets:
            name = str(r["file"]).format(subset=s)
            want = expected_ids(cohort, str(r["cohort"]), SUBSETS.get(s))
            status, detail = check_file(a.imm_dir / name, str(r["id_column"]), want, refs)
            rows.append({"no": r["no"], "file": name, "status": status, "detail": detail, "replaces": r["replaces"]})
    out = pd.DataFrame(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out / "export_status.tsv", sep="\t", index=False)
    print(out.groupby("status").size().to_string())
    print(out[out.status.isin(["INCOMPLETE", "EXTRA"])].to_string(index=False))
    bad = out.status.isin(["MISSING", "INCOMPLETE", "EXTRA"]).any()
    return 1 if (a.strict and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
