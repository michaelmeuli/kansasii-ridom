"""Shared sample and subset definitions for the Ridom exports."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from kansasii_paths import COHORT_SPECIES, OUT_RIDOM, SCREENING_MAP_RESULTS

# one MST / distance matrix per subset; the minor species are exported together
SUBSETS: dict[str, list[str] | None] = {
    "complex": None,
    "kansasii": ["kansasii"],
    "pseudokansasii": ["pseudokansasii"],
    "persicum": ["persicum"],
    "minor": ["attenuatum", "ostraviense", "innocens", "gastri"],
}
# distance matrices are sliced in code, so the minor species also get one matrix each (MST/tree exports keep SUBSETS)
DISTANCE_SUBSETS: dict[str, list[str] | None] = {
    **SUBSETS, "attenuatum": ["attenuatum"], "ostraviense": ["ostraviense"], "innocens": ["innocens"], "gastri": ["gastri"],
}
REFERENCE_PREFIX = "GCF_"  # the 7 public GTDB references are named by accession in Ridom
RIDOM_TARBALL = OUT_RIDOM / "ridom_assemblies.tar.gz"  # input of SeqSphere+ (01_export_ridom.py)
RIDOM_ZIP = OUT_RIDOM / "ridom_assemblies.zip"
# public GTDB representatives (accession, species) sent to Ridom as references; the first is the ATCC 12478 type strain
REFERENCES = {
    "GCF_000157895.3": "kansasii", "GCF_002086675.1": "persicum", "GCF_900566075.1": "pseudokansasii",
    "GCF_900566085.1": "attenuatum", "GCF_900566055.1": "innocens", "GCF_002705925.1": "ostraviense",
    "GCF_002102175.1": "gastri",
}


def ridom_id(pnr: str, label: str) -> str:
    """Sample ID in Ridom (kansasii-cgmlst 10_export_ridom.py): Mkan329-081 -> M081, controls -> safe label."""
    if label:
        return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", label)).strip("_")
    return pnr.replace("Mkan329-", "M")


def load_cohort(results: Path = SCREENING_MAP_RESULTS) -> pd.DataFrame:
    """Complex isolates (not contaminated) with columns sample, ridom_id, label, species, keep (not `excluded`).

    `keep` rows are the cohort that was sent to Ridom (174); the rest is the pool (175, adds M055, excluded by coverage only).
    """
    m = pd.read_csv(results, dtype=str, keep_default_na=False,
                    usecols=["PROBENNUMMER", "LABEL", "species", "contamination_flag", "excluded"])
    m = m[m.species.isin(COHORT_SPECIES) & (m.contamination_flag != "contaminated")]
    return pd.DataFrame({"sample": m.PROBENNUMMER,
                         "ridom_id": [ridom_id(p, lab) for p, lab in zip(m.PROBENNUMMER, m.LABEL)],
                         "label": m.LABEL,
                         "species": m.species,
                         "keep": m.excluded.str.lower() != "true"}).reset_index(drop=True)


def meg_label(raw: str) -> str:
    """Sample ID of a MEGA label: Ridom writes the accessions as GCF_000157895_3 for GCF_000157895.3."""
    return re.sub(r"^(GCF_\d+)_(\d+)$", r"\1.\2", raw)


def read_meg(path: Path) -> pd.DataFrame:
    """Square distance matrix (int, index = columns = Sample ID) from a SeqSphere+ MEGA export (`DataFormat=LowerLeft`).

    Labels (`#id`) come first, then the lower triangle without the diagonal. Raises ValueError if the file does not
    match its own header (NTaxa, number of values) or has duplicate labels.
    """
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if not lines or lines[0].strip().lower() != "#mega":
        raise ValueError(f"{path}: not a MEGA file (first line must be #mega)")
    fmt = " ".join(ln for ln in lines if ln.lower().startswith("!format"))
    m = re.search(r"NTaxa=(\d+)", fmt, re.I)
    if not m or not re.search(r"DataFormat=LowerLeft", fmt, re.I) or not re.search(r"DataType=distance", fmt, re.I):
        raise ValueError(f"{path}: expected '!Format DataType=distance DataFormat=LowerLeft NTaxa=n', got '{fmt}'")
    n = int(m.group(1))
    labels = [meg_label(ln[1:].strip()) for ln in lines[1:] if ln.startswith("#")]
    if len(labels) != n or len(set(labels)) != n:
        raise ValueError(f"{path}: NTaxa={n} but {len(labels)} labels ({len(set(labels))} unique)")
    last = max(i for i, ln in enumerate(lines) if ln.startswith("#"))
    values = [float(t) for ln in lines[last + 1:] for t in ln.split()]
    if len(values) != n * (n - 1) // 2:
        raise ValueError(f"{path}: {len(values)} values, lower triangle of {n} taxa needs {n * (n - 1) // 2}")
    mat = np.zeros((n, n))
    mat[np.tril_indices(n, -1)] = values  # row-major: row 1 has 1 value, row 2 has 2, ...
    mat = mat + mat.T
    if not np.array_equal(mat, np.rint(mat)):
        raise ValueError(f"{path}: allele distances must be integers")
    return pd.DataFrame(mat.astype(int), index=labels, columns=labels)


COMPARISON_META = {"#Missing values in Distance Columns": "missing_ridom", "Perc. Good Targets": "good_pct",
                   "Approximated Genome Size (Mbases)": "genome_mbases", "Species": "species_ridom"}


def sniff_sep(path: Path) -> str:
    """Ridom writes `;`-separated CSV; accept `,` too."""
    header = path.read_text(encoding="utf-8-sig").splitlines()[0]
    return ";" if header.count(";") > header.count(",") else ","


def read_comparison(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Comparison table of SeqSphere+ -> (meta, alleles), both indexed by Sample ID.

    meta: missing_ridom (int), good_pct, genome_mbases (float), species_ridom. alleles: one column per target, allele number
    as text, `""` = no allele called (missing / failed; Ridom does not export the reason).
    """
    t = pd.read_csv(path, sep=sniff_sep(path), encoding="utf-8-sig", dtype=str, keep_default_na=False)
    t = t.loc[:, [c for c in t.columns if c and not c.startswith("Unnamed")]]
    absent = [c for c in [*COMPARISON_META, "Sample ID"] if c not in t.columns]
    if absent:
        raise ValueError(f"{path}: columns missing: {absent}")
    if t["Sample ID"].duplicated().any():
        raise ValueError(f"{path}: duplicate Sample ID")
    t = t.set_index("Sample ID")
    meta = t[list(COMPARISON_META)].rename(columns=COMPARISON_META)
    meta = meta.astype({"missing_ridom": int, "good_pct": float, "genome_mbases": float})
    return meta, t.drop(columns=list(COMPARISON_META))


def read_clusters(path: Path) -> pd.Series[str]:
    """Sample ID -> MST cluster label of the Ridom cluster export (`""` = not in a cluster)."""
    t = pd.read_csv(path, sep=sniff_sep(path), encoding="utf-8-sig", dtype=str, keep_default_na=False)
    return pd.Series(t["-Group-"].to_numpy(), index=t["Sample ID"].to_numpy(), dtype=str)


def kind_of(sample: str, cohort: pd.DataFrame) -> str:
    """`reference` (GCF_*), `control` (strain with a LABEL) or `isolate`."""
    if sample.startswith(REFERENCE_PREFIX):
        return "reference"
    row = cohort[cohort.ridom_id == sample]
    return "control" if len(row) and row.label.iloc[0] else "isolate"


def species_of(cohort: pd.DataFrame) -> dict[str, str]:
    """Sample ID -> species for the cohort and the GTDB references."""
    return {**dict(zip(cohort.ridom_id, cohort.species)), **REFERENCES}
