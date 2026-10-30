"""Shared sample and subset definitions for the Ridom exports."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from kansasii_paths import COMPLEX, SCREENING_MAP_RESULTS

# one MST / distance matrix per subset; the minor species are exported together
SUBSETS: dict[str, list[str] | None] = {
    "complex": None,
    "kansasii": ["kansasii"],
    "pseudokansasii": ["pseudokansasii"],
    "persicum": ["persicum"],
    "minor": ["attenuatum", "ostraviense", "innocens", "gastri"],
}
REFERENCE_PREFIX = "GCF_"  # the 7 public GTDB references are named by accession in Ridom


def ridom_id(pnr: str, label: str) -> str:
    """Sample ID in Ridom (kansasii-cgmlst 10_export_ridom.py): Mkan329-081 -> M081, controls -> safe label."""
    if label:
        return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", label)).strip("_")
    return pnr.replace("Mkan329-", "M")


def load_cohort(results: Path = SCREENING_MAP_RESULTS) -> pd.DataFrame:
    """Complex isolates (not contaminated) with columns sample, ridom_id, species, keep (not `excluded`).

    `keep` rows are the cohort that was sent to Ridom (174); the rest is the pool (175, adds M055, excluded by coverage only).
    """
    m = pd.read_csv(results, dtype=str, keep_default_na=False,
                    usecols=["PROBENNUMMER", "LABEL", "species", "contamination_flag", "excluded"])
    m = m[m.species.isin(COMPLEX) & (m.contamination_flag != "contaminated")]
    return pd.DataFrame({"sample": m.PROBENNUMMER,
                         "ridom_id": [ridom_id(p, lab) for p, lab in zip(m.PROBENNUMMER, m.LABEL)],
                         "species": m.species,
                         "keep": m.excluded.str.lower() != "true"}).reset_index(drop=True)
