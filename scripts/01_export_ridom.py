"""Build output/kansasii-ridom/ridom_assemblies.tar.gz and .zip (same content) for Ridom SeqSphere+: the complex isolates that are not
contaminated (the pool, 175 incl. M055 which is `excluded` for low chewBBACA coverage only; assemblies as real files named M001, the
ATCC/DSM control strains by their LABEL) + the 7 public GTDB references (named by accession, unchanged: Ridom removes plasmids), and
metadata.tsv for the isolate-data import. Same layout and metadata as the earlier kansasii-cgmlst/scripts/10_export_ridom.py,
which left out M055."""
from __future__ import annotations

import argparse
import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

from common import REFERENCES, RIDOM_TARBALL, RIDOM_ZIP, load_cohort
from kansasii_paths import ASSEMBLIES, GTDB_REPRESENTATIVES, SCREENING_MAP_RESULTS

N_POOL, N_REFERENCES, N_CONTROLS = 175, 7, 3  # expected sizes (docs/EXCLUSIONS.md)


def species_label(sp: str) -> str:
    return f"M. {sp}"


def mlst_st(assemblies: Path, sample: str) -> str:
    """Sequence type of the pipeline's mlst run (scheme mycobacteria_2); '-' when mlst assigns none."""
    f = assemblies / sample / "5_typing" / "mlst" / f"mlst_output_{sample}.tsv"
    return f.read_text(encoding="utf-8").splitlines()[1].split("\t")[2]


def build_table(cohort: pd.DataFrame, checkm: pd.Series, assemblies: Path, references: Path) -> pd.DataFrame:
    """One row per genome sent to Ridom: name, sample, species, mlst_ST, checkm_contamination, is_control, is_reference, src."""
    rows: list[dict[str, object]] = []
    for sample, name, label, species in zip(cohort["sample"], cohort["ridom_id"], cohort["label"], cohort["species"]):
        rows.append({"name": name, "sample": sample, "species": species_label(species), "mlst_ST": mlst_st(assemblies, sample),
                     "checkm_contamination": checkm.at[sample], "is_control": bool(label), "is_reference": False,
                     "src": assemblies / sample / "1_unicycler" / f"{sample}.fasta"})
    for acc, sp in REFERENCES.items():
        rows.append({"name": acc, "sample": acc, "species": species_label(sp), "mlst_ST": "", "checkm_contamination": "",
                     "is_control": False, "is_reference": True, "src": references / f"{acc}_query.fasta"})
    df = pd.DataFrame(rows)
    assert df.name.is_unique, df.name[df.name.duplicated()].tolist()
    return df


def write_archives(df: pd.DataFrame, tarball: Path, zipfile_: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="ridom_export_") as tmp:
        top = Path(tmp) / "ridom"
        (top / "assemblies").mkdir(parents=True)
        for name, src in zip(df["name"].astype(str), df["src"]):
            dst = top / "assemblies" / f"{name}.fasta"
            shutil.copyfile(str(src), dst)  # follows symlinks: the archives hold real files
            assert dst.stat().st_size > 0, src
        df.drop(columns="src").to_csv(top / "metadata.tsv", sep="\t", index=False)
        tarball.parent.mkdir(parents=True, exist_ok=True)
        tmp_tar = tarball.with_name(tarball.name + ".part")  # renamed only when complete
        with tarfile.open(tmp_tar, "w:gz") as tf:
            tf.add(top, arcname="ridom", recursive=False)
            tf.add(top / "assemblies", arcname="ridom/assemblies", recursive=False)
            for f in sorted((top / "assemblies").iterdir()):
                tf.add(f, arcname=f"ridom/assemblies/{f.name}")
            tf.add(top / "metadata.tsv", arcname="ridom/metadata.tsv")
        tmp_tar.replace(tarball)
        tmp_zip = zipfile_.with_name(zipfile_.name + ".part")
        with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for f in sorted((top / "assemblies").iterdir()):
                zf.write(f, arcname=f"ridom/assemblies/{f.name}")
            zf.write(top / "metadata.tsv", arcname="ridom/metadata.tsv")
        tmp_zip.replace(zipfile_)


def verify(df: pd.DataFrame, tarball: Path, zipfile_: Path) -> None:
    with tarfile.open(tarball) as tf:
        fasta = [m for m in tf.getmembers() if m.name.endswith(".fasta")]
        assert all(m.isfile() and m.size > 0 for m in fasta), "link or empty file in tarball"
        assert len(fasta) == len(df), (len(fasta), len(df))
    with zipfile.ZipFile(zipfile_) as zf:
        zfasta = [i for i in zf.infolist() if i.filename.endswith(".fasta")]
        assert all(i.file_size > 0 for i in zfasta), "empty file in zip"
        assert sorted(i.filename for i in zfasta) == sorted(m.name for m in fasta), "zip and tarball differ"
    assert int(df.is_reference.sum()) == N_REFERENCES and int(df.is_control.sum()) == N_CONTROLS, "references/controls"
    assert len(df) == N_POOL + N_REFERENCES and "M055" in set(df.name), "pool incl. M055 expected"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES, help="<sample>/1_unicycler/ and 5_typing/mlst/ per sample")
    ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS, help="pool and CheckM contamination")
    ap.add_argument("--gtdb-representatives", type=Path, default=GTDB_REPRESENTATIVES, help="<accession>_query.fasta of the 7 references")
    ap.add_argument("--tarball", type=Path, default=RIDOM_TARBALL)
    ap.add_argument("--zip", type=Path, default=RIDOM_ZIP)
    a = ap.parse_args()

    chk = pd.read_csv(a.screening_map_results, dtype=str).fillna("").set_index("PROBENNUMMER")["checkm_contamination"]
    df = build_table(load_cohort(a.screening_map_results), chk, a.assemblies, a.gtdb_representatives)
    write_archives(df, a.tarball, a.zip)
    verify(df, a.tarball, a.zip)
    print(f"{a.tarball} + {a.zip.name}: {len(df)} genomes ({int(df.is_reference.sum())} references, {int(df.is_control.sum())} controls), "
          f"{a.tarball.stat().st_size / 1e6:.1f} MB; species:\n{df.species.value_counts().to_string()}")


if __name__ == "__main__":
    main()
