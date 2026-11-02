"""Scheme definition from the Ridom task template (export 6 of RIDOM.md): core targets, Target Definer settings and filters.

Reads data/imm/ridom-scheme-mkan329.tasktemplate (gzip XML, SeqSphere+ "Export task template"); the Target Definer report is
embedded in it, so no separate parameter export is needed. Writes to output/kansasii-ridom/results/:
  scheme_targets.tsv           core targets: locus tag, gene, product, protein id, length, strand, start/stop
  scheme_loci_classes.tsv      every locus the Target Definer saw: core / accessory / discarded and the filter that decided
  scheme_filters.tsv           per filter: fate, reported and listed number of loci
  scheme_funnel.tsv            candidate targets after each Target Definer step
  scheme_genomes.tsv           seed genomes, penetration query genomes, exclude sequences
  scheme_summary.tsv           key / value table for the methods (versions, counts, thresholds)
  scheme_parameters.txt        Target Definer header, settings, progress, coverage (for the supplement)
  scheme_target_definer_report.txt   the full report incl. locus lists
Optionally checks that the loci of the comparison table are exactly the core targets.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from common import read_comparison
from kansasii_paths import IMM_DIR, OUT_RIDOM
from scheme import LOCUS, parse_filters, parse_funnel, parse_genomes, read_tasktemplate, settings_part


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imm-dir", type=Path, default=IMM_DIR, help="folder with the ridom-* exports")
    ap.add_argument("--tasktemplate", type=Path, help="default: <imm-dir>/ridom-scheme-mkan329.tasktemplate")
    ap.add_argument("--comparison-table", type=Path, help="default: <imm-dir>/ridom-comparison-table.csv")
    ap.add_argument("--out", type=Path, default=OUT_RIDOM / "results")
    a = ap.parse_args()

    s = read_tasktemplate(a.tasktemplate or a.imm_dir / "ridom-scheme-mkan329.tasktemplate")
    a.out.mkdir(parents=True, exist_ok=True)
    targets = pd.DataFrame(s.targets)
    core = list(targets.locus_tag)
    filters = parse_filters(s.report)
    problems: list[str] = []

    classes = [{"locus_tag": t, "class": "core", "filter": ""} for t in core]
    for f in filters:
        loci = [str(x) for x in f["loci"]]  # type: ignore[attr-defined]
        if len(loci) != f["n_reported"]:
            problems.append(f"filter '{f['filter']}': {f['n_reported']} reported, {len(loci)} listed")
        cls = "accessory" if str(f["fate"]).endswith("Accessory") else "discarded"
        classes += [{"locus_tag": t, "class": cls, "filter": f["filter"]} for t in loci]
    cl = pd.DataFrame(classes)
    if cl.locus_tag.duplicated().any():
        problems.append(f"{int(cl.locus_tag.duplicated().sum())} loci in more than one class")
    cl.to_csv(a.out / "scheme_loci_classes.tsv", sep="\t", index=False)
    targets.drop(columns=["reference_sequence"]).to_csv(a.out / "scheme_targets.tsv", sep="\t", index=False)
    pd.DataFrame([{"filter": f["filter"], "fate": f["fate"], "n_reported": f["n_reported"], "n_listed": len(f["loci"])}  # type: ignore[arg-type]
                  for f in filters]).to_csv(a.out / "scheme_filters.tsv", sep="\t", index=False)
    funnel = parse_funnel(s.report)
    pd.DataFrame([*[{"step": k, "targets": v} for k, v in funnel],
                  {"step": "core targets after the stop codon percentage filter", "targets": len(core)}]
                 ).to_csv(a.out / "scheme_funnel.tsv", sep="\t", index=False)
    genomes = pd.DataFrame(parse_genomes(s.report))
    genomes.to_csv(a.out / "scheme_genomes.tsv", sep="\t", index=False)

    rep = s.report
    bases = re.search(r"(\d+) targets were defined for cgMLST \((\d+) bases\)", rep)
    acc = re.search(r"(\d+) targets were used as Accessory targets \((\d+) bases\)", rep)
    disc = re.search(r"(\d+) targets were discarded", rep)
    td = re.search(r"Target Definer version: ([^\n]+)", rep)
    typer = re.search(r"Ridom Typer version: ([^\n]+)", rep)
    date = re.search(r"Date: ([^\n]+)", rep)
    ident = re.search(r"overlap=(\d+)%, identity>=([\d.]+)%", rep)
    stop_pct = re.search(r"in >(\d+)% penetration", rep)
    min_len = re.search(r"Minimum Length Filter \(requires >=(\d+) bases\)", rep)
    cov = re.search(r"\* ([\d.]+)% of Seed genome", rep)
    n_core_bases = int(bases[2]) if bases else 0
    summary = {
        "scheme_name": s.name, "bundle": s.additional.get("Bundle", ""), "uuid": s.additional.get("UUID", ""),
        "template_created": s.created, "template_creator": s.creator,
        "target_definer_version": td[1] if td else "", "ridom_typer_version": typer[1] if typer else "",
        "target_definer_date": date[1] if date else "",
        "seed_genome": next((g["accession"] for g in genomes.to_dict("records") if g["role"] == "seed"), ""),
        "n_query_genomes": int((genomes.role == "query").sum()), "n_exclude_sequences": int((genomes.role == "exclude").sum()),
        "min_target_length_bp": min_len[1] if min_len else "",
        "blast_identity_min_pct": ident[2] if ident else "", "blast_overlap_pct": ident[1] if ident else "",
        "stop_codon_penetration_pct": stop_pct[1] if stop_pct else "",
        "n_core_targets": len(core), "core_bases": n_core_bases,
        "n_accessory_targets": acc[1] if acc else "", "accessory_bases": acc[2] if acc else "",
        "n_discarded_targets": disc[1] if disc else "",
        "core_target_length_median_bp": int(pd.to_numeric(targets.length_bp).median()),
        "core_target_length_min_bp": int(pd.to_numeric(targets.length_bp).min()),
        "core_target_length_max_bp": int(pd.to_numeric(targets.length_bp).max()),
        "seed_genome_coverage_pct": cov[1] if cov else "",
        **{f"ngs_{k}": v for k, v in s.ngs_params.items()},
    }
    if bases and int(bases[1]) != len(core):
        problems.append(f"report says {bases[1]} core targets, template holds {len(core)}")
    if n_core_bases != int(pd.to_numeric(targets.length_bp).sum()):
        problems.append(f"core bases: report {n_core_bases}, sum of target lengths {int(pd.to_numeric(targets.length_bp).sum())}")
    n_acc = int((cl["class"] == "accessory").sum())
    n_dis = int((cl["class"] == "discarded").sum())
    if acc and int(acc[1]) != n_acc or disc and int(disc[1]) != n_dis:
        problems.append(f"accessory/discarded listed {n_acc}/{n_dis}, report {acc[1] if acc else '?'}/{disc[1] if disc else '?'}")
    table = a.comparison_table or a.imm_dir / "ridom-comparison-table.csv"
    if table.exists():
        in_table = set(read_comparison(table)[1].columns)
        summary["comparison_table_loci_equal_core"] = str(in_table == set(core)).lower()
        if in_table != set(core):
            problems.append(f"comparison table has {len(in_table - set(core))} loci not in the scheme, lacks {len(set(core) - in_table)}")
    pd.Series(summary, name="value").rename_axis("key").reset_index().to_csv(a.out / "scheme_summary.tsv", sep="\t", index=False)

    head = f"Scheme: {s.name}\nTask template created: {s.created} by {s.creator}\nBundle: {s.additional.get('Bundle', '')}\n\n"
    (a.out / "scheme_parameters.txt").write_text(head + settings_part(rep))
    (a.out / "scheme_target_definer_report.txt").write_text(head + rep)
    print(pd.Series(summary).to_string())
    for p in problems:
        print(f"ERROR: {p}", file=sys.stderr)
    assert all(LOCUS.fullmatch(t) for t in core)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
