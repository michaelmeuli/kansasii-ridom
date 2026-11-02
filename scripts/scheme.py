"""Read a SeqSphere+ task template (`*.tasktemplate`): the cgMLST scheme and the embedded Target Definer report.

The file is gzip-compressed XML. One `PCRBankEntry type=group` holds the scheme; each core target is a child
`PCRBankEntry type=template`. The Target Definer report is the `<Comment>` that starts with `cgMLST Target Definer Results`.
"""
from __future__ import annotations

import gzip
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

REPORT_START = "cgMLST Target Definer Results"
LOCUS = re.compile(r"MKAN_RS\d+")


@dataclass
class Scheme:
    name: str
    created: str
    creator: str
    ngs_params: dict[str, str]
    additional: dict[str, str]
    targets: list[dict[str, str]]
    report: str = field(repr=False)


def _text(e: ET.Element, path: str) -> str:
    x = e.find(path)
    return (x.text or "").strip() if x is not None else ""


def read_tasktemplate(path: Path) -> Scheme:
    """Parse the scheme: raises ValueError if the file is not a cgMLST task template or has no Target Definer report."""
    with gzip.open(path, "rb") as fh:
        root = ET.parse(fh).getroot()
    groups = [e for e in root.iter("PCRBankEntry") if e.get("type") == "group"]
    if root.tag != "TaskTemplate" or len(groups) != 1:
        raise ValueError(f"{path}: expected one TaskTemplate with one group of targets, found {len(groups)} groups")
    g = groups[0]
    targets = []
    for e in g.findall("PCRBankEntry"):
        props = {k.text or "": (v.text or "") for k, v in zip(e.iter("Key"), e.iter("Value"))}
        area = e.find(".//RefAreas/RefArea")
        seq = _text(e, ".//ReferenceSequence/Bases")
        targets.append({
            "locus_tag": props.get("GenBank locus_tag") or _text(e, "Name"),
            "gene": props.get("GenBank gene", ""),
            "product": props.get("GenBank product", ""),
            "protein_id": props.get("GenBank protein_id", ""),
            "gene_id": props.get("GenBank GeneID", ""),
            "length_bp": str(len(seq)),
            "strand": ("+" if _text(area, "Forward") == "1" else "-") if area is not None else "",
            "codon_table": _text(e, ".//CodonTable"),
            "has_start": _text(area, "ContainsStart") if area is not None else "",
            "has_stop": _text(area, "ContainsStop") if area is not None else "",
            "reference_sequence": _text(e, ".//ReferenceSequence/Name"),
        })
    report = next((c.text or "" for c in root.iter("Comment") if (c.text or "").startswith(REPORT_START)), "")
    if not report:
        raise ValueError(f"{path}: no '{REPORT_START}' comment (export the scheme as a task template from the Target Definer result)")
    additional = {i.get("key", ""): i.get("value", "") for i in root.iter("item")}
    return Scheme(name=_text(g, "Name"), created=g.get("created", ""), creator=g.get("creator", ""),
                  ngs_params={c.tag: (c.text or "") for c in (g.find("NgsParams") or [])},
                  additional=additional, targets=targets, report=report.replace("\r\n", "\n").replace("\r", "\n"))


def settings_part(report: str) -> str:
    """The report up to the per-target lists (header, settings, progress, coverage)."""
    i = report.find("TARGET DETAILS")
    return report if i < 0 else report[:i].rstrip() + "\n"


def parse_genomes(report: str) -> list[dict[str, str]]:
    """Seed genomes, penetration query genomes and exclude sequences of the SETTINGS block."""
    role = ""
    rows: list[dict[str, str]] = []
    pat = re.compile(r"^ \* (?P<acc>.+?), (?P<bases>\d+) bases(?:, (?P<n>\d+) (?P<unit>contigs|genes with CDS))? \((?P<desc>.*)\)$")
    for ln in settings_part(report).split("DETAILS")[0].splitlines():
        if ln.startswith("Seed Genome Sources"):
            role = "seed"
        elif ln.startswith("Penetration Query Genomes ("):
            role = "query"
        elif ln.startswith("Exclude Sequences ("):
            role = "exclude"
        elif ln.startswith("Seed Genome Filters") or ln.startswith("Query Genome BLAST"):
            role = ""
        m = pat.match(ln)
        if role and m:
            rows.append({"role": role, "accession": m["acc"], "bases": m["bases"], "n": m["n"] or "", "unit": m["unit"] or "",
                         "description": m["desc"]})
    return rows


def parse_funnel(report: str) -> list[tuple[str, int]]:
    """Number of candidate targets after each Target Definer step (the `Progress:` block)."""
    out: list[tuple[str, int]] = []
    on = False
    for ln in report.splitlines():
        if ln.startswith("Progress:"):
            on = True
        elif on and not ln.strip():
            break
        elif on:
            m = re.match(r"^ \* (\d+) targets (.*)$", ln)
            if m:
                out.append((m[2], int(m[1])))
    return out


def parse_filters(report: str) -> list[dict[str, object]]:
    """Per filter: name, fate (`discarded` or `moved to Accessory`), reported count and the listed locus tags."""
    i = report.find("TARGET DETAILS")
    if i < 0:
        return []
    out: list[dict[str, object]] = []
    for block in re.split(r"\n\s*\n", report[i:]):
        m = re.match(r"^\s*(?P<name>[^\n]+?): filtered out (?P<n>\d+) targets \((?P<fate>[^)]+)\)\n", block)
        if m:
            loci: list[str] = []
            for ln in block[m.end():].splitlines():  # the list ends at the first line that is not locus tags (e.g. `Messages for ...`)
                if not re.fullmatch(r"(MKAN_RS\d+,?\s*)+", ln):
                    break
                loci += LOCUS.findall(ln)
            out.append({"filter": m["name"], "fate": m["fate"], "n_reported": int(m["n"]), "loci": loci})
    return out
