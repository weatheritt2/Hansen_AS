#!/usr/bin/env python3
"""
interface_overlap_fixed_columns.py

Strict fixed-column overlap annotator for interface resources.

Query file (TSV/CSV) must contain columns:
    uniprot_canonical   start_aa_position   end_aa_position
(plus optional query_id)

Per-source fixed AA positions (1-based indexing):
- IntAct:      aa_start=36, aa_end=37
- 3DID:        pair1=(11,12), pair2=(18,19)
- DOMINO:      aa_start=32, aa_end=33
- Interactome3D:
      protein columns = 5 and 6
      AA columns provided by --i3d-aa-cols (default: 11,12,18,19)
- ELM:         pair1=(4,5), pair2=(6,7)
- PDB:         user provides protein and AA columns via CLI (no guessing)

Important:
- This script does NOT infer columns from headers for these core fields.
- It uses fixed positional columns as requested.

Data sources (download pages):
  - IntAct:          https://www.ebi.ac.uk/intact/download
  - 3DID:            https://3did.irbbarcelona.org/download.php
  - Interactome3D:   https://interactome3d.irbbarcelona.org/download.php
  - ELM:             http://elm.eu.org/
  - DOMINO:          http://mint.bio.uniroma2.it/domino/
  - PDB:             https://www.rcsb.org/downloads

You can run the program with ANY single database file (e.g. only IntAct).
Only the provided sources are loaded/queried.
"""

from __future__ import annotations
import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict


# ----------------------------
# Helpers
# ----------------------------

def sniff_delimiter(path: str) -> str:
    with open(path, "r", encoding="utf-8") as fh:
        sample = fh.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except Exception:
        return "\t"

def norm_uniprot(acc: str) -> str:
    if acc is None:
        return ""
    acc = str(acc).strip().upper()
    return acc.split("-")[0]

def to_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None
    # handle ranges like "123-145" elsewhere; here pure int only
    try:
        return int(float(s))
    except Exception:
        return None

def parse_range_token(tok: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Accept:
      "123" -> (123,123)
      "123-145" -> (123,145)
      "123..145" -> (123,145)
    """
    if tok is None:
        return (None, None)
    s = str(tok).strip()
    if s == "":
        return (None, None)

    s = s.replace("..", "-")
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", s)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        return (min(a, b), max(a, b))
    v = to_int(s)
    if v is None:
        return (None, None)
    return (v, v)

def colval(row: List[str], col_1based: int) -> str:
    idx = col_1based - 1
    if idx < 0 or idx >= len(row):
        return ""
    return row[idx].strip()

def overlap_len(a1: int, a2: int, b1: int, b2: int) -> int:
    s = max(a1, b1)
    e = min(a2, b2)
    return max(0, e - s + 1)

def summarize(vals: List[float]) -> str:
    if not vals:
        return "NA"
    vs = sorted(vals)
    n = len(vs)
    def q(p):
        k = int(round((n - 1) * p))
        return vs[k]
    return f"min={vs[0]:.4f};q25={q(0.25):.4f};median={q(0.50):.4f};q75={q(0.75):.4f};max={vs[-1]:.4f};n={n}"


@dataclass(frozen=True)
class InterfaceRec:
    source: str
    protein: str
    start: int
    end: int
    partner: str = ""
    interface_id: str = ""
    record_id: str = ""
    evidence: str = ""
    taxon: str = ""


def parse_queries(path: str) -> List[Dict[str, Any]]:
    delim = sniff_delimiter(path)
    req = ["uniprot_canonical", "start_aa_position", "end_aa_position"]
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        cols = rdr.fieldnames or []
        missing = [c for c in req if c not in cols]
        if missing:
            raise ValueError(f"Query file missing required columns: {missing}. Found: {cols}")
        for i, r in enumerate(rdr, start=1):
            p = norm_uniprot(r["uniprot_canonical"])
            s = to_int(r["start_aa_position"])
            e = to_int(r["end_aa_position"])
            if not p or s is None or e is None or s <= 0 or e <= 0:
                continue
            if s > e:
                s, e = e, s
            qid = r.get("query_id", f"Q{i}")
            out.append({
                "query_id": qid,
                "uniprot_canonical": p,
                "start_aa_position": s,
                "end_aa_position": e
            })
    return out


def read_table_rows(path: str) -> List[List[str]]:
    delim = sniff_delimiter(path)
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        for row in rdr:
            rows.append(row)
    return rows


def maybe_skip_header(rows: List[List[str]], aa_cols: List[int]) -> List[List[str]]:
    if not rows:
        return rows
    first = rows[0]
    # If all AA cols in first row fail parse as int/range, assume header
    all_fail = True
    for c in aa_cols:
        tok = colval(first, c)
        s, e = parse_range_token(tok)
        if s is not None and e is not None:
            all_fail = False
            break
    return rows[1:] if all_fail else rows


def add_rec(out: List[InterfaceRec], rec: InterfaceRec, dedup: set):
    if not rec.protein or rec.start <= 0 or rec.end <= 0:
        return
    key = (rec.source, rec.protein, rec.start, rec.end, rec.interface_id or rec.record_id, rec.partner)
    if key in dedup:
        return
    dedup.add(key)
    out.append(rec)


def parse_single_side_source(
    source: str,
    path: str,
    protein_col: int,
    aa_start_col: int,
    aa_end_col: int,
    partner_col: Optional[int] = None,
    interface_id_col: Optional[int] = None,
    record_id_col: Optional[int] = None,
) -> List[InterfaceRec]:
    rows = read_table_rows(path)
    rows = maybe_skip_header(rows, [aa_start_col, aa_end_col])
    out, dedup = [], set()

    for r in rows:
        prot = norm_uniprot(colval(r, protein_col))
        if not prot:
            continue

        s = to_int(colval(r, aa_start_col))
        e = to_int(colval(r, aa_end_col))
        # if not plain int, try range token in start col
        if s is None or e is None:
            rs, re_ = parse_range_token(colval(r, aa_start_col))
            if rs is not None and re_ is not None:
                s, e = rs, re_
        if s is None or e is None:
            continue
        if s > e:
            s, e = e, s

        partner = norm_uniprot(colval(r, partner_col)) if partner_col else ""
        iid = colval(r, interface_id_col) if interface_id_col else ""
        rid = colval(r, record_id_col) if record_id_col else ""

        add_rec(out, InterfaceRec(source, prot, s, e, partner=partner, interface_id=iid, record_id=rid), dedup)

    return out


def parse_dual_partner_source(
    source: str,
    path: str,
    protA_col: int,
    protB_col: int,
    A_start_col: int,
    A_end_col: int,
    B_start_col: int,
    B_end_col: int,
    interface_id_col: Optional[int] = None,
    record_id_col: Optional[int] = None,
) -> List[InterfaceRec]:
    rows = read_table_rows(path)
    rows = maybe_skip_header(rows, [A_start_col, A_end_col, B_start_col, B_end_col])
    out, dedup = [], set()

    for r in rows:
        pa = norm_uniprot(colval(r, protA_col))
        pb = norm_uniprot(colval(r, protB_col))
        iid = colval(r, interface_id_col) if interface_id_col else ""
        rid = colval(r, record_id_col) if record_id_col else ""

        # side A
        sa, ea = to_int(colval(r, A_start_col)), to_int(colval(r, A_end_col))
        if sa is None or ea is None:
            rsa, rea = parse_range_token(colval(r, A_start_col))
            if rsa is not None and rea is not None:
                sa, ea = rsa, rea
        if pa and sa is not None and ea is not None:
            if sa > ea:
                sa, ea = ea, sa
            add_rec(out, InterfaceRec(source, pa, sa, ea, partner=pb, interface_id=iid, record_id=rid), dedup)

        # side B
        sb, eb = to_int(colval(r, B_start_col)), to_int(colval(r, B_end_col))
        if sb is None or eb is None:
            rsb, reb = parse_range_token(colval(r, B_start_col))
            if rsb is not None and reb is not None:
                sb, eb = rsb, reb
        if pb and sb is not None and eb is not None:
            if sb > eb:
                sb, eb = eb, sb
            add_rec(out, InterfaceRec(source, pb, sb, eb, partner=pa, interface_id=iid, record_id=rid), dedup)

    return out


def index_by_protein(recs: List[InterfaceRec]) -> Dict[str, List[InterfaceRec]]:
    d = defaultdict(list)
    for x in recs:
        d[x.protein].append(x)
    return d


def annotate(queries: List[Dict[str, Any]], sources: Dict[str, List[InterfaceRec]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    idx = {s: index_by_protein(v) for s, v in sources.items()}
    source_names = list(sources.keys())

    summary = []
    longrows = []

    for q in queries:
        qid = q["query_id"]
        prot = q["uniprot_canonical"]
        qs = q["start_aa_position"]
        qe = q["end_aa_position"]
        qlen = qe - qs + 1

        row = dict(q)
        all_frac_iface = []
        total_hits = 0

        for src in source_names:
            hits = []
            for rec in idx[src].get(prot, []):
                ol = overlap_len(qs, qe, rec.start, rec.end)
                if ol > 0:
                    ilen = rec.end - rec.start + 1
                    fq = ol / qlen
                    fi = ol / ilen
                    hits.append((rec, ol, fq, fi))
                    all_frac_iface.append(fi)

                    longrows.append({
                        "query_id": qid,
                        "uniprot_canonical": prot,
                        "query_start": qs,
                        "query_end": qe,
                        "source": src,
                        "partner_uniprot": rec.partner,
                        "interface_id": rec.interface_id,
                        "record_id": rec.record_id,
                        "interface_start": rec.start,
                        "interface_end": rec.end,
                        "interface_len": ilen,
                        "overlap_residue_count": ol,
                        "overlap_fraction_of_query_interval": round(fq, 6),
                        "overlap_fraction_of_interface": round(fi, 6),
                    })

            row[f"{src}_overlap_binary"] = 1 if hits else 0
            row[f"{src}_matched_interface_count"] = len(hits)
            row[f"{src}_matched_interface_ids"] = ";".join([(h[0].interface_id or h[0].record_id or f"{h[0].start}-{h[0].end}") for h in hits])
            row[f"{src}_matched_ranges"] = ";".join([f"{h[0].start}-{h[0].end}" for h in hits])
            total_hits += len(hits)

        row["total_matched_interfaces_all_sources"] = total_hits
        row["any_overlap_binary"] = 1 if total_hits > 0 else 0
        row["overlap_fraction_of_interface_distribution_all_sources"] = summarize(all_frac_iface)
        summary.append(row)

    return summary, longrows


def write_tsv(path: str, rows: List[Dict[str, Any]]):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        return
    cols = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(
        description="Fixed-column interface overlap annotator (any subset of sources).",
        epilog=(
            "Download pages:\n"
            "  IntAct:        https://www.ebi.ac.uk/intact/download\n"
            "  3DID:          https://3did.irbbarcelona.org/download.php\n"
            "  Interactome3D: https://interactome3d.irbbarcelona.org/download.php\n"
            "  ELM:           http://elm.eu.org/\n"
            "  DOMINO:        http://mint.bio.uniroma2.it/domino/\n"
            "  PDB:           https://www.rcsb.org/downloads\n\n"
            "Provide one or more --<source> files. Only provided sources are used."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--queries", required=True, help="TSV/CSV with uniprot_canonical,start_aa_position,end_aa_position")

    # Source files (all optional; provide any subset)
    ap.add_argument("--intact", default=None, help="Path to IntAct export (uses AA cols 36,37)")
    ap.add_argument("--domino", default=None, help="Path to DOMINO export (uses AA cols 32,33)")
    ap.add_argument("--interactome3d", default=None, help="Path to Interactome3D export")
    ap.add_argument("--threedid", default=None, help="Path to 3DID export")
    ap.add_argument("--elm", default=None, help="Path to ELM export")
    ap.add_argument("--pdb", default=None, help="Path to a PDB-derived interface table")

    # Optional ID columns (1-based)
    ap.add_argument("--intact-protein-col", type=int, default=1, help="IntAct protein accession column")
    ap.add_argument("--domino-protein-col", type=int, default=1, help="DOMINO protein accession column")
    ap.add_argument("--elm-protein-col", type=int, default=1, help="ELM protein accession column")

    # Interactome3D: protein cols fixed per your note, AA cols configurable
    ap.add_argument("--i3d-protA-col", type=int, default=5)
    ap.add_argument("--i3d-protB-col", type=int, default=6)
    ap.add_argument("--i3d-aa-cols", default="11,12,18,19", help="A_start,A_end,B_start,B_end")

    # 3DID protein cols (must be provided correctly for your file)
    ap.add_argument("--threedid-protA-col", type=int, default=1)
    ap.add_argument("--threedid-protB-col", type=int, default=2)

    # PDB columns must be explicit because formats vary
    ap.add_argument("--pdb-protA-col", type=int, default=1)
    ap.add_argument("--pdb-protB-col", type=int, default=2)
    ap.add_argument("--pdb-aa-cols", default="3,4,5,6", help="A_start,A_end,B_start,B_end")

    ap.add_argument("--out-prefix", required=True, help="Output prefix, writes .query_summary.tsv and .overlap_long.tsv")
    args = ap.parse_args()

    # Ensure at least one source file is provided
    if not any([args.intact, args.domino, args.interactome3d, args.threedid, args.elm, args.pdb]):
        ap.error("You must provide at least one source file (e.g. --intact intact.tsv).")

    queries = parse_queries(args.queries)

    sources: Dict[str, List[InterfaceRec]] = {}

    # IntAct: AA 36,37
    if args.intact:
        sources["intact"] = parse_single_side_source(
            source="intact",
            path=args.intact,
            protein_col=args.intact_protein_col,
            aa_start_col=36,
            aa_end_col=37,
        )

    # DOMINO: AA 32,33
    if args.domino:
        sources["domino"] = parse_single_side_source(
            source="domino",
            path=args.domino,
            protein_col=args.domino_protein_col,
            aa_start_col=32,
            aa_end_col=33,
        )

    # 3DID: AA pairs 11,12 and 18,19
    if args.threedid:
        sources["threedid"] = parse_dual_partner_source(
            source="threedid",
            path=args.threedid,
            protA_col=args.threedid_protA_col,
            protB_col=args.threedid_protB_col,
            A_start_col=11,
            A_end_col=12,
            B_start_col=18,
            B_end_col=19,
        )

    # ELM: AA pairs 4,5 and 6,7
    # many ELM exports are single-protein; we still support dual-parse if protein B is absent (set same col by default)
    if args.elm:
        sources["elm"] = parse_dual_partner_source(
            source="elm",
            path=args.elm,
            protA_col=args.elm_protein_col,
            protB_col=args.elm_protein_col,
            A_start_col=4,
            A_end_col=5,
            B_start_col=6,
            B_end_col=7,
        )

    # Interactome3D: proteins 5,6 and AA by argument
    if args.interactome3d:
        a1, a2, b1, b2 = [int(x.strip()) for x in args.i3d_aa_cols.split(",")]
        sources["interactome3d"] = parse_dual_partner_source(
            source="interactome3d",
            path=args.interactome3d,
            protA_col=args.i3d_protA_col,
            protB_col=args.i3d_protB_col,
            A_start_col=a1,
            A_end_col=a2,
            B_start_col=b1,
            B_end_col=b2,
        )

    # PDB: fully configurable
    if args.pdb:
        p1, p2, p3, p4 = [int(x.strip()) for x in args.pdb_aa_cols.split(",")]
        sources["pdb"] = parse_dual_partner_source(
            source="pdb",
            path=args.pdb,
            protA_col=args.pdb_protA_col,
            protB_col=args.pdb_protB_col,
            A_start_col=p1,
            A_end_col=p2,
            B_start_col=p3,
            B_end_col=p4,
        )

    summary, longrows = annotate(queries, sources)

    write_tsv(f"{args.out_prefix}.query_summary.tsv", summary)
    write_tsv(f"{args.out_prefix}.overlap_long.tsv", longrows)

    # small run report
    print("Loaded records:")
    for k, v in sources.items():
        print(f"  {k}: {len(v)}")
    print(f"Queries: {len(queries)}")
    print(f"Wrote: {args.out_prefix}.query_summary.tsv")
    print(f"Wrote: {args.out_prefix}.overlap_long.tsv")


if __name__ == "__main__":
    main()