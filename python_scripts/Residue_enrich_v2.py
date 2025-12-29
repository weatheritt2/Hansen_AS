#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
microexon_residue_enrichment_v2.py

Compare amino-acid composition between MICRO and LONG exons in:
  - Surrounding flanks (upstream/downstream windows)
  - Exon body
Optionally restrict to a structure-filtered subset (e.g. helix-like).

Outputs:
  1) Residue-level long-format table (log2 enrichment + p + q + -log10p)
  2) Residue-level wide table (counts + frequencies for MICRO/LONG)
  3) Feature-group table (charged/helix-breaking/etc) in long format

"""

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA20)


# ------------------------ stats helpers ------------------------

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR; returns q-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def safe_log2_ratio(a: float, b: float, pseudocount: float = 0.5) -> float:
    """log2((a+p)/(b+p)) with pseudocount to avoid div-by-zero."""
    return math.log((a + pseudocount) / (b + pseudocount), 2)


def neglog10(p: float, floor: float = 1e-300) -> float:
    p = max(float(p), floor)
    return -math.log10(p)


# ------------------------ parsing ------------------------

@dataclass
class Record:
    klass: str            # "MICRO" or "LONG"
    up: str               # upstream flank
    exon: str             # exon sequence
    down: str             # downstream flank
    struct: Optional[str] # structure string/token (optional)


def infer_class(line: str) -> Optional[str]:
    """Detect MICRO/LONG from line content."""
    u = line.upper()
    if "MICRO" in u:
        return "MICRO"
    if "LONG" in u:
        return "LONG"
    return None


def clean_seq(seq: str) -> str:
    """Uppercase and keep only canonical AA letters; stop on delimiters | or _."""
    out = []
    for ch in seq.upper():
        if ch in ("|", "_"):
            break
        if ch in AA_SET:
            out.append(ch)
        # ignore other characters
    return "".join(out)


def parse_line(line: str) -> Optional[Record]:
    """
    Expected minimal tokenization:
      tokens[0]=upstream, tokens[1]=exon, tokens[2]=downstream, tokens[3]=structure (optional)
    Additional tokens are ignored.
    """
    toks = line.strip().split()
    if len(toks) < 3:
        return None
    klass = infer_class(line)
    if klass is None:
        return None

    up = clean_seq(toks[0])
    exon = clean_seq(toks[1])
    down = clean_seq(toks[2])
    struct = toks[3] if len(toks) >= 4 else None
    return Record(klass=klass, up=up, exon=exon, down=down, struct=struct)


def load_records(path: Path) -> List[Record]:
    recs = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            r = parse_line(line)
            if r is not None:
                recs.append(r)
    return recs


# ------------------------ counting ------------------------

def init_counts() -> Dict[str, int]:
    return {aa: 0 for aa in AA20}


def add_counts(counts: Dict[str, int], seq: str) -> None:
    for aa in seq:
        if aa in counts:
            counts[aa] += 1


def summarize_counts(micro: Dict[str, int], long: Dict[str, int]) -> pd.DataFrame:
    """Return a residue table with counts and frequencies for MICRO vs LONG."""
    tm = sum(micro.values())
    tl = sum(long.values())
    rows = []
    for aa in AA20:
        cm = micro.get(aa, 0)
        cl = long.get(aa, 0)
        pm = (cm / tm) if tm else np.nan
        pl = (cl / tl) if tl else np.nan
        rows.append((aa, cm, cl, tm, tl, pm, pl))
    return pd.DataFrame(
        rows,
        columns=["residue", "micro_count", "long_count", "micro_total", "long_total", "micro_freq", "long_freq"],
    )


def fisher_residue_table(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Add Fisher p-values and log2 enrichment for each residue."""
    pvals = []
    log2fc = []
    for _, r in df_counts.iterrows():
        cm = int(r["micro_count"])
        cl = int(r["long_count"])
        tm = int(r["micro_total"])
        tl = int(r["long_total"])
        # 2x2: residue vs other residues
        # [[cm, tm-cm], [cl, tl-cl]]
        if tm <= 0 or tl <= 0:
            p = np.nan
        else:
            p = fisher_exact([[cm, tm - cm], [cl, tl - cl]], alternative="two-sided")[1]
        pvals.append(p)
        log2fc.append(safe_log2_ratio(r["micro_freq"], r["long_freq"]))
    out = df_counts.copy()
    out["p"] = pvals
    out["log2_enrich_micro_vs_long"] = log2fc
    out["neglog10p"] = out["p"].apply(lambda x: neglog10(x) if pd.notna(x) else np.nan)
    return out


def apply_fdr_by_group(df: pd.DataFrame, group_cols: List[str], p_col: str = "p") -> pd.DataFrame:
    """BH-FDR within each group (e.g. region×subset)."""
    df = df.copy()
    df["q_fdr"] = np.nan
    for _, sub in df.groupby(group_cols, dropna=False, sort=False):
        idx = sub.index
        p = sub[p_col].to_numpy()
        ok = np.isfinite(p)
        q = np.full_like(p, np.nan, dtype=float)
        if ok.sum() > 0:
            q[ok] = bh_fdr(p[ok])
        df.loc[idx, "q_fdr"] = q
    return df


# ------------------------ feature groups ------------------------

FEATURE_GROUPS = {
    "Charged": list("EKRD"),
    "Helix_Breaking": list("PG"),
    "Helix_Forming": list("LIFEYWM"),
    "Aliphatic": list("LIV"),
    "Aromatic": list("FYWH"),
}

def group_enrichment(group_aas: List[str], micro_counts: Dict[str, int], long_counts: Dict[str, int]) -> Tuple[float, float, int, int, int, int]:
    """Return log2 enrichment, p-value, and raw totals."""
    tm = sum(micro_counts.values())
    tl = sum(long_counts.values())
    cm = sum(micro_counts.get(a, 0) for a in group_aas)
    cl = sum(long_counts.get(a, 0) for a in group_aas)
    if tm <= 0 or tl <= 0:
        return np.nan, np.nan, cm, tm, cl, tl
    p = fisher_exact([[cm, tm - cm], [cl, tl - cl]], alternative="two-sided")[1]
    pm = cm / tm
    pl = cl / tl
    return safe_log2_ratio(pm, pl), p, cm, tm, cl, tl


# ------------------------ main workflow ------------------------

def filter_records(
    recs: List[Record],
    min_long_exon_len: int,
    skip_exon_start_m: bool,
    require_struct_regex: Optional[str],
) -> Tuple[List[Record], List[Record]]:
    micro = []
    long = []
    rx = re.compile(require_struct_regex) if require_struct_regex else None

    for r in recs:
        if skip_exon_start_m and r.exon.startswith("M"):
            continue

        if r.klass == "LONG" and len(r.exon) < min_long_exon_len:
            continue

        # For the structure-filtered subset, we’ll decide later; here we only filter
        # by schema-level requirements.
        if r.klass == "MICRO":
            micro.append(r)
        elif r.klass == "LONG":
            long.append(r)

    return micro, long


def count_regions(records: List[Record], flank: int) -> Dict[str, Dict[str, int]]:
    """
    Count residues in regions:
      - surround: up[-flank:] + down[:flank]
      - exon: exon
      - first: up[-flank:]
      - last: down[:flank]
    """
    out = {
        "surround": init_counts(),
        "exon": init_counts(),
        "first": init_counts(),
        "last": init_counts(),
    }
    for r in records:
        up = r.up[-flank:] if flank > 0 else ""
        down = r.down[:flank] if flank > 0 else ""
        add_counts(out["first"], up)
        add_counts(out["last"], down)
        add_counts(out["surround"], up + down)
        add_counts(out["exon"], r.exon)
    return out


def subset_by_structure(records: List[Record], require_struct_regex: Optional[str]) -> List[Record]:
    if not require_struct_regex:
        return records
    rx = re.compile(require_struct_regex)
    keep = []
    for r in records:
        if r.struct is None:
            continue
        if rx.search(r.struct):
            keep.append(r)
    return keep


def main():
    ap = argparse.ArgumentParser(description="Residue/feature enrichment: MICRO vs LONG (with optional structure subset + FDR).")
    ap.add_argument("--in", dest="inp", required=True, help="microexon_seqs.txt-like input file")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--prefix", default="microexon_enrich", help="Output file prefix")
    ap.add_argument("--flank", type=int, default=50, help="Flank window size for surround/first/last regions")
    ap.add_argument("--min-long-exon-len", type=int, default=30, help="Skip LONG records with exon length < this")
    ap.add_argument("--skip-exon-start-m", action="store_true", help="Skip records whose exon starts with 'M'")
    ap.add_argument("--struct-regex", default=None, help="Regex applied to structure token to define 'subset' (e.g. 'H|G')")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    recs = load_records(Path(args.inp))
    if not recs:
        raise SystemExit(f"[error] No parseable records found in {args.inp}")

    micro_all, long_all = filter_records(
        recs,
        min_long_exon_len=args.min_long_exon_len,
        skip_exon_start_m=args.skip_exon_start_m,
        require_struct_regex=args.struct_regex,
    )

    # Define two subsets:
    #   - ALL: all records after basic filters
    #   - SUBSET: those matching --struct-regex (if provided)
    subsets = [("ALL", micro_all, long_all)]
    if args.struct_regex:
        micro_sub = subset_by_structure(micro_all, args.struct_regex)
        long_sub = subset_by_structure(long_all, args.struct_regex)
        subsets.append(("SUBSET", micro_sub, long_sub))

    # Collect residue-level outputs in both wide and tidy formats
    tidy_rows = []
    wide_rows = []

    for subset_name, micro_recs, long_recs in subsets:
        micro_regions = count_regions(micro_recs, flank=args.flank)
        long_regions = count_regions(long_recs, flank=args.flank)

        for region in ["surround", "exon", "first", "last"]:
            df_counts = summarize_counts(micro_regions[region], long_regions[region])
            df_stats = fisher_residue_table(df_counts)

            df_stats["subset"] = subset_name
            df_stats["region"] = region

            wide_rows.append(df_stats)

            # tidy (long) version for plotting
            tidy = df_stats[[
                "residue", "subset", "region",
                "micro_count", "long_count", "micro_total", "long_total",
                "micro_freq", "long_freq",
                "log2_enrich_micro_vs_long", "p", "neglog10p"
            ]].copy()
            tidy_rows.append(tidy)

    wide = pd.concat(wide_rows, ignore_index=True)
    tidy = pd.concat(tidy_rows, ignore_index=True)

    # FDR correction within each subset×region
    tidy = apply_fdr_by_group(tidy, ["subset", "region"], p_col="p")
    wide = apply_fdr_by_group(wide, ["subset", "region"], p_col="p")
    wide["q_fdr"] = tidy["q_fdr"].values  # align (same ordering)

    # Write residue outputs
    wide_path = outdir / f"{args.prefix}_residue_wide.tsv"
    tidy_path = outdir / f"{args.prefix}_residue_tidy.tsv"
    wide.to_csv(wide_path, sep="\t", index=False)
    tidy.to_csv(tidy_path, sep="\t", index=False)

    # Feature-group enrichments
    feat_rows = []
    for subset_name, micro_recs, long_recs in subsets:
        micro_regions = count_regions(micro_recs, flank=args.flank)
        long_regions = count_regions(long_recs, flank=args.flank)
        for region in ["surround", "exon"]:
            for feat_name, aas in FEATURE_GROUPS.items():
                log2fc, p, cm, tm, cl, tl = group_enrichment(aas, micro_regions[region], long_regions[region])
                feat_rows.append({
                    "feature": feat_name,
                    "subset": subset_name,
                    "region": region,
                    "micro_count": cm, "micro_total": tm,
                    "long_count": cl, "long_total": tl,
                    "log2_enrich_micro_vs_long": log2fc,
                    "p": p,
                    "neglog10p": neglog10(p) if pd.notna(p) else np.nan,
                })

    feat_df = pd.DataFrame(feat_rows)
    if not feat_df.empty:
        feat_df = apply_fdr_by_group(feat_df, ["subset", "region"], p_col="p")
    feat_path = outdir / f"{args.prefix}_features.tsv"
    feat_df.to_csv(feat_path, sep="\t", index=False)

    # Run summary
    summary_path = outdir / f"{args.prefix}_run_summary.txt"
    with summary_path.open("w") as fh:
        fh.write(f"Input: {args.inp}\n")
        fh.write(f"Total parsed records: {len(recs)}\n")
        fh.write(f"MICRO kept (ALL): {len(micro_all)}\n")
        fh.write(f"LONG  kept (ALL): {len(long_all)}\n")
        if args.struct_regex:
            fh.write(f"Structure subset regex: {args.struct_regex}\n")
            fh.write(f"MICRO kept (SUBSET): {len(subsets[1][1])}\n")
            fh.write(f"LONG  kept (SUBSET): {len(subsets[1][2])}\n")
        fh.write(f"Flank size: {args.flank}\n")
        fh.write(f"Min LONG exon len: {args.min_long_exon_len}\n")
        fh.write(f"Skip exon starting M: {args.skip_exon_start_m}\n")
        fh.write("\nOutputs:\n")
        fh.write(f"- {wide_path}\n")
        fh.write(f"- {tidy_path}\n")
        fh.write(f"- {feat_path}\n")

    print("[OK] Wrote:")
    print("  ", wide_path)
    print("  ", tidy_path)
    print("  ", feat_path)
    print("  ", summary_path)


if __name__ == "__main__":
    main()