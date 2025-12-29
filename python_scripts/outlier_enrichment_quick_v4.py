#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-input workflow on node_impact_with_features.tsv

Inputs:
  node_impact_with_features.tsv  (columns like: gene, deltae, deltas, regulation, pfam_modules, ptms, disorder, complex, ...)

Outputs (under --outdir):
  <metric>_bin_stats.tsv
  <metric>_tail_tests.tsv
  <metric>_perm_bin_z.tsv
  (if --quantreg) <metric>_quantreg.tsv
  plots/  (bin curves, tail boxplots, perm Z, and quantreg lines)
  report.pdf  (selected summary plots)

Usage (example):
  python node_enrichment_qreg_single.py \
    --in node_impact_with_features.tsv \
    --outdir q_enrich_out \
    --metrics deltae deltas \
    --features regulation pfam_modules ptms disorder complex pathways domain_variety subcellular \
    --nbins 10 --tail 0.10 --n-perm 200 --quantreg
"""

import os, argparse, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages

# Optional: statsmodels for quantile regression
try:
    import statsmodels.formula.api as smf
except Exception:
    smf = None

plt.rcParams.update({"figure.dpi": 130, "axes.spines.right": False, "axes.spines.top": False})

# ---------------- helpers ----------------
def read_table(path: str) -> pd.DataFrame:
    # Try TSV, then CSV, then whitespace
    for sep in ["\t", ",", r"\s+"]:
        try:
            return pd.read_csv(path, sep=sep, engine="python")
        except Exception:
            continue
    raise ValueError(f"Could not parse {path} — check delimiter.")

def is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series)

def safe_sem(x: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    return np.nan if x.size <= 1 else x.std(ddof=1) / np.sqrt(x.size)

def quantile_bins(x: pd.Series, nbins: int) -> pd.Series:
    ranks = x.rank(method="average")
    return pd.qcut(ranks, nbins, labels=False, duplicates="drop")

def agg_by_bins(df: pd.DataFrame, bins: pd.Series, feature: str, numeric: bool) -> pd.DataFrame:
    out = []
    for b in sorted(bins.dropna().unique()):
        sub = df.loc[bins == b, feature]
        if numeric:
            mu, se = np.nanmean(sub), safe_sem(sub.to_numpy())
            n = sub.notna().sum()
            out.append((b, n, mu, se))
        else:
            counts = sub.value_counts(dropna=True)
            n = counts.sum()
            if n == 0:
                out.append((b, 0, np.nan, np.nan))
            else:
                frac = counts.iloc[0] / n
                out.append((b, n, float(frac), np.nan))
    return pd.DataFrame(out, columns=["bin", "n", "mean", "sem"])

def tail_test_numeric(top: pd.Series, bot: pd.Series):
    top, bot = top.dropna(), bot.dropna()
    if len(top) < 3 or len(bot) < 3:
        return np.nan, np.nan, np.nan, np.nan
    u, p = stats.mannwhitneyu(top, bot, alternative="two-sided")
    cliff = (2*u)/(len(top)*len(bot)) - 1
    diff = top.mean() - bot.mean()
    return p, cliff, diff, u

def permutation_bin_null(df, feature, metric_vals, nbins, n_perm, numeric):
    bins_true = quantile_bins(metric_vals, nbins)
    real = agg_by_bins(df, bins_true, feature, numeric)
    real_means = real["mean"].to_numpy()
    perm_means = np.zeros((n_perm, len(real_means)))
    for i in range(n_perm):
        shuffled = metric_vals.sample(frac=1.0, replace=False, random_state=np.random.randint(1e9)).reset_index(drop=True)
        shuffled.index = metric_vals.index
        b = quantile_bins(shuffled, nbins)
        pm = agg_by_bins(df, b, feature, numeric)["mean"].to_numpy()
        if len(pm) != len(real_means):
            tmp = np.full_like(real_means, np.nan, dtype=float)
            tmp[:len(pm)] = pm
            pm = tmp
        perm_means[i, :] = pm
    mu = np.nanmean(perm_means, axis=0)
    sd = np.nanstd(perm_means, axis=0, ddof=1)
    z = (real_means - mu) / sd
    return pd.DataFrame({"bin": np.arange(len(real_means)),
                         "real_mean": real_means, "perm_mean": mu, "perm_sd": sd, "z": z})

def save_bins_plot(ax, metric_name, feature, stats_df):
    ax.errorbar(stats_df["bin"], stats_df["mean"], yerr=stats_df["sem"], fmt="-o", capsize=3)
    ax.set_xlabel(f"{metric_name} quantile bin (0=lowest)")
    ax.set_ylabel(f"{feature} (mean ± SEM)")
    ax.set_title(f"{feature} across {metric_name} quantiles")

def save_tail_boxplot(ax, metric_name, feature, df, metric_vals, tail):
    n = len(metric_vals)
    k = max(1, int(np.floor(n * tail)))
    order = metric_vals.sort_values(kind="mergesort").index
    bot_ix = order[:k]
    mid_ix = order[k:-k] if (n - 2*k) > 0 else order[k:]
    top_ix = order[-k:]

    data = [df.loc[top_ix, feature].values,
            df.loc[mid_ix, feature].values,
            df.loc[bot_ix, feature].values]
    labels = ["Top tail", "Middle", "Bottom tail"]
    ax.boxplot(data, tick_labels=labels, showfliers=True)
    ax.set_title(f"{feature}: tails vs middle ({metric_name})")

def plot_perm_z(ax, metric_name, feature, zdf):
    ax.axhline(0, color="gray", lw=1)
    ax.plot(zdf["bin"], zdf["z"], marker="o")
    ax.set_xlabel(f"{metric_name} quantile bin")
    ax.set_ylabel("Z (real vs permuted mean)")
    ax.set_title(f"{feature}: per-bin Z ({metric_name})")

def run_quantreg(df, metric, features, outdir, quantiles=(0.25, 0.5, 0.75)):
    if smf is None:
        print("[warn] statsmodels not installed; skipping quantile regression.")
        return None
    qdir = os.path.join(outdir, metric, "quantreg")
    os.makedirs(qdir, exist_ok=True)
    rows = []
    for feat in features:
        if not is_numeric(df[feat]): 
            continue
        # drop rows missing either variable
        sub = df[[metric, feat]].dropna()
        if sub.shape[0] < 20:
            continue
        for q in quantiles:
            try:
                mod = smf.quantreg(f"{feat} ~ {metric}", sub)
                res = mod.fit(q=q)
                slope = res.params.get(metric, np.nan)
                p = res.pvalues.get(metric, np.nan)
                ci = res.conf_int().loc[metric].tolist() if metric in res.params.index else [np.nan, np.nan]
                rows.append({"feature": feat, "quantile": q, "slope": slope, "p": p, "ci_low": ci[0], "ci_high": ci[1]})
                # quick line overlay plot
                fig, ax = plt.subplots(figsize=(5,4))
                ax.scatter(sub[metric], sub[feat], alpha=0.35, s=10)
                xpred = np.linspace(sub[metric].min(), sub[metric].max(), 120)
                ypred = res.params["Intercept"] + slope * xpred
                ax.plot(xpred, ypred, label=f"q={q:.2f}", lw=2)
                ax.set_xlabel(metric); ax.set_ylabel(feat); ax.legend()
                fig.tight_layout()
                fig.savefig(os.path.join(qdir, f"{feat}_q{int(q*100)}.pdf"), dpi=140)
                plt.close(fig)
            except Exception as e:
                print(f"[quantreg fail] {feat} q={q}: {e}")
    if rows:
        out = pd.DataFrame(rows).sort_values(["p","feature","quantile"])
        out.to_csv(os.path.join(outdir, f"{metric}_quantreg.tsv"), sep="\t", index=False)
        return out
    return None

# ---------------- main analysis ----------------
def analyze_metric(df, metric, features, nbins, tail, n_perm, outdir, pdf=None, do_quantreg=False):
    metric_vals = pd.to_numeric(df[metric], errors="coerce")
    keep = metric_vals.notna()
    df = df.loc[keep].copy()
    metric_vals = metric_vals.loc[keep]

    bins_all, tails_all, perms_all = [], [], []
    metric_dir = os.path.join(outdir, metric)
    os.makedirs(metric_dir, exist_ok=True)
    plots_dir = os.path.join(metric_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # per-feature analyses
    for feat in features:
        if feat == "gene":  # skip ID column if present
            continue
        series = df[feat]
        numeric = is_numeric(series)

        # 1) Quantile-bin curve
        bins = quantile_bins(metric_vals, nbins)
        stats_df = agg_by_bins(df, bins, feat, numeric)
        stats_df["feature"] = feat
        bins_all.append(stats_df)

        # plot bins
        fig, ax = plt.subplots(figsize=(7,4))
        save_bins_plot(ax, metric, feat, stats_df)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, f"{metric}_feature_bins_{feat}.pdf"), dpi=140)
        if pdf: pdf.savefig(fig)
        plt.close(fig)

        # 2) Tail test (top vs bottom)
        n = len(metric_vals)
        k = max(1, int(np.floor(n * tail)))
        order = metric_vals.sort_values(kind="mergesort").index
        bot_ix = order[:k]
        top_ix = order[-k:]
        if numeric:
            p, cliff, diff, _ = tail_test_numeric(series.loc[top_ix], series.loc[bot_ix])
            tails_all.append({"feature": feat, "type": "numeric",
                              "p_mannwhitney": p, "cliffs_delta": cliff,
                              "mean_top_minus_bottom": diff,
                              "n_top": series.loc[top_ix].notna().sum(),
                              "n_bottom": series.loc[bot_ix].notna().sum()})

        # tail boxplot (top/mid/bottom)
        fig, ax = plt.subplots(figsize=(7,4))
        save_tail_boxplot(ax, metric, feat, df, metric_vals, tail)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, f"{metric}_feature_tail_box_{feat}.pdf"), dpi=140)
        if pdf: pdf.savefig(fig)
        plt.close(fig)

        # 3) Permutation null for bin curves
        zdf = permutation_bin_null(df, feat, metric_vals, nbins, n_perm, numeric)
        zdf["feature"] = feat
        perms_all.append(zdf)

        # plot Z per bin
        fig, ax = plt.subplots(figsize=(7,4))
        plot_perm_z(ax, metric, feat, zdf)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, f"{metric}_feature_perm_z_{feat}.pdf"), dpi=140)
        if pdf: pdf.savefig(fig)
        plt.close(fig)

    # write TSVs
    pd.concat(bins_all, ignore_index=True).to_csv(os.path.join(outdir, f"{metric}_bin_stats.tsv"), sep="\t", index=False)
    pd.DataFrame(tails_all).to_csv(os.path.join(outdir, f"{metric}_tail_tests.tsv"), sep="\t", index=False)
    pd.concat(perms_all, ignore_index=True).to_csv(os.path.join(outdir, f"{metric}_perm_bin_z.tsv"), sep="\t", index=False)

    # 4) optional quantile regression
    if do_quantreg:
        _ = run_quantreg(df, metric, features, outdir)

def main():
    ap = argparse.ArgumentParser(description="Single-input enrichment + optional quantile regression")
    ap.add_argument("--in", dest="inp", required=True, help="node_impact_with_features.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--metrics", nargs="+", default=["deltae","deltas"], help="metric columns (default: deltae deltas)")
    ap.add_argument("--features", nargs="*", default=[],
                    help="feature columns to test; if empty, autodetect numerics except metrics/ID")
    ap.add_argument("--nbins", type=int, default=10)
    ap.add_argument("--tail", type=float, default=0.10)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--quantreg", action="store_true", help="run quantile regression (q=0.25,0.5,0.75)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = read_table(args.inp)
    # Normalise headers
    df.columns = [c.strip().lower() for c in df.columns]

    # Ensure metric columns exist (case-insensitive handled by lowercasing)
    missing = [m for m in args.metrics if m not in df.columns]
    if missing:
        sys.exit(f"[error] Missing metric columns: {missing}. Found: {list(df.columns)}")

    # Autodetect features if not given
    feats = args.features
    if not feats:
        blocked = set(args.metrics) | {"gene", "id"}
        feats = [c for c in df.columns if c not in blocked]
    # Coerce obvious numeric features
    for c in feats:
        # leave categorical as-is; numeric detection happens per-feature
        if c in df.columns and df[c].dtype == object:
            # try numeric conversion without clobbering non-numerics
            maybe = pd.to_numeric(df[c], errors="coerce")
            # if we get a decent amount of numerics, keep them
            if maybe.notna().sum() >= 0.3 * len(maybe):
                df[c] = maybe

    # Summary
    with open(os.path.join(args.outdir, "run_summary.txt"), "w") as fh:
        fh.write(f"Rows: {len(df)}\n")
        for m in args.metrics:
            fh.write(f"{m}: non-null {df[m].notna().sum()}\n")
        for f in feats:
            fh.write(f"{f}: type={'numeric' if is_numeric(df[f]) else 'categorical'}, non-null={df[f].notna().sum()}\n")

    # PDF report
    report_path = os.path.join(args.outdir, "report.pdf")
    with PdfPages(report_path) as pdf:
        for m in args.metrics:
            analyze_metric(df, m, feats, args.nbins, args.tail, args.n_perm, args.outdir, pdf=pdf, do_quantreg=args.quantreg)

    print("[done] Outputs written to:", args.outdir)
    print("[done] PDF report:", report_path)

if __name__ == "__main__":
    main()
