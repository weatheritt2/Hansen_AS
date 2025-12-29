#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Percolation with paired Δ, degree-matched nulls, κ(f), and inter-community enrichment
(with replicate-level null distributions + per-f p-values).

Outputs (prefix = --out-prefix):
  - <prefix>_<metric>_curves.tsv                 # AS curve, null mean, null CI, paired Δ mean/CI
  - <prefix>_<metric>_auc.tsv                    # AUC summary (AS vs null), p-value, z
  - <prefix>_kappa_curves.tsv                    # κ(f) for AS and null
  - <prefix>_intercomm.tsv                       # Inter-community edges removed: AS, null mean/CI, p-values per f
  - <prefix>_<metric>.png                        # Plot S/E with null band and Δ CI
  - <prefix>_<metric>_delta.png                  # Paired Δ plot

Example:
  python percolation_auc_plus_v2.py \
    --edges interface_interactome.tab \
    --sep tab \
    --col-a A --col-b B --col-as as_prob \
    --metric both \
    --n-nulls 1000 \
    --n-jobs 32 \
    --null-mode degree \
    --deg-bins 0,1,2,3,5,8,13,21,34,55,89 \
    --out-prefix run_plus_v2
"""

import argparse, os, sys, random
import numpy as np
import pandas as pd
import networkx as nx
from multiprocessing import Pool, cpu_count
from functools import partial
from networkx.algorithms.community import greedy_modularity_communities

# ------------------------ Robust reading ------------------------

def _resolve_sep(sep_flag: str):
    if sep_flag is None or sep_flag == "auto":
        return None
    s = sep_flag.lower()
    if s in ["tab", "\\t"]: return "\t"
    if s in ["csv", ","]:   return ","
    if s in ["space", "ws", "whitespace"]: return "space"
    return None

def _read_sniff(path, comment=None):
    return pd.read_csv(path, sep=None, engine="python", comment=comment)

def _read_explicit(path, sep, comment=None):
    return pd.read_csv(path, sep=sep, engine="python", comment=comment)

def _read_whitespace(path, comment=None):
    return pd.read_csv(path, delim_whitespace=True, engine="python", comment=comment)

def read_edges(path, col_a, col_b, col_as, sep_flag="auto", comment=None):
    df = None
    sep = _resolve_sep(sep_flag)
    if sep == "space":
        try: df = _read_whitespace(path, comment=comment)
        except: pass
    if df is None and sep is None:
        try: df = _read_sniff(path, comment=comment)
        except: pass
    if df is None and sep not in (None, "space"):
        try: df = _read_explicit(path, sep=sep, comment=comment)
        except: pass
    if df is None:
        df = pd.read_csv(path, sep=r"[,\t ]+", engine="python", comment=comment)

    colmap_lower = {c.lower(): c for c in df.columns}
    def _find(name):
        if name in df.columns: return name
        key = name.lower()
        if key in colmap_lower: return colmap_lower[key]
        raise ValueError(f"Required column '{name}' not found. Available: {list(df.columns)}")

    A = _find(col_a)
    B = _find(col_b)
    AS = _find(col_as)

    df[A] = df[A].astype(str)
    df[B] = df[B].astype(str)
    df[AS] = pd.to_numeric(df[AS], errors="coerce").fillna(0.0)

    before = len(df)
    df = df.dropna(subset=[A, B])
    if len(df) < before:
        print(f"[read] Dropped {before - len(df)} rows with missing {A}/{B}", file=sys.stderr)

    return df, A, B, AS

# ------------------------ Graph & metrics ------------------------

def build_graph_from_df(df, A, B):
    G = nx.Graph()
    G.add_edges_from(zip(df[A], df[B]))
    return G

def S_fraction(G):
    if G.number_of_nodes() == 0:
        return 0.0
    comps = [len(c) for c in nx.connected_components(G)]
    return (max(comps) / G.number_of_nodes()) if comps else 0.0

def kappa_components(G):
    return nx.number_connected_components(G)

def approx_global_efficiency(G, sources=None, max_sources=200, seed=0):
    n = G.number_of_nodes()
    if n <= 1: return 0.0
    rng = random.Random(seed)
    nodes = list(G.nodes())
    if sources is None:
        sources = nodes if len(nodes) <= max_sources else rng.sample(nodes, k=max_sources)
    total = 0.0; count = 0
    for s in sources:
        dists = nx.single_source_shortest_path_length(G, s)
        for v, d in dists.items():
            if v != s and d > 0:
                total += 1.0 / d
                count += 1
    return (total / count) if count else 0.0

def auc_trapz(y, x):
    return float(np.trapz(np.asarray(y, float), np.asarray(x, float)))

# ------------------------ Communities ------------------------

def community_labels(G):
    comms = list(greedy_modularity_communities(G))
    label = {}
    for i, cset in enumerate(comms):
        for u in cset: label[u] = i
    return label

def intercomm_mask(edge_pairs, labels):
    # boolean array: True if edge is inter-community
    mask = np.array([labels.get(u, -1) != labels.get(v, -1) for (u, v) in edge_pairs], dtype=bool)
    return mask

# ------------------------ Null constructions ------------------------

def degree_bins(G, edge_pairs, bins):
    deg = dict(G.degree())
    mins = [min(deg[u], deg[v]) for u, v in edge_pairs]
    idx = np.digitize(mins, bins, right=False) - 1
    bucket = {}
    for i, bi in enumerate(idx):
        bucket.setdefault(int(bi), []).append(i)
    return bucket  # {bin_id: [edge_idx,...]}

def null_global_worker(seed, edge_pairs, k_per_f, metric, sources_k, mask_inter):
    rng = random.Random(seed)
    m = len(edge_pairs)
    order = list(range(m))
    rng.shuffle(order)

    G = nx.Graph(); G.add_edges_from(edge_pairs)
    out_metric = []; out_kappa = []; out_inter = []

    for k in k_per_f:
        k = min(k, m)
        rm_idx = order[:k]
        to_remove = [edge_pairs[i] for i in rm_idx]
        H = G.copy()
        H.remove_edges_from(to_remove)
        if metric == "S":
            out_metric.append(S_fraction(H))
        else:
            out_metric.append(approx_global_efficiency(H, max_sources=sources_k, seed=seed))
        out_kappa.append(kappa_components(H))
        out_inter.append(int(np.sum(mask_inter[rm_idx])))
    return np.asarray(out_metric), np.asarray(out_kappa), np.asarray(out_inter)

def null_degree_worker(seed, edge_pairs, k_per_f, buckets, metric, sources_k, mask_inter):
    rng = random.Random(seed)
    # fixed random order inside each bucket for this replicate
    bucket_order = {b: rng.sample(idxs, k=len(idxs)) for b, idxs in buckets.items()}
    bucket_sizes = {b: len(idxs) for b, idxs in buckets.items()}
    total_m = len(edge_pairs)

    G = nx.Graph(); G.add_edges_from(edge_pairs)
    out_metric = []; out_kappa = []; out_inter = []

    for k in k_per_f:
        k = min(k, total_m)
        # proportional target per bucket
        to_take = {b: int(round(k * (bucket_sizes[b]/total_m))) for b in buckets}
        # adjust rounding to hit exactly k
        diff = k - sum(to_take.values())
        if diff != 0:
            order_b = sorted(bucket_sizes, key=lambda b: bucket_sizes[b], reverse=True)
            i = 0
            while diff != 0 and i < len(order_b):
                to_take[order_b[i]] += 1 if diff > 0 else -1
                diff += -1 if diff > 0 else 1
                i = (i+1) % len(order_b)

        rm_idx = []
        for b, need in to_take.items():
            need = max(0, min(need, bucket_sizes[b]))
            rm_idx.extend(bucket_order[b][:need])

        H = G.copy()
        H.remove_edges_from([edge_pairs[i] for i in rm_idx])
        if metric == "S":
            out_metric.append(S_fraction(H))
        else:
            out_metric.append(approx_global_efficiency(H, max_sources=sources_k, seed=seed))
        out_kappa.append(kappa_components(H))
        out_inter.append(int(np.sum(mask_inter[rm_idx])))

    return np.asarray(out_metric), np.asarray(out_kappa), np.asarray(out_inter)

# ------------------------ AS curve & κ ------------------------

def as_curve_and_kappa(G, edge_pairs, as_order_idx, k_per_f, metric, sources_k, mask_inter):
    out_metric = []; out_kappa = []; out_inter = []
    for k in k_per_f:
        k = max(0, min(k, len(as_order_idx)))
        rm_idx = as_order_idx[:k]
        H = G.copy()
        H.remove_edges_from([edge_pairs[i] for i in rm_idx])
        if metric == "S":
            out_metric.append(S_fraction(H))
        else:
            out_metric.append(approx_global_efficiency(H, max_sources=sources_k, seed=0))
        out_kappa.append(kappa_components(H))
        out_inter.append(int(np.sum(mask_inter[rm_idx])))
    return np.asarray(out_metric), np.asarray(out_kappa), np.asarray(out_inter)

# ------------------------ Main ------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Percolation with paired Δ, degree-matched nulls, κ(f), and inter-community enrichment (with CI and p-values).")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--sep", default="auto", help='Delimiter: auto/tab/csv/space')
    ap.add_argument("--comment", default=None)
    ap.add_argument("--col-a", default="A")
    ap.add_argument("--col-b", default="B")
    ap.add_argument("--col-as", default="as_prob", help="AS edge priority; removed DESC")
    ap.add_argument("--metric", choices=["S","E","both"], default="S")
    ap.add_argument("--fractions", default="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    ap.add_argument("--n-nulls", type=int, default=1000)
    ap.add_argument("--null-mode", choices=["global","degree"], default="global")
    ap.add_argument("--deg-bins", default="0,1,2,3,5,8,13,21,34,55,89",
                    help="Comma list of min-degree bin edges for degree-matched nulls")
    ap.add_argument("--n-jobs", type=int, default=max(1, cpu_count()//2))
    ap.add_argument("--sources-k", type=int, default=200, help="E approximation BFS sources")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", default="percolation_plus_v2")
    return ap.parse_args()

def main():
    args = parse_args()
    df, A, B, AS = read_edges(args.edges, args.col_a, args.col_b, args.col_as, sep_flag=args.sep, comment=args.comment)

    G = build_graph_from_df(df, A, B)
    edge_pairs = list(zip(df[A].tolist(), df[B].tolist()))
    m = len(edge_pairs)

    # AS order: descending by AS score (stable)
    order = np.argsort(-df[AS].to_numpy(), kind="mergesort")
    as_order_idx = order.tolist()

    # Fractions -> exact removals
    f_grid = np.array([float(x) for x in args.fractions.split(",")], dtype=float)
    f_grid = np.clip(f_grid, 0.0, 1.0)
    k_per_f = [int(round(f * len(as_order_idx))) for f in f_grid]

    # Community labels and inter-community mask (fixed on base graph)
    labels = community_labels(G)
    mask_inter = intercomm_mask(edge_pairs, labels)

    # Degree buckets for matched nulls
    if args.null_mode == "degree":
        bins = [int(b) for b in args.deg_bins.split(",")]
        buckets = degree_bins(G, edge_pairs, bins)
    else:
        buckets = None

    # Metrics to run
    metrics = ["S","E"] if args.metric == "both" else [args.metric]

    for metric in metrics:
        # AS curves and inter-community counts
        as_curve, as_kappa, inter_as = as_curve_and_kappa(
            G, edge_pairs, as_order_idx, k_per_f, metric, args.sources_k, mask_inter
        )

        # Null workers
        seeds = [args.seed + i + 1 for i in range(args.n_nulls)]
        if args.null_mode == "global":
            worker = partial(null_global_worker,
                             edge_pairs=edge_pairs,
                             k_per_f=k_per_f,
                             metric=metric,
                             sources_k=args.sources_k,
                             mask_inter=mask_inter)
        else:
            worker = partial(null_degree_worker,
                             edge_pairs=edge_pairs,
                             k_per_f=k_per_f,
                             buckets=buckets,
                             metric=metric,
                             sources_k=args.sources_k,
                             mask_inter=mask_inter)

        with Pool(processes=args.n_jobs) as pool:
            res = pool.map(worker, seeds)

        null_curves = np.vstack([r[0] for r in res])   # (R, |f|)
        null_kappa  = np.vstack([r[1] for r in res])   # (R, |f|)
        null_inter  = np.vstack([r[2] for r in res])   # (R, |f|)

        # Paired Δ = AS - null_r
        delta_curves = as_curve[None, :] - null_curves
        delta_mean = delta_curves.mean(axis=0)
        delta_lo   = np.quantile(delta_curves, 0.025, axis=0)
        delta_hi   = np.quantile(delta_curves, 0.975, axis=0)

        # κ(f) bands
        kappa_null_mean = null_kappa.mean(axis=0)
        kappa_null_lo   = np.quantile(null_kappa, 0.025, axis=0)
        kappa_null_hi   = np.quantile(null_kappa, 0.975, axis=0)

        # Inter-community bands + per-f p-values (two-sided by extremeness wrt null mean)
        inter_null_mean = null_inter.mean(axis=0)
        inter_null_lo   = np.quantile(null_inter, 0.025, axis=0)
        inter_null_hi   = np.quantile(null_inter, 0.975, axis=0)
        pvals = []
        for j in range(len(f_grid)):
            col = null_inter[:, j]
            mu = np.mean(col)
            # two-sided empirical p: extremeness around mu
            if inter_as[j] <= mu:
                p = (np.sum(col <= inter_as[j]) + 1.0) / (len(col) + 1.0)
            else:
                p = (np.sum(col >= inter_as[j]) + 1.0) / (len(col) + 1.0)
            p = min(1.0, 2.0 * p)  # two-sided
            pvals.append(p)
        pvals = np.array(pvals, dtype=float)

        # Null mean/CI for metric
        null_mean = null_curves.mean(axis=0)
        null_lo   = np.quantile(null_curves, 0.025, axis=0)
        null_hi   = np.quantile(null_curves, 0.975, axis=0)

        # AUC test (one-sided: AS more disruptive => smaller AUC)
        auc_as = auc_trapz(as_curve, f_grid)
        auc_nulls = np.array([auc_trapz(cur, f_grid) for cur in null_curves], dtype=float)
        mu_auc = float(np.mean(auc_nulls))
        sd_auc = float(np.std(auc_nulls, ddof=1)) if len(auc_nulls) > 1 else float("nan")
        z_auc  = (auc_as - mu_auc) / sd_auc if (sd_auc and sd_auc > 0) else float("nan")
        p_auc  = (np.sum(auc_nulls <= auc_as) + 1.0) / (len(auc_nulls) + 1.0)

        # ---- Save tables ----
        curves_path = f"{args.out_prefix}_{metric}_curves.tsv"
        pd.DataFrame({
            "f": f_grid,
            "AS": as_curve,
            "null_mean": null_mean,
            "null_lo": null_lo,
            "null_hi": null_hi,
            "Delta_mean": delta_mean,
            "Delta_lo": delta_lo,
            "Delta_hi": delta_hi
        }).to_csv(curves_path, sep="\t", index=False)

        auc_path = f"{args.out_prefix}_{metric}_auc.tsv"
        pd.DataFrame([{
            "metric": metric,
            "auc_as": auc_as,
            "auc_null_mean": mu_auc,
            "auc_null_sd": sd_auc,
            "n_nulls": len(auc_nulls),
            "p_empirical_lower": p_auc,
            "z_vs_null": z_auc,
            "null_mode": args.null_mode
        }]).to_csv(auc_path, sep="\t", index=False)

        kappa_path = f"{args.out_prefix}_kappa_curves.tsv"
        pd.DataFrame({
            "f": f_grid,
            "kappa_as": as_kappa,
            "kappa_null_mean": kappa_null_mean,
            "kappa_null_lo": kappa_null_lo,
            "kappa_null_hi": kappa_null_hi
        }).to_csv(kappa_path, sep="\t", index=False)

        inter_path = f"{args.out_prefix}_intercomm.tsv"
        pd.DataFrame({
            "f": f_grid,
            "inter_removed_AS": inter_as,
            "inter_removed_null_mean": inter_null_mean,
            "inter_removed_null_lo": inter_null_lo,
            "inter_removed_null_hi": inter_null_hi,
            "p_empirical_two_sided": pvals
        }).to_csv(inter_path, sep="\t", index=False)

        # ---- Plots ----
        try:
            import matplotlib.pyplot as plt
            # S/E curves
            plt.figure(figsize=(8,5))
            plt.fill_between(f_grid, null_lo, null_hi, alpha=0.2, label="Null 95% CI")
            plt.plot(f_grid, null_mean, color="C0", lw=2, label="Null mean")
            plt.plot(f_grid, as_curve, color="C1", lw=2, marker="o", label="AS")
            ttl = f"{metric}(f): AUC p={p_auc:.3g}, z={z_auc:.2f}, null={args.null_mode}"
            plt.title(ttl)
            plt.xlabel("Fraction of AS edges removed (f)")
            plt.ylabel("S(f)" if metric=="S" else "E(f) (approx)")
            plt.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_{metric}.png", dpi=160)
            plt.close()

            # Δ plot
            plt.figure(figsize=(8,5))
            plt.axhline(0, ls="--", color="gray")
            plt.fill_between(f_grid, delta_lo, delta_hi, alpha=0.2, label="Δ 95% CI")
            plt.plot(f_grid, delta_mean, color="C3", lw=2, marker="o", label="Δ mean (AS - null)")
            plt.title(f"Paired Δ {metric}(f) (null={args.null_mode})")
            plt.xlabel("f")
            plt.ylabel("Δ")
            plt.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_{metric}_delta.png", dpi=160)
            plt.close()

            # Inter-community with CI + p
            plt.figure(figsize=(8,5))
            plt.fill_between(f_grid, inter_null_lo, inter_null_hi, alpha=0.2, label="Null inter-comm 95% CI")
            plt.plot(f_grid, inter_null_mean, color="C0", lw=2, label="Null mean")
            plt.plot(f_grid, inter_as, color="C2", lw=2, marker="o", label="AS inter-comm")
            # annotate p-values
            for x, p in zip(f_grid, pvals):
                plt.text(x, max(inter_null_hi.max(), max(inter_as)) * 0.02 + inter_null_hi.min(), 
                         f"p={p:.2g}", rotation=90, va="bottom", ha="center", fontsize=8)
            plt.xlabel("f")
            plt.ylabel("# inter-community edges removed")
            plt.title(f"Inter-community removals (null={args.null_mode})")
            plt.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(f"{args.out_prefix}_intercomm.png", dpi=160)
            plt.close()

        except Exception as e:
            print(f"[plot] skipped: {e}", file=sys.stderr)

        print(f"[done] {metric}: curves={curves_path}, auc={auc_path}, kappa={kappa_path}, inter={inter_path}")

if __name__ == "__main__":
    main()
