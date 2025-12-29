#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
node_feature_models_v1.py

End-to-end analysis of ΔS/ΔE vs. molecular features on edge- OR node-level tables.

"""

import argparse, os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from statsmodels.api import OLS, add_constant
from statsmodels.stats.multitest import multipletests

from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.multioutput import MultiOutputRegressor
from sklearn.metrics import r2_score

# Optional SHAP
try:
    import shap
    HAS_SHAP = True
except Exception:
    HAS_SHAP = False

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ------------------------- helpers -------------------------

def read_table(path):
    # robust flexible separator
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path, sep=r"[,\t]", engine="python")

def colname(df, name):
    for c in df.columns:
        if c.lower() == name.lower():
            return c
    return None

def is_binary_series(s):
    """
    Strict, safe binary detector:
    - Numeric: unique non-nan subset of {0,1}
    - Strings: map yes/no/true/false/y/n/t/f/0/1 for >=95% of non-nan
    """
    sv = s.dropna()
    if sv.empty:
        return False
    num = pd.to_numeric(sv, errors="coerce")
    if num.notna().all():
        vals = set(np.unique(num.values))
        return vals.issubset({0,1,0.0,1.0})
    maps = {"yes":1,"no":0,"true":1,"false":0,"y":1,"n":0,"t":1,"f":0,"1":1,"0":0}
    mapped = sv.astype(str).str.lower().map(maps)
    if mapped.notna().mean() >= 0.95:
        vals = set(np.unique(mapped.dropna().values))
        return vals.issubset({0,1})
    return False

def bh_fdr(pvals):
    pvals = np.asarray(pvals, float)
    if pvals.size == 0:
        return pvals
    _, q, _, _ = multipletests(pvals, method="fdr_bh")
    return q

def _dedup_df_columns(df_in):
    """Ensure unique column names by dropping duplicates (keeps first)."""
    return df_in.loc[:, ~df_in.columns.duplicated()].copy()

# --------------------- univariate tests ---------------------

def univariate_suite(df, y_col, feat_cols):
    rows = []
    y = pd.to_numeric(df[y_col], errors="coerce")
    for c in feat_cols:
        s = df[c]
        if s.isna().all():
            continue
        try:
            if is_binary_series(s):
                maps = {"yes":1,"no":0,"true":1,"false":0,"y":1,"n":0,"t":1,"f":0,"1":1,"0":0}
                s_num = pd.to_numeric(s, errors="coerce")
                if s_num.notna().all():
                    grp = s_num.clip(0,1).round().astype(int)
                    yy = y
                else:
                    grp_map = s.astype(str).str.lower().map(maps)
                    keep = grp_map.notna()
                    grp = grp_map[keep].astype(int)
                    yy = y[keep]
                if grp.nunique() < 2:
                    continue
                x1 = yy[grp==1]; x0 = yy[grp==0]
                if len(x1) < 10 or len(x0) < 10:
                    continue
                U, p = stats.mannwhitneyu(x1, x0, alternative="two-sided")
                eff = float(np.median(x1) - np.median(x0))
                rows.append({"y": y_col, "feature": c, "type": "binary",
                             "stat": U, "p": p, "effect_median_diff": eff,
                             "n1": int(len(x1)), "n0": int(len(x0))})
            else:
                s_num = pd.to_numeric(s, errors="coerce")
                mm = (~s_num.isna())
                if mm.sum() < 20:
                    continue
                rho, p = stats.spearmanr(y[mm], s_num[mm])
                rows.append({"y": y_col, "feature": c, "type": "numeric",
                             "stat": rho, "p": p, "effect_rho": float(rho),
                             "n": int(mm.sum())})
        except Exception as e:
            print(f"[warn] univariate {c}: {e}", file=sys.stderr)
            continue

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    out = []
    for yname, sub in res.groupby("y", sort=False):
        sub = sub.copy()
        sub["q_fdr"] = bh_fdr(sub["p"].values)
        out.append(sub)
    return pd.concat(out, ignore_index=True).sort_values(["y","q_fdr","p"])

# ------------------------- OLS --------------------------

def ols_model(df, y_col, X_cols):
    sub = df[[y_col] + X_cols].copy()
    X = sub[X_cols].apply(pd.to_numeric, errors="coerce")
    keep = ~X.isna().any(axis=1) & ~pd.to_numeric(sub[y_col], errors="coerce").isna()
    X = X[keep]; y = pd.to_numeric(sub[y_col], errors="coerce")[keep]
    if len(X) < 50:
        raise RuntimeError(f"Too few rows for OLS of {y_col}")

    Xz = (X - X.mean(0)) / X.std(0).replace(0, 1)
    Xz = add_constant(Xz, has_constant="add")
    model = OLS(y, Xz).fit()
    est = model.params.drop("const", errors="ignore")
    pvals = model.pvalues.drop("const", errors="ignore")
    out = pd.DataFrame({"predictor": est.index, "coef_std": est.values, "p": pvals.values})
    out["q_fdr"] = bh_fdr(out["p"].values)
    out = out.sort_values("p")
    adj_r2 = model.rsquared_adj
    return out, model, adj_r2

# ---------------- Random Forest + permutation + SHAP ----------------

def fit_rf(df, y_col, num_cols, cat_cols, random_state=1, n_jobs=8, n_estimators=500, max_depth=None, log_path=None):
    """
    Robust RF fitter:
      - Intersects requested columns with df.columns
      - Drops rows with NA in any used column (y + features)
      - Logs shapes and lists to log_path (if provided)
    """
    num_cols = [c for c in num_cols if c in df.columns]
    cat_cols = [c for c in cat_cols if c in df.columns]

    use_cols = num_cols + cat_cols
    needed = [y_col] + use_cols
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing columns for RF {y_col}: {missing}")

    sub = df[needed].copy().dropna(axis=0, how="any")
    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"[fit_rf:{y_col}] rows(after dropna)={len(sub)}, num={len(num_cols)}, cat={len(cat_cols)}\n")
            fh.write(f"[fit_rf:{y_col}] num_cols={num_cols}\n")
            fh.write(f"[fit_rf:{y_col}] cat_cols={cat_cols}\n")

    if len(sub) < 200:
        raise RuntimeError(f"Too few rows for RF of {y_col} (after dropna: {len(sub)})")

    X = sub[use_cols]
    y = pd.to_numeric(sub[y_col], errors="coerce")

    numeric_transformer = StandardScaler(with_mean=True, with_std=True)
    categorical_transformer = OneHotEncoder(handle_unknown="ignore", sparse_output=False)

    pre = ColumnTransformer(
        transformers=[
            ("num", numeric_transformer, num_cols),
            ("cat", categorical_transformer, cat_cols)
        ],
        remainder="drop"
    )

    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=n_jobs,
        max_depth=max_depth
    )

    pipe = Pipeline(steps=[("pre", pre), ("rf", rf)])
    pipe.fit(X, y)

    y_pred = pipe.predict(X)
    r2 = r2_score(y, y_pred)

        # Permutation importance (uses same sub)
    pi = permutation_importance(pipe, X, y, n_repeats=20, random_state=random_state, n_jobs=n_jobs)
    importances = pi.importances_mean

    # === Robust feature names from the preprocessor ===
    try:
        pre = pipe.named_steps["pre"]
        try:
            # sklearn >= 1.0
            feat_names = pre.get_feature_names_out()
        except Exception:
            # Fallback: build a simple index if names are not available
            feat_names = np.array([f"feat_{i}" for i in range(importances.shape[0])])
    except Exception:
        feat_names = np.array([f"feat_{i}" for i in range(importances.shape[0])])

    # Final safety: align lengths if anything is off
    if len(feat_names) != importances.shape[0]:
        # Trim/pad names to match importances length
        if len(feat_names) > importances.shape[0]:
            feat_names = feat_names[:importances.shape[0]]
        else:
            feat_names = np.r_[feat_names, [f"feat_{i}" for i in range(len(feat_names), importances.shape[0])]]

    imp_df = pd.DataFrame({
        "feature": feat_names,
        "perm_importance": importances
    })
    imp_df = imp_df.sort_values("perm_importance", ascending=False)
    
    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"[fit_rf:{y_col}] transformed_features={len(feat_names)}; importances={importances.shape[0]}\n")
            fh.write(f"[fit_rf:{y_col}] example names: {list(feat_names[:10])}\n")

    return pipe, r2, imp_df

def shap_summary(pipe, df, num_cols, cat_cols, y_col, out_png):
    if not HAS_SHAP:
        return False, "shap not installed"
    use_cols = [c for c in (num_cols + cat_cols) if c in df.columns]
    sub = df[[y_col] + use_cols].copy().dropna()
    if len(sub) < 200:
        return False, f"not enough rows for SHAP ({len(sub)})"
    X = sub[use_cols]
    rf = pipe.named_steps["rf"]
    pre = pipe.named_steps["pre"]
    Xt = pre.transform(X)
    try:
        explainer = shap.TreeExplainer(rf)
        shap_values = explainer.shap_values(Xt)
        plt.figure(figsize=(8,6))
        shap.summary_plot(shap_values, Xt, show=False, plot_size=(8,6))
        plt.tight_layout()
        plt.savefig(out_png, dpi=160)
        plt.close()
        return True, None
    except Exception as e:
        return False, str(e)

# -------------------- main orchestration --------------------

def main():
    ap = argparse.ArgumentParser(description="ΔS/ΔE feature modeling: univariate, OLS, RF (+SHAP), joint RF.")
    ap.add_argument("--table", required=True, help="impact.tsv or node_impact_with_features.tsv")
    ap.add_argument("--out", default="feature_models")
    ap.add_argument("--random-state", type=int, default=1)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--n-estimators", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=None)
    ap.add_argument("--binary-as-categorical", action="store_true",
                    help="Treat binary flags as categorical (default: numeric 0/1).")
    ap.add_argument("--col-deltas", default="deltaS", help="Column for ΔS (default: deltaS)")
    ap.add_argument("--col-deltae", default="deltaE", help="Column for ΔE (default: deltaE)")
    ap.add_argument("--feature-whitelist", default="",
                    help="Comma-separated list of feature columns to use; if empty, auto-detect.")
    args = ap.parse_args()

    df = read_table(args.table)

    # If no 'gene' but edge endpoints exist, create surrogate ID
    if "gene" not in df.columns and {"A","B"}.issubset(df.columns):
        df["gene"] = df["A"].astype(str) + "_" + df["B"].astype(str)

    # Targets
    dS_col = colname(df, args.col_deltas) or args.col_deltas
    dE_col = colname(df, args.col_deltae) or args.col_deltae
    for need in [dS_col, dE_col]:
        if need not in df.columns:
            sys.exit(f"Missing required column in --table: {need}. Have: {list(df.columns)}")

    # Build feature list
    if args.feature_whitelist.strip():
        feature_cols = [f.strip() for f in args.feature_whitelist.split(",") if f.strip() in df.columns]
    else:
        # exclude targets and pure identifiers
        exclude = {dS_col, dE_col, "edge_id", "A", "B", "gene", "AS_score"}
        feature_cols = [c for c in df.columns if c not in exclude]

    # Classify features
    numeric_feats, binary_feats, categorical_feats = [], [], []
    for c in feature_cols:
        s = df[c]
        if is_binary_series(s):
            if args.binary_as_categorical:
                categorical_feats.append(c)
            else:
                binary_feats.append(c)
            continue
        s_num = pd.to_numeric(s, errors="coerce")
        if s_num.notna().mean() > 0.7:
            numeric_feats.append(c)
        else:
            categorical_feats.append(c)

    os.makedirs(args.out, exist_ok=True)

    # Preflight log
    with open(f"{args.out}_preflight.txt","w") as fh:
        fh.write(f"Targets: dS={dS_col}, dE={dE_col}\n")
        fh.write(f"numeric_feats ({len(numeric_feats)}): {numeric_feats}\n")
        fh.write(f"binary_feats ({len(binary_feats)}): {binary_feats}\n")
        fh.write(f"categorical_feats ({len(categorical_feats)}): {categorical_feats}\n")

    # Correlation heatmap among numeric features (and targets)
    corr_cols = [dS_col, dE_col] + numeric_feats
    C = df[corr_cols].apply(pd.to_numeric, errors="coerce").dropna()
    if len(C) >= 50:
        corr = C.corr(method="spearman")
        plt.figure(figsize=(8,7))
        im = plt.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.xticks(range(len(corr.columns)), corr.columns, rotation=60, ha="right", fontsize=8)
        plt.yticks(range(len(corr.index)), corr.index, fontsize=8)
        plt.title("Spearman correlation heatmap (numeric)")
        plt.tight_layout()
        plt.savefig(f"{args.out}_corr_heatmap.png", dpi=180)
        plt.close()

    # Scatter ΔE vs ΔS
    plt.figure(figsize=(5,5))
    plt.scatter(pd.to_numeric(df[dS_col], errors="coerce"),
                pd.to_numeric(df[dE_col], errors="coerce"), s=8, alpha=0.4)
    plt.xlabel("ΔS"); plt.ylabel("ΔE"); plt.title("ΔE vs ΔS")
    plt.tight_layout(); plt.savefig(f"{args.out}_scatter_dE_vs_dS.png", dpi=160); plt.close()

    # ---------------- Univariate ----------------
    uni = pd.concat([
        univariate_suite(df, dS_col, feature_cols),
        univariate_suite(df, dE_col, feature_cols)
    ], ignore_index=True)
    uni.to_csv(f"{args.out}_univariate.tsv", sep="\t", index=False)

    # ---------------- Prepare matrices for multivariate ----------------
    X_base = df[numeric_feats].apply(pd.to_numeric, errors="coerce")
    for b in binary_feats:
        X_base[b] = pd.to_numeric(df[b], errors="coerce").round().clip(0,1)
    YdS = pd.to_numeric(df[dS_col], errors="coerce")
    YdE = pd.to_numeric(df[dE_col], errors="coerce")

    keep_rows = (~YdS.isna()) & (~YdE.isna())
    for c in X_base.columns:
        keep_rows &= ~X_base[c].isna()
    X_base = X_base[keep_rows].copy()
    YdS = YdS[keep_rows]; YdE = YdE[keep_rows]
    cats = [c for c in categorical_feats if c in df.columns]
    X_cats = df.loc[keep_rows, cats].copy()

    # ---------------- OLS (numeric/binary only) ----------------
    X_ols = X_base.copy()
    try:
        olsS_df, olsS_model, adjR2_S = ols_model(pd.concat([YdS.rename("dS"), X_ols], axis=1), "dS", list(X_ols.columns))
        olsS_df.to_csv(f"{args.out}_ols_deltaS.tsv", sep="\t", index=False)
        with open(f"{args.out}_ols_deltaS_model.txt","w") as fh:
            fh.write(olsS_model.summary().as_text() + f"\n\nAdj R2: {adjR2_S:.4f}\n")
    except Exception as e:
        with open(f"{args.out}_ols_deltaS_model.txt","w") as fh:
            fh.write(f"OLS ΔS failed: {e}\n")

    try:
        olsE_df, olsE_model, adjR2_E = ols_model(pd.concat([YdE.rename("dE"), X_ols], axis=1), "dE", list(X_ols.columns))
        olsE_df.to_csv(f"{args.out}_ols_deltaE.tsv", sep="\t", index=False)
        with open(f"{args.out}_ols_deltaE_model.txt","w") as fh:
            fh.write(olsE_model.summary().as_text() + f"\n\nAdj R2: {adjR2_E:.4f}\n")
    except Exception as e:
        with open(f"{args.out}_ols_deltaE_model.txt","w") as fh:
            fh.write(f"OLS ΔE failed: {e}\n")

    # ---------------- RF (separate models) ----------------
    rf_num = list(X_base.columns) + binary_feats
    rf_cat = cats
    rf_num = list(dict.fromkeys(rf_num))
    rf_cat = list(dict.fromkeys(rf_cat))
    overlap = set(rf_num) & set(rf_cat)
    if overlap:
        rf_cat = [c for c in rf_cat if c not in overlap]

    with open(f"{args.out}_preflight.txt","a") as fh:
        fh.write(f"\nRF num ({len(rf_num)}): {rf_num}\n")
        fh.write(f"RF cat ({len(rf_cat)}): {rf_cat}\n")

    # ΔS RF
    try:
        pipeS, rfS_r2, rfS_imp = fit_rf(
            df, dS_col, rf_num, rf_cat,
            random_state=args.random_state, n_jobs=args.n_jobs,
            n_estimators=args.n_estimators, max_depth=args.max_depth,
            log_path=f"{args.out}_rf_debug.log"
        )
        rfS_imp.to_csv(f"{args.out}_rf_deltaS_importance.tsv", sep="\t", index=False)
        plt.figure(figsize=(7,6))
        top = rfS_imp.head(20)
        plt.barh(range(len(top)), top["perm_importance"], alpha=0.8)
        plt.yticks(range(len(top)), top["feature"])
        plt.gca().invert_yaxis()
        plt.xlabel("Permutation importance"); plt.title(f"RF ΔS (R²={rfS_r2:.3f})")
        plt.tight_layout(); plt.savefig(f"{args.out}_rf_deltaS_importance.png", dpi=180); plt.close()

        ok, msg = shap_summary(pipeS, df, rf_num, rf_cat, dS_col, f"{args.out}_shap_deltaS.png")
        if not ok and msg:
            with open(f"{args.out}_shap_log.txt","a") as fh: fh.write(f"ΔS SHAP: {msg}\n")
    except Exception as e:
        with open(f"{args.out}_rf_deltaS_log.txt","w") as fh:
            fh.write(f"RF ΔS failed: {e}\n")

    # ΔE RF
    try:
        pipeE, rfE_r2, rfE_imp = fit_rf(
            df, dE_col, rf_num, rf_cat,
            random_state=args.random_state, n_jobs=args.n_jobs,
            n_estimators=args.n_estimators, max_depth=args.max_depth,
            log_path=f"{args.out}_rf_debug.log"
        )
        rfE_imp.to_csv(f"{args.out}_rf_deltaE_importance.tsv", sep="\t", index=False)
        plt.figure(figsize=(7,6))
        top = rfE_imp.head(20)
        plt.barh(range(len(top)), top["perm_importance"], alpha=0.8)
        plt.yticks(range(len(top)), top["feature"])
        plt.gca().invert_yaxis()
        plt.xlabel("Permutation importance"); plt.title(f"RF ΔE (R²={rfE_r2:.3f})")
        plt.tight_layout(); plt.savefig(f"{args.out}_rf_deltaE_importance.png", dpi=180); plt.close()

        ok, msg = shap_summary(pipeE, df, rf_num, rf_cat, dE_col, f"{args.out}_shap_deltaE.png")
        if not ok and msg:
            with open(f"{args.out}_shap_log.txt","a") as fh: fh.write(f"ΔE SHAP: {msg}\n")
    except Exception as e:
        with open(f"{args.out}_rf_deltaE_log.txt","w") as fh:
            fh.write(f"RF ΔE failed: {e}\n")

    # ---------------- Joint multi-output RF ----------------
    try:
        rf_num = list(dict.fromkeys(rf_num))
        rf_cat = list(dict.fromkeys(rf_cat))
        overlap = set(rf_num) & set(rf_cat)
        if overlap:
            rf_cat = [c for c in rf_cat if c not in overlap]

        rf_num = [c for c in rf_num if c in df.columns]
        rf_cat = [c for c in rf_cat if c in df.columns]

        cols_joint = list(dict.fromkeys(rf_num + rf_cat))
        needed = [dS_col, dE_col] + cols_joint
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise RuntimeError(f"Missing columns for joint RF: {missing}")

        sub = df[needed].copy().dropna(axis=0, how="any")
        if len(sub) < 200:
            raise RuntimeError(f"Too few rows for joint RF: {len(sub)}")

        sub = _dedup_df_columns(sub)
        X_joint = sub[cols_joint]
        Y_joint = sub[[dS_col, dE_col]].values

        numeric_transformer = StandardScaler(with_mean=True, with_std=True)
        categorical_transformer = OneHotEncoder(handle_unknown="ignore", sparse_output=False)

        pre = ColumnTransformer(
            transformers=[
                ("num", numeric_transformer, [c for c in cols_joint if c in rf_num]),
                ("cat", categorical_transformer, [c for c in cols_joint if c in rf_cat]),
            ],
            remainder="drop"
        )

        base_rf = RandomForestRegressor(
            n_estimators=args.n_estimators,
            random_state=args.random_state,
            n_jobs=args.n_jobs,
            max_depth=args.max_depth
        )
        morf = MultiOutputRegressor(base_rf, n_jobs=args.n_jobs)
        pipe_joint = Pipeline(steps=[("pre", pre), ("morf", morf)])
        pipe_joint.fit(X_joint, Y_joint)
        Y_pred = pipe_joint.predict(X_joint)

        r2S = r2_score(Y_joint[:,0], Y_pred[:,0])
        r2E = r2_score(Y_joint[:,1], Y_pred[:,1])
        pd.DataFrame({"target":["deltaS","deltaE"], "R2":[r2S, r2E]}).to_csv(
            f"{args.out}_rf_joint_r2.tsv", sep="\t", index=False
        )
    except Exception as e:
        with open(f"{args.out}_rf_joint_log.txt","w") as fh:
            fh.write(f"Joint RF failed: {e}\n")

    print("[done] Outputs in:", args.out)

if __name__ == "__main__":
    main()
