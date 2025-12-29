#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
node_feature_models_catboost_with_interactions_v3.py

CatBoost + SHAP (values + interactions) for AS impact modeling, with safe OLS.
- Categorical handling fixed for CatBoost (strings, 'NA' for missing).
- OLS dummy-codes categoricals and uses only numeric matrix.
- SHAP interaction values computed for ALL features (heatmap + curated pairs).
- Optional quantile-binning for large-scale count-like features.

Inputs:
  --table <TSV>           Impact/feature table (e.g., impact.tsv)
  --target {combined,deltaS,deltaE}  Target to model (default: combined)
  --feature-cols <csv>    Explicit feature list (optional; else default panel)
  --cat-cols <csv>        Extra categorical columns (if present)
  --bin-quantiles "NAME[:k],NAME2[:k2]"  Quantile-bin selected numeric features
  --outdir <dir>          Output directory
  --prefix <str>          Output filename prefix (no ext)
  --model-out <path>      Path for saved CatBoost model (.cbm)
  --test-size <float>     Holdout fraction (default 0.2)
  --random-state <int>    Seed (default 13)
  --iterations <int>      CatBoost trees (default 1500)
  --depth <int>           CatBoost depth (default 6)
  --learning-rate <float> CatBoost LR (default 0.03)
  --l2 <float>            CatBoost L2 regularization (default 3.0)
  --shap-n <int>          Rows to use for SHAP (0 = all)
  --shap-approx           Approximate SHAP for speed

Outputs (in --outdir, prefix = --prefix):
  *_preflight.txt
  *_holdout_r2.txt
  *_perm_importance.tsv
  *_ols_model.txt
  *_ols.tsv
  *_shap_summary.tsv
  *_shap_summary.png
  *_shap_interactions_mean.tsv
  *_shap_interactions_heatmap.png
  *_shap_dependence_<A>__x__<B>.png
  *_model.cbm

Requires: catboost, shap, pandas, numpy, matplotlib, seaborn, scikit-learn, statsmodels
"""

import argparse, os, sys, json, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score

import statsmodels.api as sm

from catboost import CatBoostRegressor, Pool
import shap


# --------------------------- helpers ---------------------------

def parse_kv_list(s):
    """Parse NAME or NAME:k entries into dict of {NAME: k or None}."""
    out = {}
    if not s:
        return out
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            kname, kval = tok.split(":", 1)
            out[kname.strip()] = int(kval.strip())
        else:
            out[tok] = None
    return out

def bin_quantiles(df, name, k=None):
    """Create quantile bins as categorical labels for feature 'name'."""
    x = pd.to_numeric(df[name], errors="coerce")
    k = 3 if (k is None) else k
    # Rank-based qcut is more stable with ties
    b = pd.qcut(x.rank(method="first"), q=min(k, x.nunique()), labels=False, duplicates="drop")
    return b.astype("Int64").astype(str).fillna("NA")  # explicit string tokens

def safe_exists(cols, name):
    return name in cols


# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--table", required=True, help="TSV with features & targets")
    ap.add_argument("--target", default="combined", choices=["combined","deltaS","deltaE"],
                    help="Which target to model")
    ap.add_argument("--feature-cols", default="",
                    help="Comma-separated feature names to use. If empty, use default panel.")
    ap.add_argument("--cat-cols", default="",
                    help="Comma-separated categorical features (if present).")
    ap.add_argument("--bin-quantiles", default="",
                    help="Optional quantile binning: NAME[:k],NAME2[:k2] ...")
    ap.add_argument("--model-out", default="", help="Where to save CatBoost model (.cbm). If empty, uses <outdir>/<prefix>_model.cbm")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--prefix", default="feature_models_cat_all",
                    help="Output filename prefix (no extension)")
    ap.add_argument("--test-size", type=float, default=0.2)
    ap.add_argument("--random-state", type=int, default=13)
    ap.add_argument("--iterations", type=int, default=1500)
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--learning-rate", type=float, default=0.03)
    ap.add_argument("--l2", type=float, default=3.0)
    ap.add_argument("--shap-n", type=int, default=0, help="0 = use all rows; else subsample for SHAP")
    ap.add_argument("--shap-approx", action="store_true", help="Use approximate SHAP for speed")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.prefix)

    # ---------- read ----------
    df = pd.read_csv(args.table, sep="\t", low_memory=False)

    # Targets
    if args.target == "combined":
        if "impact_combined" in df.columns:
            y = pd.to_numeric(df["impact_combined"], errors="coerce")
        else:
            if not {"deltaS","deltaE"} <= set(df.columns):
                sys.exit("Need deltaS and deltaE or impact_combined in the table.")
            y = pd.to_numeric(df["deltaS"], errors="coerce") + pd.to_numeric(df["deltaE"], errors="coerce")
    elif args.target == "deltaS":
        y = pd.to_numeric(df["deltaS"], errors="coerce")
    else:
        y = pd.to_numeric(df["deltaE"], errors="coerce")

    # Features: default panel if not specified
    default_feats = ['is_as_controlled', 'exon_len', 'as_score_agg', 'tissue', 'length', 'REGULATION', 'PFAM_MODULES', 'DEGREE', 'SUBCELLULAR', 'DOMAIN_VARIETY', 'PATHWAYS', 'DISORDER', 'PTMS', 'COMPLEX', 'gini', 'EVENTS', 'degree_A', 'degree_B', 'degree_diff', 'clustering_A', 'clustering_B', 'shared_neighbors', 'jaccard_coeff', 'community_A', 'community_B', 'community_diff', 'intercomm', 'edge_betweenness', 'log_length', 'length_sq', 'micro', 'gini_sq', 'len_micro_interact', 'A_emb_1', 'A_emb_2', 'A_emb_3', 'A_emb_4', 'A_emb_5', 'A_emb_6', 'A_emb_7', 'A_emb_8', 'A_emb_9', 'A_emb_10', 'A_emb_11', 'A_emb_12', 'A_emb_13', 'A_emb_14', 'A_emb_15', 'A_emb_16', 'A_emb_17', 'A_emb_18', 'A_emb_19', 'A_emb_20', 'A_emb_21', 'A_emb_22', 'A_emb_23', 'A_emb_24', 'A_emb_25', 'A_emb_26', 'A_emb_27', 'A_emb_28', 'A_emb_29', 'A_emb_30', 'A_emb_31', 'A_emb_32', 'A_emb_33', 'A_emb_34', 'A_emb_35', 'A_emb_36', 'A_emb_37', 'A_emb_38', 'A_emb_39', 'A_emb_40', 'A_emb_41', 'A_emb_42', 'A_emb_43', 'A_emb_44', 'A_emb_45', 'A_emb_46', 'A_emb_47', 'A_emb_48', 'A_emb_49', 'A_emb_50', 'A_emb_51', 'A_emb_52', 'A_emb_53', 'A_emb_54', 'A_emb_55', 'A_emb_56', 'A_emb_57', 'A_emb_58', 'A_emb_59', 'A_emb_60', 'A_emb_61', 'A_emb_62', 'A_emb_63', 'A_emb_64', 'A_emb_65', 'A_emb_66', 'A_emb_67', 'A_emb_68', 'A_emb_69', 'A_emb_70', 'A_emb_71', 'A_emb_72', 'A_emb_73', 'A_emb_74', 'A_emb_75', 'A_emb_76', 'A_emb_77', 'A_emb_78', 'A_emb_79', 'A_emb_80', 'A_emb_81', 'A_emb_82', 'A_emb_83', 'A_emb_84', 'A_emb_85', 'A_emb_86', 'A_emb_87', 'A_emb_88', 'A_emb_89', 'A_emb_90', 'A_emb_91', 'A_emb_92', 'A_emb_93', 'A_emb_94', 'A_emb_95', 'A_emb_96', 'A_emb_97', 'A_emb_98', 'A_emb_99', 'A_emb_100', 'A_emb_101', 'A_emb_102', 'A_emb_103', 'A_emb_104', 'A_emb_105', 'A_emb_106', 'A_emb_107', 'A_emb_108', 'A_emb_109', 'A_emb_110', 'A_emb_111', 'A_emb_112', 'A_emb_113', 'A_emb_114', 'A_emb_115', 'A_emb_116', 'A_emb_117', 'A_emb_118', 'A_emb_119', 'A_emb_120', 'A_emb_121', 'A_emb_122', 'A_emb_123', 'A_emb_124', 'A_emb_125', 'A_emb_126', 'A_emb_127', 'A_emb_128', 'B_emb_1', 'B_emb_2', 'B_emb_3', 'B_emb_4', 'B_emb_5', 'B_emb_6', 'B_emb_7', 'B_emb_8', 'B_emb_9', 'B_emb_10', 'B_emb_11', 'B_emb_12', 'B_emb_13', 'B_emb_14', 'B_emb_15', 'B_emb_16', 'B_emb_17', 'B_emb_18', 'B_emb_19', 'B_emb_20', 'B_emb_21', 'B_emb_22', 'B_emb_23', 'B_emb_24', 'B_emb_25', 'B_emb_26', 'B_emb_27', 'B_emb_28', 'B_emb_29', 'B_emb_30', 'B_emb_31', 'B_emb_32', 'B_emb_33', 'B_emb_34', 'B_emb_35', 'B_emb_36', 'B_emb_37', 'B_emb_38', 'B_emb_39', 'B_emb_40', 'B_emb_41', 'B_emb_42', 'B_emb_43', 'B_emb_44', 'B_emb_45', 'B_emb_46', 'B_emb_47', 'B_emb_48', 'B_emb_49', 'B_emb_50', 'B_emb_51', 'B_emb_52', 'B_emb_53', 'B_emb_54', 'B_emb_55', 'B_emb_56', 'B_emb_57', 'B_emb_58', 'B_emb_59', 'B_emb_60', 'B_emb_61', 'B_emb_62', 'B_emb_63', 'B_emb_64', 'B_emb_65', 'B_emb_66', 'B_emb_67', 'B_emb_68', 'B_emb_69', 'B_emb_70', 'B_emb_71', 'B_emb_72', 'B_emb_73', 'B_emb_74', 'B_emb_75', 'B_emb_76', 'B_emb_77', 'B_emb_78', 'B_emb_79', 'B_emb_80', 'B_emb_81', 'B_emb_82', 'B_emb_83', 'B_emb_84', 'B_emb_85', 'B_emb_86', 'B_emb_87', 'B_emb_88', 'B_emb_89', 'B_emb_90', 'B_emb_91', 'B_emb_92', 'B_emb_93', 'B_emb_94', 'B_emb_95', 'B_emb_96', 'B_emb_97', 'B_emb_98', 'B_emb_99', 'B_emb_100', 'B_emb_101', 'B_emb_102', 'B_emb_103', 'B_emb_104', 'B_emb_105', 'B_emb_106', 'B_emb_107', 'B_emb_108', 'B_emb_109', 'B_emb_110', 'B_emb_111', 'B_emb_112', 'B_emb_113', 'B_emb_114', 'B_emb_115', 'B_emb_116', 'B_emb_117', 'B_emb_118', 'B_emb_119', 'B_emb_120', 'B_emb_121', 'B_emb_122', 'B_emb_123', 'B_emb_124', 'B_emb_125', 'B_emb_126', 'B_emb_127', 'B_emb_128', 'emb_dot', 'emb_cos', 'emb_euclid', 'emb_normed', 'emb_sum_1', 'emb_diff_1', 'emb_sum_2', 'emb_diff_2', 'emb_sum_3', 'emb_diff_3', 'emb_sum_4', 'emb_diff_4', 'emb_sum_5', 'emb_diff_5', 'emb_sum_6', 'emb_diff_6', 'emb_sum_7', 'emb_diff_7', 'emb_sum_8', 'emb_diff_8', 'emb_sum_9', 'emb_diff_9', 'emb_sum_10', 'emb_diff_10', 'emb_sum_11', 'emb_diff_11', 'emb_sum_12', 'emb_diff_12', 'emb_sum_13', 'emb_diff_13', 'emb_sum_14', 'emb_diff_14', 'emb_sum_15', 'emb_diff_15', 'emb_sum_16', 'emb_diff_16', 'emb_sum_17', 'emb_diff_17', 'emb_sum_18', 'emb_diff_18', 'emb_sum_19', 'emb_diff_19', 'emb_sum_20', 'emb_diff_20', 'emb_sum_21', 'emb_diff_21', 'emb_sum_22', 'emb_diff_22', 'emb_sum_23', 'emb_diff_23', 'emb_sum_24', 'emb_diff_24', 'emb_sum_25', 'emb_diff_25', 'emb_sum_26', 'emb_diff_26', 'emb_sum_27', 'emb_diff_27', 'emb_sum_28', 'emb_diff_28', 'emb_sum_29', 'emb_diff_29', 'emb_sum_30', 'emb_diff_30', 'emb_sum_31', 'emb_diff_31', 'emb_sum_32', 'emb_diff_32', 'emb_sum_33', 'emb_diff_33', 'emb_sum_34', 'emb_diff_34', 'emb_sum_35', 'emb_diff_35', 'emb_sum_36', 'emb_diff_36', 'emb_sum_37', 'emb_diff_37', 'emb_sum_38', 'emb_diff_38', 'emb_sum_39', 'emb_diff_39', 'emb_sum_40', 'emb_diff_40', 'emb_sum_41', 'emb_diff_41', 'emb_sum_42', 'emb_diff_42', 'emb_sum_43', 'emb_diff_43', 'emb_sum_44', 'emb_diff_44', 'emb_sum_45', 'emb_diff_45', 'emb_sum_46', 'emb_diff_46', 'emb_sum_47', 'emb_diff_47', 'emb_sum_48', 'emb_diff_48', 'emb_sum_49', 'emb_diff_49', 'emb_sum_50', 'emb_diff_50', 'emb_sum_51', 'emb_diff_51', 'emb_sum_52', 'emb_diff_52', 'emb_sum_53', 'emb_diff_53', 'emb_sum_54', 'emb_diff_54', 'emb_sum_55', 'emb_diff_55', 'emb_sum_56', 'emb_diff_56', 'emb_sum_57', 'emb_diff_57', 'emb_sum_58', 'emb_diff_58', 'emb_sum_59', 'emb_diff_59', 'emb_sum_60', 'emb_diff_60', 'emb_sum_61', 'emb_diff_61', 'emb_sum_62', 'emb_diff_62', 'emb_sum_63', 'emb_diff_63', 'emb_sum_64', 'emb_diff_64', 'emb_sum_65', 'emb_diff_65', 'emb_sum_66', 'emb_diff_66', 'emb_sum_67', 'emb_diff_67', 'emb_sum_68', 'emb_diff_68', 'emb_sum_69', 'emb_diff_69', 'emb_sum_70', 'emb_diff_70', 'emb_sum_71', 'emb_diff_71', 'emb_sum_72', 'emb_diff_72', 'emb_sum_73', 'emb_diff_73', 'emb_sum_74', 'emb_diff_74', 'emb_sum_75', 'emb_diff_75', 'emb_sum_76', 'emb_diff_76', 'emb_sum_77', 'emb_diff_77', 'emb_sum_78', 'emb_diff_78', 'emb_sum_79', 'emb_diff_79', 'emb_sum_80', 'emb_diff_80', 'emb_sum_81', 'emb_diff_81', 'emb_sum_82', 'emb_diff_82', 'emb_sum_83', 'emb_diff_83', 'emb_sum_84', 'emb_diff_84', 'emb_sum_85', 'emb_diff_85', 'emb_sum_86', 'emb_diff_86', 'emb_sum_87', 'emb_diff_87', 'emb_sum_88', 'emb_diff_88', 'emb_sum_89', 'emb_diff_89', 'emb_sum_90', 'emb_diff_90', 'emb_sum_91', 'emb_diff_91', 'emb_sum_92', 'emb_diff_92', 'emb_sum_93', 'emb_diff_93', 'emb_sum_94', 'emb_diff_94', 'emb_sum_95', 'emb_diff_95', 'emb_sum_96', 'emb_diff_96', 'emb_sum_97', 'emb_diff_97', 'emb_sum_98', 'emb_diff_98', 'emb_sum_99', 'emb_diff_99', 'emb_sum_100', 'emb_diff_100', 'emb_sum_101', 'emb_diff_101', 'emb_sum_102', 'emb_diff_102', 'emb_sum_103', 'emb_diff_103', 'emb_sum_104', 'emb_diff_104', 'emb_sum_105', 'emb_diff_105', 'emb_sum_106', 'emb_diff_106', 'emb_sum_107', 'emb_diff_107', 'emb_sum_108', 'emb_diff_108', 'emb_sum_109', 'emb_diff_109', 'emb_sum_110', 'emb_diff_110', 'emb_sum_111', 'emb_diff_111', 'emb_sum_112', 'emb_diff_112', 'emb_sum_113', 'emb_diff_113', 'emb_sum_114', 'emb_diff_114', 'emb_sum_115', 'emb_diff_115', 'emb_sum_116', 'emb_diff_116', 'emb_sum_117', 'emb_diff_117', 'emb_sum_118', 'emb_diff_118', 'emb_sum_119', 'emb_diff_119', 'emb_sum_120', 'emb_diff_120', 'emb_sum_121', 'emb_diff_121', 'emb_sum_122', 'emb_diff_122', 'emb_sum_123', 'emb_diff_123', 'emb_sum_124', 'emb_diff_124', 'emb_sum_125', 'emb_diff_125', 'emb_sum_126', 'emb_diff_126', 'emb_sum_127', 'emb_diff_127', 'emb_sum_128', 'emb_diff_128', 'degree_mean', 'degree_min', 'degree_max', 'degree_ratio', 'same_comm', 'local_redundancy', 'log_betweenness', 'zscore_betweenness', 'log_gini', 'log_degree_mean']
    #default_feats = ['is_as_controlled','exon_len','as_score_agg','details','tissue','length','REGULATION','PFAM_MODULES','DEGREE','SUBCELLULAR','DOMAIN_VARIETY','PATHWAYS','DISORDER','PTMS','COMPLEX','GINI','EVENTS','degree_A','degree_B','degree_diff','clustering_A','clustering_B','shared_neighbors','jaccard_coeff','community_A','community_B','community_diff','edge_betweenness','log_length','length_sq','micro','gini_sq','len_micro_interact','A_length','A_degree_x','A_degree_y','A_gini','A_domain_variety','A_subcellular','A_pathways','B_length','B_degree_x','B_degree_y','B_gini','B_domain_variety','B_subcellular','B_pathways','edge_degree_x_mean','edge_degree_x_max','edge_degree_x_absdiff','edge_degree_y_same','edge_degree_y_pair','edge_domain_variety_same','edge_domain_variety_pair','edge_gini_same','edge_gini_pair','edge_length_same','edge_length_pair','edge_pathways_same','edge_pathways_pair','edge_subcellular_same','edge_subcellular_pair']
    #default_feats = ['length','log_length','micro','tissue_flag','degree_x','regulation','gini','pfam_modules','degree_y','subcellular','domain_variety','pathways','disorder','ptms', 'complex']
    #default_feats = [
    #    "REGULATION_min","PTMS_min","DISORDER_min","REGULATION_mean","PTMS_mean","DISORDER_mean",
    #    "PFAM_MODULES_mean","PFAM_MODULES_min","DOMAIN_VARIETY_mean","DOMAIN_VARIETY_min",
    #    "SUBCELLULAR_mean","SUBCELLULAR_min","COMPLEX_mean","COMPLEX_min",
    #    "mean_degree","min_degree","length","log_length","intercomm","tissue_flag","micro","PATHWAYS_mean"
    #]
    feats = [c for c in (args.feature_cols.split(",") if args.feature_cols else default_feats) if c]
    feats = [f for f in feats if f in df.columns]  # keep existing only

    # Optional: quantile bin selected features
    bin_map = parse_kv_list(args.bin_quantiles)
    for name, kval in bin_map.items():
        if name in df.columns:
            df[name + "_qbin"] = bin_quantiles(df, name, kval)
            feats.append(name + "_qbin")

    # Categorical columns
    user_cats = [c for c in (args.cat_cols.split(",") if args.cat_cols else []) if c]
    auto_cats = [c for c in ["tissue_flag","intercomm","micro"] if c in df.columns]
    qbin_cats = [c for c in df.columns if c.endswith("_qbin")]
    cat_cols = list(dict.fromkeys(user_cats + auto_cats + qbin_cats))

    # Build X
    X = df[feats].copy()

    # ---- Cast categoricals to STRING tokens for CatBoost ----
    cat_cols_in_X = [c for c in cat_cols if c in X.columns]
    for c in cat_cols_in_X:
        X[c] = (
            X[c]
            .astype(object)            # leave numeric dtype
            .where(~X[c].isna(), 'NA') # explicit NA
            .astype(str)               # cat tokens
        )
    # CatBoost needs categorical *indices* in current X
    cat_ix = [i for i, col in enumerate(X.columns) if col in cat_cols_in_X]

    # Split
    mask_valid = y.notna() & (~X.isna().any(axis=1))
    X = X.loc[mask_valid].copy()
    y = y.loc[mask_valid].astype(float)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.random_state
    )

    train_pool = Pool(X_train, y_train, cat_features=cat_ix)
    test_pool  = Pool(X_test,  y_test,  cat_features=cat_ix)

    # ---- CatBoost model ----
    model = CatBoostRegressor(
        loss_function="RMSE",
        eval_metric="RMSE",
        depth=args.depth,
        learning_rate=args.learning_rate,
        l2_leaf_reg=args.l2,
        iterations=args.iterations,
        random_seed=args.random_state,
        verbose=False
    )
    model.fit(train_pool, eval_set=test_pool, use_best_model=True)

    # Save model
    model_out = args.model_out if args.model_out else (prefix + "_model.cbm")
    model.save_model(model_out)

    # Holdout R2
    y_pred = model.predict(test_pool)
    r2 = r2_score(y_test, y_pred)
    with open(prefix + "_holdout_r2.txt", "w") as fh:
        fh.write(f"R2_test\t{r2:.6f}\n")

    # Permutation importance
    p_imp = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=args.random_state)
    imp_df = pd.DataFrame({
        "feature": X.columns,
        "perm_importance_mean": p_imp.importances_mean,
        "perm_importance_sd": p_imp.importances_std
    }).sort_values("perm_importance_mean", ascending=False)
    imp_df.to_csv(prefix + "_perm_importance.tsv", sep="\t", index=False)

    # ---- OLS (optional sanity check, numeric-only design) ----
    try:
        # Split features by type
        num_cols = [c for c in X.columns if c not in cat_cols_in_X]
        Xn = X[num_cols].apply(pd.to_numeric, errors="coerce")

        # One-hot for categoricals
        Xc = pd.get_dummies(X[cat_cols_in_X].astype("category"), drop_first=True, dtype=float)

        # Combine, clean, align
        X_ols = pd.concat([Xn, Xc], axis=1).replace([np.inf, -np.inf], np.nan)
        mask = X_ols.notna().all(axis=1) & y.notna()
        X_ols = X_ols.loc[mask]
        y_ols = y.loc[mask].astype(float)

        X_ols = sm.add_constant(X_ols, has_constant="add")
        ols_model = sm.OLS(y_ols, X_ols).fit()
        with open(prefix + "_ols_model.txt", "w") as fh:
            fh.write(ols_model.summary().as_text())
        pd.DataFrame({
            "term": ols_model.params.index,
            "coef": ols_model.params.values,
            "pval": ols_model.pvalues.values
        }).to_csv(prefix + "_ols.tsv", sep="\t", index=False)
    except Exception as e:
        with open(prefix + "_ols_model.txt", "w") as fh:
            fh.write(f"OLS skipped/failed: {e}\n")

    # ---- SHAP values (optionally subsample) ----
    if args.shap_n and args.shap_n > 0 and len(X) > args.shap_n:
        X_shap = X.sample(args.shap_n, random_state=args.random_state)
    else:
        X_shap = X

    explainer = shap.TreeExplainer(
        model,
        feature_perturbation="interventional",
        approximate=args.shap_approx
    )
    shap_vals = explainer.shap_values(X_shap)

    # SHAP summary (TSV + beeswarm)
    mean_abs = np.abs(shap_vals).mean(axis=0)
    shap_sum = pd.DataFrame({"feature": X_shap.columns, "mean_abs_shap": mean_abs})
    shap_sum = shap_sum.sort_values("mean_abs_shap", ascending=False)
    shap_sum.to_csv(prefix + "_shap_summary.tsv", sep="\t", index=False)

    plt.figure(figsize=(8, min(12, 0.35*len(X_shap.columns)+2)))
    shap.summary_plot(shap_vals, X_shap, show=False, max_display=len(X_shap.columns))
    plt.tight_layout()
    plt.savefig(prefix + "_shap_summary.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ---- SHAP interaction values (ALL features) ----
    shap_int = explainer.shap_interaction_values(X_shap)
    M = np.abs(shap_int).mean(axis=0)  # [feat x feat]
    np.fill_diagonal(M, 0.0)
    M_df = pd.DataFrame(M, index=X_shap.columns, columns=X_shap.columns)
    M_df.to_csv(prefix + "_shap_interactions_mean.tsv", sep="\t")

    plt.figure(figsize=(12, 10))
    sns.heatmap(M_df, cmap="mako", square=True, cbar_kws={"label": "mean |SHAP interaction|"})
    plt.title("Mean |SHAP interaction| (all features)")
    plt.tight_layout()
    plt.savefig(prefix + "_shap_interactions_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ---- Curated dependence plots for biologically interesting pairs ----
    pairs = [
        ("REGULATION_mean","DISORDER_min"),
        ("REGULATION_min","PTMS_min"),
        ("PTMS_min","PFAM_MODULES_mean"),
        ("log_length","micro"),
        ("length","micro"),
        ("mean_degree","REGULATION_min"),
        ("intercomm","REGULATION_mean"),
        ("tissue_flag","REGULATION_min"),
        ("PTMS_mean","DISORDER_mean"),
        ("PFAM_MODULES_mean","DOMAIN_VARIETY_mean"),
    ]
    for a, b in pairs:
        if (a in X_shap.columns) and (b in X_shap.columns):
            plt.figure(figsize=(6.4, 5.0))
            shap.dependence_plot(a, shap_vals, X_shap, interaction_index=b, show=False, dot_size=8)
            plt.tight_layout()
            plt.savefig(prefix + f"_shap_dependence_{a}__x__{b}.png",
                        dpi=200, bbox_inches="tight")
            plt.close()

    # ---- Preflight summary ----
    with open(prefix + "_preflight.txt","w") as fh:
        fh.write(textwrap.dedent(f"""
        table: {args.table}
        target: {args.target}
        n_rows: {len(df)}
        model_path: {model_out}
        R2_test: {r2:.6f}
        features_used ({len(X.columns)}): {', '.join(X.columns)}
        categorical ({len(cat_cols_in_X)}): {', '.join(cat_cols_in_X)}
        bin_quantiles: {json.dumps(parse_kv_list(args.bin_quantiles))}
        shap_rows: {len(X_shap)} (approx={args.shap_approx})
        """).strip()+"\n")

    print(f"[done] Outputs written under: {prefix}*")
    print(f"[done] Model saved to: {model_out}")

if __name__ == "__main__":
    main()
