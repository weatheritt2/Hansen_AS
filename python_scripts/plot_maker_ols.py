# plot_coefficients_from_tsv.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

IN_TSV = "intercomm_ols_HC3.tsv"   
OUTDIR = Path("intercomm_plots")
OUTDIR.mkdir(exist_ok=True, parents=True)

df = pd.read_csv(IN_TSV, sep=None, engine="python")

if "term" in df.columns:
    df = df[df["term"].str.lower() != "const"]

# Basic sanity
need = {"term","coef","se","p","ci_low","ci_high"}
missing = need - set(df.columns)
if missing:
    raise ValueError(f"Coefficient TSV is missing columns: {missing}")

# Clean up term labels (optional: order important terms first)
# Move intercept to the top or drop it
df = df.copy()
df = df[df["term"] != "Intercept"]
# Nice ordering: sort by p ascending
df = df.sort_values("p", ascending=True)

# ---------- Forest plot (coef ± 95% CI) ----------
fig, ax = plt.subplots(figsize=(6, max(2.5, 0.3*len(df))), dpi=150)

ypos = np.arange(len(df))
ax.errorbar(
    x=df["coef"], y=ypos, xerr=[df["coef"]-df["ci_low"], df["ci_high"]-df["coef"]],
    fmt='o', capsize=3, lw=1
)
ax.axvline(0, color="k", ls="--", lw=0.8)
ax.set_yticks(ypos)
ax.set_yticklabels(df["term"])
ax.set_xlabel("Coefficient (OLS, HC3 ± 95% CI)")
ax.set_title("Inter-community OLS: coefficients with 95% CI")
plt.tight_layout()
plt.savefig(OUTDIR / "intercomm_ols_coef_forest.pdf")
plt.close()

# ---------- Bubble plot: -log10(p) vs term; size = |coef|; color = sign(coef) ----------
df["nlog10p"] = -np.log10(df["p"].clip(lower=1e-300))
df["size"] = np.log10(np.abs(df["coef"]).replace(0, np.nan)).replace(-np.inf, np.nan)
smin, smax = 50, 800
size_scaled = (df["size"] - df["size"].min()) / (df["size"].max() - df["size"].min() + 1e-9)
size_scaled = (smin + size_scaled * (smax - smin)).fillna(smin)

# Color and legend
df["sign"] = np.where(df["coef"] >= 0, "Positive", "Negative")
color_map = {"Positive": "tab:red", "Negative": "tab:blue"}
colors = df["sign"].map(color_map)

fig, ax = plt.subplots(figsize=(7, max(2.5, 0.3 * len(df))), dpi=150)
ypos = np.arange(len(df))
sc = ax.scatter(df["nlog10p"], ypos, s=size_scaled, c=colors, alpha=0.7, edgecolors="none")

ax.set_yticks(ypos)
ax.set_yticklabels(df["term"])
ax.set_xlabel(r"$-\log_{10}(p)$")
plt.axvline(x=-np.log10(0.05), color="red", linestyle="--", lw=1.2, label="p = 0.05")
ax.set_title("Inter-community OLS: term significance\n(size ≈ |coef|, color = sign)")

# Legend for colors
import matplotlib.patches as mpatches
legend_patches = [
    mpatches.Patch(color="tab:red", label="Positive coefficient"),
    mpatches.Patch(color="tab:blue", label="Negative coefficient")
]
ax.legend(handles=legend_patches, title="Coefficient sign", loc="upper right", frameon=True)

plt.tight_layout()
plt.savefig(OUTDIR / "intercomm_ols_bubble_nlog10p.pdf")
plt.close()
print(f"[OK] Wrote plots to {OUTDIR.resolve()}")