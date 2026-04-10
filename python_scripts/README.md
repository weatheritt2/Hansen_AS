# README for Kjer-Hansen et al. Alternative splicing fine-tunes modular protein interfaces to control protein localisation


# Script 1: Percolation analysis with degree-matched nulls, κ(f), AUC tests, and inter-community enrichment (v2)

`percolation_auc_plus_v2.py` runs a percolation-style robustness analysis on an interaction network by removing edges in **descending AS priority** (e.g. `as_prob`) and comparing the resulting fragmentation/efficiency curves to **replicate null distributions**. It supports two metrics—largest component fraction **S(f)** and approximate global efficiency **E(f)**—and produces paired **Δ(f) = AS − null** curves with confidence intervals, **AUC-based significance tests**, **κ(f)** (number of connected components), and **inter-community edge removal enrichment** with per-f empirical p-values.

## What this script does

Given an edge list with columns `A`, `B`, and an AS priority score:

1. **Builds a graph** (`networkx.Graph`) from the edge list.
2. **Ranks edges by AS score** (descending; stable sort) and removes the top `k` edges at each fraction `f`.
3. Computes percolation curves for one or both metrics:
   - **S(f):** fraction of nodes in the largest connected component after removing `k` edges
   - **E(f):** approximate global efficiency via sampled BFS sources (`--sources-k`) for speed
4. Computes community structure once on the base graph using greedy modularity communities, then tracks:
   - **# inter-community edges removed** as `f` increases
5. Generates **null percolation distributions** over `--n-nulls` replicates using one of two null modes:
   - `--null-mode global`: random edge removal
   - `--null-mode degree`: degree-matched removal using bins over *min(deg(u), deg(v))* per edge
6. Produces:
   - **Null mean + 95% CI bands** for each curve
   - **Paired Δ curves**: `Δ(f) = AS_curve(f) − null_curve_r(f)` across replicates, with mean and 95% CI
   - **AUC test**: compares AS AUC to the null AUC distribution (empirical one-sided p-value)
   - **κ(f)** curves (AS and null mean/CI)
   - **Inter-community enrichment**: empirical two-sided p-values per `f`

## Requirements

Python packages:
- `numpy`, `pandas`
- `networkx`
- `matplotlib` (for plots)
- Standard library multiprocessing

Install example:
```bash
pip install numpy pandas networkx matplotlib
```

---

Run example:

```bash
python percolation_auc_plus_v2.py \
  --edges interface_interactome.tab \
  --sep tab \
  --col-a A --col-b B --col-as as_prob \
  --metric S \
  --null-mode global \
  --n-nulls 1000 \
  --n-jobs 16 \
  --out-prefix run_plus_v2
```

---



# Script 2: CatBoost + SHAP interaction modeling for AS impact (v3)

`node_feature_models_catboost_with_interactions_v3.py` trains a **CatBoost regression model** to predict alternative splicing (AS) “impact” from a combined feature table, and produces **feature attributions** via **SHAP values** and **SHAP interaction values** (global heatmap + selected dependence plots). It also runs a “safe” **OLS sanity-check** by dummy-coding categorical variables into a numeric design matrix.

## What this script does

Given a TSV containing targets (e.g. `impact_combined`, `deltaS`, `deltaE`) and feature columns, the script:

1. Builds a feature matrix (default panel or user-supplied list)
2. Optionally quantile-bins selected numeric features into categorical bins
3. Trains a CatBoost regressor with correct categorical handling (string tokens, explicit `NA`)
4. Evaluates holdout performance (R² on a test split)
5. Computes permutation importance on the holdout set
6. Runs OLS (dummy-coded categoricals; numeric-only matrix) as an interpretability/sanity check
7. Computes SHAP values and SHAP interaction values and writes summary tables + figures

## Requirements

Python packages:
- `catboost`
- `shap`
- `pandas`, `numpy`
- `matplotlib`, `seaborn`
- `scikit-learn`
- `statsmodels`

Install example:
```bash
pip install catboost shap pandas numpy matplotlib seaborn scikit-learn statsmodels
```

---
Run example:

```bash
python node_feature_models_catboost_with_interactions_v3.py \
  --table impact.tsv \
  --outdir results \
  --prefix as_impact_cb

```

---  

# Script 3: ΔS/ΔE feature modeling pipeline (v1)

`node_feature_models_v3.py` is an end-to-end analysis script for relating **ΔS** and **ΔE** (splicing / expression impact metrics) to molecular and network features from either **edge-level** or **node-level** tables (e.g. `impact.tsv`). It runs **univariate association tests**, **multivariate OLS** (numeric/binary only), and **Random Forest** models (with permutation importance and optional SHAP), plus an optional **joint multi-output RF** predicting ΔS and ΔE simultaneously.

## What this script does

Given an input table containing targets (`deltaS`, `deltaE`) and a set of features, the script produces:

1. **Preflight feature typing**
   - Auto-detects features as **numeric**, **binary**, or **categorical**
   - Optionally treats binary flags as categorical (`--binary-as-categorical`)
2. **Exploratory plots**
   - Spearman correlation heatmap for numeric features + targets
   - Scatter plot of ΔE vs ΔS
3. **Univariate association testing**
   - Numeric features: Spearman correlation vs ΔS / ΔE
   - Binary features: Mann–Whitney U test comparing ΔS / ΔE between groups
   - Multiple testing correction: Benjamini–Hochberg FDR (per target)
4. **OLS regression (sanity check / interpretable baseline)**
   - Uses only numeric + binary features (no categoricals)
   - Z-scores predictors, fits OLS for ΔS and ΔE separately
   - Reports standardized coefficients and FDR-adjusted p-values
5. **Random Forest regression**
   - Separate RF models for ΔS and ΔE
   - Preprocessing:
     - Standardize numeric features
     - One-hot encode categoricals (`handle_unknown="ignore"`)
   - Reports in-sample R² and permutation importance
   - Optional SHAP summary plot (if `shap` is installed and enough rows)
6. **Joint multi-output Random Forest**
   - Fits a MultiOutputRegressor predicting [ΔS, ΔE] jointly
   - Reports per-target R²

## Requirements

Python packages:
- `pandas`, `numpy`, `matplotlib`
- `scipy`
- `statsmodels`
- `scikit-learn`
- Optional: `shap` (for SHAP summary plots)

Install example:
```bash
pip install pandas numpy matplotlib scipy statsmodels scikit-learn shap
```

---
Run example:
```bash
python node_feature_models_v1.py \
  --table impact.tsv \
  --out feature_models
```

---

# Script 4: Single-input enrichment workflow + optional quantile regression

`node_enrichment_qreg_single.py` runs a single-table workflow on `node_impact_with_features.tsv` to test whether molecular/network features are enriched across **ΔE / ΔS** (or other metric) values. For each metric and feature, it produces (i) **quantile-bin trends**, (ii) **tail comparisons** (top vs bottom extremes), and (iii) **permutation-based null Z-scores** for bin-wise enrichment. Optionally, it also fits **quantile regression** (25th/50th/75th percentiles) for numeric features.

## What this script does

For each metric (e.g. `deltae`, `deltas`) and each feature:

1. **Quantile-bin curve**
   - Splits genes/nodes into `--nbins` quantile bins of the metric (rank-based `qcut`)
   - Aggregates the feature within each bin:
     - Numeric features: mean ± SEM
     - Categorical features: fraction of the most common category per bin (simple “dominance” metric)
   - Saves a per-feature bin-curve plot and a TSV summary.

2. **Tail enrichment tests**
   - Defines top and bottom tails using `--tail` (default 10% each)
   - For **numeric** features only:
     - Mann–Whitney U test (two-sided)
     - Cliff’s delta effect size
     - Mean(top) − Mean(bottom)
   - Saves tail boxplots (top/middle/bottom) for all features.

3. **Permutation null for bin curves**
   - Shuffles metric values `--n-perm` times, recomputes bin means, and builds a null distribution
   - Computes per-bin Z-scores:
     - `Z = (real_bin_mean − perm_mean) / perm_sd`
   - Saves per-feature Z-by-bin plots and a TSV.

4. **Optional quantile regression (`--quantreg`)**
   - For numeric features only, fits `feature ~ metric` using quantile regression at q = 0.25, 0.5, 0.75
   - Outputs a summary TSV plus per-feature scatter+fit-line PDFs.

Finally, the script collates “selected” plots into a single `report.pdf`.

## Requirements

Python packages:
- `pandas`, `numpy`, `matplotlib`
- `scipy`
- Optional (for `--quantreg`): `statsmodels`

Install example:
```bash
pip install pandas numpy matplotlib scipy statsmodels
```

---
Run example:

```bash
python node_enrichment_qreg_single.py \
  --in node_impact_with_features.tsv \
  --outdir q_enrich_out
```

---

# Script 5: Plot OLS coefficient summaries from a TSV

`plot_coefficients_from_tsv.py` reads an OLS coefficient summary table (TSV/CSV) and generates two publication-ready visualizations:

1. A **forest plot** showing each term’s coefficient with its **95% confidence interval** (HC3 robust CI assumed).
2. A **bubble plot** showing each term’s **significance** as `-log10(p)` with point **size scaled by |coef|** and **color indicating coefficient sign**.

## What this script does

- Loads a coefficient table from `IN_TSV` using automatic delimiter detection.
- Drops intercept/constant terms (`const` and `Intercept`) if present.
- Validates that required columns are present:
  - `term`, `coef`, `se`, `p`, `ci_low`, `ci_high`
- Sorts terms by ascending p-value (most significant at the top).
- Writes two PDF plots to `OUTDIR`:
  - `intercomm_ols_coef_forest.pdf`
  - `intercomm_ols_bubble_nlog10p.pdf`

## Requirements

Python packages:
- `pandas`
- `numpy`
- `matplotlib`

Install example:
```bash
pip install pandas numpy matplotlib
```

---
Run example:

```bash
python plot_coefficients_from_tsv.py
```

---

# Script 6: Microexon vs long-exon residue and feature enrichment (v2)

`microexon_residue_enrichment_v2.py` performs a **comparative amino-acid composition analysis** between **microexons (MICRO)** and **long exons (LONG)** across exon bodies and flanking regions. It reports **per-residue** and **per-feature (AA class)** enrichments using Fisher’s exact tests, provides **log2 enrichment effect sizes**, and applies **Benjamini–Hochberg FDR correction** within biologically meaningful strata.

---

## What this script does

Given a sequence-context file describing microexon and long-exon events, the script:

1. **Parses exon contexts** into:
   - upstream flank
   - exon body
   - downstream flank  
   with automatic detection of `MICRO` vs `LONG` classes.

2. **Counts amino acids (20 canonical AAs)** in four regions:
   - `surround` : upstream + downstream flanks
   - `first`    : upstream flank only
   - `last`     : downstream flank only
   - `exon`     : exon body

3. **Optionally restricts analysis to a structure-defined subset**  
   (e.g. helix-like regions via a regex such as `H|G`).

4. **Computes residue-level enrichment**:
   - raw counts and frequencies (MICRO vs LONG)
   - `log2(micro / long)` enrichment (with pseudocounts)
   - Fisher’s exact test p-values
   - `-log10(p)` for visualization
   - BH-FDR–adjusted q-values **per (region × subset)**

5. **Computes grouped feature enrichments** for biologically meaningful AA classes:
   - Charged
   - Helix-breaking
   - Helix-forming
   - Aliphatic
   - Aromatic

6. Writes **tidy, plot-ready tables** and a run summary.

---

## Input format

### Required input
A text file (e.g. `microexon_seqs.txt`) where each non-comment line contains:
<upstream_seq> <exon_seq> <downstream_seq> <structure_token> … MICRO|LONG …

## Usage

Basic run example:

```bash
python microexon_residue_enrichment_v2.py \
  --in microexon_seqs.txt \
  --outdir enrich_out
```

---
Script 7: protGPS localisation scoring for AS events + exon/IDR summary outputs (portable)

AS_events_protGPS.py is a takes an alternative splicing (AS) event table and runs protGPS to quantify predicted localisation changes. 

## Requirements

Python packages (minimum practical set):
	•	pandas, numpy
	•	torch
	•	scikit-learn
	•	matplotlib, seaborn
	•	POT (imported as ot)
	•	protGPS (must be importable; script exits cleanly with guidance if missing)

## Inputs

Required:
	•	Ensembl ↔ UniProt mapping table 
	•	UniProt FASTA exported from Ensembl BioMart 
	•	Event table 

## Usage

Basic run example:

```bash
python AS_events_protGPS.py \
  --event-table Gini_data4peptide_full-Hsa_peptides_vCombined.tab \
  --ensembl-uniprot-map GRCh37.p12_ensembl_uniprot.txt \
  --uniprot-fasta uniprot_GRCh38.p14_protsequence.fa.txt \
  --args-pkl checkpoints/protgps/model.args \
  --model-ckpt checkpoints/protgps/model.ckpt \
  --pretrained-hub-dir esm2/esm \
  --output-dir protgps_outputs```

# Script 8: Residue-level interface overlap annotation from multiple PPI interface resources 

interface_interactome_creator.py is a fixed-column, residue-level overlap annotator. Given a query file of canonical UniProt proteins and AA start/end positions, it identifies overlaps against interface annotations from IntAct, 3DID, Interactome3D, ELM, DOMINO, and/or PDB. 

## Data sources (download pages)
	•	IntAct: https://www.ebi.ac.uk/intact/download
	•	3DID: https://3did.irbbarcelona.org/download.php
	•	Interactome3D: https://interactome3d.irbbarcelona.org/download.php
	•	ELM: http://elm.eu.org/
	•	DOMINO: http://mint.bio.uniroma2.it/domino/
	•	PDB: https://www.rcsb.org/downloads

## Requirements

Python packages:
	•	None beyond the standard library (Python ≥3.8 recommended)

## Usage

Basic run example:

```bash
python interface_interactome_creator.py \
  --queries queries.tsv \
  --intact intact.tsv \
  --out-prefix results/intact_only
```