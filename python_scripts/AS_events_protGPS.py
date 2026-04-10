#!/usr/bin/env python3
"""
Portable version of AS_events_GPS_v4_augmented.py.

Required input files
--------------------
1. Ensembl ↔ UniProt mapping table
   Example: GRCh37.p12_ensembl_uniprot.txt

   The script expects tab-delimited columns where:
   - column 1 (0-based index 0): Ensembl protein/transcript identifier used in the event file
   - column 4 (0-based index 3): UniProt accession or alias field; must be non-empty
   - column 5 (0-based index 4): UniProt accession matching the FASTA headers

2. UniProt FASTA exported from Ensembl BioMart
   Example: uniprot_GRCh38.p14_protsequence.fa.txt

   FASTA headers are expected to contain a pipe-delimited UniProt accession, e.g.:
   >sp|P12345|...

3. Event table
   Example: Gini_data4peptide_full-Hsa_peptides_vCombined.tab

   Minimum required columns by 0-based index:
   - 0  : gene name
   - 3  : exon nucleotide length
   - 8  : inclusion flag (expects "TRUE" to keep the event)
   - 9  : Ensembl identifier that can be mapped via the mapping table
   - 17 : amino-acid position string, first entry formatted like "start-end"
   - 27 : IDR annotation (optional semantically, but expected by downstream summary)

protGPS requirement
-------------------
This script requires protGPS to actually score sequences. If protGPS is missing,
the script exits cleanly and explains what to install or where to point PYTHONPATH.

Example
-------
python AS_events_GPS_v4_augmented_portable.py \
  --event-table Gini_data4peptide_full-Hsa_peptides_vCombined.tab \
  --ensembl-uniprot-map GRCh37.p12_ensembl_uniprot.txt \
  --uniprot-fasta uniprot_GRCh38.p14_protsequence.fa.txt \
  --args-pkl checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729.args \
  --model-ckpt checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729epoch=26.ckpt \
  --pretrained-hub-dir esm2/esm \
  --output-dir protgps_outputs
"""

from __future__ import annotations

import argparse
import gc
import os
import pickle
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import seaborn as sns
import torch
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import r2_score, roc_auc_score

try:
    import statsmodels.api as sm
    HAVE_STATSMODELS = True
except Exception:
    HAVE_STATSMODELS = False

# Safe protGPS import. Do not hard-fail during import of this script.
PROTGPS_IMPORT_ERROR = None
try:
    from protgps.utils.loading import get_object
    HAVE_PROTGPS = True
except Exception as exc:
    PROTGPS_IMPORT_ERROR = exc
    HAVE_PROTGPS = False


COMPARTMENT_CLASSES = [
    "nuclear_speckle",
    "p-body",
    "pml-body",
    "post_synaptic_density",
    "stress_granule",
    "chromosome",
    "nucleolus",
    "nuclear_pore_complex",
    "cajal_body",
    "rna_granule",
    "cell_junction",
    "transcriptional",
]

COORDS = np.array([
    [1, 1],
    [5, 1],
    [1, 2],
    [8, 5],
    [5, 2],
    [2, 5],
    [1, 3],
    [3, 3],
    [1, 4],
    [5, 3],
    [8, 6],
    [0, 0],
], dtype=float)

M = ot.dist(COORDS, COORDS, metric="euclidean")
M /= M.max()
M = np.asarray(M, dtype=np.float64)

DEVICE = torch.device("cpu")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict localisation shifts for exon-removal isoforms using protGPS."
    )

    parser.add_argument("--event-table", required=True, help="Path to Gini_data4peptide_full-Hsa_peptides_vCombined.tab or equivalent TSV.")
    parser.add_argument("--ensembl-uniprot-map", required=True, help="Path to GRCh37.p12_ensembl_uniprot.txt or equivalent mapping file.")
    parser.add_argument("--uniprot-fasta", required=True, help="Path to uniprot_GRCh38.p14_protsequence.fa.txt or equivalent FASTA.")
    parser.add_argument("--args-pkl", required=True, help="Path to protGPS .args pickle file.")
    parser.add_argument("--model-ckpt", required=True, help="Path to protGPS checkpoint (.ckpt).")
    parser.add_argument("--pretrained-hub-dir", required=True, help="Path to protGPS/ESM pretrained hub directory.")
    parser.add_argument("--output-dir", required=True, help="Directory for plots, checkpoints, and TSV outputs.")

    parser.add_argument("--thresh", type=float, default=0.83, help="Score threshold used to keep a sequence pair. Default: 0.83")
    parser.add_argument("--pair-chunk", type=int, default=8, help="Number of event pairs to process together. Default: 8")
    parser.add_argument("--model-batch", type=int, default=2, help="Model batch size. Default: 2")
    parser.add_argument("--start-n", type=int, default=16403, help="Start event index. Default preserved from original script.")
    parser.add_argument("--end-n", type=int, default=1000000, help="End event index. Default: 1000000")
    parser.add_argument("--change-threshold", type=float, default=0.10, help="WD threshold defining a localisation-changing event. Default: 0.10")
    parser.add_argument("--idr-col", type=int, default=27, help="0-based index for IDR annotation column in the event table. Default: 27")
    parser.add_argument("--force-overwrite", action="store_true", help="Delete prior outputs/checkpoint in output-dir before running.")

    return parser.parse_args()


def fail(msg: str, exit_code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(exit_code)


def check_file_exists(path_str: str, label: str) -> Path:
    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        fail(f"{label} not found: {path}")
    if not path.is_file():
        fail(f"{label} is not a file: {path}")
    return path


def prepare_output_dir(path_str: str, force_overwrite: bool = False) -> Path:
    outdir = Path(path_str).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if force_overwrite:
        for name in [
            "wd_checkpoint.pkl",
            "wd_impact_table.tsv",
            "wd_scores.tsv",
            "event_level_localisation_idr_position.tsv",
            "localisation_change_idr_summary.tsv",
            "localisation_change_location_summary.tsv",
            "localisation_change_predictive_models.tsv",
        ]:
            p = outdir / name
            if p.exists():
                p.unlink()
    return outdir


def validate_protgps_available() -> None:
    if HAVE_PROTGPS:
        return

    msg = (
        "protGPS could not be imported.\n\n"
        "To run this script you need a working protGPS installation available on PYTHONPATH.\n"
        "Possible fixes:\n"
        "  1. Activate the environment where protGPS is installed.\n"
        "  2. Install protGPS in the current environment.\n"
        "  3. Export PYTHONPATH to include the protGPS package directory.\n\n"
        f"Original import error: {PROTGPS_IMPORT_ERROR}"
    )
    fail(msg)


def validate_event_table_schema(event_table: Path, idr_col: int) -> None:
    required_cols = {0, 3, 8, 9, 17, idr_col}
    with event_table.open() as fh:
        first = fh.readline().rstrip("\n").split("\t")
    ncols = len(first)
    missing = sorted(c for c in required_cols if c >= ncols)
    if missing:
        fail(
            "Event table does not contain the required columns. "
            f"Found {ncols} columns, but need at least up to column index {max(required_cols)}. "
            f"Missing required indices: {missing}."
        )


def load_uniprot_aliases(mapping_file: Path) -> Dict[str, List[str]]:
    uni_alias: Dict[str, List[str]] = {}
    with mapping_file.open() as fh:
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            if len(cols[3]) > 1:
                uni_alias[cols[0]] = [cols[3], cols[4]]

    if not uni_alias:
        fail(
            "No usable Ensembl→UniProt mappings were loaded. "
            "Check that the mapping file is tab-delimited and contains at least 5 columns."
        )
    return uni_alias


def load_uniprot_fasta(fasta_file: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    seq = ""
    current = None

    with fasta_file.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and seq:
                    seqs[current] = seq
                parts = line.split("|")
                if len(parts) < 2:
                    fail(
                        "FASTA header format not recognised. Expected pipe-delimited headers like '>sp|P12345|...'. "
                        f"Problematic header: {line}"
                    )
                current = parts[1]
                seq = ""
            else:
                seq += line.strip()

    if current is not None and seq:
        seqs[current] = seq

    if not seqs:
        fail("No sequences were loaded from the UniProt FASTA file.")
    return seqs


def load_model(snargs: argparse.Namespace):
    model = get_object(snargs.lightning_name, "lightning")(snargs)
    model = model.load_from_checkpoint(
        checkpoint_path=snargs.model_path,
        strict=not snargs.relax_checkpoint_matching,
        **{"args": snargs},
    )
    return model


@torch.no_grad()
def predict_condensates_fast(model, sequences: List[str], batch_size: int = 16, to_cpu: bool = True):
    outs = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        out = model.model({"x": batch})
        s = torch.sigmoid(out["logit"])
        outs.append(s)
    scores = torch.cat(outs, dim=0)
    return scores.cpu() if to_cpu else scores


def get_bio_wasserstein_np(b: np.ndarray, a: np.ndarray, cost_matrix: np.ndarray) -> float:
    bsum = b.sum()
    asum = a.sum()
    if bsum <= 0 or asum <= 0:
        return np.nan
    b = b / bsum
    a = a / asum
    return ot.emd2(b, a, cost_matrix)


def save_checkpoint(path: Path, wd_scores: List[float], transition_matrix: torch.Tensor, last_event_n: int) -> None:
    payload = {
        "wd_scores": wd_scores,
        "transition_matrix": transition_matrix.cpu().numpy(),
        "last_event_n": last_event_n,
    }
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as fh:
        pickle.dump(payload, fh)
    os.replace(tmp, path)


def load_checkpoint(path: Path) -> Tuple[List[float], torch.Tensor, int]:
    if not path.exists():
        return [], torch.zeros((12, 12), dtype=torch.float32), 0
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    return (
        payload.get("wd_scores", []),
        torch.tensor(payload.get("transition_matrix", np.zeros((12, 12))), dtype=torch.float32),
        payload.get("last_event_n", 0),
    )


def append_rows_tsv(path: Path, rows: List[dict]) -> None:
    if not rows:
        return
    df = pd.DataFrame(rows)
    header = not path.exists()
    df.to_csv(path, sep="\t", index=False, mode="a", header=header)


def parse_idr_value(val):
    s = str(val).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "idr", "in_idr", "inside_idr"}:
        return True
    if s in {"0", "false", "f", "no", "n", "non-idr", "ordered", "outside_idr"}:
        return False
    if "idr" in s and "non" not in s and "not" not in s:
        return True
    return np.nan


def norm_location_bin(x):
    if pd.isna(x):
        return np.nan
    x = min(max(float(x), 0.0), 1.0)
    lower = int(np.floor(x * 10.0))
    if lower == 10:
        lower = 9
    upper = lower + 1
    return f"{lower*10}-{upper*10}%"


def safe_div(a, b):
    return np.nan if b == 0 else a / b


def summarise_events(df, change_threshold=0.10):
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    x = df.copy()
    x["changed"] = x["WD"] > change_threshold

    idr_summary = (
        x.groupby("is_idr", dropna=False)
         .agg(
            n_events=("event_n", "count"),
            mean_wd=("WD", "mean"),
            median_wd=("WD", "median"),
            fraction_changed=("changed", "mean"),
            mean_exon_aa_len=("exon_aa_len", "mean"),
            mean_norm_midpoint=("normalized_midpoint", "mean")
         )
         .reset_index()
    )

    loc_summary = (
        x.groupby("location_bin", dropna=False)
         .agg(
            n_events=("event_n", "count"),
            mean_wd=("WD", "mean"),
            median_wd=("WD", "median"),
            fraction_changed=("changed", "mean"),
            idr_fraction=("is_idr_numeric", "mean")
         )
         .reset_index()
    )

    desired_order = [f"{i*10}-{(i+1)*10}%" for i in range(10)]
    loc_summary["location_bin"] = pd.Categorical(loc_summary["location_bin"], categories=desired_order, ordered=True)
    loc_summary = loc_summary.sort_values("location_bin").reset_index(drop=True)

    return idr_summary, loc_summary


def fit_predictive_models(df, out_tsv: Path, change_threshold=0.10):
    if df.empty:
        pd.DataFrame(columns=["model", "term", "coefficient", "pvalue", "metric", "metric_value", "n"]).to_csv(out_tsv, sep="\t", index=False)
        return

    x = df.copy()
    x = x.dropna(subset=["WD", "normalized_midpoint", "is_idr_numeric"])
    if x.empty:
        pd.DataFrame(columns=["model", "term", "coefficient", "pvalue", "metric", "metric_value", "n"]).to_csv(out_tsv, sep="\t", index=False)
        return

    x["changed"] = (x["WD"] > change_threshold).astype(int)
    X = x[["normalized_midpoint", "is_idr_numeric"]].astype(float)
    X = X.rename(columns={"normalized_midpoint": "norm_pos", "is_idr_numeric": "is_idr"})
    y_bin = x["changed"].astype(int)
    y_cont = x["WD"].astype(float)

    rows = []

    if HAVE_STATSMODELS and len(np.unique(y_bin)) > 1:
        try:
            X_sm = sm.add_constant(X, has_constant="add")
            logit = sm.Logit(y_bin, X_sm).fit(disp=False)
            pred = logit.predict(X_sm)
            auc = roc_auc_score(y_bin, pred) if len(np.unique(y_bin)) > 1 else np.nan
            for term in logit.params.index:
                rows.append({
                    "model": "logistic_changed",
                    "term": term,
                    "coefficient": logit.params[term],
                    "pvalue": logit.pvalues.get(term, np.nan),
                    "metric": "roc_auc",
                    "metric_value": auc,
                    "n": len(x),
                })
        except Exception:
            pass

        try:
            X_sm = sm.add_constant(X, has_constant="add")
            ols = sm.OLS(y_cont, X_sm).fit()
            for term in ols.params.index:
                rows.append({
                    "model": "linear_wd",
                    "term": term,
                    "coefficient": ols.params[term],
                    "pvalue": ols.pvalues.get(term, np.nan),
                    "metric": "r_squared",
                    "metric_value": ols.rsquared,
                    "n": len(x),
                })
        except Exception:
            pass

    if not rows:
        if len(np.unique(y_bin)) > 1:
            clf = LogisticRegression(max_iter=1000)
            clf.fit(X, y_bin)
            pred = clf.predict_proba(X)[:, 1]
            auc = roc_auc_score(y_bin, pred)
            rows.append({
                "model": "logistic_changed",
                "term": "intercept",
                "coefficient": float(clf.intercept_[0]),
                "pvalue": np.nan,
                "metric": "roc_auc",
                "metric_value": auc,
                "n": len(x),
            })
            for term, coef in zip(X.columns, clf.coef_[0]):
                rows.append({
                    "model": "logistic_changed",
                    "term": term,
                    "coefficient": float(coef),
                    "pvalue": np.nan,
                    "metric": "roc_auc",
                    "metric_value": auc,
                    "n": len(x),
                })

        reg = LinearRegression()
        reg.fit(X, y_cont)
        pred_cont = reg.predict(X)
        r2 = r2_score(y_cont, pred_cont)
        rows.append({
            "model": "linear_wd",
            "term": "intercept",
            "coefficient": float(reg.intercept_),
            "pvalue": np.nan,
            "metric": "r_squared",
            "metric_value": r2,
            "n": len(x),
        })
        for term, coef in zip(X.columns, reg.coef_):
            rows.append({
                "model": "linear_wd",
                "term": term,
                "coefficient": float(coef),
                "pvalue": np.nan,
                "metric": "r_squared",
                "metric_value": r2,
                "n": len(x),
            })

    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)


def plot_augmented_outputs(df_events, idr_summary, loc_summary, plot_dir: Path, change_threshold=0.10):
    if df_events.empty:
        return

    x = df_events.copy()
    x["changed_label"] = np.where(x["WD"] > change_threshold, "Changed", "Unchanged")

    plt.figure(figsize=(6, 5))
    box_df = x.dropna(subset=["is_idr", "WD"]).copy()
    if not box_df.empty:
        box_df["IDR_state"] = np.where(box_df["is_idr"], "IDR", "Non-IDR")
        sns.boxplot(data=box_df, x="IDR_state", y="WD")
        plt.title("Localisation shift by IDR status")
        plt.xlabel("")
        plt.ylabel("Biological Wasserstein Distance")
        plt.tight_layout()
        plt.savefig(plot_dir / "4_WD_by_IDR_boxplot.pdf")
        plt.close()
    else:
        plt.close()

    plt.figure(figsize=(7, 5))
    idr_bar = idr_summary.copy()
    if not idr_bar.empty:
        idr_bar["IDR_state"] = idr_bar["is_idr"].map({True: "IDR", False: "Non-IDR"}).fillna("Unknown")
        sns.barplot(data=idr_bar, x="IDR_state", y="fraction_changed")
        plt.ylim(0, 1)
        plt.ylabel(f"Fraction changed (WD > {change_threshold})")
        plt.xlabel("")
        plt.title("Fraction of localisation-changing exons by IDR status")
        plt.tight_layout()
        plt.savefig(plot_dir / "5_fraction_changed_by_IDR.pdf")
        plt.close()
    else:
        plt.close()

    plt.figure(figsize=(10, 5))
    loc_box = x.dropna(subset=["location_bin", "WD"]).copy()
    if not loc_box.empty:
        desired_order = [f"{i*10}-{(i+1)*10}%" for i in range(10)]
        sns.boxplot(data=loc_box, x="location_bin", y="WD", order=desired_order)
        plt.xticks(rotation=45, ha="right")
        plt.xlabel("Normalised exon midpoint in protein")
        plt.ylabel("Biological Wasserstein Distance")
        plt.title("Localisation shift by exon position along protein")
        plt.tight_layout()
        plt.savefig(plot_dir / "6_WD_by_location_bin_boxplot.pdf")
        plt.close()
    else:
        plt.close()

    plt.figure(figsize=(10, 5))
    if not loc_summary.empty:
        desired_order = [f"{i*10}-{(i+1)*10}%" for i in range(10)]
        sns.barplot(data=loc_summary, x="location_bin", y="fraction_changed", order=desired_order)
        plt.ylim(0, 1)
        plt.xticks(rotation=45, ha="right")
        plt.xlabel("Normalised exon midpoint in protein")
        plt.ylabel(f"Fraction changed (WD > {change_threshold})")
        plt.title("Fraction of localisation-changing exons by protein position")
        plt.tight_layout()
        plt.savefig(plot_dir / "7_fraction_changed_by_location_bin.pdf")
        plt.close()
    else:
        plt.close()

    plt.figure(figsize=(7, 5))
    scatter_df = x.dropna(subset=["normalized_midpoint", "WD"]).copy()
    if not scatter_df.empty:
        sns.scatterplot(data=scatter_df, x="normalized_midpoint", y="WD", hue="changed_label", alpha=0.7)
        plt.xlabel("Normalised exon midpoint in protein")
        plt.ylabel("Biological Wasserstein Distance")
        plt.title("Continuous relationship between exon position and localisation shift")
        plt.legend(title="")
        plt.tight_layout()
        plt.savefig(plot_dir / "8_WD_vs_normalized_position_scatter.pdf")
        plt.close()
    else:
        plt.close()


def iter_candidate_pairs(
    infile: Path,
    uni_alias: Dict[str, List[str]],
    uni_seqs: Dict[str, str],
    resume_after_n: int = 0,
    start_n: int = 16403,
    end_n: int = 1000000,
    idr_col: int = 27,
) -> Iterator[dict]:
    n = 0.0
    with infile.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n").split("\t")

            if len(line) <= max(17, 9, 8, 3, idr_col):
                continue

            if line[8] != "TRUE":
                continue

            raw_pos = line[17].split(",")[0]
            if "-" not in raw_pos:
                continue

            start_s, end_s = raw_pos.split("-", 1)

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            try:
                ensembl_id = line[9].rstrip()
                uniprot_acc = uni_alias[ensembl_id][1]
                full_length = uni_seqs[uniprot_acc]
            except KeyError:
                continue

            try:
                exon_nt_len = int(line[3])
            except ValueError:
                continue

            if (exon_nt_len / 3) * 0.75 > abs(end - start):
                continue

            n += 1.0
            if n % 100 == 0:
                print(f"screened {int(n)}")

            if n < start_n:
                continue
            if n > end_n:
                break
            if n <= resume_after_n:
                continue

            isoform = full_length[:start] + full_length[end:]
            protein_length = len(full_length)
            exon_aa_len = abs(end - start)
            midpoint_aa = (start + end) / 2.0
            normalized_midpoint = safe_div(midpoint_aa, protein_length)
            is_idr = parse_idr_value(line[idr_col])

            yield {
                "full_length": full_length,
                "isoform": isoform,
                "event_n": int(n),
                "gene_name": line[0],
                "ensembl_id": ensembl_id,
                "raw_position": raw_pos,
                "exon_nt_len": pd.to_numeric(line[3], errors="coerce"),
                "exon_aa_len": exon_aa_len,
                "protein_length": protein_length,
                "exon_start_aa": start,
                "exon_end_aa": end,
                "midpoint_aa": midpoint_aa,
                "normalized_midpoint": normalized_midpoint,
                "location_bin": norm_location_bin(normalized_midpoint),
                "is_idr": is_idr,
                "idr_raw": line[idr_col],
            }


def chunked(iterator: Iterable[dict], size: int) -> Iterator[List[dict]]:
    chunk = []
    for item in iterator:
        chunk.append(item)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def main() -> None:
    args = parse_args()

    validate_protgps_available()

    event_table = check_file_exists(args.event_table, "Event table")
    mapping_file = check_file_exists(args.ensembl_uniprot_map, "Ensembl-UniProt mapping file")
    fasta_file = check_file_exists(args.uniprot_fasta, "UniProt FASTA file")
    args_pkl = check_file_exists(args.args_pkl, "protGPS args pickle")
    model_ckpt = check_file_exists(args.model_ckpt, "protGPS checkpoint")
    pretrained_hub_dir = Path(args.pretrained_hub_dir).expanduser().resolve()
    if not pretrained_hub_dir.exists():
        fail(f"pretrained-hub-dir not found: {pretrained_hub_dir}")

    validate_event_table_schema(event_table, args.idr_col)
    outdir = prepare_output_dir(args.output_dir, force_overwrite=args.force_overwrite)

    checkpoint_path = outdir / "wd_checkpoint.pkl"
    impact_tsv = outdir / "wd_impact_table.tsv"
    wd_tsv = outdir / "wd_scores.tsv"
    event_tsv = outdir / "event_level_localisation_idr_position.tsv"
    idr_summary_tsv = outdir / "localisation_change_idr_summary.tsv"
    location_summary_tsv = outdir / "localisation_change_location_summary.tsv"
    model_tsv = outdir / "localisation_change_predictive_models.tsv"

    print("device:", DEVICE)
    print("Loading mapping file...")
    uni_alias = load_uniprot_aliases(mapping_file)
    print(f"Loaded {len(uni_alias):,} Ensembl→UniProt mappings")

    print("Loading FASTA sequences...")
    uni_seqs = load_uniprot_fasta(fasta_file)
    print(f"Loaded {len(uni_seqs):,} UniProt sequences")

    print("Loading protGPS args/checkpoint...")
    snargs = argparse.Namespace(**pickle.load(args_pkl.open("rb")))
    snargs.model_path = str(model_ckpt)
    snargs.pretrained_hub_dir = str(pretrained_hub_dir)

    model = load_model(snargs)
    model.eval()
    model = model.to(DEVICE)

    wd_scores, transition_matrix, resume_after_n = load_checkpoint(checkpoint_path)
    print("Resuming after event_n =", resume_after_n)

    labels = COMPARTMENT_CLASSES

    for pair_chunk in chunked(
        iter_candidate_pairs(
            infile=event_table,
            uni_alias=uni_alias,
            uni_seqs=uni_seqs,
            resume_after_n=resume_after_n,
            start_n=args.start_n,
            end_n=args.end_n,
            idr_col=args.idr_col,
        ),
        args.pair_chunk,
    ):
        flat_sequences = []
        for item in pair_chunk:
            flat_sequences.extend([item["full_length"], item["isoform"]])

        scores = predict_condensates_fast(model, flat_sequences, batch_size=args.model_batch, to_cpu=True)
        scores_np = scores.numpy().astype(np.float64)

        chunk_impact_rows = []
        chunk_wd_rows = []
        chunk_event_rows = []

        for i, item in enumerate(pair_chunk):
            b = scores_np[2 * i]
            a = scores_np[2 * i + 1]

            keep_pair = (b > args.thresh).any() or (a > args.thresh).any()
            if not keep_pair:
                continue

            wd = get_bio_wasserstein_np(b, a, M)
            if np.isnan(wd):
                continue

            wd_scores.append(wd)
            chunk_wd_rows.append({
                "event_n": item["event_n"],
                "WD": wd,
            })

            origin_idx = int(np.argmax(b))
            origin_label = labels[origin_idx]
            isoform_idx = int(np.argmax(a))
            isoform_label = labels[isoform_idx]

            chunk_impact_rows.append({
                "Compartment": origin_label,
                "WD": wd,
                "event_n": item["event_n"],
            })

            chunk_event_rows.append({
                "event_n": item["event_n"],
                "gene_name": item["gene_name"],
                "ensembl_id": item["ensembl_id"],
                "WD": wd,
                "origin_compartment": origin_label,
                "isoform_top_compartment": isoform_label,
                "max_full_score": float(np.max(b)),
                "max_isoform_score": float(np.max(a)),
                "raw_position": item["raw_position"],
                "exon_nt_len": item["exon_nt_len"],
                "exon_aa_len": item["exon_aa_len"],
                "protein_length": item["protein_length"],
                "exon_start_aa": item["exon_start_aa"],
                "exon_end_aa": item["exon_end_aa"],
                "midpoint_aa": item["midpoint_aa"],
                "normalized_midpoint": item["normalized_midpoint"],
                "location_bin": item["location_bin"],
                "is_idr": item["is_idr"],
                "is_idr_numeric": np.nan if pd.isna(item["is_idr"]) else int(bool(item["is_idr"])),
                "idr_raw": item["idr_raw"],
                "full_length_seq": item["full_length"],
                "isoform_seq": item["isoform"],
            })

            diff = a - b
            transition_matrix[origin_idx] += torch.from_numpy(diff).float()

        append_rows_tsv(impact_tsv, chunk_impact_rows)
        append_rows_tsv(wd_tsv, chunk_wd_rows)
        append_rows_tsv(event_tsv, chunk_event_rows)

        last_event_n = pair_chunk[-1]["event_n"]
        save_checkpoint(checkpoint_path, wd_scores, transition_matrix, last_event_n)
        print(f"checkpoint saved at event {last_event_n}")

        del scores, scores_np, chunk_impact_rows, chunk_wd_rows, chunk_event_rows, flat_sequences
        gc.collect()

    wd_scores_np = np.array(wd_scores)

    if impact_tsv.exists():
        df_impact = pd.read_csv(impact_tsv, sep="\t")
    else:
        df_impact = pd.DataFrame(columns=["Compartment", "WD", "event_n"])

    if event_tsv.exists():
        df_events = pd.read_csv(event_tsv, sep="\t")
    else:
        df_events = pd.DataFrame()

    if len(wd_scores_np) == 0 or df_impact.empty:
        print("No results passed filtering; no plots generated.")
        return

    counts = df_impact["Compartment"].value_counts().reindex(labels).fillna(0).values
    avg_transition = transition_matrix / (torch.tensor(counts).view(-1, 1) + 1e-9)

    threshold = args.change_threshold
    significant_exons = wd_scores_np > threshold
    fraction_changed = np.mean(significant_exons)

    print(f"Fraction of exons causing localisation change: {fraction_changed:.2%}")

    plt.figure(figsize=(8, 6))
    sns.ecdfplot(wd_scores_np, complementary=True, linewidth=2)
    plt.axvline(x=threshold, color="r", linestyle="--", label=f"Threshold ({threshold})")
    plt.text(threshold + 0.02, 0.5, f"{fraction_changed:.1%} of Exons", color="r", fontweight="bold")
    plt.title("Distribution of Localisation Shift Magnitudes (WD)")
    plt.xlabel("Biological Wasserstein Distance")
    plt.ylabel("Fraction of Exons with Change ≥ X")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "exon_change_fraction_cdf.pdf")
    plt.close()

    plt.figure(figsize=(10, 6))
    order = df_impact.groupby("Compartment")["WD"].mean().sort_values(ascending=False).index
    sns.barplot(data=df_impact, x="WD", y="Compartment", order=order, palette="magma")
    plt.title("Sensitivity to Exon Removal by Compartment")
    plt.tight_layout()
    plt.savefig(outdir / "2_Compartment_Sensitivity.pdf")
    plt.close()

    plt.figure(figsize=(12, 10))
    df_heatmap = pd.DataFrame(avg_transition.numpy(), index=labels, columns=labels)
    sns.heatmap(df_heatmap, cmap="PiYG", center=0, annot=False)
    plt.title("Average Score Redistribution per Compartment")
    plt.tight_layout()
    plt.savefig(outdir / "3_Shift_Direction_Heatmap.pdf")
    plt.close()

    if not df_events.empty:
        df_events["is_idr"] = df_events["is_idr"].replace({"True": True, "False": False})
        if "is_idr_numeric" not in df_events.columns:
            df_events["is_idr_numeric"] = df_events["is_idr"].map({True: 1, False: 0})
        idr_summary, loc_summary = summarise_events(df_events, change_threshold=args.change_threshold)
        idr_summary.to_csv(idr_summary_tsv, sep="\t", index=False)
        loc_summary.to_csv(location_summary_tsv, sep="\t", index=False)
        fit_predictive_models(df_events, model_tsv, change_threshold=args.change_threshold)
        plot_augmented_outputs(df_events, idr_summary, loc_summary, outdir, change_threshold=args.change_threshold)

    print("Done.")
    print(f"Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
