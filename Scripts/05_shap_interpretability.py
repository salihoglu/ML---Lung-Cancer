#!/usr/bin/env python3
"""
RS
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import sys
import traceback
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import importlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from scipy.special import expit
from scipy.stats import chi2, kendalltau, rankdata, spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_ipcw
from sksurv.util import Surv
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")

SEED = 42
FIG_DPI = 600
RNG = np.random.default_rng(SEED)
np.random.seed(SEED)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": True,
    "legend.framealpha": 0.95,
    "legend.edgecolor": "0.8",
    "legend.fontsize": 9,
    "figure.dpi": 120,
    "savefig.dpi": FIG_DPI,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

PAL = {
    "Driver": "#D85A30",
    "Suppressor": "#378ADD",
    "Marker/Other": "#888780",
    "high": "#B2182B",
    "low": "#2166AC",
    "neutral": "#4D4D4D",
}

FERRO_DRIVERS = {
    "ACSL4", "LPCAT3", "NOX1", "NOX4", "HMOX1", "TFRC", "SLC40A1", "SLC11A2",
    "FTH1", "FTL", "STEAP3", "PCBP1", "ALOX5", "ALOX12", "ALOX15", "PTGS2", "TP53",
    "LHFPL2", "RIPK2", "FYN", "FAS", "NCOA4", "SLC1A5",
}
FERRO_SUPPRESSORS = {
    "GPX4", "SLC7A11", "SLC3A2", "FSP1", "DHODH", "GCH1", "NFE2L2", "KEAP1",
    "HSPA5", "PRDX6", "B4GALT1", "OGT", "TXNRD1", "AIFM2",
}


@dataclass
class CohortResult:
    cohort: str
    n_samples: int
    n_events: int
    genes_used: List[str]
    final_score: np.ndarray
    final_source: str
    risk_column: Optional[str]
    Xs: pd.DataFrame
    surv: pd.DataFrame
    tau: float
    perm_df: pd.DataFrame
    shap_df: pd.DataFrame
    combined_df: pd.DataFrame
    cox_df: pd.DataFrame
    ipcw_cindex: float
    cindex_ci_low: float
    cindex_ci_high: float
    surrogate_r2: float
    surrogate_rho: float
    surrogate_r2_insample: float
    surrogate_rho_insample: float
    surrogate_r2_cv: float
    surrogate_rho_cv: float
    shap_direction_df: pd.DataFrame


def msg(text: str) -> None:
    print(text, flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Q1-grade multi-cohort external SHAP interpretability analysis")
    p.add_argument("--manifest", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--cohorts", nargs="*", default=["GSE30219", "GSE50081", "GSE68465"])
    p.add_argument("--n-perm", type=int, default=200)
    p.add_argument("--n-boot", type=int, default=300)
    p.add_argument("--n-estimators", type=int, default=500)
    p.add_argument("--surrogate-estimators", type=int, default=400)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--min-overlap-samples", type=int, default=30)
    p.add_argument("--min-events", type=int, default=10)
    p.add_argument("--min-overlap-genes", type=int, default=5)
    p.add_argument("--cox-top-n", type=int, default=30)
    p.add_argument("--consensus-rank-weight-perm", type=float, default=0.6)
    p.add_argument("--consensus-rank-weight-shap", type=float, default=0.4)
    p.add_argument("--manuscript-top-n", type=int, default=20, help="Number of genes shown in manuscript-ready interpretability plots")
    p.add_argument("--support-rho-threshold", type=float, default=0.15, help="Minimum |Spearman rho(feature, SHAP)| considered directional SHAP support")
    return p.parse_args()


def gene_cat(g: str) -> str:
    if g in FERRO_DRIVERS:
        return "Driver"
    if g in FERRO_SUPPRESSORS:
        return "Suppressor"
    return "Marker/Other"


def sanitize_gene(x: object) -> str:
    x = str(x).strip().replace(" ", "")
    if "__" in x:
        x = x.split("__")[-1]
    return x.replace(".", "-").upper()


def clean_id(vals: Iterable[object]) -> pd.Index:
    s = pd.Series(vals, dtype="string").astype(str).str.strip()
    s = s.str.replace("_", "-", regex=False)
    s = s.str.replace(r"\.", "-", regex=True)
    s = s.str.replace(r"\s+", "", regex=True)
    s = s.str.replace(r"-01A.*$|-01B.*$|-10A.*$|-10B.*$|-11A.*$|-11B.*$", "", regex=True)
    s = s.str.replace(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$", r"\1", regex=True)
    s = s.str.replace(r"^(GSM\d+).*$", r"\1", regex=True)
    s = s.str.replace(r"^(\d+)[A-Za-z]$", r"\1", regex=True)
    return pd.Index(s)


def read_flex_csv(path: str | Path, index_col=0) -> pd.DataFrame:
    try:
        return pd.read_csv(path, index_col=index_col)
    except Exception:
        try:
            return pd.read_csv(path, index_col=index_col, engine="python")
        except Exception:
            return pd.read_csv(path, index_col=index_col, sep=None, engine="python")


def read_manifest(path: str | Path) -> List[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def resolve_from_manifest(manifest_paths: Sequence[str], patterns: Sequence[str]) -> Optional[str]:
    for pat in patterns:
        rgx = re.compile(pat, flags=re.IGNORECASE)
        hits = [p for p in manifest_paths if rgx.search(p)]
        if hits:
            return hits[0]
    return None


def read_gene_list(path: str | Path) -> List[str]:
    df = read_flex_csv(path, index_col=None)
    if df is None or df.shape[0] == 0:
        return []
    if df.shape[1] == 1:
        genes = [sanitize_gene(x) for x in df.iloc[:, 0].dropna().astype(str)]
        return pd.Series([g for g in genes if g and g != "NAN"]).drop_duplicates().tolist()
    cols_low = {str(c).lower(): c for c in df.columns}
    for key in ["gene", "symbol", "gene_symbol", "genesymbol", "hgnc_symbol", "feature"]:
        if key in cols_low:
            gene_col = cols_low[key]
            break
    else:
        gene_col = next((c for c in df.columns if any(k in str(c).lower() for k in ["gene", "symbol", "feature"])), df.columns[0])
    genes = [sanitize_gene(x) for x in df[gene_col].dropna().astype(str)]
    return pd.Series([g for g in genes if g and g != "NAN"]).drop_duplicates().tolist()



def _first_matching_col(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    low = {str(c).lower().strip(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    for c in df.columns:
        lc = str(c).lower()
        if any(cand.lower() in lc for cand in candidates):
            return c
    return None


def infer_ml_direction_from_row(row: pd.Series) -> str:
    """Infer ML-derived direction as Risk / Protective / Unknown from common gene-list columns."""
    for key in ["direction", "Direction", "effect_direction", "risk_direction", "ml_direction"]:
        if key in row.index and pd.notna(row[key]):
            val = str(row[key]).strip().lower()
            if any(x in val for x in ["risk", "high", "positive", "hazard", "up"]):
                return "Risk"
            if any(x in val for x in ["protect", "low", "negative", "benefit", "down"]):
                return "Protective"
    for key in ["zmeta", "z_meta", "meta_z", "z", "coef", "beta", "ml_beta"]:
        if key in row.index:
            val = pd.to_numeric(pd.Series([row[key]]), errors="coerce").iloc[0]
            if np.isfinite(val):
                return "Risk" if val > 0 else "Protective"
    return "Unknown"


def read_gene_metadata(path: str | Path) -> pd.DataFrame:
    """Read gene list plus optional ML metadata used for manuscript-facing support tables."""
    df = read_flex_csv(path, index_col=None)
    if df is None or df.shape[0] == 0:
        return pd.DataFrame(columns=["gene", "ml_direction", "ml_z", "ml_bag_frac", "ml_source"])
    if df.shape[1] == 1:
        out = pd.DataFrame({"gene": [sanitize_gene(x) for x in df.iloc[:, 0].dropna().astype(str)]})
    else:
        gene_col = _first_matching_col(df, ["gene", "symbol", "gene_symbol", "genesymbol", "hgnc_symbol", "feature"])
        if gene_col is None:
            gene_col = df.columns[0]
        out = pd.DataFrame({"gene": [sanitize_gene(x) for x in df[gene_col].astype(str)]})
        out["ml_direction"] = df.apply(infer_ml_direction_from_row, axis=1).values
        z_col = _first_matching_col(df, ["zmeta", "z_meta", "meta_z", "z"])
        beta_col = _first_matching_col(df, ["coef", "beta", "ml_beta"])
        if z_col is not None:
            out["ml_z"] = pd.to_numeric(df[z_col], errors="coerce").values
        elif beta_col is not None:
            out["ml_z"] = pd.to_numeric(df[beta_col], errors="coerce").values
        else:
            out["ml_z"] = np.nan
        bag_col = _first_matching_col(df, ["bag_frac", "bagging_fraction", "selection_frequency", "stability", "loco_stability"])
        out["ml_bag_frac"] = pd.to_numeric(df[bag_col], errors="coerce").values if bag_col is not None else np.nan
        source_col = _first_matching_col(df, ["source", "Source", "gene_source", "annotation"])
        out["ml_source"] = df[source_col].astype(str).values if source_col is not None else ""
    out = out.loc[out["gene"].notna() & (out["gene"] != "") & (out["gene"] != "NAN")].copy()
    if "ml_direction" not in out.columns:
        out["ml_direction"] = "Unknown"
    if "ml_z" not in out.columns:
        out["ml_z"] = np.nan
    if "ml_bag_frac" not in out.columns:
        out["ml_bag_frac"] = np.nan
    if "ml_source" not in out.columns:
        out["ml_source"] = ""
    out["category"] = out["gene"].map(gene_cat)
    out = out.drop_duplicates("gene", keep="first").reset_index(drop=True)
    return out[["gene", "ml_direction", "ml_z", "ml_bag_frac", "ml_source", "category"]]


def ml_direction_sign(direction: object, z: object = np.nan) -> float:
    direction = str(direction).strip().lower()
    if any(x in direction for x in ["risk", "positive", "high"]):
        return 1.0
    if any(x in direction for x in ["protect", "negative", "low"]):
        return -1.0
    val = pd.to_numeric(pd.Series([z]), errors="coerce").iloc[0]
    if np.isfinite(val):
        return 1.0 if val > 0 else -1.0
    return np.nan


def classify_directional_support(shap_rho: float, expected_sign: float, threshold: float) -> str:
    if not np.isfinite(shap_rho) or not np.isfinite(expected_sign):
        return "Not assessable"
    if abs(shap_rho) < threshold:
        return "Weak/neutral"
    return "Concordant" if np.sign(shap_rho) == np.sign(expected_sign) else "Discordant"


def compute_shap_direction_table(shap_values: np.ndarray, X_df: pd.DataFrame, cohort: str,
                                 gene_meta: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Summarise whether high expression shifts surrogate SHAP in the expected ML direction."""
    rows = []
    meta = gene_meta.set_index("gene") if gene_meta is not None and not gene_meta.empty else pd.DataFrame().set_index(pd.Index([]))
    for j, gene in enumerate(X_df.columns):
        fv = np.asarray(X_df[gene].values, dtype=float)
        sv = np.asarray(shap_values[:, j], dtype=float)
        ok = np.isfinite(fv) & np.isfinite(sv)
        if ok.sum() < 10:
            rho = slope = delta_q4_q1 = np.nan
        else:
            rho = safe_spearman(fv[ok], sv[ok])
            try:
                slope = float(np.polyfit(fv[ok], sv[ok], 1)[0])
            except Exception:
                slope = np.nan
            q1 = np.nanpercentile(fv[ok], 25)
            q3 = np.nanpercentile(fv[ok], 75)
            low = sv[ok & (fv <= q1)]
            high = sv[ok & (fv >= q3)]
            delta_q4_q1 = float(np.nanmean(high) - np.nanmean(low)) if len(low) and len(high) else np.nan
        if gene in meta.index:
            md = meta.loc[gene, "ml_direction"]
            mz = meta.loc[gene, "ml_z"]
            mb = meta.loc[gene, "ml_bag_frac"]
            ms = meta.loc[gene, "ml_source"]
        else:
            md, mz, mb, ms = "Unknown", np.nan, np.nan, ""
        exp_sign = ml_direction_sign(md, mz)
        rows.append({
            "cohort": cohort,
            "gene": gene,
            "ml_direction": md,
            "ml_z": mz,
            "ml_bag_frac": mb,
            "ml_source": ms,
            "category": gene_cat(gene),
            "expected_risk_sign": exp_sign,
            "shap_feature_rho": rho,
            "shap_feature_slope": slope,
            "shap_delta_top_bottom_quartile": delta_q4_q1,
            "shap_direction": "Risk-shifting" if np.isfinite(rho) and rho > threshold else ("Protective-shifting" if np.isfinite(rho) and rho < -threshold else "Weak/neutral"),
            "directional_support": classify_directional_support(rho, exp_sign, threshold),
        })
    return pd.DataFrame(rows)


def aggregate_directional_support(direction_all: pd.DataFrame, gene_meta: pd.DataFrame) -> pd.DataFrame:
    if direction_all is None or direction_all.empty:
        return pd.DataFrame()
    def frac_eq(x, val):
        x = pd.Series(x).dropna()
        return float(np.mean(x == val)) if len(x) else np.nan
    agg = direction_all.groupby("gene", as_index=False).agg(
        n_cohorts_direction=("cohort", "nunique"),
        mean_shap_feature_rho=("shap_feature_rho", "mean"),
        median_shap_feature_rho=("shap_feature_rho", "median"),
        mean_shap_delta_top_bottom_quartile=("shap_delta_top_bottom_quartile", "mean"),
        concordant_fraction=("directional_support", lambda x: frac_eq(x, "Concordant")),
        discordant_fraction=("directional_support", lambda x: frac_eq(x, "Discordant")),
        weak_fraction=("directional_support", lambda x: frac_eq(x, "Weak/neutral")),
    )
    meta = gene_meta.drop_duplicates("gene") if gene_meta is not None and not gene_meta.empty else pd.DataFrame({"gene": agg["gene"]})
    out = meta.merge(agg, on="gene", how="right")
    out["category"] = out["gene"].map(gene_cat)
    return out


def load_expr(path: str | Path) -> pd.DataFrame:
    df = read_flex_csv(path, index_col=0)
    df.index = clean_id(df.index.to_series())
    df.columns = [c.split("__")[-1] if "__" in c else c for c in df.columns]
    df.columns = [sanitize_gene(c) for c in df.columns]
    if pd.Index(df.columns).duplicated().any():
        df = df.T.groupby(level=0, sort=False).mean(numeric_only=False).T
    for c in list(df.columns):
        col = df.loc[:, c]
        if isinstance(col, pd.DataFrame):
            col = col.iloc[:, 0]
        df.loc[:, c] = pd.to_numeric(col, errors="coerce")
    df = df.loc[:, ~df.isna().all(axis=0)]
    if pd.Index(df.index).duplicated().any():
        df = df.groupby(level=0, sort=False).mean(numeric_only=True)
    try:
        arr = df.to_numpy(dtype=float, copy=False)
        q99 = np.nanquantile(arr, 0.99)
        if np.isfinite(q99) and q99 > 30:
            df = np.log2(df.astype(float) + 1.0)
            msg(f"[INFO] Applied log2(x+1) transform to {Path(path).parent.name} expression")
    except Exception:
        pass
    return df.astype(float)


def detect_time_event_cols(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
    low = {c.lower(): c for c in df.columns}
    time_pats = ["os_time", "os.time", "survival_time", "futime", "time"]
    event_pats = ["os_event", "os_status", "event", "status", "fustat", "os"]
    time_col = next((low[c] for c in time_pats if c in low), None)
    event_col = next((low[c] for c in event_pats if c in low), None)
    if time_col is None:
        time_col = next((c for c in df.columns if "time" in c.lower()), None)
    if event_col is None:
        event_col = next((c for c in df.columns if "event" in c.lower() or "status" in c.lower()), None)
    return time_col, event_col


def normalize_event(x: pd.Series) -> np.ndarray:
    if pd.api.types.is_numeric_dtype(x):
        z = pd.to_numeric(x, errors="coerce")
        return np.where(pd.isna(z), np.nan, np.where(z > 0, 1, 0)).astype(float)
    z = x.astype(str).str.strip().str.upper()
    out = np.where(z.str.contains("DEAD|EVENT|TRUE|YES|^1$", regex=True), 1,
                   np.where(z.str.contains("ALIVE|FALSE|NO|CENSOR|^0$", regex=True), 0, np.nan))
    return out.astype(float)


def normalize_time(x: pd.Series) -> np.ndarray:
    z = pd.to_numeric(x, errors="coerce").astype(float)
    med = np.nanmedian(z)
    if np.isfinite(med) and med > 100:
        return z / 365.25
    if np.isfinite(med) and med > 15:
        return z / 12.0
    return z


def load_surv(path: str | Path) -> pd.DataFrame:
    df = read_flex_csv(path, index_col=0)
    df.index = clean_id(df.index.to_series())
    time_col, event_col = detect_time_event_cols(df)
    if time_col is None or event_col is None:
        raise RuntimeError(f"Could not detect time/event columns in {path}")
    out = pd.DataFrame(index=df.index)
    out["T"] = normalize_time(df[time_col])
    out["E"] = normalize_event(df[event_col])
    mappings = [("age", "age"), ("sex", "sex"), ("gender", "sex"), ("stage", "stage"), ("ajcc_pathologic_stage", "stage")]
    for src, dest in mappings:
        if src in df.columns and dest not in out.columns:
            out[dest] = df[src].values
    if "sex" in out.columns:
        out["sex"] = out["sex"].astype(str).str.strip().str.title()
    if "stage" in out.columns:
        out["stage"] = out["stage"].astype(str)
    out = out.loc[np.isfinite(out["T"]) & pd.notna(out["E"])]
    return out


def resolve_risk_path(manifest_paths: Sequence[str], cohort: str) -> Optional[str]:
    pats = [
        rf"/OUT_v44/risks/test_{cohort}__risk__STACKING\.csv$",
        rf"/OUT_v44/risks/test_{cohort}__risk__.*\.csv$",
        rf"/nomogram_results_v3/{cohort}/analysis_df\.csv$",
        rf"/{cohort}/analysis_df\.csv$",
        rf"/{cohort}/.*risk.*\.csv$",
    ]
    return resolve_from_manifest(manifest_paths, pats)


def detect_risk_col(df: pd.DataFrame) -> Optional[str]:
    for c in ["risk_score", "RiskScore", "score", "linear_predictor", "lp"]:
        if c in df.columns:
            return c
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    exclude = {"t", "e", "time", "event", "status", "surv", "os"}
    cands = [c for c in num_cols if c.lower() not in exclude and not any(e in c.lower() for e in exclude)]
    return cands[0] if cands else (num_cols[0] if num_cols else None)


def load_risk(path: str | Path) -> Tuple[Optional[pd.Series], Optional[str]]:
    df = read_flex_csv(path, index_col=0)
    df.index = clean_id(df.index.to_series())
    rc = detect_risk_col(df)
    if rc is None:
        return None, None
    s = pd.to_numeric(df[rc], errors="coerce")
    return s, rc


def safe_div(a: float, b: float) -> float:
    return float(a / b) if b not in (0, 0.0) and np.isfinite(b) else np.nan


def bh_adjust(p: Sequence[float]) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty_like(adj)
    restored[order] = adj
    out[ok] = restored
    return out


def ipcw_concordance(T: np.ndarray, E: np.ndarray, scores: np.ndarray, tau: Optional[float] = None) -> float:
    if tau is None:
        tau = np.percentile(T[E == 1], 75) if (E == 1).sum() > 0 else np.nanmax(T) * 0.9
    tau = min(float(tau), float(np.nanmax(T)) * 0.99)
    y_struct = Surv.from_arrays(event=E.astype(bool), time=T)
    try:
        return float(concordance_index_ipcw(y_struct, y_struct, estimate=-scores, tau=tau)[0])
    except Exception:
        return float(concordance_index(T, -scores, E))


def bootstrap_metric(T: np.ndarray, E: np.ndarray, scores: np.ndarray, n_boot: int = 300, tau: Optional[float] = None,
                     seed: int = SEED) -> Tuple[float, float, np.ndarray]:
    rng = np.random.default_rng(seed)
    vals = []
    n = len(T)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            vals.append(ipcw_concordance(T[idx], E[idx], scores[idx], tau=tau))
        except Exception:
            continue
    arr = np.asarray(vals, dtype=float)
    if arr.size < 20:
        return np.nan, np.nan, arr
    return float(np.nanpercentile(arr, 2.5)), float(np.nanpercentile(arr, 97.5)), arr


def permutation_importance_ipcw(X_df: pd.DataFrame, T: np.ndarray, E: np.ndarray, model, n_perm: int = 200,
                                tau: Optional[float] = None, seed: int = SEED) -> Tuple[pd.DataFrame, float]:
    rng = np.random.default_rng(seed)
    base_scores = model.predict(X_df.values.astype(np.float32))
    base_ci = ipcw_concordance(T, E, base_scores, tau=tau)
    Xp = X_df.copy()
    rows = []
    for gene in X_df.columns:
        orig = X_df[gene].values.copy()
        drops = []
        for _ in range(n_perm):
            Xp[gene] = rng.permutation(orig)
            perm_scores = model.predict(Xp.values.astype(np.float32))
            perm_ci = ipcw_concordance(T, E, perm_scores, tau=tau)
            drops.append(base_ci - perm_ci)
        Xp[gene] = orig
        arr = np.asarray(drops, dtype=float)
        rows.append({
            "gene": gene,
            "ci_drop": float(np.nanmean(arr)),
            "ci_drop_std": float(np.nanstd(arr, ddof=1)) if arr.size > 1 else np.nan,
            "ci_drop_ci_low": float(np.nanpercentile(arr, 2.5)) if arr.size > 5 else np.nan,
            "ci_drop_ci_high": float(np.nanpercentile(arr, 97.5)) if arr.size > 5 else np.nan,
        })
    out = pd.DataFrame(rows).sort_values("ci_drop", ascending=False).reset_index(drop=True)
    out["rank_perm"] = np.arange(1, len(out) + 1)
    out["category"] = out["gene"].map(gene_cat)
    return out, base_ci



def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho with stable NA handling."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def fit_risk_surrogate_with_cv(X_df: pd.DataFrame, final_score: np.ndarray, n_estimators: int,
                               seed: int = SEED, n_splits: int = 5) -> Tuple[ExtraTreesRegressor, Dict[str, float], np.ndarray]:
    """Fit the risk-score surrogate and compute in-sample and CV approximation diagnostics.

    The CV diagnostics are the reviewer-facing fidelity estimates. In-sample metrics are retained
    only to document surrogate capacity and should not be interpreted as independent validation.
    """
    X = X_df.values.astype(np.float32)
    y = np.asarray(final_score, dtype=float)
    n = X.shape[0]
    base = ExtraTreesRegressor(n_estimators=n_estimators, n_jobs=-1, random_state=seed)
    full_model = base.fit(X, y)
    pred_in = full_model.predict(X)
    metrics = {
        "r2_insample": float(r2_score(y, pred_in)) if np.isfinite(y).all() else np.nan,
        "rho_insample": safe_spearman(y, pred_in),
        "r2_cv": np.nan,
        "rho_cv": np.nan,
    }
    if n >= 30:
        k = min(n_splits, max(2, n // 15))
        cv = KFold(n_splits=k, shuffle=True, random_state=seed)
        try:
            pred_cv = cross_val_predict(base, X, y, cv=cv, n_jobs=-1)
            metrics["r2_cv"] = float(r2_score(y, pred_cv))
            metrics["rho_cv"] = safe_spearman(y, pred_cv)
        except Exception:
            pred_cv = np.full_like(y, np.nan, dtype=float)
    else:
        pred_cv = np.full_like(y, np.nan, dtype=float)
    return full_model, metrics, pred_cv


def permutation_importance_risk_surrogate(X_df: pd.DataFrame, final_score: np.ndarray, model,
                                          n_perm: int = 200, seed: int = SEED) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Permutation importance for the same target used by SHAP: risk-score surrogate output.

    Importance is quantified as the loss of surrogate approximation quality after permuting one
    feature. This makes the permutation ranking target-aligned with SHAP, unlike the separate
    survival-refit IPCW permutation sensitivity analysis.
    """
    rng = np.random.default_rng(seed)
    X = X_df.values.astype(np.float32)
    y = np.asarray(final_score, dtype=float)
    base_pred = model.predict(X)
    base_r2 = float(r2_score(y, base_pred))
    base_rho = safe_spearman(y, base_pred)
    Xp = X_df.copy()
    rows = []
    for gene in X_df.columns:
        orig = Xp[gene].values.copy()
        r2_drops, rho_drops = [], []
        for _ in range(n_perm):
            Xp[gene] = rng.permutation(orig)
            pred = model.predict(Xp.values.astype(np.float32))
            r2_drops.append(base_r2 - float(r2_score(y, pred)))
            rho_perm = safe_spearman(y, pred)
            rho_drops.append(base_rho - rho_perm if np.isfinite(base_rho) and np.isfinite(rho_perm) else np.nan)
        Xp[gene] = orig
        arr = np.asarray(r2_drops, dtype=float)
        arr_rho = np.asarray(rho_drops, dtype=float)
        rows.append({
            "gene": gene,
            "r2_drop": float(np.nanmean(arr)),
            "r2_drop_std": float(np.nanstd(arr, ddof=1)) if arr.size > 1 else np.nan,
            "r2_drop_ci_low": float(np.nanpercentile(arr, 2.5)) if arr.size > 5 else np.nan,
            "r2_drop_ci_high": float(np.nanpercentile(arr, 97.5)) if arr.size > 5 else np.nan,
            "rho_drop": float(np.nanmean(arr_rho)),
            "rho_drop_ci_low": float(np.nanpercentile(arr_rho, 2.5)) if arr_rho.size > 5 else np.nan,
            "rho_drop_ci_high": float(np.nanpercentile(arr_rho, 97.5)) if arr_rho.size > 5 else np.nan,
        })
    out = pd.DataFrame(rows).sort_values("r2_drop", ascending=False).reset_index(drop=True)
    out["rank_perm"] = np.arange(1, len(out) + 1)
    out["category"] = out["gene"].map(gene_cat)
    return out, {"base_r2_insample": base_r2, "base_rho_insample": base_rho}


def multivariate_cox_pergene(X_df: pd.DataFrame, surv_df: pd.DataFrame, genes: Sequence[str], top_n: int = 30) -> pd.DataFrame:
    records = []
    covariates = [c for c in ["stage", "age", "sex"] if c in surv_df.columns]
    for gene in list(genes)[:top_n]:
        try:
            df = pd.DataFrame({"T": surv_df["T"], "E": surv_df["E"].astype(int), gene: X_df[gene]})
            for c in covariates:
                df[c] = surv_df[c].values
            df = df.dropna()
            if df.shape[0] < 30 or df["E"].sum() < 5:
                continue
            cph = CoxPHFitter(penalizer=0.1)
            cph.fit(df, duration_col="T", event_col="E", show_progress=False)
            row = cph.summary.loc[gene]
            records.append({
                "gene": gene,
                "HR": float(np.exp(row["coef"])),
                "HR_lower": float(np.exp(row["coef lower 95%"])),
                "HR_upper": float(np.exp(row["coef upper 95%"])),
                "coef": float(row["coef"]),
                "coef_se": float(row["se(coef)"]),
                "p": float(row["p"]),
                "z": float(row["z"]),
                "covariates_adjusted": ";".join(covariates) if covariates else "none",
                "category": gene_cat(gene),
            })
        except Exception:
            continue
    out = pd.DataFrame(records)
    if out.empty:
        return out
    out = out.sort_values("p").reset_index(drop=True)
    out["p_adj"] = bh_adjust(out["p"].values)
    return out


def random_effects_meta(effects: Sequence[float], ses: Sequence[float]) -> Dict[str, float]:
    y = np.asarray(effects, dtype=float)
    se = np.asarray(ses, dtype=float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0)
    y = y[ok]
    se = se[ok]
    if y.size == 0:
        return {k: np.nan for k in ["k", "fixed_beta", "fixed_se", "random_beta", "random_se", "z", "p", "Q", "Q_p", "I2", "tau2"]}
    v = se ** 2
    w = 1.0 / v
    fixed_beta = np.sum(w * y) / np.sum(w)
    fixed_se = math.sqrt(1.0 / np.sum(w))
    Q = float(np.sum(w * (y - fixed_beta) ** 2))
    df = y.size - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / c) if y.size > 1 and c > 0 else 0.0
    wr = 1.0 / (v + tau2)
    random_beta = np.sum(wr * y) / np.sum(wr)
    random_se = math.sqrt(1.0 / np.sum(wr))
    z = random_beta / random_se if random_se > 0 else np.nan
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2)))) if np.isfinite(z) else np.nan
    Q_p = 1 - chi2.cdf(Q, df=df) if df > 0 else np.nan
    I2 = max(0.0, (Q - df) / Q) if Q > 0 and df > 0 else 0.0
    return {
        "k": int(y.size),
        "fixed_beta": float(fixed_beta),
        "fixed_se": float(fixed_se),
        "random_beta": float(random_beta),
        "random_se": float(random_se),
        "z": float(z) if np.isfinite(z) else np.nan,
        "p": float(p) if np.isfinite(p) else np.nan,
        "Q": float(Q),
        "Q_p": float(Q_p) if np.isfinite(Q_p) else np.nan,
        "I2": float(I2),
        "tau2": float(tau2),
    }


def kendalls_w(rankings_matrix: np.ndarray) -> Tuple[float, float]:
    m, n = rankings_matrix.shape
    if m < 2 or n < 2:
        return np.nan, np.nan
    ranked = np.apply_along_axis(rankdata, 1, rankings_matrix)
    Ri = ranked.sum(axis=0)
    S = np.sum((Ri - Ri.mean()) ** 2)
    W = 12 * S / (m ** 2 * (n ** 3 - n))
    chi2_stat = m * (n - 1) * W
    p = 1 - chi2.cdf(chi2_stat, df=n - 1)
    return float(W), float(p)


def save_fig(fig, outdir: Path, stem: str) -> None:
    outdir.mkdir(exist_ok=True, parents=True)
    for ext in ["pdf", "tiff", "png"]:
        try:
            fig.savefig(outdir / f"{stem}.{ext}", dpi=FIG_DPI, bbox_inches="tight")
        except Exception:
            pass


def plot_permutation_bar(perm_df: pd.DataFrame, outdir: Path, cohort: str, cindex: float,
                         ci_low: float, ci_high: float, top_n: int = 20) -> None:
    top = perm_df.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(9, max(5, top_n * 0.38)))
    y = np.arange(len(top))
    colors = top["category"].map(PAL).fillna(PAL["Marker/Other"])
    err_low = top["ci_drop"] - top["ci_drop_ci_low"]
    err_high = top["ci_drop_ci_high"] - top["ci_drop"]
    ax.barh(y, top["ci_drop"], xerr=np.vstack([err_low.fillna(0), err_high.fillna(0)]),
            color=colors, alpha=0.88, capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Permutation IPCW C-index drop")
    ttl = f"{cohort}: permutation importance\nIPCW C-index {cindex:.3f}"
    if np.isfinite(ci_low) and np.isfinite(ci_high):
        ttl += f" [{ci_low:.3f}–{ci_high:.3f}]"
    ax.set_title(ttl)
    ax.legend(handles=[mpatches.Patch(color=PAL[k], label=k) for k in ["Driver", "Suppressor", "Marker/Other"]], loc="lower right")
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_permutation_importance_q1")
    plt.close(fig)



def plot_risk_surrogate_permutation_bar(perm_df: pd.DataFrame, outdir: Path, cohort: str,
                                        r2_cv: float, rho_cv: float, top_n: int = 20) -> None:
    top = perm_df.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(9, max(5, top_n * 0.38)))
    y = np.arange(len(top))
    colors = top["category"].map(PAL).fillna(PAL["Marker/Other"])
    err_low = top["r2_drop"] - top["r2_drop_ci_low"]
    err_high = top["r2_drop_ci_high"] - top["r2_drop"]
    ax.barh(y, top["r2_drop"], xerr=np.vstack([err_low.fillna(0), err_high.fillna(0)]),
            color=colors, alpha=0.88, capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Permutation drop in surrogate R²")
    ttl = f"{cohort}: risk-score surrogate permutation importance"
    if np.isfinite(r2_cv) or np.isfinite(rho_cv):
        ttl += f"\nCV fidelity: R²={r2_cv:.3f}, ρ={rho_cv:.3f}"
    ax.set_title(ttl)
    ax.legend(handles=[mpatches.Patch(color=PAL[k], label=k) for k in ["Driver", "Suppressor", "Marker/Other"]], loc="lower right")
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_risk_surrogate_permutation_importance_q1")
    plt.close(fig)


def plot_shap_beeswarm(shap_values: np.ndarray, X_df: pd.DataFrame, top_features: Sequence[str], outdir: Path, cohort: str) -> None:
    ordered = [g for g in top_features if g in X_df.columns]
    if not ordered:
        return
    fig, ax = plt.subplots(figsize=(10, max(6, len(ordered) * 0.42)))
    for i, gene in enumerate(reversed(ordered)):
        j = list(X_df.columns).index(gene)
        sv = np.asarray(shap_values[:, j], dtype=float)
        fv = np.asarray(X_df[gene].values, dtype=float)
        if np.nanmax(np.abs(sv)) == 0:
            continue
        jitter = RNG.normal(0, 0.08, size=len(sv))
        cmap = plt.cm.RdBu_r
        norm = Normalize(vmin=np.nanpercentile(fv, 5), vmax=np.nanpercentile(fv, 95))
        ax.scatter(sv, np.full_like(sv, i, dtype=float) + jitter, c=cmap(norm(fv)), s=10, alpha=0.6, rasterized=True)
        x_text = np.nanmin(sv) - 0.03 * (np.nanmax(sv) - np.nanmin(sv) + 1e-8)
        ax.text(x_text, i, gene, ha="right", va="center", fontsize=8,
                color=PAL.get(gene_cat(gene), "black"), fontweight="bold")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels([""] * len(ordered))
    ax.set_xlabel("Surrogate SHAP value")
    ax.set_title(f"{cohort}: surrogate-model SHAP")
    sm = ScalarMappable(cmap=plt.cm.RdBu_r, norm=Normalize(0, 1))
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, pad=0.01, aspect=40)
    cb.set_label("Feature value\n(low → high)")
    ax.legend(handles=[mpatches.Patch(color=PAL[k], label=k) for k in ["Driver", "Suppressor", "Marker/Other"]], loc="lower right")
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_surrogate_shap_beeswarm_q1")
    plt.close(fig)


def plot_cox_forest(cox_df: pd.DataFrame, outdir: Path, cohort: str, top_n: int = 20) -> None:
    if cox_df is None or cox_df.empty:
        return
    df = cox_df.head(top_n).sort_values("HR")
    fig, ax = plt.subplots(figsize=(8, max(5, len(df) * 0.4)))
    y = np.arange(len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        color = PAL[row["category"]]
        ax.scatter(row["HR"], i, color=color, s=60, zorder=3)
        ax.plot([row["HR_lower"], row["HR_upper"]], [i, i], color=color, lw=2)
        if row.get("p_adj", 1) < 0.05:
            ax.text(row["HR_upper"] + 0.02, i, "*", va="center", color=color, fontsize=12)
    ax.axvline(1, color="black", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title(f"{cohort}: adjusted per-gene Cox")
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_multivariable_cox_forest_q1")
    plt.close(fig)


def fit_score_recalibration_cox(T: np.ndarray, E: np.ndarray, risk_scores: np.ndarray):
    df = pd.DataFrame({"T": T, "E": E.astype(int), "score": risk_scores}).replace([np.inf, -np.inf], np.nan).dropna()
    if df.shape[0] < 30 or df["E"].sum() < 10:
        return None
    try:
        cph = CoxPHFitter(penalizer=0.05)
        cph.fit(df, duration_col="T", event_col="E", show_progress=False)
        return cph
    except Exception:
        return None


def predict_survival_from_recalibration(cph, risk_scores: np.ndarray, times: Sequence[float]) -> Optional[pd.DataFrame]:
    if cph is None:
        return None
    X = pd.DataFrame({"score": risk_scores})
    try:
        sf = cph.predict_survival_function(X, times=list(times))
        # lifelines returns times x individuals
        if isinstance(sf, pd.DataFrame):
            return sf.T
    except Exception:
        return None
    return None


def bootstrap_calibration_slope(prob_event: np.ndarray, y_event: np.ndarray, n_boot: int = 200, seed: int = SEED):
    rng = np.random.default_rng(seed)
    vals = []
    n = len(prob_event)
    prob_event = np.clip(prob_event, 1e-6, 1 - 1e-6)
    logit = np.log(prob_event / (1 - prob_event))
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yy = y_event[idx]
        if len(np.unique(yy)) < 2:
            continue
        try:
            fit = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
            fit.fit(logit[idx].reshape(-1,1), yy)
            vals.append(float(fit.coef_[0][0]))
        except Exception:
            continue
    if len(vals) < 20:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


def plot_calibration(T: np.ndarray, E: np.ndarray, risk_scores: np.ndarray, outdir: Path, cohort: str, years=(1, 3, 5)) -> Optional[pd.DataFrame]:
    cph = fit_score_recalibration_cox(T, E, risk_scores)
    if cph is None:
        return None
    times = [yr * 365.25 if np.nanmax(T) > 500 else yr for yr in years]
    pred_surv_df = predict_survival_from_recalibration(cph, risk_scores, times)
    if pred_surv_df is None:
        return None

    fig, axes = plt.subplots(1, len(years), figsize=(5 * len(years), 5))
    if len(years) == 1:
        axes = [axes]
    rows = []
    for ax, yr, t_eval in zip(axes, years, times):
        pred_surv = pred_surv_df[t_eval].to_numpy(dtype=float)
        pred_event = np.clip(1 - pred_surv, 1e-6, 1 - 1e-6)
        observed_event = ((T <= t_eval) & (E == 1)).astype(int)
        try:
            deciles = pd.qcut(pred_event, q=10, labels=False, duplicates="drop")
        except Exception:
            ax.set_title(f"{yr}-year (insufficient spread)")
            continue
        xs, ys, ns = [], [], []
        for dec in sorted(pd.Series(deciles).dropna().unique()):
            mask = np.asarray(deciles == dec)
            if np.sum(mask) < 8:
                continue
            kmf = KaplanMeierFitter()
            kmf.fit(T[mask], E[mask], label="_")
            obs_surv = float(kmf.survival_function_at_times([t_eval]).values[0])
            xs.append(float(np.mean(pred_surv[mask])))
            ys.append(obs_surv)
            ns.append(int(np.sum(mask)))
        if len(xs) < 3:
            ax.set_title(f"{yr}-year (insufficient data)")
            continue
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        logit = np.log(pred_event / (1 - pred_event))
        slope_ci = (np.nan, np.nan)
        slope = np.nan
        intercept = np.nan
        if len(np.unique(observed_event)) >= 2:
            try:
                lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
                lr.fit(logit.reshape(-1,1), observed_event)
                slope = float(lr.coef_[0][0])
                intercept = float(lr.intercept_[0])
                slope_ci = bootstrap_calibration_slope(pred_event, observed_event, n_boot=200, seed=SEED)
            except Exception:
                pass
        brier = float(np.mean((pred_event - observed_event) ** 2))
        ax.scatter(xs, ys, color=PAL["high"], s=np.asarray(ns) * 6, alpha=0.8)
        lim = [max(0, min(xs.min(), ys.min()) - 0.05), min(1.05, max(xs.max(), ys.max()) + 0.05)]
        ax.plot(lim, lim, "k--", lw=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel(f"Predicted {yr}-year survival")
        ax.set_ylabel(f"Observed {yr}-year survival")
        ax.set_title(f"{yr}-year calibration\nBrier={brier:.3f}; slope={slope:.2f}")
        rows.append({
            "cohort": cohort, "year": yr, "time_eval": t_eval,
            "brier": brier, "calibration_slope": slope, "calibration_intercept": intercept,
            "slope_ci_lower": slope_ci[0], "slope_ci_upper": slope_ci[1],
            "n_bins": len(xs), "n_samples": len(T), "n_events_by_time": int(observed_event.sum())
        })
    plt.suptitle(f"{cohort}: score-recalibrated calibration curves", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_calibration_curves_q1")
    plt.close(fig)
    return pd.DataFrame(rows)


def plot_dca(T: np.ndarray, E: np.ndarray, risk_scores: np.ndarray, outdir: Path, cohort: str, t_yr: int = 5) -> Optional[pd.DataFrame]:
    t_threshold = t_yr * 365.25 if np.nanmax(T) > 500 else t_yr
    died_before = ((T <= t_threshold) & (E == 1)).astype(int)
    if died_before.sum() < 10:
        return None
    cph = fit_score_recalibration_cox(T, E, risk_scores)
    if cph is None:
        return None
    pred_surv_df = predict_survival_from_recalibration(cph, risk_scores, [t_threshold])
    if pred_surv_df is None:
        return None
    probs = np.clip(1 - pred_surv_df[t_threshold].to_numpy(dtype=float), 1e-6, 1 - 1e-6)
    thresholds = np.linspace(0.05, 0.80, 76)
    n = len(T)
    prev = died_before.mean()
    nb_model, nb_all = [], []
    for thr in thresholds:
        pos = probs >= thr
        tp = np.sum(pos & (died_before == 1))
        fp = np.sum(pos & (died_before == 0))
        nb_model.append(tp / n - fp / n * (thr / (1 - thr + 1e-8)))
        nb_all.append(prev - (1 - prev) * (thr / (1 - thr + 1e-8)))
    dca_df = pd.DataFrame({"threshold": thresholds, "net_benefit_model": nb_model, "net_benefit_all": nb_all})
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(thresholds, nb_model, color=PAL["high"], lw=2, label="Risk model")
    ax.plot(thresholds, nb_all, color="grey", lw=1.5, ls="--", label="Treat all")
    ax.axhline(0, color="black", lw=1, ls=":", label="Treat none")
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.set_title(f"{cohort}: score-recalibrated decision curve ({t_yr}-year mortality)")
    ax.set_xlim(0.05, 0.80)
    ax.legend()
    plt.tight_layout()
    save_fig(fig, outdir, f"{cohort}_decision_curve_{t_yr}yr_q1")
    plt.close(fig)
    dca_df.insert(0, "cohort", cohort)
    dca_df.insert(1, "year", t_yr)
    return dca_df


def plot_rank_heatmap(rank_df: pd.DataFrame, cohorts: Sequence[str], outdir: Path, stem: str, title: str, top_n: int = 20) -> None:
    top = rank_df.head(top_n).copy()
    if top.empty:
        return
    heat = np.vstack([top[f"rank_{c}"].values for c in cohorts])
    fig, ax = plt.subplots(figsize=(max(8, top_n // 2), len(cohorts) + 2))
    im = ax.imshow(heat, aspect="auto", cmap="RdYlGn_r")
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top["gene"], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels(cohorts)
    plt.colorbar(im, ax=ax, label="Rank (lower = more important)")
    ax.set_title(title)
    plt.tight_layout()
    save_fig(fig, outdir, stem)
    plt.close(fig)


def plot_consensus_bar(consensus_df: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    top = consensus_df.head(top_n).copy()
    if top.empty:
        return
    if "consensus_importance_score" not in top.columns:
        max_rank = float(np.nanmax(consensus_df["combined_mean_rank"])) if "combined_mean_rank" in consensus_df.columns else float(top_n)
        top["consensus_importance_score"] = max_rank - top["combined_mean_rank"] + 1.0
    fig, ax = plt.subplots(figsize=(9, max(5, top_n * 0.35)))
    y = np.arange(len(top))
    colors = top["category"].map(PAL).fillna(PAL["Marker/Other"])
    ax.barh(y, top["consensus_importance_score"], color=colors, alpha=0.88)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Consensus importance score (higher = more stable/important)")
    ax.set_title("Consensus ML-gene support across external cohorts")
    ax.legend(handles=[mpatches.Patch(color=PAL[k], label=k) for k in ["Driver", "Suppressor", "Marker/Other"]], loc="lower right")
    plt.tight_layout()
    save_fig(fig, outdir, "consensus_top_genes_bar_q1")
    plt.close(fig)



def plot_manuscript_ml_gene_support(summary_df: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    """Manuscript-facing figure: ML-prioritised genes supported by cross-cohort SHAP/rank evidence."""
    if summary_df is None or summary_df.empty:
        return
    df = summary_df.sort_values(["support_tier_order", "combined_mean_rank"], ascending=[True, True]).head(top_n).copy()
    if df.empty:
        return
    if "consensus_importance_score" not in df.columns:
        max_rank = float(np.nanmax(summary_df["combined_mean_rank"]))
        df["consensus_importance_score"] = max_rank - df["combined_mean_rank"] + 1.0
    df = df.sort_values("consensus_importance_score")
    fig, ax = plt.subplots(figsize=(9, max(5.5, len(df) * 0.38)))
    y = np.arange(len(df))
    dir_colors = {"Risk": "#D85A30", "Protective": "#2166AC", "Unknown": "#888780"}
    colors = df["ml_direction"].map(dir_colors).fillna("#888780")
    ax.barh(y, df["consensus_importance_score"], color=colors, alpha=0.88)
    for i, (_, row) in enumerate(df.iterrows()):
        mark = "✓" if row.get("shap_support_tier", "") in ["Strong SHAP-concordant", "Moderate SHAP-concordant"] else ("!" if row.get("shap_support_tier", "") == "Discordant" else "")
        if mark:
            ax.text(row["consensus_importance_score"] + 0.15, i, mark, va="center", ha="left", fontsize=11, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Consensus importance score (higher = stronger cross-cohort support)")
    ax.set_title("Manuscript-ready support for ML-prioritised ferroptosis-related genes")
    handles = [mpatches.Patch(color=dir_colors[k], label=k) for k in ["Risk", "Protective", "Unknown"]]
    ax.legend(handles=handles, loc="lower right", title="ML direction")
    plt.tight_layout()
    save_fig(fig, outdir, "manuscript_ml_gene_support_summary_q1")
    plt.close(fig)


def plot_directional_shap_support(summary_df: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    """Plot mean feature-SHAP correlation by ML direction; positive supports risk genes, negative supports protective genes."""
    if summary_df is None or summary_df.empty or "mean_shap_feature_rho" not in summary_df.columns:
        return
    df = summary_df.sort_values(["support_tier_order", "combined_mean_rank"], ascending=[True, True]).head(top_n).copy()
    df = df.sort_values("mean_shap_feature_rho")
    fig, ax = plt.subplots(figsize=(8.5, max(5, len(df) * 0.38)))
    y = np.arange(len(df))
    dir_colors = {"Risk": "#D85A30", "Protective": "#2166AC", "Unknown": "#888780"}
    colors = df["ml_direction"].map(dir_colors).fillna("#888780")
    ax.scatter(df["mean_shap_feature_rho"], y, c=colors, s=70, zorder=3)
    ax.hlines(y, 0, df["mean_shap_feature_rho"], color=colors, lw=2, alpha=0.75)
    ax.axvline(0, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Mean Spearman correlation: feature value vs surrogate SHAP")
    ax.set_title("Directional SHAP support for ML-derived risk/protective direction")
    ax.legend(handles=[mpatches.Patch(color=dir_colors[k], label=k) for k in ["Risk", "Protective", "Unknown"]], loc="lower right")
    plt.tight_layout()
    save_fig(fig, outdir, "manuscript_directional_shap_support_q1")
    plt.close(fig)


def write_manuscript_interpretability_note(summary_df: pd.DataFrame, fidelity_df: pd.DataFrame,
                                           kendall_W: float, kendall_p: float, outdir: Path,
                                           top_n: int = 10) -> None:
    if summary_df is None or summary_df.empty:
        return
    top = summary_df.sort_values(["support_tier_order", "combined_mean_rank"], ascending=[True, True]).head(top_n)
    genes = ", ".join(top["gene"].astype(str).tolist())
    cv_r2 = fidelity_df["surrogate_r2_cv"].dropna()
    cv_rho = fidelity_df["surrogate_spearman_rho_cv"].dropna()
    n_strong = int((summary_df["shap_support_tier"] == "Strong SHAP-concordant").sum()) if "shap_support_tier" in summary_df else 0
    n_mod = int((summary_df["shap_support_tier"] == "Moderate SHAP-concordant").sum()) if "shap_support_tier" in summary_df else 0
    lines = [
        "Suggested manuscript text for interpretability section:",
        "",
        (
            f"Cross-cohort interpretability analysis showed moderate feature-ranking concordance "
            f"across the external cohorts (Kendall's W={kendall_W:.3f}, p={kendall_p:.4g}), "
            f"indicating that the ML-derived risk structure was reproducible but not identical across datasets."
        ),
        (
            f"The most consistently supported ML-prioritised genes were {genes}. "
            f"Directional SHAP analysis classified {n_strong} genes as strongly and {n_mod} genes as moderately "
            f"concordant with their ML-derived risk/protective direction."
        ),
        (
            "SHAP values were computed from an ExtraTrees surrogate trained to approximate the final cohort-specific "
            "risk score; therefore, they explain the risk-score approximation rather than causal gene effects or the "
            "original survival learner directly."
        ),
    ]
    if len(cv_r2) > 0 and len(cv_rho) > 0:
        lines.insert(2, f"Cross-validated surrogate fidelity: median R2={cv_r2.median():.3f}, median Spearman rho={cv_rho.median():.3f}.")
    (outdir / "manuscript_interpretability_results_text_q1.txt").write_text("\n".join(lines), encoding="utf-8")


def plot_meta_importance(meta_df: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    top = meta_df.head(top_n).copy()
    if top.empty:
        return
    top = top.sort_values("random_beta")
    fig, ax = plt.subplots(figsize=(8.5, max(5, len(top) * 0.38)))
    y = np.arange(len(top))
    colors = top["category"].map(PAL).fillna(PAL["Marker/Other"])
    ax.scatter(top["random_beta"], y, c=colors, s=60)
    ax.hlines(y, top["random_beta"] - 1.96 * top["random_se"], top["random_beta"] + 1.96 * top["random_se"], color=colors)
    ax.axvline(0, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlabel("Random-effects meta importance (rank-z)")
    ax.set_title("Cross-cohort meta-analysis of feature importance")
    plt.tight_layout()
    save_fig(fig, outdir, "meta_feature_importance_forest_q1")
    plt.close(fig)


def plot_meta_cox(meta_df: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    top = meta_df.head(top_n).copy()
    if top.empty:
        return
    top = top.sort_values("random_beta")
    fig, ax = plt.subplots(figsize=(8.5, max(5, len(top) * 0.38)))
    y = np.arange(len(top))
    colors = top["category"].map(PAL).fillna(PAL["Marker/Other"])
    hr = np.exp(top["random_beta"])
    low = np.exp(top["random_beta"] - 1.96 * top["random_se"])
    high = np.exp(top["random_beta"] + 1.96 * top["random_se"])
    ax.scatter(hr, y, c=colors, s=60)
    ax.hlines(y, low, high, color=colors)
    ax.axvline(1, color="black", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlabel("Random-effects hazard ratio")
    ax.set_title("Cross-cohort meta-analysis of adjusted Cox effects")
    plt.tight_layout()
    save_fig(fig, outdir, "meta_adjusted_cox_forest_q1")
    plt.close(fig)


def z_from_rank(rank: pd.Series) -> pd.Series:
    x = pd.to_numeric(rank, errors="coerce")
    if x.notna().sum() < 2:
        return pd.Series(np.nan, index=rank.index)
    mu = x.mean()
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=rank.index)
    return -(x - mu) / sd


def leave_one_out_meta(rank_table: pd.DataFrame, cohorts: Sequence[str]) -> pd.DataFrame:
    rows = []
    for leave in cohorts:
        keep = [c for c in cohorts if c != leave]
        if len(keep) < 2:
            continue
        rank_cols = [f"rank_{c}" for c in keep]
        tmp = rank_table[["gene"] + rank_cols].copy()
        zcols = []
        for c in keep:
            zc = f"z_{c}"
            tmp[zc] = z_from_rank(tmp[f"rank_{c}"])
            zcols.append(zc)
        tmp["loo_mean_z"] = tmp[zcols].mean(axis=1)
        tmp["left_out_cohort"] = leave
        rows.append(tmp[["gene", "left_out_cohort", "loo_mean_z"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["gene", "left_out_cohort", "loo_mean_z"])


def log_packages(outdir: Path, script_name: str) -> None:
    pkgs = ["numpy", "pandas", "scipy", "sklearn", "lifelines", "sksurv", "matplotlib", "shap"]
    lines = [
        f"Python: {sys.version}",
        f"Platform: {platform.platform()}",
        f"Script: {script_name}",
        f"Random seed: {SEED}",
    ]
    for pkg in pkgs:
        try:
            mod = importlib.import_module(pkg)
            lines.append(f"{pkg}: {getattr(mod, '__version__', 'unknown')}")
        except Exception:
            lines.append(f"{pkg}: not installed")
    outdir.mkdir(exist_ok=True, parents=True)
    (outdir / "package_versions.txt").write_text("\n".join(lines), encoding="utf-8")


def cohort_rank_summary(df: pd.DataFrame, rank_col: str, cohorts: Sequence[str], fill_value: Optional[float] = None) -> pd.DataFrame:
    all_genes = sorted(set(df["gene"]))
    out = pd.DataFrame({"gene": all_genes})
    for cohort in cohorts:
        sub = df.loc[df["cohort"] == cohort, ["gene", rank_col]].copy()
        sub = sub.rename(columns={rank_col: f"rank_{cohort}"})
        out = out.merge(sub, on="gene", how="left")
    if fill_value is None:
        fill_value = len(all_genes) + 1
    for cohort in cohorts:
        out[f"rank_{cohort}"] = out[f"rank_{cohort}"].fillna(fill_value)
    out["category"] = out["gene"].map(gene_cat)
    return out


def meta_from_rank_table(rank_table: pd.DataFrame, cohorts: Sequence[str]) -> pd.DataFrame:
    rows = []
    for _, row in rank_table.iterrows():
        zs = []
        ses = []
        k = 0
        for cohort in cohorts:
            r = row.get(f"rank_{cohort}", np.nan)
            if not np.isfinite(r):
                continue
            cohort_ranks = rank_table[f"rank_{cohort}"]
            zc = z_from_rank(cohort_ranks)
            z = zc.loc[rank_table["gene"] == row["gene"]].iloc[0]
            if np.isfinite(z):
                zs.append(float(z))
                ses.append(1.0)
                k += 1
        meta = random_effects_meta(zs, ses)
        meta["gene"] = row["gene"]
        meta["category"] = row["category"]
        rows.append(meta)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p_adj"] = bh_adjust(out["p"].values)
    return out.sort_values(["p_adj", "random_beta"], ascending=[True, False]).reset_index(drop=True)


def meta_from_cox(cox_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, sub in cox_all.groupby("gene"):
        meta = random_effects_meta(sub["coef"].values, sub["coef_se"].values)
        meta["gene"] = gene
        meta["category"] = gene_cat(gene)
        meta["k_with_cox"] = int(sub.shape[0])
        rows.append(meta)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p_adj"] = bh_adjust(out["p"].values)
    return out.sort_values(["p_adj", "random_beta"], ascending=[True, False]).reset_index(drop=True)


def analyze_one_cohort(cohort: str, manifest_paths: Sequence[str], genes: Sequence[str], args: argparse.Namespace,
                       outdir: Path, notes: List[Dict[str, object]], gene_meta: Optional[pd.DataFrame] = None) -> Optional[CohortResult]:
    expr_path = resolve_from_manifest(manifest_paths, [rf"/{cohort}/X_rna_symbol\.csv$", rf"/{cohort}/X_rna_.*\.csv$"])
    surv_path = resolve_from_manifest(manifest_paths, [rf"/{cohort}/survival\.csv$"])
    risk_path = resolve_risk_path(manifest_paths, cohort)
    if expr_path is None or surv_path is None:
        notes.append({"cohort": cohort, "note": "Skipped: expression or survival file not resolved from manifest."})
        return None

    msg(f"[INFO] {cohort}: expr={expr_path}")
    msg(f"[INFO] {cohort}: surv={surv_path}")
    msg(f"[INFO] {cohort}: risk={risk_path if risk_path else 'not found; RSF fallback'}")

    expr = load_expr(expr_path)
    surv = load_surv(surv_path)
    common = sorted(set(expr.index) & set(surv.index))
    if len(common) < args.min_overlap_samples:
        notes.append({"cohort": cohort, "note": f"Skipped: too few overlapping samples ({len(common)})."})
        return None
    expr = expr.loc[common]
    surv = surv.loc[common]
    avail = [g for g in genes if g in expr.columns]
    if len(avail) < args.min_overlap_genes:
        notes.append({"cohort": cohort, "note": f"Skipped: too few overlapping genes ({len(avail)})."})
        return None

    X = expr[avail].copy()
    scaler = StandardScaler()
    Xs = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)
    T = surv["T"].astype(float).values
    E = surv["E"].astype(int).values

    final_score = None
    final_source = None
    risk_col = None
    if risk_path:
        risk_s, risk_col = load_risk(risk_path)
        if risk_s is not None:
            common2 = sorted(set(Xs.index) & set(risk_s.index))
            if len(common2) >= args.min_overlap_samples:
                Xs = Xs.loc[common2]
                surv = surv.loc[common2]
                T = surv["T"].astype(float).values
                E = surv["E"].astype(int).values
                final_score = risk_s.loc[common2].astype(float).values
                final_source = f"manifest/user risk file: {Path(risk_path).name} | column: {risk_col}"
            else:
                notes.append({"cohort": cohort, "note": f"Risk file overlap too small ({len(common2)}); RSF fallback used."})
    if final_score is None:
        y_surv = Surv.from_arrays(event=E.astype(bool), time=T)
        rsf_tmp = RandomSurvivalForest(
            n_estimators=args.n_estimators,
            min_samples_leaf=5,
            max_features="sqrt",
            n_jobs=-1,
            random_state=SEED,
        )
        rsf_tmp.fit(Xs.values.astype(np.float32), y_surv)
        final_score = rsf_tmp.predict(Xs.values.astype(np.float32))
        final_source = "RSF-generated score (fallback; no risk file found)"

    if len(Xs) < args.min_overlap_samples or np.nansum(E) < args.min_events:
        notes.append({"cohort": cohort, "note": f"Skipped: insufficient samples/events after harmonization (n={len(Xs)}, events={int(np.nansum(E))})."})
        return None

    tau = np.percentile(T[E == 1], 75) if (E == 1).sum() > 5 else np.nanmax(T) * 0.9
    c_low, c_high, _ = bootstrap_metric(T, E, final_score, n_boot=args.n_boot, tau=tau)
    ipcw_ci = ipcw_concordance(T, E, final_score, tau=tau)

    # Secondary sensitivity only: survival-refit permutation importance.
    # This is saved separately and is NOT merged into the primary SHAP/risk-surrogate consensus,
    # because it explains a cohort-refit RSF survival model rather than the final risk-score function.
    y_surv = Surv.from_arrays(event=E.astype(bool), time=T)
    rsf = RandomSurvivalForest(
        n_estimators=args.n_estimators,
        min_samples_leaf=5,
        max_features="sqrt",
        n_jobs=-1,
        random_state=SEED,
    )
    rsf.fit(Xs.values.astype(np.float32), y_surv)
    survival_perm_df, _ = permutation_importance_ipcw(Xs, T, E, rsf, n_perm=args.n_perm, tau=tau, seed=SEED)
    survival_perm_df["cohort"] = cohort
    survival_perm_df.to_csv(outdir / f"{cohort}_survival_refit_permutation_cindex_importance_q1.csv", index=False)
    plot_permutation_bar(survival_perm_df, outdir, cohort, ipcw_ci, c_low, c_high, top_n=args.top_n)

    shap_df = pd.DataFrame(columns=["gene", "mean_abs_shap", "shap_ci_low", "shap_ci_high", "rank_shap", "category", "cohort"])
    shap_direction_df = pd.DataFrame()
    surrogate_r2 = np.nan
    surrogate_rho = np.nan
    surrogate_r2_insample = np.nan
    surrogate_rho_insample = np.nan
    surrogate_r2_cv = np.nan
    surrogate_rho_cv = np.nan
    perm_df = pd.DataFrame({"gene": Xs.columns.tolist(), "rank_perm": np.arange(1, len(Xs.columns) + 1), "category": [gene_cat(g) for g in Xs.columns]})
    try:
        import shap
        et, fid, pred_cv = fit_risk_surrogate_with_cv(
            Xs, final_score, n_estimators=args.surrogate_estimators, seed=SEED, n_splits=5
        )
        surrogate_r2_insample = fid["r2_insample"]
        surrogate_rho_insample = fid["rho_insample"]
        surrogate_r2_cv = fid["r2_cv"]
        surrogate_rho_cv = fid["rho_cv"]
        surrogate_r2 = surrogate_r2_cv
        surrogate_rho = surrogate_rho_cv

        # Primary permutation importance: aligned to SHAP target, i.e. risk-score surrogate.
        perm_df, perm_base = permutation_importance_risk_surrogate(
            Xs, final_score, et, n_perm=args.n_perm, seed=SEED
        )
        perm_df["cohort"] = cohort
        perm_df.to_csv(outdir / f"{cohort}_risk_surrogate_permutation_importance_q1.csv", index=False)
        plot_risk_surrogate_permutation_bar(perm_df, outdir, cohort, surrogate_r2_cv, surrogate_rho_cv, top_n=args.top_n)

        explainer = shap.TreeExplainer(et)
        shap_values = explainer.shap_values(Xs.values.astype(np.float32))
        if getattr(shap_values, "ndim", 2) == 3:
            shap_values = shap_values[..., 0]
        abs_sv = np.abs(shap_values)
        means = np.nanmean(abs_sv, axis=0)
        ci_low = np.nanpercentile(abs_sv, 2.5, axis=0)
        ci_high = np.nanpercentile(abs_sv, 97.5, axis=0)
        shap_df = pd.DataFrame({
            "gene": Xs.columns.tolist(),
            "mean_abs_shap": means,
            "shap_ci_low": ci_low,
            "shap_ci_high": ci_high,
        }).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
        shap_df["rank_shap"] = np.arange(1, len(shap_df) + 1)
        shap_df["category"] = shap_df["gene"].map(gene_cat)
        shap_df["cohort"] = cohort
        shap_df.to_csv(outdir / f"{cohort}_surrogate_SHAP_mean_abs_importance_q1.csv", index=False)
        plot_shap_beeswarm(shap_values, Xs, shap_df["gene"].head(args.top_n).tolist(), outdir, cohort)
        shap_direction_df = compute_shap_direction_table(shap_values, Xs, cohort, gene_meta if gene_meta is not None else pd.DataFrame(), args.support_rho_threshold)
        shap_direction_df.to_csv(outdir / f"{cohort}_directional_SHAP_support_for_ML_genes_q1.csv", index=False)
    except Exception as e:
        notes.append({"cohort": cohort, "note": f"SHAP/risk-surrogate importance failed: {e}"})

    comb = pd.DataFrame({"gene": Xs.columns.tolist()})
    comb = comb.merge(perm_df[["gene", "rank_perm"]], on="gene", how="left")
    if not shap_df.empty:
        comb = comb.merge(shap_df[["gene", "rank_shap"]], on="gene", how="left")
    else:
        comb["rank_shap"] = np.nan
    fill_rank = len(comb) + 1
    comb["rank_perm"] = comb["rank_perm"].fillna(fill_rank)
    comb["rank_shap"] = comb["rank_shap"].fillna(fill_rank)
    comb["combined_mean_rank"] = (
        args.consensus_rank_weight_perm * comb["rank_perm"] +
        args.consensus_rank_weight_shap * comb["rank_shap"]
    )
    comb["category"] = comb["gene"].map(gene_cat)
    if gene_meta is not None and not gene_meta.empty:
        comb = comb.merge(gene_meta[["gene", "ml_direction", "ml_z", "ml_bag_frac", "ml_source"]], on="gene", how="left")
    comb["cohort"] = cohort
    comb = comb.sort_values("combined_mean_rank").reset_index(drop=True)
    comb.to_csv(outdir / f"{cohort}_combined_gene_ranking_q1.csv", index=False)

    cox_df = multivariate_cox_pergene(Xs, surv, comb["gene"].tolist(), top_n=max(args.cox_top_n, args.top_n))
    if not cox_df.empty:
        cox_df.insert(0, "cohort", cohort)
        cox_df.to_csv(outdir / f"{cohort}_multivariate_cox_per_gene_q1.csv", index=False)
        plot_cox_forest(cox_df, outdir, cohort, top_n=min(args.top_n, len(cox_df)))

    plot_calibration(T, E, final_score, outdir, cohort, years=(1, 3, 5))
    plot_dca(T, E, final_score, outdir, cohort, t_yr=3)
    plot_dca(T, E, final_score, outdir, cohort, t_yr=5)

    notes.append({"cohort": cohort, "note": f"Completed successfully. Final score source: {final_source}"})
    return CohortResult(
        cohort=cohort,
        n_samples=len(Xs),
        n_events=int(np.nansum(E)),
        genes_used=list(Xs.columns),
        final_score=np.asarray(final_score, dtype=float),
        final_source=final_source,
        risk_column=risk_col,
        Xs=Xs,
        surv=surv,
        tau=float(tau),
        perm_df=perm_df,
        shap_df=shap_df,
        combined_df=comb,
        cox_df=cox_df,
        ipcw_cindex=float(ipcw_ci),
        cindex_ci_low=float(c_low) if np.isfinite(c_low) else np.nan,
        cindex_ci_high=float(c_high) if np.isfinite(c_high) else np.nan,
        surrogate_r2=float(surrogate_r2) if np.isfinite(surrogate_r2) else np.nan,
        surrogate_rho=float(surrogate_rho) if np.isfinite(surrogate_rho) else np.nan,
        surrogate_r2_insample=float(surrogate_r2_insample) if np.isfinite(surrogate_r2_insample) else np.nan,
        surrogate_rho_insample=float(surrogate_rho_insample) if np.isfinite(surrogate_rho_insample) else np.nan,
        surrogate_r2_cv=float(surrogate_r2_cv) if np.isfinite(surrogate_r2_cv) else np.nan,
        surrogate_rho_cv=float(surrogate_rho_cv) if np.isfinite(surrogate_rho_cv) else np.nan,
        shap_direction_df=shap_direction_df,
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)
    script_name = Path(__file__).name if "__file__" in globals() else "05_shap_interpretability_multicohort_external_q1_v5_manuscript_ready.py"
    log_packages(outdir, script_name)

    manifest_paths = read_manifest(args.manifest)
    gene_meta = read_gene_metadata(args.genes)
    genes = gene_meta["gene"].dropna().astype(str).drop_duplicates().tolist()
    if len(genes) == 0:
        raise RuntimeError("Gene list is empty after parsing.")
    gene_meta.to_csv(outdir / "input_ml_gene_metadata_q1.csv", index=False)

    notes: List[Dict[str, object]] = []
    results: List[CohortResult] = []

    msg("[0/9] Running manuscript-ready cohort-wise interpretability analysis...")
    for cohort in args.cohorts:
        try:
            res = analyze_one_cohort(cohort, manifest_paths, genes, args, outdir, notes, gene_meta)
            if res is not None:
                results.append(res)
        except Exception as e:
            notes.append({"cohort": cohort, "note": f"Unhandled failure: {e}", "traceback": traceback.format_exc()})

    if not results:
        pd.DataFrame(notes).to_csv(outdir / "analysis_notes_for_manuscript_q1.csv", index=False)
        raise RuntimeError("No cohort completed successfully.")

    cohorts_run = [r.cohort for r in results]
    msg(f"[1/9] Completed cohorts: {', '.join(cohorts_run)}")

    perm_all = pd.concat([r.perm_df for r in results], ignore_index=True)
    perm_all.to_csv(outdir / "risk_surrogate_permutation_importance_all_cohorts_q1.csv", index=False)
    shap_nonempty = [r.shap_df for r in results if not r.shap_df.empty]
    if shap_nonempty:
        shap_all = pd.concat(shap_nonempty, ignore_index=True)
        shap_all.to_csv(outdir / "surrogate_SHAP_all_cohorts_q1.csv", index=False)
    else:
        shap_all = pd.DataFrame()

    combined_all = pd.concat([r.combined_df for r in results], ignore_index=True)
    combined_all.to_csv(outdir / "combined_gene_rankings_all_cohorts_q1.csv", index=False)

    direction_nonempty = [r.shap_direction_df for r in results if hasattr(r, "shap_direction_df") and r.shap_direction_df is not None and not r.shap_direction_df.empty]
    if direction_nonempty:
        direction_all = pd.concat(direction_nonempty, ignore_index=True)
        direction_all.to_csv(outdir / "directional_SHAP_support_all_cohorts_q1.csv", index=False)
        direction_summary = aggregate_directional_support(direction_all, gene_meta)
        direction_summary.to_csv(outdir / "directional_SHAP_support_summary_q1.csv", index=False)
    else:
        direction_all = pd.DataFrame()
        direction_summary = pd.DataFrame()

    fidelity_df = pd.DataFrame([{
        "cohort": r.cohort,
        "n_samples": r.n_samples,
        "n_events": r.n_events,
        "n_genes": len(r.genes_used),
        "score_source": r.final_source,
        "risk_column": r.risk_column,
        "ipcw_cindex": r.ipcw_cindex,
        "cindex_ci_lower": r.cindex_ci_low,
        "cindex_ci_upper": r.cindex_ci_high,
        "surrogate_r2_cv": r.surrogate_r2_cv,
        "surrogate_spearman_rho_cv": r.surrogate_rho_cv,
        "surrogate_r2_insample": r.surrogate_r2_insample,
        "surrogate_spearman_rho_insample": r.surrogate_rho_insample,
        "manuscript_use": "Report CV fidelity; in-sample fidelity is diagnostic only.",
    } for r in results])
    fidelity_df.to_csv(outdir / "surrogate_fidelity_summary_q1.csv", index=False)

    msg("[2/9] Cross-cohort stability...")
    rank_table = cohort_rank_summary(combined_all, "combined_mean_rank", cohorts_run)
    fill_value = len(rank_table) + 1
    for c in cohorts_run:
        rank_table[f"rank_{c}"] = rank_table[f"rank_{c}"].fillna(fill_value)
    rank_matrix = np.vstack([rank_table[f"rank_{c}"].values for c in cohorts_run])
    kendall_W, kendall_p = kendalls_w(rank_matrix)
    rank_table["combined_mean_rank"] = rank_table[[f"rank_{c}" for c in cohorts_run]].mean(axis=1)
    rank_table = rank_table.sort_values("combined_mean_rank").reset_index(drop=True)
    rank_table["kendall_W"] = kendall_W
    rank_table["kendall_p"] = kendall_p
    rank_table.to_csv(outdir / "cross_cohort_stability_summary_q1.csv", index=False)
    plot_rank_heatmap(rank_table, cohorts_run, outdir, "cross_cohort_stability_heatmap_q1",
                      f"Cross-cohort feature stability (Kendall's W={kendall_W:.3f}, p={kendall_p:.3g})",
                      top_n=args.top_n)
    plot_consensus_bar(rank_table.rename(columns={"combined_mean_rank": "combined_mean_rank"}), outdir, top_n=args.top_n)

    msg("[3/9] Meta-analysis of feature importance...")
    meta_importance = meta_from_rank_table(rank_table, cohorts_run)
    if not meta_importance.empty:
        meta_importance.to_csv(outdir / "meta_feature_importance_random_effects_q1.csv", index=False)
        plot_meta_importance(meta_importance, outdir, top_n=args.top_n)

    msg("[4/9] Leave-one-cohort-out stability...")
    loo_df = leave_one_out_meta(rank_table, cohorts_run)
    if not loo_df.empty:
        loo_summary = loo_df.groupby("gene", as_index=False).agg(
            loo_mean_z=("loo_mean_z", "mean"),
            loo_sd_z=("loo_mean_z", "std"),
            n_leaveouts=("loo_mean_z", "size")
        ).sort_values("loo_mean_z", ascending=False)
        loo_summary["category"] = loo_summary["gene"].map(gene_cat)
        loo_summary.to_csv(outdir / "leave_one_cohort_out_stability_q1.csv", index=False)

    msg("[5/9] Meta-analysis of adjusted Cox effects...")
    cox_nonempty = [r.cox_df for r in results if not r.cox_df.empty]
    if cox_nonempty:
        cox_all = pd.concat(cox_nonempty, ignore_index=True)
        cox_all.to_csv(outdir / "multivariable_cox_all_cohorts_q1.csv", index=False)
        meta_cox = meta_from_cox(cox_all)
        if not meta_cox.empty:
            meta_cox.to_csv(outdir / "meta_adjusted_cox_random_effects_q1.csv", index=False)
            plot_meta_cox(meta_cox, outdir, top_n=args.top_n)
    else:
        cox_all = pd.DataFrame()
        meta_cox = pd.DataFrame()

    msg("[6/9] Manuscript-facing summary tables...")
    max_rank_for_score = float(np.nanmax(rank_table["combined_mean_rank"])) if "combined_mean_rank" in rank_table.columns else float(len(rank_table))
    rank_table["consensus_importance_score"] = max_rank_for_score - rank_table["combined_mean_rank"] + 1.0
    manuscript_top = rank_table.copy()
    if gene_meta is not None and not gene_meta.empty:
        manuscript_top = manuscript_top.merge(gene_meta[["gene", "ml_direction", "ml_z", "ml_bag_frac", "ml_source"]], on="gene", how="left")
    if not direction_summary.empty:
        manuscript_top = manuscript_top.merge(direction_summary.drop(columns=[c for c in ["category"] if c in direction_summary.columns], errors="ignore"), on="gene", how="left")
    if not meta_importance.empty:
        manuscript_top = manuscript_top.merge(meta_importance[["gene", "random_beta", "random_se", "I2", "p_adj"]].rename(columns={
            "random_beta": "importance_random_beta",
            "random_se": "importance_random_se",
            "p_adj": "importance_p_adj",
        }), on="gene", how="left")
    if not meta_cox.empty:
        manuscript_top = manuscript_top.merge(meta_cox[["gene", "random_beta", "random_se", "p_adj"]].rename(columns={
            "random_beta": "cox_random_beta",
            "random_se": "cox_random_se",
            "p_adj": "cox_p_adj",
        }), on="gene", how="left")
    manuscript_top["category"] = manuscript_top["gene"].map(gene_cat)
    # Robust defaults: DataFrame.get(..., scalar) returns a scalar, which has no .fillna().
    if "ml_direction" in manuscript_top.columns:
        manuscript_top["ml_direction"] = manuscript_top["ml_direction"].fillna("Unknown")
    else:
        manuscript_top["ml_direction"] = pd.Series("Unknown", index=manuscript_top.index)

    ml_z_series = manuscript_top["ml_z"] if "ml_z" in manuscript_top.columns else pd.Series(np.nan, index=manuscript_top.index)
    manuscript_top["expected_risk_sign"] = [ml_direction_sign(d, z) for d, z in zip(manuscript_top["ml_direction"], ml_z_series)]

    concordant_fraction = manuscript_top["concordant_fraction"].fillna(0) if "concordant_fraction" in manuscript_top.columns else pd.Series(0.0, index=manuscript_top.index)
    discordant_fraction = manuscript_top["discordant_fraction"].fillna(0) if "discordant_fraction" in manuscript_top.columns else pd.Series(0.0, index=manuscript_top.index)
    manuscript_top["shap_concordance_label"] = concordant_fraction
    manuscript_top["shap_support_tier"] = np.where(
        concordant_fraction >= 0.67, "Strong SHAP-concordant",
        np.where(concordant_fraction >= 0.34, "Moderate SHAP-concordant",
                 np.where(discordant_fraction >= 0.67, "Discordant", "Weak/heterogeneous"))
    )
    order_map = {"Strong SHAP-concordant": 1, "Moderate SHAP-concordant": 2, "Weak/heterogeneous": 3, "Discordant": 4}
    manuscript_top["support_tier_order"] = manuscript_top["shap_support_tier"].map(order_map).fillna(5)
    manuscript_top = manuscript_top.sort_values(["support_tier_order", "combined_mean_rank"], ascending=[True, True]).reset_index(drop=True)
    manuscript_top.to_csv(outdir / "manuscript_top_genes_summary_q1.csv", index=False)
    manuscript_top.head(args.manuscript_top_n).to_csv(outdir / "manuscript_top_genes_for_main_text_q1.csv", index=False)
    plot_consensus_bar(manuscript_top.sort_values("combined_mean_rank"), outdir, top_n=args.manuscript_top_n)
    plot_manuscript_ml_gene_support(manuscript_top, outdir, top_n=args.manuscript_top_n)
    plot_directional_shap_support(manuscript_top, outdir, top_n=args.manuscript_top_n)
    write_manuscript_interpretability_note(manuscript_top, fidelity_df, kendall_W, kendall_p, outdir, top_n=min(10, args.manuscript_top_n))

    category_rows = []
    for r in results:
        top = r.combined_df.head(args.top_n)
        cnt = top["category"].value_counts()
        for cat in ["Driver", "Suppressor", "Marker/Other"]:
            category_rows.append({"cohort": r.cohort, "category": cat, "top_n_count": int(cnt.get(cat, 0))})
    pd.DataFrame(category_rows).to_csv(outdir / "feature_category_summary_q1.csv", index=False)

    pd.DataFrame(notes).to_csv(outdir / "analysis_notes_for_manuscript_q1.csv", index=False)

    msg("[7/9] Output manifest and run summary...")
    files = []
    for pth in sorted(outdir.rglob("*")):
        if pth.is_file():
            files.append({"relative_path": str(pth.relative_to(outdir)), "bytes": pth.stat().st_size})
    pd.DataFrame(files).to_csv(outdir / "output_manifest_q1.csv", index=False)

    run_summary = {
        "cohorts_requested": args.cohorts,
        "cohorts_completed": cohorts_run,
        "n_signature_genes_input": len(genes),
        "kendall_W": None if pd.isna(kendall_W) else float(kendall_W),
        "kendall_p": None if pd.isna(kendall_p) else float(kendall_p),
        "seed": SEED,
        "n_perm": args.n_perm,
        "n_boot": args.n_boot,
        "n_estimators": args.n_estimators,
        "surrogate_estimators": args.surrogate_estimators,
        "top_n": args.top_n,
        "script": script_name,
        "interpretation_note": "SHAP values and primary permutation ranks are derived from / aligned to an ExtraTrees surrogate trained on the final cohort-specific risk score. Directional SHAP support is compared with ML-derived risk/protective gene direction; these outputs support model interpretation but do not imply causal gene effects.",
        "q1_note": "This version reports cross-validated surrogate fidelity for manuscript use and retains in-sample fidelity as a diagnostic only. It also keeps RSF survival-refit permutation as a separate sensitivity output rather than merging it with SHAP.",
    }
    (outdir / "run_summary_q1.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")

    msg("[8/9] Plain-text manuscript note...")
    txt = [
        "Q1-oriented external-cohort interpretability analysis summary",
        f"Completed cohorts: {', '.join(cohorts_run)}",
        f"Input signature genes: {len(genes)}",
        f"Cross-cohort concordance (Kendall's W): {kendall_W:.4f}" if np.isfinite(kendall_W) else "Cross-cohort concordance: NA",
        f"Cross-cohort concordance p-value: {kendall_p:.4g}" if np.isfinite(kendall_p) else "Cross-cohort concordance p-value: NA",
        "",
        "Reviewer-facing upgrades:",
        "- IPCW C-index bootstrap intervals were added.",
        "- Primary permutation importance is target-aligned with SHAP and explains the risk-score surrogate.",
        "- Directional SHAP support tables compare feature-SHAP direction with ML-derived risk/protective direction.",
        "- Manuscript-ready support figures and text snippets are exported for the main paper.",
        "- RSF survival-refit permutation is retained separately as a sensitivity output.",
        "- SHAP remains explicitly framed as surrogate-model attribution, with CV fidelity reported for manuscript use.",
        "- Cross-cohort consensus, random-effects meta-analysis, heterogeneity, leave-one-cohort-out stability, score-recalibrated calibration, and score-recalibrated DCA were added.",
        "- Adjusted per-gene Cox models and cross-cohort meta-Cox summaries were added where estimable.",
        "- Output provenance and package versions were recorded for reproducibility.",
    ]
    (outdir / "manuscript_notes_q1.txt").write_text("\n".join(txt), encoding="utf-8")

    msg("[9/9] DONE")
    print(f"  Output directory: {outdir}")
    print(f"  Cohorts completed: {', '.join(cohorts_run)}")
    if np.isfinite(kendall_W):
        print(f"  Cross-cohort stability: Kendall's W={kendall_W:.4f}, p={kendall_p:.4g}")


if __name__ == "__main__":
    main()
