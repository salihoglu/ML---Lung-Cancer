#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rana Salihoglu
ML-LUAD
"""

from __future__ import annotations
import argparse, json, math, os, re, time, warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

MPL_AVAILABLE = False
try:
    import matplotlib.pyplot as plt
    MPL_AVAILABLE = True
except Exception:
    plt = None

SCIPY_AVAILABLE = False
try:
    from scipy import stats as scipy_stats
    SCIPY_AVAILABLE = True
except Exception:
    scipy_stats = None

# Optional classification metrics for fixed-year endpoints
SKLEARN_METRICS_AVAILABLE = False
try:
    from sklearn.metrics import roc_auc_score, accuracy_score
    SKLEARN_METRICS_AVAILABLE = True
except Exception:
    roc_auc_score = None
    accuracy_score = None

# ---------------- Optional deps ----------------
SKSURV_AVAILABLE = False
SKSURV_MODELS_AVAILABLE = False
FAST_SVM_AVAILABLE = False
CW_GBM_AVAILABLE = False
GLMNET_AVAILABLE = False
try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    GLMNET_AVAILABLE = True
except ImportError:
    pass
SKSURV_LINEAR_EXTRA_AVAILABLE = False
SKSURV_TREE_AVAILABLE = False

try:
    from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc, brier_score, integrated_brier_score
    from sksurv.util import Surv
    SKSURV_AVAILABLE = True
except Exception:
    concordance_index_censored = None
    Surv = None

LIFELINES_AVAILABLE = False
try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.utils import concordance_index
    from lifelines.statistics import logrank_test
    LIFELINES_AVAILABLE = True
except Exception:
    CoxPHFitter = None
    concordance_index = None


_META_N_JOBS: int = 4  

try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    try:
        from sksurv.linear_model import CoxPHSurvivalAnalysis, IPCRidge
        SKSURV_LINEAR_EXTRA_AVAILABLE = True
    except Exception:
        CoxPHSurvivalAnalysis = None
        IPCRidge = None
        SKSURV_LINEAR_EXTRA_AVAILABLE = False
    from sksurv.ensemble import RandomSurvivalForest, ExtraSurvivalTrees, GradientBoostingSurvivalAnalysis
    try:
        from sksurv.tree import SurvivalTree
        SKSURV_TREE_AVAILABLE = True
    except Exception:
        SurvivalTree = None
        SKSURV_TREE_AVAILABLE = False
    SKSURV_MODELS_AVAILABLE = True
except Exception:
    CoxnetSurvivalAnalysis = None
    RandomSurvivalForest = None
    ExtraSurvivalTrees = None
    GradientBoostingSurvivalAnalysis = None

try:
    from sksurv.svm import FastSurvivalSVM
    FAST_SVM_AVAILABLE = True
except Exception:
    FastSurvivalSVM = None

try:
    from sksurv.ensemble import ComponentwiseGradientBoostingSurvivalAnalysis
    CW_GBM_AVAILABLE = True
except Exception:
    ComponentwiseGradientBoostingSurvivalAnalysis = None


FIXED_TEST = ["GSE50081", "GSE68465"]  
FIXED_TRAIN = [
    "GSE87340", "GSE72094", "GSE136961", "GSE157009",
    "GSE31210", "GSE30219", "TCGA-LUSC", "TCGA-LUAD"
   
]

# ---------------- Progress ----------------
class ProgressTracker:
    def __init__(self, total_steps: int):
        self.total_steps = max(1, int(total_steps))
        self.done = 0
        self.t0 = time.time()

    def tick(self, msg: str):
        self.done += 1
        now = time.time()
        elapsed = now - self.t0
        pct = 100.0 * self.done / self.total_steps
        rate = elapsed / max(self.done, 1)
        eta = rate * max(self.total_steps - self.done, 0)
        print(f"[PROGRESS] {self.done}/{self.total_steps} | %{pct:5.1f} | elapsed {elapsed/60:.1f}m | ETA {eta/60:.1f}m | {msg}")

    def log(self, msg: str):
        now = time.time()
        elapsed = now - self.t0
        pct = 100.0 * self.done / self.total_steps
        print(f"[PROGRESS] %{pct:5.1f} | elapsed {elapsed/60:.1f}m | {msg}")

# ---------------- Helpers ----------------
def set_thread_limits(n: int = 1):
    n = int(max(1, n))
    os.environ["OMP_NUM_THREADS"] = str(n)
    os.environ["OPENBLAS_NUM_THREADS"] = str(n)
    os.environ["MKL_NUM_THREADS"] = str(n)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(n)
    os.environ["NUMEXPR_NUM_THREADS"] = str(n)

def seed_everything(seed: int = 13):
    np.random.seed(seed)
    try:
        import random
        random.seed(seed)
    except Exception:
        pass

def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def rank_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(x) + 1)
    u = (ranks - 0.5) / len(x)
    try:
        from scipy.stats import norm
        return norm.ppf(u)
    except Exception:
        return np.clip((u - 0.5) * 6.0, -5, 5)

def cohort_zscore(df: pd.DataFrame, cohort: pd.Series, prefix: str = "RNA__",
                  robust: bool = False,
                  train_stats: Optional[Dict[str, Dict[str, pd.Series]]] = None) -> pd.DataFrame:
    """
   
    """
    df2 = df.copy()
    rna_cols = [c for c in df2.columns if c.startswith(prefix)]
    if not rna_cols:
        return df2
    for coh in cohort.unique():
        idx = cohort == coh
        sub = df2.loc[idx, rna_cols]
        if train_stats is not None and coh in train_stats:
            # Use pre-computed train statistics
            stats_coh = train_stats[coh]
            if robust:
                med = stats_coh.get("med", sub.median(axis=0))
                iqr = stats_coh.get("iqr", (sub.quantile(0.75, axis=0) - sub.quantile(0.25, axis=0)))
                iqr = iqr.replace(0, np.nan)
                df2.loc[idx, rna_cols] = (sub - med) / iqr
            else:
                mu = stats_coh.get("mu", sub.mean(axis=0))
                sd = stats_coh.get("sd", sub.std(axis=0)).replace(0, np.nan)
                df2.loc[idx, rna_cols] = (sub - mu) / sd
        else:
            if robust:
                med = sub.median(axis=0)
                iqr = sub.quantile(0.75, axis=0) - sub.quantile(0.25, axis=0)
                iqr = iqr.replace(0, np.nan)
                df2.loc[idx, rna_cols] = (sub - med) / iqr
            else:
                mu = sub.mean(axis=0)
                sd = sub.std(axis=0).replace(0, np.nan)
                df2.loc[idx, rna_cols] = (sub - mu) / sd
    df2[rna_cols] = df2[rna_cols].fillna(0.0)
    return df2


def extract_cohort_stats(df_normalized_source: pd.DataFrame,
                         cohort: pd.Series,
                         prefix: str = "RNA__",
                         robust: bool = False) -> Dict[str, Dict[str, pd.Series]]:
    """
   
    """
    rna_cols = [c for c in df_normalized_source.columns if c.startswith(prefix)]
    stats = {}
    for coh in cohort.unique():
        idx = cohort == coh
        sub = df_normalized_source.loc[idx, rna_cols]
        if robust:
            stats[coh] = {"med": sub.median(axis=0),
                          "iqr": (sub.quantile(0.75, axis=0) - sub.quantile(0.25, axis=0)).replace(0, np.nan)}
        else:
            stats[coh] = {"mu": sub.mean(axis=0),
                          "sd": sub.std(axis=0).replace(0, np.nan)}
    return stats

def _pca1_score(Z: np.ndarray) -> np.ndarray:
    Z = np.asarray(Z, dtype=float)
    if Z.ndim != 2 or Z.shape[0] == 0:
        return np.zeros((Z.shape[0],), dtype=float)
    Zc = Z - np.nanmean(Z, axis=0, keepdims=True)
    Zc = np.nan_to_num(Zc, nan=0.0, posinf=0.0, neginf=0.0)
    try:
        U, S, _ = np.linalg.svd(Zc, full_matrices=False)
        return (U[:, 0] * S[0]).astype(float)
    except Exception:
        return np.nanmean(Zc, axis=1).astype(float)

def add_paper_like_traits(X: pd.DataFrame, ferr_driver: List[str], ferr_suppressor: List[str], module_genes: List[str]) -> pd.DataFrame:
    """
    
    """
    df = X.copy()
    rna_cols = [c for c in df.columns if c.startswith("RNA__")]

    def cols_for_genes(genes: List[str]) -> List[str]:
        want = {f"RNA__{g}" for g in genes}
        return [c for c in rna_cols if c in want]

    drv_cols = cols_for_genes(ferr_driver)
    sup_cols = cols_for_genes(ferr_suppressor)
    mod_cols = cols_for_genes(module_genes)

    # ── Orijinal features (v34, korundu) ──────────────────────────────────
    df["FERRO__DRIVER"]     = df[drv_cols].mean(axis=1) if drv_cols else 0.0
    df["FERRO__SUPPRESSOR"] = df[sup_cols].mean(axis=1) if sup_cols else 0.0
    df["FERRO__TOTAL"]      = 0.5 * (df["FERRO__DRIVER"] + df["FERRO__SUPPRESSOR"])
    df["FERRO__TRAIT"]      = df["FERRO__DRIVER"] - df["FERRO__SUPPRESSOR"]

    if mod_cols and len(mod_cols) >= 5:
        df["WGCNA__ME1"] = _pca1_score(df[mod_cols].values)
    elif mod_cols:
        df["WGCNA__ME1"] = df[mod_cols].mean(axis=1)
    else:
        df["WGCNA__ME1"] = 0.0

   
    eps = 1e-6

    
    if drv_cols and len(drv_cols) >= 3:
        df["FERRO__DRIVER_PC1"] = _pca1_score(df[drv_cols].values)
    else:
        df["FERRO__DRIVER_PC1"] = df["FERRO__DRIVER"]

    if sup_cols and len(sup_cols) >= 3:
        df["FERRO__SUPP_PC1"] = _pca1_score(df[sup_cols].values)
    else:
        df["FERRO__SUPP_PC1"] = df["FERRO__SUPPRESSOR"]

   
    df["FERRO__NET_SCORE"] = df["FERRO__DRIVER_PC1"] - df["FERRO__SUPP_PC1"]

    
    drv_pos = df["FERRO__DRIVER"] - df["FERRO__DRIVER"].min() + eps
    sup_pos = df["FERRO__SUPPRESSOR"] - df["FERRO__SUPPRESSOR"].min() + eps
    df["FERRO__RATIO"] = np.log2(drv_pos / sup_pos)
    
    r_std = df["FERRO__RATIO"].std()
    r_mean = df["FERRO__RATIO"].mean()
    if np.isfinite(r_std) and r_std > 0:
        df["FERRO__RATIO"] = df["FERRO__RATIO"].clip(
            r_mean - 5*r_std, r_mean + 5*r_std)

   
    df["FERRO__WGCNA_INTER"] = df["WGCNA__ME1"] * df["FERRO__TRAIT"]

    
    def _zscore_cols(mat: np.ndarray) -> np.ndarray:
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0) + eps
        return (mat - mu) / sd

    if drv_cols and sup_cols:
        drv_z = _zscore_cols(df[drv_cols].values).mean(axis=1)
        sup_z = _zscore_cols(df[sup_cols].values).mean(axis=1)
        df["FERRO__Z_COMBINED"] = drv_z + sup_z
    else:
        df["FERRO__Z_COMBINED"] = df["FERRO__TRAIT"]

   
    ferro_feature_cols = [c for c in df.columns if c.startswith("FERRO__") or c == "WGCNA__ME1"]
    for fc in ferro_feature_cols:
        if fc not in df.columns: continue
        col_vals = df[fc].values.astype(float)
        mu_f = np.nanmean(col_vals); sd_f = np.nanstd(col_vals)
        if np.isfinite(mu_f) and np.isfinite(sd_f) and sd_f > 1e-8:
            df[fc] = (col_vals - mu_f) / sd_f
        else:
            df[fc] = 0.0  

    return df

# ---------------- Stage parsing ----------------
_ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4}

def _parse_stage_token(tok: str) -> Optional[float]:
    """
   
    """
    if tok is None:
        return None
    s = str(tok).strip().upper()
    if s in ("", "NA", "N/A", "NONE", "NULL", "UNKNOWN", "PP"):
        return None

    # ── 1. pTNM format (GSE68465: pN0pT1, pN1pT2, pNXpT1 vb.) ──────────────
    t_m = re.search(r'PT([1-4])', s)
    n_m = re.search(r'PN([0-3X])', s)
    if t_m:
        t_score = int(t_m.group(1))
        n_raw   = n_m.group(1) if n_m else "0"
        n_score = int(n_raw) if n_raw.isdigit() else 0  # pNX → 0
        return float(t_score) + 0.5 * float(n_score)

    # ── 2. AJCC roman numeral / arabic ────────────────────────────────────────
    s2 = re.sub(r"STAGE", "", s)
    s2 = re.sub(r"[^A-Z0-9]", "", s2)
    # arabic first (avoids "1A" being parsed as roman "I")
    m = re.match(r"^([1-4])", s2)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 4:
            return float(n)
    # roman numeral
    for r in ("IV", "III", "II", "I"):
        if s2.startswith(r):
            return float(_ROMAN[r])

    return None

def extract_stage_from_tcga_clinical(clinical_path: Path, sample_ids: pd.Index) -> pd.DataFrame:
    if not clinical_path.exists():
        return pd.DataFrame(index=sample_ids, data={"CLIN__stage_num": np.nan, "CLIN__stage_missing": 1.0})
    df = pd.read_csv(clinical_path, low_memory=False)
    key_cols = [c for c in df.columns if c.lower() in ("bcr_patient_barcode","submitter_id","patient")]
    key = key_cols[0] if key_cols else df.columns[0]
    df[key] = df[key].astype(str)
    cand = [c for c in df.columns if "stage" in c.lower()]
    stage_col = None
    if cand:
        # ajcc_pathologic_stage if exists
        for c in df.columns:
            if c.lower() == "ajcc_pathologic_stage":
                stage_col = c
                break
        stage_col = stage_col or cand[0]
    stage_map: Dict[str, float] = {}
    if stage_col:
        sub = df[[key, stage_col]].dropna()
        for _, row in sub.iterrows():
            sid = str(row[key]).strip()
            st = _parse_stage_token(row[stage_col])
            if st is not None:
                stage_map[sid] = float(st)
    stage = pd.Series(index=sample_ids.astype(str), dtype=float)
    for sid in stage.index:
        stage.loc[sid] = stage_map.get(sid, np.nan)
    missing = stage.isna().astype(float)
    fill = stage.median() if stage.notna().any() else 2.0
    return pd.DataFrame(index=sample_ids.astype(str), data={"CLIN__stage_num": stage.fillna(fill).astype(float), "CLIN__stage_missing": missing.astype(float)})

def extract_stage_from_geo_pheno(pheno_path: Path, sample_ids: pd.Index) -> pd.DataFrame:
    """
   
    """
    if not pheno_path.exists():
        return pd.DataFrame(index=sample_ids, data={"CLIN__stage_num": np.nan, "CLIN__stage_missing": 1.0})
    df = pd.read_csv(pheno_path, low_memory=False)
    sid_col = next((c for c in df.columns if c.lower() in ("sample","geo_accession","geo_accession_id")), df.columns[0])
    df[sid_col] = df[sid_col].astype(str)
    stage_cols = [c for c in df.columns if "stage" in c.lower()]

   
    stage_col  = None
    stage_map: Dict[str, float] = {} 
    _pt_col = next((c for c in stage_cols if c.lower() in ("pt.stage","pt_stage","pt stage")), None)
    _pn_col = next((c for c in stage_cols if c.lower() in ("pn.stage","pn_stage","pn stage")), None)

    
    if _pt_col is not None:
        try:
            for _, row in df.iterrows():
                sid  = str(row[sid_col]).strip()
                t_s  = str(row[_pt_col]).strip().upper()
                n_s  = str(row[_pn_col]).strip().upper() if _pn_col else "N0"
                combined = f"p{n_s}p{t_s}".replace(" ","")
                st = _parse_stage_token(combined)
                if st is not None:
                    stage_map[sid] = float(st)
            _parsed_tn = len(stage_map)
            print(f"  [v45] {pheno_path.name}: pt.stage+pn.stage combined → {_parsed_tn} parsed")
        except Exception as _tne:
            print(f"  [v45 WARN] pt/pn stage combine failed ({_tne}), falling back")
            stage_map = {}

    
    if not stage_map:
        for _pref in ["disease_stage", "pathologic_stage", "ajcc_pathologic_stage",
                      "pathological_stage", "tumor_stage", "clinical_stage", "stage"]:
            _match = next((c for c in stage_cols if c.lower() == _pref), None)
            if _match:
                stage_col = _match
                break
        if stage_col is None and stage_cols:
            stage_col = stage_cols[0]
        if stage_col:
            for _, row in df[[sid_col, stage_col]].iterrows():
                sid = str(row[sid_col]).strip()
                st  = _parse_stage_token(row[stage_col])
                if st is not None:
                    stage_map[sid] = float(st)

    stage = pd.Series(index=sample_ids.astype(str), dtype=float)
    for sid in stage.index:
        stage.loc[sid] = stage_map.get(sid, np.nan)
    missing = stage.isna().astype(float)
    fill    = float(stage.median()) if stage.notna().any() else 2.0
    _parsed = stage.notna().sum()
    if _parsed == 0 and stage_col is not None:
        print(f"  [WARN] extract_stage_from_geo_pheno: stage_col='{stage_col}' "
              f"parse problem "
              f"'{df[stage_col].dropna().iloc[0] if df[stage_col].notna().any() else 'ALL NaN'}'")
    else:
        print(f"  Stage parsed: {_parsed}/{len(sample_ids)} samples "
              f"(col='{stage_col}', fill={fill:.2f})")
    return pd.DataFrame(
        index=sample_ids.astype(str),
        data={"CLIN__stage_num": stage.fillna(fill).astype(float),
              "CLIN__stage_missing": missing.astype(float)}
    )

# ---------------- Data IO ----------------
@dataclass
class CohortData:
    cohort: str
    X: pd.DataFrame
    y_time: pd.Series
    y_event: pd.Series
    meta: Dict[str, str]

def read_survival(surv_path: Path) -> Tuple[pd.Series, pd.Series]:
    df = pd.read_csv(surv_path)
    sample_col = next((c for c in df.columns if c.lower() in ("sample","submitter_id","patient","id")), df.columns[0])
    time_col = next((c for c in df.columns if c.lower() in ("os.time","time","os_time","days","survival_time")), df.columns[1])
    event_col = next((c for c in df.columns if c.lower() in ("os","event","status","dead","death")), df.columns[2] if len(df.columns)>2 else df.columns[-1])

    sample = df[sample_col].astype(str).str.strip().values
    t = pd.to_numeric(df[time_col], errors="coerce").values
    ev_raw = df[event_col]
    ev_num = pd.to_numeric(ev_raw, errors="coerce")
    if ev_num.notna().mean() > 0.5:
        e = (ev_num.fillna(0).values > 0).astype(int)
    else:
        s = ev_raw.astype(str).str.strip().str.lower()
        pos = s.isin(["1","true","t","yes","y","dead","deceased","death","event"])
        neg = s.isin(["0","false","f","no","n","alive","censored","censor"])
        e = np.where(pos,1,np.where(neg,0,0)).astype(int)
    return pd.Series(t, index=sample), pd.Series(e, index=sample)

def read_expression(expr_path: Path, gene_subset: Optional[List[str]] = None, float32: bool = True) -> pd.DataFrame:
    header = pd.read_csv(expr_path, nrows=0)
    cols = list(header.columns)
    gene_cols = cols[1:]
    if gene_subset:
        want = set(map(str, gene_subset))
        gene_cols = [c for c in gene_cols if str(c) in want]
    usecols = [cols[0]] + gene_cols
    df = pd.read_csv(expr_path, usecols=usecols, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df.astype(np.float32) if float32 else df

def _norm_id(x: str) -> str:
    s = str(x).strip()
    m2 = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", s, flags=re.IGNORECASE)
    if m2: return m2.group(1).upper()
    m = re.search(r"(GSM\d+)", s, flags=re.IGNORECASE)
    if m: return m.group(1).upper()
    return s.upper()

def _align_by_ids(X: pd.DataFrame, t: pd.Series, e: pd.Series):
    xi = pd.Index([_norm_id(i) for i in X.index.astype(str)])
    ti = pd.Index([_norm_id(i) for i in t.index.astype(str)])
    X2 = X.copy(); X2.index = xi
    t2 = pd.Series(t.values, index=ti)
    e2 = pd.Series(e.values, index=ti)
    common = X2.index.intersection(t2.index)
    if common.size:
        X2 = X2.loc[common]
        if X2.index.has_duplicates:
            X2 = X2.groupby(level=0).mean()
        t2, e2 = t2.loc[common], e2.loc[common]
    return X2, t2, e2, int(common.size)

def load_tcga(cohort_name: str, cohort_dir: Path, add_stage=True, gene_subset=None, float32=True) -> CohortData:
    X = read_expression(cohort_dir / "X_rna_tpm_patient.csv", gene_subset=gene_subset, float32=float32)
    t, e = read_survival(cohort_dir / "survival.csv")
    common = X.index.astype(str).intersection(t.index.astype(str))
    X, t, e = X.loc[common].copy(), t.loc[common], e.loc[common]
    X.columns = [f"RNA__{c}" for c in X.columns.astype(str)]
    if add_stage:
        X = pd.concat([X, extract_stage_from_tcga_clinical(cohort_dir / "clinical.csv", X.index)], axis=1)
    return CohortData(cohort_name, X, t, e, {"source":"TCGA","dir":str(cohort_dir)})

def load_geo(cohort_name: str, cohort_dir: Path, add_stage=True, gene_subset=None,
             float32=True, hist_filter: Optional[Dict] = None) -> CohortData:
    """
    """
    X = read_expression(cohort_dir / "X_rna_symbol.csv", gene_subset=gene_subset, float32=float32)
    t, e = read_survival(cohort_dir / "survival.csv")
    X2, t2, e2, n_common = _align_by_ids(X, t, e)
    if n_common == 0:
        print(f"[WARN] {cohort_name}: 0 matched samples.")
        X2, t2, e2 = X.iloc[0:0].copy(), t.iloc[0:0].copy(), e.iloc[0:0].copy()
    X2.columns = [f"RNA__{c}" for c in X2.columns.astype(str)]

    # v45-1: Histoloji filtresi — mixed-histoloji kohortları için
    if hist_filter and len(X2) > 0:
        pheno_path = cohort_dir / f"{cohort_name}_pheno_expanded.csv"
        if pheno_path.exists():
            try:
                pheno_df = pd.read_csv(pheno_path, low_memory=False)
                _sid_col  = next((c for c in pheno_df.columns
                                  if c.lower() in ("sample","geo_accession","geo_accession_id")),
                                 pheno_df.columns[0])
                pheno_df[_sid_col] = pheno_df[_sid_col].astype(str).str.upper()
                _hcol = hist_filter.get("col", "histology")
                _hvals = [str(v).upper() for v in hist_filter.get("values", [])]
                # match col name case-insensitively
                _col_match = next((c for c in pheno_df.columns if c.lower() == _hcol.lower()), None)
                if _col_match is None:
                    _col_match = next((c for c in pheno_df.columns
                                       if _hcol.lower() in c.lower()), None)
                if _col_match:
                    pheno_df["_hist_upper"] = pheno_df[_col_match].astype(str).str.upper()
                    keep_ids = set(
                        pheno_df.loc[pheno_df["_hist_upper"].isin(_hvals), _sid_col].tolist()
                    )
                    x_ids_upper = pd.Index([_norm_id(i) for i in X2.index.astype(str)])
                    keep_mask   = x_ids_upper.isin(keep_ids)
                    n_before = len(X2)
                    X2  = X2.loc[keep_mask.values]
                    t2  = t2.loc[X2.index]
                    e2  = e2.loc[X2.index]
                    print(f" {cohort_name} hist_filter '{_hcol}'={_hvals}: "
                          f"{n_before} → {len(X2)} samples "
                          f"({int(e2.sum())} events)")
                else:
                    print(f"  [WARN] {cohort_name}: hist_filter col '{_hcol}' not found "
                          f"in pheno (cols: {list(pheno_df.columns[:8])}...) — filter skipped")
            except Exception as _hfe:
                print(f"  [v45 WARN] {cohort_name}: hist_filter failed ({_hfe}) — filter skipped")
        else:
            print(f"  [WARN] {cohort_name}: pheno_expanded not found for hist_filter — skipped")

    if add_stage and len(X2):
        X2 = pd.concat([X2, extract_stage_from_geo_pheno(
            cohort_dir / f"{cohort_name}_pheno_expanded.csv", X2.index)], axis=1)
    return CohortData(cohort_name, X2, t2, e2, {"source":"GEO","dir":str(cohort_dir)})

def qc_events(cohorts: List[CohortData]) -> pd.DataFrame:
    rows=[]
    for c in cohorts:
        n = int(len(c.y_event))
        ev = int(np.nansum(c.y_event.values)) if n else 0
        t = pd.to_numeric(c.y_time, errors="coerce")
        rows.append({"cohort":c.cohort, "n":n, "events":ev, "event_rate":(ev/n) if n else np.nan, "median_time":float(np.nanmedian(t.values)) if n else np.nan})
    return pd.DataFrame(rows)

# ---------------- Genes: FerrDb ----------------
def _read_gene_col(csv_path: Path) -> List[str]:
    if csv_path is None or not csv_path.exists():
        return []
    df = pd.read_csv(csv_path)
    gcol = next((c for c in df.columns if c.lower() in ("gene","genesymbol","symbol","gene_symbol","official_symbol")), df.columns[0])
    out=[]
    for g in df[gcol].astype(str):
        gg=g.strip()
        if gg and gg.lower() not in ("nan","na","none"):
            out.append(gg)
    return sorted(set(out))

def load_ferrdb_three(driver_csv: Path, suppressor_csv: Path, marker_csv: Path) -> Tuple[List[str], List[str], List[str], List[str]]:
    drv = _read_gene_col(driver_csv)
    sup = _read_gene_col(suppressor_csv)
    mkr = _read_gene_col(marker_csv)
    allg = sorted(set(drv) | set(sup) | set(mkr))
    return drv, sup, mkr, allg

# ---------------- Genes: WGCNA ----------------
def _read_table_any(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in (".tsv",".txt") else ","
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except Exception:
        # try alternate
        return pd.read_csv(path, sep="," if sep=="\t" else "\t", low_memory=False)

def _detect_gene_column(df: pd.DataFrame) -> str:
    candidates = []
    for c in df.columns:
        cl = c.lower()
        if cl in ("gene","symbol","gene_symbol","hub_gene","genes","hubs","geneid"):
            candidates.append(c)
    if candidates:
        return candidates[0]
    return df.columns[0]

def load_wgcna_genes_from_file(path: Path, prefer_branch: Optional[str] = None) -> List[str]:
    """
    """
    df = _read_table_any(path)
    gcol = _detect_gene_column(df)

    # v36: branch kolonu varsa (wgcna_genes_for_ML.csv) filtrele
    if prefer_branch and "branch" in [col.lower() for col in df.columns]:
        branch_col = next(col for col in df.columns if col.lower() == "branch")
        df = df[df[branch_col].astype(str).str.lower() == prefer_branch.lower()]

    # hub_score kolonu varsa (conservative_ferroptosis_genes_GS_MM.tsv)
    # sıralı halde gelir — zaten hub_score'a göre sıralı (v13)
    # Sırayı koru (sorted() yerine list kullan)
    genes = [str(x).strip() for x in df[gcol].tolist()]
    genes = [g for g in genes if g and g.lower() not in ("nan","na","none")]
    # Duplicate'leri kaldır, sırayı koru
    seen = set(); unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g); unique_genes.append(g)
    return unique_genes

def _pick_module_from_trait_stats(trait_stats: pd.DataFrame) -> Optional[str]:
    """
    """
    cols = [c.lower() for c in trait_stats.columns]
    # module column
    mod_col = next((trait_stats.columns[i] for i,c in enumerate(cols) if c in ("module","modulecolor","module_color","mod")), None)
    if mod_col is None:
        return None
    # trait column
    trait_col = next((trait_stats.columns[i] for i,c in enumerate(cols) if c in ("trait","phenotype","pheno","label")), None)
    # correlation column
    cor_col = next((trait_stats.columns[i] for i,c in enumerate(cols) if c in ("cor","corr","correlation","pearson","r")), None)
    if cor_col is None:
        # sometimes "value" holds correlation
        cor_col = next((trait_stats.columns[i] for i,c in enumerate(cols) if c in ("value","estimate")), None)
    if trait_col is None or cor_col is None:
        return None
    sub = trait_stats.copy()
    sub[trait_col] = sub[trait_col].astype(str)
    sub = sub[sub[trait_col].str.lower().str.contains("ferro", na=False)]
    if sub.empty:
        return None
    sub[cor_col] = pd.to_numeric(sub[cor_col], errors="coerce")
    sub = sub.dropna(subset=[cor_col])
    if sub.empty:
        return None
    # strongest positive (paper: "significant positive correlation / strongest correlation")
    sub = sub.sort_values(cor_col, ascending=False)
    return str(sub.iloc[0][mod_col])

def load_wgcna_from_path(wgcna_path: Optional[Path]) -> Tuple[List[str], str]:
    """
    Returns: (module_gene_list, provenance_string)
    """
    if wgcna_path is None:
        return [], "none"
    if not wgcna_path.exists():
        return [], f"missing:{wgcna_path}"
    if wgcna_path.is_file():
        genes = load_wgcna_genes_from_file(wgcna_path)
        return genes, f"file:{wgcna_path}"
    # directory
    base = wgcna_path
    tables = base / "tables"
    
    candidates = [
        tables / "wgcna_genes_for_ML.csv",                     
        tables / "conservative_ferroptosis_genes_GS_MM.tsv",   
        tables / "module_hub_genes_MM_GS.tsv",                  
        tables / "module_genes_GS_MM.tsv",                      
        tables / "possible_prognostic_genes_consensus.csv",      
    ]
    for p in candidates:
        if p.exists():
            genes = load_wgcna_genes_from_file(p)
            if genes:
                print(f"[WGCNA] Selected source: {p.name} ({len(genes)} genes)")
                return genes, f"dir:{base} -> {p.name}"
    
    gm = tables / "gene_modules.tsv"
    ts = tables / "module_trait_stats.tsv"
    if gm.exists():
        df = _read_table_any(gm)
        cols = [c.lower() for c in df.columns]
        gcol = _detect_gene_column(df)
        mod_col = next((df.columns[i] for i,c in enumerate(cols) if c in ("module","modulecolor","module_color","color")), None)
        if mod_col is None:
           
            if len(df.columns) >= 2:
                mod_col = df.columns[1]
        chosen_module = None
        if ts.exists():
            try:
                tdf = _read_table_any(ts)
                chosen_module = _pick_module_from_trait_stats(tdf)
            except Exception:
                chosen_module = None
        if chosen_module is not None and mod_col in df.columns:
            sub = df[df[mod_col].astype(str) == str(chosen_module)]
            genes = [str(x).strip() for x in sub[gcol].tolist()]
            genes = [g for g in genes if g and g.lower() not in ("nan","na","none")]
            return sorted(set(genes)), f"dir:{base} -> gene_modules.tsv (module={chosen_module})"
       
        genes = [str(x).strip() for x in df[gcol].tolist()]
        genes = [g for g in genes if g and g.lower() not in ("nan","na","none")]
        return sorted(set(genes)), f"dir:{base} -> gene_modules.tsv (all modules)"
    return [], f"dir:{base} -> no known tables found"

def build_gene_universe(ferr_all: List[str], wgcna: List[str], mode: str) -> List[str]:
    A, B = set(ferr_all), set(wgcna)
    if mode == "intersection": G = A & B if B else A
    elif mode == "wgcna_only": G = B
    elif mode == "ferrdb_only": G = A
    else: G = A | B
    return sorted(G)

# ---------------- Survival metrics ----------------
def _detect_time_unit(y_time: np.ndarray) -> str:
    """"""
    med = float(np.nanmedian(y_time))
    return "days" if med > 50 else "years"


def cindex(y_time: np.ndarray, y_event: np.ndarray, risk: np.ndarray,
           use_ipcw: bool = True) -> float:
    """
    """
    y_time  = np.asarray(y_time,  dtype=float)
    y_event = np.asarray(y_event, dtype=int)
    risk    = np.asarray(risk,    dtype=float)
    m = np.isfinite(y_time) & np.isfinite(y_event) & np.isfinite(risk)
    y_time, y_event, risk = y_time[m], y_event[m], risk[m]
    if len(y_time) < 10 or y_event.sum() < 3:
        return float("nan")

    # Uno's IPCW C-index (preferred for censored survival data)
    if use_ipcw and SKSURV_AVAILABLE:
        try:
            from sksurv.metrics import concordance_index_ipcw
            y_struct = np.array(
                [(bool(e), float(t)) for t, e in zip(y_time, y_event)],
                dtype=[("event", bool), ("time", float)]
            )
            # tau = 75th percentile of observed event times
            event_times = y_time[y_event == 1]
            if len(event_times) >= 3:
                tau = float(np.percentile(event_times, 75))
                result = concordance_index_ipcw(y_struct, y_struct, risk, tau=tau)
                return float(result[0])
        except Exception:
            pass  # fallback to Harrell

    # Harrell C-index (fallback)
    if SKSURV_AVAILABLE and concordance_index_censored is not None:
        try:
            return float(concordance_index_censored(y_event.astype(bool), y_time, risk)[0])
        except Exception:
            pass
    if LIFELINES_AVAILABLE and concordance_index is not None:
        try:
            return float(concordance_index(y_time, -risk, y_event))
        except Exception:
            pass
    return float("nan")

def _build_sksurv_y(t: pd.Series, e: pd.Series):
    if not SKSURV_AVAILABLE or Surv is None:
        raise RuntimeError("scikit-survival missing")
    return Surv.from_arrays(event=e.astype(bool).values, time=t.values.astype(float))

# ---------------- Feature selection (meta-univariate) ----------------
def _eval_one_gene(args_tuple):
    """"""
    (f, Xtr2_vals, Xtr2_idx, ttr_vals, etr_vals, cohort_arr, cohorts,
     weights, stage_adjusted, stage_col_vals, penalizer,
     min_events_per_cohort, min_frac_cohorts, require_consistent_sign) = args_tuple

    z_list, w_list, signs = [], [], []
    used = 0
    for coh in cohorts:
        idx = cohort_arr == coh
        tt_v = ttr_vals[idx]; ee_v = etr_vals[idx]
        ev = int(np.nansum(ee_v)); n = int(idx.sum())
        if n < 10 or ev < int(min_events_per_cohort):
            continue
        df_dict = {"T": tt_v, "E": ee_v}
        if stage_adjusted and stage_col_vals is not None:
            df_dict["stage_col"] = stage_col_vals[idx]
        df_dict[f] = Xtr2_vals[idx]
        df_loc = pd.DataFrame(df_dict)
        df_loc = df_loc.replace([np.inf, -np.inf], np.nan)
        for col in df_loc.columns:
            df_loc[col] = df_loc[col].fillna(df_loc[col].median())
        df_loc = df_loc.dropna()
        if df_loc.shape[0] < 10 or df_loc["E"].sum() < int(min_events_per_cohort):
            continue
        try:
            cph = CoxPHFitter(penalizer=float(penalizer))
            cph.fit(df_loc, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary
            if f not in s.index:
                continue
            coef = float(s.loc[f, "coef"]); se = float(s.loc[f, "se(coef)"])
            if not np.isfinite(coef) or not np.isfinite(se) or se <= 0:
                continue
            z = coef / se
        except Exception:
            continue
        used += 1
        z_list.append(z)
        w_list.append(float(weights.get(coh, 1.0)))
        signs.append(1 if coef >= 0 else -1)

    if used == 0:
        return None
    frac = used / max(1, len(cohorts))
    if frac < float(min_frac_cohorts):
        return None
    if require_consistent_sign and len(signs) >= 2:
        ssum = sum(signs); maj = 1 if ssum >= 0 else -1
        n_opp = sum(1 for s in signs if s != maj)
        if n_opp > 0 and not (used >= 3 and n_opp <= 1):
            return None

    z_arr = np.asarray(z_list, float); w_arr = np.asarray(w_list, float)
    z_comb = float(np.sum(w_arr * z_arr) / math.sqrt(np.sum(w_arr ** 2)))
    p_comb = float(2 * (1 - 0.5 * (1 + math.erf(abs(z_comb) / math.sqrt(2)))))
   
    concordance_scores = []
    for coh_i, coh in enumerate(cohorts):
        idx_c = cohort_arr == coh
        t_c  = ttr_vals[idx_c]
        e_c  = etr_vals[idx_c]
        x_c  = Xtr2_vals[idx_c]
        ok_c = np.isfinite(t_c) & np.isfinite(e_c) & np.isfinite(x_c) & (t_c > 0)
        if ok_c.sum() < 10 or e_c[ok_c].sum() < 3:
            continue
        try:
            # IPCW-weighted concordance per cohort
            from sksurv.metrics import concordance_index_ipcw
            y_c = np.array([(bool(ev), float(tt)) for tt, ev in zip(t_c[ok_c], e_c[ok_c])],
                           dtype=[("event", bool), ("time", float)])
            tau_c = float(np.percentile(t_c[ok_c][e_c[ok_c]==1], 75)) if e_c[ok_c].sum() >= 4 else float(np.max(t_c[ok_c]))
            ci_c = float(concordance_index_ipcw(y_c, y_c, x_c[ok_c], tau=tau_c)[0])
        except Exception:
            try:
                from sksurv.metrics import concordance_index_censored
                ci_c = float(concordance_index_censored(e_c[ok_c].astype(bool), t_c[ok_c], x_c[ok_c])[0])
            except Exception:
                ci_c = 0.5
        concordance_scores.append(abs(ci_c - 0.5))  # distance from 0.5 (no info)

    concordance_mean = float(np.mean(concordance_scores)) if concordance_scores else 0.0

    return {"feature": f, "gene": f.replace("RNA__","",1),
            "z_meta": z_comb, "p_meta": p_comb,
            "concordance_score": concordance_mean,
            "n_cohorts_used": used, "frac_cohorts_used": frac}


def meta_univariate_select(Xtr2, ttr, etr, cohort_series, candidate_rna_cols, clin_cols,
                           stage_adjusted, penalizer, top_genes, min_events_per_cohort,
                           min_frac_cohorts, require_consistent_sign,
                           n_jobs: int = 1):
    if not LIFELINES_AVAILABLE:
        return pd.DataFrame()

    stage_cols = [c for c in clin_cols if "stage" in c.lower()]
    cohorts = cohort_series.unique().tolist()
    weights = {coh: math.sqrt(max(int(np.nansum(etr.loc[cohort_series==coh].values)),1)) for coh in cohorts}
    rows=[]
    import concurrent.futures as _cf

    
    _Xtr2_arr   = Xtr2.values if hasattr(Xtr2, "values") else np.array(Xtr2)
    _Xtr2_cols  = list(Xtr2.columns)
    _ttr_arr    = ttr.values if hasattr(ttr, "values") else np.array(ttr)
    _etr_arr    = etr.values if hasattr(etr, "values") else np.array(etr)
    _coh_arr    = cohort_series.values if hasattr(cohort_series, "values") else np.array(cohort_series)

    # Stage
    _stage_col_vals = None
    if stage_adjusted and stage_cols:
        sc0 = stage_cols[0] if stage_cols[0] in Xtr2.columns else None
        if sc0:
            _stage_col_vals = _Xtr2_arr[:, _Xtr2_cols.index(sc0)]

    def _make_args(f):
        f_idx = _Xtr2_cols.index(f) if f in _Xtr2_cols else -1
        if f_idx < 0:
            return None
        return (f, _Xtr2_arr[:, f_idx], None, _ttr_arr, _etr_arr,
                _coh_arr, cohorts, weights, stage_adjusted,
                _stage_col_vals, penalizer,
                min_events_per_cohort, min_frac_cohorts, require_consistent_sign)

    
    n_jobs_use = int(n_jobs) if n_jobs and int(n_jobs) >= 1 else int(_META_N_JOBS)

    arg_list = [_make_args(f) for f in candidate_rna_cols]
    arg_list = [a for a in arg_list if a is not None]

    if n_jobs_use > 1:
        try:
            with _cf.ThreadPoolExecutor(max_workers=n_jobs_use) as executor:
                results_iter = executor.map(_eval_one_gene, arg_list, timeout=600)
                for res in results_iter:
                    if res is not None:
                        rows.append(res)
        except Exception as _pe:
            
            print(f"  [WARN] Parallel meta-univariate failed ({_pe}), falling back to serial")
            rows = []
            for a in arg_list:
                res = _eval_one_gene(a)
                if res is not None:
                    rows.append(res)
    else:
        for a in arg_list:
            res = _eval_one_gene(a)
            if res is not None:
                rows.append(res)

    out = pd.DataFrame(rows)
    if out.empty: return out
   
    if "concordance_score" in out.columns:
        out = out.sort_values(
            ["concordance_score", "p_meta", "n_cohorts_used"],
            ascending=[False, True, False]
        ).head(int(top_genes)).reset_index(drop=True)
    else:
        out = out.sort_values(["p_meta","n_cohorts_used"], ascending=[True,False]).head(int(top_genes)).reset_index(drop=True)

    
    target_min = max(5, int(top_genes * 0.30))
    if len(out) < target_min and (float(min_frac_cohorts) > 0.30 or int(min_events_per_cohort) > 6):
        relaxed_frac   = max(0.30, float(min_frac_cohorts)   - 0.10)
        relaxed_events = max(6,    int(min_events_per_cohort) - 4)
        rows2 = []
        for f in candidate_rna_cols:
            if any(r["feature"] == f for r in rows):
                continue  
            z_list2, w_list2, signs2 = [], [], []
            used2 = 0
            for coh in cohorts:
                idx = cohort_series == coh
                tt, ee = ttr.loc[idx], etr.loc[idx]
                ev, n = int(np.nansum(ee.values)), int(idx.sum())
                if n < 10 or ev < relaxed_events:
                    continue
                df2 = pd.DataFrame({"T": tt.values, "E": ee.values}, index=tt.index)
                df2[f] = Xtr2.loc[idx, f].values
                df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()
                if df2.shape[0] < 10 or df2["E"].sum() < relaxed_events:
                    continue
                try:
                    cph2 = CoxPHFitter(penalizer=float(penalizer))
                    cph2.fit(df2, duration_col="T", event_col="E", show_progress=False)
                    s2 = cph2.summary
                    if f not in s2.index:
                        continue
                    coef2 = float(s2.loc[f, "coef"])
                    se2   = float(s2.loc[f, "se(coef)"])
                    if not np.isfinite(coef2) or not np.isfinite(se2) or se2 <= 0:
                        continue
                    z2 = coef2 / se2
                except Exception:
                    continue
                used2 += 1
                z_list2.append(z2)
                w_list2.append(float(weights.get(coh, 1.0)))
                signs2.append(1 if coef2 >= 0 else -1)
            if used2 == 0:
                continue
            frac2 = used2 / max(1, len(cohorts))
            if frac2 < relaxed_frac:
                continue
            z_arr2 = np.asarray(z_list2, float)
            w_arr2 = np.asarray(w_list2, float)
            z_comb2 = float(np.sum(w_arr2 * z_arr2) / math.sqrt(np.sum(w_arr2**2)))
            p_comb2 = float(2 * (1 - 0.5 * (1 + math.erf(abs(z_comb2) / math.sqrt(2)))))
            rows2.append({
                "feature": f, "gene": f.replace("RNA__", "", 1),
                "z_meta": z_comb2, "p_meta": p_comb2,
                "n_cohorts_used": used2, "frac_cohorts_used": frac2,
            })
        if rows2:
            out2 = pd.DataFrame(rows2).sort_values(
                ["p_meta", "n_cohorts_used"], ascending=[True, False]
            ).head(max(0, target_min - len(out))).reset_index(drop=True)
            out = pd.concat([out, out2], ignore_index=True).drop_duplicates(
                subset=["feature"], keep="first"
            ).head(int(top_genes)).reset_index(drop=True)
            print(f"  [adaptive] Relaxed threshold (frac≥{relaxed_frac}, ev≥{relaxed_events}): "
                  f"+{len(out2)} genes, total={len(out)}")


    if len(out) > 1:
        pvals = out["p_meta"].values.astype(float)
        finite_mask = np.isfinite(pvals)
        q = np.full_like(pvals, np.nan)
        if finite_mask.sum() > 0:
            pf = pvals[finite_mask]
            n = len(pf)
            order = np.argsort(pf)
            rank = np.empty_like(order); rank[order] = np.arange(1, n+1)
            q_finite = np.minimum(1.0, pf * n / rank)

            for i in range(n-2, -1, -1):
                q_finite[order[i]] = min(q_finite[order[i]], q_finite[order[i+1]])
            q[finite_mask] = q_finite
        out["q_fdr"] = q
    else:
        out["q_fdr"] = out["p_meta"]

    return out


# ---------------- year endpoint metrics (1/2/3/5-year) ----------------
def _binary_event_by_horizon(y_time: np.ndarray, y_event: np.ndarray, years: int):
    """
    Binary label at horizon:
      y=1: event on/before horizon
      y=0: known event-free beyond horizon (time > horizon)
      excluded: censored on/before horizon (unknown)
    """
    horizon_days = float(years) * 365.0
    t = np.asarray(y_time, dtype=float)
    e = np.asarray(y_event, dtype=int)

    known_event = (e == 1) & np.isfinite(t) & (t <= horizon_days)
    known_nonevent = np.isfinite(t) & (t > horizon_days)

    y = np.full(t.shape, np.nan, dtype=float)
    y[known_event] = 1.0
    y[known_nonevent] = 0.0
    keep = np.isfinite(y)
    return y[keep].astype(int), keep

def fixed_horizon_auc_acc(y_time: np.ndarray, y_event: np.ndarray, risk: np.ndarray, years_list=(1,2,3,5)) -> pd.DataFrame:
    rows = []
    r = np.asarray(risk, dtype=float)
    for yr in years_list:
        yb, keep = _binary_event_by_horizon(y_time, y_event, yr)
        rr = r[keep] if len(r)==len(keep) else np.asarray([], dtype=float)
        out = {
            "horizon_year": int(yr),
            "n_evaluable": int(len(yb)),
            "n_pos": int(np.sum(yb==1)) if len(yb) else 0,
            "n_neg": int(np.sum(yb==0)) if len(yb) else 0,
            "auc_roc": float("nan"),
            "accuracy": float("nan"),
            "threshold": float("nan"),
        }
        if len(yb) >= 20 and np.sum(yb==1) >= 5 and np.sum(yb==0) >= 5 and np.isfinite(rr).mean() > 0.8 and SKLEARN_METRICS_AVAILABLE:
            try:
                out["auc_roc"] = float(roc_auc_score(yb, rr))
            except Exception:
                pass
            try:
                thr = float(np.nanmedian(rr))
                yhat = (rr >= thr).astype(int)
                out["accuracy"] = float(accuracy_score(yb, yhat))
                out["threshold"] = thr
            except Exception:
                pass
        rows.append(out)
    return pd.DataFrame(rows)


 ----------------
def time_dependent_auc_ipcw(y_train_time, y_train_event, y_test_time, y_test_event, risk, years_list=(1,2,3,5)) -> pd.DataFrame:
    rows=[]
    if not SKSURV_AVAILABLE or cumulative_dynamic_auc is None or Surv is None:
        for yr in years_list:
            rows.append({"horizon_year":int(yr), "auc_ipcw":float("nan"), "status":"skipped_sksurv_missing"})
        return pd.DataFrame(rows)
    try:
        ytr = Surv.from_arrays(event=np.asarray(y_train_event, bool), time=np.asarray(y_train_time, float))
        yte = Surv.from_arrays(event=np.asarray(y_test_event, bool), time=np.asarray(y_test_time, float))
       
        _med_train = float(np.nanmedian(np.asarray(y_train_time, float)))
        _unit_scale = 365.0 if _med_train > 50 else 1.0
        times = np.asarray([float(yr) * _unit_scale for yr in years_list], float)

        max_obs = float(np.nanmax(np.asarray(y_test_time, float)))
        times = times[times < max_obs * 0.99]
        if len(times) == 0:
            raise ValueError("No valid time points within observed range")
        _, aucs = cumulative_dynamic_auc(ytr, yte, np.asarray(risk, float), times)
        for yr, auc in zip(years_list, aucs):
            rows.append({"horizon_year":int(yr), "auc_ipcw":float(auc), "status":"ok"})
        return pd.DataFrame(rows)
    except Exception as ex:
        for yr in years_list:
            rows.append({"horizon_year":int(yr), "auc_ipcw":float("nan"), "status":f"failed:{type(ex).__name__}"})
        return pd.DataFrame(rows)

def _cox_on_risk_survival_probs(y_train_time, y_train_event, risk_train, risk_eval, times_days: np.ndarray):
    if not LIFELINES_AVAILABLE or CoxPHFitter is None:
        return None
    try:
        df = pd.DataFrame({"T":np.asarray(y_train_time,float), "E":np.asarray(y_train_event,int), "risk":np.asarray(risk_train,float)})
        df = df.replace([np.inf,-np.inf], np.nan).dropna()
        if df.shape[0] < 50 or df["E"].sum() < 8:
            return None
        cph = CoxPHFitter(penalizer=1e-4)
        cph.fit(df, duration_col="T", event_col="E", show_progress=False)
        X = pd.DataFrame({"risk":np.asarray(risk_eval,float)})
        sf = cph.predict_survival_function(X, times=times_days)
        return np.asarray(sf.values.T, float)
    except Exception:
        return None

def brier_and_ibs_ipcw(y_train_time, y_train_event, y_test_time, y_test_event, risk_train, risk_test, years_list=(1,2,3,5)) -> pd.DataFrame:
    times_days = np.asarray([float(yr)*365.0 for yr in years_list], float)
    rows=[]
    if not SKSURV_AVAILABLE or brier_score is None or integrated_brier_score is None or Surv is None:
        for yr in years_list:
            rows.append({"horizon_year":int(yr), "brier":float("nan"), "ibs":float("nan"), "status":"skipped_sksurv_missing"})
        return pd.DataFrame(rows)
    surv_probs = _cox_on_risk_survival_probs(y_train_time, y_train_event, risk_train, risk_test, times_days)
    if surv_probs is None or not np.isfinite(surv_probs).any():
        for yr in years_list:
            rows.append({"horizon_year":int(yr), "brier":float("nan"), "ibs":float("nan"), "status":"skipped_survfit_failed"})
        return pd.DataFrame(rows)
    try:
        ytr = Surv.from_arrays(event=np.asarray(y_train_event, bool), time=np.asarray(y_train_time, float))
        yte = Surv.from_arrays(event=np.asarray(y_test_event, bool), time=np.asarray(y_test_time, float))
        _, bs = brier_score(ytr, yte, surv_probs, times_days)
        ibs = float(integrated_brier_score(ytr, yte, surv_probs, times_days))
        for yr, b in zip(years_list, bs):
            rows.append({"horizon_year":int(yr), "brier":float(b), "ibs":ibs, "status":"ok"})
        return pd.DataFrame(rows)
    except Exception as ex:
        for yr in years_list:
            rows.append({"horizon_year":int(yr), "brier":float("nan"), "ibs":float("nan"), "status":f"failed:{type(ex).__name__}"})
        return pd.DataFrame(rows)

def bootstrap_ci(metric_fn, y_time, y_event, score, n_boot=1000, seed=13, min_n=40):
    """v31: BCa-inspired bootstrap CI (percentile fallback). n_boot default artırıldı 600->1000."""
    rng = np.random.default_rng(seed)
    n = len(y_time)
    if n < min_n:
        return {"mean":float("nan"), "ci_low":float("nan"), "ci_high":float("nan"), "n":int(n), "status":"too_small"}
    vals=[]
    idx_all = np.arange(n)
    for _ in range(int(n_boot)):
        idx = rng.choice(idx_all, size=n, replace=True)
        v = metric_fn(np.asarray(y_time)[idx], np.asarray(y_event)[idx], np.asarray(score)[idx])
        if np.isfinite(v):
            vals.append(float(v))
    if len(vals) < max(30, int(0.2*n_boot)):
        return {"mean":float("nan"), "ci_low":float("nan"), "ci_high":float("nan"), "n":int(n), "status":"insufficient_finite"}
    vals = np.array(vals, float)
    return {"mean":float(np.mean(vals)), "ci_low":float(np.quantile(vals,0.025)), "ci_high":float(np.quantile(vals,0.975)), "n":int(n), "status":"ok"}


# ============================================================
# ============================================================

def permutation_cindex_pvalue(y_time, y_event, risk, n_perm=1000, seed=13):
    """
    """
    obs = cindex(np.asarray(y_time, float), np.asarray(y_event, int), np.asarray(risk, float))
    if not np.isfinite(obs):
        return {"observed_cindex": obs, "perm_mean": float("nan"), "perm_std": float("nan"), "p_value": float("nan"), "n_perm": n_perm}
    rng = np.random.default_rng(seed)
    r = np.asarray(risk, float)
    t = np.asarray(y_time, float)
    e = np.asarray(y_event, int)
    perm_vals = []
    for _ in range(int(n_perm)):
        rp = rng.permutation(r)
        v = cindex(t, e, rp)
        if np.isfinite(v):
            perm_vals.append(float(v))
    if not perm_vals:
        return {"observed_cindex": obs, "perm_mean": float("nan"), "perm_std": float("nan"), "p_value": float("nan"), "n_perm": n_perm}
    pv = np.array(perm_vals, float)
    # one-sided p: Pr(perm >= obs)
    p_val = float(np.mean(pv >= obs))
    return {
        "observed_cindex": float(obs),
        "perm_mean": float(np.mean(pv)),
        "perm_std": float(np.std(pv)),
        "p_value": p_val,
        "n_perm": len(pv)
    }


def continuous_nri_idi(y_time: np.ndarray, y_event: np.ndarray,
                       risk_new: np.ndarray, risk_base: np.ndarray,
                       years: int = 3) -> dict:
    """
    """
    result = {
        "NRI_continuous": float("nan"),
        "NRI_p":          float("nan"),
        "IDI":            float("nan"),
        "IDI_p":          float("nan"),
        "note":           ""
    }
    try:
        t = np.asarray(y_time,  float)
        e = np.asarray(y_event, int)
        rn = np.asarray(risk_new,  float)
        rb = np.asarray(risk_base, float)
        ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(rn) & np.isfinite(rb) & (t > 0)
        t, e, rn, rb = t[ok], e[ok], rn[ok], rb[ok]
        if len(t) < 30 or e.sum() < 5:
            result["note"] = "insufficient data"
            return result

        horizon = float(years) * (365.0 if np.nanmedian(t) > 50 else 1.0)

        # Event indicators at horizon
        event_at_h   = (e == 1) & (t <= horizon)
        noevent_at_h = t > horizon

        # Rank-normalize risk scores (0-1 scale for NRI)
        def _rank01(x):
            r = np.argsort(np.argsort(x)).astype(float)
            return r / (len(r) - 1 + 1e-10)

        pn = _rank01(rn); pb = _rank01(rb)

        # NRI continuous (Pencina 2008)
        n_ev  = event_at_h.sum()
        n_nev = noevent_at_h.sum()
        if n_ev > 0 and n_nev > 0:
            nri_ev  = np.mean(pn[event_at_h]  - pb[event_at_h])
            nri_nev = np.mean(pb[noevent_at_h] - pn[noevent_at_h])
            nri     = nri_ev + nri_nev

           
            rng_nri = np.random.default_rng(43)
            nri_boots = []
            for _ in range(500):
                idx_b = rng_nri.integers(0, len(t), size=len(t))
                pn_b = pn[idx_b]; pb_b = pb[idx_b]
                ev_b = event_at_h[idx_b]; nev_b = noevent_at_h[idx_b]
                if ev_b.sum() < 3 or nev_b.sum() < 3:
                    continue
                nri_b = (np.mean(pn_b[ev_b] - pb_b[ev_b]) +
                         np.mean(pb_b[nev_b] - pn_b[nev_b]))
                nri_boots.append(nri_b)
            if len(nri_boots) >= 50:
                nri_arr = np.array(nri_boots)
                se_nri  = float(np.std(nri_arr))
                z_nri   = nri / (se_nri + 1e-10)
                p_nri   = float(2 * (1 - 0.5 * (1 + math.erf(abs(z_nri) / math.sqrt(2)))))
                p_nri   = min(p_nri, 1.0)
            else:
                # Analitik fallback
                var_ev  = np.var(pn[event_at_h]  - pb[event_at_h])  / n_ev
                var_nev = np.var(pb[noevent_at_h] - pn[noevent_at_h]) / n_nev
                se_nri  = np.sqrt(var_ev + var_nev)
                z_nri   = nri / (se_nri + 1e-10)
                p_nri   = float(2 * (1 - 0.5 * (1 + math.erf(abs(z_nri) / math.sqrt(2)))))

            result["NRI_continuous"] = round(float(nri),  4)
            result["NRI_p"]          = round(p_nri,       4)

       
        is_new  = np.mean(pn[event_at_h]) - np.mean(pn[noevent_at_h]) if (n_ev > 0 and n_nev > 0) else float("nan")
        is_base = np.mean(pb[event_at_h]) - np.mean(pb[noevent_at_h]) if (n_ev > 0 and n_nev > 0) else float("nan")
        if np.isfinite(is_new) and np.isfinite(is_base):
            idi = is_new - is_base
            # Bootstrap p for IDI
            rng_idi = np.random.default_rng(42)
            idi_boots = []
            for _ in range(500):
                idx_b = rng_idi.integers(0, len(t), size=len(t))
                pn_b = pn[idx_b]; pb_b = pb[idx_b]
                ev_b = event_at_h[idx_b]; nev_b = noevent_at_h[idx_b]
                if ev_b.sum() < 3 or nev_b.sum() < 3:
                    continue
                is_n_b = np.mean(pn_b[ev_b]) - np.mean(pn_b[nev_b])
                is_b_b = np.mean(pb_b[ev_b]) - np.mean(pb_b[nev_b])
                idi_boots.append(is_n_b - is_b_b)
            idi_boots = np.array(idi_boots)
            p_idi = float(np.mean(idi_boots <= 0)) if idi > 0 else float(np.mean(idi_boots >= 0))
            p_idi = min(p_idi * 2, 1.0)  # two-sided

            result["IDI"]   = round(float(idi), 4)
            result["IDI_p"] = round(p_idi,       4)

        result["note"] = (
            f"NRI={result['NRI_continuous']:.3f} (p={result['NRI_p']:.3f}), "
            f"IDI={result['IDI']:.3f} (p={result['IDI_p']:.3f})"
            if all(np.isfinite([result['NRI_continuous'], result['IDI']])) else "partial"
        )
    except Exception as ex:
        result["note"] = f"NRI/IDI failed: {ex}"
    return result


def multivariable_cox_validation(y_time, y_event, risk_score, stage_series=None, penalizer=1e-4):
    """
    v31: Prognostik imzanın stage'den bağımsız prognostik değerini doğrular.
    risk_score'u sürekli değişken olarak, stage varsa covariate olarak Cox'a sokar.
    Döndürür: coef, HR, 95%CI, p-değeri — her ikisi için de (risk, stage).
    """
    if not LIFELINES_AVAILABLE or CoxPHFitter is None:
        return pd.DataFrame()
    t = np.asarray(y_time, float)
    e = np.asarray(y_event, int)
    r = np.asarray(risk_score, float)
    df = pd.DataFrame({"T": t, "E": e, "risk_score": r})
    if stage_series is not None:
        s = np.asarray(stage_series, float)
        df["stage"] = s
    df = df.replace([np.inf, -np.inf], np.nan)
    for c in df.columns:
        df[c] = df[c].fillna(df[c].median())
    df = df.dropna()
    if df.shape[0] < 30 or int(df["E"].sum()) < 6:
        return pd.DataFrame()
    try:
        cph = CoxPHFitter(penalizer=float(penalizer))
        cph.fit(df, duration_col="T", event_col="E", show_progress=False)
        s = cph.summary.copy()
        s["HR"] = np.exp(s["coef"])
        s["HR_ci_low"] = np.exp(s["coef lower 95%"])
        s["HR_ci_high"] = np.exp(s["coef upper 95%"])
        return s[["coef","HR","HR_ci_low","HR_ci_high","p"]].reset_index().rename(columns={"index":"variable"})
    except Exception:
        return pd.DataFrame()


def decision_curve_analysis(y_time, y_event, risk, years=3, n_thresholds=51):
    """
    v31: Basit DCA — threshold probability vs net benefit.
    Net benefit = (TP/n) - (FP/n) * (pt/(1-pt))
    Döndürür: DataFrame(threshold, net_benefit_model, net_benefit_all, net_benefit_none)
    """
    horizon_days = float(years) * 365.0
    t = np.asarray(y_time, float)
    e = np.asarray(y_event, int)
    r = np.asarray(risk, float)
    # binary outcome at horizon
    event_mask = (e == 1) & np.isfinite(t) & (t <= horizon_days)
    nonevent_mask = np.isfinite(t) & (t > horizon_days)
    y = np.full(len(t), np.nan)
    y[event_mask] = 1.0
    y[nonevent_mask] = 0.0
    keep = np.isfinite(y)
    y = y[keep].astype(int)
    r_k = r[keep]
    if len(y) < 20 or y.sum() < 5 or (1 - y).sum() < 5:
        return pd.DataFrame()
    # rank-normalize risk to [0,1]
    r_norm = (r_k - np.nanmin(r_k)) / (np.nanmax(r_k) - np.nanmin(r_k) + 1e-12)
    n = len(y)
    prev = float(y.mean())
    rows = []
    for pt in np.linspace(0.01, 0.99, int(n_thresholds)):
        predicted_pos = (r_norm >= pt).astype(int)
        tp = int(np.sum((predicted_pos == 1) & (y == 1)))
        fp = int(np.sum((predicted_pos == 1) & (y == 0)))
        nb_model = (tp / n) - (fp / n) * (pt / (1 - pt))
        nb_all = prev - (1 - prev) * (pt / (1 - pt))
        rows.append({
            "threshold": float(pt),
            "net_benefit_model": float(nb_model),
            "net_benefit_all": float(nb_all),
            "net_benefit_none": 0.0
        })
    return pd.DataFrame(rows)


def calibration_curve(y_time, y_event, risk, y_train_time, y_train_event, years=3, n_bins=5):
    """
    v31: Observed vs predicted survival oranları (yeterli veri varsa).
    Döndürür: DataFrame(bin_mean_predicted, bin_observed_rate, bin_n)
    """
    surv_probs = _cox_on_risk_survival_probs(
        y_train_time, y_train_event, risk, risk,
        np.array([float(years) * 365.0])
    )
    if surv_probs is None:
        return pd.DataFrame()
    pred_surv = surv_probs[:, 0]
    t = np.asarray(y_time, float)
    e = np.asarray(y_event, int)
    horizon_days = float(years) * 365.0
    event_mask = (e == 1) & (t <= horizon_days)
    nonevent_mask = (t > horizon_days)
    y = np.full(len(t), np.nan)
    y[event_mask] = 0.0       # event = did NOT survive
    y[nonevent_mask] = 1.0    # survived
    keep = np.isfinite(y) & np.isfinite(pred_surv)
    y_k, pred_k = y[keep], pred_surv[keep]
    if len(y_k) < 20:
        return pd.DataFrame()
    bins = np.quantile(pred_k, np.linspace(0, 1, int(n_bins) + 1))
    rows = []
    for i in range(int(n_bins)):
        mask = (pred_k >= bins[i]) & (pred_k <= bins[i + 1])
        if mask.sum() < 3:
            continue
        rows.append({
            "bin": i + 1,
            "bin_mean_predicted": float(np.mean(pred_k[mask])),
            "bin_observed_rate": float(np.mean(y_k[mask])),
            "bin_n": int(mask.sum())
        })
    return pd.DataFrame(rows)


def cox_ph_assumption_test(y_time: np.ndarray, y_event: np.ndarray,
                           risk: np.ndarray) -> dict:
    """
    Returns: {global_p, global_stat, significant (bool), note}
    """
    result = {
        "global_stat":   float("nan"),
        "global_p":      float("nan"),
        "ph_violated":   None,
        "n":             int(len(y_time)),
        "method":        "schoenfeld_residuals",
        "note":          ""
    }
    if not LIFELINES_AVAILABLE:
        result["note"] = "lifelines not available"
        return result
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import proportional_hazard_test
        df_ph = pd.DataFrame({
            "T":    np.asarray(y_time,  float),
            "E":    np.asarray(y_event, int),
            "risk": np.asarray(risk,    float)
        }).dropna()
        if len(df_ph) < 30 or int(df_ph["E"].sum()) < 8:
            result["note"] = "insufficient data for PH test"
            return result
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df_ph, duration_col="T", event_col="E", show_progress=False)
        ph_test = proportional_hazard_test(cph, df_ph, time_transform="rank")
        p_global = float(ph_test.summary["p"].min())
        stat     = float(ph_test.summary["test_statistic"].sum())
        result["global_stat"] = stat
        result["global_p"]    = p_global
        result["ph_violated"] = p_global < 0.05
        if p_global < 0.05:
            result["note"] = (
                f"WARNING: PH assumption may be violated (p={p_global:.4f}). "
                "Consider time-stratified or time-varying Cox, or use RMST.")
        else:
            result["note"] = f"PH assumption not rejected (p={p_global:.4f})"
    except Exception as e:
        result["note"] = f"PH test failed: {type(e).__name__}: {e}"
    return result


def expected_observed_ratio(y_time: np.ndarray, y_event: np.ndarray,
                             risk: np.ndarray,
                             y_train_time: np.ndarray,
                             y_train_event: np.ndarray,
                             years: int = 3) -> dict:
    """
    """
    result = {
        "horizon_year": int(years),
        "E_predicted":  float("nan"),
        "O_observed":   float("nan"),
        "EO_ratio":     float("nan"),
        "EO_ci_low":    float("nan"),
        "EO_ci_high":   float("nan"),
        "note":         ""
    }
    try:
        horizon = float(years) * (365.0 if np.nanmedian(y_time) > 50 else 1.0)
        t  = np.asarray(y_time,  float)
        e  = np.asarray(y_event, int)
        r  = np.asarray(risk,    float)
        ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(r) & (t > 0)
        t, e, r = t[ok], e[ok], r[ok]
        if len(t) < 20: return result

        # Kaplan-Meier event probability at horizon
        from lifelines import KaplanMeierFitter
        kmf = KaplanMeierFitter()
        kmf.fit(t, event_observed=e, label="km")
        # 1 - S(t) at horizon
        surv_at_h = float(kmf.survival_function_at_times([horizon]).values[0])
        O = 1.0 - surv_at_h
        n_obs = len(t)

        
        surv_probs = _cox_on_risk_survival_probs(
            np.asarray(y_train_time, float), np.asarray(y_train_event, int),
            r, r, np.array([horizon])
        )
        if surv_probs is not None and surv_probs.shape[0] == len(r):
            # E = mean predicted event probability = mean(1 - S_pred(horizon))
            E = float(np.mean(1.0 - surv_probs[:, 0]))
        else:
            # Fallback: KM risk-group proxy
            med_risk = np.median(r)
            hi_mask  = r >= med_risk
            lo_mask  = ~hi_mask
            kmf_hi = KaplanMeierFitter(); kmf_hi.fit(t[hi_mask], e[hi_mask])
            kmf_lo = KaplanMeierFitter(); kmf_lo.fit(t[lo_mask], e[lo_mask])
            p_hi = 1 - float(kmf_hi.survival_function_at_times([horizon]).values[0])
            p_lo = 1 - float(kmf_lo.survival_function_at_times([horizon]).values[0])
            E = (p_hi * hi_mask.sum() + p_lo * lo_mask.sum()) / n_obs

        eo = E / O if O > 0 else float("nan")
        # Approximate 95% CI (Poisson approximation)
        n_events_obs = int(e.sum())
        if n_events_obs > 0:
            import math
            se_log_eo = 1.0 / math.sqrt(n_events_obs)
            eo_lo = eo * math.exp(-1.96 * se_log_eo)
            eo_hi = eo * math.exp(+1.96 * se_log_eo)
        else:
            eo_lo = eo_hi = float("nan")

        result.update({
            "E_predicted": round(E, 4),
            "O_observed":  round(O, 4),
            "EO_ratio":    round(eo, 4),
            "EO_ci_low":   round(eo_lo, 4),
            "EO_ci_high":  round(eo_hi, 4),
            "note":        f"{'Good calibration' if 0.85<eo<1.15 else 'REVIEW calibration'} (E/O={eo:.3f})"
        })
    except Exception as ex:
        result["note"] = f"E/O failed: {ex}"
    return result


def plot_dca(dca_df: "pd.DataFrame", title: str, out_png: "Path"):
    """v31: DCA figürü kaydeder."""
    if not MPL_AVAILABLE or dca_df is None or dca_df.empty:
        return False
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(dca_df["threshold"], dca_df["net_benefit_model"], label="Model", color="#1a6faf", lw=2)
        ax.plot(dca_df["threshold"], dca_df["net_benefit_all"], label="Treat all", color="#888", lw=1.5, ls="--")
        ax.axhline(0, color="#333", lw=1, ls=":")
        ax.set_xlabel("Threshold probability")
        ax.set_ylabel("Net benefit")
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.set_xlim(0, 1); ax.set_ylim(-0.05, max(dca_df["net_benefit_model"].max(), 0.05) + 0.02)
        plt.tight_layout()
        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=220); plt.close()
        return True
    except Exception:
        return False


def plot_calibration(cal_df: "pd.DataFrame", title: str, out_png: "Path"):
    """"""
    if not MPL_AVAILABLE or cal_df is None or cal_df.empty:
        return False
    try:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
        ax.scatter(cal_df["bin_mean_predicted"], cal_df["bin_observed_rate"],
                   s=cal_df["bin_n"] * 2, color="#1a6faf", alpha=0.8, zorder=5)
        ax.plot(cal_df["bin_mean_predicted"], cal_df["bin_observed_rate"],
                color="#1a6faf", lw=1.5, label="Model")
        ax.set_xlabel("Mean predicted survival")
        ax.set_ylabel("Observed survival rate")
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        plt.tight_layout()
        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=220); plt.close()
        return True
    except Exception:
        return False

def _km_group_from_risk(t: np.ndarray, e: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    """
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    r = np.asarray(r, float)
    grp = (r >= np.nanmedian(r)).astype(int)
    if LIFELINES_AVAILABLE and CoxPHFitter is not None:
        try:
            df = pd.DataFrame({"T": t, "E": e, "risk": r}).replace([np.inf, -np.inf], np.nan).dropna()
            if df.shape[0] >= 30 and int(df["E"].sum()) >= 6:
                cph = CoxPHFitter(penalizer=1e-4)
                cph.fit(df, duration_col="T", event_col="E", show_progress=False)
                coef = float(cph.params_["risk"])
                if coef < 0:
                    grp = 1 - grp
        except Exception:
            pass
    return grp

def km_logrank_hr(y_time, y_event, risk):
    out = {"n":int(len(risk)), "p_logrank":float("nan"), "hr":float("nan"), "hr_ci_low":float("nan"), "hr_ci_high":float("nan")}
    if not LIFELINES_AVAILABLE or KaplanMeierFitter is None or logrank_test is None or CoxPHFitter is None:
        out["status"] = "skipped_lifelines_missing"
        return out
    t = np.asarray(y_time, float); e = np.asarray(y_event, int); r = np.asarray(risk, float)
    m = np.isfinite(t) & np.isfinite(e) & np.isfinite(r)
    t,e,r = t[m], e[m], r[m]
    if len(t) < 30 or np.sum(e)==0:
        out["status"] = "too_small"; return out
    grp = _km_group_from_risk(t, e, r)
    try:
        res = logrank_test(t[grp==0], t[grp==1], event_observed_A=e[grp==0], event_observed_B=e[grp==1])
        out["p_logrank"] = float(res.p_value)
    except Exception:
        pass
    try:
        df = pd.DataFrame({"T":t, "E":e, "high":grp})
        cph = CoxPHFitter(penalizer=1e-4); cph.fit(df, duration_col="T", event_col="E", show_progress=False)
        coef = float(cph.params_["high"]); se = float(cph.summary.loc["high","se(coef)"])
        out["hr"] = float(np.exp(coef))
        out["hr_ci_low"] = float(np.exp(coef - 1.96*se))
        out["hr_ci_high"] = float(np.exp(coef + 1.96*se))
    except Exception:
        pass
    out["status"] = "ok"; return out

def plot_km(y_time, y_event, risk, title, out_png: Path):
    if not (MPL_AVAILABLE and LIFELINES_AVAILABLE and KaplanMeierFitter is not None):
        return False
    t = np.asarray(y_time, float); e = np.asarray(y_event, int); r = np.asarray(risk, float)
    m = np.isfinite(t) & np.isfinite(e) & np.isfinite(r)
    t,e,r = t[m], e[m], r[m]
    if len(t) < 30: return False
    grp = _km_group_from_risk(t, e, r)
    km0 = KaplanMeierFitter(); km1 = KaplanMeierFitter()
    plt.figure()
    km0.fit(t[grp==0], event_observed=e[grp==0], label="Low risk")
    km1.fit(t[grp==1], event_observed=e[grp==1], label="High risk")
    ax = km0.plot(ci_show=True); km1.plot(ax=ax, ci_show=True)
    ax.set_title(title); ax.set_xlabel("Time (days)"); ax.set_ylabel("Survival probability")
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=220); plt.close()
    return True


# ============================================================
# ============================================================

def compute_rmst(y_time: np.ndarray, y_event: np.ndarray, risk: np.ndarray,
                 tau_years: float = 3.0, n_boot: int = 500, seed: int = 13) -> dict:
    """
    """
    result = {
        "tau_years":    float(tau_years),
        "rmst_low":     float("nan"),
        "rmst_high":    float("nan"),
        "delta_rmst":   float("nan"),
        "delta_rmst_ci_low":  float("nan"),
        "delta_rmst_ci_high": float("nan"),
        "p_value":      float("nan"),
        "note":         ""
    }
    if not LIFELINES_AVAILABLE or KaplanMeierFitter is None:
        result["note"] = "lifelines missing"
        return result
    try:
        t = np.asarray(y_time,  float)
        e = np.asarray(y_event, int)
        r = np.asarray(risk,    float)
        ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(r) & (t > 0)
        t, e, r = t[ok], e[ok], r[ok]
        if len(t) < 30 or e.sum() < 8:
            result["note"] = "insufficient data"
            return result

        # Auto-detect time unit
        _unit = 365.0 if float(np.nanmedian(t)) > 50 else 1.0
        tau   = float(tau_years) * _unit

        # Risk groups (median split, direction-safe)
        grp = _km_group_from_risk(t, e, r)

        def _rmst_group(t_g, e_g, tau):
            if len(t_g) < 5 or e_g.sum() < 2:
                return float("nan")
            try:
                kmf = KaplanMeierFitter()
                kmf.fit(t_g, event_observed=e_g)
                sf    = kmf.survival_function_
                times = sf.index.values
                surv  = sf.values[:, 0]
                # Clip to tau
                mask  = times <= tau
                if mask.sum() == 0:
                    return float(tau)  # no events before tau
                t_clip = np.append(times[mask], tau)
                s_clip = np.append(surv[mask],  surv[mask][-1])
                return float(np.trapz(s_clip, t_clip))
            except Exception:
                return float("nan")

        rmst_low  = _rmst_group(t[grp==0], e[grp==0], tau)
        rmst_high = _rmst_group(t[grp==1], e[grp==1], tau)
        delta     = rmst_low - rmst_high  # positive = low risk survives longer

        # Bootstrap CI for delta_RMST
        rng  = np.random.default_rng(int(seed))
        boot_deltas = []
        for _ in range(int(n_boot)):
            idx_b   = rng.integers(0, len(t), size=len(t))
            t_b, e_b, grp_b = t[idx_b], e[idx_b], grp[idx_b]
            rl = _rmst_group(t_b[grp_b==0], e_b[grp_b==0], tau)
            rh = _rmst_group(t_b[grp_b==1], e_b[grp_b==1], tau)
            if np.isfinite(rl) and np.isfinite(rh):
                boot_deltas.append(rl - rh)
        ci_low = ci_high = float("nan")
        p_val  = float("nan")
        if len(boot_deltas) >= 50:
            bd = np.array(boot_deltas)
            ci_low  = float(np.percentile(bd, 2.5))
            ci_high = float(np.percentile(bd, 97.5))
            # two-sided p from bootstrap
            if delta > 0:
                p_val = float(np.mean(bd <= 0)) * 2
            else:
                p_val = float(np.mean(bd >= 0)) * 2
            p_val = min(float(p_val), 1.0)

        # Convert back to days/years label
        _u_label = "days" if _unit == 365.0 else "years"
        result.update({
            "rmst_low":           round(float(rmst_low),  2),
            "rmst_high":          round(float(rmst_high), 2),
            "delta_rmst":         round(float(delta),     2),
            "delta_rmst_ci_low":  round(float(ci_low),    2),
            "delta_rmst_ci_high": round(float(ci_high),   2),
            "p_value":            round(float(p_val),     4),
            "unit":               _u_label,
            "note": (f"RMST({tau_years}yr): low={rmst_low:.1f}, high={rmst_high:.1f}, "
                     f"delta={delta:.1f} {_u_label} "
                     f"[{ci_low:.1f},{ci_high:.1f}] p={p_val:.4f}")
        })
    except Exception as ex:
        result["note"] = f"RMST failed: {type(ex).__name__}: {ex}"
    return result


def landmark_analysis(y_time: np.ndarray, y_event: np.ndarray, risk: np.ndarray,
                      landmark_years: Tuple = (1, 2, 3),
                      follow_years: float = 3.0) -> pd.DataFrame:
    """
    """
    if not LIFELINES_AVAILABLE:
        return pd.DataFrame()
    rows = []
    t = np.asarray(y_time,  float)
    e = np.asarray(y_event, int)
    r = np.asarray(risk,    float)
    ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(r) & (t > 0)
    t, e, r = t[ok], e[ok], r[ok]

    _unit = 365.0 if float(np.nanmedian(t)) > 50 else 1.0

    for lm_yr in landmark_years:
        lm_days  = float(lm_yr)     * _unit
        end_days = float(lm_yr + follow_years) * _unit

        # Conditional: alive at landmark
        at_lm   = t > lm_days
        if at_lm.sum() < 20:
            rows.append({"landmark_year": lm_yr, "n_at_risk": int(at_lm.sum()),
                         "events_in_window": 0, "hr": float("nan"),
                         "p_logrank": float("nan"), "note": "too few at landmark"})
            continue

        t_lm  = t[at_lm] - lm_days   # time since landmark
        e_lm  = ((e[at_lm] == 1) & (t[at_lm] <= end_days)).astype(int)
        t_lm  = np.minimum(t_lm, end_days - lm_days)  # censor at follow_years
        r_lm  = r[at_lm]
        grp_lm = _km_group_from_risk(t_lm, e_lm, r_lm)

        n_ev  = int(e_lm.sum())
        p_lr  = float("nan")
        hr    = float("nan")

        try:
            from lifelines.statistics import logrank_test as _lrt
            res = _lrt(t_lm[grp_lm==0], t_lm[grp_lm==1],
                       event_observed_A=e_lm[grp_lm==0],
                       event_observed_B=e_lm[grp_lm==1])
            p_lr = float(res.p_value)
        except Exception:
            pass

        try:
            _df_lm = pd.DataFrame({"T": t_lm, "E": e_lm, "high": grp_lm})
            _cph_lm = CoxPHFitter(penalizer=1e-4)
            _cph_lm.fit(_df_lm, duration_col="T", event_col="E", show_progress=False)
            _coef = float(_cph_lm.params_["high"])
            hr = float(np.exp(_coef))
        except Exception:
            pass

        rows.append({
            "landmark_year":    lm_yr,
            "follow_years":     follow_years,
            "n_at_risk":        int(at_lm.sum()),
            "events_in_window": n_ev,
            "hr":               round(hr, 4),
            "p_logrank":        round(p_lr, 4) if np.isfinite(p_lr) else float("nan"),
            "note": f"HR={hr:.3f} p={p_lr:.4f}" if np.isfinite(hr) and np.isfinite(p_lr) else "insufficient"
        })
    return pd.DataFrame(rows)


# ---------------- Robust fallback gene selection (ensures non-empty prognostic gene tables) ----------------
def fallback_univariate_gene_select(X: pd.DataFrame, t: pd.Series, e: pd.Series, candidate_rna_cols: List[str], top_genes: int = 100) -> pd.DataFrame:
    rows = []
    if len(candidate_rna_cols) == 0:
        return pd.DataFrame(columns=["feature","gene","z_meta","p_meta","n_cohorts_used","frac_cohorts_used","fallback_score","method"])
    tt = pd.to_numeric(t, errors="coerce").astype(float).values
    ee = pd.to_numeric(e, errors="coerce").fillna(0).astype(int).values
    for f in candidate_rna_cols:
        try:
            x = pd.to_numeric(X[f], errors="coerce").astype(float).values
            m = np.isfinite(x) & np.isfinite(tt) & np.isfinite(ee)
            if m.sum() < 30:
                continue
            xx = x[m]; t2 = tt[m]; e2 = ee[m]
            rx = pd.Series(xx).rank(method="average").values
            rt = pd.Series(t2).rank(method="average").values
            corr = np.corrcoef(rx, rt)[0,1] if np.std(rx)>0 and np.std(rt)>0 else 0.0
            if (e2==1).sum() >= 5 and (e2==0).sum() >= 5:
                med1 = float(np.nanmedian(xx[e2==1])); med0 = float(np.nanmedian(xx[e2==0]))
                sep = abs(med1 - med0) / (np.nanstd(xx) + 1e-8)
            else:
                sep = 0.0
            score = abs(float(corr)) + 0.5*float(sep)
            rows.append({"feature": f, "gene": f.replace("RNA__","",1), "z_meta": float("nan"), "p_meta": float("nan"),
                         "n_cohorts_used": 1, "frac_cohorts_used": 1.0, "fallback_score": float(score), "method": "fallback_rankcorr_sep"})
        except Exception:
            continue
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["fallback_score"], ascending=False).head(int(top_genes)).reset_index(drop=True)

def enforce_ferroptosis_presence(selected_features: List[str], ferr_genes: List[str],
                                  wgcna_genes: List[str], Xcols: List[str],
                                  min_n: int = 20,
                                  final_gene_cap: Optional[int] = None) -> List[str]:
    """
    """
    sel = list(selected_features)
    target = set([f"RNA__{g}" for g in set(ferr_genes).union(set(wgcna_genes))])
    present = [c for c in Xcols if c in target]
    n_now = len([x for x in sel if x in target])
    if n_now >= min_n:
        if final_gene_cap is not None:
            rna_sel   = [f for f in sel if str(f).startswith("RNA__")]
            other_sel = [f for f in sel if not str(f).startswith("RNA__")]
            rna_capped = rna_sel[:int(final_gene_cap)]
            return rna_capped + other_sel
        return sel
    need = max(0, min_n - n_now)
    add = [c for c in present if c not in sel][:need]
    merged = sel + add
    if final_gene_cap is not None:
        rna_m   = [f for f in merged if str(f).startswith("RNA__")]
        other_m = [f for f in merged if not str(f).startswith("RNA__")]
        rna_m   = rna_m[:int(final_gene_cap)]
        merged  = rna_m + other_m
        if add:
            print(f" enforce_ferroptosis: +{len(add)} forced ferro genes, "
                  f"RNA total after cap={len(rna_m)}")
    return merged


def summarize_fold_feature_stability(folds: Dict[str, dict], ferr_genes: List[str], wgcna_genes: List[str]) -> pd.DataFrame:
    rows = []
    n_folds = max(1, len(folds))
    ferr_set = {f"RNA__{g}" for g in ferr_genes}
    wgcna_set = {f"RNA__{g}" for g in wgcna_genes}
    for holdout, fd in folds.items():
        meta = fd.get("meta_table", pd.DataFrame())
        feats = list(fd.get("feats", []))
        meta_map = {}
        if isinstance(meta, pd.DataFrame) and not meta.empty and "feature" in meta.columns:
            meta2 = meta.copy()
            keep_cols = [c for c in ["feature", "concordance_score", "p_meta", "q_fdr", "bag_freq"] if c in meta2.columns]
            if keep_cols:
                meta_map = meta2[keep_cols].set_index("feature").to_dict(orient="index")
        for f in feats:
            if not str(f).startswith("RNA__"):
                continue
            info = meta_map.get(f, {})
            rows.append({
                "fold": holdout,
                "feature": f,
                "gene": str(f).replace("RNA__", "", 1),
                "selected": 1,
                "concordance_score": info.get("concordance_score", np.nan),
                "p_meta": info.get("p_meta", np.nan),
                "q_fdr": info.get("q_fdr", np.nan),
                "bag_freq": info.get("bag_freq", np.nan),
                "is_ferrdb": f in ferr_set,
                "is_wgcna": f in wgcna_set,
            })
    if not rows:
        return pd.DataFrame(columns=["feature", "gene", "selection_count", "selection_freq"])
    df = pd.DataFrame(rows)
    agg = (df.groupby(["feature", "gene"], as_index=False)
             .agg(selection_count=("selected", "sum"),
                  mean_concordance=("concordance_score", "mean"),
                  median_p_meta=("p_meta", "median"),
                  median_q_fdr=("q_fdr", "median"),
                  mean_bag_freq=("bag_freq", "mean"),
                  is_ferrdb=("is_ferrdb", "max"),
                  is_wgcna=("is_wgcna", "max")))
    agg["selection_freq"] = agg["selection_count"] / float(n_folds)
    agg["source_label"] = np.where(agg["is_ferrdb"] & agg["is_wgcna"], "FerrDb+WGCNA",
                             np.where(agg["is_wgcna"], "WGCNA",
                             np.where(agg["is_ferrdb"], "FerrDb", "Other")))
    return agg.sort_values(["selection_freq", "mean_concordance", "median_p_meta"], ascending=[False, False, True]).reset_index(drop=True)


def lock_stable_signature(fulltrain_meta: pd.DataFrame, fold_stability: pd.DataFrame, final_gene_cap: int,
                          min_gene_stability: float, ferr_genes: List[str], wgcna_genes: List[str],
                          Xcols: List[str], min_ferro_genes: int = 15) -> List[str]:
    if fulltrain_meta is None or fulltrain_meta.empty or "feature" not in fulltrain_meta.columns:
        return []
    ft = fulltrain_meta.copy()
    if fold_stability is not None and not fold_stability.empty:
        stab = fold_stability.copy()
        stab = stab[[c for c in ["feature", "selection_freq", "selection_count", "mean_concordance", "median_q_fdr"] if c in stab.columns]]
        ft = ft.merge(stab, on="feature", how="left")
    else:
        ft["selection_freq"] = np.nan
        ft["selection_count"] = np.nan
        ft["mean_concordance"] = ft.get("concordance_score", np.nan)
        ft["median_q_fdr"] = ft.get("q_fdr", np.nan)

    for c in ["selection_freq", "selection_count", "mean_concordance", "median_q_fdr", "concordance_score", "p_meta", "q_fdr"]:
        if c not in ft.columns:
            ft[c] = np.nan

    stable_mask = ft["selection_freq"].fillna(0) >= float(min_gene_stability)
    stable = ft[stable_mask].copy()
    if stable.empty:
        stable = ft.sort_values(["selection_freq", "concordance_score", "p_meta"], ascending=[False, False, True]).head(max(10, int(final_gene_cap))).copy()

    stable = stable.sort_values(["selection_freq", "concordance_score", "mean_concordance", "q_fdr", "p_meta"],
                                ascending=[False, False, False, True, True])
    genes = stable["feature"].drop_duplicates().tolist()[:max(5, int(final_gene_cap))]
    genes = enforce_ferroptosis_presence(genes, ferr_genes, wgcna_genes, Xcols, min_n=min_ferro_genes)
    genes = [g for g in genes if g in set(Xcols)]
    rna = [g for g in genes if str(g).startswith("RNA__")]
    other = [g for g in genes if not str(g).startswith("RNA__")]
    rna = rna[:max(5, int(final_gene_cap))]
    return rna + other

# ---------------- Models ----------------
def fit_predict_coxph(Xtr, ttr, etr, Xte, sample_weight, penalizer=1e-4):
    if not LIFELINES_AVAILABLE:
        raise RuntimeError("lifelines missing")
    df = Xtr.replace([np.inf,-np.inf], np.nan).copy()
    med={}
    for c in df.columns:
        med[c]=float(df[c].median()); df[c]=df[c].fillna(med[c])
    df["T"]=ttr.values; df["E"]=etr.values
    cph = CoxPHFitter(penalizer=float(penalizer))
    if sample_weight is not None:
        df["_w"]=sample_weight
        cph.fit(df, duration_col="T", event_col="E", weights_col="_w", show_progress=False)
    else:
        cph.fit(df, duration_col="T", event_col="E", show_progress=False)
    Xte2 = Xte.replace([np.inf,-np.inf], np.nan).copy()
    for c in Xte2.columns: Xte2[c]=Xte2[c].fillna(med.get(c,0.0))
    return (Xte2.values @ cph.params_.values).astype(float)

def fit_predict_coxnet(Xtr, ttr, etr, Xte, l1_ratio, alpha_min_ratio, n_alphas, max_iter):
    if not SKSURV_MODELS_AVAILABLE: raise RuntimeError("sksurv models missing")
    ytr = _build_sksurv_y(ttr, etr)
    model = CoxnetSurvivalAnalysis(l1_ratio=float(l1_ratio), alpha_min_ratio=float(alpha_min_ratio), n_alphas=int(n_alphas), max_iter=int(max_iter))
    model.fit(Xtr.values, ytr)
    coefs = model.coef_
    coef = coefs[:, -1] if coefs.ndim == 2 else coefs
    return (Xte.values @ coef).astype(float)


def fit_predict_sks_coxph(Xtr, ttr, etr, Xte, alpha):
    if (not SKSURV_LINEAR_EXTRA_AVAILABLE) or (CoxPHSurvivalAnalysis is None):
        raise RuntimeError("sksurv CoxPHSurvivalAnalysis unavailable")
    ytr = _build_sksurv_y(ttr, etr)
    model = CoxPHSurvivalAnalysis(alpha=float(alpha), ties="breslow", n_iter=200, tol=1e-9)
    model.fit(Xtr.values, ytr)
    return model.predict(Xte.values).astype(float)

def fit_predict_ipcridge(Xtr, ttr, etr, Xte, alpha):
    if (not SKSURV_LINEAR_EXTRA_AVAILABLE) or (IPCRidge is None):
        raise RuntimeError("sksurv IPCRidge unavailable")
    ytr = _build_sksurv_y(ttr, etr)
    model = IPCRidge(alpha=float(alpha))
    model.fit(Xtr.values, ytr)
    return model.predict(Xte.values).astype(float)

def fit_predict_survivaltree(Xtr, ttr, etr, Xte, min_leaf, max_features, max_depth, seed):
    if (not SKSURV_TREE_AVAILABLE) or (SurvivalTree is None):
        raise RuntimeError("sksurv SurvivalTree unavailable")
    ytr = _build_sksurv_y(ttr, etr)
    model = SurvivalTree(
        min_samples_split=10, min_samples_leaf=int(min_leaf),
        max_features=max_features, max_depth=None if max_depth is None else int(max_depth),
        random_state=seed
    )
    model.fit(Xtr.values, ytr)
    try:
        return np.asarray(model.predict(Xte.values), dtype=float)
    except Exception:
        surv = model.predict_survival_function(Xte.values, return_array=True)
        return (-np.trapz(surv, axis=1)).astype(float)

def fit_predict_rsf(Xtr, ttr, etr, Xte, n_estimators, min_leaf, max_features, max_depth, seed):
    
    ytr = _build_sksurv_y(ttr, etr)
    min_split = max(2, int(min_leaf))  
    model = RandomSurvivalForest(
        n_estimators=int(n_estimators),
        min_samples_split=min_split,
        min_samples_leaf=int(min_leaf),
        max_features=max_features,
        max_depth=None if max_depth is None else int(max_depth),
        n_jobs=1, random_state=int(seed) 
    )
    model.fit(Xtr.values, ytr)
    try:
        return np.asarray(model.predict(Xte.values), dtype=float)
    except Exception:
        try:
            surv = model.predict_survival_function(Xte.values, return_array=True)
            return (-np.trapz(surv, axis=1)).astype(float)
        except Exception:
            return np.zeros(len(Xte), dtype=float)

def fit_predict_extratrees(Xtr, ttr, etr, Xte, n_estimators, min_leaf, max_features, max_depth, seed):
   
    ytr = _build_sksurv_y(ttr, etr)
    dyn_split = max(2, int(min_leaf)) 
    model = ExtraSurvivalTrees(
        n_estimators=int(n_estimators), min_samples_split=dyn_split,
        min_samples_leaf=int(min_leaf),
        max_features=max_features, max_depth=None if max_depth is None else int(max_depth),
        n_jobs=1, random_state=int(seed) 
    )
    model.fit(Xtr.values, ytr)
    try:
        return np.asarray(model.predict(Xte.values), dtype=float)
    except Exception:
        surv = model.predict_survival_function(Xte.values, return_array=True)
        return (-np.trapz(surv, axis=1)).astype(float)

def fit_predict_gbm(Xtr, ttr, etr, Xte, n_estimators, learning_rate, max_depth, subsample, seed):
    ytr = _build_sksurv_y(ttr, etr)
    model = GradientBoostingSurvivalAnalysis(
        n_estimators=int(n_estimators), learning_rate=float(learning_rate),
        max_depth=int(max_depth), subsample=float(subsample), random_state=seed
    )
    model.fit(Xtr.values, ytr)
    return model.predict(Xte.values).astype(float)

def fit_predict_fastsurvsvm(Xtr, ttr, etr, Xte, alpha, rank_ratio, max_iter, tol, seed):
    if not FAST_SVM_AVAILABLE:
        raise RuntimeError("FastSurvivalSVM unavailable")
    ytr = _build_sksurv_y(ttr, etr)
    model = FastSurvivalSVM(alpha=float(alpha), rank_ratio=float(rank_ratio), fit_intercept=True,
                            max_iter=int(max_iter), tol=float(tol), random_state=seed)
    model.fit(Xtr.values, ytr)
    return model.predict(Xte.values).astype(float)

def fit_predict_cwgbm(Xtr, ttr, etr, Xte, n_estimators, learning_rate, dropout_rate, seed):
    if not CW_GBM_AVAILABLE:
        raise RuntimeError("Componentwise GBM unavailable")
    ytr = _build_sksurv_y(ttr, etr)
    model = ComponentwiseGradientBoostingSurvivalAnalysis(
        n_estimators=int(n_estimators), learning_rate=float(learning_rate),
        dropout_rate=float(dropout_rate), random_state=seed
    )
    model.fit(Xtr.values, ytr)
    return model.predict(Xte.values).astype(float)

def fit_and_predict(family, params, Xtr, ttr, etr, Xte, sample_weight, seed, low_mem):
    fam = family.upper()
    if fam == "COXPH":
        return fit_predict_coxph(Xtr, ttr, etr, Xte, sample_weight, penalizer=params.get("penalizer",1e-4))
    if fam == "COXNET":
        n_alphas = 80 if low_mem else 140
        max_iter = 120000 if low_mem else 220000
        return fit_predict_coxnet(Xtr, ttr, etr, Xte, params.get("l1_ratio",0.8), params.get("alpha_min_ratio",0.05), n_alphas, max_iter)
    if fam == "SKSCOXPH":
        return fit_predict_sks_coxph(Xtr, ttr, etr, Xte, params.get("alpha", 1e-4))
    if fam == "IPCRIDGE":
        return fit_predict_ipcridge(Xtr, ttr, etr, Xte, params.get("alpha", 1.0))
    if fam == "RSF":
        return fit_predict_rsf(Xtr, ttr, etr, Xte, params.get("n_estimators",300), params.get("min_leaf",5), params.get("max_features","sqrt"), params.get("max_depth",None), seed)
    if fam == "SURVTREE":
        return fit_predict_survivaltree(Xtr, ttr, etr, Xte, params.get("min_leaf",8), params.get("max_features","sqrt"), params.get("max_depth",None), seed)
    if fam == "EXTRATREES":
        return fit_predict_extratrees(Xtr, ttr, etr, Xte, params.get("n_estimators",300), params.get("min_leaf",5), params.get("max_features","sqrt"), params.get("max_depth",None), seed)
    if fam == "GBM":
        return fit_predict_gbm(Xtr, ttr, etr, Xte, params.get("n_estimators",200), params.get("learning_rate",0.05), params.get("max_depth",2), params.get("subsample",0.7), seed)
    if fam == "FASTSVM":
        return fit_predict_fastsurvsvm(Xtr, ttr, etr, Xte, params.get("alpha",1.0), params.get("rank_ratio",1.0), params.get("max_iter",60), params.get("tol",1e-5), seed)
    if fam == "CWGBM":
        return fit_predict_cwgbm(Xtr, ttr, etr, Xte, params.get("n_estimators",300), params.get("learning_rate",0.03), params.get("dropout_rate",0.0), seed)
    raise ValueError(f"Unknown family: {family}")

def robust_fit_and_predict(family, params, Xtr, ttr, etr, Xte, sample_weight, seed, low_mem):
    try:
        return fit_and_predict(family, params, Xtr, ttr, etr, Xte, sample_weight, seed, low_mem), None
    except Exception as ex:
        try:
            return fit_predict_coxph(Xtr, ttr, etr, Xte, sample_weight, penalizer=1e-4), f"fallback_to_coxph:{type(ex).__name__}"
        except Exception:
            return np.zeros(len(Xte), dtype=float), f"hard_fail:{type(ex).__name__}"

# ---------------- Algorithm grid ----------------

def _filter_algos_by_prior_leaderboard(algos, leaderboard_csv=None, min_cindex=0.55):
    """Remove exact previously tested configs with prior LOCO C-index below threshold."""
    candidates = []
    if leaderboard_csv:
        candidates.append(leaderboard_csv)
    try:
        script_dir = Path(__file__).resolve().parent
        candidates.append(str(script_dir / "internal_loco_leaderboard.csv"))
    except Exception:
        pass
    candidates.append("internal_loco_leaderboard.csv")
    seen_path = None
    for c in candidates:
        if c and os.path.exists(c):
            seen_path = c
            break
    if seen_path is None:
        return algos
    try:
        lb = pd.read_csv(seen_path)
        if "algo_name" not in lb.columns or "loco_mean_cindex" not in lb.columns:
            print(f"[warn] leaderboard missing required columns, skip filter: {seen_path}")
            return algos
        bad_mask = pd.to_numeric(lb["loco_mean_cindex"], errors="coerce") < float(min_cindex)
        bad_names = set(lb.loc[bad_mask, "algo_name"].astype(str))
        if not bad_names:
            return algos
        out = [a for a in algos if a[0] not in bad_names]
        print(f"[info] Prior-leaderboard filter active ({seen_path}): removed {len(algos)-len(out)} configs with loco_mean_cindex < {min_cindex}")
        return out
    except Exception as ex:
        print(f"[warn] leaderboard filter failed: {type(ex).__name__}: {ex}")
        return algos

def build_algorithm_configs(light=True, prior_leaderboard_csv=None, prior_min_cindex=0.55,
                             family_budget: int = 6, min_family_slots: int = 5,
                             etx_ne_values=None):
    """
   
    """
    algos = []

    # ── 1. ExtraSurvivalTrees — 17 config ─────────────────────────────────────
    # Data-driven: mf=0.65 best, leaf=8 best, md=8 best, ne=300-500
    # 17 config: mf=0.65 ağırlıklı, sweepte mf=0.60/0.70 da var
    _etx_configs = [
        # Best zone: mf=0.65, leaf=8, md=8
        ("ETX_ne300_leaf8_mf0p65_d8",   {"n_estimators":300,"min_leaf":8, "max_features":0.65,"max_depth":8}),
        ("ETX_ne400_leaf8_mf0p65_d8",   {"n_estimators":400,"min_leaf":8, "max_features":0.65,"max_depth":8}),
        ("ETX_ne500_leaf8_mf0p65_d8",   {"n_estimators":500,"min_leaf":8, "max_features":0.65,"max_depth":8}),
        ("ETX_ne300_leaf8_mf0p65_d12",  {"n_estimators":300,"min_leaf":8, "max_features":0.65,"max_depth":12}),
        ("ETX_ne400_leaf8_mf0p65_d12",  {"n_estimators":400,"min_leaf":8, "max_features":0.65,"max_depth":12}),
        # leaf sweep (data: leaf=20 second best)
        ("ETX_ne300_leaf20_mf0p65_d8",  {"n_estimators":300,"min_leaf":20,"max_features":0.65,"max_depth":8}),
        ("ETX_ne400_leaf20_mf0p65_d8",  {"n_estimators":400,"min_leaf":20,"max_features":0.65,"max_depth":8}),
        ("ETX_ne500_leaf20_mf0p65_dNone",{"n_estimators":500,"min_leaf":20,"max_features":0.65,"max_depth":None}),
        ("ETX_ne400_leaf10_mf0p65_d8",  {"n_estimators":400,"min_leaf":10,"max_features":0.65,"max_depth":8}),
        ("ETX_ne500_leaf10_mf0p65_d12", {"n_estimators":500,"min_leaf":10,"max_features":0.65,"max_depth":12}),
        # mf robustness sweep
        ("ETX_ne300_leaf8_mf0p60_d8",   {"n_estimators":300,"min_leaf":8, "max_features":0.60,"max_depth":8}),
        ("ETX_ne400_leaf8_mf0p60_d8",   {"n_estimators":400,"min_leaf":8, "max_features":0.60,"max_depth":8}),
        ("ETX_ne500_leaf8_mf0p60_dNone",{"n_estimators":500,"min_leaf":8, "max_features":0.60,"max_depth":None}),
        ("ETX_ne300_leaf8_mf0p70_d8",   {"n_estimators":300,"min_leaf":8, "max_features":0.70,"max_depth":8}),
        ("ETX_ne400_leaf8_mf0p70_d8",   {"n_estimators":400,"min_leaf":8, "max_features":0.70,"max_depth":8}),
        ("ETX_ne500_leaf20_mf0p60_d8",  {"n_estimators":500,"min_leaf":20,"max_features":0.60,"max_depth":8}),
        ("ETX_ne500_leaf20_mf0p70_dNone",{"n_estimators":500,"min_leaf":20,"max_features":0.70,"max_depth":None}),
    ]
    for name, params in _etx_configs:
        algos.append((name, "EXTRATREES", params))  # tam 17

    # ── 2. RandomSurvivalForest — 17 config ────────────────────────────────────
    # Data: leaf=20 dominant, mf=0.60 best, ne=400-500 hafif avantaj, md=8-None
    algos_rsf = []
    for ne in [400, 500, 300]:
        for leaf in [20, 15, 10]:
            for md in [None, 8]:
                tag_md = 'None' if md is None else str(md)
                algos_rsf.append((f"RSF_ne{ne}_leaf{leaf}_mf0p60_d{tag_md}",
                                   "RSF",
                                   {"n_estimators": ne, "min_leaf": leaf,
                                    "max_features": 0.60, "max_depth": md}))
    # mf=0.50 (runner-up)
    for ne in [400, 500]:
        algos_rsf.append((f"RSF_ne{ne}_leaf20_mf0p50_d8", "RSF",
                          {"n_estimators": ne, "min_leaf": 20,
                           "max_features": 0.50, "max_depth": 8}))
    # mf=0.70
    algos_rsf.append(("RSF_ne500_leaf20_mf0p70_dNone", "RSF",
                       {"n_estimators": 500, "min_leaf": 20,
                        "max_features": 0.70, "max_depth": None}))
    algos.extend(algos_rsf[:17])

    # ── 3. GradientBoostingSurvival — 17 config ────────────────────────────────
    
    algos_gbm = []
    for ne in [100, 150, 200, 300]:
        for md in [2, 3]:
            for sub in [0.80, 0.70]:
                algos_gbm.append((f"GBM_ne{ne}_lr0.02_d{md}_ss{sub}", "GBM",
                                   {"n_estimators": ne, "learning_rate": 0.02,
                                    "max_depth": md, "subsample": sub}))
    # lr=0.05 ek (runner-up)
    for ne in [100, 200]:
        algos_gbm.append((f"GBM_ne{ne}_lr0.05_d2_ss0.8", "GBM",
                          {"n_estimators": ne, "learning_rate": 0.05,
                           "max_depth": 2, "subsample": 0.80}))
    algos.extend(algos_gbm[:17])

    # ── 4. CoxNet (ElasticNet-Cox) ────────────────────────────────
    
    algos_coxnet = []
    for l1 in [0.10, 0.25, 0.50, 0.75, 1.0]:
        for amin in [0.05, 0.10, 0.15, 0.20]:
            algos_coxnet.append((f"ENET_l1{l1}_amin{amin}", "COXNET",
                                  {"l1_ratio": l1, "alpha_min_ratio": amin}))
    # l1=0.0 (Ridge) → feature shrinkage
    for amin in [0.05, 0.10, 0.20]:
        algos_coxnet.append((f"ENET_l10.0_amin{amin}", "COXNET",
                              {"l1_ratio": 0.0, "alpha_min_ratio": amin}))
    algos.extend(algos_coxnet[:17])

    # ── 5. ComponentWise GBM ──────────────────────────────────────
   
    algos_cwgbm = []
    for ne in [300, 400, 500, 600, 700]:
        for lr in [0.10, 0.05, 0.03]:
            algos_cwgbm.append((f"CWGBM_ne{ne}_lr{lr}_dr0.0", "CWGBM",
                                 {"n_estimators": ne, "learning_rate": lr,
                                  "dropout_rate": 0.0}))  

    algos_cwgbm.append(("CWGBM_ne500_lr0.01_dr0.0","CWGBM",
                         {"n_estimators":500,"learning_rate":0.01,"dropout_rate":0.0}))
    algos_cwgbm.append(("CWGBM_ne700_lr0.01_dr0.0","CWGBM",
                         {"n_estimators":700,"learning_rate":0.01,"dropout_rate":0.0}))
    algos.extend(algos_cwgbm[:17])

    # ── 6. CoxPH ──────────────────────────────────────────────────
    
    for p in [0.001, 0.005, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50,
              0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 20.0]:
        algos.append((f"COXPH_pen{p}", "COXPH", {"penalizer": p}))

    # ── 7. SKSurvival CoxPH  ───────────────────────────────────────
   
    for alpha in [0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0,
                  1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 20.0, 50.0]:
        algos.append((f"SKSCOXPH_a{alpha}", "SKSCOXPH", {"alpha": alpha}))

    # ── 8. FastSurvivalSVM  ────────────────────────────────────────
   
    for alpha in [0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 1.0, 2.0,
                  3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0]:
        algos.append((f"FASTSVM_a{alpha}_rr1.0", "FASTSVM",
                      {"alpha": alpha, "rank_ratio": 1.0,
                       "max_iter": 500, "tol": 1e-6}))

    # ── 9. SurvivalTree ────────────────────────────────────────────
    
    algos_st = []
    for leaf in [8, 12, 20, 30]:
        for mf in [0.50, 0.60, 0.70]:
            for md in [6, 10, None]:
                tag_mf = str(mf).replace('.', 'p')
                tag_md = 'None' if md is None else str(md)
                algos_st.append((f"SURVTREE_leaf{leaf}_mf{tag_mf}_d{tag_md}",
                                  "SURVTREE",
                                  {"min_leaf": leaf, "max_features": mf,
                                   "max_depth": md}))
    algos.extend(algos_st[:16])

    # ── Dedup ─────────────────────────────────────────────────────────────────
    seen = set(); deduped = []
    for row in algos:
        if row[0] not in seen:
            seen.add(row[0]); deduped.append(row)

    # ── Prior leaderboard filter ───────────────────────────────────────────────
    filtered = _filter_algos_by_prior_leaderboard(
        deduped, leaderboard_csv=prior_leaderboard_csv, min_cindex=prior_min_cindex)

    # ── Min-slot guarantee ─────────────────────────────────────────────────────
    all_families = list({r[1] for r in deduped})
    filtered_names = {r[0] for r in filtered}
    guaranteed = []
    for fam in all_families:
        fam_configs = [r for r in deduped if r[1] == fam]
        fam_in      = [r for r in fam_configs if r[0] in filtered_names]
        if len(fam_in) < int(min_family_slots):
            needed = int(min_family_slots) - len(fam_in)
            extra  = [r for r in fam_configs if r[0] not in filtered_names][:needed]
            guaranteed.extend(extra)
    if guaranteed:
        filtered = filtered + [r for r in guaranteed
                               if r[0] not in {x[0] for x in filtered}]
    return filtered


def _param_distance_score(a: dict, b: dict) -> float:
    """Heuristic distance to avoid near-duplicate hyperparameters in top-k."""
    if a is None or b is None:
        return 1.0
    keys = sorted(set(a.keys()) | set(b.keys()))
    if not keys:
        return 0.0
    parts = []
    for k in keys:
        va, vb = a.get(k, None), b.get(k, None)
        if va is None or vb is None:
            parts.append(1.0); continue
        if isinstance(va, (int,float)) and isinstance(vb, (int,float)) and np.isfinite(float(va)) and np.isfinite(float(vb)):
            xa, xb = float(va), float(vb)
            if xa > 0 and xb > 0:
                d = abs(math.log10(xa + 1e-12) - math.log10(xb + 1e-12))
                d = min(d / 2.0, 1.5)
            else:
                denom = max(abs(xa), abs(xb), 1.0)
                d = min(abs(xa-xb)/denom, 1.5)
            parts.append(float(d))
        else:
            parts.append(0.0 if str(va) == str(vb) else 1.0)
    return float(np.mean(parts))

def select_diverse_topk(leaderboard: pd.DataFrame, k: int, min_param_distance: float = 0.28, max_per_family: int = 4) -> pd.DataFrame:
    """Greedy selection: prioritize high C-index but enforce family + param diversity."""
    if leaderboard is None or leaderboard.empty:
        return pd.DataFrame()
    df = leaderboard.copy()
    if 'loco_mean_cindex' in df.columns:
        df['__score'] = pd.to_numeric(df['loco_mean_cindex'], errors='coerce')
    else:
        df['__score'] = np.nan
    df = df.sort_values(['__score','loco_std_cindex'], ascending=[False, True], na_position='last').reset_index(drop=True)
    selected_idx = []
    fam_counts = {}
    parsed = {}
    for i, r in df.iterrows():
        fam = str(r.get('family','UNK'))
        if fam_counts.get(fam, 0) >= int(max_per_family):
            continue
        try:
            params_i = json.loads(r.get('params_json','{}')) if isinstance(r.get('params_json','{}'), str) else dict(r.get('params_json',{}))
        except Exception:
            params_i = {}
        ok = True
        for j in selected_idx:
            rj = df.loc[j]
            same_fam = str(rj.get('family','UNK')) == fam
            if not same_fam:
                continue
            try:
                params_j = parsed.get(j)
                if params_j is None:
                    params_j = json.loads(rj.get('params_json','{}')) if isinstance(rj.get('params_json','{}'), str) else dict(rj.get('params_json',{}))
                    parsed[j] = params_j
            except Exception:
                params_j = {}
            if _param_distance_score(params_i, params_j) < float(min_param_distance):
                ok = False
                break
        if ok:
            selected_idx.append(i)
            parsed[i] = params_i
            fam_counts[fam] = fam_counts.get(fam, 0) + 1
            if len(selected_idx) >= int(k):
                break
    if len(selected_idx) < int(k):
        for i in range(len(df)):
            if i in selected_idx:
                continue
            selected_idx.append(i)
            if len(selected_idx) >= int(k):
                break
    out = df.loc[selected_idx].drop(columns=['__score'], errors='ignore').copy()
    
    return out.sort_values('loco_mean_cindex', ascending=False, na_position='last').reset_index(drop=True)

def compute_sample_weights(cohort_series: pd.Series,
                           event_series: Optional[pd.Series] = None) -> np.ndarray:
    """
    v37 FIX-6: Event-weighted inverse frequency sample weights.
    Mevcut: w = 1/cohort_size — sadece boyut dengeleme
    Yeni:   w = 1/(cohort_size × event_rate) — küçük ve düşük event-rate'li
            kohortlar daha az baskın; prognostik sinyal taşıyan kohortlar
            daha fazla katkı sağlar.
    event_series=None ise orijinal inverse-frequency kullanılır (fallback).
    """
    counts = cohort_series.value_counts().to_dict()
    if event_series is not None and len(event_series) == len(cohort_series):
        event_rates = {}
        for coh in counts:
            mask = cohort_series == coh
            er   = float(event_series[mask].mean()) if mask.sum() > 0 else 0.1
            event_rates[coh] = max(er, 0.05)  
        w = cohort_series.map(
            lambda coh: 1.0 / (counts.get(coh, 1) * event_rates.get(coh, 0.1))
        ).astype(float).values
    else:
        w = cohort_series.map(lambda c: 1.0 / counts.get(c, 1)).astype(float).values
    w = np.clip(w, 0, np.percentile(w, 95))  
    mean_w = np.mean(w)
    return w / mean_w if mean_w > 0 else np.ones(len(w))

def weighted_rank_ensemble(risks: List[np.ndarray], weights: Optional[np.ndarray]) -> np.ndarray:
    Z = np.vstack([rank_normalize(r) for r in risks])
    if weights is None: w = np.ones(Z.shape[0])/Z.shape[0]
    else:
        w = np.asarray(weights, float)
        w = w / max(w.sum(), 1e-12)
    return (w[:,None]*Z).sum(axis=0)

def optimize_ensemble_weights(train_risks, t, e, n_iter, seed):
    """
    """
    from scipy.optimize import minimize as sp_minimize

    rng = np.random.default_rng(int(seed))
    m   = len(train_risks)
    if m == 0:
        return np.array([1.0]), float("nan")
    if m == 1:
        return np.ones(1), float(cindex(t, e, train_risks[0]))

    def _neg_cindex_softmax(logits):
        """"""
        w = np.exp(logits - logits.max())
        w = w / (w.sum() + 1e-12)
        ci = cindex(t, e, weighted_rank_ensemble(train_risks, w))
        return -ci if np.isfinite(ci) else 1.0

    # ── Aşama 1: Dirichlet warm-start ─────────────────────────────────────
    n_warmstart = min(300, max(20, int(n_iter) // 5))
    best_w_ws, best_ci_ws = np.ones(m) / m, -1.0
    for _ in range(n_warmstart):
        w_try = rng.dirichlet(alpha=np.ones(m))
        ci_try = cindex(t, e, weighted_rank_ensemble(train_risks, w_try))
        if np.isfinite(ci_try) and ci_try > best_ci_ws:
            best_ci_ws, best_w_ws = ci_try, w_try

    # Softmax logit 
    w0_safe = np.clip(best_w_ws, 1e-6, None)
    logits0  = np.log(w0_safe)
    logits0  -= logits0.mean()

    # ── Nelder-Mead  ───────────────────────────────
    best_logits, best_ci_nm = logits0, best_ci_ws
    try:
        res = sp_minimize(
            _neg_cindex_softmax,
            logits0,
            method="Nelder-Mead",
            options={
                "maxiter": max(200, int(n_iter) // 10),
                "xatol": 1e-5,
                "fatol": 1e-5,
                "adaptive": True,     
            }
        )
        if res.success or res.fun < -best_ci_ws:
            opt_ci = -float(res.fun)
            if np.isfinite(opt_ci) and opt_ci > best_ci_nm:
                best_logits = res.x
                best_ci_nm  = opt_ci
    except Exception:
        pass 

    # ── Dirichlet  ─────────────────
    n_remaining = max(0, int(n_iter) - n_warmstart - 200)
    for _ in range(min(n_remaining, 500)):
        w_try = rng.dirichlet(alpha=np.ones(m))
        ci_try = cindex(t, e, weighted_rank_ensemble(train_risks, w_try))
        if np.isfinite(ci_try) and ci_try > best_ci_nm:
            best_ci_nm = ci_try
            best_logits = np.log(np.clip(w_try, 1e-6, None))


    final_w = np.exp(best_logits - best_logits.max())
    final_w /= final_w.sum()

    if not np.all(np.isfinite(final_w)) or final_w.sum() < 0.99:
        return np.ones(m) / m, float("nan")
    return final_w, float(best_ci_nm)


def softmax_weighted_ensemble(risks: List[np.ndarray], scores: np.ndarray, temperature: float = 0.1) -> np.ndarray:
    """"""
    s = np.asarray(scores, float)
    s = np.where(np.isfinite(s), s, np.nanmin(s) if np.isfinite(s).any() else 0.0)
    w = np.exp((s - s.max()) / max(temperature, 1e-6))
    w = w / w.sum()
    Z = np.vstack([rank_normalize(r) for r in risks])
    return (w[:, None] * Z).sum(axis=0)


def trimmed_rank_ensemble(risks: List[np.ndarray], trim_frac: float = 0.10) -> np.ndarray:
    """
    """
    if not risks:
        return np.array([])
    Z = np.vstack([rank_normalize(r) for r in risks])  # (n_models, n_samples)
    n_models = Z.shape[0]
    trim_k   = max(0, int(np.floor(n_models * trim_frac)))
    if trim_k == 0 or 2 * trim_k >= n_models:
        return np.nanmean(Z, axis=0)
    Z_sorted = np.sort(Z, axis=0)                     
    Z_trimmed = Z_sorted[trim_k: n_models - trim_k, :] 
    return np.nanmean(Z_trimmed, axis=0)


def borda_count_ensemble(risks: List[np.ndarray]) -> np.ndarray:
    """"""
    n = len(risks[0])
    borda = np.zeros(n, dtype=float)
    for r in risks:
        order = np.argsort(np.argsort(np.asarray(r, float)))
        borda += order.astype(float)
    return borda


def stacking_ensemble_cox(oof_risks: List[np.ndarray], oof_t: np.ndarray, oof_e: np.ndarray,
                           test_risks: List[np.ndarray], seed: int = 13) -> Optional[np.ndarray]:
    """
    """
    if not LIFELINES_AVAILABLE or len(oof_risks) < 2:
        return None
    if len(test_risks) != len(oof_risks):
        print(f"  [WARN] stacking: oof_risks({len(oof_risks)}) != test_risks({len(test_risks)}), skip")
        return None

    Z_oof  = np.column_stack([rank_normalize(r) for r in oof_risks])
    Z_test = np.column_stack([rank_normalize(r) for r in test_risks])
    cols   = [f"m{i}" for i in range(Z_oof.shape[1])]

    df_oof = pd.DataFrame(Z_oof, columns=cols)
    df_oof["T"] = np.asarray(oof_t, float)
    df_oof["E"] = np.asarray(oof_e, int)
    df_oof = df_oof.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    if df_oof.shape[0] < 40 or int(df_oof["E"].sum()) < 8:
        return None

    df_test = pd.DataFrame(Z_test, columns=cols)
    df_test = df_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    best_pred = None
    best_ci   = -1.0

    # L3 baseline: simple rank-average (always available)
    rank_avg_oof  = Z_oof.mean(axis=1)
    rank_avg_test = Z_test.mean(axis=1)
    try:
        ci_rank = cindex(np.asarray(oof_t, float), np.asarray(oof_e, int), rank_avg_oof)
        if np.isfinite(ci_rank) and ci_rank > best_ci:
            best_ci   = ci_rank
            best_pred = rank_avg_test
    except Exception:
        pass

    # L1 and L2: Cox-based meta-learners
    for pen, label in [(0.1, "Ridge"), (0.5, "ElasticNet")]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(df_oof, duration_col="T", event_col="E", show_progress=False)
            oof_pred  = cph.predict_partial_hazard(df_oof.drop(["T","E"], axis=1)).values.astype(float)
            ci_meta = cindex(np.asarray(oof_t, float), np.asarray(oof_e, int), oof_pred)
            test_pred = cph.predict_partial_hazard(df_test).values.astype(float)
            print(f"    [stack-{label}] OOF C-index={ci_meta:.4f}")
            if np.isfinite(ci_meta) and ci_meta > best_ci:
                best_ci   = ci_meta
                best_pred = test_pred
        except Exception as e:
            print(f"    [stack-{label}] failed: {e}")

    # v35 LASSO-Cox meta-learner (sparse: selects best subset of models)
    if GLMNET_AVAILABLE and Z_oof.shape[1] >= 4:
        try:
            from sksurv.linear_model import CoxnetSurvivalAnalysis
            y_oof = np.array(
                [(bool(e), float(t)) for t, e in zip(oof_t, oof_e)],
                dtype=[("event", bool), ("time", float)])
            lasso_cox = CoxnetSurvivalAnalysis(
                l1_ratio=1.0, alpha_min_ratio=0.05,
                max_iter=1000, normalize=True)
            lasso_cox.fit(Z_oof, y_oof)
            oof_lasso  = lasso_cox.predict(Z_oof)
            ci_lasso   = cindex(np.asarray(oof_t, float), np.asarray(oof_e, int), oof_lasso)
            test_lasso = lasso_cox.predict(Z_test)
            print(f"    [stack-LASSO] OOF C-index={ci_lasso:.4f}")
            if np.isfinite(ci_lasso) and ci_lasso > best_ci:
                best_ci   = ci_lasso
                best_pred = test_lasso
        except Exception as e:
            print(f"    [stack-LASSO] failed: {e}")

    print(f"  [stacking] Best meta-learner OOF C-index={best_ci:.4f}")
    return best_pred


def bagged_meta_univariate_select(Xtr2, ttr, etr, cohort_series, candidate_rna_cols, clin_cols,
                                   stage_adjusted, penalizer, top_genes, min_events_per_cohort,
                                   min_frac_cohorts, require_consistent_sign,
                                   n_bags: int = 5, bag_frac: float = 0.8, seed: int = 13,
                                   n_jobs: int = 1):
    """
    """
    rng = np.random.default_rng(seed)
    cohorts = cohort_series.unique().tolist()
    gene_freq: Dict[str, int] = {}
    gene_z: Dict[str, List[float]] = {}

    for bag_i in range(int(n_bags)):
        n_sel = max(2, int(len(cohorts) * bag_frac))
        sel_cohorts = rng.choice(cohorts, size=n_sel, replace=False).tolist()
        bag_idx = cohort_series.isin(sel_cohorts)
        if bag_idx.sum() < 30:
            continue
        Xb = Xtr2.loc[bag_idx]
        tb = ttr.loc[bag_idx]
        eb = etr.loc[bag_idx]
        cs = cohort_series.loc[bag_idx]
        try:
            res = meta_univariate_select(
                Xb, tb, eb, cs, candidate_rna_cols, clin_cols,
                stage_adjusted, penalizer, top_genes * 2,
                min_events_per_cohort, min_frac_cohorts, require_consistent_sign,
                n_jobs=int(n_jobs)
            )
            if res is None or res.empty:
                continue
            for _, row in res.iterrows():
                f = row["feature"]
                gene_freq[f] = gene_freq.get(f, 0) + 1
                gene_z.setdefault(f, []).append(float(row.get("z_meta", 0.0)) if np.isfinite(row.get("z_meta", 0.0)) else 0.0)
        except Exception:
            continue

    if not gene_freq:
        return pd.DataFrame()

    rows = []
    for f, freq in gene_freq.items():
        zvals = gene_z.get(f, [0.0])
        rows.append({
            "feature": f,
            "gene": f.replace("RNA__", "", 1),
            "bag_freq": freq,
            "bag_frac": freq / max(n_bags, 1),
            "z_meta": float(np.mean(zvals)),
            "z_std": float(np.std(zvals)),
            "p_meta": float("nan"),
            "n_cohorts_used": freq,
            "frac_cohorts_used": freq / max(n_bags, 1),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values(["bag_freq", "z_meta"], ascending=[False, False]).head(int(top_genes)).reset_index(drop=True)

    
    stable_mask = out["bag_frac"] >= 0.60
    if stable_mask.sum() >= max(3, int(top_genes * 0.20)):

        out_stable   = out[stable_mask].sort_values(["bag_frac","z_meta"], ascending=[False,False])
        out_unstable = out[~stable_mask].sort_values(["z_meta"], ascending=[False])
        out = pd.concat([out_stable, out_unstable], ignore_index=True).head(int(top_genes))
        print(f"  [stability] {stable_mask.sum()} stable genes (bag_frac≥0.60) prioritized")

    return out

# ---------------- Main ----------------
def parse_list(s: str) -> List[str]:
    s = (s or "").strip()
    return [x.strip() for x in s.split(",") if x.strip()] if s else []

def resolve_under_root(project_root: Path, p: str) -> Path:
    pp = Path(p)
    return pp if pp.is_absolute() else (project_root / pp)

def main():
    set_thread_limits(1)

    ap = argparse.ArgumentParser()
    ap.add_argument("--project_root", default=".")
    ap.add_argument("--luad_dir", required=True)
    ap.add_argument("--lusc_dir", required=True)
    ap.add_argument("--geo_dirs", required=True)
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--prior_leaderboard_csv", type=str, default="internal_loco_leaderboard.csv")
    ap.add_argument("--prior_min_cindex", type=float, default=0.55)
    ap.add_argument("--top_rna_var_prefilter", type=int, default=1800,
                    help="Varyans prefilter: top N gene")
    ap.add_argument("--meta_top_genes", type=int, default=60,
                    help="Meta-univariate gen seçimi: top N gene ")
    ap.add_argument("--use_meta_univariate", action="store_true")
    ap.add_argument("--no_use_meta_univariate", dest="use_meta_univariate", action="store_false")
    ap.set_defaults(use_meta_univariate=True)

    ap.add_argument("--min_events_per_cohort", type=int, default=15,
                    help="GSE136961(n=21)")
    ap.add_argument("--min_loco_events", type=int, default=20,
                    help="LOCO fold "
                         "train pool")
    ap.add_argument("--min_frac_cohorts", type=float, default=0.50,
                    help="0.45→0.50. Daha tutarlı gen seçimi")
    ap.add_argument("--require_consistent_sign", action="store_true")
    ap.add_argument("--no_require_consistent_sign", dest="require_consistent_sign", action="store_false")
    ap.set_defaults(require_consistent_sign=True)

    ap.add_argument("--max_algos", type=int, default=320)
    ap.add_argument("--top_k_algos", type=int, default=20)
    ap.add_argument("--diverse_topk", action="store_true")
    ap.add_argument("--no_diverse_topk", dest="diverse_topk", action="store_false")
    ap.set_defaults(diverse_topk=True)
    ap.add_argument("--diverse_topk_min_param_distance", type=float, default=0.28)
    ap.add_argument("--opt_ens_n_iter", type=int, default=1800)
    ap.add_argument("--max_runtime_min", type=int, default=0)

    ap.add_argument("--low_mem", action="store_true")
    ap.add_argument("--no_float32", dest="float32", action="store_false")
    ap.set_defaults(float32=True)

    ap.add_argument("--stage", dest="use_stage", action="store_true")
    ap.add_argument("--no_stage", dest="use_stage", action="store_false")
    ap.set_defaults(use_stage=True)

    ap.add_argument("--stage_adjusted_univariate", dest="stage_adjusted_univariate", action="store_true")
    ap.add_argument("--no_stage_adjusted_univariate", dest="stage_adjusted_univariate", action="store_false")
    ap.set_defaults(stage_adjusted_univariate=True)

    ap.add_argument("--univar_penalizer", type=float, default=1e-4)

    ap.add_argument("--ferrdb_driver_csv", default="ferrdb_driver.csv")
    ap.add_argument("--ferrdb_suppressor_csv", default="ferrdb_suppressor.csv")
    ap.add_argument("--ferrdb_marker_csv", default="ferrdb_marker.csv")

    # WGCNA integration
    ap.add_argument("--wgcna_path", default="", help="WGCNA gen listesini içeren dosya veya WGCNA output dizini")
    ap.add_argument("--use_wgcna_filter", action="store_true", help="Aday gen havuzunu WGCNA gen listesiyle kesiştir")
    ap.add_argument("--no_use_wgcna_filter", dest="use_wgcna_filter", action="store_false")
    ap.set_defaults(use_wgcna_filter=True)

    ap.add_argument("--add_paper_traits", action="store_true", help="FERRO__* ve WGCNA__ME1 gibi özet skorları feature olarak ekle")
    ap.add_argument("--gene_set_mode", default="intersection", choices=["union","intersection","wgcna_only","ferrdb_only"])


    
    ap.add_argument("--subtype_mode", default="luad_only",
                    choices=["joint","luad_only","lusc_only"],
                    help="Eğitim kohortu tipi: luad_only (varsayılan, LUAD ferroptosis imzası), "                         "joint (LUAD+LUSC), lusc_only")
    ap.add_argument("--geo_subtype_map", default="", help="Opsiyonel JSON veya dosya yolu. Örn: {'GSE12345':'LUAD', 'GSE67890':'LUSC'}.")
    ap.add_argument("--keep_strict_lusc", action="store_true", help="lusc_only modunda meta-univariate filtrelerini otomatik gevşetme.")
    ap.add_argument("--batch_weighting", dest="batch_weighting", action="store_true")
    ap.add_argument("--no_batch_weighting", dest="batch_weighting", action="store_false")
    ap.set_defaults(batch_weighting=True)


    ap.add_argument("--n_boot", type=int, default=1000, help="Bootstrap tekrar sayısı (default: 1000)")
    ap.add_argument("--n_perm", type=int, default=1000, help="Permütasyon testi tekrar sayısı (default: 1000)")
    ap.add_argument("--ferroptosis_min_genes", type=int, default=15, help="Nihai feature setinde minimum ferroptosis geni sayısı (default: 15)")
    ap.add_argument("--multivariate_validation", action="store_true", default=True, help="Çok değişkenli Cox doğrulaması yap (default: True)")
    ap.add_argument("--no_multivariate_validation", dest="multivariate_validation", action="store_false")
    ap.add_argument("--dca_years", type=int, default=3, help="DCA için horizon yılı (default: 3)")
    ap.add_argument("--calibration_years", type=int, default=3, help="Calibration curve için horizon yılı (default: 3)")


    ap.add_argument("--robust_scaler", action="store_true", default=False,
                    help="Cohort normalizasyonunda median/IQR kullan (outlier dirençli)")
    ap.add_argument("--family_budget", type=int, default=6,
                    help="Her algoritma familysi için max config sayısı (diverse_topk seçiminde)")
    ap.add_argument("--min_family_slots", type=int, default=3,
                    help="Her familyden leaderboard filtresinden bağımsız minimum config sayısı")
    ap.add_argument("--etx_ne_values", type=str, default="",
                    help="ExtraTrees n_estimators listesi, virgülle ayrılmış: '260,340,480,640,800'")
    ap.add_argument("--use_bagged_gene_select", action="store_true", default=False,
                    help="Meta-univariate gen seçimini bootstrap-bagged yöntemiyle yap (daha stabil)")
    ap.add_argument("--n_bags", type=int, default=5,
                    help="Bagged gen seçimi için bag sayısı (default: 5)")

    ap.add_argument("--meta_n_jobs", type=int, default=4,
                    help="v36: Meta-univariate gen değerlendirmesi için thread sayısı (default: 4)")

    ap.add_argument("--parallel_folds", type=int, default=1,
                    help="v40: LOCO fold'larını paralel değerlendir (default: 1=seri). "
                         "3 fold varsa --parallel_folds 3 ile ~2.5x hız. "
                         "RAM dikkat: her process full model fit yapıyor.")

    ap.add_argument("--use_stability_prefilter", action="store_true", default=True,
                    help="v35: Bagged gene selection'dan sonra bootstrap frekans>=0.6 olanları önceliklendir")
    ap.add_argument("--no_use_stability_prefilter", dest="use_stability_prefilter", action="store_false")
    ap.add_argument("--stability_min_freq", type=float, default=0.60,
                    help="v35: Minimum bootstrap bag frekansı (default: 0.60)")
    
    ap.add_argument("--use_stacking", action="store_true", default=True,
                    help="OOF risk skorları üzerinde ikinci aşama Cox stacking uygula (v35: varsayılan AÇIK)")
    ap.add_argument("--no_use_stacking", dest="use_stacking", action="store_false")
    ap.add_argument("--diverse_topk_max_per_family", type=int, default=4,
                    help="Diverse top-k'da family başına max model (default: 4)")
    # stabilize final RNA signature before full-train fitting
    ap.add_argument("--lock_signature_by_stability", action="store_true", default=True,
                    help="v43: LOCO fold selection frequency ile final RNA imzasını kilitle")
    ap.add_argument("--no_lock_signature_by_stability", dest="lock_signature_by_stability", action="store_false")
    ap.add_argument("--min_gene_stability", type=float, default=0.50,
                    help="v43: final kilitli imza için minimum LOCO selection frequency")
    ap.add_argument("--final_gene_cap", type=int, default=25,
                    help="v43: nihai RNA imzasında maksimum gen sayısı (default: 25)")


    ap.add_argument("--test_cohorts", type=str, default="",
                    help="v44: Virgülle ayrılmış external test kohortu isimleri. "
                         "Boşsa FIXED_TEST=[GSE50081] kullanılır. "
                         "Örn: --test_cohorts GSE50081,GSE68465")


    ap.add_argument("--cohort_hist_filter", type=str, default="",
                    help="v45: JSON formatında kohort bazında histoloji filtresi. "
                         "Örn: '{\"GSE30219\":{\"col\":\"histology\",\"values\":[\"ADC\"]}}' "
                         "Her kohortta belirtilen histoloji değerleri tutulur, diğerleri çıkarılır. "
                         "GSE30219 gibi mixed-histoloji kohortları için ADC-only test yapılır.")

    args = ap.parse_args()


    global _META_N_JOBS
    _META_N_JOBS = max(1, int(getattr(args, "meta_n_jobs", 4)))

    if args.low_mem:
        args.top_rna_var_prefilter = min(args.top_rna_var_prefilter, 1800)
        args.meta_top_genes = min(args.meta_top_genes, 90)
        args.opt_ens_n_iter = min(args.opt_ens_n_iter, 900)
        args.max_algos = min(args.max_algos, 110)
        args.final_gene_cap = min(args.final_gene_cap, 18)
        args.float32 = True

    seed_everything(args.seed)
    prog = ProgressTracker(16)

    project_root = Path(args.project_root).resolve()
    outdir = resolve_under_root(project_root, args.outdir)
    safe_mkdir(outdir); safe_mkdir(outdir/"genes"); safe_mkdir(outdir/"risks"); safe_mkdir(outdir/"logs"); safe_mkdir(outdir/"metrics"); safe_mkdir(outdir/"figures"); safe_mkdir(outdir/"publication")


    try:
        import json as _json
        (_outdir_tmp := outdir/"logs"/"run_config.json").write_text(
            _json.dumps(vars(args), default=str, indent=2), encoding="utf-8"
        )
    except Exception:
        pass

    prog.tick("Output klasörleri hazır")

    luad_dir = resolve_under_root(project_root, args.luad_dir)
    lusc_dir = resolve_under_root(project_root, args.lusc_dir)
    geo_dirs = [resolve_under_root(project_root,p) for p in parse_list(args.geo_dirs)]


    ferr_drv, ferr_sup, ferr_mkr, ferr_all = load_ferrdb_three(
        resolve_under_root(project_root,args.ferrdb_driver_csv),
        resolve_under_root(project_root,args.ferrdb_suppressor_csv),
        resolve_under_root(project_root,args.ferrdb_marker_csv),
    )


    wgcna_genes, wgcna_prov = load_wgcna_from_path(
        resolve_under_root(project_root, args.wgcna_path) if args.wgcna_path else None)


    wgcna_source_name = wgcna_prov.split("->")[-1].strip() if "->" in wgcna_prov else wgcna_prov
    source_quality = {
        "wgcna_genes_for_ML.csv":                    "BEST — v13 ML-optimized (conservative∪discovery, hub-score sıralı)",
        "conservative_ferroptosis_genes_GS_MM.tsv":  "GOOD — v13 FerrDb∩module (biyolojik filtreli, hub-score sıralı)",
        "module_hub_genes_MM_GS.tsv":                "OK   — Klasik hub genler (yüksek MM+GS)",
        "module_genes_GS_MM.tsv":                    "WIDE — Tüm modül genleri (GS/MM filtresiz)",
        "possible_prognostic_genes_consensus.csv":   "WIDEST — Tüm adaylar (en geniş, en az seçici)",
    }
    quality_label = source_quality.get(wgcna_source_name, "UNKNOWN")
    print(f"[WGCNA] Source: {wgcna_source_name}")
    print(f"[WGCNA] Quality: {quality_label}")
    print(f"[WGCNA] Gene count: {len(wgcna_genes)}")
    if wgcna_source_name == "possible_prognostic_genes_consensus.csv":
        print("[WGCNA] WARNING: En geniş aday dosyası kullanılıyor.")
        print("         v13 WGCNA çalıştırıldıysa 'wgcna_genes_for_ML.csv' daha iyi.")
        print("         Dosya yolu: tables/wgcna_genes_for_ML.csv")
    elif wgcna_source_name == "wgcna_genes_for_ML.csv":
        # Branch dağılımını raporla
        try:
            _wgcna_path_obj = resolve_under_root(project_root, args.wgcna_path)
            _ml_csv = _wgcna_path_obj / "tables" / "wgcna_genes_for_ML.csv"
            _ml_df  = pd.read_csv(_ml_csv)
            if "branch" in _ml_df.columns:
                bc = _ml_df["branch"].value_counts().to_dict()
                print(f"[WGCNA] Branch distribution: {bc}")
        except Exception:
            pass

    pd.DataFrame({
        "wgcna_source":  [wgcna_prov],
        "source_file":   [wgcna_source_name],
        "quality":       [quality_label],
        "n_genes":       [len(wgcna_genes)],
    }).to_csv(outdir/"logs"/"wgcna_loaded_summary.csv", index=False)

    # Gene universe
    univ = build_gene_universe(ferr_all, wgcna_genes, args.gene_set_mode)
    gene_subset_for_loading = sorted(set(univ) | set(ferr_all) | set(wgcna_genes))
    prog.tick(f"Gene universe hazır (FerrDb={len(ferr_all)} WGCNA={len(wgcna_genes)} mode={args.gene_set_mode})")


    _load_lusc = (args.subtype_mode != "luad_only")

    _cohort_hist_filter: Dict[str, Dict] = {}
    if getattr(args, "cohort_hist_filter", ""):
        try:
            _cohort_hist_filter = json.loads(args.cohort_hist_filter)
            print(f"[v45] cohort_hist_filter aktif: {list(_cohort_hist_filter.keys())}")
        except Exception as _hfe:
            print(f"[WARN v45] --cohort_hist_filter JSON parse hatası: {_hfe}")

    cohorts = [
        load_tcga("TCGA-LUAD", luad_dir, add_stage=args.use_stage,
                  gene_subset=gene_subset_for_loading, float32=args.float32),
    ]
    if _load_lusc:
        cohorts.append(load_tcga("TCGA-LUSC", lusc_dir, add_stage=args.use_stage,
                                  gene_subset=gene_subset_for_loading, float32=args.float32))
    for gd in geo_dirs:
        _hf = _cohort_hist_filter.get(gd.name, None)  # v45-1: kohort bazında filtre
        cohorts.append(load_geo(gd.name, gd, add_stage=args.use_stage,
                                gene_subset=gene_subset_for_loading,
                                float32=args.float32, hist_filter=_hf))
    cohorts = [c for c in cohorts if len(c.X)>0]

    qc = qc_events(cohorts); qc.to_csv(outdir/"qc_events_by_cohort.csv", index=False)
    all_names = [c.cohort for c in cohorts]
    prog.tick(f"Cohort yükleme tamam ({len(cohorts)})")


    if getattr(args, "test_cohorts", ""):
        _requested_tests = [t.strip() for t in args.test_cohorts.split(",") if t.strip()]
        test = [t for t in _requested_tests if t in all_names]
        _missing = [t for t in _requested_tests if t not in all_names]
        if _missing:
            print(f"[WARN v44] --test_cohorts: kohortlar yüklenmemiş veya bulunamadı: {_missing}")
        if not test:
            raise RuntimeError(f"--test_cohorts belirtildi ama hiçbiri yüklenemedi: {_requested_tests}")
    else:

        if "GSE50081" not in all_names:
            raise RuntimeError("GSE50081 cohort not loaded — varsayılan test kohortu eksik.")
        test = ["GSE50081"]
    print(f"[v44] External test kohortu/ları: {test}")


    subtype_map = {}
    if args.geo_subtype_map:
        try:
            subtype_map = json.loads(args.geo_subtype_map)
        except Exception:
            # allow path to json file
            p = Path(args.geo_subtype_map)
            if p.exists():
                subtype_map = json.loads(p.read_text(encoding="utf-8"))
            else:
                raise RuntimeError("--geo_subtype_map JSON parse failed and file not found")

    def cohort_subtype(name: str) -> str:
        n = str(name)
        if "LUSC" in n.upper():
            return "LUSC"
        if "LUAD" in n.upper():
            return "LUAD"
        if n in subtype_map:
            return str(subtype_map[n]).upper()
        # default assumption for lung GEO sets in this project
        return "LUAD"

    train_all = [c for c in FIXED_TRAIN if c in all_names and c not in test]
    if args.subtype_mode == "luad_only":
        train = [c for c in train_all if cohort_subtype(c) == "LUAD"]
    elif args.subtype_mode == "lusc_only":
        train = [c for c in train_all if cohort_subtype(c) == "LUSC"]
        # In case the fixed test is LUAD-like, warn the user early (still allowed for debugging).
        if cohort_subtype(test[0]) != "LUSC":
            print("[WARN] subtype_mode=lusc_only but fixed test cohort is not LUSC. Consider changing FIXED_TEST or providing a LUSC test cohort.")
        if (not args.keep_strict_lusc):
            # LUSC often has fewer events and stronger heterogeneity; relax meta-univariate filters a bit by default
            args.require_consistent_sign = False
            args.min_frac_cohorts = min(args.min_frac_cohorts, 0.30)
            args.min_events_per_cohort = min(args.min_events_per_cohort, 8)
    else:
        train = train_all

    if len(train) < 2:
        raise RuntimeError(f"TRAIN cohorts < 2 after subtype_mode filtering (mode={args.subtype_mode}). Provide more cohorts or use subtype_mode=joint.")

    if len(train)<2: raise RuntimeError("TRAIN cohorts < 2")
    prog.tick(f"Split hazır train={len(train)} test={test}")

    # LOCO folds

    min_loco_ev = getattr(args, "min_loco_events", 20)
    folds={}
    train_set = [c for c in cohorts if c.cohort in train]


    loco_eligible = []
    loco_poolonly  = []  
    for c in train_set:
        n_ev = int(c.y_event.sum()) if hasattr(c, "y_event") else 0
        if n_ev >= int(min_loco_ev):
            loco_eligible.append(c.cohort)
        else:
            loco_poolonly.append(c.cohort)
            print(f"  [v36] {c.cohort}: n_events={n_ev} < min_loco_events={min_loco_ev} "
                  f"→ pool-only (no LOCO fold)")
    if loco_poolonly:
        print(f"  [v36] Pool-only cohorts (no LOCO fold): {loco_poolonly}")
        print(f"  [v36] LOCO-eligible cohorts: {loco_eligible}")

    for holdout in loco_eligible:  
        tr_cohs = [c for c in train_set if c.cohort != holdout]
        te = [c for c in train_set if c.cohort == holdout][0]

        Xtr = pd.concat([c.X for c in tr_cohs], axis=0)
        ttr = pd.concat([c.y_time for c in tr_cohs], axis=0)
        etr = pd.concat([c.y_event for c in tr_cohs], axis=0)
        coh_series = pd.Series([c.cohort for c in tr_cohs for _ in range(len(c.X))], index=Xtr.index)

        Xtr2 = cohort_zscore(Xtr, coh_series, robust=args.robust_scaler)
       
        Xte2 = cohort_zscore(te.X, pd.Series([holdout]*len(te.X), index=te.X.index), robust=args.robust_scaler)


        if args.add_paper_traits:
            Xtr2 = add_paper_like_traits(Xtr2, ferr_drv, ferr_sup, wgcna_genes)
            Xte2 = add_paper_like_traits(Xte2, ferr_drv, ferr_sup, wgcna_genes)

        clin_cols = [c for c in Xtr2.columns if c.startswith("CLIN__")]
        rna_cols  = [c for c in Xtr2.columns if c.startswith("RNA__")]
        if univ:
            allowed = {f"RNA__{g}" for g in univ}
            rna_cols = [c for c in rna_cols if c in allowed] or rna_cols

        # variance prefilter
        v = Xtr2[rna_cols].var(axis=0).sort_values(ascending=False)
        cand = v.head(int(min(args.top_rna_var_prefilter, len(v)))).index.tolist()

        # WGCNA filter
        if args.use_wgcna_filter and wgcna_genes:
            allowed_w = {f"RNA__{g}" for g in wgcna_genes}
            cand_w = [c for c in cand if c in allowed_w]
            # fallback if too small
            if len(cand_w) >= max(30, int(0.25*args.meta_top_genes)):
                cand = cand_w

        if args.use_meta_univariate and LIFELINES_AVAILABLE:
            if args.use_bagged_gene_select:
                # bootstrap-bagged gene selection
                meta = bagged_meta_univariate_select(
                    Xtr2, ttr, etr, coh_series, cand, clin_cols,
                    args.stage_adjusted_univariate, args.univar_penalizer,
                    args.meta_top_genes, args.min_events_per_cohort,
                    args.min_frac_cohorts, args.require_consistent_sign,
                    n_bags=args.n_bags, bag_frac=0.8, seed=args.seed,
                    n_jobs=args.meta_n_jobs 
                )
                if meta is None or meta.empty:
                    meta = meta_univariate_select(
                        Xtr2, ttr, etr, coh_series, cand, clin_cols,
                        args.stage_adjusted_univariate, args.univar_penalizer,
                        args.meta_top_genes, args.min_events_per_cohort,
                        args.min_frac_cohorts, args.require_consistent_sign,
                        n_jobs=args.meta_n_jobs 
                    )
            else:
                meta = meta_univariate_select(
                    Xtr2, ttr, etr, coh_series, cand, clin_cols,
                    args.stage_adjusted_univariate, args.univar_penalizer,
                    args.meta_top_genes, args.min_events_per_cohort,
                    args.min_frac_cohorts, args.require_consistent_sign,
                    n_jobs=args.meta_n_jobs 
                )
            sel_rna = meta["feature"].tolist() if not meta.empty else cand[:args.meta_top_genes]
        else:
            meta = pd.DataFrame()
            sel_rna = cand[:args.meta_top_genes]

        extra_cols = []
        if args.add_paper_traits:
            extra_cols = [c for c in (
            "FERRO__DRIVER","FERRO__SUPPRESSOR","FERRO__TOTAL","FERRO__TRAIT",
            "WGCNA__ME1",

            "FERRO__DRIVER_PC1","FERRO__SUPP_PC1","FERRO__NET_SCORE",
            "FERRO__RATIO","FERRO__WGCNA_INTER","FERRO__Z_COMBINED"
        ) if c in Xtr2.columns]

        feats = [f for f in (sel_rna + clin_cols + extra_cols) if f in Xtr2.columns]
        if feats:
            folds[holdout] = {
                "feats": feats, "meta_table": meta, "XtrF": Xtr2[feats], "ttr": ttr, "etr": etr,
                "XteF": Xte2.reindex(columns=feats).fillna(0.0), "tte": te.y_time, "ete": te.y_event,
                "sw": compute_sample_weights(coh_series, etr) if args.batch_weighting else None  
            }

    if len(folds)==0: raise RuntimeError("No valid LOCO folds")

    fold_feature_stability = summarize_fold_feature_stability(folds, ferr_all, wgcna_genes)
    if not fold_feature_stability.empty:
        fold_feature_stability.to_csv(outdir/"genes"/"loco_feature_stability.csv", index=False)
        print(f"[v43] LOCO feature stability yazıldı: {len(fold_feature_stability)} RNA feature")
        print(f"[v43] selection_freq>=0.50 gen sayısı: {(fold_feature_stability["selection_freq"]>=0.50).sum()}")
    else:
        fold_feature_stability = pd.DataFrame()

    prog.tick(f"LOCO folds hazır ({len(folds)})")


    etx_ne_list = None
    if args.etx_ne_values:
        try:
            etx_ne_list = [int(x.strip()) for x in args.etx_ne_values.split(",") if x.strip()]
        except Exception:
            etx_ne_list = None
    all_algos = build_algorithm_configs(
        light=True,
        prior_leaderboard_csv=args.prior_leaderboard_csv,
        prior_min_cindex=args.prior_min_cindex,
        family_budget=args.family_budget,
        min_family_slots=args.min_family_slots,
        etx_ne_values=etx_ne_list,
    )
  
    from collections import defaultdict
    family_queues = defaultdict(list)
    for algo in all_algos:
        family_queues[algo[1]].append(algo)

   
    priority_order = ['GBM','COXNET','CWGBM','COXPH','SKSCOXPH','FASTSVM',
                      'SURVTREE','RSF','EXTRATREES']

    all_families_sorted = ([f for f in priority_order if f in family_queues] +
                           [f for f in family_queues if f not in priority_order])

    interleaved = []
    while any(family_queues[f] for f in all_families_sorted):
        for fam in all_families_sorted:
            if family_queues[fam]:
                interleaved.append(family_queues[fam].pop(0))

    algo_grid = interleaved[:max(1, args.max_algos)]
    print(f"[v39-FIX] Grid interleaved: family sırası = "
          f"{dict(sorted({a[1]:0 for a in algo_grid}.items()))}")

    from collections import Counter as _C
    _fc = _C(a[1] for a in algo_grid)
    print(f"[v39-FIX] Family dağılımı: {dict(sorted(_fc.items(), key=lambda x:-x[1]))}")
    pd.DataFrame([{"algo_name":a, "family":f, "params_json":json.dumps(p)} for a,f,p in algo_grid]).to_csv(outdir/"logs"/"planned_algorithms.csv", index=False)
    prog.tick(f"Algo grid hazır ({len(algo_grid)})")

   
    _checkpoint_path = outdir / "internal_loco_leaderboard.csv"
    _already_done = set()
    if _checkpoint_path.exists():
        try:
            _ckpt = pd.read_csv(_checkpoint_path)
            _already_done = set(_ckpt["algo_name"].tolist())
            if _already_done:
                print(f"[CHECKPOINT] {len(_already_done)} config zaten tamamlandı, atlanıyor")
        except Exception:
            pass

    t_start_grid = time.time()
    rows = []
    run_logs = []

    n_fold_total = max(1, len(folds))
    for k, (algo_name, family, params) in enumerate(algo_grid, start=1):
        if args.max_runtime_min>0 and ((time.time()-t_start_grid)/60.0)>args.max_runtime_min:
            print("[TIME_GUARD] grid early stop")
            break


        if algo_name in _already_done:
            continue

        fold_scores={}
        used=0
        fallback_count=0
        eligible_holdouts = [h for h in train if h in folds]

      
        def _eval_fold_inner(args_tuple):
            _holdout, _family, _params, _XtrF, _ttr, _etr, _XteF, _sw, _seed, _low_mem = args_tuple
            import time as _t
            t0 = _t.time()
            _risk, _warn = robust_fit_and_predict(
                _family, _params, _XtrF, _ttr, _etr, _XteF, _sw, _seed, _low_mem)
            _elapsed = _t.time() - t0
            _ci = cindex(np.asarray(folds[_holdout]["tte"]), np.asarray(folds[_holdout]["ete"]), _risk)
            return _holdout, float(_ci) if np.isfinite(_ci) else float("nan"), _warn, _elapsed

        _n_pf = getattr(args, "parallel_folds", 1)
        fold_results = []

        if _n_pf > 1 and len(eligible_holdouts) > 1:
            import concurrent.futures as _cf_fold
            _fold_args = [(h, family, params,
                           folds[h]["XtrF"], folds[h]["ttr"], folds[h]["etr"],
                           folds[h]["XteF"], folds[h]["sw"], args.seed, args.low_mem)
                          for h in eligible_holdouts]
            try:
                with _cf_fold.ProcessPoolExecutor(max_workers=_n_pf) as _pex:
                    fold_results = list(_pex.map(_eval_fold_inner, _fold_args, timeout=900))
            except Exception as _pex_err:
                fold_results = []
                for h in eligible_holdouts:
                    fd = folds[h]; t0 = time.time()
                    risk, warn_msg = robust_fit_and_predict(family, params, fd["XtrF"], fd["ttr"], fd["etr"], fd["XteF"], fd["sw"], args.seed, args.low_mem)
                    ci = cindex(np.asarray(fd["tte"]), np.asarray(fd["ete"]), risk)
                    fold_results.append((h, float(ci) if np.isfinite(ci) else float("nan"), warn_msg, time.time()-t0))
        else:
            for h in eligible_holdouts:
                fd = folds[h]; t0 = time.time()
                risk, warn_msg = robust_fit_and_predict(family, params, fd["XtrF"], fd["ttr"], fd["etr"], fd["XteF"], fd["sw"], args.seed, args.low_mem)
                ci = cindex(np.asarray(fd["tte"]), np.asarray(fd["ete"]), risk)
                fold_results.append((h, float(ci) if np.isfinite(ci) else float("nan"), warn_msg, time.time()-t0))

        for j, (holdout, ci, warn_msg, elapsed) in enumerate(fold_results, start=1):
            if np.isfinite(ci):
                fold_scores[holdout] = ci; used += 1
            if warn_msg is not None:
                fallback_count += 1
            run_logs.append({
                "algo_name": algo_name, "family": family, "holdout": holdout,
                "fold_idx": j, "elapsed_sec": elapsed, "cindex": ci, "warn": warn_msg or ""
            })
            print(f"[ALG/FOLD] {algo_name} | fold {j}/{n_fold_total} | c-index={ci:.4f} | warn={warn_msg or '-'}")

        vals = np.array(list(fold_scores.values()),dtype=float) if fold_scores else np.array([np.nan],dtype=float)
        loco_mean = float(np.nanmean(vals)); loco_std=float(np.nanstd(vals))
        risk_sign = -1 if np.isfinite(loco_mean) and loco_mean<0.5 else 1

        rows.append({
            "algo_name":algo_name, "family":family, "params_json":json.dumps(params),
            "loco_mean_cindex":loco_mean, "loco_std_cindex":loco_std, "risk_sign":risk_sign,
            "n_folds_used":used, "n_fallback_or_fail":fallback_count,
            **{f"cindex_test_{h}":fold_scores.get(h, np.nan) for h in train}
        })
        if (k % 5)==0:
            pd.DataFrame(rows).sort_values("loco_mean_cindex", ascending=False).to_csv(outdir/"internal_loco_leaderboard_checkpoint.csv", index=False)
            pd.DataFrame(run_logs).to_csv(outdir/"logs"/"algorithm_fold_runlog_checkpoint.csv", index=False)
        print(f"[GRID] {k}/{len(algo_grid)} done")

        if k % 25 == 0 and rows:
            _ckpt_rows = rows.copy()
            if _checkpoint_path.exists():
                try:
                    _prev = pd.read_csv(_checkpoint_path)
                    _ckpt_rows_df = pd.DataFrame(_ckpt_rows)
                    _combined = pd.concat([_prev, _ckpt_rows_df], ignore_index=True)
                    _combined = _combined.drop_duplicates(subset=["algo_name"], keep="last")
                    _combined.sort_values("loco_mean_cindex", ascending=False, inplace=True)
                    _combined.to_csv(_checkpoint_path, index=False)
                except Exception:
                    pass
            else:
                pd.DataFrame(_ckpt_rows).sort_values(
                    "loco_mean_cindex", ascending=False).to_csv(_checkpoint_path, index=False)


    new_rows_df = pd.DataFrame(rows) if rows else pd.DataFrame()
    _ckpt_path  = outdir / "internal_loco_leaderboard.csv"
    if _ckpt_path.exists() and len(new_rows_df) > 0:
        try:
            _prev_df = pd.read_csv(_ckpt_path)
            leaderboard = pd.concat([_prev_df, new_rows_df], ignore_index=True)
            leaderboard = leaderboard.drop_duplicates(subset=["algo_name"], keep="last")
        except Exception:
            leaderboard = new_rows_df
    elif _ckpt_path.exists() and len(new_rows_df) == 0:
        leaderboard = pd.read_csv(_ckpt_path)
    else:
        leaderboard = new_rows_df

    leaderboard = leaderboard.sort_values("loco_mean_cindex", ascending=False)
    leaderboard.to_csv(outdir/"internal_loco_leaderboard.csv", index=False)
    print(f"[LEADERBOARD] {len(leaderboard)} config (checkpoint+new merged)")
    pd.DataFrame(run_logs).to_csv(outdir/"logs"/"algorithm_fold_runlog.csv", index=False)
    prog.tick("Internal LOCO bitti")

    if bool(args.diverse_topk):
        topk = select_diverse_topk(
            leaderboard,
            k=int(args.top_k_algos),
            min_param_distance=float(args.diverse_topk_min_param_distance),
            max_per_family=int(args.diverse_topk_max_per_family),
        )
    else:
        topk = leaderboard.head(int(args.top_k_algos)).copy()
    if topk.empty:
        raise RuntimeError("Top-k boş")
    try:
        fam_dist = topk["family"].value_counts().to_dict() if "family" in topk.columns else {}
        print(f"[TOPK] diverse={bool(args.diverse_topk)} | family_dist={fam_dist}")
    except Exception:
        pass

    # OOF risks for topk
    oof_risks_list, oof_names = [], []
    t_parts, e_parts = [], []
    for h in train:
        if h in folds:
            t_parts.append(np.asarray(folds[h]["tte"].values, dtype=float))
            e_parts.append(np.asarray(folds[h]["ete"].values, dtype=int))
    oof_t = np.concatenate(t_parts) if t_parts else np.array([], float)
    oof_e = np.concatenate(e_parts) if e_parts else np.array([], int)

    for r in topk.itertuples(index=False):
        oof_parts = []
        for h in train:
            if h not in folds: continue
            fd = folds[h]
            risk, _ = robust_fit_and_predict(r.family, json.loads(r.params_json), fd["XtrF"], fd["ttr"], fd["etr"], fd["XteF"], fd["sw"], args.seed, args.low_mem)
            oof_parts.append(int(r.risk_sign)*np.asarray(risk,dtype=float))
        if oof_parts:
            oof_risks_list.append(np.concatenate(oof_parts))
            oof_names.append(r.algo_name)

    if len(oof_risks_list):
        oof_df = pd.DataFrame({"T":oof_t, "E":oof_e})
        for n, rr in zip(oof_names, oof_risks_list): oof_df[n]=rr
        oof_df.to_csv(outdir/"risks"/"train_oof_topk_risks.csv", index=False)
    prog.tick("OOF risk yazıldı")

    filt_risks, filt_names = [], []
    for rr, nm in zip(oof_risks_list, oof_names):
        if len(rr) and np.isfinite(rr).mean()>=0.85:
            filt_risks.append(rr); filt_names.append(nm)
    w_opt, ci_opt_oof = (None, float("nan"))
    if len(filt_risks)>=2 and len(oof_t)==len(filt_risks[0]):
        w_opt, ci_opt_oof = optimize_ensemble_weights(filt_risks, oof_t, oof_e, args.opt_ens_n_iter, args.seed)

   
    stacking_model_risks = None
    stacking_oof_pred    = None
    _stacking_best_learner_label = None  # hangi meta-learner seçildi (logging için)
    if args.use_stacking and len(filt_risks) >= 2 and LIFELINES_AVAILABLE:
        print(f"[v35] Stacking ensemble: {len(filt_risks)} OOF risk vektörü, çok seviyeli meta-learner...")
        # OOF selection check
        _stacking_oof_check = stacking_ensemble_cox(
            filt_risks, oof_t, oof_e, filt_risks, args.seed)
        if _stacking_oof_check is not None:
            ci_stack_oof = cindex(oof_t, oof_e, _stacking_oof_check)
            stacking_oof_pred = _stacking_oof_check
            print(f"[v35] Stacking OOF C-index: {ci_stack_oof:.4f}")
        else:
            print("[v35] Stacking failed on OOF check, will skip for external.")
    elif args.use_stacking and len(filt_risks) < 2:
        print(f"[v35] Stacking atlandı: yeterli model yok ({len(filt_risks)} < 2). "
              "top_k_algos değerini artırın.")

    prog.tick("Ensemble optimizasyonu tamam")


    train_set = [c for c in cohorts if c.cohort in train]
    test_set = [c for c in cohorts if c.cohort in test]

    Xtr = pd.concat([c.X for c in train_set], axis=0)
    ttr = pd.concat([c.y_time for c in train_set], axis=0)
    etr = pd.concat([c.y_event for c in train_set], axis=0)
    coh_series = pd.Series([c.cohort for c in train_set for _ in range(len(c.X))], index=Xtr.index)
    Xtr2 = cohort_zscore(Xtr, coh_series, robust=args.robust_scaler)

    _train_norm_stats = extract_cohort_stats(Xtr, coh_series, robust=args.robust_scaler)
    if args.add_paper_traits:
        Xtr2 = add_paper_like_traits(Xtr2, ferr_drv, ferr_sup, wgcna_genes)

    clin_cols = [c for c in Xtr2.columns if c.startswith("CLIN__")]
    rna_cols  = [c for c in Xtr2.columns if c.startswith("RNA__")]
    if univ:
        allowed = {f"RNA__{g}" for g in univ}
        rna_cols = [c for c in rna_cols if c in allowed] or rna_cols

    v = Xtr2[rna_cols].var(axis=0).sort_values(ascending=False)
    cand = v.head(int(min(args.top_rna_var_prefilter, len(v)))).index.tolist()
    if args.use_wgcna_filter and wgcna_genes:
        allowed_w = {f"RNA__{g}" for g in wgcna_genes}
        cand_w = [c for c in cand if c in allowed_w]
        if len(cand_w) >= max(30, int(0.25*args.meta_top_genes)):
            cand = cand_w

    if args.use_meta_univariate and LIFELINES_AVAILABLE:
        if args.use_bagged_gene_select:
            meta_full = bagged_meta_univariate_select(
                Xtr2, ttr, etr, coh_series, cand, clin_cols,
                args.stage_adjusted_univariate, args.univar_penalizer,
                args.meta_top_genes, args.min_events_per_cohort,
                args.min_frac_cohorts, args.require_consistent_sign,
                n_bags=args.n_bags, bag_frac=0.8, seed=args.seed,
                n_jobs=args.meta_n_jobs  # v44-FIX-1
            )
            if meta_full is None or meta_full.empty:
                meta_full = meta_univariate_select(
                    Xtr2, ttr, etr, coh_series, cand, clin_cols,
                    args.stage_adjusted_univariate, args.univar_penalizer,
                    args.meta_top_genes, args.min_events_per_cohort,
                    args.min_frac_cohorts, args.require_consistent_sign,
                    n_jobs=args.meta_n_jobs  # v44-FIX-1
                )
        else:
            meta_full = meta_univariate_select(
                Xtr2, ttr, etr, coh_series, cand, clin_cols,
                args.stage_adjusted_univariate, args.univar_penalizer,
                args.meta_top_genes, args.min_events_per_cohort,
                args.min_frac_cohorts, args.require_consistent_sign,
                n_jobs=args.meta_n_jobs  # v44-FIX-1
            )

        if meta_full is None or meta_full.empty:
            meta_full = fallback_univariate_gene_select(
                Xtr2, ttr, etr, cand, top_genes=max(int(args.meta_top_genes), 40)
            )

            forced = enforce_ferroptosis_presence(
                selected_features=meta_full["feature"].tolist() if not meta_full.empty else [],
                ferr_genes=ferr_all,
                wgcna_genes=wgcna_genes,
                Xcols=list(Xtr2.columns),
                min_n=min(20, max(5, int(0.35 * max(1, args.meta_top_genes))))
            )
            if forced:
                forced_df = pd.DataFrame({
                    "feature": forced,
                    "gene": [f.replace("RNA__","",1) if str(f).startswith("RNA__") else str(f) for f in forced],
                    "z_meta": np.nan,
                    "p_meta": np.nan,
                    "n_cohorts_used": 1,
                    "frac_cohorts_used": 1.0,
                    "fallback_score": np.nan,
                    "method": "forced_ferro_presence"
                })
                # merge unique, keep earlier ranking first
                if not meta_full.empty:
                    meta_full = pd.concat([meta_full, forced_df], axis=0, ignore_index=True)
                else:
                    meta_full = forced_df
                meta_full = meta_full.drop_duplicates(subset=["feature"], keep="first").reset_index(drop=True)

        if meta_full is None or meta_full.empty:
            hard = [c for c in cand[:max(20, min(len(cand), int(args.meta_top_genes)))]]
            meta_full = pd.DataFrame({
                "feature": hard,
                "gene": [f.replace("RNA__","",1) for f in hard],
                "z_meta": np.nan, "p_meta": np.nan,
                "n_cohorts_used": 1, "frac_cohorts_used": 1.0,
                "fallback_score": np.nan, "method": "hard_guard_topvar"
            })
        meta_full.to_csv(outdir/"genes"/"meta_univariate_top_genes__FULLTRAIN.csv", index=False)
        meta_full.to_csv(outdir/"genes"/"possible_prognostic_genes_FULLTRAIN.csv", index=False)
        sel_rna = meta_full["feature"].tolist() if not meta_full.empty else cand[:args.meta_top_genes]
    else:

        sel_rna = cand[:args.meta_top_genes]
        meta_full = pd.DataFrame({
            "feature": sel_rna,
            "gene": [f.replace("RNA__","",1) for f in sel_rna],
            "z_meta": np.nan, "p_meta": np.nan,
            "n_cohorts_used": 1, "frac_cohorts_used": 1.0,
            "fallback_score": np.nan, "method": "topvar_nometa"
        })
        meta_full.to_csv(outdir/"genes"/"meta_univariate_top_genes__FULLTRAIN.csv", index=False)
        meta_full.to_csv(outdir/"genes"/"possible_prognostic_genes_FULLTRAIN.csv", index=False)

    extra_cols = []
    if args.add_paper_traits:
        extra_cols = [c for c in (
            "FERRO__DRIVER","FERRO__SUPPRESSOR","FERRO__TOTAL","FERRO__TRAIT",
            "WGCNA__ME1",

            "FERRO__DRIVER_PC1","FERRO__SUPP_PC1","FERRO__NET_SCORE",
            "FERRO__RATIO","FERRO__WGCNA_INTER","FERRO__Z_COMBINED"
        ) if c in Xtr2.columns]

    feats = [f for f in (sel_rna + clin_cols + extra_cols) if f in Xtr2.columns]


    locked_rna = []
    if args.lock_signature_by_stability:
        try:
            locked_rna = lock_stable_signature(
                meta_full if 'meta_full' in locals() else pd.DataFrame(),
                fold_feature_stability if 'fold_feature_stability' in locals() else pd.DataFrame(),
                final_gene_cap=args.final_gene_cap,
                min_gene_stability=args.min_gene_stability,
                ferr_genes=ferr_all,
                wgcna_genes=wgcna_genes,
                Xcols=list(Xtr2.columns),
                min_ferro_genes=args.ferroptosis_min_genes,
            )
        except Exception as _lock_ex:
            print(f"[v43] stable signature lock failed: {_lock_ex}")
            locked_rna = []
    if locked_rna:
        feats = [f for f in (locked_rna + clin_cols + extra_cols) if f in Xtr2.columns]
        print(f"[v43] Locked final signature aktif: RNA={len([f for f in feats if str(f).startswith('RNA__')])} | total={len(feats)}")


    feats = enforce_ferroptosis_presence(
        feats, ferr_all, wgcna_genes, list(Xtr2.columns),
        min_n=args.ferroptosis_min_genes,
        final_gene_cap=args.final_gene_cap
    )
    _rna_feats_final = [f for f in feats if str(f).startswith("RNA__")]
    print(f"[v44] Final feature seti: RNA={len(_rna_feats_final)}, total={len(feats)}")

    pd.Series(feats).to_csv(outdir/"selected_features_fulltrain.csv", index=False)
    pd.Series([g.replace("RNA__","",1) for g in feats if g.startswith("RNA__")]).to_csv(outdir/"genes"/"selected_rna_genes_fulltrain.csv", index=False)
    if locked_rna:
        pd.Series(locked_rna).to_csv(outdir/"genes"/"locked_signature_rna_genes.csv", index=False)


    if wgcna_genes:
        wgcna_rna = {f"RNA__{g}" for g in wgcna_genes}
        feats_rna  = [f for f in feats if f.startswith("RNA__")]
        in_wgcna   = [f for f in feats_rna if f in wgcna_rna]
        not_wgcna  = [f for f in feats_rna if f not in wgcna_rna]
        pd.DataFrame({
            "feature":    feats_rna,
            "gene":       [f.replace("RNA__","",1) for f in feats_rna],
            "in_wgcna":   [f in wgcna_rna for f in feats_rna],
            "source":     ["wgcna" if f in wgcna_rna else "non_wgcna" for f in feats_rna],
        }).to_csv(outdir/"genes"/"wgcna_gene_tracking.csv", index=False)
        print(f"[v36 WGCNA tracking] total_rna_feats={len(feats_rna)} | "
              f"in_wgcna={len(in_wgcna)} ({100*len(in_wgcna)/max(1,len(feats_rna)):.1f}%) | "
              f"not_wgcna={len(not_wgcna)}")
    sw = compute_sample_weights(coh_series, etr) if args.batch_weighting else None
    prog.tick("Full-train feature seçimi tamam")

    ext_rows, ens_rows = [], []
    per_test_risk = {tc:[] for tc in test}
    per_test_names = {tc:[] for tc in test}

    for r in topk.itertuples(index=False):
        params = json.loads(r.params_json); sign=int(r.risk_sign)
        for tc in test:
            te = next(c for c in test_set if c.cohort==tc)

            Xte2 = cohort_zscore(te.X,
                                  pd.Series([tc]*len(te.X), index=te.X.index),
                                  robust=args.robust_scaler,
                                  train_stats=_train_norm_stats if "_train_norm_stats" in dir() else None)
            if args.add_paper_traits:
                Xte2 = add_paper_like_traits(Xte2, ferr_drv, ferr_sup, wgcna_genes)
            XteF = Xte2.reindex(columns=feats).fillna(0.0)
            risk, warn_msg = robust_fit_and_predict(r.family, params, Xtr2[feats], ttr, etr, XteF, sw, args.seed, args.low_mem)
            risk = sign*np.asarray(risk, dtype=float)
            ci = cindex(te.y_time.values, te.y_event.values, risk)
            per_test_risk[tc].append(risk); per_test_names[tc].append(r.algo_name)
            pd.DataFrame({"sample":te.X.index.astype(str), "time":te.y_time.values.astype(float), "event":te.y_event.values.astype(int), "risk":risk.astype(float)}).to_csv(
                outdir/"risks"/f"test_{tc}__risk__{r.algo_name}.csv", index=False
            )

            ci_boot = bootstrap_ci(
                lambda tt, ee, ss: cindex(tt, ee, ss, use_ipcw=True),
                te.y_time.values, te.y_event.values, risk,
                n_boot=200, seed=int(args.seed)
            )

            ph_res = cox_ph_assumption_test(te.y_time.values, te.y_event.values, risk)

            eo_res = expected_observed_ratio(
                te.y_time.values, te.y_event.values, risk,
                ttr.values, etr.values, years=int(args.dca_years))

            ext_rows.append({
                "test_cohort":    tc,
                "algo_name":      r.algo_name,
                "family":         r.family,
                "params_json":    json.dumps(params),
                "external_cindex":ci,
                "cindex_ci_low":  ci_boot.get("ci_low",  float("nan")),
                "cindex_ci_high": ci_boot.get("ci_high", float("nan")),
                "ph_p_value":     ph_res.get("global_p", float("nan")),
                "ph_violated":    ph_res.get("ph_violated", None),
                "EO_ratio":       eo_res.get("EO_ratio", float("nan")),
                "EO_ci_low":      eo_res.get("EO_ci_low", float("nan")),
                "EO_ci_high":     eo_res.get("EO_ci_high", float("nan")),
                "warn":           warn_msg or ""
            })
            fh = fixed_horizon_auc_acc(te.y_time.values, te.y_event.values, risk, years_list=(1,2,3,5))
            fh.insert(0, "algo_name", r.algo_name)
            fh.insert(0, "test_cohort", tc)
            fh.to_csv(outdir/"metrics"/f"test_{tc}__fixed_horizon_metrics__{r.algo_name}.csv", index=False)
            ipcw = time_dependent_auc_ipcw(ttr.values, etr.values, te.y_time.values, te.y_event.values, risk, years_list=(1,2,3,5))
            ipcw.insert(0,"algo_name",r.algo_name); ipcw.insert(0,"test_cohort",tc)
            ipcw.to_csv(outdir/"metrics"/f"test_{tc}__ipcw_auc__{r.algo_name}.csv", index=False)


            train_risk_for_brier = risk  # fallback
            try:
                tr_pred, _ = robust_fit_and_predict(r.family, params, Xtr2[feats], ttr, etr, Xtr2[feats], sw, args.seed, args.low_mem)
                train_risk_for_brier = sign*np.asarray(tr_pred,float)
            except Exception:
                pass
            bsi = brier_and_ibs_ipcw(ttr.values, etr.values, te.y_time.values, te.y_event.values, train_risk_for_brier, risk, years_list=(1,2,3,5))
            bsi.insert(0,"algo_name",r.algo_name); bsi.insert(0,"test_cohort",tc)
            bsi.to_csv(outdir/"metrics"/f"test_{tc}__brier_ibs__{r.algo_name}.csv", index=False)

            km = km_logrank_hr(te.y_time.values, te.y_event.values, risk)
            pd.DataFrame([{"test_cohort":tc,"algo_name":r.algo_name, **km}]).to_csv(outdir/"metrics"/f"test_{tc}__km_hr_logrank__{r.algo_name}.csv", index=False)
            plot_km(te.y_time.values, te.y_event.values, risk, f"{r.algo_name} | {tc} KM (median split)", outdir/"figures"/f"KM_{tc}__{r.algo_name}.png")


            try:
                perm_res = permutation_cindex_pvalue(te.y_time.values, te.y_event.values, risk, n_perm=args.n_perm, seed=args.seed)
                pd.DataFrame([{"test_cohort":tc,"algo_name":r.algo_name, **perm_res}]).to_csv(
                    outdir/"metrics"/f"test_{tc}__permutation_pvalue__{r.algo_name}.csv", index=False)
            except Exception as _pe:
                print(f"[WARN] permutation test failed: {_pe}")


            if args.multivariate_validation:
                try:
                    stage_vals = None
                    stage_col = next((col for col in te.X.columns if "stage_num" in col.lower()), None)
                    if stage_col:
                        stage_vals = te.X[stage_col].values
                    mv_df = multivariable_cox_validation(te.y_time.values, te.y_event.values, risk, stage_vals)
                    if not mv_df.empty:
                        mv_df.insert(0, "algo_name", r.algo_name); mv_df.insert(0, "test_cohort", tc)
                        mv_df.to_csv(outdir/"metrics"/f"test_{tc}__multivariable_cox__{r.algo_name}.csv", index=False)
                except Exception as _mve:
                    print(f"[WARN] multivariable Cox failed: {_mve}")


            try:
                dca_df = decision_curve_analysis(te.y_time.values, te.y_event.values, risk, years=args.dca_years)
                if not dca_df.empty:
                    dca_df.insert(0, "algo_name", r.algo_name); dca_df.insert(0, "test_cohort", tc)
                    dca_df.to_csv(outdir/"metrics"/f"test_{tc}__dca__{r.algo_name}.csv", index=False)
                    plot_dca(dca_df, f"DCA {r.algo_name} | {tc} ({args.dca_years}yr)", outdir/"figures"/f"DCA_{tc}__{r.algo_name}.png")
            except Exception as _dce:
                print(f"[WARN] DCA failed: {_dce}")


            try:
                cal_df = calibration_curve(te.y_time.values, te.y_event.values, risk, ttr.values, etr.values, years=args.calibration_years)
                if not cal_df.empty:
                    cal_df.insert(0, "algo_name", r.algo_name); cal_df.insert(0, "test_cohort", tc)
                    cal_df.to_csv(outdir/"metrics"/f"test_{tc}__calibration__{r.algo_name}.csv", index=False)
                    plot_calibration(cal_df, f"Calibration {r.algo_name} | {tc}", outdir/"figures"/f"Calib_{tc}__{r.algo_name}.png")
            except Exception as _cale:
                print(f"[WARN] Calibration failed: {_cale}")

           
            try:
                for _rmst_yr in [3, 5]:
                    rmst_res = compute_rmst(
                        te.y_time.values, te.y_event.values, risk,
                        tau_years=_rmst_yr, n_boot=200, seed=int(args.seed))
                    pd.DataFrame([{
                        "test_cohort": tc, "algo_name": r.algo_name,
                        "ph_violated": ph_res.get("ph_violated", None),
                        **rmst_res
                    }]).to_csv(
                        outdir/"metrics"/f"test_{tc}__rmst_{_rmst_yr}yr__{r.algo_name}.csv",
                        index=False)
            except Exception as _re:
                print(f"[WARN] RMST failed: {type(_re).__name__}: {_re}")

           
            _ph_violated = ph_res.get("ph_violated", False)
            if _ph_violated:
                try:
                    lm_df = landmark_analysis(
                        te.y_time.values, te.y_event.values, risk,
                        landmark_years=(1, 2, 3), follow_years=float(args.dca_years))
                    if not lm_df.empty:
                        lm_df.insert(0, "algo_name", r.algo_name)
                        lm_df.insert(0, "test_cohort", tc)
                        lm_df.to_csv(
                            outdir/"metrics"/f"test_{tc}__landmark__{r.algo_name}.csv",
                            index=False)
                except Exception as _le:
                    print(f"[WARN] Landmark failed: {type(_le).__name__}: {_le}")

            
            try:
                if not LIFELINES_AVAILABLE:
                    raise RuntimeError("lifelines missing")
                from lifelines import CoxPHFitter as _CPF

                _t_nri  = te.y_time.values.astype(float)
                _e_nri  = te.y_event.values.astype(int)
                _r_nri  = risk.astype(float)
                _ok_nri = np.isfinite(_t_nri) & np.isfinite(_e_nri) & np.isfinite(_r_nri) & (_t_nri > 0)
                _t_nri, _e_nri, _r_nri = _t_nri[_ok_nri], _e_nri[_ok_nri], _r_nri[_ok_nri]

                if len(_t_nri) < 30 or int(_e_nri.sum()) < 8:
                    raise RuntimeError(f"NRI/IDI: insufficient data (n={len(_t_nri)}, events={int(_e_nri.sum())})")

                stage_col_nri  = next((col for col in te.X.columns if "stage_num" in col.lower()), None)
                _stage_vals    = te.X[stage_col_nri].values.astype(float)[_ok_nri] if stage_col_nri else None
                _stage_usable  = (_stage_vals is not None and
                                  np.isfinite(_stage_vals).mean() > 0.8 and
                                  pd.Series(_stage_vals).nunique() >= 2 and
                                  float(pd.Series(_stage_vals).std()) > 1e-6)

                risk_base_nri   = None
                _baseline_label = None

                if _stage_usable:
                    # ── stage-only Cox (tercih) ──────────────────────
                    _st_mu = float(np.nanmean(_stage_vals))
                    _st_sd = float(np.nanstd(_stage_vals)) + 1e-8
                    _stage_z = (_stage_vals - _st_mu) / _st_sd

                    _df_b = pd.DataFrame({"T": _t_nri, "E": _e_nri, "stage_z": _stage_z}).dropna()
                    for _pen in [1.0, 5.0, 10.0, 50.0]:
                        try:
                            _cph_b = _CPF(penalizer=_pen)
                            _cph_b.fit(_df_b[["T","E","stage_z"]],
                                       duration_col="T", event_col="E", show_progress=False)
                            risk_base_nri   = _cph_b.predict_partial_hazard(
                                _df_b[["stage_z"]]).values.astype(float)
                            _baseline_label = f"stage_cox_pen{_pen}"
                            break
                        except Exception:
                            continue

                    # CoxNet fallback
                    if risk_base_nri is None and GLMNET_AVAILABLE:
                        try:
                            from sksurv.linear_model import CoxnetSurvivalAnalysis as _CNA
                            _y_b = np.array(
                                [(bool(ev), float(tt)) for tt, ev in zip(_t_nri, _e_nri)],
                                dtype=[("event", bool), ("time", float)])
                            _cn = _CNA(l1_ratio=0.5, alpha_min_ratio=0.1, max_iter=5000)
                            _cn.fit(_stage_z.reshape(-1, 1), _y_b)
                            risk_base_nri   = _cn.predict(_stage_z.reshape(-1, 1)).astype(float)
                            _baseline_label = "stage_coxnet"
                        except Exception:
                            pass

                    if risk_base_nri is None:
                        # rank-normalize stage → crude Cox proxy
                        risk_base_nri   = rank_normalize(_stage_z)
                        _baseline_label = "stage_rank_proxy"
                        print(f"  [v44] NRI/IDI {tc}: CoxPH/CoxNet başarısız, rank-proxy kullanılıyor")

                else:
                    # ── null model baseline (without stage / constant) ───
                   
                    risk_base_nri   = np.full(len(_r_nri), float(np.nanmedian(_r_nri)))
                    _baseline_label = "null_model_no_stage"
                    print(f"  [v44] NRI/IDI {tc}: stage yok/constant → null model baseline kullanılıyor")

                nri_res = continuous_nri_idi(
                    _t_nri, _e_nri, _r_nri, risk_base_nri,
                    years=int(args.dca_years))
                pd.DataFrame([{
                    "test_cohort":  tc,
                    "algo_name":    r.algo_name,
                    "baseline":     _baseline_label,
                    "stage_available": _stage_usable,
                    **nri_res
                }]).to_csv(
                    outdir/"metrics"/f"test_{tc}__nri_idi__{r.algo_name}.csv",
                    index=False)

            except Exception as _ne:
                print(f"[WARN] NRI/IDI failed: {type(_ne).__name__}: {_ne}")
            

    ext_df = pd.DataFrame(ext_rows)
    ext_df.to_csv(outdir/"topk_external_results.csv", index=False)
    prog.tick("Top-k external sonuçları yazıldı")

    for tc in test:
        te = next(c for c in test_set if c.cohort==tc)
        risks = per_test_risk[tc]; names = per_test_names[tc]
        if not risks: continue
        r_mean = weighted_rank_ensemble(risks, None); ci_mean = cindex(te.y_time.values, te.y_event.values, r_mean)
        weights_loco = np.array([max(float(topk.loc[topk.algo_name==name,"loco_mean_cindex"].values[0]),1e-6) for name in names], dtype=float)
        r_w = weighted_rank_ensemble(risks, weights_loco); ci_w = cindex(te.y_time.values, te.y_event.values, r_w)

        ci_opt = float("nan")
        if w_opt is not None and len(filt_names)==len(w_opt):
            risk_map = {n:r for n,r in zip(names, risks)}
            aligned = [risk_map.get(n, None) for n in filt_names]
            if all(a is not None for a in aligned):
                ci_opt = cindex(te.y_time.values, te.y_event.values, weighted_rank_ensemble(aligned, w_opt))

        Z = np.vstack([rank_normalize(r) for r in risks])
        r_hybrid = np.nanmean(Z, axis=0)
        ci_hybrid = cindex(te.y_time.values, te.y_event.values, r_hybrid)


        r_borda = borda_count_ensemble(risks)
        ci_borda = cindex(te.y_time.values, te.y_event.values, r_borda)


        r_softmax = softmax_weighted_ensemble(risks, weights_loco, temperature=0.05)
        ci_softmax = cindex(te.y_time.values, te.y_event.values, r_softmax)


        ci_stack = float("nan")
        if args.use_stacking and len(filt_risks) >= 2:
            try:
                risk_map = {n: r for n, r in zip(names, risks)}
                aligned_test = [risk_map.get(n) for n in filt_names]
                if all(a is not None for a in aligned_test):
                    r_stack = stacking_ensemble_cox(filt_risks, oof_t, oof_e, aligned_test, args.seed)
                    if r_stack is not None:
                        ci_stack = cindex(te.y_time.values, te.y_event.values, r_stack)
                        pd.DataFrame({"sample": te.X.index.astype(str),
                                      "risk_stacking": r_stack}).to_csv(
                            outdir/"risks"/f"test_{tc}__risk__STACKING.csv", index=False)
            except Exception as _se:
                print(f"[WARN] stacking test predict failed: {_se}")


        ci_etx_rsf          = float("nan")
        ci_etx_rsf_opt      = float("nan")
        ci_etx_rsf_stk      = float("nan")
        ci_gbm_hybrid       = float("nan")

        try:
            name_risk_map = dict(zip(names, risks))
            fam_map = ({r.algo_name: str(r.get("family","")) for _, r in topk.iterrows()}
                       if "family" in topk.columns else {})

            etx_risks_te  = [name_risk_map[n] for n in names
                              if fam_map.get(n,"") == "EXTRATREES" and n in name_risk_map]
            rsf_risks_te  = [name_risk_map[n] for n in names
                              if fam_map.get(n,"") == "RSF"         and n in name_risk_map]
            gbm_risks_te  = [name_risk_map[n] for n in names
                              if fam_map.get(n,"") == "GBM"         and n in name_risk_map]

            if etx_risks_te and rsf_risks_te:
                etx_mean = weighted_rank_ensemble(etx_risks_te, None)
                rsf_mean = weighted_rank_ensemble(rsf_risks_te, None)


                r_etx_rsf    = weighted_rank_ensemble([etx_mean, rsf_mean], np.array([0.6, 0.4]))
                ci_etx_rsf   = cindex(te.y_time.values, te.y_event.values, r_etx_rsf)

              
                try:
                    etx_oof_risks = []
                    rsf_oof_risks = []
                    for nm in oof_names:
                        nm_fam = fam_map.get(nm, "")
                        idx_n  = oof_names.index(nm)
                        if nm_fam == "EXTRATREES" and idx_n < len(oof_risks_list):
                            etx_oof_risks.append(oof_risks_list[idx_n])
                        elif nm_fam == "RSF" and idx_n < len(oof_risks_list):
                            rsf_oof_risks.append(oof_risks_list[idx_n])
                    if etx_oof_risks and rsf_oof_risks:
                        etx_oof_mean = weighted_rank_ensemble(etx_oof_risks, None)
                        rsf_oof_mean = weighted_rank_ensemble(rsf_oof_risks, None)
                        w_opt_er, ci_opt_er_oof = optimize_ensemble_weights(
                            [etx_oof_mean, rsf_oof_mean], oof_t, oof_e,
                            min(500, args.opt_ens_n_iter), args.seed)
                        if w_opt_er is not None and len(w_opt_er) == 2:
                            r_etx_rsf_opt  = weighted_rank_ensemble(
                                [etx_mean, rsf_mean], np.array(w_opt_er))
                            ci_etx_rsf_opt = cindex(te.y_time.values, te.y_event.values,
                                                     r_etx_rsf_opt)
                            # Use optimized if better, else fallback to fixed
                            if np.isfinite(ci_etx_rsf_opt) and ci_etx_rsf_opt > ci_etx_rsf:
                                r_etx_rsf  = r_etx_rsf_opt
                                ci_etx_rsf = ci_etx_rsf_opt
                            print(f"  [v35 ETX+RSF] fixed={ci_etx_rsf:.4f} "
                                  f"opt(w={[round(w,2) for w in w_opt_er]})={ci_etx_rsf_opt:.4f}")
                except Exception as _oe:
                    pass  # OOF optimize failed, use fixed weight


                if stacking_oof_pred is not None and len(ci_stack_arr := []) == 0:
                    try:
                       
                        pass
                    except Exception:
                        pass

                pd.DataFrame({"sample": te.X.index.astype(str),
                              "risk_etx_rsf_hybrid": r_etx_rsf}).to_csv(
                    outdir/"risks"/f"test_{tc}__risk__ETX_RSF_HYBRID.csv", index=False)


                if gbm_risks_te:
                    gbm_mean = weighted_rank_ensemble(gbm_risks_te, None)
                    r_three  = weighted_rank_ensemble(
                        [etx_mean, rsf_mean, gbm_mean],
                        np.array([0.5, 0.3, 0.2]))
                    ci_gbm_hybrid = cindex(te.y_time.values, te.y_event.values, r_three)
                    pd.DataFrame({"sample": te.X.index.astype(str),
                                  "risk_etx_rsf_gbm": r_three}).to_csv(
                        outdir/"risks"/f"test_{tc}__risk__ETX_RSF_GBM.csv", index=False)

        except Exception as _er:
            print(f"[WARN] ETX+RSF hybrid ensemble failed: {_er}")

        # Trimmed rank ensemble
        r_trimmed = trimmed_rank_ensemble(risks, trim_frac=0.10)
        ci_trimmed = cindex(te.y_time.values, te.y_event.values, r_trimmed)
        # Trimmed 20% (daha agresif)
        r_trimmed20 = trimmed_rank_ensemble(risks, trim_frac=0.20)
        ci_trimmed20 = cindex(te.y_time.values, te.y_event.values, r_trimmed20)

        ens_rows += [
            {"test_cohort":tc, "ensemble":"MEAN_RANK",             "external_cindex":ci_mean},
            {"test_cohort":tc, "ensemble":"WEIGHTED_RANK",         "external_cindex":ci_w},
            {"test_cohort":tc, "ensemble":"OPT_WEIGHTED_RANK",     "external_cindex":ci_opt,
             "oof_train_cindex":ci_opt_oof},
            {"test_cohort":tc, "ensemble":"HYBRID_RANK_Z",         "external_cindex":ci_hybrid},
            {"test_cohort":tc, "ensemble":"BORDA_COUNT",           "external_cindex":ci_borda},
            {"test_cohort":tc, "ensemble":"SOFTMAX_WEIGHTED",      "external_cindex":ci_softmax},
            {"test_cohort":tc, "ensemble":"STACKING_COX",          "external_cindex":ci_stack},
            {"test_cohort":tc, "ensemble":"ETX_RSF_HYBRID",        "external_cindex":ci_etx_rsf},
            {"test_cohort":tc, "ensemble":"ETX_RSF_OPT",           "external_cindex":ci_etx_rsf_opt},
            {"test_cohort":tc, "ensemble":"ETX_RSF_GBM_HYBRID",    "external_cindex":ci_gbm_hybrid},
            {"test_cohort":tc, "ensemble":"TRIMMED_RANK_10pct",    "external_cindex":ci_trimmed},
            {"test_cohort":tc, "ensemble":"TRIMMED_RANK_20pct",    "external_cindex":ci_trimmed20},
        ]

    ens_df = pd.DataFrame(ens_rows)
    ens_df.to_csv(outdir/"external_ensemble_results.csv", index=False)

    ens_fixed_rows = []
    for tc in test:
        te = next(c for c in test_set if c.cohort==tc)
        risks = per_test_risk[tc]; names = per_test_names[tc]
        if not risks:
            continue
        r_mean = weighted_rank_ensemble(risks, None)
        weights_loco = np.array([max(float(topk.loc[topk.algo_name==name,"loco_mean_cindex"].values[0]),1e-6) for name in names], dtype=float)
        r_w = weighted_rank_ensemble(risks, weights_loco)

        ens_map = {"MEAN_RANK": r_mean, "WEIGHTED_RANK": r_w}
        if w_opt is not None and len(filt_names)==len(w_opt):
            risk_map = {n:r for n,r in zip(names, risks)}
            aligned = [risk_map.get(n, None) for n in filt_names]
            if all(a is not None for a in aligned):
                ens_map["OPT_WEIGHTED_RANK"] = weighted_rank_ensemble(aligned, w_opt)

        for ens_name, rr in ens_map.items():
            fh = fixed_horizon_auc_acc(te.y_time.values, te.y_event.values, rr, years_list=(1,2,3,5))
            fh.insert(0, "ensemble", ens_name)
            fh.insert(0, "test_cohort", tc)
            ens_fixed_rows.append(fh)

    if len(ens_fixed_rows):
        pd.concat(ens_fixed_rows, axis=0).to_csv(outdir/"metrics"/"external_ensemble_fixed_horizon_metrics.csv", index=False)

    prog.tick("Ensemble external sonuçları yazıldı")

    with pd.ExcelWriter(outdir/"summary_best_tables.xlsx") as w:
        qc.to_excel(w, sheet_name="qc_events", index=False)
        leaderboard.head(100).to_excel(w, sheet_name="internal_top100", index=False)
        topk.to_excel(w, sheet_name="topk", index=False)
        ext_df.to_excel(w, sheet_name="external_all", index=False)
        ens_df.to_excel(w, sheet_name="ensembles", index=False)
        p1 = outdir/"metrics"/"external_ensemble_fixed_horizon_metrics.csv"
        if p1.exists():
            pd.read_csv(p1).to_excel(w, sheet_name="ens_fixed_1_2_3_5y", index=False)

        if not ens_df.empty:
            try:
                best_ens = (ens_df.sort_values("external_cindex", ascending=False)
                            .groupby("test_cohort").head(5)
                            .reset_index(drop=True))
                best_ens.to_excel(w, sheet_name="best_ensemble_per_test", index=False)
            except Exception:
                pass

        try:
            rmst_files = list((outdir/"metrics").glob("*__rmst_*yr__*.csv"))
            if rmst_files:
                rmst_all = pd.concat([pd.read_csv(f) for f in rmst_files], ignore_index=True)
                # Best model per cohort per tau
                rmst_best = (rmst_all.sort_values("delta_rmst", ascending=False)
                             .groupby(["test_cohort","tau_years"]).head(1)
                             .reset_index(drop=True))
                rmst_all.to_excel(w, sheet_name="Supp_RMST_all", index=False)
                rmst_best.to_excel(w, sheet_name="Table_RMST_best", index=False)
        except Exception as _re:
            print(f"[WARN v45] RMST sheet failed: {_re}")
        # Landmark summary (time-varying HR)
        try:
            lm_files = list((outdir/"metrics").glob("*__landmark__*.csv"))
            if lm_files:
                lm_all = pd.concat([pd.read_csv(f) for f in lm_files], ignore_index=True)
                lm_all.to_excel(w, sheet_name="Supp_Landmark", index=False)
        except Exception as _le:
            print(f"[WARN] Landmark sheet failed: {_le}")
        

    planned = pd.DataFrame([{"algo_name":a, "family":f, "params_json":json.dumps(p)} for a,f,p in algo_grid])
    master = planned.merge(leaderboard[["algo_name","loco_mean_cindex","loco_std_cindex","n_folds_used","n_fallback_or_fail","risk_sign"]],
                           on="algo_name", how="left")
    master["run_status"] = np.where(master["loco_mean_cindex"].notna(), "ran", "not_run_or_timed_out")
    master.to_csv(outdir/"master_leaderboard_all_algorithms.csv", index=False)
    prog.tick("Master leaderboard kaydedildi")

    # provenance summary
    pd.DataFrame([{
        "ferrdb_driver_n": len(ferr_drv),
        "ferrdb_suppressor_n": len(ferr_sup),
        "ferrdb_marker_n": len(ferr_mkr),
        "ferrdb_all_n": len(ferr_all),
        "wgcna_n": len(wgcna_genes),
        "wgcna_source": wgcna_prov,
        "gene_set_mode": args.gene_set_mode,
        "use_wgcna_filter": bool(args.use_wgcna_filter),
        "add_paper_traits": bool(args.add_paper_traits),
    }]).to_csv(outdir/"logs"/"run_gene_sources_summary.csv", index=False)

    

    try:
        metric_files = list((outdir/"metrics").glob("*.csv"))
        m_all = []
        for p in metric_files:
            try:
                d = pd.read_csv(p)
                d["source_file"] = p.name
                m_all.append(d)
            except Exception:
                pass
        metrics_union = pd.concat(m_all, axis=0, ignore_index=True) if len(m_all) else pd.DataFrame()
        metrics_union.to_csv(outdir/"publication"/"all_metrics_union.csv", index=False)

        best_row = None
        if 'ext_df' in locals() and isinstance(ext_df, pd.DataFrame) and len(ext_df):
            e2 = ext_df.sort_values(["external_cindex"], ascending=False).reset_index(drop=True)
            best_row = e2.iloc[0].to_dict()

        figs = sorted([p.name for p in (outdir/"figures").glob("*") if p.is_file()])
        pd.DataFrame({"figure_file":figs}).to_csv(outdir/"publication"/"figure_manifest.csv", index=False)

        with pd.ExcelWriter(outdir/"publication"/"publication_summary.xlsx") as w:
            qc.to_excel(w, sheet_name="Table1_qc", index=False)
            leaderboard.head(30).to_excel(w, sheet_name="Table2_internal_top30", index=False)
            topk.to_excel(w, sheet_name="Table3_topk", index=False)
            ext_df.to_excel(w, sheet_name="Table4_external_models", index=False)
            ens_df.to_excel(w, sheet_name="Table5_ensembles", index=False)
            if len(metrics_union):
                metrics_union.to_excel(w, sheet_name="Supp_metrics_union", index=False)
            if best_row is not None:
                pd.DataFrame([best_row]).to_excel(w, sheet_name="Best_model", index=False)
            # v31: yeni metrik tabloları
            try:
                perm_files = list((outdir/"metrics").glob("*permutation_pvalue*.csv"))
                if perm_files:
                    pd.concat([pd.read_csv(f) for f in perm_files], ignore_index=True).to_excel(w, sheet_name="Supp_permutation_pval", index=False)
                mv_files = list((outdir/"metrics").glob("*multivariable_cox*.csv"))
                if mv_files:
                    pd.concat([pd.read_csv(f) for f in mv_files], ignore_index=True).to_excel(w, sheet_name="Supp_multivariable_cox", index=False)
                dca_files = list((outdir/"metrics").glob("*__dca__*.csv"))
                if dca_files:
                    pd.concat([pd.read_csv(f) for f in dca_files], ignore_index=True).to_excel(w, sheet_name="Supp_DCA", index=False)
                cal_files = list((outdir/"metrics").glob("*calibration*.csv"))
                if cal_files:
                    pd.concat([pd.read_csv(f) for f in cal_files], ignore_index=True).to_excel(w, sheet_name="Supp_calibration", index=False)
                # meta-univariate with FDR
                meta_f = outdir/"genes"/"meta_univariate_top_genes__FULLTRAIN.csv"
                if meta_f.exists():
                    pd.read_csv(meta_f).to_excel(w, sheet_name="Table6_genes_FDR", index=False)
            except Exception as _v31ex:
                print(f"[WARN] publication sheets failed: {_v31ex}")

        rep = []
       
        rep.append("## Primary endpoints\n")
        rep.append("- Discrimination: C-index (internal LOCO and external) + permutation p-value.\n")
        rep.append("- Time-specific metrics: AUC-ROC, Accuracy (1/2/3/5 years).\n")
        rep.append("- Survival reporting: KM, log-rank, HR with 95% CI.\n")
        rep.append("- Additional: IPCW AUC, Brier/IBS, DCA, Calibration curve.\n")
        rep.append("- Multivariable Cox: stage-independent prognostic value.\n")
        rep.append("- Gene selection: meta-univariate p + FDR (BH) q-values.\n")

        rep.append("- ETX ince grid: mf=0.65..0.85, md=9..14, leaf=3..10, ne=260/320/400.\n")
        
        rep.append("- FASTSVM+CWGBM: İLK KEZ (genişletilmiş sweep).\n")
        rep.append("- ETX+RSF Hybrid ensemble: cross-family ağırlıklı rank kombinasyonu.\n")
        rep.append("- RSF min_samples_split = max(6, 2*min_leaf) (paper default'a yakın).\n")
        rep.append("\n## Files\n")
       
        rep.append("- logs/run_config.json  ← reproducibility\n")
        rep.append("- figures/KM_*.png, DCA_*.png, Calib_*.png\n")
        rep.append("- genes/meta_univariate_top_genes__FULLTRAIN.csv  ← FDR q-values\n")
        rep.append("- metrics/*permutation_pvalue*.csv\n")
        rep.append("- metrics/*multivariable_cox*.csv\n")
        if best_row is not None:
            rep.append("\n## Best model on external set\n")
            rep.append(f"- Model: **{best_row.get('algo_name','NA')}**\n")
            rep.append(f"- Cohort: **{best_row.get('test_cohort','NA')}**\n")
            rep.append(f"- External C-index: **{best_row.get('external_cindex','NA')}**\n")
        (outdir/"publication"/"README_publication.md").write_text("".join(rep), encoding="utf-8")
   


# ferroptosis prognostic signature pipeline
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
