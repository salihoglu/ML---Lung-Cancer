#!/usr/bin/env python3
"""
Rana Salihoglu

FerrDb-integrated target druggability and target-compound mapping pipeline
for ML-derived LUAD ferroptosis-related signature genes.

"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator

SEED = 20260418
np.random.seed(SEED)

DEFAULT_ROOT = "/home/rana/Desktop/PROJECT_16/test"
DEFAULT_MANIFESTS = [
    "/mnt/data/full_paths.txt",
    "/mnt/data/full_paths(2).txt",
    "full_paths.txt",
    "full_paths(2).txt",
]
PPI_TREE_RELATIVE = os.path.join(
    "bioinf_scripts", "scripts", "ppi_network_Q1_v2", "tables", "node_centrality.csv"
)

OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


# -----------------------------------------------------------------------------
# CLI and helpers
# -----------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Advanced Q1-grade FerrDb-integrated target druggability / compound mapping analysis"
    )
    p.add_argument("--project_root", default=DEFAULT_ROOT)
    p.add_argument("--manifest", default="")
    p.add_argument("--gene_file", default="possible_prognostic_genes_FULLTRAIN.csv")
    p.add_argument("--ppi_table", default="")
    p.add_argument("--ferrdb_driver_csv", default="ferrdb_driver.csv")
    p.add_argument("--ferrdb_suppressor_csv", default="ferrdb_suppressor.csv")
    p.add_argument("--ferrdb_marker_csv", default="ferrdb_marker.csv")
    p.add_argument("--outdir", default="")
    p.add_argument("--top_n", type=int, default=50)
    p.add_argument("--sleep_sec", type=float, default=0.10)
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--max_retries", type=int, default=3)
    p.add_argument("--max_chembl_activities_per_target", type=int, default=400)
    p.add_argument("--max_molecule_phase_queries_per_target", type=int, default=120)
    p.add_argument("--min_pchembl_for_phase_support", type=float, default=5.0)
    p.add_argument("--max_ot_known_drugs_per_target", type=int, default=150)
    p.add_argument("--min_priority_display", type=float, default=0.0)
    p.add_argument("--include_marker_genes", default="true")
    p.add_argument("--cache_dir", default="")
    return p.parse_args()


def msg(*x) -> None:
    print(*x, flush=True)


def as_bool(x: str) -> bool:
    return str(x).strip().lower() in {"1", "true", "yes", "y", "t"}


def normalize_path(x: str) -> str:
    try:
        return str(Path(os.path.expanduser(x)).resolve())
    except Exception:
        return os.path.expanduser(x)


def read_manifest(path: str) -> List[str]:
    if not path or not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def pick_manifest(user_value: str) -> str:
    cands = [user_value] if user_value else []
    cands.extend(DEFAULT_MANIFESTS)
    for c in cands:
        if c and os.path.exists(os.path.expanduser(c)):
            return normalize_path(c)
    return ""


def infer_root(project_root: str, manifest_paths: List[str]) -> str:
    direct = normalize_path(project_root)
    if os.path.exists(direct):
        return direct
    for p in manifest_paths:
        m = re.search(r"/(TCGA_GDC_download|GEO_PREP|OUT_v44|bioinf_scripts|luad_r_analysis_bundle)/", p)
        if m:
            return p[:m.start()]
    return direct


def resolve_existing_file(target: str, candidates: Iterable[str]) -> str:
    ordered: List[str] = []
    seen = set()
    for x in [target, *candidates]:
        if not x:
            continue
        xn = normalize_path(x)
        if xn not in seen:
            seen.add(xn)
            ordered.append(xn)
    for x in ordered:
        if os.path.isfile(x):
            return x
    base = os.path.basename(target)
    for x in ordered:
        alt = os.path.join(os.path.dirname(x), base)
        if os.path.isfile(alt):
            return normalize_path(alt)
    return ""


def resolve_from_manifest(manifest_paths: List[str], patterns: Sequence[str]) -> str:
    for pat in patterns:
        rgx = re.compile(pat, flags=re.IGNORECASE)
        hits = [p for p in manifest_paths if rgx.search(p)]
        if hits:
            return hits[0]
    return ""


def ensure_dirs(outdir: Path) -> Dict[str, Path]:
    paths = {
        "root": outdir,
        "tables": outdir / "tables",
        "figures": outdir / "figures",
        "logs": outdir / "logs",
        "cache": outdir / "cache",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def safe_read_table(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    ext = Path(path).suffix.lower()
    if ext in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def sanitize_gene(x: str) -> str:
    x = str(x).strip()
    if "__" in x:
        x = x.split("__")[-1]
    return x.replace(".", "-").upper()


def scaled_series(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce").astype(float)
    if s.notna().sum() == 0:
        return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    mn = float(s.min(skipna=True))
    mx = float(s.max(skipna=True))
    if (not np.isfinite(mn)) or (not np.isfinite(mx)) or mx <= mn:
        return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    return ((s - mn) / (mx - mn)).astype(float)


def neglog10(x: pd.Series, floor: float = 1e-300) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    return -np.log10(np.clip(x.astype(float), floor, None))


# -----------------------------------------------------------------------------
# Requests session with retries and caching
# -----------------------------------------------------------------------------
def get_requests_session(max_retries: int = 3):
    import requests
    from requests.adapters import HTTPAdapter
    try:
        from urllib3.util.retry import Retry
        retry = Retry(
            total=max_retries,
            read=max_retries,
            connect=max_retries,
            status=max_retries,
            backoff_factor=0.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET", "POST"]),
            raise_on_status=False,
        )
        adapter = HTTPAdapter(max_retries=retry)
    except Exception:
        adapter = HTTPAdapter()
    s = requests.Session()
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "Q1-target-druggability-pipeline/2026-04-18"})
    return s


def _cache_key(prefix: str, value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return f"{prefix}_{safe}.json"


def read_json_cache(cache_dir: Path, prefix: str, key: str) -> Optional[dict]:
    fp = cache_dir / _cache_key(prefix, key)
    if not fp.exists():
        return None
    try:
        return json.loads(fp.read_text(encoding="utf-8"))
    except Exception:
        return None


def write_json_cache(cache_dir: Path, prefix: str, key: str, payload: dict) -> None:
    fp = cache_dir / _cache_key(prefix, key)
    try:
        fp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception:
        pass


# -----------------------------------------------------------------------------
# FerrDb readers
# -----------------------------------------------------------------------------
def _pick_gene_col(df: pd.DataFrame) -> str:
    low = {str(c).lower(): c for c in df.columns}
    priority = ["gene", "symbol", "gene_symbol", "genesymbol", "hgnc_symbol", "official_symbol", "marker", "feature"]
    for k in priority:
        if k in low:
            return low[k]
    for c in df.columns:
        cl = str(c).lower()
        if "gene" in cl or "symbol" in cl:
            return c
    return df.columns[0]


def read_ferrdb_gene_set(path: str) -> Set[str]:
    df = safe_read_table(path)
    if df.empty:
        return set()
    gene_col = _pick_gene_col(df)
    genes = {sanitize_gene(x) for x in df[gene_col].dropna().astype(str)}
    return {g for g in genes if g and g not in {"NAN", "NONE"}}


def resolve_ferrdb_file(root: str, manifest_paths: List[str], user_value: str) -> str:
    base = os.path.basename(user_value)
    return resolve_existing_file(
        user_value,
        [
            os.path.join(root, base),
            os.path.join(root, "bioinf_scripts", "scripts", base),
            os.path.join(root, "OUT_v44", base),
            resolve_from_manifest(manifest_paths, [rf"/{re.escape(base)}$"]),
        ],
    )


def build_ferrdb_sets(driver_path: str, suppressor_path: str, marker_path: str) -> Dict[str, Set[str]]:
    drivers = read_ferrdb_gene_set(driver_path)
    suppressors = read_ferrdb_gene_set(suppressor_path)
    markers = read_ferrdb_gene_set(marker_path) if marker_path and os.path.exists(marker_path) else set()
    known = set().union(drivers, suppressors, markers)
    return {"drivers": drivers, "suppressors": suppressors, "markers": markers, "known": known}


def ferrdb_category(gene: str, ferrdb_sets: Dict[str, Set[str]]) -> str:
    if gene in ferrdb_sets["drivers"]:
        return "Driver"
    if gene in ferrdb_sets["suppressors"]:
        return "Suppressor"
    if gene in ferrdb_sets["markers"]:
        return "Marker"
    if gene in ferrdb_sets["known"]:
        return "FerrDb-related"
    return "Other"


# -----------------------------------------------------------------------------
# Input readers
# -----------------------------------------------------------------------------
def read_gene_table(path: str, ferrdb_sets: Dict[str, Set[str]]) -> pd.DataFrame:
    df = safe_read_table(path)
    if df.empty:
        raise FileNotFoundError(f"Gene file not found or empty: {path}")

    cols_low = {str(c).lower(): c for c in df.columns}
    gene_col = next((cols_low[k] for k in ["gene", "symbol", "gene_symbol", "genesymbol", "feature"] if k in cols_low), df.columns[0])
    z_col = next((cols_low[k] for k in ["z_meta", "z", "meta_z"] if k in cols_low), None)
    bag_col = next((cols_low[k] for k in ["bag_frac", "bag_freq", "freq"] if k in cols_low), None)
    p_col = next((cols_low[k] for k in ["p_meta", "pvalue", "p_val"] if k in cols_low), None)
    q_col = next((cols_low[k] for k in ["q_fdr", "fdr", "qvalue"] if k in cols_low), None)

    out = pd.DataFrame({
        "gene": [sanitize_gene(x) for x in df[gene_col].astype(str)],
        "z_meta": pd.to_numeric(df[z_col], errors="coerce") if z_col else 0.0,
        "bag_frac": pd.to_numeric(df[bag_col], errors="coerce") if bag_col else 1.0,
        "p_meta": pd.to_numeric(df[p_col], errors="coerce") if p_col else np.nan,
        "q_fdr": pd.to_numeric(df[q_col], errors="coerce") if q_col else np.nan,
    }).replace([np.inf, -np.inf], np.nan)

    out["z_meta"] = out["z_meta"].fillna(0.0).astype(float)
    out["bag_frac"] = out["bag_frac"].fillna(1.0).astype(float)
    out = out[out["gene"].ne("") & out["gene"].notna()].copy()
    out["direction"] = np.where(out["z_meta"] > 0, "Risk-promoting", "Protective")
    out["abs_z"] = out["z_meta"].abs()
    out["neglog10_q"] = neglog10(out["q_fdr"])
    out["gene_category"] = out["gene"].map(lambda g: ferrdb_category(g, ferrdb_sets))
    out["is_ferrdb_gene"] = out["gene"].isin(ferrdb_sets["known"])
    out["is_driver"] = out["gene"].isin(ferrdb_sets["drivers"])
    out["is_suppressor"] = out["gene"].isin(ferrdb_sets["suppressors"])
    out["is_marker"] = out["gene"].isin(ferrdb_sets["markers"])

    grouped = out.groupby("gene", as_index=False).agg({
        "z_meta": "max",
        "bag_frac": "max",
        "p_meta": "min",
        "q_fdr": "min",
        "neglog10_q": "max",
        "direction": "first",
        "abs_z": "max",
        "gene_category": "first",
        "is_ferrdb_gene": "max",
        "is_driver": "max",
        "is_suppressor": "max",
        "is_marker": "max",
    })
    return grouped.sort_values(["bag_frac", "abs_z"], ascending=[False, False]).reset_index(drop=True)


def read_ppi_table(path: str) -> pd.DataFrame:
    empty = pd.DataFrame(columns=["gene", "degree", "betweenness", "closeness", "eigenvector", "pagerank", "hub_score", "community"])
    if not path or not os.path.exists(path):
        return empty
    df = safe_read_table(path)
    if df.empty:
        return empty

    cols_low = {str(c).lower(): c for c in df.columns}
    gene_col = next((cols_low[k] for k in ["gene", "symbol"] if k in cols_low), df.columns[0])
    degree_col = next((cols_low[k] for k in ["degree", "node_degree"] if k in cols_low), None)
    bet_col = next((cols_low[k] for k in ["betweenness", "betweenness_centrality"] if k in cols_low), None)
    close_col = next((cols_low[k] for k in ["closeness", "closeness_centrality"] if k in cols_low), None)
    eig_col = next((cols_low[k] for k in ["eigenvector", "eigenvector_centrality"] if k in cols_low), None)
    page_col = next((cols_low[k] for k in ["pagerank", "page_rank"] if k in cols_low), None)
    hub_col = next((cols_low[k] for k in ["hub_score", "composite_hub_score", "priority_score"] if k in cols_low), None)
    comm_col = next((cols_low[k] for k in ["community", "module", "cluster"] if k in cols_low), None)

    out = pd.DataFrame({
        "gene": [sanitize_gene(x) for x in df[gene_col].astype(str)],
        "degree": pd.to_numeric(df[degree_col], errors="coerce") if degree_col else np.nan,
        "betweenness": pd.to_numeric(df[bet_col], errors="coerce") if bet_col else np.nan,
        "closeness": pd.to_numeric(df[close_col], errors="coerce") if close_col else np.nan,
        "eigenvector": pd.to_numeric(df[eig_col], errors="coerce") if eig_col else np.nan,
        "pagerank": pd.to_numeric(df[page_col], errors="coerce") if page_col else np.nan,
        "hub_score": pd.to_numeric(df[hub_col], errors="coerce") if hub_col else np.nan,
        "community": df[comm_col].astype(str) if comm_col else np.nan,
    })
    numeric_cols = ["degree", "betweenness", "closeness", "eigenvector", "pagerank", "hub_score"]
    out[numeric_cols] = out[numeric_cols].apply(pd.to_numeric, errors="coerce")
    out = out.groupby("gene", as_index=False).agg({
        "degree": "max",
        "betweenness": "max",
        "closeness": "max",
        "eigenvector": "max",
        "pagerank": "max",
        "hub_score": "max",
        "community": "first",
    })
    return out


# -----------------------------------------------------------------------------
# External API helpers
# -----------------------------------------------------------------------------
def safe_get_json(
    session,
    url: str,
    timeout: int,
    cache_dir: Path,
    cache_prefix: str,
    cache_key: str,
    params: Optional[dict] = None,
) -> dict:
    cached = read_json_cache(cache_dir, cache_prefix, cache_key)
    if cached is not None:
        return cached
    try:
        res = session.get(url, timeout=timeout, params=params)
        payload = res.json()
        write_json_cache(cache_dir, cache_prefix, cache_key, payload)
        return payload
    except Exception:
        return {}


def safe_post_json(
    session,
    url: str,
    timeout: int,
    cache_dir: Path,
    cache_prefix: str,
    cache_key: str,
    payload: dict,
) -> dict:
    cached = read_json_cache(cache_dir, cache_prefix, cache_key)
    if cached is not None:
        return cached
    try:
        res = session.post(url, timeout=timeout, json=payload)
        js = res.json()
        write_json_cache(cache_dir, cache_prefix, cache_key, js)
        return js
    except Exception:
        return {}


def query_dgidb(session, symbol: str, timeout: int, cache_dir: Path) -> List[Dict[str, object]]:
    js = safe_get_json(
        session=session,
        url="https://dgidb.org/api/v2/interactions.json",
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="dgidb",
        cache_key=symbol,
        params={"genes": symbol},
    )
    matched = js.get("matchedTerms", []) or []
    rows: List[Dict[str, object]] = []
    for mt in matched:
        for it in mt.get("interactions", []) or []:
            rows.append({
                "gene": symbol,
                "source": "DGIdb",
                "drug_name": it.get("drugName") or it.get("drug_name"),
                "interaction_type": ", ".join(it.get("interactionTypes", []) or []),
                "score": pd.to_numeric(it.get("score"), errors="coerce"),
                "pmids": ";".join(str(x) for x in (it.get("pmids") or [])),
                "source_dbs": ";".join(str(x) for x in (it.get("sources") or [])),
                "mechanism_of_action": np.nan,
                "max_phase": np.nan,
                "target_class": np.nan,
                "approval_status": np.nan,
                "source_priority": 0.70,
                "evidence_kind": "interaction_db",
            })
    return rows


def query_chembl_target(session, symbol: str, timeout: int, cache_dir: Path) -> Optional[dict]:
    js = safe_get_json(
        session=session,
        url=f"{CHEMBL_BASE}/target/search.json",
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="chembl_target_search",
        cache_key=symbol,
        params={"q": symbol},
    )
    targets = js.get("targets", []) or []
    chosen = None
    for t in targets:
        pref = str(t.get("pref_name", ""))
        if symbol.upper() == pref.upper():
            chosen = t
            break
    if chosen is None and targets:
        chosen = targets[0]
    return chosen


def query_chembl_mechanisms(session, mol_id: str, timeout: int, cache_dir: Path) -> List[dict]:
    js = safe_get_json(
        session=session,
        url=f"{CHEMBL_BASE}/mechanism.json",
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="chembl_mechanism",
        cache_key=mol_id,
        params={"molecule_chembl_id": mol_id, "limit": 100},
    )
    return js.get("mechanisms", []) or []


def query_chembl_molecule(session, mol_id: str, timeout: int, cache_dir: Path) -> dict:
    return safe_get_json(
        session=session,
        url=f"{CHEMBL_BASE}/molecule/{mol_id}.json",
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="chembl_molecule",
        cache_key=mol_id,
    )


def query_chembl(
    session,
    symbol: str,
    timeout: int,
    cache_dir: Path,
    max_activities: int,
    max_phase_queries: int,
    min_pchembl_for_phase_support: float,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    chosen = query_chembl_target(session, symbol, timeout, cache_dir)
    if chosen is None:
        return [], {}

    chembl_target = chosen.get("target_chembl_id")
    target_meta = {
        "gene": symbol,
        "chembl_target_id": chembl_target,
        "chembl_pref_name": chosen.get("pref_name"),
        "target_type": chosen.get("target_type"),
        "organism": chosen.get("organism"),
    }
    if not chembl_target:
        return [], target_meta

    act_js = safe_get_json(
        session=session,
        url=f"{CHEMBL_BASE}/activity.json",
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="chembl_activity",
        cache_key=f"{symbol}_{max_activities}",
        params={"target_chembl_id": chembl_target, "limit": max_activities},
    )
    acts = act_js.get("activities", []) or []
    if not acts:
        return [], target_meta

    activity_rows: List[Dict[str, object]] = []
    mol_ids: List[str] = []

    for a in acts:
        mol = a.get("molecule_chembl_id")
        if not mol:
            continue
        pchembl = pd.to_numeric(a.get("pchembl_value"), errors="coerce")
        std_value = pd.to_numeric(a.get("standard_value"), errors="coerce")
        std_units = a.get("standard_units")
        activity_rows.append({
            "gene": symbol,
            "source": "ChEMBL",
            "drug_name": a.get("molecule_pref_name") or mol,
            "drug_id": mol,
            "interaction_type": a.get("standard_type"),
            "score": pchembl,
            "activity_value": std_value,
            "activity_units": std_units,
            "assay_type": a.get("assay_type"),
            "pmids": a.get("document_chembl_id"),
            "source_dbs": chembl_target,
            "mechanism_of_action": np.nan,
            "max_phase": np.nan,
            "target_class": chosen.get("target_type"),
            "approval_status": np.nan,
            "source_priority": 0.85,
            "evidence_kind": "activity",
        })
        if pd.notna(pchembl) and float(pchembl) >= min_pchembl_for_phase_support:
            mol_ids.append(mol)

    mol_ids = list(dict.fromkeys(mol_ids))[:max_phase_queries]

    phase_cache: Dict[str, float] = {}
    pref_name_cache: Dict[str, str] = {}
    moa_cache: Dict[str, str] = {}
    approval_cache: Dict[str, str] = {}

    for mol in mol_ids:
        mol_js = query_chembl_molecule(session, mol, timeout, cache_dir)
        phase_cache[mol] = pd.to_numeric(mol_js.get("max_phase"), errors="coerce")
        pref_name_cache[mol] = mol_js.get("pref_name") or mol
        approval_cache[mol] = (
            "Approved_or_clinical" if pd.to_numeric(mol_js.get("max_phase"), errors="coerce") >= 1 else "Preclinical_or_unknown"
        )
        mechs = query_chembl_mechanisms(session, mol, timeout, cache_dir)
        if mechs:
            moa_cache[mol] = "; ".join(sorted({
                str(m.get("mechanism_of_action")).strip()
                for m in mechs if m.get("mechanism_of_action")
            }))
        else:
            moa_cache[mol] = ""

    for row in activity_rows:
        mol = row.get("drug_id")
        row["drug_name"] = pref_name_cache.get(mol, row["drug_name"])
        row["max_phase"] = phase_cache.get(mol, np.nan)
        row["mechanism_of_action"] = moa_cache.get(mol, np.nan)
        row["approval_status"] = approval_cache.get(mol, np.nan)

    return activity_rows, target_meta


def query_opentargets_search(session, symbol: str, timeout: int, cache_dir: Path) -> List[dict]:
    payload = {
        "query": """
        query SearchTarget($q: String!) {
          search(queryString: $q, entityNames: ["target"]) {
            hits { id entity name description }
          }
        }
        """,
        "variables": {"q": symbol},
    }
    js = safe_post_json(
        session=session,
        url=OPEN_TARGETS_GRAPHQL,
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="ot_search",
        cache_key=symbol,
        payload=payload,
    )
    return (((js.get("data") or {}).get("search") or {}).get("hits") or [])


def query_opentargets_known_drugs(
    session,
    symbol: str,
    timeout: int,
    cache_dir: Path,
    max_known_drugs: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    hits = query_opentargets_search(session, symbol, timeout, cache_dir)
    target_id = None
    description = None
    for h in hits:
        if str(h.get("name", "")).upper() == symbol.upper():
            target_id = h.get("id")
            description = h.get("description")
            break
    if target_id is None and hits:
        target_id = hits[0].get("id")
        description = hits[0].get("description")
    if target_id is None:
        return [], {}

    payload = {
        "query": """
        query TargetKnownDrugs($ensemblId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            approvedName
            biotype
            tractability {
              label
              modality
              value
            }
            knownDrugs(size: $size) {
              rows {
                drug {
                  id
                  name
                  drugType
                  maximumClinicalTrialPhase
                  mechanismsOfAction {
                    rows {
                      actionType
                      mechanismOfAction
                      targets {
                        approvedSymbol
                        id
                      }
                    }
                  }
                }
                phase
                status
                disease {
                  id
                  name
                }
                targetClass {
                  id
                  label
                }
              }
            }
          }
        }
        """,
        "variables": {"ensemblId": target_id, "size": int(max_known_drugs)},
    }

    js = safe_post_json(
        session=session,
        url=OPEN_TARGETS_GRAPHQL,
        timeout=timeout,
        cache_dir=cache_dir,
        cache_prefix="ot_known_drugs",
        cache_key=f"{symbol}_{max_known_drugs}",
        payload=payload,
    )
    tgt = ((js.get("data") or {}).get("target") or {})
    rows = (((tgt.get("knownDrugs") or {}).get("rows")) or [])
    tract = tgt.get("tractability") or []

    out: List[Dict[str, object]] = []
    for r in rows:
        drug = r.get("drug") or {}
        moa_rows = (((drug.get("mechanismsOfAction") or {}).get("rows")) or [])
        moa_strings = sorted({
            str(m.get("mechanismOfAction")).strip()
            for m in moa_rows if m.get("mechanismOfAction")
        })
        action_types = sorted({
            str(m.get("actionType")).strip()
            for m in moa_rows if m.get("actionType")
        })
        tc = r.get("targetClass") or {}
        out.append({
            "gene": symbol,
            "source": "OpenTargetsKnownDrugs",
            "drug_name": drug.get("name"),
            "drug_id": drug.get("id"),
            "interaction_type": "; ".join(action_types) if action_types else np.nan,
            "score": np.nan,
            "activity_value": np.nan,
            "activity_units": np.nan,
            "assay_type": np.nan,
            "pmids": np.nan,
            "source_dbs": target_id,
            "mechanism_of_action": "; ".join(moa_strings) if moa_strings else np.nan,
            "max_phase": pd.to_numeric(drug.get("maximumClinicalTrialPhase"), errors="coerce"),
            "target_class": tc.get("label"),
            "approval_status": r.get("status"),
            "source_priority": 0.95,
            "evidence_kind": "known_drug",
            "disease_context": ((r.get("disease") or {}).get("name")),
        })

    tract_rows = []
    for t in tract:
        tract_rows.append({
            "gene": symbol,
            "target_id": target_id,
            "approved_symbol": tgt.get("approvedSymbol"),
            "approved_name": tgt.get("approvedName"),
            "biotype": tgt.get("biotype"),
            "tractability_label": t.get("label"),
            "tractability_modality": t.get("modality"),
            "tractability_value": pd.to_numeric(t.get("value"), errors="coerce"),
            "description": description,
        })
    tract_df = pd.DataFrame(tract_rows)
    tract_meta = tract_df.to_dict(orient="list") if not tract_df.empty else {
        "gene": [symbol],
        "target_id": [target_id],
        "approved_symbol": [tgt.get("approvedSymbol")],
        "approved_name": [tgt.get("approvedName")],
        "biotype": [tgt.get("biotype")],
        "description": [description],
    }
    return out, tract_meta


# -----------------------------------------------------------------------------
# Normalization and aggregation
# -----------------------------------------------------------------------------
LINK_BASE_COLS = [
    "gene", "source", "drug_name", "drug_id", "interaction_type", "score",
    "activity_value", "activity_units", "assay_type", "pmids", "source_dbs",
    "mechanism_of_action", "max_phase", "target_class", "approval_status",
    "source_priority", "evidence_kind", "disease_context"
]


def normalize_links(links: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        return pd.DataFrame(columns=LINK_BASE_COLS)
    links = links.copy()
    for c in LINK_BASE_COLS:
        if c not in links.columns:
            links[c] = np.nan
    links["gene"] = links["gene"].astype(str).map(sanitize_gene)
    for c in ["drug_name", "drug_id", "interaction_type", "pmids", "source_dbs", "mechanism_of_action", "target_class", "approval_status", "evidence_kind", "disease_context"]:
        links[c] = links[c].astype(str).replace({"nan": np.nan, "None": np.nan, "": np.nan})
    links["source"] = links["source"].astype(str)
    links["score"] = pd.to_numeric(links["score"], errors="coerce")
    links["activity_value"] = pd.to_numeric(links["activity_value"], errors="coerce")
    links["max_phase"] = pd.to_numeric(links["max_phase"], errors="coerce")
    links["source_priority"] = pd.to_numeric(links["source_priority"], errors="coerce").fillna(0.0)
    links = links.drop_duplicates().sort_values(["gene", "source", "drug_name"], na_position="last").reset_index(drop=True)
    return links[LINK_BASE_COLS]


def normalize_target_metadata(chembl_meta_rows: List[dict], tractability_meta_rows: List[dict]) -> pd.DataFrame:
    rows = []
    for x in chembl_meta_rows:
        rows.append({
            "gene": sanitize_gene(x.get("gene", "")),
            "chembl_target_id": x.get("chembl_target_id"),
            "chembl_pref_name": x.get("chembl_pref_name"),
            "chembl_target_type": x.get("target_type"),
            "chembl_organism": x.get("organism"),
        })
    tract_rows = []
    for x in tractability_meta_rows:
        if isinstance(x, dict) and all(isinstance(v, list) for v in x.values()):
            df = pd.DataFrame(x)
            tract_rows.append(df)
        elif isinstance(x, dict):
            tract_rows.append(pd.DataFrame([x]))
    tract_df = pd.concat(tract_rows, ignore_index=True) if tract_rows else pd.DataFrame()
    tract_expected = [
        "gene", "target_id", "approved_symbol", "approved_name", "biotype", "description",
        "tractability_label", "tractability_modality", "tractability_value"
    ]

    if not tract_df.empty:
        for col in tract_expected:
            if col not in tract_df.columns:
                tract_df[col] = np.nan
        tract_df["gene"] = tract_df["gene"].astype(str).map(sanitize_gene)
        tract_df["tractability_value"] = pd.to_numeric(tract_df["tractability_value"], errors="coerce")
        tract_agg = tract_df.groupby("gene", as_index=False).agg({
            "target_id": "first",
            "approved_symbol": "first",
            "approved_name": "first",
            "biotype": "first",
            "description": "first",
            "tractability_label": lambda x: "; ".join(sorted({str(v) for v in x if pd.notna(v) and str(v).strip() not in {"", "nan", "None"}})),
            "tractability_modality": lambda x: "; ".join(sorted({str(v) for v in x if pd.notna(v) and str(v).strip() not in {"", "nan", "None"}})),
            "tractability_value": "max",
        })
    else:
        tract_agg = pd.DataFrame(columns=tract_expected)

    chembl_df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "gene", "chembl_target_id", "chembl_pref_name", "chembl_target_type", "chembl_organism"
    ])
    if not chembl_df.empty:
        chembl_df = chembl_df.groupby("gene", as_index=False).agg({
            "chembl_target_id": "first",
            "chembl_pref_name": "first",
            "chembl_target_type": "first",
            "chembl_organism": "first",
        })

    if chembl_df.empty and tract_agg.empty:
        return pd.DataFrame(columns=["gene"])
    if chembl_df.empty:
        return tract_agg
    if tract_agg.empty:
        return chembl_df
    return chembl_df.merge(tract_agg, on="gene", how="outer")


def summarise_top_compounds_per_gene(links: pd.DataFrame, top_k: int = 8) -> pd.DataFrame:
    if links.empty:
        return pd.DataFrame(columns=["gene", "top_compounds"])
    tmp = links.copy()
    tmp = tmp[tmp["drug_name"].notna()].copy()
    if tmp.empty:
        return pd.DataFrame(columns=["gene", "top_compounds"])
    tmp["compound_score_proxy"] = (
        0.45 * scaled_series(tmp["max_phase"].fillna(0)) +
        0.35 * scaled_series(tmp["score"].fillna(0)) +
        0.20 * scaled_series(tmp["source_priority"].fillna(0))
    )
    agg = (tmp.groupby(["gene", "drug_name"], as_index=False)
           .agg(best_phase=("max_phase", "max"),
                best_score=("score", "max"),
                best_source_priority=("source_priority", "max"),
                n_sources=("source", "nunique")))
    agg["compound_rank_score"] = (
        0.40 * scaled_series(agg["best_phase"].fillna(0)) +
        0.30 * scaled_series(agg["best_score"].fillna(0)) +
        0.20 * scaled_series(agg["n_sources"].fillna(0)) +
        0.10 * scaled_series(agg["best_source_priority"].fillna(0))
    )
    agg = agg.sort_values(["gene", "compound_rank_score", "best_phase", "best_score"], ascending=[True, False, False, False])
    top = agg.groupby("gene").head(top_k)
    out = (top.groupby("gene", as_index=False)
           .agg(top_compounds=("drug_name", lambda x: "; ".join(list(x)))))
    return out


def build_summary(genes: pd.DataFrame, ppi_df: pd.DataFrame, links: pd.DataFrame, target_meta: pd.DataFrame) -> pd.DataFrame:
    out = genes.merge(ppi_df, on="gene", how="left").merge(target_meta, on="gene", how="left")

    if links.empty:
        links = pd.DataFrame(columns=LINK_BASE_COLS)

    out["n_compounds"] = out["gene"].map(links.groupby("gene")["drug_name"].nunique(dropna=True)).fillna(0).astype(int)
    out["n_records"] = out["gene"].map(links.groupby("gene").size()).fillna(0).astype(int)
    out["n_sources"] = out["gene"].map(links.groupby("gene")["source"].nunique(dropna=True)).fillna(0).astype(int)

    out["n_dgidb_interactions"] = out["gene"].map(
        links[links["source"].eq("DGIdb")].groupby("gene").size()
    ).fillna(0).astype(int)
    out["n_chembl_activities"] = out["gene"].map(
        links[links["source"].eq("ChEMBL")].groupby("gene").size()
    ).fillna(0).astype(int)
    out["n_ot_known_drugs"] = out["gene"].map(
        links[links["source"].eq("OpenTargetsKnownDrugs")].groupby("gene")["drug_name"].nunique(dropna=True)
    ).fillna(0).astype(int)
    out["n_clinical_compounds"] = out["gene"].map(
        links.loc[links["max_phase"].fillna(-1) >= 1].groupby("gene")["drug_name"].nunique(dropna=True)
    ).fillna(0).astype(int)
    out["n_late_phase_compounds"] = out["gene"].map(
        links.loc[links["max_phase"].fillna(-1) >= 3].groupby("gene")["drug_name"].nunique(dropna=True)
    ).fillna(0).astype(int)
    out["best_chembl_pchembl"] = out["gene"].map(
        links[links["source"].eq("ChEMBL")].groupby("gene")["score"].max()
    )
    out["best_phase"] = out["gene"].map(
        links.groupby("gene")["max_phase"].max()
    )
    out["has_dgidb"] = out["n_dgidb_interactions"] > 0
    out["has_chembl"] = out["n_chembl_activities"] > 0
    out["has_ot_known_drugs"] = out["n_ot_known_drugs"] > 0
    out["has_clinical_compound"] = out["n_clinical_compounds"] > 0
    out["has_late_phase_compound"] = out["n_late_phase_compounds"] > 0
    out["has_multisource_support"] = out["n_sources"] >= 2
    out["best_target_class"] = out["gene"].map(
        links.groupby("gene")["target_class"].agg(lambda x: next((v for v in x if pd.notna(v)), np.nan))
    )
    out["best_approval_status"] = out["gene"].map(
        links.groupby("gene")["approval_status"].agg(lambda x: next((v for v in x if pd.notna(v)), np.nan))
    )

    top_compounds = summarise_top_compounds_per_gene(links, top_k=8)
    out = out.merge(top_compounds, on="gene", how="left")

    for c in ["degree", "betweenness", "closeness", "eigenvector", "pagerank", "hub_score", "best_chembl_pchembl", "best_phase", "tractability_value"]:
        if c not in out.columns:
            out[c] = np.nan
        out[c] = pd.to_numeric(out[c], errors="coerce")

    ferrdb_weight = (
        out["is_driver"].astype(float) * 1.00 +
        out["is_suppressor"].astype(float) * 0.90 +
        out["is_marker"].astype(float) * 0.60
    )
    net_fallback = out["hub_score"].combine_first(out["betweenness"]).combine_first(out["degree"])

    out["score_expr"] = (
        0.35 * scaled_series(out["bag_frac"]) +
        0.40 * scaled_series(out["abs_z"]) +
        0.25 * scaled_series(out["neglog10_q"].fillna(0))
    )
    out["score_network"] = (
        0.30 * scaled_series(out["degree"]) +
        0.20 * scaled_series(out["betweenness"]) +
        0.10 * scaled_series(out["closeness"]) +
        0.15 * scaled_series(out["eigenvector"]) +
        0.25 * scaled_series(net_fallback)
    )
    out["score_compound"] = (
        0.22 * scaled_series(out["n_compounds"]) +
        0.16 * scaled_series(out["n_dgidb_interactions"]) +
        0.20 * scaled_series(out["n_chembl_activities"]) +
        0.16 * scaled_series(out["n_ot_known_drugs"]) +
        0.16 * scaled_series(out["best_chembl_pchembl"].fillna(0.0)) +
        0.10 * scaled_series(out["n_sources"])
    )
    out["score_clinical"] = (
        0.35 * scaled_series(out["n_clinical_compounds"]) +
        0.35 * scaled_series(out["best_phase"].fillna(0.0)) +
        0.30 * scaled_series(out["n_late_phase_compounds"])
    )
    out["score_ferrdb"] = np.clip(ferrdb_weight, 0, 1)
    out["score_tractability"] = (
        0.60 * scaled_series(out["tractability_value"].fillna(0)) +
        0.40 * scaled_series(out["n_ot_known_drugs"].fillna(0))
    )
    out["score_multisource"] = (
        0.60 * scaled_series(out["n_sources"].fillna(0)) +
        0.40 * np.where(out["has_multisource_support"], 1.0, 0.0)
    )

    out["priority_score"] = (
        0.20 * out["score_expr"] +
        0.17 * out["score_network"] +
        0.18 * out["score_compound"] +
        0.15 * out["score_clinical"] +
        0.15 * out["score_ferrdb"] +
        0.08 * out["score_tractability"] +
        0.07 * out["score_multisource"]
    )

    out["priority_tier"] = pd.cut(
        out["priority_score"],
        bins=[-np.inf, 0.25, 0.50, 0.75, np.inf],
        labels=["Exploratory", "Moderate", "High", "Very High"],
    ).astype(str)

    out["therapeutic_hypothesis"] = np.where(
        out["direction"].eq("Risk-promoting"),
        "Prefer inhibition / suppression of target activity",
        "Prefer activation / restoration or context-aware modulation",
    )

    evidence_parts = []
    evidence_parts.append(np.where(out["is_driver"], "FerrDb driver", ""))
    evidence_parts.append(np.where(out["is_suppressor"], "FerrDb suppressor", ""))
    evidence_parts.append(np.where(out["is_marker"], "FerrDb marker", ""))
    evidence_parts.append(np.where(out["has_dgidb"], "DGIdb interaction", ""))
    evidence_parts.append(np.where(out["has_chembl"], "ChEMBL activity", ""))
    evidence_parts.append(np.where(out["has_ot_known_drugs"], "OpenTargets known drug", ""))
    evidence_parts.append(np.where(out["has_late_phase_compound"], "late-phase compound", ""))
    evidence_parts.append(np.where(out["has_multisource_support"], "multisource support", ""))
    evidence_parts.append(np.where(out["degree"].fillna(0) > out["degree"].fillna(0).median(), "network-central", ""))
    evidence_df = pd.DataFrame({f"e{i}": v for i, v in enumerate(evidence_parts, start=1)})
    out["evidence_summary"] = evidence_df.apply(lambda r: "; ".join([x for x in r if x]), axis=1)

    out = out.sort_values(
        ["priority_score", "score_clinical", "score_compound", "bag_frac", "abs_z"],
        ascending=[False, False, False, False, False]
    ).reset_index(drop=True)
    return out


def build_actionable_table(summary: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "gene", "gene_category", "direction", "therapeutic_hypothesis", "z_meta", "bag_frac", "p_meta", "q_fdr",
        "degree", "betweenness", "closeness", "eigenvector", "pagerank", "hub_score", "community",
        "chembl_target_id", "chembl_target_type", "target_id", "approved_name", "biotype",
        "tractability_label", "tractability_modality", "tractability_value",
        "n_compounds", "n_records", "n_sources", "n_dgidb_interactions", "n_chembl_activities",
        "n_ot_known_drugs", "n_clinical_compounds", "n_late_phase_compounds",
        "best_chembl_pchembl", "best_phase", "best_target_class", "best_approval_status",
        "has_dgidb", "has_chembl", "has_ot_known_drugs", "has_clinical_compound", "has_late_phase_compound",
        "has_multisource_support", "is_ferrdb_gene", "is_driver", "is_suppressor", "is_marker",
        "score_expr", "score_network", "score_compound", "score_clinical", "score_ferrdb",
        "score_tractability", "score_multisource", "priority_score", "priority_tier",
        "top_compounds", "evidence_summary"
    ]
    out = summary.copy()
    for c in keep:
        if c not in out.columns:
            out[c] = np.nan
    return out[keep].copy()


def _safe_three_level_bin(
    values: pd.Series,
    low_threshold: float,
    mid_threshold: float,
    labels: Sequence[str],
    right: bool = True,
) -> pd.Series:
    """Robust three-level binning that avoids duplicated bin edges."""
    vals = pd.to_numeric(values, errors="coerce")
    low_threshold = float(low_threshold)
    mid_threshold = float(mid_threshold)

    if not np.isfinite(mid_threshold) or mid_threshold <= low_threshold:
        out = pd.Series(labels[2], index=vals.index, dtype="object")
        out = out.mask(vals.isna(), np.nan)
        out = out.mask(vals <= low_threshold, labels[0])
        out = out.mask((vals > low_threshold) & vals.notna(), labels[2])
        return out.astype(str)

    return pd.cut(
        vals,
        bins=[-np.inf, low_threshold, mid_threshold, np.inf],
        labels=list(labels),
        right=right,
        duplicates="drop",
    ).astype(str)


def build_evidence_matrix(summary: pd.DataFrame) -> pd.DataFrame:
    ev = summary[[
        "gene", "gene_category", "direction", "bag_frac", "abs_z", "neglog10_q", "degree", "hub_score",
        "n_dgidb_interactions", "n_chembl_activities", "n_ot_known_drugs", "n_clinical_compounds", "best_phase",
        "tractability_value", "is_driver", "is_suppressor", "is_marker", "priority_score"
    ]].copy()
    ev["expr_stability_bin"] = pd.cut(
        ev["bag_frac"],
        bins=[-np.inf, 0.74, 0.99, np.inf],
        labels=["Variable", "Stable", "Core"],
    ).astype(str)

    degree_nonneg = pd.to_numeric(ev["degree"], errors="coerce").fillna(0)
    degree_median = float(degree_nonneg.median()) if len(degree_nonneg) else 0.0
    ev["network_bin"] = _safe_three_level_bin(
        degree_nonneg,
        low_threshold=0.0,
        mid_threshold=degree_median,
        labels=["Low", "Intermediate", "High"],
    )

    ev["clinical_bin"] = pd.cut(
        ev["n_clinical_compounds"],
        bins=[-np.inf, 0, 1, np.inf],
        labels=["None", "Sparse", "Rich"],
    ).astype(str)
    return ev


def build_source_summary(links: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        return pd.DataFrame(columns=["source", "n_rows", "n_genes", "n_compounds", "median_score", "median_phase", "max_phase"])
    tmp = links.copy()
    return (tmp.groupby("source", as_index=False)
            .agg(
                n_rows=("gene", "size"),
                n_genes=("gene", "nunique"),
                n_compounds=("drug_name", lambda x: pd.Series(x).nunique(dropna=True)),
                median_score=("score", "median"),
                median_phase=("max_phase", "median"),
                max_phase=("max_phase", "max"),
            )
            .sort_values(["n_genes", "n_rows"], ascending=False))


def build_mechanism_summary(links: pd.DataFrame) -> pd.DataFrame:
    if links.empty:
        return pd.DataFrame(columns=["interaction_type", "n_rows", "n_genes", "n_compounds"])
    tmp = links.copy()
    tmp["interaction_type"] = tmp["interaction_type"].fillna("Unspecified").astype(str).str.strip()
    return (tmp.groupby("interaction_type", as_index=False)
            .agg(
                n_rows=("gene", "size"),
                n_genes=("gene", "nunique"),
                n_compounds=("drug_name", lambda x: pd.Series(x).nunique(dropna=True)),
            )
            .sort_values(["n_rows", "n_genes"], ascending=False))


def build_ferrdb_overlap_summary(summary: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "category": ["FerrDb overall", "FerrDb drivers", "FerrDb suppressors", "FerrDb markers"],
        "n_targets": [
            int(summary["is_ferrdb_gene"].sum()),
            int(summary["is_driver"].sum()),
            int(summary["is_suppressor"].sum()),
            int(summary["is_marker"].sum()),
        ],
        "median_priority_score": [
            float(summary.loc[summary["is_ferrdb_gene"], "priority_score"].median()) if summary["is_ferrdb_gene"].any() else np.nan,
            float(summary.loc[summary["is_driver"], "priority_score"].median()) if summary["is_driver"].any() else np.nan,
            float(summary.loc[summary["is_suppressor"], "priority_score"].median()) if summary["is_suppressor"].any() else np.nan,
            float(summary.loc[summary["is_marker"], "priority_score"].median()) if summary["is_marker"].any() else np.nan,
        ]
    })


def build_tier_breakdown(summary: pd.DataFrame) -> pd.DataFrame:
    return (summary.groupby("priority_tier", as_index=False)
            .agg(
                n_targets=("gene", "nunique"),
                median_priority=("priority_score", "median"),
                median_clinical_compounds=("n_clinical_compounds", "median"),
                median_best_phase=("best_phase", "median"),
                median_tractability=("tractability_value", "median"),
                n_ferrdb=("is_ferrdb_gene", "sum"),
                n_drivers=("is_driver", "sum"),
                n_suppressors=("is_suppressor", "sum"),
            )
            .sort_values("median_priority", ascending=False))


def build_compound_candidate_table(summary: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    if links.empty or summary.empty:
        return pd.DataFrame(columns=[
            "drug_name", "genes_supported", "n_genes_supported", "n_sources", "best_phase",
            "best_pchembl", "target_classes", "mechanism_of_action", "approval_statuses",
            "compound_priority_score", "supported_genes", "gene_directions", "rationale"
        ])

    tmp = links.copy()
    tmp = tmp[tmp["drug_name"].notna()].copy()
    if tmp.empty:
        return pd.DataFrame(columns=[
            "drug_name", "genes_supported", "n_genes_supported", "n_sources", "best_phase",
            "best_pchembl", "target_classes", "mechanism_of_action", "approval_statuses",
            "compound_priority_score", "supported_genes", "gene_directions", "rationale"
        ])

    gene_priority = summary.set_index("gene")["priority_score"].to_dict()
    gene_direction = summary.set_index("gene")["direction"].to_dict()

    agg = (tmp.groupby("drug_name", as_index=False)
           .agg(
               n_link_rows=("gene", "size"),
               n_genes_supported=("gene", "nunique"),
               n_sources=("source", "nunique"),
               best_phase=("max_phase", "max"),
               best_pchembl=("score", "max"),
               supported_genes=("gene", lambda x: "; ".join(sorted(set([str(v) for v in x if pd.notna(v)])))),
               gene_directions=("gene", lambda x: "; ".join(sorted({f"{g}:{gene_direction.get(g, 'NA')}" for g in set(x)}))),
               genes_supported=("gene", lambda x: "; ".join(sorted(set([str(v) for v in x if pd.notna(v)])))),
               target_classes=("target_class", lambda x: "; ".join(sorted(set([str(v) for v in x if pd.notna(v)])))),
               mechanism_of_action=("mechanism_of_action", lambda x: "; ".join(sorted(set([str(v) for v in x if pd.notna(v)]))[:8])),
               approval_statuses=("approval_status", lambda x: "; ".join(sorted(set([str(v) for v in x if pd.notna(v)])))),
           ))

    agg["mean_gene_priority"] = agg["genes_supported"].map(
        lambda s: np.mean([gene_priority.get(g, np.nan) for g in str(s).split("; ") if g]) if str(s).strip() else np.nan
    )
    agg["compound_priority_score"] = (
        0.30 * scaled_series(agg["n_genes_supported"]) +
        0.20 * scaled_series(agg["n_sources"]) +
        0.20 * scaled_series(agg["best_phase"].fillna(0)) +
        0.15 * scaled_series(agg["best_pchembl"].fillna(0)) +
        0.15 * scaled_series(agg["mean_gene_priority"].fillna(0))
    )
    agg["rationale"] = agg.apply(
        lambda r: (
            f"Supports {int(r['n_genes_supported'])} ML-derived target(s); "
            f"best phase={int(r['best_phase']) if pd.notna(r['best_phase']) else 'NA'}; "
            f"n_sources={int(r['n_sources'])}; "
            f"best_pchembl={f'{r['best_pchembl']:.2f}' if pd.notna(r['best_pchembl']) else 'NA'}"
        ),
        axis=1,
    )
    agg = agg.sort_values(
        ["compound_priority_score", "n_genes_supported", "best_phase", "best_pchembl"],
        ascending=[False, False, False, False]
    ).reset_index(drop=True)
    return agg


def build_gene_compound_network_edges(summary: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    if links.empty or summary.empty:
        return pd.DataFrame(columns=[
            "gene", "drug_name", "source", "interaction_type", "score", "max_phase",
            "gene_priority_score", "gene_direction", "edge_support_score"
        ])
    gene_priority = summary.set_index("gene")["priority_score"].to_dict()
    gene_direction = summary.set_index("gene")["direction"].to_dict()

    edges = links.copy()
    edges = edges[edges["drug_name"].notna()].copy()
    edges["gene_priority_score"] = edges["gene"].map(gene_priority)
    edges["gene_direction"] = edges["gene"].map(gene_direction)
    edges["edge_support_score"] = (
        0.35 * scaled_series(edges["gene_priority_score"].fillna(0)) +
        0.25 * scaled_series(edges["max_phase"].fillna(0)) +
        0.20 * scaled_series(edges["score"].fillna(0)) +
        0.20 * scaled_series(edges["source_priority"].fillna(0))
    )
    return edges.sort_values("edge_support_score", ascending=False).reset_index(drop=True)


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
FIGURE_DPI = 600
PAPER_COLORS = {
    "ink": "#1F2937",
    "muted": "#667085",
    "grid": "#E4E7EC",
    "risk": "#C4321C",
    "protective": "#1D4ED8",
    "neutral": "#98A2B3",
    "accent": "#7C3AED",
    "success": "#039855",
    "teal": "#0F766E",
    "gold": "#B7791F",
    "rose": "#C11574",
    "background": "#FFFFFF",
}


def _paper_seq_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        "paper_seq", ["#F8FAFC", "#DCE7F7", "#A8C4EA", "#5B8CCB", "#1D3557"]
    )


def _paper_div_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        "paper_div", ["#175CD3", "#F8FAFC", "#B42318"]
    )


def apply_q1_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 16,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.labelsize": 12,
        "axes.labelcolor": PAPER_COLORS["ink"],
        "axes.labelweight": "bold",
        "axes.edgecolor": PAPER_COLORS["muted"],
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": PAPER_COLORS["ink"],
        "ytick.color": PAPER_COLORS["ink"],
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "legend.frameon": False,
        "legend.fontsize": 10,
        "legend.title_fontsize": 11,
        "figure.facecolor": PAPER_COLORS["background"],
        "axes.facecolor": PAPER_COLORS["background"],
        "savefig.facecolor": PAPER_COLORS["background"],
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.dpi": FIGURE_DPI,
        "savefig.dpi": FIGURE_DPI,
    })


def _style_axis(ax, add_x_grid: bool = True, add_y_grid: bool = False) -> None:
    ax.set_axisbelow(True)
    if add_x_grid:
        ax.grid(axis="x", color=PAPER_COLORS["grid"], linewidth=0.85, alpha=0.95)
    if add_y_grid:
        ax.grid(axis="y", color=PAPER_COLORS["grid"], linewidth=0.85, alpha=0.95)
    ax.spines["left"].set_color(PAPER_COLORS["muted"])
    ax.spines["bottom"].set_color(PAPER_COLORS["muted"])
    ax.tick_params(length=4)


def _truncate_label(x: object, max_len: int = 28) -> str:
    s = str(x)
    return s if len(s) <= max_len else s[: max_len - 1] + "…"


def _size_scale(series: pd.Series, min_size: float = 80.0, max_size: float = 900.0) -> np.ndarray:
    vals = scaled_series(series.fillna(0)).to_numpy(dtype=float)
    return min_size + (max_size - min_size) * vals


def _clean_colorbar(cbar, label: str) -> None:
    cbar.set_label(label, fontsize=12, weight="bold", color=PAPER_COLORS["ink"])
    cbar.outline.set_linewidth(0.7)
    cbar.outline.set_edgecolor(PAPER_COLORS["muted"])
    cbar.ax.tick_params(labelsize=10, colors=PAPER_COLORS["ink"])


def _set_xmargin(ax, max_x: float, min_x: float = 0.0, right_extra_frac: float = 0.18, min_pad: float = 1.0) -> Tuple[float, float, float]:
    max_x = float(max_x) if np.isfinite(max_x) else 0.0
    min_x = float(min_x) if np.isfinite(min_x) else 0.0
    x_range = max(max_x - min_x, 1.0)
    pad = max(min_pad, x_range * right_extra_frac)
    ax.set_xlim(min_x - 0.05 * x_range, max_x + pad)
    return min_x, max_x, pad


def _annotate_aligned_value_labels(
    ax,
    x_values: Sequence[float],
    y_values: Sequence[float],
    labels: Sequence[str],
    text_x: Optional[float] = None,
    text_color: str = PAPER_COLORS["muted"],
    line_color: str = PAPER_COLORS["grid"],
) -> None:
    x_arr = np.asarray(pd.to_numeric(pd.Series(list(x_values)), errors="coerce"), dtype=float)
    y_arr = np.asarray(pd.to_numeric(pd.Series(list(y_values)), errors="coerce"), dtype=float)
    labels = [str(v) for v in labels]
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if not np.any(mask):
        return
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]
    labels = [lab for lab, keep in zip(labels, mask) if keep]

    x0, x1 = ax.get_xlim()
    xr = max(x1 - x0, 1.0)
    if text_x is None:
        text_x = x1 - 0.30 * xr

    connector_end = text_x - 0.018 * xr
    for xv, yv, lab in zip(x_arr, y_arr, labels):
        ax.plot([xv + 0.012 * xr, connector_end], [yv, yv], color=line_color, linewidth=0.8, zorder=2)
        txt = ax.text(
            text_x,
            yv,
            lab,
            va="center",
            ha="left",
            fontsize=9.5,
            color=text_color,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.92),
            zorder=4,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])


def _spread_close_points(
    x: Sequence[float],
    y: Sequence[float],
    min_dx: float,
    min_dy: float,
    max_iter: int = 25,
) -> Tuple[np.ndarray, np.ndarray]:
    x_adj = np.asarray(pd.to_numeric(pd.Series(list(x)), errors="coerce"), dtype=float).copy()
    y_adj = np.asarray(pd.to_numeric(pd.Series(list(y)), errors="coerce"), dtype=float).copy()
    if len(x_adj) <= 1:
        return x_adj, y_adj

    for i in range(len(x_adj)):
        if not (np.isfinite(x_adj[i]) and np.isfinite(y_adj[i])):
            continue
        for _ in range(max_iter):
            moved = False
            for j in range(i):
                if not (np.isfinite(x_adj[j]) and np.isfinite(y_adj[j])):
                    continue
                if abs(x_adj[i] - x_adj[j]) < min_dx and abs(y_adj[i] - y_adj[j]) < min_dy:
                    angle = ((i + 1) * 2.399963229728653) % (2 * np.pi)
                    x_adj[i] += 0.55 * min_dx * np.cos(angle)
                    y_adj[i] += 0.55 * min_dy * np.sin(angle)
                    moved = True
            if not moved:
                break
    return x_adj, y_adj


def _annotate_scatter_right_column(
    ax,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    label_col: str,
    n: int = 8,
    rank_by: Optional[str] = None,
    min_sep_frac: float = 0.06,
) -> None:
    if df.empty or any(c not in df.columns for c in [x_col, y_col, label_col]):
        return

    work = df.copy()
    rank_source = rank_by if rank_by in work.columns else x_col
    work["__rank__"] = pd.to_numeric(work[rank_source], errors="coerce")
    work[x_col] = pd.to_numeric(work[x_col], errors="coerce")
    work[y_col] = pd.to_numeric(work[y_col], errors="coerce")
    work = work.dropna(subset=[x_col, y_col, "__rank__"]).sort_values("__rank__", ascending=False).head(min(n, len(work)))
    if work.empty:
        return

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xr = max(x1 - x0, 1e-6)
    yr = max(y1 - y0, 1e-6)

    label_x = x1 - 0.22 * xr
    label_y = work[y_col].to_numpy(dtype=float).copy()
    order = np.argsort(label_y)
    min_sep = min_sep_frac * abs(yr)
    if min_sep <= 0:
        min_sep = 0.03

    for k in range(1, len(order)):
        cur = order[k]
        prev = order[k - 1]
        if label_y[cur] - label_y[prev] < min_sep:
            label_y[cur] = label_y[prev] + min_sep

    upper_bound = y1 - 0.03 * abs(yr)
    lower_bound = y0 + 0.03 * abs(yr)
    if len(order) > 0 and label_y[order[-1]] > upper_bound:
        label_y[order[-1]] = upper_bound
        for k in range(len(order) - 2, -1, -1):
            cur = order[k]
            nxt = order[k + 1]
            label_y[cur] = min(label_y[cur], label_y[nxt] - min_sep)
    if len(order) > 0 and label_y[order[0]] < lower_bound:
        shift = lower_bound - label_y[order[0]]
        label_y = label_y + shift

    for (_, row), ly in zip(work.iterrows(), label_y):
        ax.plot([row[x_col] + 0.01 * xr, label_x - 0.015 * xr], [row[y_col], ly], color=PAPER_COLORS["grid"], linewidth=0.8, zorder=2)
        txt = ax.text(
            label_x,
            ly,
            str(row[label_col]),
            va="center",
            ha="left",
            fontsize=9.5,
            color=PAPER_COLORS["ink"],
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=PAPER_COLORS["grid"], linewidth=0.5, alpha=0.96),
            zorder=5,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])


def savefig_both(fig: plt.Figure, png_path: Path, pdf_path: Path) -> None:
    fig.savefig(png_path, dpi=FIGURE_DPI, metadata={"Creator": "ChatGPT", "Title": png_path.stem})
    fig.savefig(pdf_path, metadata={"Creator": "ChatGPT", "Title": pdf_path.stem})
    plt.close(fig)


def fig_actionability_bubble(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    df = summary.head(min(18, len(summary))).copy()
    df = df.sort_values(["priority_score", "score_clinical", "n_compounds"], ascending=[True, True, True])
    y = np.arange(len(df), dtype=float)

    fig, ax = plt.subplots(figsize=(14.2, 9.8))
    sc = ax.scatter(
        df["n_compounds"].fillna(0),
        y,
        s=_size_scale(df["priority_score"], min_size=140, max_size=930),
        c=df["score_clinical"].fillna(0),
        cmap=_paper_seq_cmap(),
        edgecolors="white",
        linewidths=1.0,
        alpha=0.95,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Linked compounds (n)")
    ax.set_ylabel("")
    ax.set_title("Prioritized actionable targets across the ferroptosis-associated LUAD signature")
    _style_axis(ax, add_x_grid=True, add_y_grid=False)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    _set_xmargin(ax, max_x=float(df["n_compounds"].fillna(0).max()), min_x=float(min(0.0, df["n_compounds"].fillna(0).min())), right_extra_frac=0.24, min_pad=7.0)
    _annotate_aligned_value_labels(ax, df["n_compounds"].fillna(0), y, [f"Priority={v:.2f}" for v in df["priority_score"].fillna(0)])

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    _clean_colorbar(cbar, "Clinical-support score")

    size_handles = [
        plt.scatter([], [], s=s, facecolor="#AFC6E9", edgecolor="white", linewidth=0.9)
        for s in [180, 460, 820]
    ]
    size_labels = ["low", "intermediate", "high"]
    ax.legend(size_handles, size_labels, title="Priority bubble size", loc="lower right")

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T1_target_actionability_bubble.png", outdir / "Fig_T1_target_actionability_bubble.pdf")


def fig_priority_waterfall(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    df = summary.head(min(25, len(summary))).copy().sort_values("priority_score", ascending=True)
    y = np.arange(len(df))
    colors = df["direction"].map({
        "Risk-promoting": PAPER_COLORS["risk"],
        "Protective": PAPER_COLORS["protective"],
    }).fillna(PAPER_COLORS["neutral"])

    fig, ax = plt.subplots(figsize=(14.6, 10.2))
    ax.barh(y, df["priority_score"], color=colors, edgecolor="white", linewidth=0.9, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Integrated priority score")
    ax.set_ylabel("")
    ax.set_title("Integrated therapeutic ranking of signature-associated targets")
    _style_axis(ax, add_x_grid=True, add_y_grid=False)
    ax.set_xlim(0, max(1.02, float(df["priority_score"].max()) + 0.24))

    q75 = float(df["priority_score"].quantile(0.75))
    ax.axvline(q75, color=PAPER_COLORS["accent"], linestyle="--", linewidth=1.2, alpha=0.95)
    ax.text(q75 + 0.012, len(df) - 0.6, "upper quartile", fontsize=9.5, color=PAPER_COLORS["accent"], va="top")

    for yi, val in enumerate(df["priority_score"]):
        ax.text(float(val) + 0.014, yi, f"{val:.2f}", va="center", fontsize=10, color=PAPER_COLORS["ink"])

    legend_handles = [
        Patch(facecolor=PAPER_COLORS["risk"], edgecolor="none", label="Risk-promoting"),
        Patch(facecolor=PAPER_COLORS["protective"], edgecolor="none", label="Protective"),
    ]
    ax.legend(handles=legend_handles, title="Direction", loc="lower right")

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T2_priority_waterfall.png", outdir / "Fig_T2_priority_waterfall.pdf")


def fig_weighted_evidence_heatmap(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    top = summary.head(min(18, len(summary))).copy()
    mat = pd.DataFrame({
        "FerrDb prior": np.clip(top["score_ferrdb"].fillna(0).values, 0, 1),
        "PPI support": np.clip(top["score_network"].fillna(0).values, 0, 1),
        "Compound evidence": np.clip(top["score_compound"].fillna(0).values, 0, 1),
        "Clinical maturity": np.clip(top["score_clinical"].fillna(0).values, 0, 1),
        "Tractability": np.clip(top["score_tractability"].fillna(0).values, 0, 1),
        "Overall priority": np.clip(top["priority_score"].fillna(0).values, 0, 1),
    }, index=top["gene"].tolist())

    fig, ax = plt.subplots(figsize=(11.8, 10.8))
    im = ax.imshow(mat.values, aspect="auto", cmap=_paper_seq_cmap(), norm=Normalize(vmin=0, vmax=1))
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=28, ha="right")
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index)
    ax.set_title("Evidence architecture of the top-ranked actionable targets")

    ax.set_xticks(np.arange(-0.5, mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            value = float(mat.iat[i, j])
            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=9.5,
                color=("white" if value >= 0.58 else PAPER_COLORS["ink"]),
            )

    cbar = fig.colorbar(im, ax=ax, pad=0.025)
    _clean_colorbar(cbar, "Scaled evidence score")
    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T3_weighted_evidence_heatmap.png", outdir / "Fig_T3_weighted_evidence_heatmap.pdf")


def fig_phase_distribution(links: pd.DataFrame, outdir: Path) -> None:
    if links.empty:
        return
    phase_links = links[links["max_phase"].notna()].copy()
    if phase_links.empty:
        fig, ax = plt.subplots(figsize=(8.8, 4.2))
        ax.axis("off")
        ax.text(0.02, 0.62, "No non-missing clinical phase assignments were retrieved.", fontsize=12, weight="bold", transform=ax.transAxes)
        ax.text(0.02, 0.40, "The clinical-development phase panel is therefore omitted.", fontsize=11, color=PAPER_COLORS["muted"], transform=ax.transAxes)
        ax.set_title("Clinical-development phase distribution of linked compounds")
        fig.tight_layout()
        savefig_both(fig, outdir / "Fig_T4_phase_distribution.png", outdir / "Fig_T4_phase_distribution.pdf")
        return

    phase = phase_links["max_phase"].fillna(-1)
    phase_links["phase_group"] = np.select(
        [phase == 0, phase <= 2, phase >= 3],
        ["Preclinical / unknown", "Early clinical (I-II)", "Late clinical (III-IV)"],
        default="Other",
    )
    plot_df = (phase_links.groupby("phase_group", as_index=False)
               .agg(n_records=("gene", "size"),
                    n_unique_compounds=("drug_name", lambda x: pd.Series(x).nunique(dropna=True)))
               .sort_values(["n_records", "n_unique_compounds"], ascending=False))

    palette = {
        "Late clinical (III-IV)": PAPER_COLORS["success"],
        "Early clinical (I-II)": PAPER_COLORS["gold"],
        "Preclinical / unknown": PAPER_COLORS["muted"],
        "Other": PAPER_COLORS["neutral"],
    }

    fig, ax = plt.subplots(figsize=(10.4, 6.4))
    ypos = np.arange(len(plot_df))
    ax.barh(ypos, plot_df["n_records"], color=[palette.get(x, PAPER_COLORS["neutral"]) for x in plot_df["phase_group"]], edgecolor="white", linewidth=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels(plot_df["phase_group"])
    ax.invert_yaxis()
    ax.set_xlabel("Linked records (n)")
    ax.set_ylabel("")
    ax.set_title("Clinical-development phase distribution of linked compounds")
    _style_axis(ax, add_x_grid=True, add_y_grid=False)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    _set_xmargin(ax, max_x=float(plot_df["n_records"].max()), right_extra_frac=0.18, min_pad=1.2)

    for yi, (_, r) in zip(ypos, plot_df.iterrows()):
        ax.text(float(r["n_records"]) + 0.12, yi, f"{int(r['n_unique_compounds'])} unique compounds", va="center", fontsize=9.5, color=PAPER_COLORS["muted"])

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T4_phase_distribution.png", outdir / "Fig_T4_phase_distribution.pdf")


def fig_ferrdb_class_distribution(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    plot_df = (summary.groupby("gene_category", as_index=False)
               .agg(n_targets=("gene", "nunique"), median_priority=("priority_score", "median"))
               .sort_values(["n_targets", "median_priority"], ascending=False))

    fig, ax1 = plt.subplots(figsize=(10.8, 6.8))
    x = np.arange(len(plot_df))
    bars = ax1.bar(x, plot_df["n_targets"], color="#CFE1F2", edgecolor="white", linewidth=0.9, zorder=3)
    ax1.set_xticks(x)
    ax1.set_xticklabels([_truncate_label(v, 20) for v in plot_df["gene_category"]], rotation=20, ha="right")
    ax1.set_ylabel("Number of targets")
    ax1.set_xlabel("")
    ax1.set_title("FerrDb category structure of prioritized LUAD targets")
    _style_axis(ax1, add_x_grid=False, add_y_grid=True)

    ax2 = ax1.twinx()
    ax2.plot(x, plot_df["median_priority"], color=PAPER_COLORS["accent"], marker="o", linewidth=2.0, markersize=5.5)
    ax2.set_ylabel("Median priority score", color=PAPER_COLORS["accent"], weight="bold")
    ax2.tick_params(axis="y", colors=PAPER_COLORS["accent"])
    ax2.spines["right"].set_visible(False)

    for rect, val in zip(bars, plot_df["n_targets"]):
        ax1.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 0.15, f"{int(val)}", ha="center", va="bottom", fontsize=9.5)

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T5_ferrdb_class_distribution.png", outdir / "Fig_T5_ferrdb_class_distribution.pdf")


def fig_ferrdb_priority_composition(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    tmp = (summary.groupby(["priority_tier", "gene_category"], as_index=False).agg(n=("gene", "nunique")))
    tiers = ["Exploratory", "Moderate", "High", "Very High"]
    piv = tmp.pivot(index="priority_tier", columns="gene_category", values="n").reindex(tiers).fillna(0)
    piv = piv.loc[piv.sum(axis=1) > 0]
    if piv.empty:
        return

    cat_order = list(piv.columns)
    palette = ["#DDEAF7", "#B7D0EE", "#7FA9DC", "#4D81C2", "#7A5AF8", "#344054"]
    fig, ax = plt.subplots(figsize=(10.8, 7.2))
    bottom = np.zeros(len(piv), dtype=float)

    for i, col in enumerate(cat_order):
        vals = piv[col].values.astype(float)
        ax.bar(piv.index, vals, bottom=bottom, label=col, color=palette[i % len(palette)], edgecolor="white", linewidth=0.9)
        bottom += vals

    totals = piv.sum(axis=1).values.astype(float)
    for xi, total in enumerate(totals):
        ax.text(xi, total + 0.18, f"n={int(total)}", ha="center", va="bottom", fontsize=9.5, color=PAPER_COLORS["muted"])

    ax.set_ylabel("Number of targets")
    ax.set_xlabel("Priority tier")
    ax.set_title("FerrDb-aware target composition across therapeutic-priority tiers")
    _style_axis(ax, add_x_grid=False, add_y_grid=True)
    ax.legend(title="Gene category", bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T6_ferrdb_priority_composition.png", outdir / "Fig_T6_ferrdb_priority_composition.pdf")


def fig_structured_priority_panel(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    df = summary.head(min(18, len(summary))).copy()
    x_vals = df["score_network"].fillna(0)
    y_vals = df["score_compound"].fillna(0)
    x_spread, y_spread = _spread_close_points(x_vals, y_vals, min_dx=0.028, min_dy=0.032)
    df["__plot_x__"] = x_spread
    df["__plot_y__"] = y_spread

    fig, ax = plt.subplots(figsize=(12.6, 9.8))
    sc = ax.scatter(
        df["__plot_x__"],
        df["__plot_y__"],
        s=_size_scale(df["score_clinical"], min_size=120, max_size=720),
        c=df["score_ferrdb"].fillna(0),
        cmap=_paper_seq_cmap(),
        edgecolors="white",
        linewidths=1.0,
        alpha=0.88,
        zorder=3,
    )

    xmed = float(pd.to_numeric(df["score_network"], errors="coerce").median())
    ymed = float(pd.to_numeric(df["score_compound"], errors="coerce").median())
    ax.axvline(xmed, linestyle="--", linewidth=1.0, color=PAPER_COLORS["muted"], alpha=0.95)
    ax.axhline(ymed, linestyle="--", linewidth=1.0, color=PAPER_COLORS["muted"], alpha=0.95)
    ax.set_xlim(max(-0.02, float(np.nanmin(x_spread)) - 0.04), max(1.02, float(np.nanmax(x_spread)) + 0.28))
    ax.set_ylim(max(-0.02, float(np.nanmin(y_spread)) - 0.04), max(1.02, float(np.nanmax(y_spread)) + 0.06))
    xr = ax.get_xlim()[1] - ax.get_xlim()[0]
    yr = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax.text(xmed + 0.01 * xr, ax.get_ylim()[1] - 0.02 * yr, "median network", fontsize=9, color=PAPER_COLORS["muted"], va="top")
    ax.text(ax.get_xlim()[0] + 0.01 * xr, ymed + 0.012 * yr, "median compound", fontsize=9, color=PAPER_COLORS["muted"], va="bottom")

    label_df = df.copy()
    label_df["__label__"] = label_df["gene"].astype(str) + " (P=" + label_df["priority_score"].map(lambda v: f"{float(v):.2f}" if pd.notna(v) else "NA") + ")"
    _annotate_scatter_right_column(
        ax,
        label_df.sort_values("priority_score", ascending=False),
        "__plot_x__",
        "__plot_y__",
        "__label__",
        n=min(8, len(label_df)),
        rank_by="priority_score",
        min_sep_frac=0.06,
    )

    ax.set_xlabel("Network-support score")
    ax.set_ylabel("Compound-evidence score")
    ax.set_title("Network centrality versus compound support for top-ranked targets")
    _style_axis(ax, add_x_grid=True, add_y_grid=True)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    _clean_colorbar(cbar, "FerrDb-support score")

    size_handles = [
        plt.scatter([], [], s=s, facecolor="#AFC6E9", edgecolor="white", linewidth=0.9)
        for s in [160, 380, 640]
    ]
    ax.legend(size_handles, ["low", "intermediate", "high"], title="Clinical-support size", loc="lower right")

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T7_network_vs_compound_scatter.png", outdir / "Fig_T7_network_vs_compound_scatter.pdf")


def fig_compound_candidate_bubble(compounds: pd.DataFrame, outdir: Path) -> None:
    if compounds.empty:
        return
    df = compounds.head(min(18, len(compounds))).copy()
    df = df.sort_values(["compound_priority_score", "best_phase", "n_genes_supported"], ascending=[True, True, True])
    y = np.arange(len(df), dtype=float)

    fig, ax = plt.subplots(figsize=(14.6, 10.0))
    sc = ax.scatter(
        df["n_genes_supported"].fillna(0),
        y,
        s=_size_scale(df["compound_priority_score"], min_size=130, max_size=820),
        c=df["best_phase"].fillna(0),
        cmap=_paper_seq_cmap(),
        edgecolors="white",
        linewidths=1.0,
        alpha=0.95,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([_truncate_label(x, 36) for x in df["drug_name"]])
    ax.set_xlabel("Supported ML-derived genes (n)")
    ax.set_ylabel("")
    ax.set_title("Candidate compounds supported by the LUAD ferroptosis signature")
    _style_axis(ax, add_x_grid=True, add_y_grid=False)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    _set_xmargin(ax, max_x=float(df["n_genes_supported"].fillna(0).max()), min_x=float(min(0.0, df["n_genes_supported"].fillna(0).min())), right_extra_frac=0.26, min_pad=0.9)

    phase_labels = []
    for phase in df["best_phase"]:
        phase_labels.append(f"Phase {int(phase)}" if pd.notna(phase) else "Phase NA")
    _annotate_aligned_value_labels(ax, df["n_genes_supported"].fillna(0), y, phase_labels)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    _clean_colorbar(cbar, "Best clinical phase")

    size_handles = [
        plt.scatter([], [], s=s, facecolor="#AFC6E9", edgecolor="white", linewidth=0.9)
        for s in [150, 360, 700]
    ]
    ax.legend(size_handles, ["low", "intermediate", "high"], title="Compound priority size", loc="lower right")

    fig.tight_layout()
    savefig_both(fig, outdir / "Fig_T8_compound_candidate_bubble.png", outdir / "Fig_T8_compound_candidate_bubble.pdf")

# -----------------------------------------------------------------------------
# Workbook and logs
# -----------------------------------------------------------------------------
def write_excel(out_xlsx: Path, sheets: Dict[str, pd.DataFrame]) -> None:
    try:
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
            for sheet, df in sheets.items():
                safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", sheet)[:31]
                (df if not df.empty else pd.DataFrame({"note": ["no data"]})).to_excel(xw, sheet_name=safe_name, index=False)
    except Exception as e:
        msg(f"[WARN] Excel export failed: {e}")


def log_versions(log_dir: Path) -> None:
    lines = [
        f"Python: {sys.version}",
        f"Platform: {platform.platform()}",
        f"Script: {Path(__file__).name}",
        f"Seed: {SEED}",
        f"numpy: {np.__version__}",
        f"pandas: {pd.__version__}",
        f"matplotlib: {matplotlib.__version__}",
    ]
    (log_dir / "package_versions.txt").write_text("\n".join(lines), encoding="utf-8")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    apply_q1_style()
    args = parse_args()

    manifest = pick_manifest(args.manifest)
    manifest_paths = read_manifest(manifest)
    root = infer_root(args.project_root, manifest_paths)

    gene_file = resolve_existing_file(
        args.gene_file,
        [
            os.path.join(root, os.path.basename(args.gene_file)),
            os.path.join(root, "OUT_v44", "genes", os.path.basename(args.gene_file)),
            os.path.join(root, "bioinf_scripts", "scripts", os.path.basename(args.gene_file)),
            os.path.join(root, "GEO_PREP", "GSE81089", "WGCNA_GSE81089_out_v12_signature_fix", "tables", os.path.basename(args.gene_file)),
            resolve_from_manifest(manifest_paths, [rf"/{re.escape(os.path.basename(args.gene_file))}$"]),
        ],
    )

    ferrdb_driver = resolve_ferrdb_file(root, manifest_paths, args.ferrdb_driver_csv)
    ferrdb_suppressor = resolve_ferrdb_file(root, manifest_paths, args.ferrdb_suppressor_csv)
    ferrdb_marker = resolve_ferrdb_file(root, manifest_paths, args.ferrdb_marker_csv)

    if not ferrdb_driver:
        raise FileNotFoundError(f"FerrDb driver CSV not found: {args.ferrdb_driver_csv}")
    if not ferrdb_suppressor:
        raise FileNotFoundError(f"FerrDb suppressor CSV not found: {args.ferrdb_suppressor_csv}")

    ferrdb_sets = build_ferrdb_sets(ferrdb_driver, ferrdb_suppressor, ferrdb_marker)
    if not as_bool(args.include_marker_genes):
        ferrdb_sets["known"] = set().union(ferrdb_sets["drivers"], ferrdb_sets["suppressors"])
        ferrdb_sets["markers"] = set()

    ppi_hard_default = os.path.join(root, PPI_TREE_RELATIVE)
    ppi_manifest_default = resolve_from_manifest(
        manifest_paths,
        [
            r"/ppi_network_Q1_v2/tables/node_centrality\.csv$",
            r"/ppi_network_Q1_v2/tables/.*node.*central.*\.csv$",
            r"/ppi_network_Q1_v2/tables/.*central.*\.csv$",
            r"/ppi_network_Q1_v2/tables/.*hub.*\.csv$",
        ],
    )
    ppi_file = resolve_existing_file(
        args.ppi_table or PPI_TREE_RELATIVE,
        [
            ppi_hard_default,
            ppi_manifest_default,
            os.path.join(root, "ppi_network_Q1_v2", "tables", "node_centrality.csv"),
        ],
    )

    outdir = args.outdir or os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1_ferrdb_advanced")
    outpaths = ensure_dirs(Path(normalize_path(outdir)))
    cache_dir = Path(normalize_path(args.cache_dir)) if args.cache_dir else outpaths["cache"]
    cache_dir.mkdir(parents=True, exist_ok=True)
    log_versions(outpaths["logs"])

    msg(f"[INFO] project_root           : {root}")
    msg(f"[INFO] manifest               : {manifest or 'not_found'}")
    msg(f"[INFO] gene_file              : {gene_file}")
    msg(f"[INFO] ferrdb_driver_csv      : {ferrdb_driver}")
    msg(f"[INFO] ferrdb_suppressor_csv  : {ferrdb_suppressor}")
    msg(f"[INFO] ferrdb_marker_csv      : {ferrdb_marker or 'not_found'}")
    msg(f"[INFO] ppi_file               : {ppi_file or 'not_found'}")
    msg(f"[INFO] outdir                 : {outdir}")
    msg(f"[INFO] cache_dir              : {cache_dir}")

    genes = read_gene_table(gene_file, ferrdb_sets)
    ppi_df = read_ppi_table(ppi_file)
    genes = genes.head(max(1, int(args.top_n))).copy()

    session = get_requests_session(max_retries=args.max_retries)

    link_rows: List[Dict[str, object]] = []
    chembl_meta_rows: List[Dict[str, object]] = []
    tractability_meta_rows: List[Dict[str, object]] = []

    for i, symbol in enumerate(genes["gene"].tolist(), start=1):
        msg(f"[QUERY {i}/{len(genes)}] {symbol}")
        try:
            link_rows.extend(query_dgidb(session, symbol, timeout=args.timeout, cache_dir=cache_dir))
        except Exception as e:
            msg(f"[WARN] DGIdb failed for {symbol}: {e}")
        time.sleep(args.sleep_sec)

        try:
            chembl_rows, chembl_meta = query_chembl(
                session=session,
                symbol=symbol,
                timeout=args.timeout,
                cache_dir=cache_dir,
                max_activities=args.max_chembl_activities_per_target,
                max_phase_queries=args.max_molecule_phase_queries_per_target,
                min_pchembl_for_phase_support=args.min_pchembl_for_phase_support,
            )
            link_rows.extend(chembl_rows)
            if chembl_meta:
                chembl_meta_rows.append(chembl_meta)
        except Exception as e:
            msg(f"[WARN] ChEMBL failed for {symbol}: {e}")
        time.sleep(args.sleep_sec)

        try:
            ot_rows, tract_meta = query_opentargets_known_drugs(
                session=session,
                symbol=symbol,
                timeout=args.timeout,
                cache_dir=cache_dir,
                max_known_drugs=args.max_ot_known_drugs_per_target,
            )
            link_rows.extend(ot_rows)
            if tract_meta:
                tractability_meta_rows.append(tract_meta)
        except Exception as e:
            msg(f"[WARN] Open Targets failed for {symbol}: {e}")
        time.sleep(args.sleep_sec)

    links = normalize_links(pd.DataFrame(link_rows))
    target_meta = normalize_target_metadata(chembl_meta_rows, tractability_meta_rows)
    summary = build_summary(genes, ppi_df, links, target_meta)

    if args.min_priority_display > 0:
        summary = summary.loc[summary["priority_score"] >= args.min_priority_display].copy()

    actionable = build_actionable_table(summary)
    evidence = build_evidence_matrix(summary)
    source_summary = build_source_summary(links)
    mechanism_summary = build_mechanism_summary(links)
    ferrdb_overlap_summary = build_ferrdb_overlap_summary(summary)
    tier_breakdown = build_tier_breakdown(summary)
    compound_candidates = build_compound_candidate_table(summary, links)
    gene_compound_edges = build_gene_compound_network_edges(summary, links)

    actionable.to_csv(outpaths["tables"] / "table_actionable_targets.tsv", sep="\t", index=False)
    links.to_csv(outpaths["tables"] / "table_target_compound_links.tsv", sep="\t", index=False)
    summary.to_csv(outpaths["tables"] / "table_target_priority_summary.tsv", sep="\t", index=False)
    evidence.to_csv(outpaths["tables"] / "table_target_evidence_matrix.tsv", sep="\t", index=False)
    source_summary.to_csv(outpaths["tables"] / "table_source_summary.tsv", sep="\t", index=False)
    mechanism_summary.to_csv(outpaths["tables"] / "table_mechanism_summary.tsv", sep="\t", index=False)
    ferrdb_overlap_summary.to_csv(outpaths["tables"] / "table_ferrdb_overlap_summary.tsv", sep="\t", index=False)
    tier_breakdown.to_csv(outpaths["tables"] / "table_priority_tier_breakdown.tsv", sep="\t", index=False)
    compound_candidates.to_csv(outpaths["tables"] / "table_candidate_compounds.tsv", sep="\t", index=False)
    gene_compound_edges.to_csv(outpaths["tables"] / "table_gene_compound_edges.tsv", sep="\t", index=False)
    target_meta.to_csv(outpaths["tables"] / "table_target_metadata.tsv", sep="\t", index=False)

    write_excel(outpaths["tables"] / "Target_Druggability_Results_Q1_FerrDb_Advanced.xlsx", {
        "actionable_targets": actionable,
        "target_compound_links": links,
        "priority_summary": summary,
        "evidence_matrix": evidence,
        "source_summary": source_summary,
        "mechanism_summary": mechanism_summary,
        "ferrdb_overlap_summary": ferrdb_overlap_summary,
        "tier_breakdown": tier_breakdown,
        "candidate_compounds": compound_candidates,
        "gene_compound_edges": gene_compound_edges,
        "target_metadata": target_meta,
    })

    fig_actionability_bubble(summary, outpaths["figures"])
    fig_priority_waterfall(summary, outpaths["figures"])
    fig_weighted_evidence_heatmap(summary, outpaths["figures"])
    fig_phase_distribution(links, outpaths["figures"])
    fig_ferrdb_class_distribution(summary, outpaths["figures"])
    fig_ferrdb_priority_composition(summary, outpaths["figures"])
    fig_structured_priority_panel(summary, outpaths["figures"])
    fig_compound_candidate_bubble(compound_candidates, outpaths["figures"])

    metadata = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "project_root": root,
        "manifest": manifest,
        "gene_file": gene_file,
        "ppi_file": ppi_file,
        "ferrdb_driver_csv": ferrdb_driver,
        "ferrdb_suppressor_csv": ferrdb_suppressor,
        "ferrdb_marker_csv": ferrdb_marker,
        "ppi_tree_relative": PPI_TREE_RELATIVE,
        "outdir": str(outpaths["root"]),
        "cache_dir": str(cache_dir),
        "n_signature_genes": int(len(genes)),
        "n_ferrdb_known": int(len(ferrdb_sets["known"])),
        "n_ferrdb_drivers": int(len(ferrdb_sets["drivers"])),
        "n_ferrdb_suppressors": int(len(ferrdb_sets["suppressors"])),
        "n_ferrdb_markers": int(len(ferrdb_sets["markers"])),
        "n_link_rows": int(len(links)),
        "n_actionable_targets": int(len(actionable)),
        "n_candidate_compounds": int(len(compound_candidates)),
    }
    (outpaths["logs"] / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    methods_notes = [
        "Signature genes were loaded dynamically from the ML-derived prognostic gene table.",
        "FerrDb priors were loaded dynamically from driver, suppressor, and optional marker tables.",
        "PPI support was integrated from the node-centrality table in the project PPI workflow.",
        "Target-compound evidence was aggregated from DGIdb, ChEMBL, and Open Targets known-drug records.",
        "Target prioritization combined expression/stability, network centrality, FerrDb prior, tractability, compound evidence, and clinical maturity.",
        "Compound nomination aggregated evidence across supported ML-derived targets and retained source provenance.",
    ]
    (outpaths["logs"] / "methods_notes.txt").write_text("\n".join(methods_notes), encoding="utf-8")

    msg("[DONE] Advanced FerrDb-integrated target druggability / compound mapping analysis completed.")


if __name__ == "__main__":
    main()
