#!/usr/bin/env python3
"""
Rana Salihoglu
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260418
np.random.seed(SEED)

DEFAULT_ROOT = "/home/rana/Desktop/PROJECT_16/test"
DEFAULT_MANIFESTS = [
    "/mnt/data/full_paths.txt",
    "/mnt/data/full_paths(2).txt",
    "full_paths.txt",
    "full_paths(2).txt",
]


# ------------------------------ CLI ---------------------------------------- #
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Q1-grade structural target prioritization")
    p.add_argument("--project_root", default=DEFAULT_ROOT)
    p.add_argument("--manifest", default="")
    p.add_argument("--gene_file", default="possible_prognostic_genes_FULLTRAIN.csv")
    p.add_argument("--priority_table", default="")
    p.add_argument("--outdir", default="")
    p.add_argument("--top_n", type=int, default=25)
    p.add_argument("--timeout", type=int, default=25)
    p.add_argument("--sleep", type=float, default=0.15)
    p.add_argument("--max_pdb_entries", type=int, default=50)
    p.add_argument("--cache_ttl_hours", type=float, default=168.0)
    p.add_argument("--disable_cache", action="store_true")
    return p.parse_args()


# ---------------------------- small helpers -------------------------------- #
def msg(*x) -> None:
    print(*x, flush=True)


def normalize_path(x: str) -> str:
    try:
        return str(Path(os.path.expanduser(x)).resolve())
    except Exception:
        return os.path.expanduser(x)


def requests_import():
    import requests
    return requests


def safe_bool(x) -> bool:
    if pd.isna(x):
        return False
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    s = str(x).strip().lower()
    return s in {"1", "true", "t", "yes", "y", "reviewed", "swiss-prot", "swissprot"}


def safe_float(x, default=np.nan) -> float:
    try:
        y = float(x)
        if math.isfinite(y):
            return y
        return default
    except Exception:
        return default


def safe_int(x, default=0) -> int:
    try:
        y = int(float(x))
        return y
    except Exception:
        return default


def pick_manifest(user_value: str) -> str:
    cands = [user_value] if user_value else []
    cands.extend(DEFAULT_MANIFESTS)
    for c in cands:
        if c and os.path.exists(os.path.expanduser(c)):
            return normalize_path(c)
    return ""


def read_manifest(path: str) -> List[str]:
    if not path or not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def infer_root(project_root: str, manifest_paths: Sequence[str]) -> str:
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
            ordered.append(xn)
            seen.add(xn)
    for x in ordered:
        if os.path.isfile(x):
            return x
    base = os.path.basename(target)
    for x in ordered:
        alt = os.path.join(os.path.dirname(x), base)
        if os.path.isfile(alt):
            return alt
    return ""


def ensure_dirs(outdir: Path) -> Dict[str, Path]:
    out = {
        "root": outdir,
        "tables": outdir / "tables",
        "figures": outdir / "figures",
        "logs": outdir / "logs",
        "cache": outdir / "cache",
    }
    for p in out.values():
        p.mkdir(parents=True, exist_ok=True)
    return out


def safe_read_table(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    ext = Path(path).suffix.lower()
    try:
        if ext in {".tsv", ".txt"}:
            return pd.read_csv(path, sep="\t")
        return pd.read_csv(path)
    except Exception:
        try:
            return pd.read_csv(path, sep=None, engine="python")
        except Exception:
            return pd.DataFrame()


def sanitize_gene(x: str) -> str:
    s = str(x).strip()
    if "__" in s:
        s = s.split("__")[-1]
    s = re.sub(r"\s+", "", s)
    s = s.replace(".", "-")
    return s.upper()


def pick_first_column(cols: Sequence[str], patterns: Sequence[str]) -> Optional[str]:
    low = {str(c).lower(): c for c in cols}
    for p in patterns:
        for k, c in low.items():
            if re.search(p, k):
                return c
    return None


def scaled(x: pd.Series) -> pd.Series:
    s = pd.to_numeric(x, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series(0.0, index=s.index, dtype=float)
    mn = float(s.min(skipna=True))
    mx = float(s.max(skipna=True))
    if not math.isfinite(mn) or not math.isfinite(mx) or mx <= mn:
        return pd.Series(np.where(s.notna(), 0.5, 0.0), index=s.index, dtype=float)
    return ((s - mn) / (mx - mn)).fillna(0.0)


def inverse_scaled(x: pd.Series) -> pd.Series:
    return 1.0 - scaled(x)


def placeholder_figure(path: Path, title: str, note: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    ax.text(0.5, 0.65, title, ha="center", va="center", fontsize=15, fontweight="bold")
    ax.text(0.5, 0.40, note, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def savefig(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def clip01(x: pd.Series) -> pd.Series:
    s = pd.to_numeric(x, errors="coerce").fillna(0.0)
    return s.clip(lower=0.0, upper=1.0)


# ---------------------------- caching -------------------------------------- #
class HTTPCache:
    def __init__(self, cache_dir: Path, ttl_hours: float = 168.0, enabled: bool = True):
        self.cache_dir = cache_dir
        self.ttl_seconds = float(ttl_hours) * 3600.0
        self.enabled = enabled

    def _key_to_path(self, key: str) -> Path:
        h = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{h}.json"

    def get(self, key: str):
        if not self.enabled:
            return None
        p = self._key_to_path(key)
        if not p.exists():
            return None
        age = time.time() - p.stat().st_mtime
        if age > self.ttl_seconds:
            return None
        try:
            with open(p, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return None

    def set(self, key: str, obj) -> None:
        if not self.enabled:
            return
        p = self._key_to_path(key)
        try:
            with open(p, "w", encoding="utf-8") as fh:
                json.dump(obj, fh)
        except Exception:
            pass


def cached_get_json(url: str, cache: HTTPCache, timeout: int = 25):
    key = f"GET::{url}"
    hit = cache.get(key)
    if hit is not None:
        return hit
    requests = requests_import()
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    js = r.json()
    cache.set(key, js)
    return js


def cached_post_json(url: str, payload: dict, cache: HTTPCache, timeout: int = 25):
    key = f"POST::{url}::{json.dumps(payload, sort_keys=True)}"
    hit = cache.get(key)
    if hit is not None:
        return hit
    requests = requests_import()
    r = requests.post(url, json=payload, timeout=timeout)
    r.raise_for_status()
    js = r.json()
    cache.set(key, js)
    return js


# --------------------------- input parsing --------------------------------- #
def read_gene_table(path: str) -> pd.DataFrame:
    df = safe_read_table(path)
    if df.empty:
        raise FileNotFoundError(f"Gene table not found or empty: {path}")
    gcol = pick_first_column(df.columns, [r"^gene$", r"gene.*symbol", r"symbol", r"feature"])
    if gcol is None:
        gcol = df.columns[0]
    zcol = pick_first_column(df.columns, [r"^z_meta$", r"meta.*z", r"^z$"])
    bagcol = pick_first_column(df.columns, [r"bag_frac", r"bag_freq", r"freq"])
    dircol = pick_first_column(df.columns, [r"direction"])
    out = pd.DataFrame({"gene": df[gcol].map(sanitize_gene)})
    out["z_meta"] = pd.to_numeric(df[zcol], errors="coerce") if zcol else np.nan
    out["bag_frac"] = pd.to_numeric(df[bagcol], errors="coerce") if bagcol else np.nan
    if dircol:
        out["direction"] = df[dircol].astype(str)
    else:
        out["direction"] = np.where(out["z_meta"].fillna(0) > 0, "Risk-promoting", "Protective")
    out["abs_z"] = out["z_meta"].abs()
    out = out.dropna(subset=["gene"])
    out = out.loc[out["gene"].astype(str).str.len() > 0].copy()
    out = out.sort_values(["abs_z", "bag_frac"], ascending=[False, False], na_position="last")
    out = out.drop_duplicates(subset=["gene"], keep="first").reset_index(drop=True)
    return out


def harmonize_priority_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    out = df.copy()
    if "gene" not in out.columns:
        gcol = pick_first_column(out.columns, [r"^gene$", r"symbol", r"target", r"feature"])
        if gcol:
            out = out.rename(columns={gcol: "gene"})
    if "gene" not in out.columns:
        return pd.DataFrame()
    out["gene"] = out["gene"].map(sanitize_gene)

    rename_candidates = {
        "priority_score": [r"^priority_score$", r"composite_priority_score", r"translational_priority_score", r"score_total"],
        "actionability_score": [r"^actionability_score$", r"evidence_score", r"score_compound", r"score_clinical"],
        "priority_tier": [r"^priority_tier$", r"structure_priority_tier"],
        "n_compounds": [r"^n_compounds$", r"compound"],
        "n_clinical_compounds": [r"clinical"],
        "best_phase": [r"best_phase", r"max_phase"],
        "gene_category": [r"gene_category", r"ferrdb"],
    }
    for new, patterns in rename_candidates.items():
        if new in out.columns:
            continue
        col = pick_first_column(out.columns, patterns)
        if col:
            out = out.rename(columns={col: new})

    keep = [c for c in ["gene", "priority_score", "actionability_score", "priority_tier",
                        "n_compounds", "n_clinical_compounds", "best_phase", "gene_category"] if c in out.columns]
    out = out[keep].copy()
    for c in ["priority_score", "actionability_score", "n_compounds", "n_clinical_compounds", "best_phase"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.drop_duplicates(subset=["gene"], keep="first")


# ---------------------------- remote queries -------------------------------- #
def query_uniprot(symbol: str, cache: HTTPCache, timeout: int = 25) -> Dict[str, object]:
    empty = {
        "uniprot_id": "",
        "protein_name": "",
        "reviewed": np.nan,
        "sequence_length": np.nan,
        "keywords": "",
        "gene_names": "",
        "alphafold_xref": "",
        "subcellular_location_note": "",
    }
    try:
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query=(gene_exact:{symbol}+AND+organism_id:9606)"
            "&format=json&size=10"
            "&fields=accession,reviewed,id,protein_name,gene_names,length,keyword,xref_alphafolddb,cc_subcellular_location"
        )
        js = cached_get_json(url, cache=cache, timeout=timeout)
        results = js.get("results", []) or []
        if not results:
            url2 = (
                "https://rest.uniprot.org/uniprotkb/search"
                f"?query=(gene:{symbol}+AND+organism_id:9606)"
                "&format=json&size=20"
                "&fields=accession,reviewed,id,protein_name,gene_names,length,keyword,xref_alphafolddb,cc_subcellular_location"
            )
            results = (cached_get_json(url2, cache=cache, timeout=timeout) or {}).get("results", []) or []
        if not results:
            return empty

        def reviewed_flag(rec: dict) -> int:
            et = str(rec.get("entryType", "")).lower()
            rv = rec.get("reviewed", None)
            if isinstance(rv, bool):
                return int(rv)
            if "reviewed" in et or "swiss" in et:
                return 1
            return 0

        def gene_match_score(rec: dict) -> int:
            gsec = rec.get("genes", []) or []
            names = []
            for g in gsec:
                gn = g.get("geneName", {}) or {}
                if gn.get("value"):
                    names.append(str(gn.get("value")))
                for syn in g.get("synonyms", []) or []:
                    if syn.get("value"):
                        names.append(str(syn.get("value")))
            names = {sanitize_gene(x) for x in names if x}
            return int(sanitize_gene(symbol) in names)

        results = sorted(results, key=lambda rec: (gene_match_score(rec), reviewed_flag(rec)), reverse=True)
        best = results[0]

        prot = ""
        pdsc = best.get("proteinDescription", {}) or {}
        for node in [
            ((pdsc.get("recommendedName") or {}).get("fullName") or {}).get("value", ""),
            ((pdsc.get("submissionNames") or [{}])[0].get("fullName") or {}).get("value", ""),
        ]:
            if node:
                prot = node
                break

        rv = best.get("reviewed", None)
        et = str(best.get("entryType", "")).lower()
        reviewed = bool(rv) if isinstance(rv, bool) else ("reviewed" in et or "swiss" in et)

        kws = []
        for kw in best.get("keywords", []) or []:
            nm = kw.get("name", "")
            if nm:
                kws.append(str(nm))

        gene_names = []
        for g in best.get("genes", []) or []:
            gn = g.get("geneName", {}) or {}
            if gn.get("value"):
                gene_names.append(str(gn.get("value")))
            for syn in g.get("synonyms", []) or []:
                if syn.get("value"):
                    gene_names.append(str(syn.get("value")))

        alphafold_xref = ""
        for x in best.get("uniProtKBCrossReferences", []) or []:
            db = str(x.get("database", ""))
            if db.upper() == "ALPHAFOLDDB":
                alphafold_xref = str(x.get("id", ""))
                break

        subloc_text = ""
        for comment in best.get("comments", []) or []:
            if str(comment.get("commentType", "")).upper() == "SUBCELLULAR LOCATION":
                parts = []
                for sl in comment.get("subcellularLocations", []) or []:
                    loc = ((sl.get("location") or {}).get("value")) or ""
                    if loc:
                        parts.append(str(loc))
                if parts:
                    subloc_text = ";".join(sorted(set(parts)))
                    break

        return {
            "uniprot_id": best.get("primaryAccession", "") or "",
            "protein_name": prot,
            "reviewed": reviewed,
            "sequence_length": safe_float((best.get("sequence") or {}).get("length", np.nan)),
            "keywords": ";".join(sorted(set(kws))),
            "gene_names": ";".join(sorted(set(gene_names))),
            "alphafold_xref": alphafold_xref,
            "subcellular_location_note": subloc_text,
        }
    except Exception:
        return empty


def query_rcsb(symbol: str, uniprot_id: str, cache: HTTPCache, timeout: int = 25) -> Dict[str, object]:
    empty = {"pdb_ids": "", "n_pdb_entries": 0, "has_pdb": False}
    try:
        hits = set()
        queries = [
            {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": symbol,
                    },
                },
                "return_type": "entry",
                "request_options": {"return_all_hits": True},
            }
        ]
        if uniprot_id:
            queries.append(
                {
                    "query": {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": uniprot_id,
                        },
                    },
                    "return_type": "entry",
                    "request_options": {"return_all_hits": True},
                }
            )

        for q in queries:
            try:
                js = cached_post_json("https://search.rcsb.org/rcsbsearch/v2/query", q, cache=cache, timeout=timeout) or {}
                for x in js.get("result_set", []) or []:
                    ident = x.get("identifier", "")
                    if ident:
                        hits.add(str(ident))
            except Exception:
                pass

        pdb_ids = sorted(hits)
        return {"pdb_ids": ";".join(pdb_ids), "n_pdb_entries": len(pdb_ids), "has_pdb": len(pdb_ids) > 0}
    except Exception:
        return empty


def query_rcsb_entry_details(
    pdb_ids: Sequence[str], cache: HTTPCache, timeout: int = 25, max_entries: int = 50
) -> Dict[str, object]:
    empty = {
        "best_resolution": np.nan,
        "median_resolution": np.nan,
        "experimental_methods": "",
        "xray_count": 0,
        "em_count": 0,
        "nmr_count": 0,
        "ligand_bound_entries": 0,
        "total_nonpolymer_instances": 0,
        "median_polymer_entity_count": np.nan,
        "median_deposited_model_count": np.nan,
    }
    if not pdb_ids:
        return empty

    methods = []
    resolutions = []
    xray_count = em_count = nmr_count = 0
    ligand_bound_entries = 0
    nonpolymer_counts = []
    polymer_counts = []
    model_counts = []

    for pdb_id in list(pdb_ids)[:max_entries]:
        try:
            js = cached_get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", cache=cache, timeout=timeout) or {}
            exptl = js.get("exptl", []) or []
            method_values = sorted({str(x.get("method", "")) for x in exptl if x.get("method")})
            if method_values:
                methods.extend(method_values)
            joined = ";".join(method_values).upper()
            if "X-RAY" in joined:
                xray_count += 1
            if "ELECTRON MICROSCOPY" in joined or "CRYO" in joined:
                em_count += 1
            if "NMR" in joined:
                nmr_count += 1

            rr = ((js.get("rcsb_entry_info", {}) or {}).get("resolution_combined", []) or [])
            for r0 in rr:
                val = safe_float(r0)
                if math.isfinite(val):
                    resolutions.append(val)

            info = js.get("rcsb_entry_info", {}) or {}
            npoly = safe_int(info.get("nonpolymer_entity_count", 0))
            ppoly = safe_int(info.get("polymer_entity_count", 0))
            nmodels = safe_int(info.get("deposited_model_count", 0))
            nonpolymer_counts.append(npoly)
            polymer_counts.append(ppoly)
            model_counts.append(nmodels)
            if npoly > 0:
                ligand_bound_entries += 1
        except Exception:
            continue

    return {
        "best_resolution": min(resolutions) if resolutions else np.nan,
        "median_resolution": float(np.median(resolutions)) if resolutions else np.nan,
        "experimental_methods": ";".join(sorted(set(methods))),
        "xray_count": int(xray_count),
        "em_count": int(em_count),
        "nmr_count": int(nmr_count),
        "ligand_bound_entries": int(ligand_bound_entries),
        "total_nonpolymer_instances": int(np.nansum(nonpolymer_counts)) if nonpolymer_counts else 0,
        "median_polymer_entity_count": float(np.nanmedian(polymer_counts)) if polymer_counts else np.nan,
        "median_deposited_model_count": float(np.nanmedian(model_counts)) if model_counts else np.nan,
    }


def query_alphafold(uniprot_id: str, alphafold_xref: str, cache: HTTPCache, timeout: int = 25) -> Dict[str, object]:
    empty = {"has_alphafold": False, "alphafold_url": "", "alphafold_model_id": ""}
    model_id = alphafold_xref or uniprot_id
    if not model_id:
        return empty
    try:
        js = cached_get_json(f"https://alphafold.ebi.ac.uk/api/prediction/{model_id}", cache=cache, timeout=timeout)
        if isinstance(js, list) and len(js) > 0:
            return {
                "has_alphafold": True,
                "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{model_id}",
                "alphafold_model_id": model_id,
            }
        return empty
    except Exception:
        return empty


# ------------------------------- scoring ------------------------------------ #
def infer_target_class(gene: str, protein_name: str, keywords: str, subloc: str) -> str:
    txt = f"{gene} {protein_name} {keywords} {subloc}".upper()
    if any(k in txt for k in ["KINASE", "TYROSINE-PROTEIN KINASE", "SERINE/THREONINE", "NON-RECEPTOR TYROSINE"]):
        return "Kinase"
    if any(k in txt for k in ["RECEPTOR", "INTERLEUKIN", "CHEMOKINE", "CYTOKINE", "TOLL-LIKE", "IMMUNOGLOBULIN", "TRANSMEMBRANE"]):
        return "Receptor / immune surface"
    if any(k in txt for k in ["INTEGRIN", "SELECTIN", "ADHESION", "MEMBRANE", "CELL SURFACE", "EXTRACELLULAR"]):
        return "Cell adhesion / membrane"
    if any(k in txt for k in ["TRANSFERASE", "GLYCOSYL", "PARP", "ENZYME", "OXIDASE", "PHOSPHATASE", "DEHYDROGENASE", "HYDROLASE"]):
        return "Enzyme / catalytic"
    if any(k in txt for k in ["TRANSCRIPTION", "NUCLEAR", "IRF", "MEF2", "BCL11", "BTG", "CHROMATIN"]):
        return "Transcription / nuclear"
    return "Other / mixed"


def infer_subcellular_class(s: str) -> str:
    txt = str(s).upper()
    if any(k in txt for k in ["CELL MEMBRANE", "PLASMA MEMBRANE", "EXTRACELLULAR", "SECRETED", "CELL SURFACE"]):
        return "Membrane/extracellular"
    if "NUCLEUS" in txt:
        return "Nuclear"
    if any(k in txt for k in ["CYTOPLASM", "CYTOSOL"]):
        return "Cytoplasmic"
    if any(k in txt for k in ["MITOCHONDR", "ENDOPLASMIC RETICULUM", "GOLGI", "LYSOSOME"]):
        return "Organelle-associated"
    return "Unresolved"


def classify_dominant_method(methods: str) -> str:
    m = str(methods).upper()
    if "X-RAY" in m:
        return "X-ray"
    if "ELECTRON MICROSCOPY" in m or "CRYO" in m:
        return "Cryo-EM"
    if "NMR" in m:
        return "NMR"
    return "No experimental structure"


def docking_readiness_label(has_pdb: bool, best_resolution: float, has_alphafold: bool, ligand_bound_entries: int) -> str:
    if bool(has_pdb) and math.isfinite(best_resolution) and best_resolution <= 3.0 and ligand_bound_entries >= 1:
        return "Experimental-ready"
    if bool(has_pdb) and math.isfinite(best_resolution) and best_resolution <= 4.0:
        return "Structure-available"
    if bool(has_alphafold):
        return "AlphaFold-only"
    return "No structure"


def structure_source_label(has_pdb: bool, has_alphafold: bool) -> str:
    if bool(has_pdb) and bool(has_alphafold):
        return "PDB + AlphaFold"
    if bool(has_pdb):
        return "PDB only"
    if bool(has_alphafold):
        return "AlphaFold only"
    return "No structure"


def assign_tiers(score: pd.Series) -> pd.Series:
    s = pd.to_numeric(score, errors="coerce")
    valid = s.dropna()
    if valid.empty:
        return pd.Series(["Tier D"] * len(s), index=s.index)
    if valid.nunique() == 1:
        return pd.Series(["Tier B"] * len(s), index=s.index)

    q75 = valid.quantile(0.75)
    q50 = valid.quantile(0.50)
    q25 = valid.quantile(0.25)

    def one(x: float) -> str:
        if not math.isfinite(x):
            return "Tier D"
        if x >= q75:
            return "Tier A"
        if x >= q50:
            return "Tier B"
        if x >= q25:
            return "Tier C"
        return "Tier D"

    return s.map(one)


def score_structural_priority(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    def num_series(name: str, default: float = np.nan) -> pd.Series:
        if name in out.columns:
            return pd.to_numeric(out[name], errors="coerce")
        return pd.Series(default, index=out.index, dtype=float)

    def bool_series(name: str, default: bool = False) -> pd.Series:
        if name in out.columns:
            return out[name].map(safe_bool)
        return pd.Series(default, index=out.index, dtype=bool)

    out["reviewed"] = bool_series("reviewed", False)
    out["has_pdb"] = bool_series("has_pdb", False)
    out["has_alphafold"] = bool_series("has_alphafold", False)

    out["n_pdb_entries"] = num_series("n_pdb_entries", 0).fillna(0)
    out["best_resolution"] = num_series("best_resolution", np.nan)
    out["median_resolution"] = num_series("median_resolution", np.nan)
    out["xray_count"] = num_series("xray_count", 0).fillna(0)
    out["em_count"] = num_series("em_count", 0).fillna(0)
    out["nmr_count"] = num_series("nmr_count", 0).fillna(0)
    out["sequence_length"] = num_series("sequence_length", np.nan)
    out["priority_score"] = num_series("priority_score", np.nan)
    out["actionability_score"] = num_series("actionability_score", np.nan)
    out["bag_frac"] = num_series("bag_frac", np.nan)
    out["abs_z"] = num_series("abs_z", np.nan)
    out["n_compounds"] = num_series("n_compounds", np.nan)
    out["n_clinical_compounds"] = num_series("n_clinical_compounds", np.nan)
    out["best_phase"] = num_series("best_phase", np.nan)
    out["ligand_bound_entries"] = num_series("ligand_bound_entries", 0).fillna(0)
    out["total_nonpolymer_instances"] = num_series("total_nonpolymer_instances", 0).fillna(0)

    out["pdb_presence_score"] = out["has_pdb"].astype(float)
    out["alphafold_score"] = out["has_alphafold"].astype(float)
    out["reviewed_score"] = out["reviewed"].astype(float)
    out["pdb_density_score"] = scaled(out["n_pdb_entries"])
    out["resolution_quality_score"] = np.where(out["has_pdb"], inverse_scaled(out["best_resolution"]).fillna(0), 0.0)
    out["experimental_diversity_score"] = scaled(out["xray_count"] + out["em_count"] + out["nmr_count"])
    out["ligand_bound_score"] = scaled(out["ligand_bound_entries"])
    out["cofactor_density_score"] = scaled(out["total_nonpolymer_instances"])
    out["sequence_length_scaled"] = scaled(out["sequence_length"].fillna(out["sequence_length"].median()))

    target_class = out.get("target_class", pd.Series("Other / mixed", index=out.index)).fillna("Other / mixed").astype(str)
    subcell_class = out.get("subcellular_class", pd.Series("Unresolved", index=out.index)).fillna("Unresolved").astype(str)

    out["surface_target_score"] = target_class.isin(["Receptor / immune surface", "Cell adhesion / membrane"]).astype(float)
    out["kinase_target_score"] = target_class.eq("Kinase").astype(float)
    out["enzyme_target_score"] = target_class.eq("Enzyme / catalytic").astype(float)
    out["nuclear_penalty"] = target_class.eq("Transcription / nuclear").astype(float) * 0.18
    out["subcell_bonus"] = subcell_class.eq("Membrane/extracellular").astype(float) * 0.12

    pri_or_act = out["priority_score"].combine_first(out["actionability_score"]).fillna(0)
    act_or_pri = out["actionability_score"].combine_first(out["priority_score"]).fillna(0)
    compound_signal = (
        0.45 * scaled(out["n_compounds"].fillna(0)) +
        0.35 * scaled(out["n_clinical_compounds"].fillna(0)) +
        0.20 * scaled(out["best_phase"].fillna(0))
    )

    out["transcriptomic_support_score"] = clip01(
        0.40 * scaled(out["abs_z"].fillna(0)) +
        0.30 * scaled(out["bag_frac"].fillna(0)) +
        0.20 * scaled(pri_or_act) +
        0.10 * compound_signal
    )

    out["structural_evidence_score"] = clip01(
        0.22 * out["pdb_presence_score"] +
        0.14 * out["pdb_density_score"] +
        0.14 * out["resolution_quality_score"] +
        0.10 * out["experimental_diversity_score"] +
        0.10 * out["alphafold_score"] +
        0.10 * out["reviewed_score"] +
        0.12 * out["ligand_bound_score"] +
        0.08 * out["cofactor_density_score"]
    )

    out["ligandability_proxy_score"] = clip01(
        0.22 * out["kinase_target_score"] +
        0.18 * out["enzyme_target_score"] +
        0.12 * out["surface_target_score"] +
        0.12 * out["ligand_bound_score"] +
        0.10 * out["cofactor_density_score"] +
        0.14 * scaled(act_or_pri) +
        0.07 * out["pdb_presence_score"] +
        0.05 * out["resolution_quality_score"] +
        out["subcell_bonus"] -
        out["nuclear_penalty"]
    )

    out["docking_readiness_score"] = clip01(
        0.42 * np.where(out["has_pdb"], out["resolution_quality_score"].clip(lower=0.20), 0.0) +
        0.18 * out["alphafold_score"] +
        0.16 * out["pdb_density_score"] +
        0.14 * out["ligand_bound_score"] +
        0.10 * out["cofactor_density_score"]
    )

    out["translational_priority_score"] = clip01(
        0.34 * out["structural_evidence_score"] +
        0.28 * out["ligandability_proxy_score"] +
        0.20 * out["transcriptomic_support_score"] +
        0.10 * compound_signal +
        0.08 * scaled(pri_or_act)
    )

    out["docking_readiness"] = [
        docking_readiness_label(pdb, res, af, lig)
        for pdb, res, af, lig in zip(
            out["has_pdb"], out["best_resolution"].fillna(np.nan),
            out["has_alphafold"], out["ligand_bound_entries"]
        )
    ]
    out["structure_source"] = [
        structure_source_label(pdb, af)
        for pdb, af in zip(out["has_pdb"], out["has_alphafold"])
    ]
    out["dominant_method"] = out.get("experimental_methods", pd.Series("", index=out.index)).fillna("").map(classify_dominant_method)
    out["structure_priority_tier"] = assign_tiers(out["translational_priority_score"])
    out["priority_rank"] = out["translational_priority_score"].rank(ascending=False, method="dense")
    out["docking_recommendation"] = np.select(
        [
            out["docking_readiness"].eq("Experimental-ready"),
            out["docking_readiness"].eq("Structure-available"),
            out["docking_readiness"].eq("AlphaFold-only"),
        ],
        [
            "Proceed with experimental-structure-guided docking",
            "Proceed with cautious structure-based screening",
            "Use model-guided exploratory docking only",
        ],
        default="Defer docking; prioritize target biology first",
    )
    out = out.sort_values(
        ["translational_priority_score", "structural_evidence_score", "ligandability_proxy_score", "transcriptomic_support_score"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    return out


# ---------------------------- plotting -------------------------------------- #

def select_plot_subset(df: pd.DataFrame, top_n: int = 20, min_priority: float = 0.0) -> pd.DataFrame:
    tmp = df.copy()
    if tmp.empty:
        return tmp
    tmp["translational_priority_score"] = pd.to_numeric(tmp["translational_priority_score"], errors="coerce")
    tmp["structural_evidence_score"] = pd.to_numeric(tmp["structural_evidence_score"], errors="coerce")
    tmp["ligandability_proxy_score"] = pd.to_numeric(tmp["ligandability_proxy_score"], errors="coerce")
    tmp = tmp.sort_values(["translational_priority_score", "structural_evidence_score", "ligandability_proxy_score"],
                          ascending=[False, False, False]).reset_index(drop=True)
    if min_priority > 0:
        tmp = tmp.loc[tmp["translational_priority_score"].fillna(0) >= min_priority].copy()
    if len(tmp) > top_n:
        tmp = tmp.head(top_n).copy()
    return tmp

def apply_publication_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", alpha=0.22, linewidth=0.7)

def annotate_top_points(ax, df: pd.DataFrame, x: str, y: str, label_col: str = "gene", n: int = 10,
                        dx: float = 0.012, dy: float = 0.006) -> None:
    if df.empty:
        return
    tmp = df.sort_values("translational_priority_score", ascending=False).head(min(n, len(df)))
    placed = []
    for _, r in tmp.iterrows():
        xx = float(r[x]) if pd.notna(r[x]) else 0.0
        yy = float(r[y]) if pd.notna(r[y]) else 0.0
        label = str(r[label_col])
        # simple deterministic vertical repel
        yoff = dy
        for px, py in placed:
            if abs(xx - px) < 0.03 and abs((yy + yoff) - py) < 0.03:
                yoff += dy * 1.6
        ax.text(xx + dx, yy + yoff, label, fontsize=9, ha="left", va="bottom")
        placed.append((xx + dx, yy + yoff))



def plot_priority_bubble(df: pd.DataFrame, out: Path) -> None:
    top = select_plot_subset(df, top_n=20)
    if top.empty:
        placeholder_figure(out, "Structural prioritization", "No targets available")
        return
    plot_df = top.copy()
    zero_mask = plot_df["structural_evidence_score"].fillna(0) <= 1e-9
    if zero_mask.any():
        offsets = np.linspace(-0.018, 0.018, zero_mask.sum())
        plot_df.loc[zero_mask, "structural_evidence_plot"] = offsets
    plot_df.loc[~zero_mask, "structural_evidence_plot"] = plot_df.loc[~zero_mask, "structural_evidence_score"]
    plot_df = plot_df.sort_values(["structural_evidence_score", "translational_priority_score"], ascending=[True, False])
    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=plot_df["gene"].tolist(), ordered=True)

    fig, ax = plt.subplots(figsize=(11.0, 8.2))
    sizes = 220 + 1500 * scaled(plot_df["translational_priority_score"])
    sc = ax.scatter(
        plot_df["structural_evidence_plot"],
        plot_df["gene"],
        s=sizes,
        c=clip01(plot_df["ligandability_proxy_score"]),
        cmap="viridis",
        edgecolor="black",
        linewidth=0.7,
        alpha=0.94,
    )
    ax.axvline(0, color="#8F8F8F", linestyle=":", linewidth=1.0)
    ax.set_xlabel("Structural evidence score", fontweight="bold")
    ax.set_ylabel("")
    ax.set_title("Structural prioritization of LUAD ferroptosis-associated candidates", fontweight="bold")
    ax.text(0.01, 0.02, "Targets with no experimental structure are jittered around 0 for visibility.",
            transform=ax.transAxes, fontsize=9)
    apply_publication_axes(ax)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Ligandability proxy score", fontweight="bold")
    savefig(fig, out)



def plot_score_scatter(df: pd.DataFrame, out: Path) -> None:
    top = select_plot_subset(df, top_n=25)
    if top.empty:
        placeholder_figure(out, "Integrated ranking", "No targets available")
        return
    palette = {
        "Kinase": "#4C97C9",
        "Receptor / immune surface": "#B48FD3",
        "Enzyme / catalytic": "#D98BC2",
        "Cell adhesion / membrane": "#C8C84B",
        "Transcription / nuclear": "#A6D7E2",
        "Other / mixed": "#68BE63",
    }
    plot_df = top.copy()
    zero_mask = plot_df["structural_evidence_score"].fillna(0) <= 1e-9
    if zero_mask.any():
        offsets = np.linspace(-0.015, 0.015, zero_mask.sum())
        plot_df.loc[zero_mask, "structural_evidence_plot"] = offsets
    plot_df.loc[~zero_mask, "structural_evidence_plot"] = plot_df.loc[~zero_mask, "structural_evidence_score"]

    fig, ax = plt.subplots(figsize=(11.8, 8.2))
    for cls, sub in plot_df.groupby("target_class", dropna=False):
        ax.scatter(
            sub["structural_evidence_plot"],
            sub["translational_priority_score"],
            s=180 + 900 * scaled(sub["ligandability_proxy_score"]),
            color=palette.get(cls, "#888888"),
            alpha=0.86,
            edgecolor="black",
            linewidth=0.55,
            label=cls,
        )
    annotate_top_points(ax, plot_df, "structural_evidence_plot", "translational_priority_score", n=12)
    ax.axvline(0, color="#8F8F8F", linestyle=":", linewidth=1.0)
    ax.set_xlabel("Structural evidence score", fontweight="bold")
    ax.set_ylabel("Translational priority score", fontweight="bold")
    ax.set_title("Integrated structural versus translational target ranking", fontweight="bold")
    apply_publication_axes(ax)
    ax.legend(frameon=True, fontsize=9, title="Target class", loc="upper left")
    savefig(fig, out)



def plot_evidence_heatmap(df: pd.DataFrame, out: Path) -> None:
    top = select_plot_subset(df, top_n=15)
    if top.empty:
        placeholder_figure(out, "Structural evidence matrix", "No targets available")
        return
    top = top.sort_values(["translational_priority_score", "structural_evidence_score"], ascending=[False, False]).copy()
    mat = pd.DataFrame({
        "PDB": top["has_pdb"].astype(float).values,
        "AlphaFold": top["has_alphafold"].astype(float).values,
        "Ligand-bound": clip01(top["ligand_bound_score"]).values,
        "Resolution": clip01(top["resolution_quality_score"]).values,
        "Ligandability": clip01(top["ligandability_proxy_score"]).values,
        "Priority": clip01(top["translational_priority_score"]).values,
    }, index=top["gene"])
    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    im = ax.imshow(mat.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    ax.set_title("Structural evidence matrix for top-ranked targets", fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = float(mat.iat[i, j])
            txt = f"{val:.2f}" if j >= 2 else f"{int(round(val))}"
            ax.text(j, i, txt, ha="center", va="center",
                    color="white" if val < 0.55 else "black", fontsize=8.5)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Scaled score", fontweight="bold")
    savefig(fig, out)



def plot_readiness_bar(df: pd.DataFrame, out: Path) -> None:
    if df.empty or "docking_readiness" not in df.columns:
        placeholder_figure(out, "Docking readiness", "No docking-readiness data")
        return
    order = ["Experimental-ready", "Structure-available", "AlphaFold-only", "No structure"]
    cnt = df["docking_readiness"].value_counts().reindex(order).fillna(0)
    cnt = cnt.loc[cnt > 0]
    if cnt.sum() == 0:
        placeholder_figure(out, "Docking readiness", "All categories empty")
        return
    colors = {"Experimental-ready": "#EF8A62", "Structure-available": "#D9BD8D", "AlphaFold-only": "#A6D354", "No structure": "#C9CED6"}
    fig, ax = plt.subplots(figsize=(9.8, 6.0))
    bars = ax.bar(cnt.index, cnt.values, color=[colors[i] for i in cnt.index], edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Number of targets", fontweight="bold")
    ax.set_title("Docking readiness strata across prioritized structural targets", fontweight="bold")
    ax.tick_params(axis="x", rotation=22)
    means = df.groupby("docking_readiness", dropna=False)["translational_priority_score"].mean()
    for b, cat in zip(bars, cnt.index):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.08,
                f"n={int(cnt[cat])}\nmean={means.get(cat, np.nan):.2f}",
                ha="center", va="bottom", fontsize=9)
    apply_publication_axes(ax)
    ax.grid(False, axis="x")
    savefig(fig, out)



def plot_method_priority(df: pd.DataFrame, out: Path) -> None:
    tmp = df.copy()
    if tmp.empty or "dominant_method" not in tmp.columns:
        placeholder_figure(out, "Method summary", "No structural method data")
        return
    summ = tmp.groupby("dominant_method", dropna=False)["translational_priority_score"].agg(["mean", "median", "count"]).reset_index()
    summ = summ.loc[summ["count"] > 0].sort_values("mean", ascending=True)
    if summ.empty:
        placeholder_figure(out, "Method summary", "No structural method data")
        return
    colors = ["#56268E", "#AEB1CF", "#73A2C6", "#D95F5F", "#909090"]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    bars = ax.barh(summ["dominant_method"], summ["mean"], color=colors[:len(summ)], edgecolor="black")
    ax.set_xlabel("Mean translational priority score", fontweight="bold")
    ax.set_title("Priority distribution by dominant structural method", fontweight="bold")
    for b, (_, r) in zip(bars, summ.iterrows()):
        ax.text(b.get_width() + 0.01, b.get_y() + b.get_height()/2,
                f"n={int(r['count'])}; median={r['median']:.2f}", va="center", fontsize=9)
    apply_publication_axes(ax)
    savefig(fig, out)



def plot_direction_readiness_heatmap(df: pd.DataFrame, out: Path) -> None:
    tmp = df.copy()
    if tmp.empty or "direction" not in tmp.columns or "docking_readiness" not in tmp.columns:
        placeholder_figure(out, "Direction vs readiness", "No direction/readiness data")
        return
    order = ["Experimental-ready", "Structure-available", "AlphaFold-only", "No structure"]
    ctab = pd.crosstab(tmp["direction"], tmp["docking_readiness"]).reindex(columns=order).fillna(0)
    if ctab.empty:
        placeholder_figure(out, "Direction vs readiness", "Cross-tabulation empty")
        return
    rownorm = ctab.div(ctab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    im = ax.imshow(rownorm.values, cmap="GnBu", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(np.arange(rownorm.shape[1]))
    ax.set_xticklabels(rownorm.columns, rotation=28, ha="right")
    ax.set_yticks(np.arange(rownorm.shape[0]))
    ax.set_yticklabels(rownorm.index)
    ax.set_title("Docking-readiness distribution by prognostic direction", fontweight="bold")
    for i in range(rownorm.shape[0]):
        for j in range(rownorm.shape[1]):
            val = rownorm.iat[i, j]
            n = int(ctab.iat[i, j])
            ax.text(j, i, f"{val:.2f}\n(n={n})", ha="center", va="center", color="black", fontsize=8.5)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Row-normalized fraction", fontweight="bold")
    savefig(fig, out)



def plot_structure_source_summary(df: pd.DataFrame, out: Path) -> None:
    tmp = df.copy()
    if tmp.empty or "structure_source" not in tmp.columns:
        placeholder_figure(out, "Structure source summary", "No structure-source data")
        return
    order = ["PDB + AlphaFold", "PDB only", "AlphaFold only", "No structure"]
    ctab = pd.crosstab(tmp["structure_source"], tmp["direction"]).reindex(order).fillna(0)
    ctab = ctab.loc[ctab.sum(axis=1) > 0]
    if ctab.sum().sum() == 0:
        placeholder_figure(out, "Structure source summary", "No structure-source data")
        return
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    left = np.zeros(len(ctab))
    colors = {"Risk-promoting": "#C53929", "Protective": "#2C7FB8"}
    for direction in ["Risk-promoting", "Protective"]:
        vals = ctab[direction] if direction in ctab.columns else pd.Series(0, index=ctab.index)
        ax.barh(ctab.index, vals.values, left=left, label=direction, color=colors.get(direction, "#888888"), edgecolor="white")
        left += vals.values
    ax.set_xlabel("Targets", fontweight="bold")
    ax.set_ylabel("Structure source", fontweight="bold")
    ax.set_title("Experimental and predicted structure-source availability", fontweight="bold")
    for i, total in enumerate(ctab.sum(axis=1).values):
        ax.text(total + 0.08, i, f"n={int(total)}", va="center", fontsize=9)
    ax.legend(frameon=True)
    apply_publication_axes(ax)
    ax.grid(False, axis="y")
    savefig(fig, out)



def plot_reviewed_summary(df: pd.DataFrame, out: Path) -> None:
    tmp = df.copy()
    if tmp.empty or "reviewed" not in tmp.columns:
        placeholder_figure(out, "UniProt review status", "No review-status data")
        return
    tmp["review_status"] = np.where(tmp["reviewed"], "Reviewed (Swiss-Prot)", "Unreviewed")
    summ = tmp.groupby("review_status", dropna=False)["translational_priority_score"].agg(["mean", "median", "count"]).reset_index()
    if summ.empty:
        placeholder_figure(out, "UniProt review status", "No review-status data")
        return
    order = ["Reviewed (Swiss-Prot)", "Unreviewed"]
    summ["review_status"] = pd.Categorical(summ["review_status"], categories=order, ordered=True)
    summ = summ.sort_values("review_status")
    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    bars = ax.bar(summ["review_status"], summ["mean"], color=["#284B63", "#AFAFAF"][:len(summ)], edgecolor="black")
    ax.set_ylim(0, max(1.0, float(summ["mean"].max()) * 1.35))
    ax.set_ylabel("Mean translational priority", fontweight="bold")
    ax.set_xlabel("UniProt curation status", fontweight="bold")
    ax.set_title("Priority profile by UniProt review status", fontweight="bold")
    for b, (_, r) in zip(bars, summ.iterrows()):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.02,
                f"n={int(r['count'])}; median={r['median']:.2f}", ha="center", fontsize=9)
    apply_publication_axes(ax)
    ax.grid(False, axis="x")
    savefig(fig, out)


def plot_tier_distribution(df: pd.DataFrame, out: Path) -> None:
    tmp = df.copy()
    if tmp.empty or "structure_priority_tier" not in tmp.columns:
        placeholder_figure(out, "Tier distribution", "No tier data")
        return
    order = ["Tier A", "Tier B", "Tier C", "Tier D"]
    ctab = pd.crosstab(tmp["structure_priority_tier"], tmp["direction"]).reindex(order).fillna(0)
    ctab = ctab.loc[ctab.sum(axis=1) > 0]
    if ctab.sum().sum() == 0:
        placeholder_figure(out, "Tier distribution", "Tier counts are empty")
        return
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    x = np.arange(len(ctab.index))
    width = 0.38
    vals_risk = ctab["Risk-promoting"].values if "Risk-promoting" in ctab.columns else np.zeros(len(x))
    vals_prot = ctab["Protective"].values if "Protective" in ctab.columns else np.zeros(len(x))
    ax.bar(x - width/2, vals_risk, width=width, label="Risk-promoting", color="#C53929", edgecolor="black")
    ax.bar(x + width/2, vals_prot, width=width, label="Protective", color="#2C7FB8", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(ctab.index)
    ax.set_ylabel("Targets", fontweight="bold")
    ax.set_xlabel("Structure-priority tier", fontweight="bold")
    ax.set_title("Structure-priority tier distribution", fontweight="bold")
    for xi, r, p in zip(x, vals_risk, vals_prot):
        total = r + p
        ax.text(xi, total + 0.08, f"n={int(total)}", ha="center", fontsize=9)
    ax.legend(frameon=True)
    apply_publication_axes(ax)
    ax.grid(False, axis="x")
    savefig(fig, out)


def plot_shortlist_lollipop(df: pd.DataFrame, out: Path) -> None:
    top = df.head(min(12, len(df))).copy()
    if top.empty:
        placeholder_figure(out, "Shortlist", "No shortlisted targets")
        return
    top = top.sort_values("translational_priority_score", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    y = np.arange(len(top))
    ax.hlines(y, 0, top["translational_priority_score"], color="#808080", alpha=0.7, linewidth=1.5)
    sc = ax.scatter(
        top["translational_priority_score"],
        y,
        s=180 + 650 * scaled(top["structural_evidence_score"]),
        c=top["ligandability_proxy_score"],
        cmap="plasma",
        edgecolor="black",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlabel("Translational priority score", fontweight="bold")
    ax.set_title("Shortlisted structure-guided follow-up targets", fontweight="bold")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Ligandability proxy score", fontweight="bold")
    savefig(fig, out)


# ---------------------------- exports --------------------------------------- #
def write_excel_bundle(path: Path, sheets: Dict[str, pd.DataFrame]) -> None:
    try:
        with pd.ExcelWriter(path, engine="openpyxl") as xw:
            for name, df in sheets.items():
                if df is None or df.empty:
                    pd.DataFrame({"note": ["no data"]}).to_excel(xw, sheet_name=name[:31], index=False)
                else:
                    df.to_excel(xw, sheet_name=name[:31], index=False)
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


# -------------------------------- main -------------------------------------- #
def main() -> None:
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
        ],
    )
    priority_candidates = []
    if args.priority_table and os.path.isfile(normalize_path(args.priority_table)):
        priority_candidates.append(args.priority_table)
    priority_candidates.extend([
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1_ferrdb_refined", "tables", "table_target_priority_summary.tsv"),
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1_ferrdb_advanced", "tables", "table_target_priority_summary.tsv"),
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1_ferrdb_advanced_fixed", "tables", "table_target_priority_summary.tsv"),
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1", "tables", "table_target_priority_summary.tsv"),
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1_updated", "tables", "table_target_priority_summary.tsv"),
        os.path.join(root, "bioinf_scripts", "scripts", "target_druggability_q1", "tables", "table_actionable_targets.tsv"),
    ])
    priority_table = resolve_existing_file(args.priority_table, priority_candidates)

    outdir = args.outdir or os.path.join(root, "bioinf_scripts", "scripts", "structural_target_prioritization_q1_enhanced")
    outpaths = ensure_dirs(Path(normalize_path(outdir)))
    log_versions(outpaths["logs"])
    cache = HTTPCache(outpaths["cache"], ttl_hours=args.cache_ttl_hours, enabled=(not args.disable_cache))

    msg(f"[INFO] project_root: {root}")
    msg(f"[INFO] manifest    : {manifest or 'not_found'}")
    msg(f"[INFO] gene_file   : {gene_file}")
    msg(f"[INFO] priority    : {priority_table or 'not_found'}")
    msg(f"[INFO] outdir      : {outdir}")

    genes_df = read_gene_table(gene_file).head(max(1, args.top_n)).copy()
    priority_df = harmonize_priority_table(safe_read_table(priority_table)) if priority_table else pd.DataFrame()

    records = []
    for i, row in genes_df.iterrows():
        gene = row["gene"]
        msg(f"[QUERY {i+1}/{len(genes_df)}] {gene}")
        u = query_uniprot(gene, cache=cache, timeout=args.timeout)
        time.sleep(args.sleep)
        r = query_rcsb(gene, u.get("uniprot_id", ""), cache=cache, timeout=args.timeout)
        time.sleep(args.sleep)
        pdb_ids = [x for x in str(r.get("pdb_ids", "")).split(";") if x]
        rd = query_rcsb_entry_details(pdb_ids, cache=cache, timeout=args.timeout, max_entries=args.max_pdb_entries)
        time.sleep(args.sleep)
        af = query_alphafold(u.get("uniprot_id", ""), u.get("alphafold_xref", ""), cache=cache, timeout=args.timeout)
        time.sleep(args.sleep)
        rec = {
            **row.to_dict(),
            **u,
            **r,
            **rd,
            **af,
        }
        rec["target_class"] = infer_target_class(gene, rec.get("protein_name", ""), rec.get("keywords", ""), rec.get("subcellular_location_note", ""))
        rec["subcellular_class"] = infer_subcellular_class(rec.get("subcellular_location_note", ""))
        records.append(rec)

    struct_df = pd.DataFrame(records)
    merged = struct_df.merge(priority_df, on="gene", how="left", suffixes=("", "_priority"))
    ranked = score_structural_priority(merged)

    evidence_cols = [
        "gene", "direction", "target_class", "subcellular_class", "uniprot_id", "protein_name", "reviewed",
        "sequence_length", "has_pdb", "n_pdb_entries", "best_resolution", "median_resolution",
        "experimental_methods", "ligand_bound_entries", "total_nonpolymer_instances",
        "has_alphafold", "structure_source", "docking_readiness", "docking_recommendation",
        "priority_score", "actionability_score", "n_compounds", "n_clinical_compounds", "best_phase",
        "structural_evidence_score", "ligandability_proxy_score", "docking_readiness_score",
        "transcriptomic_support_score", "translational_priority_score", "structure_priority_tier", "priority_rank",
    ]
    evidence_matrix = ranked[[c for c in evidence_cols if c in ranked.columns]].copy()
    docking_ready = ranked.loc[ranked["docking_readiness"].isin(["Experimental-ready", "Structure-available"])].copy()
    shortlist = ranked.head(min(12, len(ranked))).copy()

    target_class_summary = ranked.groupby("target_class", dropna=False).agg(
        n=("gene", "count"),
        mean_priority=("translational_priority_score", "mean"),
        mean_structural=("structural_evidence_score", "mean"),
        mean_ligandability=("ligandability_proxy_score", "mean"),
    ).reset_index().sort_values("mean_priority", ascending=False)

    method_summary = ranked.groupby("dominant_method", dropna=False).agg(
        n=("gene", "count"),
        mean_priority=("translational_priority_score", "mean"),
        median_best_resolution=("best_resolution", "median"),
    ).reset_index().sort_values("mean_priority", ascending=False)

    resolution_summary = ranked.groupby("docking_readiness", dropna=False).agg(
        n=("gene", "count"),
        mean_resolution=("best_resolution", "mean"),
        mean_priority=("translational_priority_score", "mean"),
    ).reset_index()

    reviewed_summary = ranked.assign(review_status=np.where(ranked["reviewed"], "Reviewed (Swiss-Prot)", "Unreviewed")).groupby(
        "review_status", dropna=False
    ).agg(n=("gene", "count"), mean_priority=("translational_priority_score", "mean")).reset_index()

    structure_source_summary = ranked.groupby("structure_source", dropna=False).agg(
        n=("gene", "count"),
        mean_priority=("translational_priority_score", "mean"),
    ).reset_index()

    tier_summary = ranked.groupby("structure_priority_tier", dropna=False).agg(
        n=("gene", "count"),
        mean_priority=("translational_priority_score", "mean"),
    ).reset_index().sort_values("structure_priority_tier")

    direction_readiness = pd.crosstab(ranked["direction"], ranked["docking_readiness"]).reset_index()
    plot_diagnostics = pd.DataFrame({
        "metric": [
            "n_ranked_targets",
            "n_targets_without_experimental_structure",
            "n_targets_with_zero_structural_evidence",
            "n_experimental_ready",
            "n_alphafold_only",
            "n_no_structure",
            "max_translational_priority",
            "median_translational_priority",
        ],
        "value": [
            int(len(ranked)),
            int((~ranked["has_pdb"].astype(bool)).sum()),
            int((pd.to_numeric(ranked["structural_evidence_score"], errors="coerce").fillna(0) <= 1e-9).sum()),
            int((ranked["docking_readiness"] == "Experimental-ready").sum()),
            int((ranked["docking_readiness"] == "AlphaFold-only").sum()),
            int((ranked["docking_readiness"] == "No structure").sum()),
            float(pd.to_numeric(ranked["translational_priority_score"], errors="coerce").max()),
            float(pd.to_numeric(ranked["translational_priority_score"], errors="coerce").median()),
        ],
    })

    subcellular_summary = ranked.groupby("subcellular_class", dropna=False).agg(
        n=("gene", "count"),
        mean_priority=("translational_priority_score", "mean"),
    ).reset_index().sort_values("mean_priority", ascending=False)

    source_manifest = pd.DataFrame({
        "gene": ranked["gene"],
        "uniprot_id": ranked["uniprot_id"],
        "pdb_ids": ranked["pdb_ids"],
        "alphafold_model_id": ranked["alphafold_model_id"],
        "structure_source": ranked["structure_source"],
    })

    table_map = {
        "table_structural_target_prioritization.tsv": ranked,
        "table_docking_ready_targets.tsv": docking_ready,
        "table_structure_evidence_matrix.tsv": evidence_matrix,
        "table_target_class_summary.tsv": target_class_summary,
        "table_experimental_method_summary.tsv": method_summary,
        "table_priority_shortlist.tsv": shortlist,
        "table_resolution_summary.tsv": resolution_summary,
        "table_reviewed_summary.tsv": reviewed_summary,
        "table_structure_source_summary.tsv": structure_source_summary,
        "table_tier_summary.tsv": tier_summary,
        "table_direction_readiness.tsv": direction_readiness,
        "table_plot_diagnostics.tsv": plot_diagnostics,
        "table_subcellular_summary.tsv": subcellular_summary,
        "table_structure_source_manifest.tsv": source_manifest,
    }
    for fname, df in table_map.items():
        df.to_csv(outpaths["tables"] / fname, sep="\t", index=False)

    write_excel_bundle(
        outpaths["tables"] / "Structural_Target_Prioritization_Q1_Enhanced.xlsx",
        {
            "ranked_targets": ranked,
            "evidence_matrix": evidence_matrix,
            "docking_ready": docking_ready,
            "shortlist": shortlist,
            "target_class_summary": target_class_summary,
            "method_summary": method_summary,
            "resolution_summary": resolution_summary,
            "reviewed_summary": reviewed_summary,
            "source_summary": structure_source_summary,
            "tier_summary": tier_summary,
            "plot_diagnostics": plot_diagnostics,
            "subcellular_summary": subcellular_summary,
        },
    )

    plot_priority_bubble(ranked, outpaths["figures"] / "Fig_S1_priority_bubble")
    plot_score_scatter(ranked, outpaths["figures"] / "Fig_S2_priority_scatter")
    plot_evidence_heatmap(ranked, outpaths["figures"] / "Fig_S3_evidence_heatmap")
    plot_readiness_bar(ranked, outpaths["figures"] / "Fig_S4_docking_readiness")
    plot_method_priority(ranked, outpaths["figures"] / "Fig_S5_method_priority")
    plot_direction_readiness_heatmap(ranked, outpaths["figures"] / "Fig_S6_direction_readiness_heatmap")
    plot_shortlist_lollipop(ranked, outpaths["figures"] / "Fig_S7_shortlist_lollipop")
    plot_tier_distribution(ranked, outpaths["figures"] / "Fig_S8_tier_distribution")
    plot_structure_source_summary(ranked, outpaths["figures"] / "Fig_S9_structure_source_summary")
    plot_reviewed_summary(ranked, outpaths["figures"] / "Fig_S10_reviewed_summary")

    run_meta = {
        "project_root": root,
        "manifest": manifest or "",
        "gene_file": gene_file,
        "priority_table": priority_table,
        "n_targets": int(len(ranked)),
        "timestamp_unix": time.time(),
        "script": Path(__file__).name,
        "seed": SEED,
        "cache_enabled": (not args.disable_cache),
        "cache_ttl_hours": args.cache_ttl_hours,
    }
    (outpaths["logs"] / "run_metadata.json").write_text(json.dumps(run_meta, indent=2), encoding="utf-8")
    msg("[DONE] Enhanced structural target prioritization analysis completed.")


if __name__ == "__main__":
    main()
