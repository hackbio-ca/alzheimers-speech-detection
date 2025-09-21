#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a merged ML table from Tahoe-100M pseudobulk DE + SMILES features.

What it does
------------
1) Streams rows LIVE from Hugging Face (no local dataset needed):
      dataset='tahoebio/Tahoe-100M', subset='pseudobulk_differential_expression'
2) Normalizes core fields to a common schema:
      gene, drug, cell, dose, logFC, pval, qval, stat
3) Attaches SMILES to Tahoe rows in THREE passes:
      (1) explicit drug-map text file (DrugName,SMILES per line)
      (2) fallback from features CSV (drug_name → smiles)
      (3) optional PubChem lookup by drug name (cached for reproducibility)
4) Canonicalizes SMILES (RDKit if available; safe lowercase fallback)
5) Merges Tahoe rows with Noah's SMILES feature table on 'smiles'
6) Saves CSV (always) + Parquet (if pyarrow/fastparquet installed)
7) Prints helpful previews + join diagnostics

Quickstart (Windows paths shown as example)
-------------------------------------------
# Place the features CSV in your repo root OR pass the full path
python make_merged_table.py ^
  --n 1000 ^
  --features "C:\\Users\\you\\repo\\drug_descriptors_scaled.csv" ^
  --drug-map "C:\\Users\\you\\repo\\unique_drugs_tahoe100m.txt" ^
  --use-pubchem ^
  --out-prefix merged_1k

Then scale up:
  --n 10000 --out-prefix merged_10k
  --n 100000 --out-prefix merged_100k

Environment
-----------
- Recommended: conda env with:
    rdkit (conda-forge), datasets, pandas, numpy, (optional) pyarrow
- Optional enforcement: set EXPECTED_CONDA_ENV to your env name (default 'hackbio')

Notes
-----
- We intentionally avoid early dedup to keep as much signal as possible.
- We dedup only at the end using all key identifiers plus SMILES.
- PubChem is only queried for rows that STILL miss SMILES after local sources.
- A JSON cache ensures reproducibility and avoids repeated lookups.
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Iterable, Any
from urllib.parse import quote

import pandas as pd
import numpy as np

# -----------------------------
# Optional: RDKit for SMILES canonicalization
# -----------------------------
try:
    from rdkit import Chem
    RDKit_OK = True
except Exception:
    RDKit_OK = False

# -----------------------------
# Required: Hugging Face datasets (streaming)
# -----------------------------
try:
    from datasets import load_dataset
except Exception as e:
    print("=" * 60)
    print("IMPORT ERROR: Missing dependency 'datasets'.")
    print("Install inside your env:")
    print("  pip install datasets   (or: conda install -c conda-forge datasets)")
    print("=" * 60)
    print("Detail:", e)
    sys.exit(1)

# -----------------------------
# Optional: requests for PubChem lookups
# -----------------------------
try:
    import requests
    REQUESTS_OK = True
except Exception:
    REQUESTS_OK = False

# -----------------------------
# Environment banner (warning-only)
# -----------------------------
EXPECTED_ENV = os.environ.get("EXPECTED_CONDA_ENV", "hackbio")

def check_environment() -> None:
    """Warn (non-fatal) if not running in the expected conda environment."""
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    if conda_env != EXPECTED_ENV:
        print("=" * 60)
        print("ENVIRONMENT WARNING")
        print("=" * 60)
        print(f"Current env: '{conda_env or 'None'}' | Expected: '{EXPECTED_ENV}'")
        print("Continuing anyway… (set EXPECTED_CONDA_ENV to silence this)")
        print("=" * 60)

# -----------------------------
# Helpers: normalization & SMILES utilities
# -----------------------------
SALT_TOKENS = {
    "hydrochloride","hydrobromide","nitrate","nitrite","sulfate","sulphate",
    "phosphate","acetate","succinate","tartrate","fumarate","maleate",
    "mesylate","tosylate","bitartrate","citrate","oxalate","lactate",
    "sodium","potassium","calcium","magnesium","zinc","ammonium","hydrate"
}

def canonicalize_smiles(s: Optional[str]) -> Optional[str]:
    """Canonicalize SMILES via RDKit (if available); else lowercase/strip."""
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None
    if RDKit_OK:
        try:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            return None
    return s.lower()

def clean_drug_name(name: Optional[str]) -> Optional[str]:
    """
    Normalize drug names for better matching across sources:
      - lowercase
      - remove parenthetical suffixes: 'verapamil (hydrochloride)' → 'verapamil'
      - replace separators with spaces, normalize whitespace
      - drop trailing salt/counter-ion words
      - convert slashes/hyphens to spaces: 'AZD6244/MEK inhibitor' → 'azd6244 mek inhibitor'
    """
    if name is None:
        return None
    s = str(name).lower().strip()

    # remove parenthetical pieces
    while "(" in s and ")" in s and s.index("(") < s.index(")"):
        i, j = s.index("("), s.index(")")
        s = (s[:i] + " " + s[j+1:]).strip()

    # normalize separators and spaces
    for sep in ["/", "-", "_"]:
        s = s.replace(sep, " ")
    s = " ".join(s.split())

    # remove trailing salts/counter-ions
    parts = s.split()
    while parts and parts[-1] in SALT_TOKENS:
        parts = parts[:-1]
    s = " ".join(parts)

    return s or None

def first_non_null(d: Dict, keys: List[str], default=None):
    """Return the first present, non-null value among candidate keys."""
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return default

def normalize_tahoe_row(row: Dict) -> Dict:
    """
    Normalize a raw Tahoe row into our schema.
    We include multiple candidate keys per field since schema names can vary.
    """
    gene = first_non_null(row, [
        "gene","gene_symbol","gene_name","symbol","target_gene","hgnc_symbol",
        "ensembl_gene","ensembl_id","gene_id"
    ])
    drug = first_non_null(row, [
        "drug","pert_name","compound","drug_name","pert_iname","compound_name",
        "treatment","perturbation"
    ])
    cell = first_non_null(row, [
        "cell","cell_type","cellline","cell_line","cell_id","cell_name",
        "line","line_name"
    ])
    dose = first_non_null(row, [
        "dose","dose_um","dose_uM","dose_umol","dose_uM_final","concentration_uM",
        "concentration","dose_value"
    ])

    logfc = first_non_null(row, [
        "logFC","log2_fc","log2FoldChange","lfc","log2fc","log2_fc_mean"
    ])
    pval  = first_non_null(row, [
        "p_value","pval","pvalue","pval_raw","p_value_raw"
    ])
    qval  = first_non_null(row, [
        "q_value","padj","qval","adj_p_value","fdr_qvalue","qvalue"
    ])
    stat  = first_non_null(row, [
        "stat","t_stat","tstat","z_score","wald_stat","score"
    ])

    return {
        "gene": gene, "drug": drug, "cell": cell, "dose": dose,
        "logFC": logfc, "pval": pval, "qval": qval, "stat": stat,
    }

# -----------------------------
# Streaming Tahoe rows (no early dedup)
# -----------------------------
def stream_tahoe(n: int,
                 dataset: str = "tahoebio/Tahoe-100M",
                 subset: str = "pseudobulk_differential_expression",
                 split: str = "train",
                 streaming: bool = True,
                 target_unique_drugs: Optional[int] = None,
                 shuffle_buffer: int = 50_000,
                 seed: int = 42) -> pd.DataFrame:
    """
    Stream up to N rows, shuffling online. If target_unique_drugs is set,
    stop early once we have that many distinct 'drug_clean' names.
    """
    print(f"Connecting HF → {dataset}/{subset} (split={split}, streaming={streaming})")
    ds = load_dataset(dataset, subset, split=split, streaming=streaming)

    # Shuffle the iterator (supported by IterableDataset)
    try:
        ds = ds.shuffle(buffer_size=shuffle_buffer, seed=seed)
        print(f"Shuffling stream with buffer_size={shuffle_buffer}, seed={seed}")
    except Exception:
        print("Note: shuffle not available; proceeding without it.")

    rows = []
    seen = set()
    for i, row in enumerate(ds, start=1):
        norm = normalize_tahoe_row(row)
        if norm.get("drug") is None:
            continue
        norm["drug_clean"] = clean_drug_name(norm["drug"])
        if norm["drug_clean"] is None:
            continue

        rows.append(norm)
        if target_unique_drugs:
            seen.add(norm["drug_clean"])
            if len(seen) >= target_unique_drugs:
                break
        if len(rows) >= n:
            break

    df = pd.DataFrame(rows)
    if df.empty:
        print("No rows retrieved from Tahoe.")
        return df

    print(f"Collected {len(df)} rows; unique drugs: {df['drug_clean'].nunique()}")
    return df


# -----------------------------
# Load feature table & drug map
# -----------------------------
def load_features_csv(path: str) -> pd.DataFrame:
    """
    Features CSV must include 'smiles'; ideally also 'drug_name' to help fallback attach.
    We canonicalize 'smiles' and precompute 'drug_clean' if drug_name exists.
    """
    p = Path(path)
    if not p.exists():
        print(f"ERROR: features CSV not found → {path}")
        sys.exit(1)

    fdf = pd.read_csv(p)
    if "smiles" not in fdf.columns:
        print("ERROR: features CSV must include a 'smiles' column.")
        sys.exit(1)

    if "drug_name" in fdf.columns:
        fdf["drug_clean"] = fdf["drug_name"].map(clean_drug_name)
    else:
        fdf["drug_clean"] = pd.NA

    fdf["smiles"] = fdf["smiles"].apply(canonicalize_smiles)
    fdf = fdf.dropna(subset=["smiles"]).drop_duplicates(subset=["smiles"]).reset_index(drop=True)

    num_features = len([c for c in fdf.columns if c not in {"drug_name", "drug_clean", "smiles"}])
    print(f"Loaded features: {len(fdf)} rows; ~{num_features} feature columns (+ 'smiles').")
    return fdf

def load_drug_map(path: Optional[str]) -> Optional[pd.DataFrame]:
    """
    Load a simple text file with lines: DrugName,SMILES
    Returns DataFrame with ['drug_clean','smiles'] canonicalized.
    """
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        print(f"WARNING: drug map not found → {path} (continuing without it).")
        return None

    names, smiles = [], []
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if len(parts) < 2:
                continue
            names.append(clean_drug_name(parts[0]))
            smiles.append(canonicalize_smiles(parts[1]))

    dm = pd.DataFrame({"drug_clean": names, "smiles": smiles}).dropna().drop_duplicates()
    print(f"Loaded drug map: {len(dm)} unique drug→SMILES pairs.")
    return dm

# -----------------------------
# PubChem lookup (optional, cached)
# -----------------------------
def load_pubchem_cache(path: Optional[str]) -> Dict[str, Optional[str]]:
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}

def save_pubchem_cache(cache: Dict[str, Optional[str]], path: Optional[str]) -> None:
    if not path:
        return
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)
    except Exception:
        pass

def pubchem_smiles_from_name(name: str, timeout: float = 8.0) -> Optional[str]:
    """
    Query PubChem PUG REST by NAME → CanonicalSMILES.
    Tries variants (original, cleaned, quoted). Returns canonicalized SMILES or None.
    """
    if not REQUESTS_OK or not name:
        return None

    variants = [name]
    cleaned = clean_drug_name(name)
    if cleaned and cleaned != name:
        variants.append(cleaned)
    variants.extend([f'"{v}"' for v in list(variants)])  # quoted variants

    for q in variants:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(q)}/property/CanonicalSMILES/JSON"
            r = requests.get(url, timeout=timeout)
            if r.status_code != 200:
                continue
            data = r.json()
            props = data.get("PropertyTable", {}).get("Properties", [])
            if props:
                smi = props[0].get("CanonicalSMILES")
                return canonicalize_smiles(smi)
        except Exception:
            continue
    return None

# -----------------------------
# Attach SMILES to Tahoe rows: map → features → PubChem (optional)
# -----------------------------
def attach_smiles_to_tahoe(tahoe: pd.DataFrame,
                           drug_map: Optional[pd.DataFrame],
                           feats: pd.DataFrame,
                           use_pubchem: bool = False,
                           pubchem_cache_path: Optional[str] = None,
                           pubchem_delay_ms: int = 120) -> pd.DataFrame:
    df = tahoe.copy()
    df["smiles"] = pd.NA

    # Pass 1: explicit mapping on 'drug_clean'
    if drug_map is not None and not drug_map.empty:
        df = df.merge(drug_map, on="drug_clean", how="left", suffixes=("", "_map"))
        if "smiles_map" in df.columns:
            df["smiles"] = df["smiles"].fillna(df["smiles_map"])
            df = df.drop(columns=["smiles_map"])

    # Pass 2: features fallback (drug_clean → smiles)
    if df["smiles"].isna().any() and "drug_clean" in feats.columns:
        fb = feats[["drug_clean", "smiles"]].dropna().drop_duplicates()
        df = df.merge(fb, on="drug_clean", how="left", suffixes=("", "_feat"))
        if "smiles_feat" in df.columns:
            df["smiles"] = df["smiles"].fillna(df["smiles_feat"])
            df = df.drop(columns=["smiles_feat"])

    # Pass 3: PubChem lookups for remaining gaps (optional)
    if use_pubchem and df["smiles"].isna().any():
        cache = load_pubchem_cache(pubchem_cache_path)
        delay = max(0, pubchem_delay_ms) / 1000.0
        missing_idxs = list(df.index[df["smiles"].isna()])
        resolved = 0

        for idx in missing_idxs:
            raw = str(df.at[idx, "drug"] or "").strip()
            if not raw:
                continue
            key = raw.lower()
            if key in cache:
                smi = cache[key]
            else:
                smi = pubchem_smiles_from_name(raw)
                cache[key] = smi
                if delay > 0:
                    time.sleep(delay)
            if smi:
                df.at[idx, "smiles"] = smi
                resolved += 1

        if pubchem_cache_path:
            save_pubchem_cache(cache, pubchem_cache_path)
        print(f"PubChem resolved {resolved} additional SMILES.")

    # Final clean-up
    df["smiles"] = df["smiles"].apply(canonicalize_smiles)
    got = df["smiles"].notna().sum()
    print(f"Attached SMILES for {got}/{len(df)} Tahoe rows (after all passes).")

    # Keep rows with SMILES; dedup late using full identity
    df = df.dropna(subset=["smiles"]).drop_duplicates(subset=["gene","drug","cell","dose","smiles"])
    print(f"Usable Tahoe rows after attach: {len(df)}")
    return df

# -----------------------------
# Merge & Save
# -----------------------------
def merge_on_smiles(tahoe_df: pd.DataFrame, feats_df: pd.DataFrame) -> pd.DataFrame:
    merged = tahoe_df.merge(feats_df, on="smiles", how="inner")
    print(f"Merged rows: {len(merged)} (inner join on 'smiles').")
    return merged

def save_outputs(df: pd.DataFrame, out_prefix: str, save_schema: bool = True) -> None:
    out_dir = Path.cwd()
    csv_path = out_dir / f"{out_prefix}.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV     → {csv_path}")

    try:
        pq_path = out_dir / f"{out_prefix}.parquet"
        df.to_parquet(pq_path, index=False)
        print(f"Saved Parquet → {pq_path}")
    except Exception as e:
        print("Parquet save skipped (install 'pyarrow' or 'fastparquet' to enable). Detail:", e)

    if save_schema:
        schema = pd.DataFrame({
            "column": df.columns,
            "dtype": df.dtypes.astype(str),
            "non_null": df.notnull().sum().values
        })
        schema_path = out_dir / f"{out_prefix}__schema.csv"
        schema.to_csv(schema_path, index=False)
        print(f"Saved schema  → {schema_path}")

# -----------------------------
# CLI
# -----------------------------
def parse_args():
    SCRIPT_DIR = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="Merge Tahoe pseudobulk DE (streamed) with SMILES features via map→features→PubChem.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--n", type=int, default=100, help="Rows to stream from Tahoe (default: 100).")
    p.add_argument("--features", default=str(SCRIPT_DIR / "drug_descriptors_scaled.csv"),
                   help="Path to features CSV (must include 'smiles'). Default: ./drug_descriptors_scaled.csv")
    p.add_argument("--drug-map", default=None, help="Optional drug→SMILES map text file (DrugName,SMILES per line).")
    p.add_argument("--dataset", default="tahoebio/Tahoe-100M")
    p.add_argument("--subset", default="pseudobulk_differential_expression")
    p.add_argument("--split", default="train")
    p.add_argument("--no-stream", action="store_true", help="Disable streaming mode (downloads the split).")
    p.add_argument("--use-pubchem", action="store_true",
                   help="Query PubChem for missing SMILES after local sources.")
    p.add_argument("--pubchem-cache", default="pubchem_smiles_cache.json",
                   help="JSON file to cache name→SMILES lookups (default: pubchem_smiles_cache.json).")
    p.add_argument("--pubchem-delay-ms", type=int, default=120,
                   help="Delay between PubChem requests in milliseconds (default: 120).")
    p.add_argument("--out-prefix", default=None, help="Output filename prefix (default: merged_<n>).")
    p.add_argument("--debug", action="store_true", help="Print extra schema/preview info.")
    p.add_argument("--unique-drugs", type=int, default=None,
                   help="Stop streaming once this many distinct drugs are seen (in addition to --n cap).")
    p.add_argument("--shuffle-buffer", type=int, default=50000,
                   help="Buffer size for streaming shuffle (default: 50k).")
    p.add_argument("--seed", type=int, default=42, help="Shuffle seed (default: 42).")

    return p.parse_args()

# -----------------------------
# Main
# -----------------------------
def main():
    check_environment()
    args = parse_args()
    out_prefix = args.out_prefix or f"merged_{args.n}"
    streaming = not args.no_stream

    # 1) Stream Tahoe
    tahoe = stream_tahoe(
        args.n, args.dataset, args.subset, args.split, streaming,
        target_unique_drugs=args.unique_drugs,
        shuffle_buffer=args.shuffle_buffer,
        seed=args.seed
    )
    if tahoe.empty:
        print("No Tahoe rows collected. Exiting.")
        sys.exit(1)

    # (Optional) drop rows without a target BEFORE merging to avoid big drops later
    if "logFC" in tahoe.columns:
        before = len(tahoe)
        tahoe = tahoe.dropna(subset=["logFC"])
        print(f"Filtered logFC NaNs in Tahoe: dropped {before - len(tahoe)} rows; remaining {len(tahoe)}")

    if args.debug:
        # Peek at actual raw schema keys (first record)
        try:
            sample_iter = load_dataset(args.dataset, args.subset, split=args.split, streaming=True)
            first_row = next(iter(sample_iter))
            if isinstance(first_row, dict):  # makes PyCharm happy
                keys_preview = sorted(first_row.keys())
                print("\nFirst streamed record keys (schema preview, truncated to 60):")
                print(keys_preview[:60])
        except Exception:
            pass

        print("\nSample raw drug names from Tahoe:")
        print(tahoe["drug"].dropna().head(10).to_list())

        print("\nTahoe normalized preview:")
        print(tahoe.head(5).to_string(index=False))

    # 2) Load features
    print(f"\nLoading features: {args.features}")
    feats = load_features_csv(args.features)

    # 3) Load mapping + attach SMILES (map → features → PubChem)
    dm = load_drug_map(args.drug_map)
    tahoe_smiles = attach_smiles_to_tahoe(
        tahoe,
        dm,
        feats,
        use_pubchem=args.use_pubchem,
        pubchem_cache_path=args.pubchem_cache,
        pubchem_delay_ms=args.pubchem_delay_ms
    )
    if tahoe_smiles.empty:
        print("Could not attach SMILES to Tahoe rows. Check mapping/cleaning or enable --use-pubchem.")
        sys.exit(1)

    # 4) Merge on SMILES
    print("\nMerging on 'smiles' …")
    merged = merge_on_smiles(tahoe_smiles, feats)
    if merged.empty:
        print("Merge returned 0 rows. Check name normalization/mapping.")
        sys.exit(1)

    # Optional: ensure numeric target (logFC) if present
    if "logFC" in merged.columns and np.issubdtype(merged["logFC"].dtype, np.number):
        before = len(merged)
        merged = merged.dropna(subset=["logFC"])
        dropped = before - len(merged)
        if dropped > 0:
            print(f"Dropped {dropped} rows with NaN logFC.")

    # Optional: final dedup with full identity
    merged = merged.drop_duplicates(subset=["gene","drug","cell","dose","smiles"])

    # 5) Save outputs
    print("\nSaving outputs …")
    save_outputs(merged, out_prefix, save_schema=True)
    print(f"Merged shape: {merged.shape}")

    # Preview to console
    print("\nMerged preview (first 5 rows):")
    print(merged.head().to_string(index=False))

    # Small summary
    id_cols = ["gene", "drug", "cell", "dose"]
    uniq = {c: merged[c].nunique(dropna=True) for c in id_cols if c in merged.columns}
    feat_cols = [c for c in merged.columns if c not in {"gene","drug","cell","dose","smiles","logFC","pval","qval","stat","drug_name","drug_clean"}]
    print("\nSummary:")
    for k, v in uniq.items():
        print(f"  unique {k}: {v}")
    print(f"  feature columns merged: ~{len(feat_cols)}")
    print("\nDone.")

if __name__ == "__main__":
    main()
