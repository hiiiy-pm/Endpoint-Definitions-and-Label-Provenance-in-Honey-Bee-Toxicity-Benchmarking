"""Shared utilities for the revised ApisTox endpoint-decomposition analysis.

The module keeps all data definitions, split reconstruction, molecular representations,
and resampling helpers in one auditable place. No test labels are used in fitting,
feature scaling, probability calibration, or model combination.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import os
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from rdkit import Chem, DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

SEED = 20260829
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW_CSV = DATA / "raw" / "apistox.csv"
SPLIT_DIR = DATA / "official_splits"
RESULTS = ROOT / "results"
CACHE = RESULTS / "cache"
CACHE.mkdir(parents=True, exist_ok=True)

SPLIT_KEYS = ("random", "maxmin", "time")
SPLIT_NAMES = {"random": "Random", "maxmin": "MaxMin", "time": "Time"}
THRESHOLDS = (100, 11, 1)
REPRESENTATIONS = ("ECFP", "Avalon", "MACCS", "WL-HI", "Descriptors12")

DESC_NAMES = [
    "MolWt", "MolLogP", "TPSA", "HDonors", "HAcceptors", "RotBonds",
    "RingCount", "AromaticRings", "FractionCSP3", "HeavyAtoms", "MolMR",
    "FormalCharge",
]

META_CATEGORICAL = ["source", "toxicity_type"]
META_BINARY = ["herbicide", "fungicide", "insecticide", "other_agrochemical"]


@dataclass(frozen=True)
class StudyData:
    df: pd.DataFrame
    tier: np.ndarray
    y: Dict[int, np.ndarray]
    splits: Dict[str, Tuple[np.ndarray, np.ndarray]]


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            block = f.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def derive_tier(df: pd.DataFrame) -> np.ndarray:
    tier = np.select(
        [
            df["ppdb_level"].eq(0),
            df["ppdb_level"].eq(1) & df["label"].eq(0),
            df["ppdb_level"].eq(1) & df["label"].eq(1),
            df["ppdb_level"].eq(2),
        ],
        [0, 1, 2, 3],
        default=-1,
    ).astype(np.int8)
    if np.any(tier < 0):
        bad = df.loc[tier < 0, ["label", "ppdb_level"]].drop_duplicates()
        raise ValueError(f"Unmapped label combinations found:\n{bad}")
    return tier


def load_study_data(raw_csv: Path = RAW_CSV, split_dir: Path = SPLIT_DIR) -> StudyData:
    data_path = raw_csv
    df = pd.read_csv(data_path)
    required = {
        "name", "CID", "CAS", "SMILES", "source", "year", "toxicity_type",
        "herbicide", "fungicide", "insecticide", "other_agrochemical",
        "label", "ppdb_level",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if len(df) != 1035 or df[["CID", "CAS", "SMILES"]].isna().any().any():
        raise ValueError("The frozen analysis expects 1,035 rows with complete compound identities.")
    if df["SMILES"].duplicated().any():
        raise ValueError("The frozen analysis expects unique SMILES rows.")

    tier = derive_tier(df)
    y = {
        100: (tier >= 1).astype(np.int8),
        11: (tier >= 2).astype(np.int8),
        1: (tier >= 3).astype(np.int8),
    }
    if not np.all(y[1] <= y[11]) or not np.all(y[11] <= y[100]):
        raise AssertionError("Nested endpoint invariant failed.")

    sm2idx = {s: i for i, s in enumerate(df["SMILES"])}
    splits: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for key in SPLIT_KEYS:
        trf = pd.read_csv(split_dir / f"{key}_train.csv")
        tef = pd.read_csv(split_dir / f"{key}_test.csv")
        try:
            tr = np.array([sm2idx[s] for s in trf["SMILES"]], dtype=np.int32)
            te = np.array([sm2idx[s] for s in tef["SMILES"]], dtype=np.int32)
        except KeyError as e:
            raise ValueError(f"Split {key} contains an unknown SMILES: {e}") from e
        if len(tr) != 828 or len(te) != 207 or len(set(tr).intersection(te)):
            raise AssertionError(f"Unexpected official split geometry for {key}")
        if len(np.unique(tr)) != len(tr) or len(np.unique(te)) != len(te):
            raise AssertionError(f"Duplicate compound within official split {key}")
        if set(tr).union(te) != set(range(len(df))):
            raise AssertionError(f"Official split {key} does not cover the complete raw cohort")
        splits[SPLIT_NAMES[key]] = (tr, te)

    return StudyData(df=df, tier=tier, y=y, splits=splits)


def scaffold_key(mol: Chem.Mol) -> str:
    smi = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if smi:
        return smi
    # Acyclic molecules must not all collapse into one artificial empty scaffold.
    return "ACYCLIC:" + Chem.MolToSmiles(mol, canonical=True)


def _tanimoto_matrix(fps: Sequence[DataStructs.ExplicitBitVect]) -> np.ndarray:
    n = len(fps)
    out = np.empty((n, n), dtype=np.float32)
    for i, fp in enumerate(fps):
        out[i] = DataStructs.BulkTanimotoSimilarity(fp, fps)
    return out


def _wlhi_kernel(mols: Sequence[Chem.Mol]) -> np.ndarray:
    labels = [
        [f"{a.GetAtomicNum()}_{a.GetFormalCharge()}_{int(a.GetIsAromatic())}" for a in m.GetAtoms()]
        for m in mols
    ]
    molecule_features = [{} for _ in mols]
    for iteration in range(3):
        for mi, labs in enumerate(labels):
            d = molecule_features[mi]
            for lab in labs:
                key = f"{iteration}:{lab}"
                d[key] = d.get(key, 0) + 1
        if iteration < 2:
            nxt = []
            for mol, labs in zip(mols, labels):
                new_labs = []
                for atom in mol.GetAtoms():
                    neigh = []
                    for bond in atom.GetBonds():
                        other = bond.GetOtherAtom(atom).GetIdx()
                        neigh.append(f"{bond.GetBondType()}:{labs[other]}")
                    signature = labs[atom.GetIdx()] + "|" + "|".join(sorted(neigh))
                    new_labs.append(hashlib.sha1(signature.encode("utf-8")).hexdigest()[:16])
                nxt.append(new_labs)
            labels = nxt

    feature_ids = {}
    rows, cols, vals = [], [], []
    for i, d in enumerate(molecule_features):
        for key, count in d.items():
            # Histogram intersection feature expansion used by the frozen round-2 code.
            for t in range(1, int(count) + 1):
                expanded = (key, t)
                j = feature_ids.setdefault(expanded, len(feature_ids))
                rows.append(i)
                cols.append(j)
                vals.append(1.0)
    x = sparse.csr_matrix(
        (vals, (rows, cols)), shape=(len(mols), len(feature_ids)), dtype=np.float32
    )
    k = (x @ x.T).toarray().astype(np.float32)
    diag = np.sqrt(np.maximum(np.diag(k), 1e-12))
    k /= diag[:, None] * diag[None, :]
    return k


def compute_representations(study: StudyData) -> Mapping[str, object]:
    """Recompute every representation from raw ordered SMILES, without cache reads."""
    mols = [Chem.MolFromSmiles(s) for s in study.df["SMILES"]]
    if not all(mols):
        raise ValueError("At least one SMILES could not be parsed by RDKit.")

    morgan = GetMorganGenerator(radius=2, fpSize=1024)
    ecfp = [morgan.GetFingerprint(m) for m in mols]
    avalon = [pyAvalonTools.GetAvalonFP(m, nBits=1024) for m in mols]
    maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]

    kernels = {
        "ECFP": _tanimoto_matrix(ecfp),
        "Avalon": _tanimoto_matrix(avalon),
        "MACCS": _tanimoto_matrix(maccs),
        "WL-HI": _wlhi_kernel(mols),
    }
    funcs = [
        Descriptors.MolWt,
        Descriptors.MolLogP,
        Descriptors.TPSA,
        Descriptors.NumHDonors,
        Descriptors.NumHAcceptors,
        Descriptors.NumRotatableBonds,
        Descriptors.RingCount,
        Descriptors.NumAromaticRings,
        Descriptors.FractionCSP3,
        Descriptors.HeavyAtomCount,
        Descriptors.MolMR,
        lambda m: Chem.GetFormalCharge(m),
    ]
    descriptors = np.array([[f(m) for f in funcs] for m in mols], dtype=np.float64)
    scaffolds = np.array([scaffold_key(m) for m in mols], dtype=object)
    return {"study": study, "kernels": kernels, "descriptors": descriptors, "scaffolds": scaffolds}


def build_or_load_cache(force: bool = False, cache_dir: Path = CACHE) -> Mapping[str, object]:
    """Load a validated cache, or explicitly regenerate it with recorded provenance.

    Existing partial, stale or unverifiable caches fail rather than silently changing
    the analysis. Legacy caches require a separate current raw-SMILES validation
    record; their unknown original generation environment is never backfilled.
    """
    from _cache_validation import cache_paths, validate_cache, write_generation_metadata

    study = load_study_data()
    expected = cache_paths(cache_dir)
    if not force and (any(p.exists() for p in expected.values())
                      or (cache_dir / "CACHE_METADATA.json").exists()
                      or (cache_dir / "CACHE_VALIDATION.json").exists()):
        return validate_cache(study, cache_dir)

    arrays = compute_representations(study)
    cache_dir.mkdir(parents=True, exist_ok=True)
    # An interrupted rebuild must not leave an old attestation applicable to new files.
    for name in ("CACHE_METADATA.json", "CACHE_VALIDATION.json"):
        (cache_dir / name).unlink(missing_ok=True)
    for rep, kernel in arrays["kernels"].items():
        np.save(expected[rep], kernel)
    np.save(expected["Descriptors12"], arrays["descriptors"])
    scaffolds = arrays["scaffolds"]
    pd.DataFrame({"Index": np.arange(len(scaffolds)), "Scaffold": scaffolds}).to_csv(
        expected["Scaffolds"], index=False
    )

    processed = study.df.copy()
    processed["Tier"] = study.tier
    for thr in THRESHOLDS:
        processed[f"Y{thr}"] = study.y[thr]
    processed.to_csv(expected["Processed"], index=False)

    write_generation_metadata(study, cache_dir)
    return validate_cache(study, cache_dir)


def tie_aware_weighted_auc(y: np.ndarray, score: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted AUROC for one or many weight vectors with half credit for ties.

    Parameters
    ----------
    y : (n,) binary array
    score : (n,) prediction scores
    weights : (b, n) or (n,) non-negative sample weights
    """
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=float)
    w = np.asarray(weights, dtype=float)
    if w.ndim == 1:
        w = w[None, :]
    order = np.argsort(score, kind="mergesort")
    ys = y[order]
    ss = score[order]
    ws = w[:, order]

    pos_total = (ws * (ys == 1)).sum(axis=1)
    neg_total = (ws * (ys == 0)).sum(axis=1)
    den = pos_total * neg_total
    num = np.zeros(len(w), dtype=float)
    cum_neg = np.zeros(len(w), dtype=float)

    starts = np.r_[0, np.flatnonzero(np.diff(ss) != 0) + 1]
    ends = np.r_[starts[1:], len(ss)]
    for start, end in zip(starts, ends):
        block_y = ys[start:end]
        block_w = ws[:, start:end]
        pos_w = (block_w * (block_y == 1)).sum(axis=1)
        neg_w = (block_w * (block_y == 0)).sum(axis=1)
        num += pos_w * (cum_neg + 0.5 * neg_w)
        cum_neg += neg_w

    out = np.full(len(w), np.nan, dtype=float)
    valid = den > 0
    out[valid] = num[valid] / den[valid]
    return out


def multinomial_weights(n_units: int, b: int, rng: np.random.Generator) -> np.ndarray:
    if n_units <= 0:
        raise ValueError("n_units must be positive")
    return rng.multinomial(n_units, np.full(n_units, 1.0 / n_units), size=b).astype(np.int16)


def percentile_ci(values: np.ndarray, alpha: float = 0.05) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return (math.nan, math.nan)
    return (float(np.quantile(values, alpha / 2)), float(np.quantile(values, 1 - alpha / 2)))


def simultaneous_centered_ci(
    observed: np.ndarray, bootstrap: np.ndarray, alpha: float = 0.05
) -> Tuple[np.ndarray, np.ndarray]:
    """Approximate family-wise centered-bootstrap intervals for a contrast family."""
    observed = np.asarray(observed, dtype=float)
    bootstrap = np.asarray(bootstrap, dtype=float)
    deviations = np.abs(bootstrap - observed[None, :])
    max_dev = np.nanmax(deviations, axis=1)
    q = np.nanquantile(max_dev, 1 - alpha)
    return observed - q, observed + q


def entropy_from_neighbor_labels(labels: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
    """Normalized four-class Shannon entropy for one or many label permutations."""
    lab = np.asarray(labels)
    if lab.ndim == 1:
        lab = lab[None, :]
    z = lab[:, neighbors]
    counts = np.stack([(z == c).sum(axis=2) for c in range(4)], axis=2).astype(float)
    p = counts / neighbors.shape[1]
    logp = np.zeros_like(p)
    np.log(p, out=logp, where=p > 0)
    return -(p * logp).sum(axis=2) / np.log(4.0)


def permute_within_groups(
    labels: np.ndarray, groups: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    out = np.asarray(labels).copy()
    groups = np.asarray(groups)
    for g in pd.unique(groups):
        idx = np.flatnonzero(groups == g)
        if len(idx) > 1:
            out[idx] = rng.permutation(out[idx])
    return out


def cramer_v(table: pd.DataFrame) -> float:
    from scipy.stats import chi2_contingency

    arr = np.asarray(table, dtype=float)
    chi2, _, _, _ = chi2_contingency(arr)
    n = arr.sum()
    r, k = arr.shape
    denom = min(k - 1, r - 1)
    return float(np.sqrt((chi2 / n) / denom)) if denom > 0 else math.nan
