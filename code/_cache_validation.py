"""Cache integrity and current recomputation evidence, separate from generation history."""
from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import version
import hashlib
import inspect
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

FILENAMES = {
    "ECFP": "kernel_ecfp.npy", "Avalon": "kernel_avalon.npy",
    "MACCS": "kernel_maccs.npy", "WL-HI": "kernel_wlhi.npy",
    "Descriptors12": "descriptors12.npy", "Scaffolds": "scaffolds.csv",
    "Processed": "processed_data.csv",
}


class CacheValidationError(ValueError):
    """An existing cache cannot safely be used with the current inputs/protocol."""


def cache_paths(directory: Path) -> dict[str, Path]:
    return {key: directory / name for key, name in FILENAMES.items()}


def runtime_environment() -> dict[str, str]:
    return {"python": platform.python_version(),
            **{name: version(name) for name in ("numpy", "pandas", "scipy", "rdkit", "scikit-learn")}}


def _digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def current_context(study) -> dict:
    import _common as c

    inputs = [c.RAW_CSV, *(c.SPLIT_DIR / f"{key}_{part}.csv"
                          for key in c.SPLIT_KEYS for part in ("train", "test"))]
    implementation = "\n".join(inspect.getsource(f) for f in (
        c.compute_representations, c._tanimoto_matrix, c._wlhi_kernel, c.scaffold_key))
    return {
        "input_sha256": {p.relative_to(c.ROOT).as_posix(): c.sha256_file(p) for p in inputs},
        "ordered_smiles_sha256": _digest(study.df.SMILES.tolist()),
        "ordered_identity_sha256": _digest(study.df[["CID", "CAS", "SMILES"]].astype(str).values.tolist()),
        "n_molecules": len(study.df),
        "parameters": {
            "ECFP": {"radius": 2, "bits": 1024, "chirality": False, "kernel": "Tanimoto"},
            "Avalon": {"bits": 1024, "kernel": "Tanimoto"},
            "MACCS": {"bits": 167, "kernel": "Tanimoto"},
            "WL-HI": {"iterations": 3, "normalization": "self-diagonal", "dtype": "float32"},
            "Descriptors12": c.DESC_NAMES,
            "Scaffolds": "Murcko without chirality; canonical full SMILES for acyclic molecules",
        },
        "representation_implementation_sha256": hashlib.sha256(implementation.encode()).hexdigest(),
        "environment": runtime_environment(),
    }


def artifact_hashes(directory: Path) -> dict[str, str]:
    from _common import sha256_file
    return {p.name: sha256_file(p) for p in cache_paths(directory).values()}


def write_generation_metadata(study, directory: Path) -> None:
    from _common import RAW_CSV, SEED, sha256_file
    metadata = {
        "schema_version": 2, "provenance_kind": "generated_from_raw_smiles",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "seed": SEED,
        "source_sha256": sha256_file(RAW_CSV), **current_context(study),
        "artifact_sha256": artifact_hashes(directory),
        "array_shapes": {rep: list(np.load(path, mmap_mode="r").shape)
                         for rep, path in cache_paths(directory).items() if path.suffix == ".npy"},
    }
    (directory / "CACHE_METADATA.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def validate_cache(study, directory: Path, *, require_provenance: bool = True) -> dict:
    """Check current inputs, row identities, numerical invariants and artifact hashes.

    ``require_provenance=False`` is reserved for full raw-SMILES recomputation:
    numerical validity alone does not attest unknown legacy generation history.
    """
    from _common import RAW_CSV, THRESHOLDS, sha256_file

    paths = cache_paths(directory)
    metadata_path = directory / "CACHE_METADATA.json"
    missing = [p.name for p in [*paths.values(), metadata_path] if not p.is_file()]
    if missing:
        raise CacheValidationError(f"Incomplete cache: {missing}; explicitly rebuild with force=True.")
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as error:
        raise CacheValidationError("Unreadable cache metadata") from error
    if metadata.get("source_sha256") != sha256_file(RAW_CSV):
        raise CacheValidationError("Cache raw source SHA-256 differs from the current raw CSV")
    if metadata.get("n_molecules") != len(study.df):
        raise CacheValidationError("Cache molecule count differs from the raw cohort")
    processed = pd.read_csv(paths["Processed"])
    expected = study.df.copy()
    expected["Tier"] = study.tier
    for threshold in THRESHOLDS:
        expected[f"Y{threshold}"] = study.y[threshold]
    try:
        pd.testing.assert_frame_equal(processed, expected, check_dtype=False, check_exact=True)
    except AssertionError as error:
        raise CacheValidationError("Cached processed compound identities, row order or labels differ from raw inputs") from error
    scaffolds = pd.read_csv(paths["Scaffolds"], keep_default_na=False)
    n = len(study.df)
    if (list(scaffolds.columns) != ["Index", "Scaffold"] or len(scaffolds) != n
            or not np.array_equal(scaffolds.Index.to_numpy(), np.arange(n))
            or scaffolds.Scaffold.eq("").any()):
        raise CacheValidationError("Scaffold cache rows/indices are incomplete or reordered")
    kernels = {}
    for rep in ("ECFP", "Avalon", "MACCS", "WL-HI"):
        try:
            kernel = np.load(paths[rep], mmap_mode="r", allow_pickle=False)
        except (OSError, ValueError) as error:
            raise CacheValidationError(f"Unreadable kernel {rep}") from error
        if kernel.shape != (n, n) or not np.isfinite(kernel).all():
            raise CacheValidationError(f"Kernel {rep} has invalid shape or nonfinite values")
        if (np.min(kernel) < -1e-6 or np.max(kernel) > 1 + 1e-6
                or not np.allclose(kernel, kernel.T, atol=1e-7, rtol=0)
                or not np.allclose(np.diag(kernel), 1, atol=1e-6, rtol=0)):
            raise CacheValidationError(f"Kernel {rep} fails symmetry, range or self-similarity checks")
        kernels[rep] = kernel
    descriptors = np.load(paths["Descriptors12"], mmap_mode="r", allow_pickle=False)
    if descriptors.shape != (n, 12) or not np.isfinite(descriptors).all():
        raise CacheValidationError("Descriptor cache has invalid shape or nonfinite values")

    schema = metadata.get("schema_version", 1)
    if schema == 2:
        evidence = metadata
        if evidence.get("provenance_kind") != "generated_from_raw_smiles":
            raise CacheValidationError("Unknown generation metadata provenance")
        recorded_shapes = evidence.get("array_shapes", {})
        actual_shapes = {**{rep: list(k.shape) for rep, k in kernels.items()},
                         "Descriptors12": list(descriptors.shape)}
        if recorded_shapes != actual_shapes:
            raise CacheValidationError("Cache array shapes differ from generation metadata")
        provenance = "generation_metadata_v2"
    elif schema == 1:
        report_path = directory / "CACHE_VALIDATION.json"
        if require_provenance:
            if not report_path.is_file():
                raise CacheValidationError(
                    "Legacy cache generation environment is unknown. Run "
                    "python verification/verify_revision.py --cache-only to independently "
                    "recompute and validate it without inventing historical provenance.")
            evidence = json.loads(report_path.read_text(encoding="utf-8"))
            if (evidence.get("basis") != "fresh_raw_smiles_recomputation"
                    or evidence.get("status") != "passed"
                    or evidence.get("legacy_metadata_sha256") != sha256_file(metadata_path)):
                raise CacheValidationError("Legacy cache validation record is incomplete or stale")
        else:
            evidence = None
        provenance = "legacy_cache_current_recomputation" if require_provenance else "legacy_unattested"
    else:
        raise CacheValidationError(f"Unsupported cache metadata schema {schema}")
    if evidence is not None:
        for key, value in current_context(study).items():
            if evidence.get(key) != value:
                raise CacheValidationError(f"Cache validation context changed: {key}")
        if evidence.get("artifact_sha256") != artifact_hashes(directory):
            raise CacheValidationError("Cache artifact SHA-256 mismatch")
    return {"study": study, "kernels": kernels, "descriptors": descriptors,
            "scaffolds": scaffolds.Scaffold.astype(str).to_numpy(),
            "validation": {"status": "passed", "provenance": provenance,
                           "legacy_generation_environment": "not_recorded" if schema == 1 else None}}


def audit_cache_from_raw(study, directory: Path) -> tuple[dict, dict]:
    """Fresh computation uses the documented implementation, never cached arrays.

    This is independence from cache contents, not a second independent algorithm.
    A validation report certifies present agreement, not the legacy build environment.
    """
    from _common import compute_representations, sha256_file

    stored = validate_cache(study, directory, require_provenance=False)
    fresh = compute_representations(study)
    comparisons = {}
    for rep in (*fresh["kernels"], "Descriptors12"):
        a = fresh["descriptors"] if rep == "Descriptors12" else fresh["kernels"][rep]
        b = stored["descriptors"] if rep == "Descriptors12" else stored["kernels"][rep]
        atol, rtol = (1e-10, 1e-12) if rep == "Descriptors12" else (1e-7, 1e-7)
        matches = bool(np.allclose(a, b, atol=atol, rtol=rtol))
        comparisons[rep] = {"shape": list(a.shape), "allclose": matches,
                            "exact_equal": bool(np.array_equal(a, b)),
                            "max_absolute_difference": float(np.max(np.abs(a - b))),
                            "atol": atol, "rtol": rtol}
        if not matches:
            raise CacheValidationError(f"Raw-SMILES recomputation differs for {rep}: {comparisons[rep]}")
    if not np.array_equal(fresh["scaffolds"], stored["scaffolds"]):
        raise CacheValidationError("Raw-SMILES scaffold recomputation differs")
    metadata_path = directory / "CACHE_METADATA.json"
    legacy = json.loads(metadata_path.read_text(encoding="utf-8")).get("schema_version", 1) == 1
    report = {
        "schema_version": 1, "status": "passed", "basis": "fresh_raw_smiles_recomputation",
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Current numerical agreement with the documented implementation; not independent algorithm validation",
        "legacy_generation_environment": "not_recorded" if legacy else None,
        "legacy_metadata_sha256": sha256_file(metadata_path),
        **current_context(study), "artifact_sha256": artifact_hashes(directory),
        "comparisons": comparisons, "scaffolds_exact_match": True,
    }
    (directory / "CACHE_VALIDATION.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    validate_cache(study, directory)
    return report, fresh
