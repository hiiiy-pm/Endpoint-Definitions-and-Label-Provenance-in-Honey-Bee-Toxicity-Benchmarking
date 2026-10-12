"""Describe ECFP similarity to selected regulatory cohorts; no fitting.

For each benchmark compound, report its maximum Tanimoto similarity to the
selected source/route/all-three-determinate cohort. Exclude self matches for
selected compounds. This is a reference-set-specific description, not a
representativeness test or an independent experimental validation.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/record_audit"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = {"benchmark": ROOT / "data/raw/apistox.csv",
             "kernel": ROOT / "results/cache/kernel_ecfp.npy",
             "cache_validation": ROOT / "results/cache/CACHE_VALIDATION.json",
             "route_matched": OUT / "07_route_matched_external_records.csv",
             "cohort_flow": OUT / "13_concordance_cohort_flow.csv"}
    validation = json.loads(paths["cache_validation"].read_text())
    assert validation["status"] == "passed"
    assert sha(paths["kernel"]) == validation["artifact_sha256"]["kernel_ecfp.npy"]
    assert sha(paths["benchmark"]) == validation["input_sha256"]["data/raw/apistox.csv"]
    bench = pd.read_csv(paths["benchmark"])
    kernel = np.load(paths["kernel"], allow_pickle=False)
    assert kernel.shape == (len(bench), len(bench)) and np.isfinite(kernel).all()
    cohort = pd.read_csv(paths["route_matched"])
    flow = pd.read_csv(paths["cohort_flow"])
    # Independent primitive check: regenerate fingerprints from raw SMILES and
    # compare every maximum against RDKit bulk Tanimoto, not the cached maxima.
    generator = GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
    fps = [generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in bench.SMILES]
    frames, summaries, errors = [], [], {}
    for screen in ("OpenFoodTox", "OpenFoodTox_48h_screen", "EPA"):
        src = cohort[cohort.source.eq("EPA" if screen == "EPA" else "OpenFoodTox")]
        if screen == "OpenFoodTox_48h_screen":
            src = src[src.source_48h_screen]
        selected = []
        for idx, group in src.groupby("Index"):
            if all(set(group[f"external_y_{t}"]) in ({0}, {1}) for t in (1, 11, 100)):
                selected.append(int(idx))
        selected = np.asarray(sorted(selected), dtype=int)
        expected = flow[(flow.source_screen == screen) & (flow.stage == "all_three_thresholds_determinate")].benchmark_compounds.item()
        assert len(selected) == expected and len(selected) > 1
        values = np.asarray(kernel[:, selected], dtype=float).copy()
        for column, idx in enumerate(selected):
            values[idx, column] = -np.inf
        maximum = values.max(axis=1)
        nearest = selected[values.argmax(axis=1)]
        independent = []
        for idx, fp in enumerate(fps):
            exact = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, [fps[j] for j in selected]))
            exact[selected == idx] = -np.inf
            independent.append(float(exact.max()))
        error = float(np.max(np.abs(maximum - np.asarray(independent))))
        assert error <= 1e-7, (screen, error)
        errors[screen] = error
        included = np.isin(np.arange(len(bench)), selected)
        frame = pd.DataFrame({"source_screen": screen, "Index": np.arange(len(bench)),
                              "CAS": bench.CAS, "name": bench["name"],
                              "membership": np.where(included, "included", "excluded"),
                              "selected_cohort_size": len(selected),
                              "reference_compounds": len(selected) - included.astype(int),
                              "max_ecfp_tanimoto_to_selected": maximum,
                              "nearest_selected_Index": nearest,
                              "nearest_selected_CAS": bench.CAS.to_numpy()[nearest]})
        assert not ((frame.membership == "included") & (frame.Index == frame.nearest_selected_Index)).any()
        frames.append(frame)
        for membership, group in frame.groupby("membership", sort=True):
            q = group.max_ecfp_tanimoto_to_selected.to_numpy()
            summaries.append({"source_screen": screen, "membership": membership,
                              "n_compounds": len(group), "selected_cohort_size": len(selected),
                              "reference_compounds_per_query": int(group.reference_compounds.iloc[0]),
                              "median": float(np.median(q)), "q1": float(np.quantile(q, .25)),
                              "q3": float(np.quantile(q, .75)), "minimum": float(q.min()), "maximum": float(q.max()),
                              "self_matches_excluded": True,
                              "metric": "max ECFP4 Tanimoto to selected cohort, excluding self",
                              "interpretation": "descriptive similarity to this selected set; not proof of representativeness"})
    pd.DataFrame(summaries).to_csv(OUT / "16_selection_structure_summary.csv", index=False, lineterminator="\n")
    pd.concat(frames, ignore_index=True).to_csv(OUT / "17_selection_structure_per_compound.csv", index=False, lineterminator="\n")
    metadata = {"status": "passed", "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in paths.values()},
                "parameters": {"radius": 2, "bits": 1024, "chirality": False, "quantiles": "NumPy default linear interpolation"},
                "independent_raw_smiles_maximum_checks": len(bench) * 3,
                "maximum_absolute_difference_from_direct_rdkit": errors,
                "nearest_tie_rule": "smallest benchmark Index among tied cached similarities",
                "scope": "whole-benchmark selected cohorts; not split-specific frozen-score reassessment",
                "limitations": "Maxima depend on selected cohort composition and size. Included queries use n-1 references and excluded queries n; this comparison does not establish representativeness, distributional equivalence or experiment independence."}
    (OUT / "SELECTION_STRUCTURE_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(pd.DataFrame(summaries)[["source_screen", "membership", "n_compounds", "median", "q1", "q3", "minimum", "maximum"]].to_string(index=False))


if __name__ == "__main__":
    main()
