"""Export submission tables and figure data from the numerical analysis outputs.

This step fits no models and recomputes no inferential analyses. It selects,
renames and aggregates saved outputs from steps 00--08. Existing final/ and
source_data/ files are outputs only; neither is a permitted input directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_SEED = 20260829
SPLITS = ("Random", "MaxMin", "Time")
MODEL_NAMES = {
    "AgrochemicalFlags": "Agrochemical flags",
    "OriginRoute": "Source + route",
    "AllMetadata": "All metadata",
    "ECFPStructure": "ECFP structure",
    "ECFP+AllMetadata": "Structure + metadata",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(input_results: Path, baseline_dir: Path, output_results: Path,
           source_data: Path) -> dict:
    final = output_results / "final"
    final.mkdir(parents=True, exist_ok=True)
    source_data.mkdir(parents=True, exist_ok=True)
    inputs, outputs = {}, {}

    def canonical_dependencies(dependencies: list[str]) -> list[str]:
        return [d if d.startswith(("data/", "results/")) else "results/" + d
                for d in dependencies]

    def read(relative: str) -> pd.DataFrame:
        path = input_results / relative
        if relative.split("/")[0] in {"final", "source_data"}:
            raise ValueError("A generated export cannot be used as input")
        inputs[f"results/{relative}"] = sha256(path)
        return pd.read_csv(path)

    def read_baseline(relative: str) -> pd.DataFrame:
        path = baseline_dir / relative
        inputs[f"data/{relative}"] = sha256(path)
        return pd.read_csv(path)

    def csv(frame: pd.DataFrame, name: str, dependencies: list[str], operation: str,
            figure: bool = False) -> None:
        folder = source_data if figure else final
        path = folder / name
        frame.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
        key = ("source_data/" if figure else "results/final/") + name
        outputs[key] = {"rows": len(frame), "columns": list(frame.columns),
                        "inputs": canonical_dependencies(dependencies), "operation": operation,
                        "sha256": sha256(path)}

    def json_file(value: dict, name: str, dependencies: list[str], operation: str,
                  figure: bool = False) -> None:
        path = (source_data if figure else final) / name
        path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        key = ("source_data/" if figure else "results/final/") + name
        outputs[key] = {"inputs": canonical_dependencies(dependencies), "operation": operation, "sha256": sha256(path)}

    paths = {
        "perf": "primary/01_representation_performance.csv",
        "boot": "primary/03_ecfp_paired_bootstrap.csv",
        "rboot": "primary/05_representation_paired_bootstrap.csv",
        "comp": "decomposition/01_tier_composition_long.csv",
        "assoc": "decomposition/03_tier_composition_association_tests.csv",
        "model": "decomposition/04_model_decomposition_performance.csv",
        "curve": "learning_curves/01_matched_learning_curves_long.csv",
        "matched": "learning_curves/02_matched_learning_curves_summary.csv",
        "contrast": "learning_curves/03_matched_learning_curve_contrasts.csv",
        "perm": "geometry/02_conditional_neighborhood_permutation.csv",
        "scaf": "geometry/08_scaffold_sensitivity_bootstrap.csv",
        "nov": "applicability/03_boundary_specific_brier_slopes.csv",
        "joint": "applicability/04_joint_interaction_tests.csv",
    }
    data = {key: read(path) for key, path in paths.items()}
    perf, boot, rboot = data["perf"], data["boot"], data["rboot"]
    comp, model, perm, scaf, nov = (data[x] for x in ("comp", "model", "perm", "scaf", "nov"))
    matched = data["matched"].loc[lambda d: d.n_per_class.eq(125)].copy()
    primary_boot = boot[boot.PrimaryComparison.eq(True)].copy()

    # Validate the raw labels and official split membership directly.
    raw = read_baseline("raw/apistox.csv")
    tier = np.select([raw.ppdb_level.eq(0), raw.ppdb_level.eq(1) & raw.label.eq(0),
                      raw.ppdb_level.eq(1) & raw.label.eq(1), raw.ppdb_level.eq(2)],
                     [0, 1, 2, 3], default=-1)
    if (tier < 0).any() or raw.SMILES.duplicated().any():
        raise ValueError("Raw labels or molecule identities violate the frozen protocol")
    validation = {
        "seed": PROTOCOL_SEED, "n_rows": len(raw), "n_unique_smiles": int(raw.SMILES.nunique()),
        "tier_counts": {str(k): int(v) for k, v in pd.Series(tier).value_counts().sort_index().items()},
        "source_counts": {str(k): int(v) for k, v in raw.source.value_counts().items()},
        "nested_endpoint_invariant": bool(np.all((tier >= 3) <= (tier >= 2)) and np.all((tier >= 2) <= (tier >= 1))),
        "source_csv_sha256": inputs["data/raw/apistox.csv"],
    }
    for split in SPLITS:
        a = read_baseline(f"official_splits/{split.lower()}_train.csv")
        b = read_baseline(f"official_splits/{split.lower()}_test.csv")
        if not set(a.SMILES).union(b.SMILES).issubset(set(raw.SMILES)):
            raise ValueError(f"Unknown split molecule: {split}")
        validation[f"split_{split}"] = {"train_n": len(a), "test_n": len(b),
                                        "overlap_n": len(set(a.SMILES).intersection(b.SMILES))}
    raw_deps = [x for x in inputs if x.startswith("data/")]
    json_file(validation, "00_data_validation.json", raw_deps, "Validate raw label mapping and official split membership")

    rep_summary = perf.groupby(["Split", "Threshold"], as_index=False).agg(
        MeanAUROC=("AUROC", "mean"), MinAUROC=("AUROC", "min"), MaxAUROC=("AUROC", "max"))
    csv(rep_summary, "01_representation_summary.csv", [paths["perf"]], "Group AUROC by split and threshold")
    csv(primary_boot, "02_primary_bootstrap_table.csv", [paths["boot"]], "Select PrimaryComparison=True; preserve inferential columns")
    csv(comp[comp.Variable.eq("agrochemical_flag")], "03_key_tier_composition.csv", [paths["comp"]], "Select four separate agrochemical flags")
    csv(model, "04_model_decomposition_table.csv", [paths["model"]], "Preserve original decomposition estimates; preprocessing sensitivity remains separate")
    csv(matched, "05_matched_n125_table.csv", [paths["matched"]], "Select n_per_class=125; intervals describe repeated fits")
    perm_summary = perm.groupby("PermutationScheme", sort=False, as_index=False).agg(
        Significant_q05=("Significant_q05", "sum"), Total=("Significant_q05", "size"),
        Max_q=("BH_q_within_scheme", "max"), MedianObservedEffect=("ObservedEffect", "median"))
    csv(perm_summary, "06_conditional_neighborhood_summary.csv", [paths["perm"]], "Summarize saved within-scheme BH decisions")
    csv(scaf, "07_scaffold_sensitivity_table.csv", [paths["scaf"]], "Preserve all saved scaffold sensitivity contrasts")
    csv(nov, "08_novelty_slope_table.csv", [paths["nov"]], "Preserve endpoint-specific Brier slopes and unadjusted marginal intervals")

    pv = perf.pivot(index=["Split", "Representation"], columns="Threshold", values="AUROC")
    ordering = (pv[100] < pv[11]) & (pv[100] < pv[1])
    mv = model.pivot(index=["Split", "Threshold"], columns="Model", values="AUROC")
    counts = {"combined_better_than_metadata": int((mv["ECFP+AllMetadata"] > mv.AllMetadata).sum()),
              "combined_better_than_structure": int((mv["ECFP+AllMetadata"] > mv.ECFPStructure).sum())}
    c125 = data["contrast"].query("n_per_class == 125 and Comparison == '11-100'").set_index("Split")
    cv = float(data["assoc"].set_index("Variable").loc["insecticide", "CramerV"])
    novelty100 = nov[nov.Boundary.eq(100)].set_index("Split")
    joint = data["joint"].query("Outcome == 'BrierLoss'").set_index("Split")
    key_numbers = {
        "seed": PROTOCOL_SEED, "n_molecules": len(raw), "representation_ordering_pass": int(ordering.sum()),
        "insecticide_cramers_v": cv,
        "matched_n125_11_minus_100": {s: float(c125.loc[s, "MeanDeltaAUROC"]) for s in sorted(c125.index)},
        "neighborhood_significant": {r.PermutationScheme: int(r.Significant_q05) for r in perm_summary.itertuples()},
        **counts, "scaffold_positive_ci_count": int(scaf.CI_low.gt(0).sum()),
        "novelty_100_slopes": {s: {"slope_per_0.1": float(novelty100.loc[s, "Slope_per_0.1_novelty"]),
                                   "ci_low": float(novelty100.loc[s, "CI_low"]), "ci_high": float(novelty100.loc[s, "CI_high"])} for s in SPLITS},
        "brier_joint_interaction_p": {s: float(joint.loc[s, "P"]) for s in SPLITS},
    }
    json_file(key_numbers, "01_key_numbers.json", list(paths.values()), "Calculate counts and select reported estimates; no frozen-number input")

    # Figure values retain full precision. Display rounding belongs to the renderer.
    csv(perf[["Split", "Representation", "Threshold", "AUROC"]].rename(columns={"Split": "split", "Representation": "representation", "Threshold": "endpoint", "AUROC": "auroc"}),
        "figure_primary_data.csv", [paths["perf"]], "Rename four plotting columns; retain full precision", True)
    fb = primary_boot[["Split", "Bootstrap", "Comparison", "ObservedDeltaAUROC", "SimultaneousCI_low", "SimultaneousCI_high"]].copy()
    fb.columns = ["split", "level", "contrast", "estimate", "low", "high"]
    fb["level"] = fb.level.replace({"ScaffoldCluster": "Scaffold"})
    fb["contrast"] = fb.contrast.str.replace("-", "−", regex=False)
    csv(fb, "figure_bootstrap_data.csv", [paths["boot"]], "Primary ECFP contrasts with simultaneous intervals", True)
    comp_columns = [("source", "PPDB", "PPDB"), ("source", "ECOTOX", "ECOTOX"),
                    ("toxicity_type", "Contact", "Contact"), ("agrochemical_flag", "herbicide", "Herbicide"),
                    ("agrochemical_flag", "fungicide", "Fungicide"), ("agrochemical_flag", "insecticide", "Insecticide"),
                    ("source", "BPDB", "BPDB"), ("toxicity_type", "Oral", "Oral"),
                    ("toxicity_type", "Other", "OtherRoute"), ("agrochemical_flag", "other_agrochemical", "OtherAgrochemical")]
    comp_plot = pd.DataFrame({"Tier": [f"Tier{x}" for x in sorted(comp.Tier.unique())]})
    for variable, level, label in comp_columns:
        series = comp[comp.Variable.eq(variable) & comp.Level.eq(level)].set_index("Tier").Proportion
        comp_plot[label] = [float(series.get(t, 0.0)) for t in sorted(comp.Tier.unique())]
    csv(comp_plot, "figure_composition_data.csv", [paths["comp"]], "Pivot all source/route/use groups shown in current Figure3; flags are not exclusive", True)
    cm = model[["Split", "Model", "Threshold", "AUROC"]].copy()
    cm.columns = ["split", "model", "endpoint", "auroc"]
    cm["model"] = cm.model.map(MODEL_NAMES)
    if cm.model.isna().any():
        raise ValueError("Unknown decomposition model name")
    csv(cm, "figure_model_decomposition.csv", [paths["model"]], "All five designated model baselines, including source+route", True)
    curve = data["curve"][["Split", "Threshold", "n_per_class", "AUROC"]].copy()
    curve.columns = ["split", "endpoint", "n", "auroc"]
    csv(curve, "figure_learning_curve_raw.csv", [paths["curve"]], "All repeat-level AUROC rows in original order; n is compounds per class", True)
    mat_plot = matched[["Split", "Threshold", "MeanAUROC"]].copy()
    mat_plot.columns = ["split", "endpoint", "auroc"]
    csv(mat_plot, "figure_matched_n125.csv", [paths["matched"]], "Mean AUROC at n_per_class=125", True)
    nov_plot = nov[["Split", "Boundary", "Slope_per_0.1_novelty", "CI_low", "CI_high"]].copy()
    nov_plot.columns = ["split", "endpoint", "estimate", "low", "high"]
    csv(nov_plot, "figure_novelty_raw.csv", [paths["nov"]], "Brier-loss slopes per 0.1 novelty with marginal intervals", True)
    perm_plot = perm[["Representation", "k", "PermutationScheme", "ObservedEffect", "NullMean", "NullQ025", "NullQ975", "BH_q_within_scheme", "Significant_q05"]].copy()
    perm_plot.columns = ["representation", "k", "scheme", "observed_effect", "null_mean", "null_q025", "null_q975", "bh_q", "significant_q05"]
    csv(perm_plot, "figure_permutation_raw.csv", [paths["perm"]], "Permutation-null quantiles and within-scheme Benjamini-Hochberg adjusted p values", True)
    rb = rboot[rboot.PrimaryComparison.eq(True) & rboot.Bootstrap.eq("ScaffoldCluster")].copy()
    csv(rb, "figure_representation_bootstrap.csv", [paths["rboot"]], "All representation primary scaffold contrasts", True)
    sf = scaf[scaf.Analysis.isin(["ScaffoldEqualWeighted", "PairedScaffoldsContainingBothTiers"]) & scaf.Comparison.isin(["T2-T1", "T2-T3"])].copy()
    sf["MinScaffoldN"] = sf.MinScaffoldN.map(lambda n: f"n≥{n}")
    sf["Analysis"] = sf.Analysis.map({"ScaffoldEqualWeighted": "Scaffold-equal", "PairedScaffoldsContainingBothTiers": "Paired shared-scaffold"})
    sf["Comparison"] = sf.Comparison.str.replace("-", "−", regex=False)
    sf = sf[["MinScaffoldN", "Analysis", "Comparison", "ObservedDifference", "CI_low", "CI_high"]]
    sf.columns = ["minimum", "analysis", "contrast", "estimate", "low", "high"]
    csv(sf, "figure_scaffold_corrected.csv", [paths["scaf"]], "Supplementary Fig. S3 scaffold-equal/shared-scaffold T2-T1 and T2-T3 contrasts", True)

    dictionary = {
        "precision": "Exported values retain analysis precision; renderers may format to 3 decimal places.",
        "endpoint": "Cumulative label threshold in ug/bee; these are frozen benchmark labels, not re-extracted continuous measurements.",
        "figure_permutation_raw.csv": {"null_q025": "2.5% quantile of the permutation null distribution", "bh_q": "Within-permutation-scheme Benjamini-Hochberg adjusted p value"},
        "figure_composition_data.csv": {"Tier": "Benchmark toxicity-tier row identifier", "OtherRoute": "toxicity_type=Other", "OtherAgrochemical": "Separate binary use flag"},
        "figure_learning_curve_raw.csv": {"n": "Training compounds per class", "rows": "Repeat-level fits; original repeat order retained from analysis table"},
        "figure_bootstrap_data.csv": {"low/high": "Simultaneous 95% interval for the two primary endpoint contrasts within split"},
        "figure_scaffold_corrected.csv": {"low/high": "Saved percentile 95% scaffold-bootstrap interval"},
        "figure_novelty_raw.csv": {"estimate": "Brier-loss slope per +0.1 novelty", "low/high": "Unadjusted marginal 95% interval"},
    }
    json_file(dictionary, "SOURCE_DATA_DICTIONARY.json", list(paths.values()), "Column definitions and precision policy", True)
    provenance = {"generator": "code/09_export_submission_tables.py", "model_fits": 0,
                  "reads_previous_exports": False, "input_sha256": inputs, "outputs": outputs}
    (final / "EXPORT_PROVENANCE.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-results", type=Path, default=ROOT / "results")
    parser.add_argument("--baseline-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--output-results", type=Path, default=ROOT / "results")
    parser.add_argument("--source-data", type=Path, default=ROOT / "source_data")
    args = parser.parse_args()
    result = export(args.input_results.resolve(), args.baseline_dir.resolve(),
                    args.output_results.resolve(), args.source_data.resolve())
    print(f"Exported {len(result['outputs'])} submission artifacts from {len(result['input_sha256'])} raw/analysis inputs; no model fits.")


if __name__ == "__main__":
    main()
