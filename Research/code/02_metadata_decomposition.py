from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVC

from _common import (
    META_BINARY,
    META_CATEGORICAL,
    RESULTS,
    SEED,
    THRESHOLDS,
    build_or_load_cache,
    multinomial_weights,
    percentile_ci,
    simultaneous_centered_ci,
    tie_aware_weighted_auc,
)

OUT = RESULTS / "decomposition"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
df = study.df.copy()
k_ecfp = np.asarray(cache["kernels"]["ECFP"])
scaffolds = np.asarray(cache["scaffolds"], dtype=object)


def make_metadata_pipeline(categorical: list[str], binary: list[str]) -> Pipeline:
    transformers = []
    if categorical:
        transformers.append(
            (
                "cat",
                OneHotEncoder(handle_unknown="ignore", sparse_output=False),
                categorical,
            )
        )
    if binary:
        transformers.append(("bin", "passthrough", binary))
    prep = ColumnTransformer(transformers=transformers, remainder="drop")
    clf = LogisticRegression(
        C=1.0,
        solver="lbfgs",
        class_weight="balanced",
        max_iter=5000,
        random_state=SEED,
    )
    return Pipeline([("preprocess", prep), ("classifier", clf)])


def fit_structure_scores(tr: np.ndarray, te: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    model = SVC(C=1, kernel="precomputed", class_weight="balanced")
    model.fit(k_ecfp[np.ix_(tr, tr)], y[tr])
    train_score = model.decision_function(k_ecfp[np.ix_(tr, tr)])
    test_score = model.decision_function(k_ecfp[np.ix_(te, tr)])
    return train_score, test_score


def oof_structure_scores(tr: np.ndarray, y: np.ndarray, seed: int) -> np.ndarray:
    ytr = y[tr]
    local = np.arange(len(tr))
    oof = np.empty(len(tr), dtype=float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fit_local, val_local in cv.split(local, ytr):
        fit_idx = tr[fit_local]
        val_idx = tr[val_local]
        model = SVC(C=1, kernel="precomputed", class_weight="balanced")
        model.fit(k_ecfp[np.ix_(fit_idx, fit_idx)], y[fit_idx])
        oof[val_local] = model.decision_function(k_ecfp[np.ix_(val_idx, fit_idx)])
    return oof


model_specs = {
    "AgrochemicalFlags": ([], META_BINARY),
    "OriginRoute": (META_CATEGORICAL, []),
    "AllMetadata": (META_CATEGORICAL, META_BINARY),
}

performance_rows = []
pred_rows = []
pred_store: dict[tuple[str, int, str], np.ndarray] = {}

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    for thr in THRESHOLDS:
        y = study.y[thr]
        ytr, yte = y[tr], y[te]

        # Metadata-only baselines.
        for model_name, (categorical, binary) in model_specs.items():
            pipe = make_metadata_pipeline(list(categorical), list(binary))
            pipe.fit(df.iloc[tr], ytr)
            score = pipe.decision_function(df.iloc[te])
            pred_store[(split, thr, model_name)] = np.asarray(score, dtype=float)

        # Frozen ECFP structural probe.
        _, structure_test = fit_structure_scores(tr, te, y)
        pred_store[(split, thr, "ECFPStructure")] = structure_test

        # OOF stacking: structural score is out-of-fold for every training molecule;
        # the test score is generated only after refitting the structural model on all train data.
        oof_score = oof_structure_scores(tr, y, seed=SEED + split_i * 100 + thr)
        all_meta = make_metadata_pipeline(META_CATEGORICAL, META_BINARY)
        meta_train = all_meta.named_steps["preprocess"].fit_transform(df.iloc[tr])
        meta_test = all_meta.named_steps["preprocess"].transform(df.iloc[te])
        combined_train = np.column_stack([oof_score, meta_train])
        combined_test = np.column_stack([structure_test, meta_test])
        combined = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "classifier",
                    LogisticRegression(
                        C=1.0,
                        solver="lbfgs",
                        class_weight="balanced",
                        max_iter=5000,
                        random_state=SEED,
                    ),
                ),
            ]
        )
        combined.fit(combined_train, ytr)
        combined_score = combined.decision_function(combined_test)
        pred_store[(split, thr, "ECFP+AllMetadata")] = np.asarray(combined_score, dtype=float)

        for model_name in [
            "AgrochemicalFlags",
            "OriginRoute",
            "AllMetadata",
            "ECFPStructure",
            "ECFP+AllMetadata",
        ]:
            score = pred_store[(split, thr, model_name)]
            performance_rows.append(
                {
                    "Split": split,
                    "Threshold": thr,
                    "Model": model_name,
                    "AUROC": float(roc_auc_score(yte, score)),
                    "AUPRC": float(average_precision_score(yte, score)),
                    "n_test": int(len(te)),
                    "n_pos": int(yte.sum()),
                }
            )
            for local_i, global_i in enumerate(te):
                pred_rows.append(
                    {
                        "Split": split,
                        "Threshold": thr,
                        "Model": model_name,
                        "TestOrder": local_i,
                        "Index": int(global_i),
                        "CID": int(df.iloc[global_i]["CID"]),
                        "Scaffold": str(scaffolds[global_i]),
                        "Y": int(yte[local_i]),
                        "Score": float(score[local_i]),
                    }
                )

performance = pd.DataFrame(performance_rows)
performance.to_csv(OUT / "04_model_decomposition_performance.csv", index=False)
predictions = pd.DataFrame(pred_rows)
predictions.to_csv(OUT / "05_model_decomposition_predictions.csv", index=False)

# ---------- paired model-comparison bootstrap ----------
B = 5000
contrasts = [
    ("ECFPStructure", "AllMetadata"),
    ("ECFP+AllMetadata", "AllMetadata"),
    ("ECFP+AllMetadata", "ECFPStructure"),
    ("AllMetadata", "AgrochemicalFlags"),
]
main_family = contrasts[:3]
boot_rows = []
npz_payload = {}

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    rng = np.random.default_rng(SEED + 1700 + split_i)
    wm = multinomial_weights(len(te), B, rng)
    unique_scaffolds, inverse = np.unique(scaffolds[te], return_inverse=True)
    wg = multinomial_weights(len(unique_scaffolds), B, rng)
    ws = wg[:, inverse]

    for thr in THRESHOLDS:
        yy = study.y[thr][te]
        observed_auc = {
            model: float(roc_auc_score(yy, pred_store[(split, thr, model)]))
            for model in sorted({m for pair in contrasts for m in pair})
        }
        for mode, weights in [("Molecule", wm), ("ScaffoldCluster", ws)]:
            auc_boot = {
                model: tie_aware_weighted_auc(yy, pred_store[(split, thr, model)], weights)
                for model in observed_auc
            }
            family_obs = np.array([observed_auc[a] - observed_auc[b] for a, b in main_family])
            family_boot = np.column_stack([auc_boot[a] - auc_boot[b] for a, b in main_family])
            valid_family = np.all(np.isfinite(family_boot), axis=1)
            sim_low, sim_high = simultaneous_centered_ci(
                family_obs, family_boot[valid_family], alpha=0.05
            )

            for a, b in contrasts:
                d = auc_boot[a] - auc_boot[b]
                d = d[np.isfinite(d)]
                low, high = percentile_ci(d)
                obs = observed_auc[a] - observed_auc[b]
                if (a, b) in main_family:
                    j = main_family.index((a, b))
                    family_low, family_high = float(sim_low[j]), float(sim_high[j])
                else:
                    family_low, family_high = np.nan, np.nan
                comparison = f"{a}-{b}"
                boot_rows.append(
                    {
                        "Split": split,
                        "Threshold": thr,
                        "Bootstrap": mode,
                        "Comparison": comparison,
                        "PrimaryFamily": (a, b) in main_family,
                        "ObservedDeltaAUROC": float(obs),
                        "BootstrapMean": float(d.mean()),
                        "PercentileCI_low": low,
                        "PercentileCI_high": high,
                        "SimultaneousCI_low": family_low,
                        "SimultaneousCI_high": family_high,
                        "ProbDelta_gt0": float(np.mean(d > 0)),
                        "B_valid": int(len(d)),
                    }
                )
                npz_payload[f"{split}_{thr}_{mode}_{comparison}"] = d.astype(np.float32)

boot = pd.DataFrame(boot_rows)
boot.to_csv(OUT / "06_model_decomposition_bootstrap.csv", index=False)
np.savez_compressed(OUT / "06_model_decomposition_bootstrap_distributions.npz", **npz_payload)

# ---------- performance within insecticide strata (same globally trained models) ----------
strata_rows = []
for split, (tr, te) in study.splits.items():
    for thr in THRESHOLDS:
        yy = study.y[thr][te]
        for insecticide_status in [0, 1]:
            mask = df.iloc[te]["insecticide"].to_numpy() == insecticide_status
            for model in ["AllMetadata", "ECFPStructure", "ECFP+AllMetadata"]:
                ysub = yy[mask]
                score = pred_store[(split, thr, model)][mask]
                auc = np.nan
                if len(ysub) >= 10 and len(np.unique(ysub)) == 2:
                    auc = float(roc_auc_score(ysub, score))
                strata_rows.append(
                    {
                        "Split": split,
                        "Threshold": thr,
                        "Insecticide": insecticide_status,
                        "Model": model,
                        "n": int(mask.sum()),
                        "n_pos": int(ysub.sum()),
                        "AUROC": auc,
                    }
                )
pd.DataFrame(strata_rows).to_csv(OUT / "07_insecticide_stratified_performance.csv", index=False)

# ---------- leave-one-source-out sensitivity (supplementary) ----------
source_rows = []
source_prediction_rows = []
source_score_store = {}
for heldout in sorted(df["source"].unique()):
    te = np.flatnonzero(df["source"].to_numpy() == heldout)
    tr = np.flatnonzero(df["source"].to_numpy() != heldout)
    for thr in THRESHOLDS:
        y = study.y[thr]
        if len(np.unique(y[te])) < 2 or len(np.unique(y[tr])) < 2:
            continue
        structure = SVC(C=1, kernel="precomputed", class_weight="balanced")
        structure.fit(k_ecfp[np.ix_(tr, tr)], y[tr])
        structure_score = structure.decision_function(k_ecfp[np.ix_(te, tr)])

        route_flags = make_metadata_pipeline(["toxicity_type"], META_BINARY)
        route_flags.fit(df.iloc[tr], y[tr])
        metadata_score = route_flags.decision_function(df.iloc[te])
        for model, score in [
            ("ECFPStructure", structure_score),
            ("Route+Agrochemical", metadata_score),
        ]:
            source_score_store[(heldout, thr, model)] = np.asarray(score, dtype=float)
            source_rows.append(
                {
                    "HeldOutSource": heldout,
                    "Threshold": thr,
                    "Model": model,
                    "n_test": int(len(te)),
                    "n_pos": int(y[te].sum()),
                    "AUROC": float(roc_auc_score(y[te], score)),
                    "AUPRC": float(average_precision_score(y[te], score)),
                }
            )
            for local_i, global_i in enumerate(te):
                source_prediction_rows.append({
                    "HeldOutSource": heldout, "Threshold": thr, "Model": model,
                    "TestOrder": local_i, "Index": int(global_i),
                    "Scaffold": str(scaffolds[global_i]), "Y": int(y[global_i]),
                    "Score": float(score[local_i]),
                })
pd.DataFrame(source_rows).to_csv(OUT / "08_leave_one_source_out.csv", index=False)
pd.DataFrame(source_prediction_rows).to_csv(OUT / "09_leave_one_source_out_predictions.csv", index=False)

source_boot_rows = []
for source_i, heldout in enumerate(sorted(df["source"].unique())):
    te = np.flatnonzero(df["source"].to_numpy() == heldout)
    rng = np.random.default_rng(SEED + 2600 + source_i)
    wm = multinomial_weights(len(te), B, rng)
    unique_sc, inverse = np.unique(scaffolds[te], return_inverse=True)
    ws = multinomial_weights(len(unique_sc), B, rng)[:, inverse]
    for mode, weights in [("Molecule", wm), ("ScaffoldCluster", ws)]:
        auc_obs = {}
        auc_boot = {}
        for thr in THRESHOLDS:
            yy = study.y[thr][te]
            score = source_score_store[(heldout, thr, "ECFPStructure")]
            auc_obs[thr] = float(roc_auc_score(yy, score))
            auc_boot[thr] = tie_aware_weighted_auc(yy, score, weights)
        family = [(11, 100), (1, 100)]
        observed = np.array([auc_obs[a] - auc_obs[b] for a, b in family])
        boot_matrix = np.column_stack([auc_boot[a] - auc_boot[b] for a, b in family])
        valid = np.all(np.isfinite(boot_matrix), axis=1)
        sim_low, sim_high = simultaneous_centered_ci(observed, boot_matrix[valid])
        for j, (a, b) in enumerate(family):
            d = boot_matrix[:, j]
            d = d[np.isfinite(d)]
            low, high = percentile_ci(d)
            source_boot_rows.append({
                "HeldOutSource": heldout, "Bootstrap": mode,
                "Comparison": f"{a}-{b}",
                "ObservedDeltaAUROC": float(observed[j]),
                "PercentileCI_low": low, "PercentileCI_high": high,
                "SimultaneousCI_low": float(sim_low[j]),
                "SimultaneousCI_high": float(sim_high[j]),
                "ProbDelta_gt0": float(np.mean(d > 0)), "B_valid": int(len(d)),
            })
pd.DataFrame(source_boot_rows).to_csv(OUT / "10_leave_one_source_out_bootstrap.csv", index=False)

print("Model decomposition performance:")
print(
    performance.pivot_table(index=["Split", "Model"], columns="Threshold", values="AUROC")
    .round(3)
    .to_string()
)
print("\nPrimary decomposition contrasts (molecule bootstrap):")
print(
    boot[(boot["Bootstrap"] == "Molecule") & boot["PrimaryFamily"]][
        [
            "Split",
            "Threshold",
            "Comparison",
            "ObservedDeltaAUROC",
            "PercentileCI_low",
            "PercentileCI_high",
            "SimultaneousCI_low",
            "SimultaneousCI_high",
            "ProbDelta_gt0",
        ]
    ].round(4).to_string(index=False)
)
