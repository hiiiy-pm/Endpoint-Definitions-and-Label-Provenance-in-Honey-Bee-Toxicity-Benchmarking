"""Diagnostics added in response to pre-submission review.

1. Cross-endpoint evaluation: scores trained for one endpoint, evaluated on another.
2. Brier skill against a constant training-prevalence forecast.
3. Metadata versus fusion under identical preprocessing, with and without exposure route.
4. Simultaneous intervals over the wider family of ten contrasts per split.
5. Learner sensitivity: L2 logistic regression and random forest on ECFP4 bits.

Nothing here changes the frozen primary analysis; every model is fitted on training
compounds only.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
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
    tie_aware_weighted_auc,
)

OUT = RESULTS / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)
B = 2000

cache = build_or_load_cache(force=False)
study = cache["study"]
df = study.df
k_ecfp = np.asarray(cache["kernels"]["ECFP"])
scaffolds = np.asarray(cache["scaffolds"], dtype=object)


def scaffold_weights(te: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    uniq, inv = np.unique(scaffolds[te], return_inverse=True)
    return multinomial_weights(len(uniq), B, rng)[:, inv]


# ---------------------------------------------------------------- 1. cross-endpoint
pred = pd.read_csv(RESULTS / "primary" / "02_test_predictions.csv")
rows = []
for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    rng = np.random.default_rng(SEED + 7001 + split_i)
    w = scaffold_weights(te, rng)
    for rep in ("ECFP", "WL-HI"):
        for train_thr in THRESHOLDS:
            s = (
                pred[(pred.Split == split) & (pred.Representation == rep)
                     & (pred.Threshold == train_thr)]
                .sort_values("TestOrder").Score.to_numpy()
            )
            for eval_thr in THRESHOLDS:
                y = study.y[eval_thr][te]
                boot = tie_aware_weighted_auc(y, s, w)
                lo, hi = percentile_ci(boot)
                rows.append({
                    "Split": split, "Representation": rep,
                    "TrainEndpoint": train_thr, "EvalEndpoint": eval_thr,
                    "AUROC": float(roc_auc_score(y, s)), "CI_low": lo, "CI_high": hi,
                })
cross = pd.DataFrame(rows)
cross.to_csv(OUT / "08a_cross_endpoint_evaluation.csv", index=False)

# ---------------------------------------------------------------- 2. Brier skill
cal = pd.read_csv(RESULTS / "applicability" / "01_calibrated_predictions_and_novelty.csv")
brier_rows = []
for (split, thr), g in cal.groupby(["Split", "Boundary"]):
    tr, _ = study.splits[split]
    p0 = float(study.y[int(thr)][tr].mean())
    model = float(g.BrierLoss.mean())
    base = float(((g.Y - p0) ** 2).mean())
    brier_rows.append({
        "Split": split, "Endpoint": int(thr), "ModelBrier": model,
        "TrainingPrevalence": p0, "ConstantBrier": base, "BrierSkill": 1 - model / base,
    })
brier = pd.DataFrame(brier_rows)
brier.to_csv(OUT / "08b_brier_skill.csv", index=False)

# ---------------------------------------------------------------- 3. fair metadata
def meta_prep(categorical):
    parts = []
    if categorical:
        parts.append(("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False),
                      categorical))
    parts.append(("bin", "passthrough", META_BINARY))
    return ColumnTransformer(parts, remainder="drop")


def logit():
    return LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000,
                              random_state=SEED)


def oof_ecfp(tr, y, seed):
    oof = np.empty(len(tr))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for f, v in cv.split(np.arange(len(tr)), y[tr]):
        m = SVC(C=1, kernel="precomputed", class_weight="balanced")
        m.fit(k_ecfp[np.ix_(tr[f], tr[f])], y[tr[f]])
        oof[v] = m.decision_function(k_ecfp[np.ix_(tr[v], tr[f])])
    return oof


meta_rows = []
for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    for thr in THRESHOLDS:
        y = study.y[thr]
        svm = SVC(C=1, kernel="precomputed", class_weight="balanced")
        svm.fit(k_ecfp[np.ix_(tr, tr)], y[tr])
        s_te = svm.decision_function(k_ecfp[np.ix_(te, tr)])
        s_oof = oof_ecfp(tr, y, SEED + split_i * 100 + thr)
        for label, cats in (("AllMetadata", list(META_CATEGORICAL)),
                            ("MetadataNoRoute", ["source"])):
            prep = meta_prep(cats)
            m_tr = prep.fit_transform(df.iloc[tr])
            m_te = prep.transform(df.iloc[te])
            meta_only = Pipeline([("scale", StandardScaler()), ("clf", logit())])
            meta_only.fit(m_tr, y[tr])
            fused = Pipeline([("scale", StandardScaler()), ("clf", logit())])
            fused.fit(np.column_stack([s_oof, m_tr]), y[tr])
            meta_rows.append({
                "Split": split, "Endpoint": thr, "MetadataSet": label,
                "MetadataAUROC": float(roc_auc_score(y[te], meta_only.decision_function(m_te))),
                "StructureAUROC": float(roc_auc_score(y[te], s_te)),
                "FusedAUROC": float(roc_auc_score(
                    y[te], fused.decision_function(np.column_stack([s_te, m_te])))),
            })
meta = pd.DataFrame(meta_rows)
meta.to_csv(OUT / "08c_matched_preprocessing_metadata.csv", index=False)

# ---------------------------------------------------------------- 4. wider family
dist = np.load(RESULTS / "primary" / "05_representation_paired_bootstrap_distributions.npz")
obs = pd.read_csv(RESULTS / "primary" / "05_representation_paired_bootstrap.csv")
obs = obs[obs.PrimaryComparison]
reps = ["ECFP", "Avalon", "MACCS", "WL-HI", "Descriptors12"]
fam_rows = []
for mode in ("Molecule", "ScaffoldCluster"):
    for split in study.splits:
        keys, o = [], []
        for rep in reps:
            for a, b in ((11, 100), (1, 100)):
                keys.append((rep, a, b))
                o.append(float(obs[(obs.Representation == rep) & (obs.Split == split)
                                   & (obs.Bootstrap == mode)
                                   & (obs.Comparison == f"{a}-{b}")].ObservedDeltaAUROC.iloc[0]))
        o = np.array(o)
        boots = np.column_stack([dist[f"{r}_{split}_{mode}_{a}_{b}"] for r, a, b in keys])
        q = np.quantile(np.max(np.abs(boots - o[None, :]), axis=1), 0.95)
        for (rep, a, b), d in zip(keys, o):
            fam_rows.append({"Bootstrap": mode, "Split": split, "Representation": rep,
                             "Comparison": f"{a}-{b}", "Delta": d,
                             "FamilySize": len(keys), "CI_low": d - q, "CI_high": d + q})
fam = pd.DataFrame(fam_rows)
fam.to_csv(OUT / "08d_ten_contrast_family_intervals.csv", index=False)

# ---------------------------------------------------------------- 5. learners
gen = GetMorganGenerator(radius=2, fpSize=1024)
bits = np.zeros((len(df), 1024), dtype=np.uint8)
for i, smi in enumerate(df.SMILES):
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(smi)), bits[i])

learn_rows = []
for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    for thr in THRESHOLDS:
        y = study.y[thr]
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED + split_i)
        lr = GridSearchCV(
            LogisticRegression(class_weight="balanced", max_iter=5000, random_state=SEED),
            {"C": [0.01, 0.1, 1.0, 10.0]}, scoring="roc_auc", cv=cv)
        lr.fit(bits[tr], y[tr])
        rf = RandomForestClassifier(n_estimators=500, class_weight="balanced_subsample",
                                    n_jobs=-1, random_state=SEED)
        rf.fit(bits[tr], y[tr])
        for name, score in (("LogisticL2", lr.decision_function(bits[te])),
                            ("RandomForest", rf.predict_proba(bits[te])[:, 1])):
            learn_rows.append({"Split": split, "Endpoint": thr, "Learner": name,
                               "AUROC": float(roc_auc_score(y[te], score)),
                               "SelectedC": float(lr.best_params_["C"]) if name == "LogisticL2"
                               else np.nan})
learn = pd.DataFrame(learn_rows)
learn.to_csv(OUT / "08e_learner_sensitivity.csv", index=False)

# ---------------------------------------------------------------- summary
sc = fam[fam.Bootstrap == "ScaffoldCluster"]
order = learn.pivot_table(index=["Split", "Learner"], columns="Endpoint", values="AUROC")
summary = {
    "wider_family_scaffold_ci_above_zero": int((sc.CI_low > 0).sum()),
    "wider_family_by_split": {s: int((g.CI_low > 0).sum()) for s, g in sc.groupby("Split")},
    "learner_settings_with_100_lowest": int(
        ((order[100] < order[11]) & (order[100] < order[1])).sum()),
    "learner_settings_total": int(len(order)),
}
(OUT / "08_SENSITIVITY_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

pd.set_option("display.width", 160)
print("== cross-endpoint (ECFP) ==")
print(cross[cross.Representation == "ECFP"].pivot_table(
    index=["Split", "TrainEndpoint"], columns="EvalEndpoint", values="AUROC").round(3))
print("\n== Brier skill ==")
print(brier.round(4).to_string(index=False))
print("\n== matched preprocessing ==")
print(meta.round(3).to_string(index=False))
print("\n== learners ==")
print(order.round(3))
print("\n", json.dumps(summary, indent=2))
