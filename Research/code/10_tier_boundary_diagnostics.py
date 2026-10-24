"""Tier-boundary diagnostics for the broad (<=100 ug/bee) endpoint.

Three analyses locate where the broad-endpoint AUROC deficit arises:

1. Tier-pair AUROC. Every endpoint AUROC is a weighted average of tier-pair
   ranking probabilities (Supplementary Methods, Eq. S2). The saved held-out
   scores of every representation and training endpoint are evaluated on each
   pair of analytical tiers, with the composition weight of that pair.
2. Tier1 exclusion. Tier1 (the (11, 100] ug/bee annotations) is removed from
   evaluation only, and then from both training and evaluation. With Tier1
   removed, the 100 and 11 ug/bee cut-offs define the same task (Tier0 versus
   Tier2+Tier3). The broad-minus-severe gap is decomposed exactly into an
   evaluation component, a training component and a residual.
3. Source concordance by tier. The route-matched regulatory concordance records
   are broken down by the benchmark tier and by the benchmark record source.

Original models are not refitted: their saved scores are read from
results/primary. Only the Tier1-excluded models are new fits, using the same
fixed SVM protocol. Bootstrap draws reproduce the primary analysis draws
(same seeds and call order), so intervals are paired with Table 2.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import RobustScaler
from sklearn.svm import SVC

from _common import (
    REPRESENTATIONS,
    RESULTS,
    ROOT,
    SEED,
    build_or_load_cache,
    multinomial_weights,
    percentile_ci,
    tie_aware_weighted_auc,
)

OUT = RESULTS / "tier_boundary"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE_DATA = ROOT / "source_data"
B = 5000  # identical to the primary paired bootstrap

cache = build_or_load_cache(force=False)
study = cache["study"]
kernels = cache["kernels"]
descriptors = np.asarray(cache["descriptors"])
scaffolds = np.asarray(cache["scaffolds"], dtype=object)
tier = study.tier.astype(int)

saved = pd.read_csv(RESULTS / "primary" / "02_test_predictions.csv")


def saved_scores(split: str, rep: str, thr: int, te: np.ndarray) -> np.ndarray:
    d = saved[(saved.Split == split) & (saved.Representation == rep) & (saved.Threshold == thr)]
    d = d.sort_values("TestOrder")
    if not np.array_equal(d.Index.to_numpy(), te):
        raise AssertionError(f"Saved prediction order differs for {split}/{rep}/{thr}")
    return d.Score.to_numpy(dtype=float)


def fit_without_tier1(rep: str, tr: np.ndarray, te: np.ndarray) -> np.ndarray:
    """Fit Tier0 versus Tier2+Tier3 on training compounds outside Tier1."""
    trc = tr[tier[tr] != 1]
    tec = te[tier[te] != 1]
    y = (tier >= 2).astype(np.int8)
    if rep == "Descriptors12":
        scaler = RobustScaler().fit(descriptors[trc])
        model = SVC(C=1, kernel="rbf", gamma="scale", class_weight="balanced")
        model.fit(scaler.transform(descriptors[trc]), y[trc])
        return model.decision_function(scaler.transform(descriptors[tec]))
    kernel = kernels[rep]
    model = SVC(C=1, kernel="precomputed", class_weight="balanced")
    model.fit(np.asarray(kernel[np.ix_(trc, trc)]), y[trc])
    return model.decision_function(np.asarray(kernel[np.ix_(tec, trc)]))


def bootstrap_weights(split_i: int, te: np.ndarray) -> dict[str, np.ndarray]:
    """Reproduce the molecule and scaffold-cluster draws of 01_primary_endpoint_models."""
    rng = np.random.default_rng(SEED + 101 * (split_i + 1))
    molecule = multinomial_weights(len(te), B, rng)
    unique_scaffolds, inverse = np.unique(scaffolds[te], return_inverse=True)
    scaffold = multinomial_weights(len(unique_scaffolds), B, rng)[:, inverse]
    return {"Molecule": molecule, "ScaffoldCluster": scaffold}


def auc_with_draws(y: np.ndarray, s: np.ndarray, w: np.ndarray) -> tuple[float, np.ndarray]:
    return float(roc_auc_score(y, s)), tie_aware_weighted_auc(y, s, w)


pair_rows, excl_rows, contrast_rows, refit_rows = [], [], [], []
PAIRS = [(hi, lo) for hi in range(1, 4) for lo in range(hi)]
KTHR = {100: 1, 11: 2, 1: 3}

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    tt = tier[te]
    counts = np.bincount(tt, minlength=4)
    weights = bootstrap_weights(split_i, te)
    keep = tt != 1  # Tier1-excluded evaluation cohort
    for rep in REPRESENTATIONS:
        scores = {thr: saved_scores(split, rep, thr, te) for thr in (100, 11, 1)}

        # ---------------- 1. tier-pair AUROC ---------------------------------
        for trained, s in scores.items():
            for hi, lo in PAIRS:
                m = (tt == hi) | (tt == lo)
                y = (tt[m] == hi).astype(np.int8)
                est, draws = auc_with_draws(y, s[m], weights["ScaffoldCluster"][:, m])
                low, high = percentile_ci(draws)
                row = {
                    "Split": split, "Representation": rep, "TrainedFor": trained,
                    "HigherTier": hi, "LowerTier": lo, "Pair": f"T{hi}>T{lo}",
                    "n_higher": int(counts[hi]), "n_lower": int(counts[lo]),
                    "AUROC": est, "ScaffoldCI_low": low, "ScaffoldCI_high": high,
                }
                # Composition weight of this pair in each endpoint AUROC (Eq. S2).
                for thr, k in KTHR.items():
                    npos = counts[k:].sum()
                    nneg = counts[:k].sum()
                    inside = hi >= k > lo
                    row[f"Weight_le{thr}"] = counts[hi] * counts[lo] / (npos * nneg) if inside else 0.0
                pair_rows.append(row)

        # ---------------- 2. Tier1 exclusion ---------------------------------
        clean_scores = fit_without_tier1(rep, tr, te)
        y_clean = (tt[keep] >= 2).astype(np.int8)
        for idx, tier_i, score in zip(te[keep], tt[keep], clean_scores):
            refit_rows.append({"Split": split, "Representation": rep, "Index": int(idx),
                               "Tier": int(tier_i), "Y": int(tier_i >= 2), "Score": float(score)})
        quantities = {}
        for mode, w in weights.items():
            est, draws = {}, {}
            est["A100_full"], draws["A100_full"] = auc_with_draws(
                (tt >= 1).astype(np.int8), scores[100], w)
            est["A11_full"], draws["A11_full"] = auc_with_draws(
                (tt >= 2).astype(np.int8), scores[11], w)
            est["A1_full"], draws["A1_full"] = auc_with_draws(
                (tt >= 3).astype(np.int8), scores[1], w)
            est["A100_noT1_eval"], draws["A100_noT1_eval"] = auc_with_draws(
                y_clean, scores[100][keep], w[:, keep])
            est["A11_noT1_eval"], draws["A11_noT1_eval"] = auc_with_draws(
                y_clean, scores[11][keep], w[:, keep])
            est["Aclean_noT1_train_eval"], draws["Aclean_noT1_train_eval"] = auc_with_draws(
                y_clean, clean_scores, w[:, keep])
            quantities[mode] = (est, draws)
            for name, value in est.items():
                low, high = percentile_ci(draws[name])
                excl_rows.append({
                    "Split": split, "Representation": rep, "Bootstrap": mode,
                    "Quantity": name, "AUROC": value, "CI_low": low, "CI_high": high,
                    "n_test": int(keep.sum()) if "noT1" in name else int(len(te)),
                })
            # Exact decomposition: A11 - A100 = eval + train + residual.
            parts = {
                "Gap_11_minus_100": ("A11_full", "A100_full"),
                "Evaluation_component": ("A100_noT1_eval", "A100_full"),
                "Training_component": ("Aclean_noT1_train_eval", "A100_noT1_eval"),
                "Tier1_total_component": ("Aclean_noT1_train_eval", "A100_full"),
                "Residual_11_minus_clean": ("A11_full", "Aclean_noT1_train_eval"),
                "Residual_1_minus_clean": ("A1_full", "Aclean_noT1_train_eval"),
                "Tier1_share_of_gap_11": None,
            }
            for name, pair in parts.items():
                if pair is None:
                    gap = est["A11_full"] - est["A100_full"]
                    removed = est["Aclean_noT1_train_eval"] - est["A100_full"]
                    contrast_rows.append({
                        "Split": split, "Representation": rep, "Bootstrap": mode,
                        "Contrast": name, "Estimate": removed / gap if gap else math.nan,
                        "CI_low": math.nan, "CI_high": math.nan, "B_valid": math.nan,
                    })
                    continue
                a, b = pair
                values = draws[a] - draws[b]
                values = values[np.isfinite(values)]
                low, high = percentile_ci(values)
                contrast_rows.append({
                    "Split": split, "Representation": rep, "Bootstrap": mode,
                    "Contrast": name, "Estimate": est[a] - est[b],
                    "CI_low": low, "CI_high": high, "B_valid": int(len(values)),
                })

pairs = pd.DataFrame(pair_rows)
pairs.to_csv(OUT / "01_tier_pair_auroc.csv", index=False)
excl = pd.DataFrame(excl_rows)
excl.to_csv(OUT / "02_tier1_exclusion_auroc.csv", index=False)
contrasts = pd.DataFrame(contrast_rows)
contrasts.to_csv(OUT / "03_tier1_exclusion_decomposition.csv", index=False)
pd.DataFrame(refit_rows).to_csv(OUT / "02b_tier1_refit_test_predictions.csv", index=False)

# ---------------- 3. regulatory concordance by tier and record source -------
conc = pd.read_csv(ROOT / "output/external_data_audit/same_route_label_concordance.csv")
det = conc[(conc[["external_y_1", "external_y_11", "external_y_100"]] >= 0).all(axis=1)].copy()
det["tier"] = det.apistox_y_100 + det.apistox_y_11 + det.apistox_y_1
SOURCE_NAME = {"OpenFoodTox": "OFT", "OpenFoodTox_48h_screen": "OFT 48 h", "PLOS2022": "EPA"}


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (centre - half, centre + half)


tier_rows, source_rows = [], []
full_share = np.bincount(tier, minlength=4) / len(tier)
for src, g in det.groupby("source", sort=False):
    for t in range(4):
        q = g[g.tier == t]
        dis = q[q.agreement_100 == 0]
        lo, hi = wilson(int((q.agreement_100 == 1).sum()), len(q))
        tier_rows.append({
            "Source": SOURCE_NAME[src], "Tier": t, "n": int(len(q)),
            "Cohort_share": len(q) / len(g), "Benchmark_share": float(full_share[t]),
            "Agree_100": int((q.agreement_100 == 1).sum()), "Disagree_100": int(len(dis)),
            "Agreement_rate_100": float((q.agreement_100 == 1).mean()) if len(q) else math.nan,
            "Agreement_Wilson_low": lo, "Agreement_Wilson_high": hi,
            "Disagree_benchmark_pos_source_neg": int(((dis.apistox_y_100 == 1) & (dis.external_y_100 == 0)).sum()),
            "Disagree_benchmark_neg_source_pos": int(((dis.apistox_y_100 == 0) & (dis.external_y_100 == 1)).sum()),
            "Disagree_with_right_censored": int(dis.contains_right_censored_100.sum()),
        })
    for rs, q in g.groupby("apistox_source"):
        source_rows.append({
            "Source": SOURCE_NAME[src], "Benchmark_record_source": rs, "n": int(len(q)),
            "Disagree_100": int((q.agreement_100 == 0).sum()),
            "Disagree_100_Tier1": int(((q.agreement_100 == 0) & (q.tier == 1)).sum()),
            "n_Tier1": int((q.tier == 1).sum()),
        })
by_tier = pd.DataFrame(tier_rows)
by_tier.to_csv(OUT / "04_concordance_by_tier.csv", index=False)
by_source = pd.DataFrame(source_rows)
by_source.to_csv(OUT / "05_concordance_by_record_source.csv", index=False)

# ---------------- figure source data and summary ----------------------------
pairs[(pairs.Representation == "ECFP")].to_csv(SOURCE_DATA / "figure_tier_pair_auroc.csv", index=False)
excl[(excl.Representation == "ECFP") & (excl.Bootstrap == "ScaffoldCluster")].to_csv(
    SOURCE_DATA / "figure_tier1_exclusion.csv", index=False)
by_tier.to_csv(SOURCE_DATA / "figure_concordance_by_tier.csv", index=False)


def pick(df, **kw):
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    if len(q) != 1:
        raise AssertionError(kw)
    return q.iloc[0]


summary = {"bootstrap_draws": B, "resampling": "scaffold-cluster draws identical to the primary analysis",
           "new_model_fits": 15, "splits": {}}
for split in ("Random", "MaxMin", "Time"):
    e = lambda q: float(pick(excl, Split=split, Representation="ECFP", Bootstrap="ScaffoldCluster", Quantity=q).AUROC)
    c = lambda q: pick(contrasts, Split=split, Representation="ECFP", Bootstrap="ScaffoldCluster", Contrast=q)
    p = lambda hi, lo, t: pick(pairs, Split=split, Representation="ECFP", TrainedFor=t, HigherTier=hi, LowerTier=lo)
    summary["splits"][split] = {
        "A100_full": e("A100_full"), "A100_noT1_eval": e("A100_noT1_eval"),
        "Aclean_noT1_train_eval": e("Aclean_noT1_train_eval"),
        "A11_full": e("A11_full"), "A1_full": e("A1_full"),
        "gap_11_minus_100": float(c("Gap_11_minus_100").Estimate),
        "evaluation_component": [float(c("Evaluation_component").Estimate), float(c("Evaluation_component").CI_low), float(c("Evaluation_component").CI_high)],
        "training_component": [float(c("Training_component").Estimate), float(c("Training_component").CI_low), float(c("Training_component").CI_high)],
        "tier1_total_component": [float(c("Tier1_total_component").Estimate), float(c("Tier1_total_component").CI_low), float(c("Tier1_total_component").CI_high)],
        "residual_11": [float(c("Residual_11_minus_clean").Estimate), float(c("Residual_11_minus_clean").CI_low), float(c("Residual_11_minus_clean").CI_high)],
        "tier1_share_of_gap_11": float(c("Tier1_share_of_gap_11").Estimate),
        "T1_vs_T0_model100": [float(p(1, 0, 100).AUROC), float(p(1, 0, 100).ScaffoldCI_low), float(p(1, 0, 100).ScaffoldCI_high)],
        "T3_vs_T0_model100": float(p(3, 0, 100).AUROC),
        "T1_vs_T0_weight_le100": float(p(1, 0, 100).Weight_le100),
        "T1_vs_T0_model11": float(p(1, 0, 11).AUROC), "T1_vs_T0_model1": float(p(1, 0, 1).AUROC),
        "T2_vs_T1_model11": float(p(2, 1, 11).AUROC),
    }
summary["concordance_tier1"] = {
    r.Source: {"n": int(r.n), "agree": int(r.Agree_100), "rate": float(r.Agreement_rate_100),
               "wilson": [float(r.Agreement_Wilson_low), float(r.Agreement_Wilson_high)]}
    for r in by_tier[by_tier.Tier == 1].itertuples()
}
summary["disagreements_by_record_source"] = {
    s: {r.Benchmark_record_source: int(r.Disagree_100) for r in g.itertuples()}
    for s, g in by_source.groupby("Source")
}
(OUT / "TIER_BOUNDARY_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

pd.set_option("display.width", 200)
print(json.dumps(summary, indent=1))
print(by_tier.round(3).to_string(index=False))
print(by_source.to_string(index=False))
