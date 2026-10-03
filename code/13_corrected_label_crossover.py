"""Corrected-label crossover for the ECOTOX-derived compounds.

Script 12 re-derives the labels of the 441 ECOTOX-derived compounds from the
upstream records with inequality operators kept. Here those labels replace the
distributed labels of these compounds (PPDB- and BPDB-derived compounds keep
their distributed labels, whose raw records are not available), and the fixed
SVM protocol is evaluated in a 2 x 2 crossover on the official test sets:

  train labels  original | corrected
  eval labels   original | corrected

For a full cohort, original-train scores are the saved held-out scores (a refit
is checked to reproduce them). When unresolved compounds are excluded, both
original-label and corrected-label models are refitted on the same retained
training identities. Every cell therefore uses the same training and test
members within a scenario/handling/split. Label definitions (scenarios):

  primary       operators kept; interval median at 1 and 100 ug/bee, all records
                agree at 11 ug/bee (each the upstream rule with operators kept);
                minimum over the upstream-eligible routes; all records
  conservative  operators kept; all records agree at every cut-off
  adult48       primary rule on adult, 48 h, non-daily-dose records only
  regulatory    labels from the route-matched OpenFoodTox and EPA records of the
                source audit (operators kept); conflicting sources are unresolved.
                This scenario also covers PPDB- and BPDB-derived compounds.
  combined      primary ECOTOX labels and regulatory labels; disagreement between
                them is unresolved

Unresolved corrected labels are handled three ways: ``keep`` (distributed label),
``exclude`` (compound removed from training and evaluation at every cut-off)
and ``negative`` (unresolved labels at 100 ug/bee set to negative; an extreme
scenario). Intervals are scaffold-cluster percentile intervals conditional on
the fitted models (2,000 draws, identical within a split and cohort).
"""
from __future__ import annotations

import hashlib
import json

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

OUT = RESULTS / "label_audit"
B = 2000
THRESHOLDS = (100, 11, 1)
SCENARIOS = {
    "primary": {"filter": "all", "routes": "eligible", "rule": {100: "D", 11: "C", 1: "D"}},
    "conservative": {"filter": "all", "routes": "eligible", "rule": {100: "C", 11: "C", 1: "C"}},
    "adult48": {"filter": "adult48", "routes": "eligible", "rule": {100: "D", 11: "C", 1: "D"}},
    "regulatory": None,
    "combined": None,
}
HANDLING = ("keep", "exclude", "negative")

cache = build_or_load_cache(force=False)
study = cache["study"]
kernels = cache["kernels"]
descriptors = np.asarray(cache["descriptors"])
scaffolds = np.asarray(cache["scaffolds"], dtype=object)
tier = study.tier.astype(int)
n = len(tier)
original = {t: study.y[t].astype(int) for t in THRESHOLDS}
saved = pd.read_csv(RESULTS / "primary" / "02_test_predictions.csv")
audit = pd.read_csv(OUT / "02_compound_labels_by_rule.csv")
EXT = ROOT / "output" / "external_data_audit"


def regulatory_labels() -> dict[int, np.ndarray]:
    """Route-matched OpenFoodTox/EPA labels: 1, 0, -1 (unresolved) or -9 (no record)."""
    conc = pd.read_csv(EXT / "same_route_label_concordance.csv")
    conc = conc[conc.source.isin(["OpenFoodTox", "PLOS2022"])]  # the 48 h screen is a subset of OpenFoodTox
    struct = pd.read_csv(EXT / "apistox_structures.csv")
    sm2idx = {sm: i for i, sm in enumerate(study.df.SMILES)}
    assert len(struct) == n and struct.SMILES.map(sm2idx).notna().all()
    key2idx = dict(zip(struct.inchikey, struct.SMILES.map(sm2idx).astype(int)))
    assert conc.inchikey.isin(key2idx).all()
    out = {}
    for t in THRESHOLDS:
        y = np.full(n, -9)
        for key, g in conc.groupby("inchikey"):
            vals = set(g[f"external_y_{t}"].astype(int))
            det = vals - {-1}
            # Every route-matched source must be determinate and agree.
            y[key2idx[key]] = det.pop() if len(det) == 1 and -1 not in vals else -1
            assert (g[f"apistox_y_{t}"] == original[t][key2idx[key]]).all(), "regulatory row misaligned"
        out[t] = y
    return out


def corrected_labels(scenario: str) -> dict[int, np.ndarray]:
    """Corrected label per compound and cut-off: 1, 0 or -1 (unresolved)."""
    if scenario in ("regulatory", "combined"):
        reg = regulatory_labels()
        eco = corrected_labels("primary") if scenario == "combined" else None
        out = {}
        for t in THRESHOLDS:
            y = original[t].copy()
            has = reg[t] != -9
            y[has] = reg[t][has]
            if eco is not None:
                is_eco = study.df.source.eq("ECOTOX").to_numpy()
                both = is_eco & has
                y[is_eco & ~has] = eco[t][is_eco & ~has]
                y[both] = np.where(eco[t][both] == reg[t][both], reg[t][both], -1)
            out[t] = y
        return out
    spec = SCENARIOS[scenario]
    out = {}
    for t in THRESHOLDS:
        y = original[t].copy()
        q = audit[(audit["filter"] == spec["filter"]) & (audit.routes == spec["routes"])
                  & (audit.threshold == t) & (audit["rule"] == spec["rule"][t])]
        assert q.Index.is_unique and len(q) == 441
        assert (q.original.to_numpy() == original[t][q.Index.to_numpy()]).all(), "label rows misaligned"
        lab = q.label.to_numpy()
        lab = np.where(lab == -2, -1, lab)  # no eligible record counts as unresolved
        y[q.Index.to_numpy()] = lab
        out[t] = y
    return out


def apply_handling(y: dict[int, np.ndarray], how: str) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """Return labels and the mask of compounds kept in training and evaluation."""
    unresolved = np.zeros(n, bool)
    for t in THRESHOLDS:
        unresolved |= y[t] == -1
    out = {t: y[t].copy() for t in THRESHOLDS}
    if how == "keep":
        for t in THRESHOLDS:
            out[t] = np.where(y[t] == -1, original[t], y[t])
        keep = np.ones(n, bool)
    elif how == "exclude":
        keep = ~unresolved
    else:  # negative at 100; other cut-offs keep the distributed label
        out[100] = np.where(y[100] == -1, 0, y[100])
        for t in (11, 1):
            out[t] = np.where(y[t] == -1, original[t], y[t])
        keep = np.ones(n, bool)
    if how != "exclude":
        assert all((out[t] >= 0).all() for t in THRESHOLDS)
    return out, keep


def fit_scores(rep: str, tr: np.ndarray, te: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fixed primary protocol: C=1, balanced class weights."""
    if rep == "Descriptors12":
        scaler = RobustScaler().fit(descriptors[tr])
        model = SVC(C=1, kernel="rbf", gamma="scale", class_weight="balanced")
        model.fit(scaler.transform(descriptors[tr]), y[tr])
        return model.decision_function(scaler.transform(descriptors[te]))
    k = kernels[rep]
    model = SVC(C=1, kernel="precomputed", class_weight="balanced")
    model.fit(np.asarray(k[np.ix_(tr, tr)]), y[tr])
    return model.decision_function(np.asarray(k[np.ix_(te, tr)]))


def saved_scores(split: str, rep: str, t: int, te: np.ndarray) -> np.ndarray:
    d = saved[(saved.Split == split) & (saved.Representation == rep) & (saved.Threshold == t)].sort_values("TestOrder")
    assert np.array_equal(d.Index.to_numpy(), te)
    return d.Score.to_numpy(float)


# Reproduction check: refitting with the original labels returns the saved scores.
for split, (tr, te) in study.splits.items():
    for rep in REPRESENTATIONS:
        for t in THRESHOLDS:
            diff = np.max(np.abs(fit_scores(rep, tr, te, original[t]) - saved_scores(split, rep, t, te)))
            assert diff < 1e-9, (split, rep, t, diff)

auc_rows, gap_rows, change_rows, prediction_frames = [], [], [], []
for scenario in SCENARIOS:
    ycorr = corrected_labels(scenario)
    for how in HANDLING:
        ylab, keep = apply_handling(ycorr, how)
        for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
            trk, tek = tr[keep[tr]], te[keep[te]]
            membership = {
                "n_train": len(trk), "n_test": len(tek),
                "train_index_sha256": hashlib.sha256(np.asarray(trk, dtype="<i8").tobytes()).hexdigest(),
                "test_index_sha256": hashlib.sha256(np.asarray(tek, dtype="<i8").tobytes()).hexdigest(),
            }
            pos = {t: np.flatnonzero(keep[te]) for t in THRESHOLDS}
            rng = np.random.default_rng(SEED + 13000 + 97 * split_i + 7 * HANDLING.index(how))
            uniq, inv = np.unique(scaffolds[tek], return_inverse=True)
            w = multinomial_weights(len(uniq), B, rng)[:, inv]
            changes = {t: int((ylab[t][tek] != original[t][tek]).sum()) for t in THRESHOLDS}
            changes_tr = {t: int((ylab[t][trk] != original[t][trk]).sum()) for t in THRESHOLDS}
            change_rows.append({"Scenario": scenario, "Unresolved": how, "Split": split,
                                **membership,
                                **{f"train_changed_le{t}": changes_tr[t] for t in THRESHOLDS},
                                **{f"test_changed_le{t}": changes[t] for t in THRESHOLDS},
                                **{f"test_pos_orig_le{t}": int(original[t][tek].sum()) for t in THRESHOLDS},
                                **{f"test_pos_corr_le{t}": int(ylab[t][tek].sum()) for t in THRESHOLDS}})
            for rep in REPRESENTATIONS:
                est, draws = {}, {}
                for t in THRESHOLDS:
                    s_orig = (fit_scores(rep, trk, tek, original[t]) if how == "exclude"
                              else saved_scores(split, rep, t, te)[pos[t]])
                    s_corr = fit_scores(rep, trk, tek, ylab[t])
                    prediction_frames.append(pd.DataFrame({
                        "Scenario": scenario, "Unresolved": how, "Split": split,
                        "Representation": rep, "Threshold": t,
                        "Index": tek, "CAS": study.df.CAS.to_numpy()[tek],
                        "OriginalLabel": original[t][tek], "CorrectedLabel": ylab[t][tek],
                        "OriginalTrainScore": s_orig, "CorrectedTrainScore": s_corr,
                        **membership,
                    }))
                    for train_lab, s in (("original", s_orig), ("corrected", s_corr)):
                        for eval_lab, yv in (("original", original[t][tek]), ("corrected", ylab[t][tek])):
                            key = (t, train_lab, eval_lab)
                            est[key] = float(roc_auc_score(yv, s))
                            draws[key] = tie_aware_weighted_auc(yv.astype(np.int8), s, w)
                            lo, hi = percentile_ci(draws[key])
                            auc_rows.append({"Scenario": scenario, "Unresolved": how, "Split": split,
                                             "Representation": rep, "Threshold": t, "TrainLabels": train_lab,
                                             "EvalLabels": eval_lab, **membership, "n_pos": int(yv.sum()),
                                             "B_requested": B, "B_valid": int(np.isfinite(draws[key]).sum()),
                                             "B_degenerate": int((~np.isfinite(draws[key])).sum()),
                                             "AUROC": est[key], "CI_low": lo, "CI_high": hi})
                # Severe-minus-broad gaps within each label combination, and their change.
                gaps, gdraws = {}, {}
                for train_lab in ("original", "corrected"):
                    for eval_lab in ("original", "corrected"):
                        for sev in (11, 1):
                            g = est[(sev, train_lab, eval_lab)] - est[(100, train_lab, eval_lab)]
                            gd = draws[(sev, train_lab, eval_lab)] - draws[(100, train_lab, eval_lab)]
                            gaps[(sev, train_lab, eval_lab)], gdraws[(sev, train_lab, eval_lab)] = g, gd
                            lo, hi = percentile_ci(gd)
                            gap_rows.append({"Scenario": scenario, "Unresolved": how, "Split": split,
                                             "Representation": rep, "Severe": sev, "TrainLabels": train_lab,
                                             "EvalLabels": eval_lab, **membership,
                                             "B_requested": B, "B_valid": int(np.isfinite(gd).sum()),
                                             "B_degenerate": int((~np.isfinite(gd)).sum()),
                                             "Gap": g, "CI_low": lo, "CI_high": hi})
                for sev in (11, 1):
                    base = ("original", "original")
                    for cell in (("original", "corrected"), ("corrected", "corrected")):
                        d = gdraws[(sev, *cell)] - gdraws[(sev, *base)]
                        lo, hi = percentile_ci(d)
                        g0 = gaps[(sev, *base)]
                        g1 = gaps[(sev, *cell)]
                        gap_rows.append({"Scenario": scenario, "Unresolved": how, "Split": split,
                                         "Representation": rep, "Severe": sev,
                                         "TrainLabels": f"change:{cell[0]}-vs-original",
                                         "EvalLabels": f"change:{cell[1]}-vs-original",
                                         **membership,
                                         "B_requested": B, "B_valid": int(np.isfinite(d).sum()),
                                         "B_degenerate": int((~np.isfinite(d)).sum()),
                                         "Gap": g1 - g0, "CI_low": lo, "CI_high": hi,
                                         "RelativeChange": (g1 - g0) / g0 if g0 != 0 else np.nan})

auc = pd.DataFrame(auc_rows)
gap = pd.DataFrame(gap_rows)
chg = pd.DataFrame(change_rows)
auc.to_csv(OUT / "06_crossover_auroc.csv", index=False)
gap.to_csv(OUT / "07_crossover_gap.csv", index=False)
chg.to_csv(OUT / "08_crossover_label_changes.csv", index=False)
pd.concat(prediction_frames, ignore_index=True).to_csv(OUT / "09_crossover_predictions.csv.gz",
                                                      index=False, compression={"method": "gzip", "mtime": 0})

# Readable summary: ECFP, severe = 11.
q = gap[(gap.Representation == "ECFP") & (gap.Severe == 11)]
summary = {}
for (scenario, how, split), g in q.groupby(["Scenario", "Unresolved", "Split"]):
    cell = lambda a, b: g[(g.TrainLabels == a) & (g.EvalLabels == b)].iloc[0]
    o = cell("original", "original")
    e = cell("original", "corrected")
    c = cell("corrected", "corrected")
    dc = cell("change:corrected-vs-original", "change:corrected-vs-original")
    summary[f"{scenario}/{how}/{split}"] = {
        "gap_original": round(o.Gap, 4), "gap_eval_corrected": round(e.Gap, 4),
        "gap_train_eval_corrected": round(c.Gap, 4),
        "change_full": [round(dc.Gap, 4), round(dc.CI_low, 4), round(dc.CI_high, 4)],
    }
(OUT / "CROSSOVER_SUMMARY.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
print(chg.to_string(index=False))
print(json.dumps(summary, indent=1))
