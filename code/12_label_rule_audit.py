"""Separate qualifier handling from aggregation in the ECOTOX-derived labels.

The upstream ApisTox ECOTOX step (replicated in ``11``) ignores the operator of
each LD50 record, labels the 11 ug/bee cut-off by requiring all records of a
CAS x route group to lie on one side (groups that straddle 11 are removed) and
labels the 1 and 100 ug/bee cut-offs by the group median. The compound takes the
labels of the eligible route with the lowest median.

This script re-labels the 441 ECOTOX-derived benchmark compounds under a 2 x 2
design that changes one factor at a time:

  rule A  operators ignored, group median          (upstream rule at 1 and 100)
  rule B  operators ignored, all records agree      (upstream rule at 11, strict <=)
  rule C  operators kept,    all records agree      (rule of the earlier audit)
  rule D  operators kept,    interval median        (bounds of the feasible median)

for each cut-off (1, 11, 100 ug/bee), three route definitions

  fixed     the historical label-determining CAS x route group
  eligible  minimum over the routes retained by the upstream 11 ug/bee filter
  all       minimum over all routes (no cut-off-specific route removal)

and three record filters

  all        records accepted by the upstream unit and value filters
  adult48    adult life stage, observed duration of exactly 2 days and a dose unit
             that is not per day
  adult48ai  adult48 restricted to active-ingredient dose reports.

Labels are 1 (LD50 <= cut-off), 0 (> cut-off), -1 (unresolved) or -2 (no
eligible record). No model is fitted and no benchmark label is changed.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from _common import RESULTS, ROOT, load_study_data

OUT = RESULTS / "label_audit"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / "output/external_data_audit/cache"
THRESHOLDS = (100, 11, 1)
RULES = ("A", "B", "C", "D")
ROUTE_DEFS = ("fixed", "eligible", "all")
FILTERS = ("all", "adult48", "adult48ai")
EPS = 1e-9  # turns open interval ends into strict comparisons

UNIT_FACTOR = {
    "AI ug/org": 1.0, "AI ug/org/d": 1.0, "ug/bee": 1.0, "ug/org": 1.0, "ug/org/d": 1.0,
    "AI ng/org": 1e-3, "AI ng/org/d": 1e-3, "ng/org": 1e-3,
    "AI mg/org": 1e3, "mg/bee": 1e3, "mg/org": 1e3, "pg/org": 1e-6,
}
ROUTE = {
    "Diet, unspecified": "Oral", "Drinking water": "Oral", "Food": "Oral",
    "Dermal": "Contact", "Direct application": "Contact", "Topical, general": "Contact",
    "Multiple routes between application groups": "Other", "Oral via capsule": "Other",
    "Spray, unspecified": "Other", "Environmental, unspecified": "Other",
}
ROUTE_ORDER = {"Contact": 0, "Oral": 1, "Other": 2}
EXACT_OPS = ("", "~")
KNOWN_OPS = {"", "~", ">", ">=", "<", "<="}


def standard_cas(value: str) -> str:
    s = str(int(float(value)))
    return f"{s[:-3]}-{s[-3:-1]}-{s[-1]}"


# --------------------------------------------------------------------------- records
raw = pd.read_csv(CACHE / "upstream_ecotox.csv", sep="|", dtype=str)
raw.columns = [c.strip() for c in raw.columns]
s = lambda c: raw[c].fillna("").str.strip()
rec = pd.DataFrame({
    "record": np.arange(len(raw)),
    "CAS": raw["CAS Number"].map(standard_cas),
    "exposure": s("Exposure Type"),
    "op": s("Observed Response Mean Op"),
    "value_raw": s("Observed Response Mean"),
    "unit": s("Observed Response Units"),
    "stage": s("Organism Lifestage"),
    "duration": s("Observed Duration (Days)"),
    "duration_op": s("Observed Duration Op (Days)"),
    "duration_unit": s("Observed Duration Units (Days)"),
    "conc_type": s("Conc 1 Type (Author)"),
    "range_min_op": s("Observed Response Min Op"),
    "range_max_op": s("Observed Response Max Op"),
    "reference": s("Reference Number"),
    "author": s("Author"),
    "year": s("Publication Year"),
    "title": s("Title"),
})
n_export = len(rec)
assert set(rec.op) <= KNOWN_OPS, set(rec.op) - KNOWN_OPS
# Upstream filters: convertible unit, reported value, mapped exposure route.
rec = rec[rec.unit.isin(UNIT_FACTOR) & (rec.value_raw != "NR")].copy()
rec["value"] = rec.value_raw.str.replace("/", "", regex=False).astype(float) * rec.unit.map(UNIT_FACTOR)
rec["route"] = rec.exposure.map(ROUTE)
rec = rec.dropna(subset=["route"]).copy()
rec["daily_unit"] = rec.unit.str.endswith("/d")
rec["adult"] = rec.stage.eq("Adult")
rec["duration_2d"] = rec.duration.eq("2") & rec.duration_unit.eq("Day(s)") & rec.duration_op.eq("")
rec["active_ingredient"] = rec.conc_type.eq("Active ingredient")
rec["adult48"] = rec.adult & rec.duration_2d & ~rec.daily_unit
rec["adult48ai"] = rec.adult48 & rec.active_ingredient
# Interval of each record (operators kept); open ends are shifted by EPS.
lo = np.where(rec.op.isin(EXACT_OPS) | rec.op.eq(">="), rec.value,
              np.where(rec.op.eq(">"), rec.value * (1 + EPS), 0.0))
hi = np.where(rec.op.isin(EXACT_OPS) | rec.op.eq("<="), rec.value,
              np.where(rec.op.eq("<"), rec.value * (1 - EPS), np.inf))
rec["lower"], rec["upper"] = lo, hi


def upstream_binary(values: np.ndarray):
    """Upstream 11 ug/bee rule, including its treatment of values equal to 11."""
    if values.min() <= 11 and values.max() <= 11:
        return 1
    if values.min() >= 11 and values.max() >= 11:
        return 0
    return "Unspecified"


def upstream_level(median: float) -> int:
    return 0 if median > 100 else (1 if median > 1 else 2)


def group_label(g: pd.DataFrame, rule: str, t: float) -> int:
    if g.empty:
        return -2
    v = g.value.to_numpy()
    if rule == "A":
        return int(np.median(v) <= t)
    if rule == "B":
        return 1 if (v <= t).all() else (0 if (v > t).all() else -1)
    lo, hi = g.lower.to_numpy(), g.upper.to_numpy()
    if rule == "C":
        pos, neg = hi <= t, lo > t
        return 1 if pos.all() else (0 if neg.all() else -1)
    # rule D: the feasible median lies between the medians of the lower and upper bounds.
    med_lo, med_hi = np.median(lo), np.median(hi)
    return 1 if med_hi <= t else (0 if med_lo > t else -1)


def combine_routes(labels: list[int]) -> int:
    """Compound label = label of the minimum over routes."""
    labels = [x for x in labels if x != -2]
    if not labels:
        return -2
    if 1 in labels:
        return 1
    return 0 if all(x == 0 for x in labels) else -1


# ------------------------------------------------- upstream replication (hard checks)
groups = []
for (cas, route), g in rec.groupby(["CAS", "route"]):
    v = g.value.to_numpy()
    groups.append({"CAS": cas, "route": route, "n_records": len(g), "median": float(np.median(v)),
                   "upstream_label": upstream_binary(v), "upstream_level": upstream_level(float(np.median(v)))})
groups = pd.DataFrame(groups)
groups["eligible"] = groups.upstream_label.ne("Unspecified")
elig = groups[groups.eligible].assign(route_order=lambda x: x.route.map(ROUTE_ORDER))
determining = (elig.sort_values(["CAS", "median", "route_order"])
               .groupby("CAS", as_index=False).first()[["CAS", "route", "median", "upstream_label", "upstream_level"]])

study = load_study_data()
bench = study.df.assign(Tier=study.tier, Index=np.arange(len(study.df)))
eco = bench[bench.source.eq("ECOTOX")].merge(determining, on="CAS", how="left")
assert len(eco) == 441, len(eco)
assert eco["route"].notna().all(), "every ECOTOX compound must have a determining group"
assert (eco.upstream_label.astype(int) == eco.label).all(), "binary label replication failed"
assert (eco.upstream_level.astype(int) == eco.ppdb_level).all(), "ternary level replication failed"
assert (eco.route == eco.toxicity_type).all(), "determining route differs from the distributed route"
# Values exactly equal to 11 in retained groups: upstream counts them on both sides.
exact11 = []
for cas, route in zip(eco.CAS, eco.route):
    g = rec[(rec.CAS == cas) & (rec.route == route)]
    if np.isclose(g.value, 11).any() and (g.value > 11).any() and not (g.value < 11).any():
        exact11.append({"CAS": cas, "route": route, "values": ";".join(f"{o}{x:g}" for o, x in zip(g.op, g.value))})

# -------------------------------------------------------------- relabelling matrix
original = {100: (eco.Tier >= 1).astype(int), 11: (eco.Tier >= 2).astype(int), 1: (eco.Tier >= 3).astype(int)}
rows = []
by_cas = {cas: g for cas, g in rec.groupby("CAS")}
eligible_routes = groups[groups.eligible].groupby("CAS").route.apply(set).to_dict()
for i, c in enumerate(eco.itertuples(index=False)):
    g_all = by_cas[c.CAS]
    for filt in FILTERS:
        g_f = g_all if filt == "all" else g_all[g_all[filt]]
        for route_def in ROUTE_DEFS:
            if route_def == "fixed":
                route_sets = [{c.route}]
            elif route_def == "eligible":
                route_sets = [eligible_routes.get(c.CAS, set())]
            else:
                route_sets = [set(g_all.route)]
            routes = sorted(route_sets[0], key=ROUTE_ORDER.get)
            for t in THRESHOLDS:
                for rule in RULES:
                    labs = [group_label(g_f[g_f.route == r], rule, t) for r in routes]
                    lab = labs[0] if route_def == "fixed" else combine_routes(labs)
                    rows.append({"Index": c.Index, "CAS": c.CAS, "Tier": int(c.Tier), "filter": filt,
                                 "routes": route_def, "threshold": t, "rule": rule, "label": lab,
                                 "original": int(original[t].iloc[i]),
                                 "n_records": int(sum(len(g_f[g_f.route == r]) for r in routes))})
labels = pd.DataFrame(rows)
# Index must address the benchmark row of each compound.
for t in THRESHOLDS:
    q = labels[labels.threshold == t]
    assert (q.original.to_numpy() == study.y[t][q.Index.to_numpy()]).all(), "Index does not address the benchmark row"
    assert (study.df.CAS.to_numpy()[q.Index.to_numpy()] == q.CAS.to_numpy()).all()
labels.to_csv(OUT / "02_compound_labels_by_rule.csv", index=False)

# The rule-A labels on the historical group with all records are the distributed labels.
chk = labels[(labels["filter"] == "all") & (labels.routes == "fixed") & (labels["rule"] == "A") & labels.threshold.isin([1, 100])]
assert (chk.label == chk.original).all(), "rule A must reproduce the distributed 1 and 100 labels"
chk = labels[(labels["filter"] == "all") & (labels.routes == "eligible") & (labels["rule"] == "A") & labels.threshold.isin([1, 100])]
assert (chk.label == chk.original).all(), "minimum over eligible routes must reproduce the distributed labels"

# --------------------------------------------------------------------- summaries
def counts(q: pd.DataFrame) -> dict:
    return {"n": int(len(q)), "positive": int((q.label == 1).sum()), "negative": int((q.label == 0).sum()),
            "unresolved": int((q.label == -1).sum()), "no_eligible_record": int((q.label == -2).sum()),
            "original_positive": int(q.original.sum()),
            "positive_to_negative": int(((q.original == 1) & (q.label == 0)).sum()),
            "negative_to_positive": int(((q.original == 0) & (q.label == 1)).sum()),
            "to_unresolved": int((q.label == -1).sum())}


summary_rows = []
for keys, q in labels.groupby(["filter", "routes", "threshold", "rule", "Tier"]):
    summary_rows.append(dict(zip(["filter", "routes", "threshold", "rule", "Tier"], keys)) | counts(q))
for keys, q in labels.groupby(["filter", "routes", "threshold", "rule"]):
    summary_rows.append(dict(zip(["filter", "routes", "threshold", "rule"], keys)) | {"Tier": "all"} | counts(q))
summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "03_rule_transition_counts.csv", index=False)

# Record-level ledger for the historical determining groups.
det_keys = set(zip(eco.CAS, eco.route))
ledger = rec[[k in det_keys for k in zip(rec.CAS, rec.route)]].merge(
    eco[["CAS", "name", "Tier", "Index"]], on="CAS")
ledger = ledger[["Index", "name", "CAS", "Tier", "route", "exposure", "op", "value_raw", "unit", "value", "lower",
                 "upper", "stage", "duration", "duration_op", "duration_unit", "daily_unit", "conc_type",
                 "adult48", "adult48ai", "range_min_op", "range_max_op", "reference", "author", "year", "title"]]
ledger.to_csv(OUT / "01_determining_group_records.csv", index=False)

# Heterogeneity of the historical determining groups (group level).
het = []
for (t, ), q in ledger.groupby(["Tier"]):
    gg = q.groupby("CAS")
    het.append({"Tier": int(t), "groups": int(gg.ngroups),
                "all_adult": int(gg.stage.apply(lambda x: x.eq("Adult").all()).sum()),
                "any_larva_or_pupa": int(gg.stage.apply(lambda x: x.isin(["Larva", "Pupa"]).any()).sum()),
                "stage_not_reported_only": int(gg.stage.apply(lambda x: x.eq("Not reported").all()).sum()),
                "all_duration_2d": int(gg.apply(lambda x: (x.duration.eq("2") & x.duration_unit.eq("Day(s)") & x.duration_op.eq("")).all(), include_groups=False).sum()),
                "any_daily_unit": int(gg.daily_unit.any().sum()),
                "all_adult48": int(gg.adult48.all().sum()),
                "any_adult48": int(gg.adult48.any().sum()),
                "any_right_censored": int(gg.op.apply(lambda x: x.isin([">", ">="]).any()).sum()),
                "records": int(len(q))})
het = pd.DataFrame(het)
het = pd.concat([het, pd.DataFrame([{"Tier": "all", **het.drop(columns="Tier").sum().to_dict()}])], ignore_index=True)
het.to_csv(OUT / "04_determining_group_heterogeneity.csv", index=False)

# Upstream cut-off-specific route removal and numeric straddlers.
up = []
all_cas = set(groups.CAS)
removed = groups[~groups.eligible]
for t in THRESHOLDS:
    straddle = []
    for (cas, route), g in rec[rec.CAS.isin(set(eco.CAS))].groupby(["CAS", "route"]):
        if (cas, route) in det_keys:
            v = g.value.to_numpy()
            straddle.append(bool((v <= t).any() and (v > t).any()))
    up.append({"threshold": t, "rule_in_upstream": "all records agree; straddling route groups removed" if t == 11
               else "group median; straddling route groups kept",
               "determining_groups": len(det_keys),
               "determining_groups_straddling_cutoff_numerically": int(sum(straddle))})
up = pd.DataFrame(up)
up.to_csv(OUT / "05_upstream_cutoff_rules.csv", index=False)

summary_json = {
    "ecotox_export_records": n_export,
    "records_after_upstream_filters": int(len(rec)),
    "cas_route_groups": int(len(groups)),
    "groups_removed_by_upstream_11_filter": int(len(removed)),
    "cas_with_removed_route": int(removed.CAS.nunique()),
    "cas_with_removed_route_retained_by_another_route": int(len(set(removed.CAS) & set(elig.CAS))),
    "benchmark_ecotox_compounds": int(len(eco)),
    "replication_label_level_route_match": int(len(eco)),
    "groups_with_exact_11_and_values_above_only": exact11,
    "records_with_range_operator_fields": int(((rec.range_min_op != "") | (rec.range_max_op != "")).sum()),
    "determining_group_heterogeneity": het.to_dict(orient="records"),
    "upstream_cutoff_rules": up.to_dict(orient="records"),
}
(OUT / "LABEL_AUDIT_SUMMARY.json").write_text(json.dumps(summary_json, indent=2, default=int), encoding="utf-8")

show = summary[(summary.Tier == 1) | (summary.Tier == "all")]
for filt in FILTERS:
    for route_def in ROUTE_DEFS:
        q = show[(show["filter"] == filt) & (show.routes == route_def)]
        print(f"\n== filter={filt} routes={route_def}")
        print(q[["threshold", "rule", "Tier", "n", "positive", "negative", "unresolved", "no_eligible_record",
                 "positive_to_negative", "negative_to_positive"]].to_string(index=False))
print(json.dumps({k: v for k, v in summary_json.items() if k not in ("determining_group_heterogeneity",)}, indent=1, default=int))
print(het.to_string(index=False))
