"""Trace inequality qualifiers through the upstream ECOTOX labelling step.

The upstream ApisTox ECOTOX loader reads the ``Observed Response Mean`` column
but not its operator column (``Observed Response Mean Op``); the PPDB/BPDB parser
strips leading ``<``, ``>`` and ``=`` characters. A limit-test result such as
``>100`` therefore enters the ternary rule as 100.0, which the moderate class
(1 < LD50 <= 100) contains.

This script replicates the upstream ECOTOX labelling from the cached ApisTox
ECOTOX export, verifies that the replicated labels equal the distributed
benchmark labels, and then re-labels the 100 ug/bee cut-off with the operators
kept. No model is fitted and no benchmark label is changed.

Inputs (cached copies of the public ApisTox v2 release and code):
  output/external_data_audit/cache/upstream_ecotox.csv
  output/external_data_audit/cache/upstream_ecotox.py (reference only)
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from _common import RESULTS, ROOT, load_study_data

OUT = RESULTS / "tier_boundary"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / "output/external_data_audit/cache"

# Unit conversions accepted by the upstream loader; every other unit is removed there.
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


def standard_cas(value: str) -> str:
    s = str(int(float(value)))
    return f"{s[:-3]}-{s[-3:-1]}-{s[-1]}"


def ppdb_level(median: float) -> int:
    if median > 100:
        return 0
    if 1 < median <= 100:
        return 1
    return 2


def binary_label(values: np.ndarray):
    if values.min() <= 11 and values.max() <= 11:
        return 1
    if values.min() >= 11 and values.max() >= 11:
        return 0
    return "Unspecified"


def interval_label_100(ops: pd.Series, values: pd.Series) -> int:
    """Qualifier-aware label at 100 ug/bee: 1, 0, or -1 (unresolved).

    All records of the group must be determinate and agree, as in the regulatory
    source audit.
    """
    labels = []
    for op, v in zip(ops, values):
        if op in ("", "=", "~"):
            labels.append(int(v <= 100))
        elif op == ">":
            labels.append(0 if v >= 100 else -1)
        elif op == ">=":
            labels.append(0 if v > 100 else -1)
        elif op in ("<", "<="):
            labels.append(1 if v <= 100 else -1)
        else:
            labels.append(-1)
    labels = set(labels)
    return labels.pop() if len(labels) == 1 and -1 not in labels else -1


raw = pd.read_csv(CACHE / "upstream_ecotox.csv", sep="|", dtype=str)
raw.columns = [c.strip() for c in raw.columns]
d = pd.DataFrame({
    "CAS": raw["CAS Number"].map(standard_cas),
    "exposure": raw["Exposure Type"].str.strip(),
    "op": raw["Observed Response Mean Op"].fillna("").str.strip(),
    "value": raw["Observed Response Mean"].str.strip(),
    "unit": raw["Observed Response Units"].str.strip(),
})
d = d[d.unit.isin(UNIT_FACTOR)]
d = d[d.value != "NR"].copy()
d["value"] = d.value.str.replace("/", "", regex=False).astype(float) * d.unit.map(UNIT_FACTOR)
d["route"] = d.exposure.map(ROUTE)
d = d.dropna(subset=["route"])

groups = []
for (cas, route), g in d.groupby(["CAS", "route"]):
    values = g.value.to_numpy()
    groups.append({
        "CAS": cas, "route": route, "n_records": len(g),
        "median": float(np.median(values)), "label": binary_label(values),
        "ppdb_level": ppdb_level(float(np.median(values))),
        "n_right_censored": int(g.op.isin([">", ">="]).sum()),
        "n_right_censored_ge100": int((g.op.isin([">", ">="]) & (g.value >= 100)).sum()),
        "interval_label_100": interval_label_100(g.op, g.value),
    })
groups = pd.DataFrame(groups)
groups = groups[groups.label != "Unspecified"].copy()
groups["label"] = groups.label.astype(int)
# Upstream rule: the route with the lowest median defines the compound (Contact, Oral, Other order).
order = {"Contact": 0, "Oral": 1, "Other": 2}
groups["route_order"] = groups.route.map(order)
determining = (groups.sort_values(["CAS", "median", "route_order"])
               .groupby("CAS", as_index=False).first())

study = load_study_data()
bench = study.df.assign(Tier=study.tier)
eco = bench[bench.source == "ECOTOX"].merge(determining, on="CAS", how="left",
                                           suffixes=("", "_replicated"))
matched = eco.dropna(subset=["median"]).copy()
matched["replicated_label_matches"] = matched.label == matched.label_replicated
matched["replicated_level_matches"] = matched.ppdb_level == matched.ppdb_level_replicated
matched["median_exactly_100"] = np.isclose(matched["median"], 100.0)
matched.to_csv(OUT / "06_ecotox_qualifier_propagation.csv", index=False,
               columns=["name", "CAS", "Tier", "label", "ppdb_level", "route", "n_records", "median",
                        "n_right_censored", "n_right_censored_ge100", "interval_label_100",
                        "median_exactly_100", "replicated_label_matches", "replicated_level_matches"])

consistent = matched[matched.replicated_label_matches & matched.replicated_level_matches]
rows = []
for t in range(4):
    q = consistent[consistent.Tier == t]
    rows.append({
        "Tier": t, "n_ecotox_compounds": int(len(q)),
        "any_right_censored": int((q.n_right_censored > 0).sum()),
        "right_censored_ge100": int((q.n_right_censored_ge100 > 0).sum()),
        "median_exactly_100": int(q.median_exactly_100.sum()),
        "interval_label_100_negative": int((q.interval_label_100 == 0).sum()),
        "interval_label_100_unresolved": int((q.interval_label_100 == -1).sum()),
        "interval_label_100_positive": int((q.interval_label_100 == 1).sum()),
    })
by_tier = pd.DataFrame(rows)
by_tier.to_csv(OUT / "07_ecotox_qualifier_by_tier.csv", index=False)
by_tier.to_csv(ROOT / "source_data" / "figure_qualifier_propagation.csv", index=False)

summary = {
    "ecotox_export_records": int(len(raw)),
    "records_with_right_censoring_operator": int(raw["Observed Response Mean Op"].fillna("").str.strip().isin([">", ">="]).sum()),
    "benchmark_ecotox_compounds": int(len(eco)),
    "replicated_compounds": int(len(matched)),
    "replication_label_and_level_match": int(len(consistent)),
    "tier1": by_tier[by_tier.Tier == 1].iloc[0].to_dict(),
    "by_tier": by_tier.to_dict(orient="records"),
}
(OUT / "QUALIFIER_PROPAGATION_SUMMARY.json").write_text(json.dumps(summary, indent=2, default=int),
                                                      encoding="utf-8")
print(json.dumps(summary, indent=1, default=int))
