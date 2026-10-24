"""Deterministic cached-record audit; offline by default and no fitted models.

Preserve original record ordinals, canonical row/file hashes and available
identifiers. A source record replay establishes extraction behaviour; a
chemical/route match between databases does not establish the same study.
Run from any directory with the pinned analysis environment.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import re
import argparse
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/record_audit"
EXT = ROOT / "output/external_data_audit"
THRESHOLDS = (1, 11, 100)
FACTORS = {
    "AI ug/org": "1", "AI ug/org/d": "1", "ug/bee": "1", "ug/org": "1", "ug/org/d": "1",
    "AI ng/org": ".001", "AI ng/org/d": ".001", "ng/org": ".001", "AI mg/org": "1000",
    "mg/bee": "1000", "mg/org": "1000", "pg/org": ".000001",
}
ROUTES = {"Diet, unspecified": "Oral", "Drinking water": "Oral", "Food": "Oral",
          "Dermal": "Contact", "Direct application": "Contact", "Topical, general": "Contact",
          "Multiple routes between application groups": "Other", "Oral via capsule": "Other",
          "Spray, unspecified": "Other", "Environmental, unspecified": "Other"}
INPUTS = ["data/raw/apistox.csv", "output/external_data_audit/cache/upstream_ecotox.csv",
          "output/external_data_audit/apistox_structures.csv",
          "output/external_data_audit/oft_acute_endpoints_aligned.csv",
          "output/external_data_audit/oft_bee_source_rows.csv",
          "output/external_data_audit/plos_adult_acute_endpoints_aligned.csv",
          "output/external_data_audit/plos_mrid_citations.csv",
          "output/external_data_audit/same_route_label_concordance.csv",
          "results/label_audit/02_compound_labels_by_rule.csv",
          "data/external_raw/OFT3.0 export repository.xlsx",
          "data/external_raw/journal.pone.0265962.s001.xlsx",
          "data/external_raw/journal.pone.0265962.s004.xlsx"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row_sha(row):
    """UTF-8 JSON [column,value] pairs in original column order, missing=''."""
    payload = json.dumps(list(row.items()), ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read(path, sep=","):
    return pd.read_csv(ROOT / path, dtype=str, keep_default_na=False, sep=sep)


def save(df, name):
    df.to_csv(OUT / name, index=False, lineterminator="\n")


def cas_format(s):
    s = str(int(s))
    return f"{s[:-3]}-{s[-3:-1]}-{s[-1]}"


def median(values):
    v = sorted(values)
    return v[len(v)//2] if len(v) % 2 else (v[len(v)//2-1] + v[len(v)//2]) / 2


def interval(op, value):
    # Approximate ('~') values are retained as point values, matching script 12;
    # this convention is flagged explicitly in the exported record ledger.
    if op in ("", "=", "~"):
        return value, True, value, True
    if op in (">", ">="):
        return value, op == ">=", math.inf, False
    if op in ("<", "<="):
        return Fraction(0), True, value, op == "<="
    raise ValueError(op)


def interval_label(iv, t):
    lo, lc, hi, hc = iv
    return 1 if hi <= t else 0 if lo > t or (lo == t and not lc) else -1


def median_interval(intervals):
    """Exact rational bounds and attainability; no epsilon perturbations.

    At a tied lower bound, rank k is attainable when at least k records have
    a smaller lower bound or a closed lower bound at that value. Upper-bound
    attainability is the symmetric count from the top. The arithmetic median
    for even n attains its bound only if both central ranks can attain theirs.
    """
    n = len(intervals)
    ranks = [n//2 + 1] if n % 2 else [n//2, n//2 + 1]
    lowers, uppers = sorted(x[0] for x in intervals), sorted(x[2] for x in intervals)
    lo = sum((lowers[k-1] for k in ranks), Fraction(0)) / len(ranks)
    hu = [uppers[k-1] for k in ranks]
    hi = math.inf if math.inf in hu else sum(hu, Fraction(0)) / len(ranks)
    lc = all(sum(a < lowers[k-1] or (a == lowers[k-1] and ac)
                 for a, ac, b, bc in intervals) >= k for k in ranks)
    hc = hi != math.inf and all(sum(b > uppers[k-1] or (b == uppers[k-1] and bc)
                                 for a, ac, b, bc in intervals) >= n-k+1 for k in ranks)
    return lo, lc, hi, hc


def unanimous(labels):
    return labels[0] if len(set(labels)) == 1 and labels[0] in (0, 1) else -1


def fmt(v):
    return "inf" if v == math.inf else f"{float(v):g}"


def render(iv):
    lo, lc, hi, hc = iv
    return f"{'[' if lc else '('}{fmt(lo)}, {fmt(hi)}{']' if hc else ')'}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-upstream", action="store_true",
                        help="Read three files at the pinned upstream commit and save hash comparisons; never replace caches")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    # Boundary checks include even-n and tied open/closed bounds, which an
    # arbitrary epsilon need not handle correctly for arbitrary future data.
    assert interval_label(median_interval([interval(">", Fraction(100))]), 100) == 0
    assert interval_label(median_interval([interval(">=", Fraction(100))]), 100) == -1
    assert interval_label(median_interval([interval("", Fraction(90)), interval("", Fraction(110))]), 100) == 1
    assert interval_label(median_interval([interval(">", Fraction(90)), interval("", Fraction(110))]), 100) == 0
    assert median_interval([interval(">", Fraction(100)), interval("", Fraction(100)), interval("", Fraction(90))])[1]

    hashes = {p: sha(ROOT / p) for p in INPUTS}
    bench = read("data/raw/apistox.csv")
    bench["Index"] = range(len(bench))
    bench["Tier"] = [0 if int(p) == 0 else 3 if int(p) == 2 else 1 + int(y)
                     for p, y in zip(bench.ppdb_level, bench.label)]
    structs = read("output/external_data_audit/apistox_structures.csv")
    bench = bench.merge(structs[["SMILES", "inchikey"]], on="SMILES", validate="one_to_one")
    eco = bench[bench.source.eq("ECOTOX")]
    eco_by_key = {(r.CAS, r.toxicity_type): r for r in eco.itertuples(index=False)}
    bench_by_key = {r.inchikey: r for r in bench.itertuples(index=False)}
    assert len(bench_by_key) == len(bench)

    source = "output/external_data_audit/cache/upstream_ecotox.csv"
    raw = read(source, sep="|")
    records, groups = [], {}
    for ordinal, row in raw.iterrows():
        s = {k.strip(): v.strip() for k, v in row.items()}
        unit, raw_value = s["Observed Response Units"], s["Observed Response Mean"]
        route = ROUTES.get(s["Exposure Type"])
        key = (cas_format(s["CAS Number"]), route)
        if unit not in FACTORS or raw_value == "NR" or key not in eco_by_key:
            continue
        c = eco_by_key[key]
        v = Fraction(raw_value.replace("/", "")) * Fraction(FACTORS[unit])
        op = s["Observed Response Mean Op"]
        iv = interval(op, v)
        r = {"Index": c.Index, "CAS": c.CAS, "name": c.name, "Tier": c.Tier,
             "inchikey": c.inchikey, "route": route, "source_path": source,
             "source_file_sha256": hashes[source], "source_data_row_0based": ordinal,
             "source_row_sha256": row_sha(row), "ecotox_reference_number": s["Reference Number"],
             "ecotox_test_id": "not supplied in cached export", "ecotox_result_id": "not supplied in cached export",
             "MRID": "not supplied in cached export", "reference_is_study_identifier": False,
             "reference_title": s["Title"], "reference_author": s["Author"], "reference_year": s["Publication Year"],
             "exposure": s["Exposure Type"], "life_stage": s["Organism Lifestage"],
             "duration_value": s["Observed Duration (Days)"], "duration_operator": s["Observed Duration Op (Days)"],
             "duration_unit": s["Observed Duration Units (Days)"], "chemical_grade": s["Chemical Grade"],
             "purity_operator": s["Chemical Purity Mean Op"], "purity_percent": s["Chemical Purity Mean(%)"],
             "test_material_dose_basis": s["Conc 1 Type (Author)"], "operator": op, "value_raw": raw_value,
             "unit_raw": unit, "value_ug_bee_or_ug_bee_day": float(v), "daily_dose_unit": unit.endswith("/d"),
             "approximation_treated_as_point": op == "~", "record_interval": render(iv),
             "parser_replay": "same cached source row; qualifier omitted upstream",
             "record_pair_is_independent_validation": False}
        r["adult48"] = (r["life_stage"] == "Adult" and r["duration_value"] == "2"
                         and not r["duration_operator"] and r["duration_unit"] == "Day(s)" and not r["daily_dose_unit"])
        r["adult48ai"] = r["adult48"] and r["test_material_dose_basis"] == "Active ingredient"
        for t in THRESHOLDS:
            r[f"benchmark_y_{t}"] = int(c.Tier >= {1: 3, 11: 2, 100: 1}[t])
            r[f"record_operator_ignored_y_{t}"] = int(v <= t)
            r[f"record_qualifier_aware_y_{t}"] = interval_label(iv, t)
        records.append(r)
        groups.setdefault(c.Index, []).append((r, v, iv))
    ledger = pd.DataFrame(records).sort_values(["Index", "source_data_row_0based"])
    assert len(ledger) == 828 and len(groups) == 441
    save(ledger, "01_ecotox_record_ledger.csv")

    candidates = []
    code12 = read("results/label_audit/02_compound_labels_by_rule.csv")
    code12 = code12[(code12["filter"] == "all") & (code12.routes == "fixed")]
    ref = {(int(r.Index), int(r.threshold), r.rule): int(r.label) for r in code12.itertuples(index=False)}
    compared = 0
    for idx, items in sorted(groups.items()):
        rr = [x[0] for x in items]
        r = {k: rr[0][k] for k in ["Index", "CAS", "name", "Tier", "route", "inchikey"]}
        iv = median_interval([x[2] for x in items])
        r.update(n_records=len(items), median_operator_ignored=float(median([x[1] for x in items])),
                 qualifier_median_interval=render(iv), all_adult48=all(x["adult48"] for x in rr),
                 all_adult48ai=all(x["adult48ai"] for x in rr),
                 all_point=all(x["operator"] in ("", "~") for x in rr),
                 contains_exactly_gt100=any(x["operator"] == ">" and x["value_ug_bee_or_ug_bee_day"] == 100 for x in rr),
                 contains_other_right_limit=any(x["operator"] in (">", ">=") and x["value_ug_bee_or_ug_bee_day"] != 100 for x in rr),
                 source_data_rows_0based="|".join(str(x["source_data_row_0based"]) for x in rr),
                 source_row_sha256="|".join(x["source_row_sha256"] for x in rr),
                 ecotox_reference_numbers="|".join(sorted(set(x["ecotox_reference_number"] for x in rr))),
                 source_values="; ".join(f"{x['operator']}{x['value_raw']} {x['unit_raw']}" for x in rr),
                 stages="|".join(sorted(set(x["life_stage"] for x in rr))),
                 durations="|".join(sorted(set(x["duration_operator"] + x["duration_value"] + " " + x["duration_unit"] for x in rr))),
                 material_dose_bases="|".join(sorted(set(x["test_material_dose_basis"] for x in rr))))
        for t in THRESHOLDS:
            r[f"benchmark_y_{t}"] = rr[0][f"benchmark_y_{t}"]
            r[f"median_ignored_y_{t}"] = int(median([x[1] for x in items]) <= t)
            r[f"median_aware_y_{t}"] = interval_label(iv, t)
            r[f"unanimous_aware_y_{t}"] = unanimous([x[f"record_qualifier_aware_y_{t}"] for x in rr])
            if t in (1, 100):
                assert r[f"median_ignored_y_{t}"] == r[f"benchmark_y_{t}"], (idx, t)
            for rule, value in [("D", r[f"median_aware_y_{t}"]), ("C", r[f"unanimous_aware_y_{t}"])]:
                assert value == ref[(idx, t, rule)], (idx, t, rule, value, ref[(idx, t, rule)])
                compared += 1
        candidates.append(r)
    cand = pd.DataFrame(candidates)
    save(cand, "02_all_ecotox_case_candidates.csv")
    strata = [
        ("E01", "Tier1 determinate negative; includes >100", (cand.Tier == 1) & (cand.median_aware_y_100 == 0) & cand.contains_exactly_gt100),
        ("E02", "Tier1 determinate negative; other right limit, no >100", (cand.Tier == 1) & (cand.median_aware_y_100 == 0) & cand.contains_other_right_limit & ~cand.contains_exactly_gt100),
        ("E03", "Tier1 unresolved", (cand.Tier == 1) & (cand.median_aware_y_100 == -1)),
        ("E04", "Tier1 determinate positive; point records", (cand.Tier == 1) & (cand.median_aware_y_100 == 1) & cand.all_point),
        ("E05", "Tier0 determinate negative control", (cand.Tier == 0) & (cand.median_aware_y_100 == 0)),
        ("E06", "Tier3 changes to negative at 1", (cand.Tier == 3) & (cand.median_aware_y_1 == 0)),
        ("E07", "Tier3 unresolved at 1", (cand.Tier == 3) & (cand.median_aware_y_1 == -1)),
        ("E08", "Tier2 originally positive at 11; unanimous-aware unresolved", (cand.Tier == 2) & (cand.unanimous_aware_y_11 == -1)),
    ]
    selected, selection = [], []
    for case_id, desc, mask in strata:
        pool = cand[mask].copy()
        pool["cas_numeric"] = pool.CAS.str.replace("-", "", regex=False).astype(int)
        pool = pool.sort_values(["all_adult48", "n_records", "cas_numeric", "Index"], ascending=[False, True, True, True])
        selection.append({"case_id": case_id, "stratum": desc, "candidate_count": len(pool),
                          "selection_rule": "all adult48 first; fewest records; ascending numeric CAS; benchmark Index"})
        if len(pool):
            selected.append({"case_id": case_id, "stratum": desc, "candidate_count": len(pool), **pool.iloc[0].drop("cas_numeric").to_dict()})
    selected = pd.DataFrame(selected)
    save(pd.DataFrame(selection), "03_selection_rules.csv")
    save(selected, "04_selected_ecotox_cases.csv")
    # Strata may legitimately select the same compound (e.g. unresolved at
    # both 11 and 100). Each case retains every determining source row.
    save(selected[["case_id", "Index"]].merge(ledger, on="Index", validate="many_to_many"), "05_selected_ecotox_records.csv")

    # External record ledger. Preserve missingness rather than infer adult,
    # 48-hour, pure active ingredient or same underlying study from a screen.
    oft = read("output/external_data_audit/oft_acute_endpoints_aligned.csv")
    bee = read("output/external_data_audit/oft_bee_source_rows.csv")
    epa = read("output/external_data_audit/plos_adult_acute_endpoints_aligned.csv")
    cites = read("output/external_data_audit/plos_mrid_citations.csv")
    cite_map = cites.groupby("MRID").Citation.apply(lambda x: " | ".join(sorted(set(x)))).to_dict()
    lit = pd.read_excel(ROOT / "data/external_raw/OFT3.0 export repository.xlsx", sheet_name="LIT", dtype=str).fillna("")
    lit_map = lit.set_index("Document UUID").to_dict("index")
    oftraw = {k: g for k, g in bee.groupby("ResultsAndDiscussion.EffectConcentrations.UUID")}
    ext_records = []
    for source_name, df, fname in [("OpenFoodTox", oft, "oft_acute_endpoints_aligned.csv"),
                                   ("EPA", epa, "plos_adult_acute_endpoints_aligned.csv")]:
        path = "output/external_data_audit/" + fname
        for ordinal, row in df.iterrows():
            c = bench_by_key.get(row.inchikey)
            is_oft = source_name == "OpenFoodTox"
            mrid = "" if is_oft else re.sub(r"\.0$", "", row.MRID).zfill(8) if row.MRID not in ("", "-") else ""
            stage = "|".join(sorted(set(oftraw[row.result_uuid]["ResultsAndDiscussion.EffectConcentrations.LifeStage"]) - {""})) if is_oft else "adult (S1 endpoint category)"
            rs = [] if not is_oft else [lit_map[x] for x in row.literature_uuids.split("|") if x in lit_map]
            r = {"source": source_name, "source_path": path, "source_file_sha256": hashes[path],
                 "source_data_row_0based": ordinal, "source_row_sha256": row_sha(row),
                 "name": row["name"], "CAS": row.CAS, "inchikey": row.inchikey, "route": row.route,
                 "result_id": row.result_uuid if is_oft else row.source_cell,
                 "document_id": row.document_uuid if is_oft else "",
                 "literature_ids": row.literature_uuids if is_oft else mrid,
                 "source_excel_rows_or_cell": row.source_excel_rows if is_oft else row.source_cell,
                 "MRID": mrid, "study_id_status": "OFT document UUID; underlying trial ID not established" if is_oft else "MRID supplied" if mrid else "MRID missing",
                 "citation": " | ".join(sorted(set(x.get("GeneralInfo.Name", "") for x in rs))) if is_oft else cite_map.get(mrid, ""),
                 "source_doi": row.source_dois if is_oft else "10.1371/journal.pone.0265962",
                 "literature_report_numbers": " | ".join(sorted(set(x.get("GeneralInfo.ReportNo", "") for x in rs))) if is_oft else "",
                 "life_stage": stage or "not reported in aligned source fields",
                 "duration_hours": row.duration_hours if is_oft else "",
                 "test_material": row.test_material if is_oft else "",
                 "dose_basis": row.dose_basis if is_oft else "",
                 "dose_basis_other": row.dose_basis_other if is_oft else "",
                 "unit": row.unit if is_oft else "ug/bee",
                 "operator": row.lower_qualifier, "value_ug_bee": row.lower_ug_bee,
                 "upper_operator": row.upper_qualifier if is_oft else "",
                 "upper_value_ug_bee": row.upper_ug_bee if is_oft else "",
                 "material_duration_missingness_scope": "explicit source fields" if is_oft else "not supplied by S1/S4 endpoint/citation tables; not inferred from MRID",
                 "Index": c.Index if c else "", "benchmark_source": c.source if c else "",
                 "benchmark_CAS": c.CAS if c else "", "benchmark_name": c.name if c else "",
                 "benchmark_original_route": c.toxicity_type if c else "",
                 "benchmark_ecotox_data_rows_0based": "|".join(str(x[0]["source_data_row_0based"]) for x in groups.get(c.Index, [])) if c else "",
                 "Tier": c.Tier if c else "", "identity_route_match": bool(c and c.toxicity_type == row.route),
                 "source_screen": row.single_structure_screen == "True" if is_oft else True,
                 "source_48h_screen": row["48h_no_material_flag_screen"] == "True" if is_oft else False,
                 "same_study_as_benchmark": "not established",
                 "study_match_basis": "no shared trial ID supplied; chemical/route alone insufficient"}
            r["main_route_matched_cohort"] = r["identity_route_match"] and r["source_screen"]
            for t in THRESHOLDS:
                r[f"external_y_{t}"] = int(row[f"y_{t}"])
                r[f"benchmark_y_{t}"] = int(c.Tier >= {1: 3, 11: 2, 100: 1}[t]) if c else ""
            ext_records.append(r)
    ext = pd.DataFrame(ext_records)
    save(ext, "06_external_record_ledger.csv")
    cohort = ext[ext.main_route_matched_cohort].copy()
    save(cohort, "07_route_matched_external_records.csv")
    # Every original main-cohort source/key pair must be represented.
    conc = read("output/external_data_audit/same_route_label_concordance.csv")
    main_conc = conc[conc.source.isin(["OpenFoodTox", "PLOS2022"])].copy()
    main_conc["source"] = main_conc.source.replace({"PLOS2022": "EPA"})
    assert set(zip(cohort.source, cohort.inchikey)) == set(zip(main_conc.source, main_conc.inchikey))
    external_cases = []
    for source_name in ("OpenFoodTox", "EPA"):
        for benchmark_source in ("ECOTOX", "PPDB"):
            q = main_conc[(main_conc.source == source_name) & (main_conc.apistox_source == benchmark_source) & (main_conc.agreement_100 == "0")].copy()
            if q.empty:
                continue
            q = q.sort_values(["inchikey", "route"])
            r = q.iloc[0]
            evidence = cohort[(cohort.source == source_name) & (cohort.inchikey == r.inchikey)]
            for _, er in evidence.iterrows():
                external_cases.append({"case_id": f"X{len(external_cases)+1:02d}", "stratum": f"{source_name}/{benchmark_source} discordance at 100",
                                       "candidate_compounds": len(q), "selection_rule": "ascending full InChIKey, then route; include every record of selected compound", **er.to_dict()})
    save(pd.DataFrame(external_cases), "08_selected_cross_source_cases.csv")

    # Identifier groups are not all equivalent to independent experiments.
    inventory = []
    for scope, df in [("all_aligned", ext), ("route_matched_main", cohort)]:
        for src, q in df.groupby("source"):
            for field, meaning in [("MRID", "EPA report/study identifier; multiple compounds or endpoints may share it"),
                                   ("document_id", "OFT document UUID, not proven distinct physical trial"),
                                   ("result_id", "source result UUID or EPA spreadsheet cell, not independent-study proof")]:
                available = q[q[field].ne("")]
                inventory.append({"scope": scope, "source": src, "id_field": field, "records": len(q),
                                  "records_with_id": len(available), "unique_ids": available[field].nunique(),
                                  "repeated_id_records_beyond_first": len(available) - available[field].nunique(), "interpretation": meaning})
    save(pd.DataFrame(inventory), "09_identifier_inventory.csv")
    mrids = []
    for mrid, q in ext[(ext.source == "EPA") & ext.MRID.ne("")].groupby("MRID"):
        mrids.append({"MRID": mrid, "endpoint_records": len(q), "unique_chemicals": q.inchikey.nunique(),
                      "routes": "|".join(sorted(set(q.route))), "source_cells": "|".join(q.result_id),
                      "CAS": "|".join(sorted(set(q.CAS))), "names": " | ".join(sorted(set(q["name"]))),
                      "citation": cite_map.get(mrid, ""), "main_cohort_records": int(q.main_route_matched_cohort.sum())})
    save(pd.DataFrame(mrids), "10_epa_mrid_deduplication.csv")
    missingness = []
    for src, q in cohort.groupby("source"):
        missingness.append({"source": src, "records": len(q), "compounds": q.inchikey.nunique(),
                            "stage_not_reported": int(q.life_stage.eq("not reported in aligned source fields").sum()),
                            "duration_missing": int(q.duration_hours.eq("").sum()),
                            "duration_48h": int(q.duration_hours.eq("48.0").sum() + q.duration_hours.eq("48").sum()),
                            "material_field_missing": int(q.test_material.eq("").sum()),
                            "dose_basis_missing": int(q.dose_basis.eq("").sum()),
                            "dose_basis_test_material": int(q.dose_basis.eq("test mat.").sum()),
                            "dose_basis_active_ingredient": int(q.dose_basis.eq("act. ingr.").sum()),
                            "MRID_missing": int(q.MRID.eq("").sum()),
                            "same_study_as_benchmark_established": 0})
    save(pd.DataFrame(missingness), "11_route_matched_metadata_missingness.csv")
    flow, composition = [], []
    for source_name, src in ext.groupby("source"):
        screens = [(source_name, src[src.source_screen])]
        if source_name == "OpenFoodTox":
            screens.append(("OpenFoodTox_48h_screen", src[src.source_48h_screen]))
        for screen, screened in screens:
            matched = screened[screened.Index.ne("")]
            routed = matched[matched.identity_route_match]
            determined = []
            for key, g in routed.groupby("inchikey"):
                if all(set(g[f"external_y_{t}"]) in ({0}, {1}) for t in THRESHOLDS):
                    determined.append(key)
            final = routed[routed.inchikey.isin(determined)]
            for stage, q in [("source_screen", screened), ("same_full_InChIKey", matched),
                             ("same_benchmark_route", routed), ("all_three_thresholds_determinate", final)]:
                flow.append({"source_screen": screen, "stage": stage, "records": len(q),
                             "nonempty_full_keys": q.loc[q.inchikey.ne(""), "inchikey"].nunique(),
                             "benchmark_compounds": q.loc[q.Index.ne(""), "Index"].nunique(),
                             "scope": "benchmark-wide label concordance; not split-specific saved-score reassessment"})
            selected_indices = set(final.Index)
            for member, pop in [("included", bench[bench.Index.isin(selected_indices)]),
                                ("excluded", bench[~bench.Index.isin(selected_indices)])]:
                for field in ["source", "toxicity_type", "Tier", "herbicide", "fungicide", "insecticide", "other_agrochemical"]:
                    for value, count in pop[field].value_counts(dropna=False).sort_index().items():
                        composition.append({"source_screen": screen, "membership": member,
                                            "n_compounds": len(pop), "field": field,
                                            "value": value, "count": int(count),
                                            "fraction_within_membership": count / len(pop) if len(pop) else ""})
    save(pd.DataFrame(flow), "13_concordance_cohort_flow.csv")
    save(pd.DataFrame(composition), "14_included_excluded_composition.csv")
    duplicate_rows = ledger[ledger.source_row_sha256.duplicated(keep=False)].copy()
    save(duplicate_rows, "15_exact_duplicate_ecotox_rows.csv")
    a = set(main_conc.loc[main_conc.source == "OpenFoodTox", "inchikey"])
    b = set(main_conc.loc[main_conc.source == "EPA", "inchikey"])
    transitions = []
    for tier, q in cand.groupby("Tier"):
        for t in THRESHOLDS:
            transitions.append({"Tier": int(tier), "threshold": t, "compounds": len(q),
                                "median_positive": int((q[f"median_aware_y_{t}"] == 1).sum()),
                                "median_negative": int((q[f"median_aware_y_{t}"] == 0).sum()),
                                "median_unresolved": int((q[f"median_aware_y_{t}"] == -1).sum())})
    save(pd.DataFrame(transitions), "12_exact_interval_transition_counts.csv")
    summary = {"ecotox_determining_records": len(ledger), "ecotox_compounds": len(cand),
               "exact_rational_interval_vs_code12_comparisons": compared, "mismatches": 0,
               "selected_ecotox_cases": len(selected), "source_overlap": {
                   "OpenFoodTox_compounds": len(a), "EPA_compounds": len(b), "both_sources_compounds": len(a & b),
                   "union_compounds": len(a | b), "same_physical_study_established_across_sources": 0,
                   "interpretation": "zero established means unverifiable with supplied trial identifiers, not proven absence of shared studies"},
               "ecotox_export_identifier_limits": "Reference Number only; no test/result ID or MRID columns. References may identify databases, not individual experiments.",
               "OFT_MRID_keyword_hits_in_bee_export": int(bee.apply(lambda c: c.str.contains("MRID", case=False, regex=False)).any(axis=1).sum()),
               "OFT_MRID_keyword_hits_in_literature_export": int(lit.apply(lambda c: c.str.contains("MRID", case=False, regex=False)).any(axis=1).sum()),
               "ecotox_exactly_gt100_records": int(((ledger.operator == ">") & (ledger.value_ug_bee_or_ug_bee_day == 100)).sum()),
               "ecotox_other_right_limit_records": int((ledger.operator.isin([">", ">="]) & (ledger.value_ug_bee_or_ug_bee_day != 100)).sum()),
               "ecotox_exact_duplicate_rows_beyond_first": int(ledger.source_row_sha256.duplicated().sum()),
               "duplicate_handling": "retained for exact replay of upstream aggregation; not asserted to be distinct experiments",
               "metadata_missingness": missingness, "selection": selection,
               "row_hash_encoding": "SHA256 of UTF-8 compact JSON [column,value] pairs in original column order; missing fields empty strings",
               "source_data_row_0based": "parsed source data-record ordinal, header excluded; not a physical text line number",
               "labels": {"1": "LD50 <= cutoff", "0": "LD50 > cutoff", "-1": "unresolved"}}
    (OUT / "RECORD_AUDIT_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (OUT / "INPUT_SHA256.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
    if args.verify_upstream:
        pin = "57c2c05de293a0f54c119744eb3ca9cb03bc030d"
        checks = []
        for local, remote in [("upstream_ecotox.py", "ecotox.py"),
                              ("upstream_ppdb_and_bpdb.py", "ppdb_and_bpdb.py"),
                              ("upstream_processing.py", "processing.py")]:
            url = f"https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/{pin}/dataset_creation/{remote}"
            with urllib.request.urlopen(url, timeout=30) as response:
                remote_bytes = response.read()
            local_path = EXT / "cache" / local
            local_bytes = local_path.read_bytes()
            normalize = lambda b: b.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            checks.append({"commit": pin, "url": url, "cache_path": str(local_path.relative_to(ROOT)).replace("\\", "/"),
                           "cache_sha256": hashlib.sha256(local_bytes).hexdigest(),
                           "remote_sha256": hashlib.sha256(remote_bytes).hexdigest(),
                           "bytes_equal": local_bytes == remote_bytes,
                           "newline_normalized_equal": normalize(local_bytes) == normalize(remote_bytes),
                           "result": "exact match" if local_bytes == remote_bytes else "newline-only difference" if normalize(local_bytes) == normalize(remote_bytes) else "content differs"})
        assert all(x["newline_normalized_equal"] for x in checks), checks
        (OUT / "UPSTREAM_PIN_CHECK.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("selection",)}, indent=2, ensure_ascii=True))
    print(selected[["case_id", "CAS", "name", "source_values", "benchmark_y_1", "median_aware_y_1", "benchmark_y_11", "unanimous_aware_y_11", "benchmark_y_100", "median_aware_y_100"]].to_json(orient="records"))


if __name__ == "__main__":
    main()
