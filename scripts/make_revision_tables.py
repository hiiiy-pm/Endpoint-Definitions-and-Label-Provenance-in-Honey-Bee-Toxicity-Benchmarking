"""Export review-stage supplementary tables from auditable saved results.

No models are fitted. Numeric estimates and intervals are copied from their
stated analysis outputs; quantities calculated here are descriptive counts.
"""
from pathlib import Path
import json
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/latex_tables/supp_tables_revision.tex'
sys.path.insert(0, str(ROOT / 'code'))
from _common import load_study_data  # noqa: E402

lines = []


def esc(x):
    return str(x).replace('&', r'\&').replace('_', r'\_').replace('%', r'\%')


def f(x, digits=3):
    return '--' if pd.isna(x) else f'{x:.{digits}f}'


def table(number, title, spec, head, rows, note=''):
    n = head.count('&') + 1
    lines.extend([r'\Needspace{12\baselineskip}\subsection*{Supplementary Table ' + number + '. ' + title + '}',
                  note, r'\begingroup\footnotesize\setlength{\tabcolsep}{4pt}',
                  r'\begin{longtable}{' + spec + '}',
                  r'\toprule ' + head + r' \\ \midrule\endfirsthead',
                  rf'\multicolumn{{{n}}}{{l}}{{\textit{{Supplementary Table {number} (continued)}}}} \\',
                  r'\toprule ' + head + r' \\ \midrule\endhead',
                  rf'\midrule\multicolumn{{{n}}}{{r}}{{Continued on next page}}\\\endfoot\bottomrule\endlastfoot'])
    lines.extend(' & '.join(map(str, row)) + r' \\' for row in rows)
    lines.append(r'\end{longtable}\endgroup')


# Configuration details are source-code specifications, not estimated results.
config = [
    ['ECFP4', 'Morgan radius 2; 1,024 bits; chirality disabled; Tanimoto kernel.'],
    ['Avalon / MACCS', 'Avalon 1,024 bits; RDKit MACCS 167-bit vector (bit zero unused); Tanimoto kernels.'],
    ['WL-HI', 'Two refinement iterations; element, charge, aromaticity and bond types; normalized histogram intersection (Eq.~\\eqref{eq:si-wl-kernel}).'],
    ['Primary SVM', r'$C=1$; \texttt{class\_weight=balanced}; precomputed kernels; default SVC tolerance $10^{-3}$ and no probability fit.'],
    ['Descriptor SVM', r'12 descriptors; training-fitted \texttt{RobustScaler}; RBF kernel; \texttt{gamma=scale}; $C=1$, balanced class weights.'],
    ['Metadata / fusion', r'L2 logistic, $C=1$, \texttt{lbfgs}, balanced class weights, 5,000 iterations; categorical one-hot encoding. Fusion standardizes its structural score and metadata inputs.'],
    ['Calibration', r'Unweighted logistic sigmoid on five-fold out-of-fold SVM scores; $C=10^6$, \texttt{lbfgs}, 5,000 iterations. Apply that sigmoid to the final full-training model score.'],
    ['Additional logistic', r'ECFP bits; L2, \texttt{lbfgs}, balanced weights, 5,000 iterations; $C\in\{0.01,0.1,1,10\}$ selected by five-fold training AUROC.'],
    ['Random forest', r'500 trees; \texttt{balanced\_subsample}; Gini criterion; bootstrap enabled; \texttt{max\_features=sqrt}; unlimited depth; minimum split/leaf sizes 2/1.'],
    ['Seeds and folds', r'Base seed 20260829; all five-fold partitions are stratified and shuffled. Calibration uses base+$100s+e$, fusion base+$100s+t$, additional-learner CV base+$s$; $s=0,1,2$ for Random/MaxMin/Time, $e=0,1,2$ for thresholds $100,11,1$, and $t$ is the numeric threshold.'],
]
pins = dict(line.strip().split('==',1) for line in (ROOT/'requirements-analysis.txt').read_text().splitlines() if '==' in line)
config.append(['Pinned environment', '; '.join(f'{name} {pins[name]}' for name in ['numpy','pandas','scipy','scikit-learn','rdkit','statsmodels']) + r'. Exact pins are in \path{requirements-analysis.txt}; execution validation is recorded separately.'])
table('S21a', 'Reproducible learner and calibration settings', r'p{0.19\linewidth}p{0.72\linewidth}',
      'Component & Implemented settings', config,
      r'Settings are taken from \path{code/_common.py}, \path{code/01_primary_endpoint_models.py}, \path{code/02_metadata_decomposition.py}, \path{code/06_novelty_calibration.py} and \path{code/08_sensitivity_diagnostics.py}. Fixed learner specifications do not imply identical class weights across endpoints. Environment versions are retained in the repository environment lock and validation records.')

a = pd.read_csv(ROOT / 'results/decomposition/03_tier_composition_association_tests.csv')
names = {'source': 'Source', 'toxicity_type': 'Exposure route', 'insecticide': 'Insecticide',
         'herbicide': 'Herbicide', 'fungicide': 'Fungicide', 'other_agrochemical': 'Other agrochemical'}
rows = [[names.get(r.Variable, esc(r.Variable)), f(r.Chi2), int(r.df), f'{r.P:.3g}', f(r.CramerV), f'{r.BH_q:.3g}'] for r in a.itertuples()]
table('S21b', 'The six tier--metadata association tests', 'lrrrrr',
      r'Variable & $\chi^2$ & df & $p$ & Cram\'er\textquotesingle s $V$ & BH $q$', rows,
      r'All six association tests form one Benjamini--Hochberg family. The four use flags can overlap. These associations do not identify causal effects of use or data source.')

joint = pd.read_csv(ROOT / 'results/applicability/04_joint_interaction_tests.csv')
inter = pd.read_csv(ROOT / 'results/applicability/02_boundary_novelty_interactions.csv')
rows = []
for r in inter.itertuples():
    j = joint[(joint.Split == r.Split) & (joint.Outcome == r.Outcome)].iloc[0]
    rows.append([r.Split, 'Brier' if r.Outcome == 'BrierLoss' else 'Error logit',
                 r.Contrast.split(' novelty')[0], f(r.SlopeDifference), f(r.CI_low) + ' to ' + f(r.CI_high),
                 f(r.HolmP), f(j.JointInteractionChi2), f(j.P)])
table('S21c', 'Formal endpoint-by-novelty interactions', 'lllccccc',
      r'Split & Outcome & Contrast & Coefficient & 95\% CI & Holm $p$ & Joint $\chi^2$ & Joint $p$', rows,
      r'Coefficients and Wald intervals are per one-unit novelty increase; divide by ten to obtain the 0.1-unit scale used for marginal slopes in Fig.~S3c. Classification-error interactions use a logit link. Holm adjustment covers the two interactions within each split and outcome; each joint Wald test has two degrees of freedom. Nonsignificance does not establish absence of an effect.')

study = load_study_data()
pred = pd.read_csv(ROOT / 'results/primary/02_test_predictions.csv')
rows = []
for sp in ['Random', 'MaxMin', 'Time']:
    q = pred[(pred.Split == sp) & (pred.Representation == 'ECFP') & (pred.Threshold == 100)]
    sizes = q.groupby('Scaffold').size()
    rows.append([sp, len(q), len(sizes), int((sizes == 1).sum()), int(sizes.min()), f(sizes.median(), 1), int(sizes.max())])
table('S21d', 'Held-out scaffold-cluster sizes and permutation strata', 'lrrrrrr',
      'Split & Compounds & Clusters & Singletons & Minimum & Median & Maximum', rows,
      r'Bootstrap clusters are the stored Bemis--Murcko scaffold identifiers, with acyclic compounds assigned normalized-structure clusters. Singletons are retained.')
rows = []
df = study.df.assign(Tier=study.tier)
for label, cols in [('Global', []), ('Source + insecticide', ['source', 'insecticide']),
                    ('Source + route + insecticide', ['source', 'toxicity_type', 'insecticide'])]:
    groups = [df] if not cols else [q for _, q in df.groupby(cols)]
    rows.append([label, len(groups), min(len(g) for g in groups), max(len(g) for g in groups),
                 sum(len(g) == 1 for g in groups), sum(g.Tier.nunique() == 1 for g in groups)])
strata_head = r'\toprule Permutation scheme & Strata & Min. size & Max. size & Singleton & Single-tier \\ \midrule'
lines.extend([r'\begin{longtable}{lrrrrr}', strata_head + r'\endfirsthead',
              r'\multicolumn{6}{l}{\textit{Supplementary Table S21d (continued)}} \\',
              strata_head + r'\endhead', r'\bottomrule\endlastfoot'])
lines.extend(' & '.join(map(str, row)) + r' \\' for row in rows)
lines.append(r'\end{longtable}')
lines.append(r'No strata are merged. Singleton strata remain unchanged; other strata are permuted internally (a single-tier stratum has no effective label exchange). Fifteen representation-by-neighbourhood tests are BH-adjusted separately under each of the three null schemes.')

counts = pd.read_csv(ROOT / 'results/label_audit/03_rule_transition_counts.csv', dtype={'Tier': str})


def tally(q, t):
    r = q[q.threshold == t].iloc[0]
    return f'{int(r.positive)}/{int(r.negative)}/{int(r.unresolved)}/{int(r.no_eligible_record)}'


rule_note = (r'A: numerical median, operators ignored; B: numerical unanimous-record rule, operators ignored; '
             r'C: qualifier-aware unanimous-record rule; D: interval median. Each cell is positive/negative/unresolved/no eligible record, '
             r'in that order, at the stated cut-off in $\mu$g/bee. The denominator is fixed at the distributed ECOTOX-derived identities; '
             r'absence after a strict screen is kept distinct from an interval that crosses the threshold. A and D differ only in qualifier handling; '
             r'A and B differ only in aggregation. At 11, rule B uses strict $>11$ for negative records; the original code also accepted '
             r'groups containing exactly 11 and larger values as negative. All compound-level definitions and transitions are in '
             r'\path{results/label_audit/02_compound_labels_by_rule.csv} and \path{results/label_audit/03_rule_transition_counts.csv}.')
rows = []
for tier in ['0', '1', '2', '3', 'all']:
    for rule in ['A', 'B', 'C', 'D']:
        q = counts[(counts['filter'] == 'all') & (counts.routes == 'fixed') & (counts.Tier == tier) & (counts.rule == rule)]
        rows.append(['All' if tier == 'all' else 'Tier' + tier, rule, int(q.n.iloc[0])] + [tally(q,t) for t in [100,11,1]])
table('S22a', 'Qualifier handling and aggregation on the historical determining route', 'llrccc',
      r'Tier & Rule & $n$ & $\leq100$: P/N/U/none & $\leq11$: P/N/U/none & $\leq1$: P/N/U/none', rows, rule_note)

rows = []
for filt in ['all', 'adult48', 'adult48ai']:
    for routes in ['fixed', 'eligible', 'all']:
        for rule in ['A', 'B', 'C', 'D']:
            for tier in ['all', '1']:
                q = counts[(counts['filter'] == filt) & (counts.routes == routes) & (counts.Tier == tier) & (counts.rule == rule)]
                rows.append([{'all':'All records','adult48':'Adult48','adult48ai':'Adult48+AI'}[filt],
                             routes.capitalize(), rule, 'All' if tier == 'all' else 'Tier1'] + [tally(q,t) for t in [100,11,1]])
table('S22b', 'Record-filter and route-selection sensitivity at all three cut-offs', 'llllccc',
      r'Record filter & Routes & Rule & Tier & $\leq100$: P/N/U/none & $\leq11$: P/N/U/none & $\leq1$: P/N/U/none', rows,
      r'Cell order and rules follow S22a. Adult48 retains explicitly adult, exactly two-day observations without per-day dose units. Adult48+AI additionally requires an active-ingredient concentration-type field; it does not independently establish purity. Fixed retains the historical determining route. Eligible combines only routes retained by the upstream numerical 11-cut-off filter; All also restores routes removed by that filter. Records are screened before route-level aggregation. A compound is positive if any eligible route is positive, negative only if every eligible route is negative, and otherwise unresolved. This evaluates threshold labels for the minimum over routes; it does not identify a unique corrected route when interval order is ambiguous. All and Tier1 rows have original denominators 441 and 212, respectively.')

het = pd.read_csv(ROOT / 'results/label_audit/04_determining_group_heterogeneity.csv', dtype={'Tier':str})
rows = [[('All' if r.Tier == 'all' else 'Tier'+r.Tier), r.groups, r.records, r.all_adult,
         r.all_duration_2d, r.any_daily_unit, r.all_adult48, r.any_adult48] for r in het.itertuples()]
table('S22c', 'Historical record-group heterogeneity and upstream selection', 'lrrrrrrr',
      r'Tier & Groups & Records & All adult & All 2 d & Any /d & All Adult48 & Any Adult48', rows,
      r'Counts describe the original label-determining ECOTOX groups before stricter record filtering. These benchmark annotations are not harmonized adult 48-hour toxicity measurements. All Adult48 requires every record to satisfy adulthood, exactly two days and non-daily dose basis simultaneously.')
summary = json.loads((ROOT / 'results/label_audit/LABEL_AUDIT_SUMMARY.json').read_text())
lines.append(f"The 2,674-record cached export yields {summary['records_after_upstream_filters']:,} records after the upstream value, unit and route filters, in {summary['cas_route_groups']} CAS--route groups. The original 11-cut-off rule removes {summary['groups_removed_by_upstream_11_filter']} groups across {summary['cas_with_removed_route']} CAS identities; {summary['cas_with_removed_route_retained_by_another_route']} of those identities have another eligible ECOTOX route. These are export-level counts, not counts of the 441 retained benchmark identities. All 441 original binary labels, ternary levels and determining routes are reproduced.")
lines.append(r'The historical determining groups include 34 numerical threshold-straddlers at 100, three at 11 and 16 at 1. The three 11 cases contain an exact numerical value of 11 and larger values, which the upstream implementation accepts as negative. Qualifier-aware conclusions therefore cannot be inferred from the native binary filter alone.')

gap = pd.read_csv(ROOT / 'results/label_audit/07_crossover_gap.csv')
auc = pd.read_csv(ROOT / 'results/label_audit/06_crossover_auroc.csv')
changes = pd.read_csv(ROOT / 'results/label_audit/08_crossover_label_changes.csv')
rows = []
for sp in ['Random','MaxMin','Time']:
    for rep in ['ECFP','Avalon','MACCS','WL-HI','Descriptors12']:
        for t in [100,11,1]:
            q = auc[(auc.Scenario == 'primary') & (auc.Unresolved == 'keep') & (auc.Split == sp) & (auc.Representation == rep) & (auc.Threshold == t)].set_index(['TrainLabels','EvalLabels'])
            vals = [f(q.loc[(tr,ev),'AUROC']) for tr,ev in [('original','original'),('original','corrected'),('corrected','original'),('corrected','corrected')]]
            rows.append([sp, rep, t] + vals + [f"{int(q.loc[('original','original'),'n_pos'])}/{int(q.loc[('original','corrected'),'n_pos'])}"])
table('S23a', 'Fixed-membership primary label crossover for every representation', 'llrccccc',
      r'Split & Representation & Cut-off & O/O & O/C & C/O & C/C & Positive O/C', rows,
      r'Each cell is AUROC. Before/after the slash denotes training/evaluation labels: O, distributed benchmark; C, partially revised labels. The primary scenario uses qualifier-aware interval medians at 100 and 1, and unanimity at 11, combining the upstream-eligible routes. Unresolved labels retain their distributed values; PPDB/BPDB labels are unchanged. Each split retains its original 828 training and 207 test identities, and no external cohort is introduced. Endpoint class counts and 2,000-draw scaffold percentile intervals for all four cells are in \path{results/label_audit/06_crossover_auroc.csv}. Revised labels are an auditable sensitivity definition, not established ground truth.')

rows = []
for scenario in ['primary','conservative','adult48','regulatory','combined']:
    for how in ['keep','exclude','negative']:
        for sp in ['Random','MaxMin','Time']:
            q = gap[(gap.Scenario == scenario) & (gap.Unresolved == how) & (gap.Split == sp) & (gap.Representation == 'ECFP') & (gap.Severe == 11)].set_index(['TrainLabels','EvalLabels'])
            c = changes[(changes.Scenario == scenario) & (changes.Unresolved == how) & (changes.Split == sp)].iloc[0]
            d = q.loc[('change:corrected-vs-original','change:corrected-vs-original')]
            vals = [f(q.loc[(tr,ev),'Gap']) for tr,ev in [('original','original'),('original','corrected'),('corrected','original'),('corrected','corrected')]]
            rows.append([{'primary':'Primary','conservative':'Cons.','adult48':'Adult48','regulatory':'Reg.','combined':'Comb.'}[scenario],
                         {'keep':'K','exclude':'E','negative':'N'}[how], sp, f'{int(c.n_train)}/{int(c.n_test)}'] + vals + [f'{f(d.Gap)} [{f(d.CI_low)}, {f(d.CI_high)}]'])
table('S23b', 'ECFP 11-minus-100 crossover contrasts and uncertainty across every sensitivity scenario', 'llllccccc',
      r'Scenario & Unresolved & Split & Train/test & O/O & O/C & C/O & C/C & C/C minus O/O [95\% CI]', rows,
      r'Scenario and unresolved-label rules are defined in Supplementary Methods: Cons., conservative; Reg., regulatory; Comb., combined; K, keep; E, exclude; N, negative. O/C notation follows S23a. Both training-label arms use the same retained training and test members, including refitting the original-label arm after exclusion. Intervals are paired scaffold-cluster percentiles from 2,000 draws conditional on fitted models; they are not simultaneous and do not include uncertainty in label adjudication. Negative can create severe-positive/broad-negative labels: it is a nonnested extreme stress test, excluded from substantive correction inference and not a feasible corrected nested benchmark. The complete two-contrast results for all five representations and all four cells, including evaluation-only changes, are in \path{results/label_audit/07_crossover_gap.csv}. Label-change and class counts are in \path{results/label_audit/08_crossover_label_changes.csv}.')
nest = pd.read_csv(ROOT / 'results/label_audit/LABEL_NESTING_CHECK.csv')
negative = nest[nest.Unresolved == 'negative']
assert (nest[nest.Unresolved.isin(['keep','exclude'])].CorrectedLabel_nonnested == 0).all()
lines.append(f'Independent verification of exported test labels found no nesting violations under keep or exclude. All {len(negative)} negative-policy test cohorts contained violations ({int(negative.CorrectedLabel_nonnested.min())}--{int(negative.CorrectedLabel_nonnested.max())} of 207 identities). ' + r'The per-cohort counts and offending indices are in \path{results/label_audit/LABEL_NESTING_CHECK.csv}. This check concerns exported test labels; the prediction artifact does not contain training-label vectors.')

pilot = ROOT / 'output/unseen_identity_case_study'
candidates = pd.read_csv(pilot / 'data_curation/all_candidate_structures.csv')
newpred = pd.concat([pd.read_csv(pilot / f'predictions/{batch}/new_compound_predictions.csv') for batch in ['initial8','added2']], ignore_index=True)
eligible = pd.read_csv(pilot / 'data_curation/frozen_scoreable_48h_endpoints.csv')
reasons = {'inpyrfluxam':'Adult-stage evidence inadequate', 'oxazosulfyl':'Only 72/96-hour observations',
           'dimpropyridaz':'Observation time unverified', 'fluhexafon':'No eligible endpoint',
           'tyclopyrazoflor':'No verified primary endpoint', 'fluchlordiniliprole':'No verified primary endpoint',
           'cyclopyrimorate':'No verified primary endpoint', 'benzpyrimoxan':'Unresolved test-material identity'}
rows = []
for r in candidates.itertuples():
    predicted = r.candidate in set(newpred.compound_id)
    included = r.candidate in set(eligible.candidate)
    overlap = bool(r.overlap_connectivity)
    reason = 'Eligible 48-hour endpoint' if included else ('Benchmark connectivity overlap' if overlap else reasons[r.candidate])
    rows.append([r.candidate.capitalize(), '1' if r.candidate_origin == 'previous8' else '2',
                 'Yes' if overlap else 'No', 'Yes' if predicted else 'No', 'Yes' if included else 'No', reason])
table('S19c', 'Candidate flow, prediction stages and endpoint eligibility', r'lccccp{0.28\linewidth}',
      'Candidate & Batch & Overlap & Predicted & Evaluated & Status', rows,
      r'Sixteen candidates yielded 14 identities absent by connectivity, ten predicted structures and six identities with 11 accepted route endpoints. This is a retrospective, evidence-limited pilot. Batch 1 structures were frozen before prediction; the two added batch 2 structures were selected using endpoint availability. The prediction manifests record 2026-09-23 01:21:07 and 01:22:28 (UTC+8) for the two runs. No endpoint file was read by either prediction runner. This sequencing does not make candidate selection prospective or establish an independent population sample.')
lines.append(r'The structure and prediction hashes are in \path{output/unseen_identity_case_study/predictions/initial8/prediction_run_manifest.json} and \path{output/unseen_identity_case_study/predictions/added2/prediction_run_manifest.json}. The endpoint-file hashes and the curation no-prediction-file check are in \path{output/unseen_identity_case_study/data_curation/endpoint_freeze_manifest.json}; that manifest has no embedded timestamp, so a separate endpoint-freeze time is not claimed.')
rows = []
for (name, sp), q in newpred.groupby(['compound_id','split'], sort=False):
    q = q.set_index('threshold_ug_bee')
    rows.append([name.capitalize(), sp, f(q.nearest_training_tanimoto.iloc[0])] +
                [f'{f(q.loc[t,"score"])} ({int(q.loc[t,"predicted_label"])})' for t in [100,11,1]])
table('S19d', 'Individual new-structure scores, default decisions and nearest-training similarity', 'llcccc',
      r'Identity & Split & Nearest Tanimoto & $\leq100$ score (class) & $\leq11$ score (class) & $\leq1$ score (class)', rows,
      r'Continuous values are uncalibrated ECFP SVM decision scores, not probabilities. Class uses the saved zero decision threshold of the balanced SVC; it was not chosen for a regulatory decision, an external calibration target or an error-cost function. Similarity is the largest training-set ECFP Tanimoto for the corresponding split. Scores are recorded for all ten predicted identities, including the four without eligible outcomes. Source-linked outcomes and correctness at each determinate threshold are in \path{output/unseen_identity_case_study/evaluation/new_compound_case_predictions.csv}; full nearest-neighbour identities and model hashes are in the two prediction CSVs.')

cases = pd.read_csv(ROOT / 'results/record_audit/04_selected_ecotox_cases.csv')
rows = []
for r in cases.itertuples():
    vector = lambda prefix: '/'.join('U' if int(getattr(r,f'{prefix}_y_{t}')) < 0 else str(int(getattr(r,f'{prefix}_y_{t}'))) for t in [100,11,1])
    raw = esc(r.source_values).replace('>', r'$>$')
    rows.append([r.case_id, r.CAS, r.route, raw,
                 'Yes/Yes' if r.all_adult48ai else ('Yes/No' if r.all_adult48 else 'No/No'),
                 vector('benchmark'), vector('median_aware'), esc(r.source_data_rows_0based)])
table('S24a', 'Illustrative same-record ECOTOX operator replays', r'lllp{0.21\linewidth}cccl',
      'Case & CAS & Route & Original value and unit & Adult48/AI & Original & Rule D & Raw rows', rows,
      r'Label vectors are ordered 100/11/1, with U unresolved. Raw rows are zero-based data-row ordinals in the unchanged cached export. Each replay compares the same retained source row(s) with and without its operator; it is not an independent validation study. Selection used documented strata with deterministic tie handling and retains the empty E02 stratum in the selection ledger. Rule C equals D for these displayed examples; broader disagreements between the two rules are retained in S22. For E06, three $>1000$ ng/organism observations convert to $>1$ $\mu$g/bee, but stage is unreported and durations are two, three and four days. AI denotes the reported active-ingredient dose field, not independently verified purity.')
lines.append(r'Full study metadata, original units, purity fields, source-file and source-row hashes, and exact replay values are in \path{results/record_audit/05_selected_ecotox_records.csv}. The cached export does not supply ECOTOX test/result IDs or MRIDs; its reference number identifies a cited source, not necessarily an individual experiment. Cases therefore establish record-level processing differences without claiming matched underlying experiments across databases.')
cross = pd.read_csv(ROOT / 'results/record_audit/08_selected_cross_source_cases.csv')
rows = []
for r in cross.itertuples():
    cas = r.CAS if pd.notna(r.CAS) else r.benchmark_CAS
    rows.append([r.case_id, 'OFT' if r.source == 'OpenFoodTox' else 'EPA', esc(r['name']) if isinstance(r,dict) else esc(getattr(r,'name')),
                 esc(cas), r.route, '$>$'+f(r.value_ug_bee,1), '--' if pd.isna(r.duration_hours) else str(int(r.duration_hours)),
                 '--' if pd.isna(r.MRID) else str(int(r.MRID))])
table('S24b', 'Illustrative regulatory-source disagreements kept separate from parser replays', 'lllllcrc',
      r'Case & Source & Identity & CAS & Route & Value & Hours & MRID', rows,
      r'Values are in $\mu$g/bee. All four identities are benchmark-positive and source-negative at 100 on the matched route. The same underlying experiment as the benchmark is not established for any displayed case; missing duration or material information remains missing. They illustrate source disagreement and are not counted as adjudicated extraction errors. Exact identity keys, study/document fields, material and dose-basis evidence, row hashes and selection reasons are in \path{results/record_audit/08_selected_cross_source_cases.csv}.')

flow = pd.read_csv(ROOT / 'results/record_audit/13_concordance_cohort_flow.csv')
source_labels = {'OpenFoodTox':'OFT','OpenFoodTox_48h_screen':'OFT48','EPA':'EPA'}
stages = ['same_full_InChIKey','same_benchmark_route','all_three_thresholds_determinate']
rows = []
for source in source_labels:
    q = flow[flow.source_screen == source].set_index('stage')
    rows.append([source_labels[source]] + [int(q.loc[stage,'benchmark_compounds']) for stage in stages])
table('S24c', 'Benchmark-wide concordance cohort flow and source-metadata coverage', 'lrrr',
      'Source screen & Full-key matched & Same route & Three cut-offs determinate', rows,
      r'Counts are unique benchmark identities at each sequential stage, not the split-specific saved-score cohorts in S18. Retaining a common three-cut-off cohort introduces selection by source availability, route and threshold determinacy. Source screens and matching exclusions are recorded in \path{results/record_audit/13_concordance_cohort_flow.csv}.')
missing = pd.read_csv(ROOT / 'results/record_audit/11_route_matched_metadata_missingness.csv')
lines.extend([r'\begin{longtable}{lrrrrrr}',
              r'\toprule Source & Records & Time missing & 48 h & Material missing & Dose basis missing & MRID missing \\ \midrule'])
for r in missing.itertuples():
    lines.append(' & '.join(map(str,[source_labels[r.source],r.records,r.duration_missing,r.duration_48h,
                                   r.material_field_missing,r.dose_basis_missing,r.MRID_missing])) + r' \\')
lines.append(r'\bottomrule\end{longtable}')
lines.append(r'OFT has 12 further records at non-48-hour durations; all 263 matched-route records lack explicit life-stage fields, and 11 state an active-ingredient dose basis. EPA records belong to the source\textquotesingle s adult category, but individual study duration, material and dose-basis fields are not supplied in these tables. Its 128 nonmissing MRIDs are distinct. Missing metadata are not evidence that two studies match. The 216 OFT and 133 EPA same-route identities overlap in 71 identities (278 in their union); these are neither disjoint samples nor independent replications.')

composition = pd.read_csv(ROOT / 'results/record_audit/14_included_excluded_composition.csv')
rows = []
for field in ['source','toxicity_type','Tier','insecticide','herbicide','fungicide','other_agrochemical']:
    subset = composition[composition.field == field]
    for value in sorted(subset.value.unique(), key=str):
        cells = []
        for source in source_labels:
            for membership in ['included','excluded']:
                q = subset[(subset.source_screen == source) & (subset.membership == membership) & (subset.value == value)]
                cells.append('0 (0.0)' if q.empty else f'{int(q.iloc[0]["count"])} ({100*q.iloc[0].fraction_within_membership:.1f})')
        rows.append([{'toxicity_type':'Route','other_agrochemical':'Other use'}.get(field,field.capitalize()),esc(value)] + cells)
table('S24d', 'Composition of included and excluded benchmark identities', 'llcccccc',
      'Field & Value & OFT in & OFT out & OFT48 in & OFT48 out & EPA in & EPA out', rows,
      r'Cells give counts and percentages within membership groups. Included means same-route and all-three-determinate; excluded includes every other benchmark identity. Denominators are OFT 149/886, OFT48 72/963 and EPA 103/932 (included/excluded). Use flags overlap. These distributions describe selection and do not establish representativeness; no structural-selection test is inferred from them. Machine-readable counts are in \path{results/record_audit/14_included_excluded_composition.csv}.')

structure = pd.read_csv(ROOT / 'results/record_audit/16_selection_structure_summary.csv')
rows = []
for source in source_labels:
    for membership in ['included','excluded']:
        r = structure[(structure.source_screen == source) & (structure.membership == membership)].iloc[0]
        rows.append([source_labels[source], membership.capitalize(), int(r.n_compounds),
                     int(r.selected_cohort_size), int(r.reference_compounds_per_query),
                     f(r['median']), f(r.q1) + ' to ' + f(r.q3)])
table('S24e', 'Structural proximity to the selected concordance cohorts', 'llrrrrc',
      r'Source & Membership & $n$ & Selected set & References/query & Median & Q1--Q3', rows,
      r'For every benchmark molecule, the descriptive statistic is its largest ECFP4 Tanimoto similarity to the corresponding same-route, all-three-determinate selected cohort. An included molecule is compared with the other selected molecules, excluding its own identity; an excluded molecule is compared with the entire selected set. The table reports medians and interquartile ranges across query molecules, not confidence intervals. The reference-set size is shown because nearest-neighbour similarity depends on that size; the three screens are not directly comparable standardized coverage measures. These summaries quantify proximity to each selected chemical pool and do not prove representativeness or identify a causal selection effect. The saved kernel is reused without model fitting or cohort changes. Per-molecule values and cohort summaries are supplied with \path{results/record_audit/16_selection_structure_summary.csv}.')

OUT.write_text('\n\n'.join(lines) + '\n', encoding='utf-8')
print('Wrote supplementary revision tables S19c-d, S21a-d, S22a-c, S23a-b and S24a-e to', OUT)
