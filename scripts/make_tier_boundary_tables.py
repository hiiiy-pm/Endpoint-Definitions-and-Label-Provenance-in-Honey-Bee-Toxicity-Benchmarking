"""Write Supplementary Tables S11, S14, S16 and S17 from the tier-boundary outputs.

Inputs are the saved CSV files written by code/10_tier_boundary_diagnostics.py
and code/11_qualifier_propagation_audit.py. No model is fitted here.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / 'results' / 'tier_boundary'
OUT = ROOT / 'output' / 'latex_tables' / 'supp_tables_tier_boundary.tex'
SPLITS = ['Random', 'MaxMin', 'Time']
REPS = ['ECFP', 'Avalon', 'MACCS', 'WL-HI', 'Descriptors12']
PAIRS = [(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2)]
lines = []


def add(text):
    lines.append(text)


def f3(x):
    s = f'{x:.3f}'
    return '$-$' + s[1:] if s.startswith('-') else s


def interval(low, high):
    return f'{f3(low)} to {f3(high)}'


def longtable(title, spec, head, rows, note, size=r'\small', sep='4pt'):
    add(r'\Needspace{12\baselineskip}\subsection*{' + title + '}')
    if note:
        add(note)
    add(r'\begingroup' + size + r'\setlength{\tabcolsep}{' + sep + '}')
    add(r'\begin{longtable}{' + spec + '}')
    add(r'\toprule ' + head + r' \\ \midrule\endfirsthead')
    add(r'\toprule ' + head + r' \\ \midrule\endhead')
    add(r'\midrule\multicolumn{' + str(head.count('&') + 1) +
        r'}{r}{Continued on next page}\endfoot\bottomrule\endlastfoot')
    add('\n'.join(' & '.join(r) + r' \\' for r in rows))
    add(r'\end{longtable}\endgroup')


# ---------------- S11: tier-pair AUROC and composition weights ----------------
pairs = pd.read_csv(RES / '01_tier_pair_auroc.csv')
rows = []
for sp in SPLITS:
    for rep in ('ECFP', 'WL-HI'):
        for t in (100, 11, 1):
            q = pairs[(pairs.Split == sp) & (pairs.Representation == rep) & (pairs.TrainedFor == t)]
            cells = []
            for hi, lo in PAIRS:
                r = q[(q.HigherTier == hi) & (q.LowerTier == lo)].iloc[0]
                cells.append(f3(r.AUROC))
            key = q[(q.HigherTier == 1) & (q.LowerTier == 0)].iloc[0]
            rows.append([sp, rep, rf'$\leq{t}$'] + cells + [interval(key.ScaffoldCI_low, key.ScaffoldCI_high)])
longtable('Supplementary Table S11. Tier-pair AUROC of the saved held-out scores',
          'lllccccccc',
          r'Split & Rep. & Model & T1$>$T0 & T2$>$T0 & T3$>$T0 & T2$>$T1 & T3$>$T1 & T3$>$T2 & T1$>$T0 95\% CI',
          rows,
          r'Each cell is the probability that a test compound of the higher tier receives a higher score than one of the '
          r'lower tier, with half credit for ties, computed from the saved scores of the model trained for the stated '
          r'endpoint. The last column gives the scaffold-cluster 95\% percentile interval for the Tier1--Tier0 pair '
          r'(5,000 draws identical to the primary analysis). Intervals for every pair and all five representations are in '
          r'\texttt{results\_v2/tier\_boundary/01\_tier\_pair\_auroc.csv}.',
          size=r'\footnotesize', sep='3pt')
weight_rows = []
ecfp = pairs[(pairs.Representation == 'ECFP') & (pairs.TrainedFor == 100)]
for sp in SPLITS:
    q = ecfp[ecfp.Split == sp]
    for t in (100, 11, 1):
        cells = []
        for hi, lo in PAIRS:
            w = float(q[(q.HigherTier == hi) & (q.LowerTier == lo)].iloc[0][f'Weight_le{t}'])
            cells.append(f3(w) if w > 0 else '--')
        weight_rows.append([sp, rf'$\leq{t}$'] + cells)
add(r'Composition weights $w_{r\ell}^{(b)}=n_rn_\ell/(n_b^+n_b^-)$ of each tier pair in each endpoint AUROC '
    r'(Supplementary Eq.~S2); -- marks pairs that do not enter that endpoint.')
add(r'\begingroup\footnotesize\setlength{\tabcolsep}{3pt}')
add(r'\begin{longtable}{llcccccc}')
add(r'\toprule Split & Endpoint & T1$>$T0 & T2$>$T0 & T3$>$T0 & T2$>$T1 & T3$>$T1 & T3$>$T2 \\ \midrule\endfirsthead')
add(r'\toprule Split & Endpoint & T1$>$T0 & T2$>$T0 & T3$>$T0 & T2$>$T1 & T3$>$T1 & T3$>$T2 \\ \midrule\endhead')
add(r'\bottomrule\endlastfoot')
add('\n'.join(' & '.join(r) + r' \\' for r in weight_rows))
add(r'\end{longtable}\endgroup')

# ---------------- S14: Tier1 exclusion for all representations ----------------
excl = pd.read_csv(RES / '02_tier1_exclusion_auroc.csv')
comp = pd.read_csv(RES / '03_tier1_exclusion_decomposition.csv')
rows = []
for sp in SPLITS:
    for rep in REPS:
        e = excl[(excl.Split == sp) & (excl.Representation == rep) & (excl.Bootstrap == 'ScaffoldCluster')].set_index('Quantity')
        c = comp[(comp.Split == sp) & (comp.Representation == rep) & (comp.Bootstrap == 'ScaffoldCluster')].set_index('Contrast')
        rows.append([sp, rep,
                     f3(e.loc['A100_full', 'AUROC']), f3(e.loc['A100_noT1_eval', 'AUROC']),
                     f3(e.loc['Aclean_noT1_train_eval', 'AUROC']), f3(e.loc['A11_full', 'AUROC']),
                     f3(c.loc['Tier1_total_component', 'Estimate']),
                     interval(c.loc['Tier1_total_component', 'CI_low'], c.loc['Tier1_total_component', 'CI_high']),
                     f3(c.loc['Residual_11_minus_clean', 'Estimate']),
                     interval(c.loc['Residual_11_minus_clean', 'CI_low'], c.loc['Residual_11_minus_clean', 'CI_high'])])
longtable('Supplementary Table S14. Tier1 exclusion for all representations',
          'llcccccccc',
          r'Split & Rep. & $\leq100$ & No-T1 test & No-T1 refit & $\leq11$ & T1 comp. & 95\% CI & Residual & 95\% CI',
          rows,
          r'$\leq100$ and $\leq11$: original models on the full test set. No-T1 test: the $\leq100$ model evaluated after '
          r'removing Tier1 test compounds (Tier0 versus Tier2 and Tier3). No-T1 refit: the same fixed learner refitted '
          r'without Tier1 training compounds and evaluated on the Tier1-free test set. T1 comp.: No-T1 refit minus '
          r'$\leq100$. Residual: $\leq11$ minus No-T1 refit. The two components add to the $\leq11$ minus $\leq100$ '
          r'contrast. Intervals are scaffold-cluster 95\% percentile intervals over the primary bootstrap draws; '
          r'molecule-level intervals and the separate evaluation and training components are in '
          r'\texttt{results\_v2/tier\_boundary/03\_tier1\_exclusion\_decomposition.csv}.',
          size=r'\footnotesize', sep='3pt')

# ---------------- S16: concordance by benchmark tier and record source --------
tier = pd.read_csv(RES / '04_concordance_by_tier.csv')
rows = []
for r in tier.itertuples():
    rows.append([r.Source, f'Tier{r.Tier}', str(r.n), str(r.Agree_100), str(r.Disagree_100),
                 f'{100 * r.Agreement_rate_100:.1f}',
                 f'{100 * r.Agreement_Wilson_low:.1f} to {100 * r.Agreement_Wilson_high:.1f}',
                 f'{r.Disagree_benchmark_pos_source_neg}/{r.Disagree_benchmark_neg_source_pos}',
                 str(r.Disagree_with_right_censored)])
longtable('Supplementary Table S16. Route-matched concordance at 100~\\ugbee by benchmark tier and record source',
          'llrrrrcrr',
          r'Source & Tier & $n$ & Agree & Disagree & Agree (\%) & Wilson 95\% CI (\%) & B+S$-$/B$-$S+ & Censored',
          rows,
          r'All-three-determinate cohorts of Supplementary Table~S15. B+S$-$: benchmark label $\leq100$ and source '
          r'label $>100$; B$-$S+: the reverse. Censored: disagreements whose source evidence includes a right-censored '
          r'$>100$~\ugbee value. OFT: OpenFoodTox; OFT 48 h: 48-hour observations only; EPA: the 2022 retrospective.',
          size=r'\footnotesize', sep='3pt')
src = pd.read_csv(RES / '05_concordance_by_record_source.csv')
add(r'Disagreements at 100~\ugbee by the benchmark record source (the database that supplied the ApisTox label).')
add(r'\begingroup\footnotesize\setlength{\tabcolsep}{4pt}')
add(r'\begin{longtable}{llrrrr}')
add(r'\toprule Source & Benchmark record & $n$ & Disagree & Tier1 $n$ & Tier1 disagree \\ \midrule\endfirsthead')
add(r'\toprule Source & Benchmark record & $n$ & Disagree & Tier1 $n$ & Tier1 disagree \\ \midrule\endhead')
add(r'\bottomrule\endlastfoot')
add('\n'.join(' & '.join([r.Source, r.Benchmark_record_source, str(r.n), str(r.Disagree_100), str(r.n_Tier1),
                          str(r.Disagree_100_Tier1)]) + r' \\' for r in src.itertuples()))
add(r'\end{longtable}\endgroup')

# ---------------- S17: qualifier propagation through the ECOTOX labelling ------
qual = pd.read_csv(RES / '07_ecotox_qualifier_by_tier.csv')
rows = []
for r in qual.itertuples():
    rows.append([f'Tier{r.Tier}', str(r.n_ecotox_compounds), str(r.any_right_censored),
                 str(r.right_censored_ge100), str(r.median_exactly_100),
                 str(r.interval_label_100_positive), str(r.interval_label_100_unresolved),
                 str(r.interval_label_100_negative)])
longtable('Supplementary Table S17. Inequality qualifiers in the label-determining ECOTOX records',
          'lrrrrrrr',
          r'Tier & $n$ & Any $>$ & $>$ at $\geq100$ & Median $=100$ & $\leq100$ & Unresolved & $>100$',
          rows,
          r'The upstream ECOTOX labelling was replicated from the cached ApisTox ECOTOX export and reproduced the label '
          r'and ternary level of all 441 ECOTOX-derived benchmark compounds. For each compound, the label-determining '
          r'record group (the exposure route with the lowest median) was re-read with its operator column. Any $>$: '
          r'group contains at least one right-censored record; $>$ at $\geq100$: a right-censored record at or above '
          r'100~\ugbee. The last three columns give the qualifier-aware label at 100~\ugbee, which requires every record '
          r'of the group to fall determinately on one side of the cut-off. Compound-level rows are in '
          r'\texttt{results\_v2/tier\_boundary/06\_ecotox\_qualifier\_propagation.csv}.',
          size=r'\footnotesize', sep='4pt')

OUT.write_text('\n\n'.join(lines) + '\n', encoding='utf-8')
print('Supplementary Tables S11, S14, S16 and S17 written to', OUT.name)
