"""Generate main-text Tables 3 and 4 from the scientific analysis outputs.

The rows are written to output/main_tables/. External-evaluation results are generated in the supplementary tables.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output' / 'main_tables'
OUT.mkdir(parents=True, exist_ok=True)
SPLITS = ['Random', 'MaxMin', 'Time']


def num(x):
    """Three decimals with a typographic minus sign."""
    s = f'{x:.3f}'
    return '$-$' + s[1:] if s.startswith('-') else s


def table3_rows():
    tb = ROOT / 'results' / 'tier_boundary'
    e = pd.read_csv(tb / '02_tier1_exclusion_auroc.csv')
    c = pd.read_csv(tb / '03_tier1_exclusion_decomposition.csv')
    rows = []
    for sp in SPLITS:
        q = e[(e.Split == sp) & (e.Representation == 'ECFP') & (e.Bootstrap == 'ScaffoldCluster')].set_index('Quantity')
        k = c[(c.Split == sp) & (c.Representation == 'ECFP') & (c.Bootstrap == 'ScaffoldCluster')].set_index('Contrast')
        comp, res = k.loc['Tier1_total_component'], k.loc['Residual_11_minus_clean']
        rows.append(' & '.join([
            sp, num(q.loc['A100_full', 'AUROC']), num(q.loc['A100_noT1_eval', 'AUROC']),
            num(q.loc['Aclean_noT1_train_eval', 'AUROC']), num(q.loc['A11_full', 'AUROC']),
            num(k.loc['Gap_11_minus_100', 'Estimate']),
            f"{num(comp.Estimate)} ({num(comp.CI_low)} to {num(comp.CI_high)})",
            f"{num(res.Estimate)} ({num(res.CI_low)} to {num(res.CI_high)})"]) + r' \\')
    return rows


def table4_rows():
    t = pd.read_csv(ROOT / 'results' / 'tier_boundary' / '04_concordance_by_tier.csv')
    rows = []
    for tier in range(4):
        cells = [f'Tier{tier}']
        for src in ['OFT', 'OFT 48 h', 'EPA']:
            r = t[(t.Source == src) & (t.Tier == tier)].iloc[0]
            cells.append(f'{int(r.Agree_100)}/{int(r.n)} ({100 * r.Agreement_rate_100:.1f})')
        rows.append(' & '.join(cells) + r' \\')
    return rows


if __name__ == '__main__':
    for name, rows in [('table3_tier1_exclusion', table3_rows()),
                       ('table4_concordance_by_tier', table4_rows())]:
        (OUT / f'{name}.tex').write_text('\n'.join(rows) + '\n', encoding='utf-8')
        print(f'== {name}')
        print('\n'.join(rows))
