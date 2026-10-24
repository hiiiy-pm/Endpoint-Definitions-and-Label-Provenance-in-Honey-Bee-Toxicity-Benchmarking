"""Generate Figure 1 from the frozen benchmark and official split records.

Panel a shows the analytical tiers and nested labels, panel b shows the
test-set composition, and panel c connects the three research questions.
The diagram uses the 1,035-compound benchmark without filtering or resampling.
Exports are a 170 x 114 mm PDF, SVG, 600-dpi PNG, and source-data metadata.
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'jcheminf_fig1_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--research-root', type=Path, required=True,
                    help='Research directory containing the supplied data and results.')
parser.add_argument('--output-dir', type=Path, default=HERE / 'figures')
args = parser.parse_args()
src, out = args.research_root, args.output_dir
out.mkdir(parents=True, exist_ok=True)

proc_path = src / 'results/cache/processed_data.csv'
split_dir = src / 'data/official_splits'
composition_path = src / 'source_data/figure_composition_data.csv'
proc = pd.read_csv(proc_path)
composition = pd.read_csv(composition_path).set_index('Tier')
assert proc.SMILES.is_unique, 'Expected the original unique compound set.'
tier_n = [int((proc.Tier == t).sum()) for t in range(4)]
assert tier_n == [177, 562, 125, 171]
n_total = len(proc)
assert n_total == 1035
ins_share = [100 * proc.loc[proc.Tier == t, 'insecticide'].mean() for t in range(4)]
for t in range(4):
    assert np.isclose(ins_share[t] / 100, composition.loc[f'Tier{t}', 'Insecticide'])

splits = ['Random', 'MaxMin', 'Time']
endpoints = ['100', '11', '1']
expected = {'Random': [179, 59, 33], 'MaxMin': [162, 42, 23], 'Time': [167, 41, 24]}
tier_lookup = proc.set_index('SMILES').Tier
test_pos = {}
for split in splits:
    test = pd.read_csv(split_dir / f'{split.lower()}_test.csv')
    train = pd.read_csv(split_dir / f'{split.lower()}_train.csv')
    assert len(train) == 828 and len(test) == 207
    assert set(train.SMILES).isdisjoint(set(test.SMILES))
    assert set(train.SMILES) | set(test.SMILES) == set(proc.SMILES)
    tier = tier_lookup.reindex(test.SMILES)
    assert tier.notna().all()
    counts = [int((tier >= t).sum()) for t in (1, 2, 3)]
    assert counts == expected[split]
    test_pos[split] = dict(zip(endpoints, counts))

ink = '#1A1A1A'
muted = '#656565'
ep_color = {'100': '#595959', '11': '#0072B2', '1': '#D55E00'}
tier_color = ['#E3EAEF', '#A7C0D0', '#F0A868', '#BF5130']
tier_text_color = [ink, ink, ink, 'white']
tier_range = ['> 100', '(11,100]', '(1,11]', '\u2264 1']
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans'],
    'font.size': 7.7,
    'text.color': ink,
    'axes.labelcolor': ink,
    'axes.edgecolor': '#9A9A9A',
    'axes.linewidth': 0.6,
    'xtick.color': ink,
    'ytick.color': ink,
    'xtick.major.size': 2.6,
    'ytick.major.size': 2.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

fig = plt.figure(figsize=(17 / 2.54, 11.4 / 2.54), facecolor='white')

# Panel a: original population, tier intervals, compositions and nested labels.
ax_a = fig.add_axes([0.012, 0.559, 0.976, 0.424])
ax_a.set_xlim(0, 1)
ax_a.set_ylim(0, 1)
ax_a.axis('off')
ax_a.text(0.0, 1.0, 'a', fontsize=10.8, fontweight='bold', va='top')
ax_a.text(0.118, 0.996, 'Nested endpoints from one fixed compound set',
          fontsize=9.4, fontweight='bold', va='top')
ax_a.text(0.118, 0.874, f'{n_total:,} ApisTox compounds ordered by annotation tier',
          fontsize=7.7, color=muted, va='top')

x_start, x_end = 0.146, 0.994
usable = x_end - x_start
edges = np.cumsum([0] + tier_n)
xs = x_start + usable * edges / n_total
label_x = x_start - 0.017
for t in range(4):
    left, right = xs[t], xs[t+1]
    width = right - left
    ax_a.add_patch(Rectangle((left, 0.637), width, 0.115, facecolor='#F2F2F2',
                             edgecolor='white', lw=0.8))
    ax_a.add_patch(Rectangle((left, 0.637), width, 0.115 * ins_share[t] / 100,
                             facecolor=ep_color['1'], edgecolor='white', lw=0.8,
                             alpha=0.85))
    ax_a.text((left + right) / 2, 0.77, f'{ins_share[t]:.0f}%', ha='center',
              fontsize=7.6, fontweight='bold', color=ep_color['1'])
    ax_a.add_patch(Rectangle((left, 0.332), width, 0.254,
                             facecolor=tier_color[t], edgecolor='white', lw=1.0))
    for yy, label, size, weight in [
        (0.528, f'Tier{t}', 8.1, 'bold'),
        (0.450, f'n = {tier_n[t]}', 7.5, 'normal'),
        (0.373, tier_range[t], 7.5, 'normal'),
    ]:
        ax_a.text((left + right) / 2, yy, label, ha='center', va='center',
                                fontsize=size, fontweight=weight, color=tier_text_color[t])
ax_a.text(label_x, 0.697, 'Insecticide\nshare', ha='right', va='center',
          fontsize=7.3, color=muted, linespacing=1.15)
ax_a.text(label_x, 0.459, 'Tier range\n(\u00b5g/bee)', ha='right', va='center',
          fontsize=7.3, color=muted, linespacing=1.15)

positive_counts = {}
for i, (ep, first) in enumerate(zip(endpoints, (1, 2, 3))):
    yy = 0.241 - i * 0.085
    left = xs[first]
    n_positive = n_total - int(edges[first])
    positive_counts[ep] = n_positive
    ax_a.add_patch(Rectangle((left, yy), x_end - left, 0.064,
                             facecolor=ep_color[ep], edgecolor='none', alpha=0.92))
    label = f'{n_positive} positive ({n_positive / n_total * 100:.1f}%)'
    if ep != '1':
        ax_a.text((left + x_end) / 2, yy + 0.032, label, ha='center', va='center',
                  fontsize=7.5, fontweight='bold', color='white')
    else:
        ax_a.text(left - 0.012, yy + 0.032, label, ha='right', va='center',
                  fontsize=7.5, fontweight='bold', color=ep_color[ep])
    ax_a.text(label_x, yy + 0.032, f'\u2264 {ep} \u00b5g/bee', ha='right', va='center',
              fontsize=7.6, fontweight='bold', color=ep_color[ep])
ax_a.text(xs[2] - 0.012, 0.188, 'Native benchmark task', ha='right', va='center',
          fontsize=7.1, color=muted)

# Panel b: the original official split positive counts, not new model results.
fig.text(0.012, 0.522, 'b', fontsize=10.8, fontweight='bold', va='top')
fig.text(0.095, 0.518, 'Official evaluation splits', fontsize=9.1,
         fontweight='bold', va='top')
fig.text(0.095, 0.478, '828 training and 207 test compounds each',
         fontsize=7.4, color=muted, va='top')
ax_b = fig.add_axes([0.095, 0.080, 0.382, 0.361])
x = np.arange(len(splits))
bar_width = 0.255
for i, ep in enumerate(endpoints):
    vals = [test_pos[split][ep] for split in splits]
    positions = x + (i - 1) * bar_width
    ax_b.bar(positions, vals, width=bar_width, color=ep_color[ep],
             edgecolor='white', linewidth=0.5)
    for xx, v in zip(positions, vals):
        ax_b.text(xx, v + 4, str(v), ha='center', va='bottom', fontsize=7.5)
ax_b.set_xticks(x, splits, fontsize=7.9)
ax_b.set_ylim(0, 262)
ax_b.set_yticks([0, 50, 100, 150, 200])
ax_b.tick_params(axis='y', labelsize=7.5)
ax_b.set_ylabel('Positive test compounds', fontsize=7.9, labelpad=5)
ax_b.spines[['top', 'right']].set_visible(False)
ax_b.legend(handles=[Rectangle((0, 0), 1, 1, color=ep_color[e],
                                label=f'\u2264 {e}') for e in endpoints],
             ncol=3, frameon=False, loc='upper left', fontsize=7.3,
             handlelength=1.0, columnspacing=1.1, handletextpad=0.4,
             borderaxespad=0.0)
ax_b.text(1.0, 0.986, '\u00b5g/bee', transform=ax_b.transAxes, ha='right', va='top',
          fontsize=7.2, color=muted)

# Panel c: the three scientific audit questions.
fig.text(0.525, 0.522, 'c', fontsize=10.8, fontweight='bold', va='top')
fig.text(0.562, 0.518, 'Three linked audit questions', fontsize=9.1,
         fontweight='bold', va='top')
ax_c = fig.add_axes([0.553, 0.076, 0.435, 0.395])
ax_c.set_xlim(0, 1)
ax_c.set_ylim(0, 1)
ax_c.axis('off')
questions = [
    ('1', 'Localize endpoint differences',
     'Paired performance; tier pairs;\ncomposition controls', ep_color['11']),
    ('2', 'Trace label construction',
     'Upstream replay; qualifiers;\naggregation', '#00866A'),
    ('3', 'Test revision consequences',
     'Matched training/evaluation\nlabel crossover', ep_color['1']),
]
for i, (number, title, body, color) in enumerate(questions):
    top = 0.99 - i * 0.300
    bottom = top - 0.266
    ax_c.add_patch(FancyBboxPatch((0.0, bottom), 1.0, 0.266,
                                 boxstyle='round,pad=0,rounding_size=0.025',
                                 facecolor=color, edgecolor='none', alpha=0.075))
    ax_c.scatter([0.057], [top - 0.078], s=103, color=color, edgecolors='none', zorder=3)
    ax_c.text(0.057, top - 0.078, number, ha='center', va='center', color='white',
              fontsize=7.4, fontweight='bold', zorder=4)
    ax_c.text(0.115, top - 0.077, title, va='center', fontsize=8.0,
              fontweight='bold', color=color)
    ax_c.text(0.115, top - 0.173, body, va='center', fontsize=7.6,
              color=ink, linespacing=1.20)
ax_c.text(0.02, 0.004, 'Robustness analyses',
          fontsize=7.1, color=muted, va='bottom')

stem = out / 'Figure_1'
fig.savefig(stem.with_suffix('.pdf'), metadata={'Title': 'Study design: three linked audit questions'})
fig.savefig(stem.with_suffix('.svg'))
fig.savefig(stem.with_suffix('.png'), dpi=600)
plt.close(fig)

manifest = {
    'width_mm': 170, 'height_mm': 114,
    'source_files': [str(proc_path), str(composition_path)] +
                    [str(split_dir / f'{s.lower()}_test.csv') for s in splits],
    'compound_count': n_total,
    'tier_counts': tier_n,
    'insecticide_percent_exact': ins_share,
    'insecticide_percent_display': [f'{x:.0f}' for x in ins_share],
    'endpoint_positive_counts': positive_counts,
    'endpoint_positive_percent_display': {e: f'{v / n_total * 100:.1f}' for e, v in positive_counts.items()},
    'train_count_each': 828, 'test_count_each': 207,
    'test_positive_counts': test_pos,
    'new_experimental_data': False,
    'schematic_scope': 'Three audit questions with observed population, label, composition and split counts.',
    'matplotlib_version': matplotlib.__version__,
    'pandas_version': pd.__version__,
    'numpy_version': np.__version__,
}
(out / 'Figure_1_source_data.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps(manifest, indent=2))
print(f'Wrote {stem.with_suffix(".pdf")}')
