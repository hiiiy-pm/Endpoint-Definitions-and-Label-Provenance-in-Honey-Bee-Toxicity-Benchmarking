"""Generate all manuscript and supplementary figures.

Figures are sized for the IEEE Access two-column layout: double-column figures are
7.16 in wide, single-column figures 3.5 in. All numbers are read from the frozen
machine-readable source data in source_data/ and results/; nothing is hard-coded.
"""
from __future__ import annotations

from pathlib import Path
import textwrap
import argparse
import json
try:  # optional panel-alignment auditor (see output/figure_qa/ for the recorded reports)
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    require_matplotlib_panel_alignment = None

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
DATA = ROOT / 'source_data'
RES = ROOT / 'results'
OUT.mkdir(parents=True, exist_ok=True)
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--only', nargs='+', help='Export only the requested figure stems')
_SELECTED = _parser.parse_args().only
_QA = ROOT / 'output/figure_qa'
_QA.mkdir(parents=True, exist_ok=True)
MAIN_STEMS = ('Fig1_', 'Fig2_', 'Fig3_', 'Fig4_', 'Fig5_', 'Fig6_')

# IEEE Access column geometry (inches).
W2 = 7.16   # double column
W1 = 3.50   # single column

# Okabe-Ito colourblind-safe palette.
C_EP = {'100': '#595959', '11': '#0072B2', '1': '#D55E00'}
C_REP = {'ECFP': '#0072B2', 'Avalon': '#009E73', 'MACCS': '#E69F00',
         'WL-HI': '#CC79A7', 'Descriptors12': '#595959'}
C_SCHEME = {'Global': '#595959', 'Source×Insecticide': '#0072B2',
            'Source×Route×Insecticide': '#009E73'}
INK = '#1A1A1A'
MUTED = '#6B6B6B'
GRID = '#E4E4E4'

SPLITS = ['Random', 'MaxMin', 'Time']
SPLIT_SUB = {'Random': 'random holdout', 'MaxMin': 'structure-diverse holdout',
             'Time': 'publication-year proxy'}
REPS = ['ECFP', 'Avalon', 'MACCS', 'WL-HI', 'Descriptors12']
EPS = ['100', '11', '1']
# The 11 ug/bee marker is drawn smaller and on top so it stays visible when its
# value coincides with the 1 ug/bee marker.
S_EP = {'100': 30, '11': 15, '1': 30}
Z_EP = {'100': 3, '11': 5, '1': 4}
M_EP = {'100': 4.9, '11': 3.5, '1': 4.9}

matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8.0,
    'axes.titlesize': 8.8,
    'axes.labelsize': 8.0,
    'xtick.labelsize': 7.6,
    'ytick.labelsize': 7.6,
    'legend.fontsize': 7.6,
    'axes.linewidth': 0.6,
    'axes.edgecolor': '#9A9A9A',
    'axes.labelcolor': INK,
    'text.color': INK,
    'xtick.color': MUTED,
    'ytick.color': MUTED,
    'xtick.labelcolor': INK,
    'ytick.labelcolor': INK,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.6,
    'ytick.major.size': 2.6,
    'lines.solid_capstyle': 'round',
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
    'ps.fonttype': 42,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    # PDF output resamples any image at this dpi; the default of 100 blurs raster panels.
    'savefig.dpi': 600,
    'image.interpolation': 'none',
})


def panel(ax, letter, x=-0.14, y=1.14):
    ax.text(x, y, letter, transform=ax.transAxes, ha='left', va='top',
            fontsize=10.6, fontweight='bold', color=INK)


def title(ax, text, sub=None, y=1.06):
    ax.set_title(text, loc='left', pad=9 if sub else 5, fontweight='bold', color=INK)
    if sub:
        ax.text(0, y - 0.075, sub, transform=ax.transAxes, ha='left', va='top',
                fontsize=7.6, color=MUTED)


def header(ax, main, sub=None, y_main=1.20, y_sub=1.085):
    """Two-line panel header drawn in axes coordinates to avoid title collisions."""
    ax.text(0, y_main, main, transform=ax.transAxes, ha='left', va='baseline',
            fontsize=9.0, fontweight='bold', color=INK)
    if sub:
        ax.text(0, y_sub, sub, transform=ax.transAxes, ha='left', va='baseline',
                fontsize=7.4, color=MUTED)


def trim(ax, keep_left=True):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    if not keep_left:
        ax.spines['left'].set_visible(False)


def xgrid(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='x', color=GRID, lw=0.55)


def save(fig, stem, align=None, tick_gap_axis=None):
    """Export one figure, after the panel-alignment check when the auditor is available.

    ``align`` holds keyword options for the alignment gate of a multi-panel
    figure; ``None`` marks a single-panel or schematic figure.
    """
    if _SELECTED and stem not in _SELECTED:
        plt.close(fig)
        return
    if stem.startswith(('Fig1_', 'Fig3_', 'Fig4_', 'Fig6_')):
        from matplotlib.text import Text
        fig.canvas.draw()  # Materialize automatically generated tick labels.
        for axis in fig.axes:
            axis.tick_params(axis='both', labelsize=9.2)
        for label in fig.findobj(Text):
            if label.get_text():
                label.set_fontsize(max(label.get_fontsize(), 9.2))
    fig.canvas.draw()
    if align is not None and require_matplotlib_panel_alignment is not None:
        require_matplotlib_panel_alignment(
            fig, json_out=str(_QA / f'{stem}.alignment.json'),
            overlay_svg=str(_QA / f'{stem}.alignment.svg'), tolerance_pt=1.5,
            gutter_tolerance_pt=1.5, strict=True, **align)
    if tick_gap_axis is not None:
        labels = tick_gap_axis.get_xticklabels()
        boxes = [x.get_window_extent(fig.canvas.get_renderer()) for x in labels]
        gaps = [(boxes[i + 1].x0 - boxes[i].x1) * 72 / fig.dpi for i in range(len(boxes) - 1)]
        assert min(gaps) >= 4, gaps
        (_QA / f'{stem}_tick_spacing.json').write_text(json.dumps({'gaps_pt': gaps,
            'labels': [x.get_text() for x in labels], 'minimum_required_pt': 4}, indent=2),
            encoding='utf-8')
    export = {'bbox_inches': matplotlib.transforms.Bbox.from_bounds(0, 0, *fig.get_size_inches())} if stem.startswith(('Fig1_', 'Fig3_', 'Fig4_', 'Fig6_')) else {}
    fig.savefig(OUT / f'{stem}.pdf', **export)
    fig.savefig(OUT / f'{stem}.png', dpi=600, **export)
    if stem.startswith(MAIN_STEMS):
        fig.savefig(OUT / f'{stem}.svg', **export)
        fig.savefig(OUT / f'{stem}.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'}, **export)
    plt.close(fig)
    print('  wrote', stem)


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
perf = pd.read_csv(DATA / 'figure_primary_data.csv', dtype={'endpoint': str})
# Model decomposition is read from the frozen analysis output so that every
# designated baseline, including the source-and-route model, is shown.
_mdl = pd.read_csv(RES / 'decomposition' / '04_model_decomposition_performance.csv')
_MDL_NAME = {'OriginRoute': 'Source + route',
             'AgrochemicalFlags': 'Agrochemical flags',
             'AllMetadata': 'All metadata',
             'ECFPStructure': 'ECFP structure',
             'ECFP+AllMetadata': 'Structure + metadata'}
model = pd.DataFrame({
    'split': _mdl.Split,
    'model': _mdl.Model.map(_MDL_NAME),
    'endpoint': _mdl.Threshold.astype(str),
    'auroc': _mdl.AUROC,
}).dropna(subset=['model'])
curve = pd.read_csv(DATA / 'figure_learning_curve_raw.csv', dtype={'endpoint': str})
scaf = pd.read_csv(DATA / 'figure_scaffold_corrected.csv')
nov = pd.read_csv(DATA / 'figure_novelty_raw.csv', dtype={'endpoint': str})
rboot = pd.read_csv(RES / 'primary' / '05_representation_paired_bootstrap.csv')
perm = pd.read_csv(RES / 'geometry' / '02_conditional_neighborhood_permutation.csv')
comp = pd.read_csv(RES / 'decomposition' / '02_tier_composition_wide.csv')

rboot = rboot[(rboot.PrimaryComparison) & (rboot.Bootstrap == 'ScaffoldCluster')]

COMP_COLS = [
    ('source:PPDB', 'PPDB', 'Source database'),
    ('source:ECOTOX', 'ECOTOX', 'Source database'),
    ('source:BPDB', 'BPDB', 'Source database'),
    ('toxicity_type:Contact', 'Contact', 'Exposure route'),
    ('toxicity_type:Oral', 'Oral', 'Exposure route'),
    ('toxicity_type:Other', 'Other', 'Exposure route'),
    ('agrochemical_flag:herbicide', 'Herbicide', 'Agrochemical class'),
    ('agrochemical_flag:fungicide', 'Fungicide', 'Agrochemical class'),
    ('agrochemical_flag:insecticide', 'Insecticide', 'Agrochemical class'),
    ('agrochemical_flag:other_agrochemical', 'Other', 'Agrochemical class'),
]
TIER_N = [177, 562, 125, 171]
TIER_RANGE = ['> 100', '(11, 100]', '(1, 11]', '≤ 1']

print('Generating figures for IEEE Access two-column layout')

# --------------------------------------------------------------------------
# Figure 1 - study design
# --------------------------------------------------------------------------
proc = pd.read_csv(RES / 'cache' / 'processed_data.csv')
TIER_COL = ['#E3EAEF', '#A7C0D0', '#F0A868', '#BF5130']
TIER_TXT = [INK, INK, INK, 'white']
N_TOT = len(proc)
edges = np.cumsum([0] + TIER_N)
ins_share = [proc.loc[proc.Tier == t, 'insecticide'].mean() * 100 for t in range(4)]

test_pos = {}
_tier_by_smiles = proc.drop_duplicates('SMILES').set_index('SMILES').Tier
for sp in SPLITS:
    te = pd.read_csv(ROOT / 'data' / 'official_splits' / f'{sp.lower()}_test.csv')
    tier = _tier_by_smiles.reindex(te.SMILES).dropna().values
    test_pos[sp] = {'100': int((tier >= 1).sum()), '11': int((tier >= 2).sum()),
                    '1': int((tier >= 3).sum()), 'n': len(te)}

fig = plt.figure(figsize=(W2, 6.65))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.50], width_ratios=[0.94, 1.06],
                      hspace=0.42, wspace=0.27, left=0.13, right=0.98,
                      bottom=0.07, top=0.90)

# a: The horizontal scale remains proportional to the original tier counts.
ax = fig.add_subplot(gs[0, :])
ax.set_xlim(-14, N_TOT + 14)
ax.set_ylim(-2.58, 1.30)
ax.axis('off')
panel(ax, 'a', x=-0.116, y=1.26)
header(ax, f'Nested labels on {N_TOT:,} compounds',
       'Annotation ranges in µg/bee; tier widths follow compound counts',
       y_main=1.20, y_sub=1.08)
for t in range(4):
    w = TIER_N[t]
    ax.add_patch(plt.Rectangle((edges[t], 0.30), w, 0.62, facecolor='#F2F2F2',
                               edgecolor='white', lw=1.0))
    ax.add_patch(plt.Rectangle((edges[t], 0.30), w, 0.62 * ins_share[t] / 100,
                               facecolor=C_EP['1'], edgecolor='white', lw=1.0, alpha=0.85))
    ax.text(edges[t]+w/2, 1.00, f'{ins_share[t]:.0f}%', ha='center', va='bottom',
            color=C_EP['1'], fontweight='bold')
ax.text(-20, 0.61, 'Insecticide\nshare', ha='right', va='center', color=MUTED)
for t in range(4):
    w = TIER_N[t]
    ax.add_patch(plt.Rectangle((edges[t], -0.86), w, 0.96, facecolor=TIER_COL[t],
                               edgecolor='white', lw=1.2))
    for yy, txt, bold in [(-0.13, f'Tier{t}', True),
                          (-0.44, f'n = {TIER_N[t]}', False),
                          (-0.73, TIER_RANGE[t], False)]:
        ax.text(edges[t]+w/2, yy, txt, ha='center', va='center',
                fontweight='bold' if bold else 'normal', color=TIER_TXT[t])
ax.text(-20, -0.38, 'Tier', ha='right', va='center', color=MUTED)
for i, (ep, first) in enumerate([('100', 1), ('11', 2), ('1', 3)]):
    y0=-1.40-i*0.48; x0=edges[first]; npos=N_TOT-x0
    ax.add_patch(plt.Rectangle((x0,y0), npos, 0.34, facecolor=C_EP[ep], edgecolor='none', alpha=0.9))
    txt=f'{npos} positive ({npos/N_TOT*100:.1f}%)'
    if ep=='100':
        ax.text(x0+npos/2,y0+0.17,txt,ha='center',va='center',color='white',fontweight='bold')
    else:
        ax.text(x0-18,y0+0.17,txt,ha='right',va='center',color=C_EP[ep],fontweight='bold')
    ax.text(-20,y0+0.17,f'≤{ep}',ha='right',va='center',color=C_EP[ep],fontweight='bold')

# b: Official splits, unchanged endpoint counts.
ax = fig.add_subplot(gs[1,0])
panel(ax,'b',x=-0.28,y=1.16)
header(ax,'Official evaluation splits','828 train / 207 test per split',y_main=1.10,y_sub=1.03)
xs=np.arange(len(SPLITS)); wd=0.26
for i,ep in enumerate(EPS):
    vals=[test_pos[sp][ep] for sp in SPLITS]
    ax.bar(xs+(i-1)*wd,vals,width=wd,color=C_EP[ep],edgecolor='white',lw=0.5)
    for x0,v in zip(xs+(i-1)*wd,vals):
        ax.text(x0,v+4,str(v),ha='center',va='bottom')
ax.set_xticks(xs,SPLITS); ax.set_ylabel('Positive test compounds')
ax.set_ylim(0,270);ax.set_yticks([0,50,100,150,200])
ax.legend(ncol=3,handles=[plt.Rectangle((0,0),1,1,color=C_EP[e],label=f'≤{e}') for e in EPS],
          frameon=False,loc='upper center',handlelength=0.9,handletextpad=0.3,columnspacing=0.8)
trim(ax)

# c: Two-line details retain the five original analysis stages.
ax=fig.add_subplot(gs[1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
panel(ax,'c',x=-0.10,y=1.16)
header(ax,'Audit workflow','Paired baseline; diagnostic subsets',y_main=1.10,y_sub=1.03)
steps=[('1','Endpoint contrasts','Five representations; paired\nscaffold bootstrap',C_EP['11']),
       ('2','Composition and controls','Metadata; matched classes;\nsource holdouts','#B8860B'),
       ('3','Tier-boundary localization','Tier-pair AUROC;\nTier1 exclusion','#009E73'),
       ('4','Label provenance','Record rules; matched-label\ncorrection tests',C_EP['1']),
       ('5','Unseen identities','Structure-only predictions;\ndefault decisions','#6B6B6B')]
for i,(num,head_t,body_t,col) in enumerate(steps):
    y0=0.985-i*0.199
    ax.add_patch(FancyBboxPatch((0.005,y0-0.185),0.99,0.185,
        boxstyle='round,pad=0,rounding_size=0.015',facecolor=col,edgecolor='none',alpha=0.08))
    ax.scatter(0.055,y0-0.070,s=110,color=col,edgecolor='none',zorder=3)
    ax.text(0.055,y0-0.070,num,ha='center',va='center',color='white',fontweight='bold',zorder=4)
    ax.text(0.115,y0-0.037,head_t,ha='left',va='center',fontweight='bold',color=col)
    ax.text(0.115,y0-0.116,body_t,ha='left',va='center',color=INK,linespacing=1.2)

save(fig, 'Fig1_study_design', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Figure 2 - endpoint discriminability
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(W2, 5.9))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.18], hspace=0.78, wspace=0.16)

for j, sp in enumerate(SPLITS):
    ax = fig.add_subplot(gs[0, j])
    if j == 0:
        panel(ax, 'a', x=-0.30, y=1.42)
    d = perf[perf.split == sp]
    ys = np.arange(len(REPS))[::-1]
    for y0, rep in zip(ys, REPS):
        r = d[d.representation == rep].set_index('endpoint').reindex(EPS)
        v = r.auroc.values
        ax.plot([v[0], max(v[1], v[2])], [y0, y0], color='#C8C8C8', lw=2.2,
                solid_capstyle='round', zorder=1)
        for ep, val in zip(EPS, v):
            ax.scatter(val, y0, s=S_EP[ep], color=C_EP[ep], zorder=Z_EP[ep],
                       edgecolor='white', linewidth=0.6)
    ax.set_yticks(ys, REPS if j == 0 else [''] * len(REPS))
    ax.set_xlim(0.35, 1.0)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(-0.7, len(REPS) - 0.3)
    ax.set_xlabel('Test AUROC')
    xgrid(ax)
    trim(ax, keep_left=False)
    ax.tick_params(axis='y', length=0, labelcolor=INK)
    header(ax, sp, SPLIT_SUB[sp])

h = [Line2D([0], [0], marker='o', lw=0, color=C_EP[e], markersize=M_EP[e],
            markeredgecolor='white', markeredgewidth=0.6, label=f'≤ {e} µg/bee')
     for e in EPS]
fig.legend(handles=h, frameon=False, ncol=3, loc='upper center',
           bbox_to_anchor=(0.5, 1.02), handletextpad=0.35, columnspacing=1.6)

for j, sp in enumerate(SPLITS):
    ax = fig.add_subplot(gs[1, j])
    if j == 0:
        panel(ax, 'b', x=-0.30, y=1.42)
    d = rboot[rboot.Split == sp]
    ys = np.arange(len(REPS))[::-1]
    off = {'11-100': 0.19, '1-100': -0.19}
    for y0, rep in zip(ys, REPS):
        for comp_key, col in (('11-100', C_EP['11']), ('1-100', C_EP['1'])):
            r = d[(d.Representation == rep) & (d.Comparison == comp_key)]
            if r.empty:
                continue
            r = r.iloc[0]
            yy = y0 + off[comp_key]
            sig = r.SimultaneousCI_low > 0
            ax.plot([r.SimultaneousCI_low, r.SimultaneousCI_high], [yy, yy],
                    color=col, lw=1.15, alpha=0.95 if sig else 0.5, zorder=2)
            ax.scatter(r.ObservedDeltaAUROC, yy, s=22, zorder=3,
                       color=col if sig else 'white',
                       edgecolor=col, linewidth=0.9)
    ax.axvline(0, color='#9A9A9A', lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks(ys, REPS if j == 0 else [''] * len(REPS))
    ax.set_xlim(-0.16, 0.46)
    ax.set_xticks([0.0, 0.2, 0.4])
    ax.set_ylim(-0.7, len(REPS) - 0.3)
    ax.set_xlabel('Δ AUROC vs ≤ 100 µg/bee')
    xgrid(ax)
    trim(ax, keep_left=False)
    ax.tick_params(axis='y', length=0, labelcolor=INK)
    n_sig = int((d.SimultaneousCI_low > 0).sum())
    header(ax, sp, f'{n_sig} of {len(d)} intervals above zero')

h2 = [
    Line2D([0], [0], marker='o', lw=0, color=C_EP['11'], markersize=4.4, label='≤ 11 − ≤ 100'),
    Line2D([0], [0], marker='o', lw=0, color=C_EP['1'], markersize=4.4, label='≤ 1 − ≤ 100'),
    Line2D([0], [0], marker='o', lw=1.1, color='#6B6B6B', markersize=4.4,
           markerfacecolor='#6B6B6B', label='95% interval excludes zero'),
    Line2D([0], [0], marker='o', lw=1.1, color='#6B6B6B', markersize=4.4,
           markerfacecolor='white', alpha=0.6, label='interval includes zero'),
]
fig.legend(handles=h2, frameon=False, ncol=4, loc='upper center',
           bbox_to_anchor=(0.5, 0.545), handletextpad=0.35, columnspacing=1.3)

save(fig, 'Fig2_endpoint_discriminability', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Figure 3 - composition and model decomposition
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(W2, 5.65))
gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1.0], hspace=0.83, wspace=0.18,
                         left=0.23, right=0.93, top=0.86, bottom=0.11)

ax = fig.add_subplot(gs[0, :])
panel(ax, 'a', x=-0.28, y=1.23)
arr = np.array([[comp.loc[t, k] for k, _, _ in COMP_COLS] for t in range(4)])
im = ax.imshow(arr, aspect='auto', vmin=0, vmax=0.75, cmap='Blues')
ax.set_xticks(np.arange(len(COMP_COLS)), [{'Herbicide':'Herb.', 'Fungicide':'Fung.', 'Insecticide':'Insect.'}.get(lab,lab)
               for _, lab, _ in COMP_COLS])
ax.set_yticks(np.arange(4),
              [f'Tier{t}  (n = {TIER_N[t]})' for t in range(4)])
for i in range(arr.shape[0]):
    for j in range(arr.shape[1]):
        ax.text(j, i, f'{arr[i, j] * 100:.1f}', ha='center', va='center',
                fontsize=7.4, color='white' if arr[i, j] > 0.45 else INK)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
# group separators and headers
bounds, start = [], 0
for gi in range(1, len(COMP_COLS)):
    if COMP_COLS[gi][2] != COMP_COLS[gi - 1][2]:
        bounds.append((start, gi - 1, COMP_COLS[gi - 1][2]))
        start = gi
bounds.append((start, len(COMP_COLS) - 1, COMP_COLS[-1][2]))
for a0, b0, name in bounds:
    if a0 > 0:
        ax.axvline(a0 - 0.5, color='white', lw=2.4)
    ax.text((a0 + b0) / 2, -0.85, name, ha='center', va='bottom',
            fontsize=7.7, fontweight='bold', color=MUTED)
ax.set_ylim(3.5, -1.0)
ax.patch.set_visible(False)  # Group headings occupy the white band above the matrix.
# highlight the insecticide gradient
ins = [lab for _, lab, _ in COMP_COLS].index('Insecticide')
ax.add_patch(plt.Rectangle((ins - 0.5, -0.5), 1, 4, fill=False,
                           edgecolor=C_EP['1'], lw=1.5, zorder=5))
cb = fig.colorbar(im, ax=ax, fraction=0.020, pad=0.012)
cb.set_label('Within-tier proportion (%)', fontsize=7.7)
cb.set_ticks([0, 0.25, 0.50, 0.75])
cb.set_ticklabels(['0', '25', '50', '75'])
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(length=2, labelsize=6.6)
ax.set_title('Tier composition',
             loc='left', pad=16, fontweight='bold', color=INK)

MODELS = ['Source + route', 'Agrochemical flags', 'All metadata',
          'ECFP structure', 'Structure + metadata']
avail = [m for m in MODELS if m in set(model.model)]
for j, sp in enumerate(SPLITS):
    ax = fig.add_subplot(gs[1, j])
    if j == 0:
        panel(ax, 'b', x=-0.96, y=1.22)
    d = model[model.split == sp]
    ys = np.arange(len(avail))[::-1]
    for y0, m in zip(ys, avail):
        r = d[d.model == m].set_index('endpoint').reindex(EPS)
        v = r.auroc.values
        finite = v[np.isfinite(v)]
        if len(finite) > 1:
            ax.plot([finite.min(), finite.max()], [y0, y0], color='#C8C8C8',
                    lw=2.2, solid_capstyle='round', zorder=1)
        for ep, val in zip(EPS, v):
            if np.isfinite(val):
                ax.scatter(val, y0, s=S_EP[ep], color=C_EP[ep], zorder=Z_EP[ep],
                           edgecolor='white', linewidth=0.6)
    ax.set_yticks(ys, avail if j == 0 else [''] * len(avail))
    ax.set_xlim(0.40, 1.0)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(-0.7, len(avail) - 0.3)
    ax.set_xlabel('Test AUROC')
    xgrid(ax)
    trim(ax, keep_left=False)
    ax.tick_params(axis='y', length=0)
    ax.set_title(sp, loc='left', fontweight='bold', pad=6, color=INK)

fig.legend(handles=h, frameon=False, ncol=3, loc='upper center',
           bbox_to_anchor=(0.59, 0.49), handletextpad=0.35, columnspacing=1.6)
save(fig, 'Fig3_composition_and_models', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Figure 4 - class-count controls and tier-boundary localization
# --------------------------------------------------------------------------
tpair = pd.read_csv(DATA / 'figure_tier_pair_auroc.csv')
t1x = pd.read_csv(DATA / 'figure_tier1_exclusion.csv')
PAIRS = [(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2)]
# Tier pairs entering each training endpoint's own AUROC (Eq. S2).
OWN_PAIRS = {100: {(1, 0), (2, 0), (3, 0)},
             11: {(2, 0), (3, 0), (2, 1), (3, 1)},
             1: {(3, 0), (3, 1), (3, 2)}}
# Diverging map centred on chance ranking (0.5); vermillion marks inverted ranking.
AUC_DIVERGING = matplotlib.colors.LinearSegmentedColormap.from_list(
    'auc_diverging', ['#B3440E', '#F2D3BE', '#F7F7F7', '#C6DBEC', '#0B5C91'])

fig = plt.figure(figsize=(W2, 8.1))
outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.40, 1.0], hspace=0.74,
                            left=0.23, right=0.95, top=0.90, bottom=0.15)
gsa = outer[0].subgridspec(1, 3, wspace=0.16)
# Panels b and c use the complete width; only panel a has three comparable axes.

for j, sp in enumerate(SPLITS):
    ax = fig.add_subplot(gsa[0, j])
    if j == 0:
        panel(ax, 'a', x=-0.95, y=1.30)
    d = curve[curve.split == sp]
    for ep in EPS:
        q = d[d.endpoint == ep].groupby('n').auroc.quantile([.025, .5, .975]).unstack()
        ax.fill_between(q.index, q[.025], q[.975], color=C_EP[ep], alpha=0.10, lw=0)
        ax.plot(q.index, q[.5], color=C_EP[ep], lw=1.4, marker='o', ms=2.8,
                markeredgecolor='white', markeredgewidth=0.5, zorder=3)
    ax.set_xticks([25, 50, 75, 100, 125])
    ax.set_xlim(18, 132)
    ax.set_ylim(0.42, 1.0)
    ax.set_yticks([0.5, 0.7, 0.9])
    if j == 0:
        ax.set_ylabel('Test AUROC')
    else:
        ax.set_yticklabels([])
    ax.set_xlabel('Compounds per class')
    ax.set_axisbelow(True)
    ax.grid(color=GRID, lw=0.55)
    trim(ax)
    ax.set_title(sp, loc='left', fontweight='bold', pad=6, color=INK)
fig.text(0.59, 0.985, 'Class-count-matched learning curves', ha='center',
         va='top', fontsize=9.0, fontweight='bold', color=INK)
hc = [Line2D([0], [0], color=C_EP[e], lw=1.4, marker='o', ms=3.2,
             markeredgecolor='white', markeredgewidth=0.5, label=f'≤ {e} µg/bee')
      for e in EPS]
fig.legend(handles=hc, frameon=False, ncol=3, loc='upper center',
           bbox_to_anchor=(0.59, 0.968), handletextpad=0.35, columnspacing=1.8,
           handlelength=1.6)

# b: tier-pair ranking probabilities of the saved ECFP scores
ax_b = fig.add_subplot(outer[1])
panel(ax_b, 'b', x=-0.27, y=1.15)
rows_b = [(sp, t) for sp in SPLITS for t in (100, 11, 1)]
M = np.array([[float(tpair[(tpair.Split == sp) & (tpair.TrainedFor == t) &
                           (tpair.HigherTier == hi) & (tpair.LowerTier == lo)].AUROC.iloc[0])
               for hi, lo in PAIRS] for sp, t in rows_b])
im = ax_b.imshow(M, aspect='auto', cmap=AUC_DIVERGING,
                 norm=matplotlib.colors.TwoSlopeNorm(vmin=0.3, vcenter=0.5, vmax=1.0))
for i, (sp, t) in enumerate(rows_b):
    for j, (hi, lo) in enumerate(PAIRS):
        v = M[i, j]
        ax_b.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=9.2,
                  color='white' if (v > 0.86 or v < 0.36) else INK)
        if (hi, lo) in OWN_PAIRS[t]:
            ax_b.add_patch(plt.Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False,
                                         edgecolor='#3A3A3A', lw=0.75, zorder=4))
for k in (3, 6):
    ax_b.axhline(k - 0.5, color='white', lw=2.6)
ax_b.set_xticks(range(len(PAIRS)), [f'T{hi}>T{lo}' for hi, lo in PAIRS], fontsize=7.0)
ax_b.get_xticklabels()[0].set_color(C_EP['1'])
ax_b.get_xticklabels()[0].set_fontweight('bold')
ax_b.set_yticks(range(len(rows_b)), [f'{sp} · ≤{t}' for sp, t in rows_b], fontsize=7.0)
for s in ax_b.spines.values():
    s.set_visible(False)
ax_b.tick_params(length=0)
ax_b.set_title('Tier-pair AUROC of saved ECFP scores', loc='left', fontweight='bold',
               pad=7, color=INK, fontsize=8.5)
cb = fig.colorbar(im, ax=ax_b, fraction=0.045, pad=0.025)
cb.set_ticks([0.3, 0.5, 0.75, 1.0])
cb.set_label('Pairwise AUROC', fontsize=9.2)
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(length=2, labelsize=6.6)

# c: Tier1 exclusion from evaluation, then from training
ax_c = fig.add_subplot(outer[2])
panel(ax_c, 'c', x=-0.27, y=1.20)
QUANT = [('A100_full', '≤100, all tiers', C_EP['100'], 'o', True),
         ('A100_noT1_eval', '≤100, Tier1 removed from test', C_EP['100'], 'o', False),
         ('Aclean_noT1_train_eval', 'Refit and tested without Tier1', '#009E73', 'D', True),
         ('A11_full', '≤11 reference', C_EP['11'], 's', True)]
GAP = 5.4
ticks_c = []
for j, sp in enumerate(SPLITS):
    base = (len(SPLITS) - 1 - j) * GAP
    for k, (q, _lab, col, mk, filled) in enumerate(QUANT):
        y0 = base + (len(QUANT) - 1 - k)
        r = t1x[(t1x.Split == sp) & (t1x.Quantity == q)].iloc[0]
        ax_c.plot([r.CI_low, r.CI_high], [y0, y0], color=col, lw=1.0, alpha=0.9, zorder=2)
        ax_c.scatter(r.AUROC, y0, s=24, marker=mk, color=col if filled else 'white',
                     edgecolor=col, linewidth=0.9, zorder=3)
    ticks_c.append(base + (len(QUANT) - 1) / 2)
ax_c.set_yticks(ticks_c, SPLITS)
ax_c.tick_params(axis='y', length=0)
ax_c.set_ylim(-0.9, (len(SPLITS) - 1) * GAP + len(QUANT) - 0.1)
ax_c.set_xlim(0.40, 1.0)
ax_c.set_xticks([0.4, 0.6, 0.8, 1.0])
ax_c.set_xlabel('Test AUROC (scaffold-cluster 95% interval)')
xgrid(ax_c)
trim(ax_c, keep_left=False)
ax_c.set_title('Tier1 exclusion', loc='left', fontweight='bold', pad=7, color=INK,
               fontsize=8.5)
hq = [Line2D([0], [0], marker=mk, lw=0, color=col, markersize=4.4,
             markerfacecolor=col if filled else 'white', markeredgecolor=col,
             markeredgewidth=0.9, label=lab) for _q, lab, col, mk, filled in QUANT]
ax_c.legend(handles=hq, frameon=False, loc='upper left', bbox_to_anchor=(-0.22, -0.39),
            ncol=2, fontsize=6.6, handletextpad=0.3, columnspacing=1.0)

# Align panel letters in figure coordinates after the colorbar fixes plot widths.
for axis in fig.axes:
    for label in axis.texts:
        if label.get_text() in {'a', 'b', 'c'}:
            label.set_x((0.025 - axis.get_position().x0) / axis.get_position().width)

save(fig, 'Fig4_class_count_and_tier_boundary', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Supplementary Figure S3 - neighbourhood, scaffold and novelty diagnostics
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(W2, 3.30))
gss = fig.add_gridspec(1, 3, wspace=0.95)

# a: significance matrix for Tier2 excess local mixing
ax = fig.add_subplot(gss[0, 0])
panel(ax, 'a', x=-0.52, y=1.26)
schemes = ['Global', 'Source×Insecticide', 'Source×Route×Insecticide']
rows = [(rep, k) for rep in REPS for k in (5, 10, 20)]
M = np.zeros((len(rows), len(schemes)))
for i, (rep, k) in enumerate(rows):
    for j, sc in enumerate(schemes):
        r = perm[(perm.Representation == rep) & (perm.k == k) &
                 (perm.PermutationScheme == sc)]
        M[i, j] = float(r.Significant_q05.iloc[0]) if not r.empty else np.nan
ax.imshow(M, aspect='auto', cmap=matplotlib.colors.ListedColormap(['#ECECEC', '#0072B2']),
          vmin=0, vmax=1)
# cell separators so the grid reads as 15 x 3 rather than solid blocks
for i in range(len(rows) + 1):
    ax.axhline(i - 0.5, color='white', lw=0.9)
for j in range(4):
    ax.axvline(j - 0.5, color='white', lw=0.9)
ax.set_xticks(range(3), ['Global', 'Source ×\ninsect.', '+ route'], fontsize=6.8,
              rotation=0, ha='center')
# Shift the outside labels away from the two-line middle label.
ax.get_xticklabels()[0].set_ha('right')
ax.get_xticklabels()[2].set_ha('left')
ax.set_yticks(range(len(rows)), [f'{rep} · k={k}' for rep, k in rows], fontsize=7.0)
for j in range(3):
    ax.text(j, -1.4, f'{int(np.nansum(M[:, j]))}/15', ha='center', va='bottom',
            fontsize=7.8, fontweight='bold', color=INK)
for s in ax.spines.values():
    s.set_visible(False)
# The column totals sit above the matrix; drop the white axes fill behind them.
ax.patch.set_visible(False)
ax.tick_params(length=0)
ax.set_ylim(len(rows) - 0.5, -2.0)
ax.set_title('Tier2 excess local mixing', loc='left', fontweight='bold',
             pad=16, color=INK, fontsize=8.5)
ax.legend(handles=[
    Line2D([0], [0], marker='s', lw=0, markersize=5, color='#0072B2',
           label='significant at BH-FDR 0.05'),
    Line2D([0], [0], marker='s', lw=0, markersize=5, color='#ECECEC',
           markeredgecolor='#C4C4C4', markeredgewidth=0.5, label='not significant')],
    frameon=False, loc='upper left', bbox_to_anchor=(-0.02, -0.14), fontsize=7.0,
    handletextpad=0.4, ncol=1)
ax_mix = ax

# b: corrected scaffold contrasts
ax = fig.add_subplot(gss[0, 1])
panel(ax, 'b', x=-0.42, y=1.26)
rows_s = list(scaf.itertuples(index=False))
ys = np.arange(len(rows_s))[::-1]
for y0, r in zip(ys, rows_s):
    col = C_EP['11'] if r.contrast == 'T2−T1' else C_EP['1']
    sig = r.high < 0
    ax.plot([r.low, r.high], [y0, y0], color=col, lw=1.15, alpha=0.95 if sig else 0.5)
    ax.scatter(r.estimate, y0, s=20, color=col if sig else 'white',
               edgecolor=col, linewidth=0.9, zorder=3,
               marker='o' if r.analysis == 'Scaffold-equal' else 's')
ax.axvline(0, color='#9A9A9A', lw=0.7, ls=(0, (3, 2)))
ax.set_yticks(ys, [f'{r.minimum} · {"equal" if r.analysis=="Scaffold-equal" else "paired"}'
                   f' · {r.contrast}' for r in rows_s], fontsize=7.0)
ax.set_xlabel('Δ leave-one-out scaffold entropy')
xgrid(ax)
trim(ax, keep_left=False)
ax.tick_params(axis='y', length=0)
ax.set_title('Scaffold context', loc='left', fontweight='bold',
             pad=14, color=INK, fontsize=8.5)

# c: novelty slopes
ax = fig.add_subplot(gss[0, 2])
panel(ax, 'c', x=-0.42, y=1.26)
order = [(sp, ep) for sp in SPLITS for ep in EPS]
ys = np.arange(len(order))[::-1]
for y0, (sp, ep) in zip(ys, order):
    r = nov[(nov.split == sp) & (nov.endpoint == ep)]
    if r.empty:
        continue
    r = r.iloc[0]
    sig = r.low > 0
    ax.plot([r.low, r.high], [y0, y0], color=C_EP[ep], lw=1.15,
            alpha=0.95 if sig else 0.5)
    ax.scatter(r.estimate, y0, s=20, zorder=3,
               color=C_EP[ep] if sig else 'white', edgecolor=C_EP[ep], linewidth=0.9)
ax.axvline(0, color='#9A9A9A', lw=0.7, ls=(0, (3, 2)))
ax.set_yticks(ys, [f'{sp} · ≤{ep}' for sp, ep in order], fontsize=7.0)
ax.set_xlabel('Brier-loss slope per +0.1 novelty')
ax.set_xticks([0.0, 0.02, 0.04])
xgrid(ax)
trim(ax, keep_left=False)
ax.tick_params(axis='y', length=0)
ax.set_title('Novelty and calibrated error', loc='left', fontweight='bold',
             pad=14, color=INK, fontsize=8.5)

save(fig, 'FigS3_structural_diagnostics', align={'require_panel_labels': False},
     tick_gap_axis=ax_mix)

# --------------------------------------------------------------------------
# Figure 5 - label provenance at the 100 ug/bee cut-off
# --------------------------------------------------------------------------
conc = pd.read_csv(DATA / 'figure_concordance_by_tier.csv')
rules = pd.read_csv(RES / 'label_audit/03_rule_transition_counts.csv')
qual = rules[(rules['filter'] == 'all') & (rules.routes == 'fixed')
             & (rules.threshold == 100) & (rules.Tier.astype(str) == '1')].copy()
assert set(qual.rule) == {'A', 'B', 'C', 'D'} and (qual.n == 212).all()
assert (qual[['positive', 'unresolved', 'negative']].sum(axis=1) == qual.n).all()
qual.to_csv(DATA / 'figure_label_rule_comparison.csv', index=False)
SOURCES = [('OFT', 'OpenFoodTox', '#0072B2', 'o'),
           ('OFT 48 h', 'OpenFoodTox, 48 h', '#56B4E9', 's'),
           ('EPA', 'EPA review', '#E69F00', 'D')]
fig = plt.figure(figsize=(W2, 3.05))
gs5 = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.06], wspace=0.66)

ax = fig.add_subplot(gs5[0, 0])
panel(ax, 'a', x=-0.20, y=1.22)
ax.axvspan(0.5, 1.5, color='#FBE3D2', alpha=0.6, lw=0, zorder=0)
offset = {'OFT': -0.22, 'OFT 48 h': 0.0, 'EPA': 0.22}
for key, _lab, col, mk in SOURCES:
    d = conc[conc.Source == key].set_index('Tier')
    for t in range(4):
        r = d.loc[t]
        x = t + offset[key]
        ax.plot([x, x], [100 * r.Agreement_Wilson_low, 100 * r.Agreement_Wilson_high],
                color=col, lw=1.0, zorder=2)
        ax.scatter(x, 100 * r.Agreement_rate_100, s=26, marker=mk, color=col,
                   edgecolor='white', linewidth=0.5, zorder=3)
n_by = {(k, t): int(conc[(conc.Source == k) & (conc.Tier == t)].n.iloc[0])
        for k, *_ in SOURCES for t in range(4)}
ax.set_xticks(range(4), [f'Tier{t}\n' + '/'.join(str(n_by[(k, t)]) for k, *_ in SOURCES)
                         for t in range(4)], fontsize=6.8)
ax.set_xlim(-0.6, 3.6)
ax.set_ylim(0, 106)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel('Agreement at ≤100 µg/bee (%)')
ax.set_axisbelow(True)
ax.grid(axis='y', color=GRID, lw=0.55)
trim(ax)
header(ax, 'Agreement with regulatory records',
       'Wilson 95% intervals; n per source under each tier', y_main=1.13, y_sub=1.035)
ax.legend(handles=[Line2D([0], [0], marker=mk, lw=0, color=col, markersize=4.6,
                          markeredgecolor='white', markeredgewidth=0.5, label=lab)
                   for _k, lab, col, mk in SOURCES],
          frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.04), fontsize=6.8,
          handletextpad=0.3)

ax = fig.add_subplot(gs5[0, 1])
panel(ax, 'b', x=-0.30, y=1.22)
CATS = [('positive', 'Positive', '#4D4D4D', 'white'),
        ('unresolved', 'Unresolved', '#CFCFCF', INK),
        ('negative', 'Negative', C_EP['1'], 'white')]
ys = np.arange(4)[::-1]
for y0, rule in zip(ys, ['A', 'B', 'C', 'D']):
    r = qual[qual.rule == rule].iloc[0]
    n = int(r.n)
    left = 0.0
    for col_name, _lab, col, txt in CATS:
        v = 100 * r[col_name] / n
        ax.barh(y0, v, left=left, height=0.62, color=col, edgecolor='white', lw=0.6)
        if v >= 6:
            ax.text(left + v / 2, y0, str(int(r[col_name])), ha='center', va='center',
                    fontsize=6.6, color=txt)
        left += v
ax.set_yticks(ys, ['A: Numeric median', 'B: Numeric consensus',
                  'C: Interval consensus', 'D: Interval median'], fontsize=7.0)
ax.set_xlim(0, 100)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel('Tier1 compounds (%)')
ax.tick_params(axis='y', length=0)
trim(ax, keep_left=False)
header(ax, 'Qualifier and aggregation effects',
       'Same 212 Tier1 groups; counts in bars', y_main=1.13, y_sub=1.035)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col, label=lab)
                   for _c, lab, col, _t in CATS],
          frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.45, -0.24),
          fontsize=6.8, handlelength=1.0, handletextpad=0.4, columnspacing=1.2)

save(fig, 'Fig5_label_provenance', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Supplementary Figure S1 - full AUROC matrix
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(W2, 2.6), sharey=True)
for ax, sp in zip(axes, SPLITS):
    d = perf[perf.split == sp]
    A = np.array([[float(d[(d.representation == r) & (d.endpoint == e)].auroc.iloc[0])
                   for e in EPS] for r in REPS])
    im = ax.imshow(A, aspect='auto', vmin=0.35, vmax=1.0, cmap='Blues')
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            ax.text(j, i, f'{A[i, j]:.3f}', ha='center', va='center', fontsize=7.2,
                    color='white' if A[i, j] > 0.75 else INK)
    ax.set_xticks(range(3), [f'≤{e}' for e in EPS])
    ax.set_yticks(range(len(REPS)), REPS)
    ax.set_title(sp, loc='left', fontweight='bold', pad=5, color=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axes, fraction=0.020, pad=0.015)
cb.set_label('Test AUROC', fontsize=7.7)
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(length=2, labelsize=6.6)
save(fig, 'FigS1_full_performance', align={'require_panel_labels': False})

# --------------------------------------------------------------------------
# Supplementary Figure S2 - model decomposition matrix
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(W2, 2.6), sharey=True)
for ax, sp in zip(axes, SPLITS):
    d = model[model.split == sp]
    A = np.full((len(avail), 3), np.nan)
    for i, m in enumerate(avail):
        for j, e in enumerate(EPS):
            r = d[(d.model == m) & (d.endpoint == e)]
            if not r.empty:
                A[i, j] = float(r.auroc.iloc[0])
    im = ax.imshow(A, aspect='auto', vmin=0.45, vmax=1.0, cmap='Blues')
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if np.isfinite(A[i, j]):
                ax.text(j, i, f'{A[i, j]:.3f}', ha='center', va='center', fontsize=7.2,
                        color='white' if A[i, j] > 0.80 else INK)
    ax.set_xticks(range(3), [f'≤{e}' for e in EPS])
    ax.set_yticks(range(len(avail)), avail)
    ax.set_title(sp, loc='left', fontweight='bold', pad=5, color=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axes, fraction=0.020, pad=0.015)
cb.set_label('Test AUROC', fontsize=7.7)
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(length=2, labelsize=6.6)
save(fig, 'FigS2_model_decomposition', align={'require_panel_labels': False})


# --------------------------------------------------------------------------
# Figure 6 - chemical structures behind the endpoint definitions
# --------------------------------------------------------------------------
try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
    from PIL import Image
    import io as _io
    RDLogger.DisableLog('rdApp.*')
    _HAVE_RDKIT = True
except ImportError:
    _HAVE_RDKIT = False
    print('  skipped Fig5 (RDKit or Pillow not installed)')

if _HAVE_RDKIT:
    mol_df = proc.copy()
    mol_df['scaffold'] = pd.read_csv(RES / 'cache' / 'scaffolds.csv')['Scaffold'].values
    _gen = GetMorganGenerator(radius=2, fpSize=1024)
    _fps, _keep = [], []
    for _i, _smi in enumerate(mol_df.SMILES):
        _m = Chem.MolFromSmiles(_smi)
        if _m is None:
            continue
        _fps.append(_gen.GetFingerprint(_m))
        _keep.append(_i)
    mol_df = mol_df.iloc[_keep].reset_index(drop=True)

    # Database records carry systematic names; these are their common names.
    _COMMON = {
        "N'-(3,4-Dichlorophenyl)-N-methoxy-N-methylurea": 'Linuron',
        '2-[[[[(4-Chloro-6-methoxy-2-pyrimidinyl)amino]carbonyl]amino]sulfonyl]'
        'benzoic acid, Ethyl ester': 'Chlorimuron-ethyl',
        'alpha-[(Dimethoxyphosphinothioyl)thio] benzeneacetic acid, Ethyl ester': 'Phenthoate',
        'Phosphorothioic acid, O,O-Diethyl-O-(4-nitrophenyl)ester': 'Parathion',
        'Phosphoric acid, Dimethyl 4-nitrophenyl ester': 'Methyl paraoxon',
        '2,6-Dinitro-N,N-dipropyl-4-(trifluoromethyl)benzenamine': 'Trifluralin',
        'N2-Ethyl-N4-(1-methylethyl)-6-(methylthio)-1,3,5-triazine-2,4-diamine': 'Ametryn',
    }

    def _short(name, n=26):
        s = str(name).strip()
        if s in _COMMON:
            return _COMMON[s]
        s = s.split('(Ref')[0].strip().rstrip(',')
        return s if len(s) <= n else s[:n - 1].rstrip() + '\u2026'

    _SS = 4  # supersampling factor for crisp embedded structures

    def _draw(smi, display_width_pt, w=430, h=300):
        """Size atom labels and strokes for the final physical panel width.

        Bond length is a preferred maximum: RDKit can reduce it to fit a complex
        molecule, so the figure does not claim a common physical bond length.
        """
        m = Chem.MolFromSmiles(smi)
        d2d = rdMolDraw2D.MolDraw2DCairo(w * _SS, h * _SS)
        o = d2d.drawOptions()
        pixels_per_point = w * _SS / display_width_pt
        o.bondLineWidth = 0.60 * pixels_per_point
        o.scaleBondWidth = False
        o.fixedFontSize = int(round(6.2 * pixels_per_point))
        o.minFontSize = o.fixedFontSize
        o.maxFontSize = o.fixedFontSize
        o.fixedBondLength = 8.5 * pixels_per_point
        o.padding = 0.05
        o.clearBackground = True
        o.setBackgroundColour((1.0, 1.0, 1.0, 1.0))
        rdMolDraw2D.PrepareAndDrawMolecule(d2d, m)
        d2d.FinishDrawing()
        rgb = np.array(Image.open(_io.BytesIO(d2d.GetDrawingText())).convert('RGB'))
        # Make the white canvas transparent: PDF viewers colour-convert embedded images,
        # which renders an opaque white canvas one level off the page white.
        alpha = np.where((rgb == 255).all(axis=2), 0, 255).astype(np.uint8)
        return np.dstack([rgb, alpha])

    # representative compound per tier: the ECFP medoid
    medoids = []
    for _t in range(4):
        _idx = list(mol_df.index[mol_df.Tier == _t])
        _sub = [_fps[i] for i in _idx]
        _best, _bi = -1.0, _idx[0]
        for _j, _f in enumerate(_sub):
            _mean = (sum(DataStructs.BulkTanimotoSimilarity(_f, _sub)) - 1.0) / (len(_sub) - 1)
            if _mean > _best:
                _best, _bi = _mean, _idx[_j]
        medoids.append(_bi)

    # most similar cross-tier pairs that share a Bemis-Murcko scaffold
    _cand = []
    for _sc, _g in mol_df.groupby('scaffold'):
        if not _sc or len(_g) < 2:
            continue
        _lo, _hi = _g[_g.Tier <= 1], _g[_g.Tier == 3]
        if _lo.empty or _hi.empty:
            continue
        for _a in _lo.index:
            for _b in _hi.index:
                _cand.append((DataStructs.TanimotoSimilarity(_fps[_a], _fps[_b]), _a, _b))
    _cand.sort(key=lambda r: -r[0])
    _used, pairs = set(), []
    for _s, _a, _b in _cand:
        if _a in _used or _b in _used:
            continue
        pairs.append((_s, _a, _b))
        _used.update({_a, _b})
        if len(pairs) == 3:
            break

    fig = plt.figure(figsize=(W2, 5.1))
    gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 0.92], hspace=1.12, wspace=0.05,
                          top=0.82, bottom=0.12, left=0.025, right=0.975)
    fig.text(0.5, 0.965, 'Moving the cut-off reassigns the same compound',
             ha='center', va='baseline', fontsize=8.6, fontweight='bold', color=INK)

    # ---- a: one representative compound per tier, with its class at each endpoint
    for _k, _bi in enumerate(medoids):
        ax = fig.add_subplot(gs[0, _k * 3:_k * 3 + 3])
        ax.imshow(_draw(mol_df.loc[_bi, 'SMILES'],
                        ax.get_position().width * fig.get_figwidth() * 72), interpolation='none')
        ax.axis('off')
        if _k == 0:
            panel(ax, 'a', x=-0.02, y=1.63)
        _t = int(mol_df.loc[_bi, 'Tier'])
        ax.text(0.5, 1.20, f'Tier{_t}', transform=ax.transAxes, ha='center', va='baseline',
                fontsize=8.2, fontweight='bold', color=TIER_TXT[_t],
                bbox=dict(boxstyle='round,pad=0.32', facecolor=TIER_COL[_t], edgecolor='none'))
        ax.text(0.5, 1.06, _short(mol_df.loc[_bi, 'name']), transform=ax.transAxes,
                ha='center', va='baseline', fontsize=6.6, color=MUTED)
        for _j, _ep in enumerate(EPS):
            _pos = _t >= {'100': 1, '11': 2, '1': 3}[_ep]
            _x = 0.5 + (_j - 1) * 0.27
            ax.scatter(_x, -0.06, s=58, transform=ax.transAxes, clip_on=False, zorder=5,
                       color=C_EP[_ep] if _pos else 'white', edgecolor=C_EP[_ep], linewidth=1.1)
            ax.text(_x, -0.17, f'\u2264{_ep}', transform=ax.transAxes, ha='center',
                    va='top', fontsize=6.8, color=INK)
    fig.text(0.5, 0.485, 'filled marker, compound is positive at that endpoint;  '
                         'open marker, negative',
             ha='center', va='bottom', fontsize=7.2, color=MUTED)
    fig.text(0.5, 0.415, 'Similar structures carry opposite database labels '
                         'at the severe endpoints',
             ha='center', va='baseline', fontsize=8.6, fontweight='bold', color=INK)

    # ---- b: near-identical structures separated by the endpoint definition
    _pair_axes = []
    for _k, (_s, _a, _b) in enumerate(pairs):
        _row = []
        for _side, _i in enumerate((_a, _b)):
            c0 = _k * 4 + _side * 2
            ax = fig.add_subplot(gs[1, c0:c0 + 2])
            _row.append(ax)
            ax.imshow(_draw(mol_df.loc[_i, 'SMILES'],
                            ax.get_position().width * fig.get_figwidth() * 72, 380, 300),
                      interpolation='none')
            ax.axis('off')
            _t = int(mol_df.loc[_i, 'Tier'])
            ax.text(0.5, 1.04, f'Tier{_t}', transform=ax.transAxes, ha='center', va='baseline',
                    fontsize=7.4, fontweight='bold', color=TIER_TXT[_t],
                    bbox=dict(boxstyle='round,pad=0.26', facecolor=TIER_COL[_t], edgecolor='none'))
            ax.text(0.5, -0.02, textwrap.fill(_short(mol_df.loc[_i, 'name'], 24), width=14, break_long_words=False), transform=ax.transAxes,
                    ha='center', va='top', fontsize=6.4, color=MUTED)
            if _side == 0 and _k == 0:
                panel(ax, 'b', x=-0.04, y=1.46)
        _pair_axes.append((_row[0], _row[1], _s))

    # Similarity labels are drawn at figure level so they sit above the structures.
    for _axL, _axR, _s in _pair_axes:
        _pL, _pR = _axL.get_position(), _axR.get_position()
        _x = (_pL.x0 + _pR.x1) / 2
        fig.text(_x, 0.035, f'Tanimoto $T_c$ = {_s:.2f}', ha='center',
                 va='baseline', fontsize=7.2, color=INK)
    save(fig, 'Fig6_chemical_structures', align={'require_panel_labels': False})

print('Done:', len(list(OUT.glob('*.pdf'))), 'vector figures')
