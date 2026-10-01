"""Plot paired label-source endpoint contrasts (Supplementary Fig. S4).

Question: how do endpoint contrasts change when saved held-out scores are evaluated
against source-retrieved labels? All three splits and both primary contrasts
are shown for broad and route-matched cohorts. No model or label is fitted.
"""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
try:  # optional panel-alignment auditor (see output/figure_qa/ for the recorded reports)
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    require_matplotlib_panel_alignment = None

ROOT=Path(__file__).resolve().parents[1]
QA=ROOT/'output/figure_qa'
QA.mkdir(parents=True,exist_ok=True)
DATA=ROOT/'output/external_evaluation/frozen_contrasts_ECFP.csv'
source=pd.read_csv(DATA)
blocks=[
    ('OFT_single_structure_contact','all_external_contact','OFT: broad'),
    ('PLOS_resolved_contact','all_external_contact','EPA: broad'),
    ('OFT_single_structure_contact','original_contact_only','OFT: matched route'),
    ('PLOS_resolved_contact','original_contact_quality_screen','EPA: matched route + quality'),
]
rows=[]
for block,(src,cohort,label) in enumerate(blocks):
    for split in ['random','maxmin','time']:
        for contrast in ['AUROC11-AUROC100','AUROC1-AUROC100']:
            g=source[source.source.eq(src)&source.cohort.eq(cohort)&source['split'].eq(split)&source.comparison.eq(contrast)&source.labels.isin(['original','external'])]
            assert len(g)==2
            for _,r in g.iterrows():
                rows.append({**r.to_dict(),'block':block,'block_label':label})
table=pd.DataFrame(rows)
table.to_csv(ROOT/'source_data/figure_external_reassessment.csv',index=False,encoding='utf-8-sig')

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':7.5,
    'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False})
fig,ax=plt.subplots(figsize=(7.16,7.0))
fig.subplots_adjust(left=.37,right=.985,top=.925,bottom=.18)
colors={'original':'#777E87','external':'#126A8A'}
markers={'original':'o','external':'s'}
ticks=[]; ticklabels=[]
for block,(src,cohort,label) in enumerate(blocks):
    if block%2==1:
        ax.axhspan(block*7-.55,block*7+5.55,color='#F3F5F7',zorder=0)
    for si,split in enumerate(['random','maxmin','time']):
        for ci,contrast in enumerate(['AUROC11-AUROC100','AUROC1-AUROC100']):
            y=block*7+si*2+ci
            g=table[table.block.eq(block)&table['split'].eq(split)&table.comparison.eq(contrast)].set_index('labels')
            splitname={'random':'Random','maxmin':'MaxMin','time':'Time'}[split]
            comp='11 - 100' if ci==0 else '1 - 100'
            ticks.append(y);ticklabels.append(f'{splitname}, {comp} (n={int(g.n_compounds.iloc[0])})')
            values=g.estimate
            if np.isfinite(values).all():
                ax.plot([values['original'],values['external']],[y-.14,y+.14],color='#BAC3CB',lw=.65,zorder=1)
            for lab,offset in [('original',-.14),('external',.14)]:
                r=g.loc[lab]
                if not np.isfinite(r.estimate):
                    continue
                sparse=bool(r.exploratory_sparse)
                if not sparse and np.isfinite(r.ci_low) and np.isfinite(r.ci_high):
                    ax.plot([r.ci_low,r.ci_high],[y+offset,y+offset],color=colors[lab],lw=.8,zorder=2)
                    ax.plot([r.ci_low,r.ci_high],[y+offset,y+offset],linestyle='none',marker='|',color=colors[lab],ms=4,zorder=2)
                ax.plot(r.estimate,y+offset,marker=markers[lab],ms=4.2,linestyle='none',
                    markerfacecolor='white' if sparse else colors[lab],markeredgecolor=colors[lab],markeredgewidth=.85,zorder=3)
    ax.text(-.57,block*7-.94,label,transform=ax.get_yaxis_transform(),ha='left',va='bottom',fontweight='bold',fontsize=8,clip_on=False)
ax.axvline(0,color='#C0C6CC',ls=':',lw=.85,zorder=0)
ax.set_yticks(ticks,ticklabels)
ax.tick_params(axis='y',length=0,pad=7)
ax.set_xlim(-.65,1.05)
ax.set_xticks([-.5,0,.5,1])
ax.set_ylim(27.0,-1.25)
ax.set_xlabel('Difference in AUROC (severe endpoint minus 100 µg/bee)')
ax.spines[['top','right','left']].set_visible(False)
ax.spines['bottom'].set_color('#777E87')
handles=[Line2D([],[],color=colors[k],marker=markers[k],linestyle='none',ms=4.5,
    label='Original labels' if k=='original' else 'External-source labels') for k in colors]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.985,.995),ncol=2,frameon=False,fontsize=8)
fig.text(.03,.033,'Matched compounds and saved ECFP scores within each row. Open markers: a constituent class has fewer than five compounds.\n'
    'Intervals are omitted for sparse comparisons; filled-marker intervals are exploratory 95% paired bootstrap intervals.\n'
    'An absent marker denotes an undefined contrast. Complete counts, intervals and eligibility rules are supplied in the Supplementary Information.',
    fontsize=7,ha='left',va='bottom',color='#40464F')
fig.canvas.draw()
if require_matplotlib_panel_alignment is not None:
 require_matplotlib_panel_alignment(fig,json_out=str(QA/'FigS4_label_source_contrasts.alignment.json'),
    overlay_svg=str(QA/'FigS4_label_source_contrasts.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
path=ROOT/'figures/FigS4_label_source_contrasts'
fig.savefig(path.with_suffix('.pdf'))
fig.savefig(path.with_suffix('.svg'))
fig.savefig(path.with_suffix('.png'),dpi=600)
fig.savefig(path.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
(QA/'FigS4_contract.json').write_text(json.dumps({'question':'Sensitivity of saved-score endpoint contrasts to label source',
    'archetype':'single quantitative paired-contrast panel','width_inches':7.16,
    'all_splits_and_primary_contrasts_retained':True,'n_rows':24,'n_estimates':len(table),
    'uncertainty_exception':'Omit intervals when a constituent class has fewer than5 compounds; keep the estimates and full machine-readable intervals.',
    'source':str(DATA.relative_to(ROOT))},indent=2),encoding='utf-8')
print('Exported Supplementary Figure S4 with all 24 paired rows at final journal width.')
