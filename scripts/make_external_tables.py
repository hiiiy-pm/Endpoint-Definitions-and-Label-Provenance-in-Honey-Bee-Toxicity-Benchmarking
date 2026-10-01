from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path.cwd()
lines=[]
def add(s): lines.append(s)
def fmt(x,n=3): return '--' if pd.isna(x) else f'{x:.{n}f}'
def esc(x):
 return str(x).replace('&',r'\&').replace('_',r'\_').replace('%',r'\%').replace('>',r'$>$')
def table(title,spec,head,rows,note):
 add(r'\Needspace{12\baselineskip}\subsection*{'+title+'}')
 add(note)
 add(r'\begingroup\small\setlength{\tabcolsep}{4pt}')
 add(r'\begin{longtable}{'+spec+'}')
 add(r'\toprule '+head+r' \\ \midrule\endfirsthead')
 add(r'\toprule '+head+r' \\ \midrule\endhead')
 add(r'\midrule\multicolumn{'+str(len(head.split('&')))+r'}{r}{Continued on next page}\endfoot\bottomrule\endlastfoot')
 add('\n'.join(' & '.join(row)+r' \\' for row in rows))
 add(r'\end{longtable}\endgroup')
# S8: every prespecified descriptor subset; no result-guided choice of k.
d=pd.read_csv(ROOT/'results/supplementary/02_topk_descriptor_performance.csv')
rows=[]
for sp in ['Random','MaxMin','Time']:
 for k in [1,2,3,5,8,12]:
  q=d[(d.Split==sp)&(d.TopK==k)].set_index('Threshold')
  rows.append([sp,str(k)]+[fmt(q.loc[t,'TestAUROC']) for t in [100,11,1]])
table('Supplementary Table S8. Nested descriptor sensitivity','lrccc',r'Split & $k$ & $\leq100$ & $\leq11$ & $\leq1$',rows,
 r'Test AUROC for every prespecified subset size, with training-fold-specific scaling and stability selection. The official test set remains 207 compounds. Selected descriptors and training-selected penalties for all 54 fits are in \texttt{results\_v2/supplementary/02\_topk\_descriptor\_performance.csv}. This is a separate L2-logistic sensitivity analysis; its 12-descriptor fit is not the main descriptor SVM.')
# S15 fixed concordance cohorts.
d=pd.read_csv(ROOT/'output/external_data_audit/label_concordance_summary.csv')
rows=[]
for source in ['OpenFoodTox','OpenFoodTox_48h_screen','PLOS2022']:
 q=d[(d.source==source)&(d.cohort=='all_three_determinate')].set_index('threshold_ug_bee')
 rows.append([{'OpenFoodTox':'OFT','OpenFoodTox_48h_screen':'OFT, 48 h','PLOS2022':'EPA retrospective'}[source],str(int(q.loc[1,'n']))]+[str(int(q.loc[t,'disagreements'])) for t in [1,11,100]]+[str(int(q.loc[100,'discordance_with_any_right_censored_100']))])
table('Supplementary Table S15. Route-matched source-label concordance','lrcccc',r'Source & $n$ & At 1 & At 11 & At 100 & At 100, censored',rows,
 r'Counts are disagreements against the original ApisTox labels on the same compounds with all three source labels determined. The final column counts 100-cut-off discrepancies containing right-censored $>100$ evidence. These are source disagreements rather than adjudicated extraction errors. Threshold-specific cohorts, compound identities and underlying observations are retained in the concordance CSV files.')
# S18 all ECFP cohorts and labels; retain every available source/split/screen.
m=pd.read_csv(ROOT/'output/external_evaluation/frozen_metrics_ECFP.csv')
m=m[m.metric.eq('AUROC')]
source_name={'OFT_single_structure_contact':'OFT','OFT_48h_contact':'OFT48','PLOS_resolved_contact':'EPA'}
cohort_name={'all_external_contact':'Broad','original_contact_only':'Route','quality_screen':'Quality','original_contact_quality_screen':'Route+quality'}
rows=[]
for (src,co,sp,lab),q in m.groupby(['source','cohort','split','labels'],sort=False):
 q=q.set_index('threshold')
 vals=[]
 for t in [100,11,1]:
  a=q.loc[t]; sparse='*' if bool(a.exploratory_sparse) else ''
  vals.append(fmt(a.estimate)+sparse+' ('+str(int(a.n_pos))+'/'+str(int(a.n_neg))+')')
 rows.append([source_name[src],cohort_name[co],sp.title().replace("Maxmin","MaxMin"),('Orig.' if lab=='original' else 'Source'),str(int(q.n_compounds.iloc[0]))]+vals)
table('Supplementary Table S18a. Complete saved-score ECFP label reassessment','llllrccc',r'Source & Screen & Split & Labels & $n$ & $\leq100$ & $\leq11$ & $\leq1$',rows,
 r'Cells report AUROC and (positive/negative) counts. An asterisk marks fewer than five compounds in a class; -- denotes undefined AUROC. Orig.: original benchmark labels; Source: regulatory-source labels. Broad: external contact, original route unrestricted; Route: original contact required; OFT48: 48-hour observations required; EPA Quality: structure-quality and MRID evidence required. All cohorts exclude the corresponding model\textquotesingle s training identities. The cohorts and sources overlap. Exact eligibility, AP, pointwise percentile intervals and degenerate bootstrap counts remain in the machine-readable tables.')
# All contrast rows with interval; separate label types; preserve flags.
c=pd.read_csv(ROOT/'output/external_evaluation/frozen_contrasts_ECFP.csv')
rows=[]
for (src,co,sp,lab),q in c.groupby(['source','cohort','split','labels'],sort=False):
 q=q.set_index('comparison'); vals=[]
 for key in ['AUROC11-AUROC100','AUROC1-AUROC100']:
  a=q.loc[key]
  suffix='*' if bool(a.exploratory_sparse) else ''
  vals.append(fmt(a.estimate)+suffix+' ['+fmt(a.ci_low)+', '+fmt(a.ci_high)+']')
 rows.append([source_name[src],cohort_name[co],sp.title().replace("Maxmin","MaxMin"),{'original':'Orig.','external':'Source','external_minus_original':'Change'}[lab]]+vals)
table('Supplementary Table S18b. Paired endpoint contrasts and label-source changes','llllcc',r'Source & Screen & Split & Labels & $11-100$ & $1-100$',rows,
 r'Cells report paired AUROC differences with exploratory 95\% compound-bootstrap percentile intervals. Change is the source-label contrast minus the original-label contrast on identical compounds. Asterisks flag sparse classes. Conditional intervals after removing single-class bootstrap draws can be misleading, especially with one positive; they are retained for audit and do not establish replication. These intervals are not the original scaffold-cluster simultaneous intervals.')
# S19 raw new-chemical evidence and lineage locations.
e=pd.read_csv(ROOT/'output/unseen_identity_case_study/data_curation/frozen_scoreable_48h_endpoints.csv')
rows=[]
for name,q in e.groupby('candidate',sort=False):
 vals=[]
 for route in ['contact','oral']:
  a=q[q.route==route]
  vals.append('--' if not len(a) else esc(a.iloc[0].qualifier)+fmt(a.iloc[0].LD50_value,3).rstrip('0').rstrip('.'))
 adult='Explicit' if q.life_stage_evidence.str.startswith('explicit').all() else 'Context'
 rows.append([name.capitalize()]+vals+[adult])
table('Supplementary Table S19a. Eligible unseen-identity 48-hour endpoints','lccc',r'Chemical identity & Contact LD$_{50}$ & Oral LD$_{50}$ & Adult evidence',rows,
 r'Values retain their original qualifiers and are reported in $\mu$g per bee. Active-ingredient dose basis is indicated only where stated by the source; test-material and dose-basis limitations remain in the endpoint ledger. -- means no eligible endpoint, not no toxicity. Adult evidence is explicit in the source or inferred from regulatory separation of adult acute and larval studies. The latter cases are excluded from the explicit-adult sensitivity. For acynonapyr, $>98.8$ cannot determine the 100-cut-off label; for cyclobutrifluram, $>72.2$ cannot determine it. Same-connectivity stereoisomer mixtures are evaluated by the non-chiral model without asserting exact material-level stereochemical identity.')
add(r'\subsection*{Supplementary Table S19b. Source traceability and eligibility exclusions}')
for name,q in e.groupby('candidate',sort=False):
 add(r'\paragraph{'+name.capitalize()+r'.}')
 for _,a in q.iterrows():
  add(esc(a.route.capitalize())+': '+esc(a.study_id)+'. '+esc(a.source_location)+'. '+r'\url{'+a.source_url+r'}.')
add(r'All 11 accepted endpoints are in \texttt{frozen\_scoreable\_48h\_endpoints.csv}, including source material, stage evidence, full key and connectivity. The complete \texttt{frozen\_endpoint\_ledger.csv} retains rejected observations and reasons. The ten predicted candidate structures include four identities without an accepted 48-hour endpoint; they were not removed according to prediction success. Inpyrfluxam lacks adequate stage confirmation. Longer-time observations are stored separately in \texttt{extended\_times\_not\_primary.csv}. Sources not downloaded locally remain identified as full-text checked or unavailable in the source manifests; their local availability is not overstated.')
# S20 all pilot metrics, transparent confusion matrices and sizes.
pdft=pd.read_csv(ROOT/'output/unseen_identity_case_study/evaluation/pilot_performance.csv')
rows=[]
for _,a in pdft.iterrows():
 rows.append([a.route.capitalize(),('All' if a.quality_screen=='regulatory_context' else 'Explicit'),('Each' if a.cohort=='threshold_specific' else 'Common'),a.split,str(int(a.threshold)),f'{int(a.n_compounds)}/{int(a.n_positive)}',fmt(a.AUROC),'/'.join(str(int(a[k])) for k in ['TP','TN','FP','FN'])])
table('Supplementary Table S20. Complete exploratory unseen-identity evaluation','llllrrcc',r'Route & Stage & Cohort & Split & Cut-off & $n$/pos. & AUROC & TP/TN/FP/FN',rows,
 r'Stage All retains contextual adult evidence; Explicit requires explicit evidence. Each uses determinate labels at that cut-off; Common requires all three labels determined. Confusion counts use zero decision score as the positive boundary. All class configurations are sparse; single-class AUROC is --. The full CSV supplies AP, exact sensitivity/specificity intervals, exploratory bootstrap intervals and valid/degenerate draws. These are selected pilot cases, not population-level performance estimates. Contact labels at 11 and 100 are identical in the primary cohort; common oral labels are identical at all three cut-offs. Their AUROC contrasts compare trained-model rankings rather than a fresh threshold-relabeling effect.')
(ROOT/'output'/'latex_tables'/'supp_tables_external.tex').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')
print('Supplementary Tables S8, S15 and S18-S20 written.')
