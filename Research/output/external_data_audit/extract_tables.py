"""Extract source tables without changing supplied workbooks."""
from pathlib import Path
import csv
import json
import pandas as pd

BASE = (Path(__file__).resolve().parents[2] / 'data/external_raw')
OUT = Path(__file__).resolve().parent
CACHE = OUT/'cache'
CACHE.mkdir(exist_ok=True)
for name in ['REF_SUB', 'SUB', 'END_STUDY_REC.TerrestEcotox', 'LIT']:
    d = pd.read_excel(BASE/'OFT3.0 export repository.xlsx', sheet_name=name)
    d.insert(0, 'excel_row', range(2, len(d)+2))
    d.to_pickle(CACHE/f'{name}.pkl')
    print(name, d.shape)
    if name == 'END_STUDY_REC.TerrestEcotox':
        species_cols = [c for c in d.columns if 'Species' in c]
        mask = d[species_cols].fillna('').astype(str).apply(lambda s: s.str.contains('apis|bee|bombus|osmia',case=False,regex=True)).any(axis=1)
        bees = d.loc[mask].copy()
        bees.to_csv(OUT/'oft_bee_source_rows.csv', index=False, encoding='utf-8-sig')
        print('BEE ROWS', len(bees))
        for c in bees.columns:
            if bees[c].notna().any() and any(k in c for k in ['Species','DoseMethod','Duration','Endpoint','EffectLevel','EffectConc.','LifeStage','StudyType']):
                print(c, bees[c].value_counts(dropna=False).head(16).to_dict())
        print('BEE_EXAMPLES', bees.dropna(axis=1,how='all').head(4).to_json(orient='records',force_ascii=False))
rows = list(csv.reader(open(BASE/'excluded_data.csv',encoding='utf-8-sig',newline='')))
records = []
block = 0
for line, row in enumerate(rows, 1):
    if not row:
        continue
    if 'Exclusion reason' in row:
        header = row
        block += 1
        continue
    assert len(row) == len(header), (line, row, header)
    records.append(dict(zip(header,row), source_csv_record=line, source_block=block))
excluded = pd.DataFrame(records)
excluded.to_csv(OUT/'excluded_all_parsed.csv',index=False,encoding='utf-8-sig')
uncertain = excluded[excluded['Exclusion reason'].eq('Unspecifed toxicity level')]
print('EXCLUDED',len(excluded),excluded['Exclusion reason'].value_counts().to_dict())
print('UNCERTAIN',len(uncertain),uncertain['CAS'].nunique(),uncertain.groupby(['CAS','toxicity_type']).ngroups)
print('UNCERTAIN_NAMES',uncertain[['CAS','name','toxicity_type']].drop_duplicates().to_json(orient='records',force_ascii=False))
for i, header in [(1,6),(4,2),(5,2)]:
    d=pd.read_excel(BASE/f'journal.pone.0265962.s{i:03}.xlsx',header=header)
    d.insert(0,'excel_row', range(header+2,len(d)+header+2))
    d.to_pickle(CACHE/f'plos_s{i}.pkl')
    print('PLOS',i,d.shape, list(d.columns))
    if i == 1:
        print(d.head(8).to_json(orient='records',force_ascii=False))
        print('TAIL', d.tail(12).to_json(orient='records',force_ascii=False))
