"""Freeze fixed eight-compound inference panel before examining predictions."""
from pathlib import Path
import pandas as pd,json,hashlib
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
OUT=Path(__file__).resolve().parent
S=pd.read_csv(OUT/'all_candidate_structures.csv').fillna('')
S=S[S.candidate_origin.eq('previous8')].copy()
assert len(S)==8 and not S[['overlap_inchikey','overlap_connectivity','overlap_parent_connectivity']].any().any()
S['compound_id']=S.candidate
S['name']=S.candidate
S['SMILES']=S.canonical_smiles
S['material_identity_status']='connectivity_supported; full_material_stereochemical_identity_not_claimed'
S['stereochemistry_note']=S.candidate.map({
 'cyclobutrifluram':'Registered material comprises two enantiomers of the same connectivity; PubChem name represents one enantiomer.',
 'acynonapyr':'Regulator specifies3-endo; PubChem representation has partial stereochemical specification.',
 'isocycloseram':'Technical material may comprise stereoisomers; PubChem representation has unspecified stereochemistry.'}).fillna('No equivalence of registry product and exact study batch is claimed.')
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024,includeChirality=False)
ok=[]
for smi in S.SMILES:
    m=Chem.MolFromSmiles(smi);m2=Chem.Mol(m);Chem.RemoveStereochemistry(m2)
    ok.append(gen.GetFingerprint(m).ToBitString()==gen.GetFingerprint(m2).ToBitString())
assert all(ok)
S['fingerprint_equal_after_removing_stereo']=ok
p=OUT/'frozen_initial8_inference_structures.csv'
S[['compound_id','name','SMILES','inchikey','connectivity','parent_connectivity','PubChem_CID','material_identity_status','stereochemistry_note','fingerprint_equal_after_removing_stereo']].to_csv(p,index=False,encoding='utf-8-sig')
(OUT/'frozen_initial8_manifest.json').write_text(json.dumps({'frozen_before_prediction_review':True,'n_compounds':8,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'contains_endpoint_labels':False,'fingerprint_verification':{'radius':2,'fpSize':1024,'includeChirality':False},'policy':'Single chemical connectivity is the experimental unit for exploratory application of chirality-insensitive ECFP. Same-connectivity stereoisomer mixtures may be represented, with full material stereochemistry unconfirmed. Mixtures with different active chemical connectivities are excluded. Model output must not determine endpoint inclusion.'},indent=2),encoding='utf-8')
print(p)
