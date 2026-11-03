"""Extract public case metadata with verified frozen JSON hashes; no fits."""
from pathlib import Path
import csv
import hashlib
import json

from configuration import load
CFG=load()
p=Path(CFG['case_manifest']);raw=p.read_bytes();d=json.loads(raw)
assert hashlib.sha256(raw).hexdigest()=='92e9b47a3fd32fb2d85bdd87da2d2797193bbdde9f6ade9969705ad8e1e6703c'
rows=[]
for task in d['rows']:
    f=Path(task['case_path']);raw=f.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    assert digest==task['hashes']['case.json'],f
    case=json.loads(raw);pr=case['provenance'];region=case['region']
    assert task['case_id']==case['case_id'] and task['p']==case['n_variants']
    rows.append(dict(panel=task['panel'],source_index=task['index'],case_id=case['case_id'],
                     chromosome=region[0],start_1based_inclusive=region[1],end_1based_inclusive=region[2],
                     variants=case['n_variants'],analysis_n=case['study_n'],reference_n=case['reference_n'],
                     reference_population=case['population'],genome_build=case['genome_build'],
                     source=pr['source'],release=pr['release'],gene_name=pr.get('gene_name',''),
                     qc_policy=case['qc_policy'],case_json_sha256=digest,
                     z_sha256=task['hashes']['z.bin'],R_sha256=task['hashes']['R.bin'],
                     variant_order_sha256=task['hashes']['variant_ids.txt'],
                     historical_R_reference_sha256=task['frozen_reference_sha256']))
assert len(rows)==14357 and sum(r['panel']=='GWAS' for r in rows)==200
out=Path(CFG['public_manifest'])
with out.open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print('Verified metadata and wrote',len(rows),'case-coordinate rows')
