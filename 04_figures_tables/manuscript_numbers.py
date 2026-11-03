"""Numerical tables and descriptive claims from the published records alone."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parents[1] / '03_results'

def tsv(name):
    with (DATA/name).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def stats(x):
    x = np.asarray(x, float)
    return dict(zip(('minimum','q05','q25','median','q75','q95','maximum'), map(float,np.quantile(x,[0,.05,.25,.5,.75,.95,1]))))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with gzip.open(DATA/'benchmark_records.jsonl.gz','rt') as f:
        records=[json.loads(line) for line in f if line.strip()]
    summary=json.loads((DATA/'benchmark_summary.json').read_text())
    result={'definitions':{'case_ratio':'median R calls / median prusie calls', 'case_median_total_ratio':'sum of R case medians / sum of prusie case medians', 'retained_total_ratio':'sum of all retained R calls / sum of all retained prusie calls'},'panels':{},'summary_tables':summary,'derived_tables':json.loads((DATA/'derived_analyses.json').read_text()),'nested_windows':tsv('benchmark_windows.tsv')}
    for panel in ('GWAS','eQTL'):
        rows=[r for r in records if r['panel']==panel]
        members=[r for r in tsv('cs_membership.tsv') if r['panel']==panel]
        manifest=[r for r in tsv('primary_case_manifest.tsv') if r['panel']==panel]
        # The frozen TSV provides full precision differences and exactly the case order.
        primary=[r for r in tsv('benchmark_primary.tsv') if r['panel']==panel]
        p=np.array([int(r['p']) for r in primary]); n=np.array([int(r['prusie_niter']) for r in primary])
        alpha=np.array([float(r['max_abs_alpha_error']) for r in primary]); pip=np.array([float(r['max_abs_pip_error']) for r in primary])
        ratio=np.array([float(r['speedup']) for r in primary]); conv=np.array([r['prusie_converged']=='True' and r['susieR_converged']=='True' for r in primary])
        metrics={}
        for metric,field in (('complete','times'),('fitting','fitting_times')):
            times={b:[r[b][field] for r in rows] for b in ('R','prusie')}
            totals={b:math.fsum(float(np.median(v)) for v in vals) for b,vals in times.items()}
            retained={b:math.fsum(t for v in vals for t in v) for b,vals in times.items()}
            outliers=[{'case_id':r['case_id'],'seconds':t,'case_median':float(np.median(r['R'][field]))} for r in rows for t in r['R'][field] if t>2*np.median(r['R'][field])]
            metrics[metric]={'sum_case_medians_seconds':totals,'case_median_total_ratio':totals['R']/totals['prusie'],'sum_retained_seconds':retained,'retained_total_ratio':retained['R']/retained['prusie'],'R_calls_above_twice_case_median':outliers,'R_outlier_fraction':math.fsum(o['seconds'] for o in outliers)/retained['R'],'retained_ratio_replacing_R_outliers':(retained['R']-math.fsum(o['seconds']-o['case_median'] for o in outliers))/retained['prusie']}
        entry={'cases':len(rows),'variants':stats(p),'analysis_n':sorted({int(r['analysis_n']) for r in manifest}),'reference_n':sorted({int(r['reference_n']) for r in manifest}),'timing':metrics,'speedup':stats(ratio),'speedup_geometric_mean':float(np.exp(np.log(ratio).mean())),'speedup_IQR_over_median':float((np.quantile(ratio,.75)-np.quantile(ratio,.25))/np.median(ratio)),'speedup_counts_above':{str(k):int((ratio>k).sum()) for k in (10,20,30,50)},'speedup_below_10':int((ratio<10).sum()),'spearman_p_speedup':float(spearmanr(p,ratio).statistic),'p_tertile_median_speedup':[float(np.median(ratio[ix])) for ix in np.array_split(np.argsort(p,kind='stable'),3)],'alpha':stats(alpha),'pip':stats(pip),'alpha_at_most_1e8_count':int((alpha<=1e-8).sum()),'alpha_at_most_1e8_percent':float(100*(alpha<=1e-8).mean()),'alpha_exact_zero_count':int((alpha==0).sum()),'converged_alpha_max':float(alpha[conv].max()),'converged_spearman_iterations_alpha':float(spearmanr(n[conv],alpha[conv]).statistic),'median_iterations':float(np.median(n)),'alpha_by_iterations':{},'credible_sets':{b:sum(int(r[b+'_cs_count']) for r in members) for b in ('susieR','prusie')},'different_cs_membership':sum(r['identical_membership']!='True' for r in members)}
        for lo,hi in [(2,2),(3,5),(6,10),(11,20),(21,100)]:
            mask=conv & (n>=lo) & (n<=hi)
            entry['alpha_by_iterations'][f'{lo}-{hi}']={'cases':int(mask.sum()),'median':float(np.median(alpha[mask])) if mask.any() else None}
        result['panels'][panel]=entry
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Wrote numerical tables and claims:',args.output.name)

if __name__=='__main__': main()
