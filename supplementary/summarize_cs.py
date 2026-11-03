"""Merge shards into results/cs_membership.tsv and results/summary.json."""
import json, glob, csv
from collections import Counter
from pathlib import Path
from configuration import load
CFG=load()
RESULTS=Path(CFG["cs_results"])
RESULTS.mkdir(parents=True,exist_ok=True)
rows = sorted((json.loads(l) for f in glob.glob(str(RESULTS/'shard_*.jsonl')) for l in open(f)), key=lambda r: r['case_key'])
assert len(rows) == 14357 and len({r['case_key'] for r in rows}) == 14357
cols = ['case_key','case_id','panel','p','R_converged','R_niter','prusie_converged','prusie_niter',
        'R_saved_cs_count','prusie_saved_cs_count','R_derived_cs_count','prusie_derived_cs_count',
        'n_matched','n_identical_members','identical','R_cs_components','prusie_cs_components',
        'only_R_components','only_prusie_components','max_abs_alpha_error','max_abs_purity_diff',
        'max_abs_coverage_diff','R_trimmed_rows','prusie_trimmed_rows','n_raw_coverage_set_differences','membership_differences']
with open(RESULTS/'cs_membership.tsv','w',newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(cols)
    for r in rows:
        r['n_raw_coverage_set_differences'] = len(r['raw_coverage_set_differences'])
        w.writerow([json.dumps(r[c]) if isinstance(r[c], (list, dict)) else r[c] for c in cols])
summ = {}
for panel in ('GWAS', 'eQTL'):
    s = [r for r in rows if r['panel'] == panel]
    summ[panel] = dict(
        cases=len(s),
        derived_eq_saved_R=sum(r['R_derived_cs_count'] == r['R_saved_cs_count'] for r in s),
        derived_eq_saved_prusie=sum(r['prusie_derived_cs_count'] == r['prusie_saved_cs_count'] for r in s),
        saved_count_constant_across_samples=sum(len(r['R_saved_cs_count_all_samples']) == 1 and len(r['prusie_saved_cs_count_all_samples']) == 1 for r in s),
        saved_counts_equal=sum(r['R_saved_cs_count'] == r['prusie_saved_cs_count'] for r in s),
        same_cs_component_ids=sum(r['R_cs_components'] == r['prusie_cs_components'] for r in s),
        identical_membership=sum(r['identical'] for r in s),
        total_cs_R=sum(r['R_derived_cs_count'] for r in s), total_cs_prusie=sum(r['prusie_derived_cs_count'] for r in s),
        cases_with_zero_cs=sum(r['R_derived_cs_count'] == 0 for r in s),
        trimmed_row_pattern_differs=sum(r['R_trimmed_rows'] != r['prusie_trimmed_rows'] for r in s),
        cases_with_raw_coverage_set_differences=sum(bool(r['raw_coverage_set_differences']) for r in s),
        raw_diff_cases_by_niter_cap=dict(Counter(f"R_conv={r['R_converged']},prusie_conv={r['prusie_converged']}" for r in s if r['raw_coverage_set_differences'])),
        max_purity_diff=max(r['max_abs_purity_diff'] for r in s), max_coverage_diff=max(r['max_abs_coverage_diff'] for r in s),
        nonconverged_both=sum(not r['R_converged'] and not r['prusie_converged'] for r in s),
        nonconverged_both_identical=sum(r['identical'] for r in s if not r['R_converged'] and not r['prusie_converged']))
summ['non_identical_cases'] = [{k: r[k] for k in ('case_key','case_id','panel','p','R_converged','prusie_converged','R_niter','prusie_niter',
    'only_R_components','only_prusie_components','membership_differences','max_abs_alpha_error')} for r in rows if not r['identical']]
summ['count_mismatch_cases'] = [r['case_key'] for r in rows if r['R_derived_cs_count'] != r['R_saved_cs_count'] or r['prusie_derived_cs_count'] != r['prusie_saved_cs_count']]
summ['trim_differs_cases'] = [{k: r[k] for k in ('case_key','case_id','R_trimmed_rows','prusie_trimmed_rows','raw_coverage_set_differences')} for r in rows if r['R_trimmed_rows'] != r['prusie_trimmed_rows']]
summ['method'] = ('NumPy re-implementation of susieR 0.16.6 susie_get_cs (in_CS: stable decreasing order, smallest prefix with long-double cumsum >= coverage; '
                  'dedup by first occurrence; V filter via trimmed rows = constant alpha row; min |r| purity >= 0.5 over upper triangle of R[cs,cs], early exit) '
                  'applied identically to both backends; cross-checked against prusie.susie_get_cs(alpha, Xcorr=R) on 76 sampled backend fits.')
json.dump(summ, open(RESULTS/'summary.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k in ('GWAS', 'eQTL')}, indent=1))
print('non-identical', len(summ['non_identical_cases']), 'count mismatch', summ['count_mismatch_cases'][:20], 'trim differs', summ['trim_differs_cases'][:5])
