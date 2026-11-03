#!/usr/bin/env python3
"""Bind prepared cases to the fixed manuscript order before running any fits."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import yaml
from campaign import resolve
from common import read, sha, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--allow-missing', action='store_true', help='For a fixed smoke subset: preserve all panel entries but report unavailable cases')
    args = parser.parse_args()
    source = args.config.resolve()
    config = yaml.safe_load(source.read_text())
    inputs, benchmark = config['inputs'], config['benchmark']
    scratch = resolve(inputs['scratch'], source.parent)
    panel = resolve(benchmark.get('primary_case_manifest', '../03_results/primary_case_manifest.tsv'), source.parent)
    output = resolve(benchmark['cases_manifest'], source.parent)
    rows, receipts, missing = [], [], []
    with panel.open() as handle:
        fixed = list(csv.DictReader(handle, delimiter='\t'))
    assert len(fixed) == 14357
    for row in fixed:
        case_path = scratch / 'plink_ld' / row['case_id'] / 'case.json'
        item = dict(case_key=f"{row['panel']}_{int(row['source_index']):05d}", panel=row['panel'],
                    case_path=str(case_path), window=None, kind='primary')
        rows.append(item)
        if not case_path.exists():
            missing.append(item['case_key'])
            continue
        case = read(case_path)
        assert case['case_id'] == row['case_id']
        assert case['n_variants'] == int(row['variants'])
        assert case['parameters']['n'] == int(row['analysis_n'])
        assert case['parameters']['max_iter'] == 100
        assert case['parameters']['L'] == 10
        assert case['parameters']['n_purity'] == int(row['variants'])
        checks = {key: sha(case_path.parent / case[field]) == row[key] for key, field in
                  [('z_sha256', 'z_path'), ('variant_order_sha256', 'variant_ids_path')]}
        if not all(checks.values()):
            raise ValueError(f"Input digest differs from the frozen order: {item['case_key']}: {checks}")
        receipts.append(dict(case_key=item['case_key'],checks=checks,
                              regenerated_case_sha256=sha(case_path),historical_case_sha256=row['case_json_sha256'],
                              case_json_note='Portable provenance differs; signed z and variant-order bytes are checked directly.'))
    largest = max((row for row in fixed if row['panel'] == 'GWAS'), key=lambda r: (int(r['variants']), r['case_id']))
    assert largest['case_id'] == 'T2D_chr10_61277948_71054950'
    assert int(largest['variants']) == 19221
    for size in (500, 1000, 2000, 4000, 10000, 19221):
        rows.append(dict(case_key=f'window_{size}',panel='windows',
                         case_path=str(scratch / 'plink_ld' / largest['case_id'] / 'case.json'),
                         window=size,kind=f'window{size}'))
    if missing and not args.allow_missing:
        raise FileNotFoundError(f'{len(missing)} primary cases are missing; no manifest written')
    write(output, rows)
    write(output.with_suffix('.inputs.json'), dict(primary_case_manifest_sha256=sha(panel),
          total_entries=len(rows),available_primary_cases=len(receipts),missing_primary_cases=missing,inputs=receipts))
    print(f'Wrote {len(rows)} fixed entries, {len(receipts)} available primary cases, {len(missing)} missing.')


if __name__ == '__main__':
    main()
