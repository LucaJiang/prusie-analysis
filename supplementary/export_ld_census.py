"""Merge complete census JSONL shards into the seven-column public table."""
import argparse
import csv
import json
from pathlib import Path

FIELDS = ('case_key', 'case_id', 'panel', 'p', 'G', 'largest_class', 'variants_in_shared_classes')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary-table', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('shards', type=Path, nargs='+')
    args = parser.parse_args()
    with args.primary_table.open() as handle:
        primary = [r for r in csv.DictReader(handle, delimiter='\t') if r['kind'] == 'primary']
    by_key = {r['case_key']: r for r in primary}
    if len(by_key) != len(primary):
        raise ValueError('Duplicate primary case keys')
    census = {}
    for path in args.shards:
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            key = row['case_key']
            if key in census:
                raise ValueError(f'Duplicate census case: {key}')
            expected = by_key[key]
            if (row['case_id'], row['panel'], int(row['p'])) != (
                    expected['case_id'], expected['panel'], int(expected['p'])):
                raise ValueError(f'Census identity or variant count differs: {key}')
            if not (1 <= int(row['G']) <= int(row['p']) and
                    1 <= int(row['largest_class']) <= int(row['p']) and
                    0 <= int(row['variants_in_shared_classes']) <= int(row['p'])):
                raise ValueError(f'Invalid class counts: {key}')
            census[key] = {name: row[name] for name in FIELDS}
    if set(census) != set(by_key):
        raise ValueError('Census must cover all primary GWAS and eQTL cases')
    with args.output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(census[key] for key in sorted(census))
    print(f'Exported all {len(census)} primary cases to {args.output}')


if __name__ == '__main__':
    main()
