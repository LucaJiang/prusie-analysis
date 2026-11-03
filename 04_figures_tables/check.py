"""Recompute numerical derivatives and compare decoded JSON/TSV values exactly."""
from pathlib import Path
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / '03_results'
OUTPUTS = ('benchmark_primary.tsv', 'benchmark_windows.tsv', 'benchmark_summary.json', 'derived_analyses.json')

def values(path):
    if path.suffix == '.json':
        return json.loads(path.read_text())
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def main():
    build = ROOT / 'build'
    build.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=build) as tmp:
        target = Path(tmp)
        for name in ('ld_class_census.tsv', 'thread_scaling.tsv'):
            shutil.copyfile(DATA / name, target / name)
        env = dict(os.environ, ANALYSIS_OUTPUT=str(target))
        subprocess.run([sys.executable, str(Path(__file__).with_name('summarize_benchmark.py'))], env=env, check=True, stdout=subprocess.DEVNULL)
        env['ANALYSIS_INPUT'] = str(target)
        subprocess.run([sys.executable, str(Path(__file__).with_name('derived_analyses.py'))], env=env, check=True, stdout=subprocess.DEVNULL)
        for name in OUTPUTS:
            assert values(target/name) == values(DATA/name), f'Numerical mismatch: {name}'
            print('Exact values:', name)
        subprocess.run([sys.executable, str(Path(__file__).with_name('manuscript_numbers.py')), '--output', str(target/'manuscript_numbers.json')], check=True)
        expected = DATA / 'manuscript_numbers.json'
        if expected.exists():
            assert values(target/'manuscript_numbers.json') == values(expected), 'Manuscript number mismatch'
            print('Exact values: manuscript_numbers.json')
        else:
            raise RuntimeError('Missing frozen manuscript_numbers.json')
    checksums = json.loads((DATA / 'checksums.json').read_text())
    for name, digest in checksums.items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest() == digest, name
    print(f'Frozen source checksums: {len(checksums)} files')

if __name__ == '__main__':
    main()
