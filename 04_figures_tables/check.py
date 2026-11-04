"""Recompute numerical derivatives and compare them with the published values.

Non-numeric values must match exactly. Floating-point values may differ by a
relative 1e-12, because BLAS kernels and SIMD paths differ between CPUs.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
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

RTOL = 1e-12


def number(value):
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compare(got, want, where):
    if isinstance(want, dict):
        assert isinstance(got, dict) and got.keys() == want.keys(), f'{where}: keys differ'
        for key in want:
            compare(got[key], want[key], f'{where}.{key}')
    elif isinstance(want, list):
        assert isinstance(got, list) and len(got) == len(want), f'{where}: lengths differ'
        for i, (g, w) in enumerate(zip(got, want)):
            compare(g, w, f'{where}[{i}]')
    elif got != want:
        g, w = number(got), number(want)
        assert g is not None and w is not None, f'{where}: {got!r} != {want!r}'
        assert math.isclose(g, w, rel_tol=RTOL, abs_tol=0.0), f'{where}: {got!r} != {want!r}'


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
            compare(values(target / name), values(DATA / name), name)
            print('Matched:', name)
        subprocess.run([sys.executable, str(Path(__file__).with_name('manuscript_numbers.py')), '--output', str(target/'manuscript_numbers.json')], check=True)
        expected = DATA / 'manuscript_numbers.json'
        if expected.exists():
            compare(values(target / 'manuscript_numbers.json'), values(expected), 'manuscript_numbers.json')
            print('Matched: manuscript_numbers.json')
        else:
            raise RuntimeError('Missing frozen manuscript_numbers.json')
    checksums = json.loads((DATA / 'checksums.json').read_text())
    for name, digest in checksums.items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest() == digest, name
    print(f'Frozen source checksums: {len(checksums)} files')

if __name__ == '__main__':
    main()
