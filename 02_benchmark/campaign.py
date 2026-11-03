#!/usr/bin/env python3
"""Run the frozen paired complete-call experiment from configured portable cases."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import yaml
from common import THREAD_VARS, read, sha, write

HERE = Path(__file__).resolve().parent


def resolve(value, base):
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def balanced_orders(rows):
    """The exact original seed and per-panel balanced shuffle."""
    rng = random.Random(20261028)
    assignments = {}
    for panel in ('GWAS', 'eQTL', 'windows'):
        selected = sorted((row for row in rows if row['panel'] == panel), key=lambda r: r['case_key'])
        rng.shuffle(selected)
        for index, row in enumerate(selected):
            assignments[row['case_key']] = ['R', 'prusie'] if index % 2 == 0 else ['prusie', 'R']
    if len(assignments) != len(rows):
        raise ValueError('Each case needs a unique case_key and panel GWAS, eQTL or windows')
    return assignments


def environment(config, output):
    settings = dict(os.environ, **{key: '1' for key in THREAD_VARS},
                    PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
                    R_PROFILE_USER='/dev/null', R_ENVIRON_USER='/dev/null')
    for name, key in [('r_libs', 'R_LIBS'), ('r_libs_user', 'R_LIBS_USER')]:
        if config.get(name):
            values = config[name] if isinstance(config[name], list) else [config[name]]
            settings[key] = os.pathsep.join(str(resolve(value, config['_base'])) for value in values)
    settings['TMPDIR'] = str(output / 'tmp')
    (output / 'tmp').mkdir()
    python = str(resolve(config['python'], config['_base'])) if '/' in config['python'] else config['python']
    rscript = str(resolve(config['rscript'], config['_base'])) if '/' in config['rscript'] else config['rscript']
    if config.get('blas_probe'):
        probe = resolve(config['blas_probe'], config['_base'])
    else:
        build = output / 'build'
        build.mkdir()
        shutil.copy2(HERE / 'blas_probe.c', build / 'blas_probe.c')
        rbuild = config.get('r_build') or str(Path(shutil.which(rscript) or rscript).with_name('R'))
        if '/' in rbuild:
            rbuild = str(resolve(rbuild, config['_base']))
        subprocess.run([rbuild, 'CMD', 'SHLIB', 'blas_probe.c', '-ldl'], cwd=build, env=settings, check=True)
        probe = build / 'blas_probe.so'
    subprocess.run([rscript, '--vanilla', str(HERE / 'inspect_reference.R'), str(probe),
                    str(output / 'R_identity.json')], env=settings, check=True)
    native = subprocess.check_output([python, '-c',
      'import json,prusie;from prusie import _native;print(json.dumps(dict(python=__import__("sys").executable,native_file=_native.__file__,version=prusie.__version__,backend=_native.backend_version())))'], env=settings, text=True)
    identity = json.loads(native)
    if identity['version'] != '0.3.0':
        raise ValueError('The manuscript experiment requires prusie 0.3.0')
    identity['native_sha256'] = sha(identity.pop('native_file'))
    write(output / 'native_identity.json', identity)
    return settings, python, rscript, probe, identity, read(output / 'R_identity.json')


def execute(row, order, config, output, runtime):
    import numpy as np
    settings, python, rscript, probe, identity, reference = runtime
    out = output / row['case_key']
    out.mkdir()
    case_path = resolve(row['case_path'], config['_manifest_base'])
    case = read(case_path)
    task = dict(case_path=str(case_path), case_id=case['case_id'], panel=row['panel'],
                parameters=case['parameters'], p=case['n_variants'], window=row.get('window'),
                n_threads=1, overwrite=False, kind=row.get('kind', 'primary'))
    request = dict(task=task, identity=identity,
                   placement={'cpus': sorted(os.sched_getaffinity(0))},
                   protocol_path=str(HERE / 'protocol.json'), protocol_sha256=sha(HERE / 'protocol.json'),
                   R_expected=reference, blas_probe=str(probe))
    if config.get('instrumented', False):
        request['fit_timing'] = dict(native_helper_directory=str(HERE), R_helper=str(HERE / 'fit_timer.R'),
                                      R_source_file=str(HERE / 'archive/R_source/susie.R.txt'))
    request['signature_sha256'] = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    write(out / 'request.json', request)
    for backend in order:
        argv = [rscript, '--vanilla', str(HERE / 'm2_reference.R')] if backend == 'R' else [python, '-B', str(HERE / 'm2_native.py')]
        argv += ['--request', str(out / 'request.json'), '--mode', 'time', '--output', str(out / (backend + '.json'))]
        with (out / (backend + '.log')).open('w') as log:
            subprocess.run(argv, env=settings, stdout=log, stderr=subprocess.STDOUT, check=True)
        metadata = read(out / (backend + '.json'))
        if metadata['status'] != 'ok' or len(metadata['samples']) != 4:
            raise ValueError(f'{backend} returned incomplete samples for {row["case_key"]}')
    records = {backend: read(out / (backend + '.json')) for backend in ('R', 'prusie')}
    shape = records['R']['alpha_shape']
    assert shape == records['prusie']['alpha_shape']
    arrays = {backend: np.fromfile(out / (backend + '.json.alpha.bin'), dtype='<f8').reshape(shape, order='F') for backend in records}
    difference = float(np.max(abs(arrays['R'] - arrays['prusie'])))
    summary = dict(case_key=row['case_key'],case_id=case['case_id'],panel=row['panel'],backend_order=order,
                   max_abs_alpha_difference=difference,parameters=case['parameters'],
                   alpha_shape=shape,complete_call_repeats=3,warmup_calls=1,
                   instrumentation='direct_same_call' if config.get('instrumented') else 'complete_call_only',
                   inputs={'case_sha256':sha(case_path), 'z_sha256':sha(case_path.parent / case['z_path']),
                           'variant_order_sha256':sha(case_path.parent / case['variant_ids_path'])})
    write(out / 'comparison.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1, help='Independent cases; backends remain sequential')
    parser.add_argument('--case-keys', help='Optional fixed comma-separated case keys; selection occurs after order freezing')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())['benchmark']
    config['_base'] = config_path.parent
    manifest = resolve(config['cases_manifest'], config_path.parent)
    config['_manifest_base'] = manifest.parent
    rows = read(manifest)
    if isinstance(rows, dict):
        rows = rows['rows']
    orders = balanced_orders(rows)
    output = resolve(config['output_dir'], config_path.parent)
    output.mkdir(parents=True, exist_ok=False)
    write(output / 'backend_orders.json', orders)
    runtime = environment(config, output)
    if args.case_keys:
        selection = args.case_keys.split(',')
        rows = [row for row in rows if row['case_key'] in selection]
        if len(rows) != len(selection):
            raise ValueError('Unknown or duplicated selected case key')
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        summaries = list(pool.map(lambda row: execute(row, orders[row['case_key']], config, output, runtime), rows))
    write(output / 'comparisons.json', summaries)
    print(f'Completed {len(summaries)} paired cases. Results: {output}')


if __name__ == '__main__':
    main()
