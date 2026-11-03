"""Small file helpers shared by the portable benchmark drivers."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path

COMMIT = '8e56a8e038e989856d106d9ca5175cc664fea9d2'
THREAD_VARS = ('OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS BLIS_NUM_THREADS '
               'VECLIB_MAXIMUM_THREADS NUMEXPR_NUM_THREADS RAYON_NUM_THREADS PRUSIE_NUM_THREADS').split()


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)
