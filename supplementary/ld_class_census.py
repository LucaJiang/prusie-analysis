"""Count exact signed LD column classes for benchmark cases.

A class groups columns that are identical up to sign in float64. This is the
structure prusie's exact-class products can exploit; it does not depend on the
dispatch rule. Usage: ld_class_census.py PANEL SHARD NSHARDS OUT.jsonl
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from configuration import load
CFG = load()
TABLE = Path(CFG['primary_table'])
GWAS_DIR = Path(CFG['gwas_cases'])
EQTL_DIR = Path(CFG['eqtl_cases'])

EQTL_STRIDE = 20


def cases(panel):
    rows = [r for r in csv.DictReader(TABLE.open(), delimiter="\t")
            if r["kind"] == "primary" and r["panel"] == panel]
    rows.sort(key=lambda r: r["case_key"])
    if panel == "eQTL":
        rows = rows[::EQTL_STRIDE]
    base = GWAS_DIR if panel == "GWAS" else EQTL_DIR
    return [(r["case_key"], r["case_id"], int(r["p"]), base / r["case_id"]) for r in rows]


def census(path, p):
    meta = json.loads((path / "case.json").read_text())
    assert meta["n_variants"] == p
    R = np.memmap(path / meta["R_path"], dtype="<f8", mode="r", shape=(p, p))
    groups = {}
    for j in range(p):
        col = np.array(R[j])  # symmetric: row j equals column j
        first = col[np.flatnonzero(col)[0]] if col.any() else 1.0
        canon = (col if first > 0 else -col) + 0.0
        key = hashlib.blake2b(canon.tobytes(), digest_size=16).digest()
        reps = groups.setdefault(key, [])
        for rep in reps:
            if np.array_equal(rep[1], canon):
                rep[0] += 1
                break
        else:
            reps.append([1, canon])
    sizes = [s for reps in groups.values() for s, _ in reps]
    return {"G": len(sizes), "largest_class": max(sizes),
            "variants_in_shared_classes": int(sum(s for s in sizes if s > 1))}


def main():
    panel, shard, nshards, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4])
    done = set()
    if out.exists():
        done = {json.loads(line)["case_key"] for line in out.read_text().splitlines() if line}
    with out.open("a") as fh:
        for i, (key, case_id, p, path) in enumerate(cases(panel)):
            if i % nshards != shard or key in done:
                continue
            row = {"case_key": key, "case_id": case_id, "panel": panel, "p": p, **census(path, p)}
            row["G_over_p"] = row["G"] / p
            fh.write(json.dumps(row) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
