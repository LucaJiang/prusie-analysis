#!/usr/bin/env python3
"""Rebuild full signed float64 LD with the manuscript's PLINK commands."""
from __future__ import annotations
from pathlib import Path
import argparse,csv,gzip,json,os,subprocess,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config
from _core.acquire import sha256 as sha, atomic_json as write
THREADS=("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS","RAYON_NUM_THREADS","VECLIB_MAXIMUM_THREADS","PRUSIE_NUM_THREADS")
PLINK_SHA256="47ce3069a1349128bdc5193aa8af33ec20dc99ece920843c80453c34eea1a062"
def count_missing(path, n_samples, n_variants):
    """Inspect BED missing codes only; no genotype correlation is calculated here."""
    width = (n_samples + 3) // 4
    counts = np.array([sum(((byte >> shift) & 3) == 1 for shift in (0, 2, 4, 6)) for byte in range(256)], dtype=np.uint8)
    tail = n_samples % 4
    tail_counts = np.array([sum(((byte >> shift) & 3) == 1 for shift in range(0, 2*tail, 2)) for byte in range(256)], dtype=np.uint8)
    missing = 0
    with Path(path).open("rb") as stream:
        if stream.read(3) != b"\x6c\x1b\x01":
            raise ValueError("BED must be SNP-major with canonical header")
        for start in range(0, n_variants, 8192):
            size = min(8192, n_variants-start)
            data = np.frombuffer(stream.read(size * width), dtype=np.uint8).reshape(size, width)
            if tail:
                missing += int(counts[data[:, :-1]].sum()) + int(tail_counts[data[:, -1]].sum())
            else:
                missing += int(counts[data].sum())
        if stream.read(1):
            raise ValueError("Unexpected BED trailing bytes")
    return missing

def validate_matrix(path, p, tolerances):
    if path.stat().st_size != 8*p*p:
        raise ValueError("LD byte count mismatch")
    matrix = np.memmap(path, mode="r", dtype="<f8", shape=(p, p))
    diag = symmetry = maximum = 0.0
    for start in range(0, p, 128):
        stop = min(p, start + 128)
        a = matrix[start:stop]
        if not np.isfinite(a).all():
            raise ValueError("Nonfinite signed LD")
        diag = max(diag, float(np.max(np.abs(a[np.arange(stop-start), np.arange(start, stop)]-1))))
        maximum = max(maximum, float(np.max(np.abs(a))))
        symmetry = max(symmetry, float(np.max(np.abs(a-matrix[:, start:stop].T))))
    del matrix
    if diag > tolerances["diagonal_atol"] or symmetry > tolerances["symmetry_atol"] or maximum > 1+tolerances["range_atol"]:
        raise ValueError(f"LD validation failed: diagonal={diag}, symmetry={symmetry}, range={maximum}")
    return {"finite": True, "diagonal_max_abs_error": diag, "symmetry_max_abs_error": symmetry,
            "max_abs_r": maximum, "tolerances": tolerances}

def reconstruct(v,case_id):
    base=Path(v['scratch'])/'cases'/case_id
    original_case=json.loads((base/'case.json').read_text())
    prefix=Path(v['reference_root'])/original_case['population']/('chr'+str(original_case['region'][0]))
    out=Path(v['scratch'])/'plink_ld'/case_id;out.mkdir(parents=True,exist_ok=True)
    plink=v['plink1']
    case = json.loads((base / "case.json").read_text())
    p = case["n_variants"]
    ids = (base / "variant_ids.txt").read_text().splitlines()[:p]
    if len(ids) != p or len(set(ids)) != p:
        raise ValueError("Frozen variant count/order/uniqueness mismatch")
    with gzip.open(base / "variants.tsv.gz", "rt") as stream:
        variants = list(csv.DictReader(stream, delimiter="\t"))[:p]
    if [x["variant_id"] for x in variants] != ids:
        raise ValueError("Frozen allele-table order differs from variant_ids")
    if any(x["counted_allele"] != x["alt"] for x in variants):
        raise ValueError("Unexpected original counted-allele contract")
    ref_bim_path = prefix.with_suffix(".bim")
    ref_bim = tuple(tuple(line.split()) for line in ref_bim_path.open())
    source_ids = []
    for frozen in variants:
        actual = ref_bim[int(frozen["reference_index"])]
        if int(actual[3]) != int(frozen["pos"]) or {actual[4], actual[5]} != {frozen["ref"], frozen["alt"]}:
            raise ValueError("Frozen reference_index does not identify the exact position/alleles")
        source_ids.append(actual[1])
    if len(set(source_ids)) != p:
        raise ValueError("Frozen reference indices do not identify unique PLINK variants")
    (out / "variant_ids.txt").write_text("\n".join(ids)+"\n")
    (out / "extract_ids.txt").write_text("\n".join(source_ids)+"\n")
    (out / "a1.tsv").write_text("".join(f"{source_id}\t{x['counted_allele']}\n" for source_id,x in zip(source_ids,variants)))
    z = np.fromfile(base / "z.bin", dtype="<f8", count=p)
    if z.size != p or not np.isfinite(z).all():
        raise ValueError("Frozen z values invalid")
    if not np.array_equal(z, np.array([float(x["z"]) for x in variants])):
        raise ValueError("Frozen z differs from aligned variant provenance")
    z.tofile(out / "z.bin")
    plink = str(Path(plink).resolve())
    plink_hash = sha(plink)
    if plink_hash != PLINK_SHA256:
        raise ValueError("PLINK executable differs from the frozen reference-preparation binary")
    version = subprocess.check_output([plink, "--version"], text=True).strip()
    if not ("PLINK v1.9" in version or "PLINK v1.90" in version):
        raise ValueError(f"This adapter requires PLINK1.9, got {version}")
    env = dict(os.environ)
    env.update({key: "1" for key in THREADS})
    subset, ld = out / "subset", out / "plink_ld"
    common = ["--keep-allele-order", "--make-founders", "--threads", "1", "--memory", "2048"]
    commands = [[plink, "--bfile", str(prefix), "--extract", str(out/"extract_ids.txt"),
                 "--a1-allele", str(out/"a1.tsv"), "2", "1", *common,
                 "--make-bed", "--out", str(subset)],
                [plink, "--bfile", str(subset), *common, "--r", "square", "bin", "--out", str(ld)]]
    for index, command in enumerate(commands):
        with (out / f"plink_{index}.stdout.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    original_fam = [x.split()[:2] for x in prefix.with_suffix(".fam").read_text().splitlines()]
    subset_fam = [x.split()[:2] for x in subset.with_suffix(".fam").read_text().splitlines()]
    if original_fam != subset_fam or len(subset_fam) != case["reference_n"]:
        raise ValueError("PLINK sample list/order changed")
    bim = [x.split() for x in subset.with_suffix(".bim").read_text().splitlines()]
    if len(bim) != p or {x[1] for x in bim} != set(source_ids):
        raise ValueError("PLINK variant set/count changed")
    by_id = {x[1]: (i, x) for i, x in enumerate(bim)}
    permutation, signs = [], []
    for frozen, source_id in zip(variants, source_ids):
        index, actual = by_id[source_id]
        if int(actual[3]) != int(frozen["pos"]) or {actual[4], actual[5]} != {frozen["ref"], frozen["alt"]}:
            raise ValueError("PLINK variant allele/position mismatch")
        permutation.append(index)
        signs.append(1 if actual[4] == frozen["counted_allele"] else -1)
    missing = count_missing(subset.with_suffix(".bed"), len(subset_fam), p)
    if missing:
        raise ValueError("Selected genotypes contain missing calls")
    raw = ld.with_suffix(".ld.bin")
    if raw.stat().st_size != 8*p*p:
        raise ValueError("PLINK did not emit a full float64 square matrix")
    if permutation == list(range(p)) and all(x == 1 for x in signs):
        raw.replace(out / "R.bin")
    else:
        matrix = np.memmap(raw, mode="r", dtype="<f8", shape=(p, p))
        order, flips = np.array(permutation), np.array(signs)
        with (out / "R.bin").open("wb") as stream:
            for i in range(0, p, 128):
                a = np.asarray(matrix[order[i:i+128]][:, order], dtype="<f8", order="C")
                a *= flips[i:i+128, None] * flips[None, :]
                a.tofile(stream)
        del matrix
        raw.unlink()
    validation = validate_matrix(out / "R.bin", p, dict(diagonal_atol=1e-10,symmetry_atol=1e-12,range_atol=1e-10))

    case.update(R_path='R.bin',z_path='z.bin',variant_ids_path='variant_ids.txt',
                variant_provenance_path=str(base/'variants.tsv.gz'),preparation_path=str(base/'preparation.json'))
    case.pop('R_column_major_path',None);case.pop('R_gzip_path',None)
    write(out/'case.json',case)
    record=dict(case_id=case_id,R_sha256=sha(out/'R.bin'),z_sha256=sha(out/'z.bin'),
                variant_order_sha256=sha(out/'variant_ids.txt'),validation=validation,
                plink_sha256=plink_hash,plink_version=version,commands=commands,
                missing_genotypes=missing,n_reference=len(subset_fam))
    write(out/'validation.json',record)
    for suffix in ('.bed','.bim','.fam'):subset.with_suffix(suffix).unlink()
    return record

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--case-id',required=True);a=ap.parse_args()
    print(json.dumps(reconstruct(config(a.config)['inputs'],a.case_id)))
if __name__=='__main__':main()
