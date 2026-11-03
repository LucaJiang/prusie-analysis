"""Re-derive susieR 0.16.6 credible sets from saved alpha for both backends.

Usage: cs_membership.py SHARD NSHARDS OUT.jsonl [CASE_KEY ...]
"""
import csv, gzip, json, sys, time
from pathlib import Path
import numpy as np

from configuration import load
CFG = load()
TABLE = Path(CFG['primary_table'])
GWAS_DIR = Path(CFG['gwas_cases'])
EQTL_DIR = Path(CFG['eqtl_cases'])
BENCH = Path(CFG["saved_results"])

BLOCK = 256


def result_dir(r):
    for sub in ("results_v2", "results"):
        d = BENCH / sub / r["host"] / r["job_id"] / r["case_key"]
        if d.is_dir():
            return d
    raise FileNotFoundError(r["case_key"])


def load_alpha(d, backend):
    data = json.loads((d / f"{backend}.json").read_text())
    with gzip.open(d / f"{backend}.json.alpha.bin.gz", "rb") as f:
        alpha = np.frombuffer(f.read(), dtype="<f8").reshape(tuple(data["alpha_shape"]), order="F")
    return alpha, data


def coverage_set(x, coverage):
    o = np.argsort(-x, kind="stable")  # R order(decreasing=TRUE): radix, ties keep input order
    cs = np.cumsum(x[o].astype(np.longdouble))  # R cumsum accumulates in long double
    n = int(np.sum(cs < coverage)) + 1
    return np.sort(o[:n]), o, n


def min_abs_corr(R, pos, threshold):
    """Min |r| over upper triangle; early exit (returns value < threshold) once failing."""
    if len(pos) == 1:
        return 1.0, True
    m = 1.0
    for s in range(0, len(pos), BLOCK):
        rows = pos[s:s + BLOCK]
        sub = np.abs(np.asarray(R[rows][:, pos]))
        for i in range(len(rows)):
            sub[i, :s + i + 1] = np.inf
        m = min(m, float(sub.min()))
        if m < threshold:
            return m, False
    return m, True


def derive(alpha, R, coverage, min_corr):
    L, p = alpha.shape
    trimmed = [bool(np.all(alpha[l] == alpha[l, 0])) for l in range(L)]
    sets = [coverage_set(alpha[l], coverage) for l in range(L)]
    keys = [s[0].tobytes() for s in sets]
    seen, out = set(), {}
    for l in range(L):
        dup = keys[l] in seen
        seen.add(keys[l])
        if trimmed[l] or dup:
            continue
        pos = sets[l][0]
        purity, ok = min_abs_corr(R, pos, min_corr)
        if purity >= min_corr:
            out[l] = dict(members=pos, purity=purity, coverage=float(alpha[l, pos].sum()), order=sets[l][1], n=sets[l][2])
    return out, trimmed


def boundary_info(aR, aP, l, mR, mP):
    diff = np.union1d(np.setdiff1d(mR, mP), np.setdiff1d(mP, mR))
    return dict(component=l, size_R=len(mR), size_prusie=len(mP), symdiff=len(diff),
                symdiff_variants=diff.tolist()[:50],
                symdiff_alpha_R=aR[l, diff].tolist()[:50], symdiff_alpha_prusie=aP[l, diff].tolist()[:50],
                symdiff_max_abs_alpha_diff=float(np.max(np.abs(aR[l, diff] - aP[l, diff]))),
                row_max_abs_alpha_diff=float(np.max(np.abs(aR[l] - aP[l]))))


def process(r):
    d = result_dir(r)
    base = GWAS_DIR if r["panel"] == "GWAS" else EQTL_DIR
    meta = json.loads((base / r["case_id"] / "case.json").read_text())
    p = int(r["p"])
    assert meta["n_variants"] == p
    params = meta["parameters"]
    R = np.memmap(base / r["case_id"] / meta["R_path"], dtype="<f8", mode="r", shape=(p, p))
    row = dict(case_key=r["case_key"], case_id=r["case_id"], panel=r["panel"], p=p, result_dir=str(d))
    alphas, derived = {}, {}
    for b in ("R", "prusie"):
        a, data = load_alpha(d, b)
        assert a.shape[1] == p
        assert data["parameters"]["coverage"] == params["coverage"] and data["parameters"]["min_abs_corr"] == params["min_abs_corr"]
        s = data["samples"][-1]
        cs, trimmed = derive(a, R, params["coverage"], params["min_abs_corr"])
        alphas[b], derived[b] = a, cs
        row[f"{b}_niter"] = s["niter"]; row[f"{b}_converged"] = s["converged"]
        row[f"{b}_saved_cs_count"] = s["cs_count"]
        row[f"{b}_saved_cs_count_all_samples"] = sorted({x["cs_count"] for x in data["samples"]})
        row[f"{b}_derived_cs_count"] = len(cs)
        row[f"{b}_trimmed_rows"] = [l for l, t in enumerate(trimmed) if t]
        row[f"{b}_cs_components"] = sorted(cs)
        row[f"{b}_cs_sizes"] = [len(cs[l]["members"]) for l in sorted(cs)]
    aR, aP = alphas["R"], alphas["prusie"]
    row["max_abs_alpha_error"] = float(np.max(np.abs(aR - aP)))
    cR, cP = derived["R"], derived["prusie"]
    matched = sorted(set(cR) & set(cP))
    row["n_matched"] = len(matched)
    row["only_R_components"] = sorted(set(cR) - set(cP))
    row["only_prusie_components"] = sorted(set(cP) - set(cR))
    same, diffs = 0, []
    for l in matched:
        if np.array_equal(cR[l]["members"], cP[l]["members"]):
            same += 1
        else:
            diffs.append(boundary_info(aR, aP, l, cR[l]["members"], cP[l]["members"]))
    row["n_identical_members"] = same
    row["membership_differences"] = diffs
    row["max_abs_purity_diff"] = max([abs(cR[l]["purity"] - cP[l]["purity"]) for l in matched], default=0.0)
    row["max_abs_coverage_diff"] = max([abs(cR[l]["coverage"] - cP[l]["coverage"]) for l in matched], default=0.0)
    # Also compare raw coverage sets of every row (before purity), as a sensitivity check.
    raw_diff = []
    for l in range(aR.shape[0]):
        sR, sP = coverage_set(aR[l], params["coverage"])[0], coverage_set(aP[l], params["coverage"])[0]
        if not np.array_equal(sR, sP):
            raw_diff.append(dict(component=l, size_R=len(sR), size_prusie=len(sP),
                                 symdiff=int(len(np.setxor1d(sR, sP)))))
    row["raw_coverage_set_differences"] = raw_diff
    row["identical"] = (not row["only_R_components"] and not row["only_prusie_components"] and same == len(matched))
    return row


def main():
    shard, nshards, out = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
    only = set(sys.argv[4:])
    rows = sorted((r for r in csv.DictReader(TABLE.open(), delimiter="\t") if r["kind"] == "primary"),
                  key=lambda r: r["case_key"])
    done = set()
    if out.exists():
        done = {json.loads(l)["case_key"] for l in out.read_text().splitlines() if l.strip()}
    with out.open("a") as fh:
        for i, r in enumerate(rows):
            if (only and r["case_key"] not in only) or (not only and i % nshards != shard) or r["case_key"] in done:
                continue
            fh.write(json.dumps(process(r)) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
