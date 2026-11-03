"""Derived analyses: per-sweep window times, memory by region size, per-sweep
scaling exponents, the exact-column-class (G/p) comparison and thread scaling.

Reads data/benchmark_primary.tsv, data/benchmark_windows.tsv,
data/ld_class_census.tsv and data/thread_scaling.tsv; writes data/derived_analyses.json.
"""
import os
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

DATA = Path(os.environ.get("ANALYSIS_INPUT", Path(__file__).resolve().parents[1] / "03_results"))
OUT = Path(os.environ.get("ANALYSIS_OUTPUT", DATA))
MIB = 2**20
MEMORY_BINS = [0, 1000, 2000, 5000, 10000, None]


def read(name):
    with (DATA / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def windows():
    rows = sorted(read("benchmark_windows.tsv"), key=lambda r: int(r["p"]))
    out = []
    for r in rows:
        sweeps = int(r["prusie_niter"])
        assert sweeps == int(r["susieR_niter"])
        prusie = float(r["prusie_fitting_seconds"]) / sweeps
        susier = float(r["susieR_fitting_seconds"]) / sweeps
        out.append({"p": int(r["p"]), "sweeps": sweeps, "prusie_ms_per_sweep": 1e3 * prusie,
                    "susieR_ms_per_sweep": 1e3 * susier, "per_sweep_ratio": susier / prusie})
    p = np.log([w["p"] for w in out])
    slopes = {b: float(np.polyfit(p, np.log([w[f"{b}_ms_per_sweep"] for w in out]), 1)[0])
              for b in ("prusie", "susieR")}
    return {"rows": out, "log_log_slope": slopes}


def memory(primary):
    out = {}
    for panel in ("GWAS", "eQTL"):
        rows = [r for r in primary if r["panel"] == panel]
        p = np.array([int(r["p"]) for r in rows], float)
        ld = 8 * p**2
        prusie = np.array([int(r["prusie_peak_rss_bytes"]) for r in rows], float)
        susier = np.array([int(r["susieR_peak_rss_bytes"]) for r in rows], float)
        bins = []
        for lo, hi in zip(MEMORY_BINS, MEMORY_BINS[1:]):
            m = (p >= lo) & (p < (hi or np.inf))
            if not m.any():
                continue
            bins.append({"p_from": lo, "p_below": hi, "cases": int(m.sum()),
                         "ld_MiB": float(np.median(ld[m]) / MIB),
                         "prusie_MiB": float(np.median(prusie[m]) / MIB),
                         "susieR_MiB": float(np.median(susier[m]) / MIB),
                         "prusie_over_ld": float(np.median(prusie[m] / ld[m])),
                         "susieR_over_ld": float(np.median(susier[m] / ld[m])),
                         "paired_ratio": float(np.median(susier[m] / prusie[m]))})
        large = p >= 5000
        out[panel] = {"bins": bins, "p_at_least_5000": {
            "cases": int(large.sum()),
            "prusie_over_ld": float(np.median(prusie[large] / ld[large])),
            "susieR_over_ld": float(np.median(susier[large] / ld[large]))}}
    return out


def sweep_exponents(primary):
    out = {}
    for panel in ("GWAS", "eQTL"):
        rows = [r for r in primary if r["panel"] == panel and int(r["p"]) >= 1000 and int(r["prusie_niter"]) >= 3]
        p = np.log([int(r["p"]) for r in rows])
        sweeps = np.array([int(r["prusie_niter"]) for r in rows], float)
        out[panel] = {"cases": len(rows)}
        for b in ("prusie", "susieR"):
            t = np.array([float(r[f"{b}_fitting_seconds"]) for r in rows]) / sweeps
            out[panel][b] = float(np.polyfit(p, np.log(t), 1)[0])
    return out


def column_classes(primary):
    by_key = {r["case_key"]: r for r in primary}
    census = read("ld_class_census.tsv")
    out = {}
    for panel in ("GWAS", "eQTL"):
        rows = [c for c in census if c["panel"] == panel]
        for c in rows:
            assert int(c["p"]) == int(by_key[c["case_key"]]["p"])
        g = np.array([int(c["G"]) / int(c["p"]) for c in rows])
        p = np.array([int(c["p"]) for c in rows])
        ratio = np.array([float(by_key[c["case_key"]]["fitting_speedup"]) for c in rows])
        m = (g >= 0.35) & (g < 0.45) & (p >= 2000) & (p <= 6000)
        out[panel] = {"cases": len(rows), "median_G_over_p": float(np.median(g)),
                      "spearman_p_G_over_p": float(spearmanr(p, g).statistic),
                      "matched_bin_cases": int(m.sum()),
                      "matched_bin_median_fitting_ratio": float(np.median(ratio[m]))}
        if panel == "eQTL":
            order = np.argsort(g, kind="stable")
            out[panel]["tertiles"] = [
                {"cases": len(t), "G_over_p_min": float(g[t].min()), "G_over_p_max": float(g[t].max()),
                 "median_fitting_ratio": float(np.median(ratio[t]))}
                for t in np.array_split(order, 3)]
            high = g >= 0.7
            out[panel]["G_over_p_at_least_0.7"] = {
                "cases": int(high.sum()), "median_fitting_ratio": float(np.median(ratio[high])),
                "fitting_ratio_q1_q3": [float(q) for q in np.quantile(ratio[high], [0.25, 0.75])]}
    return out


def threads():
    rows = read("thread_scaling.tsv")
    out = {"cases": len({r["case_id"] for r in rows}), "settings": []}
    for k in sorted({int(r["n_threads"]) for r in rows}):
        sel = [r for r in rows if int(r["n_threads"]) == k]
        speed = np.array([1 / float(r["median_paired_round_ratio"]) for r in sel])
        large = np.array([int(r["p"]) >= 10000 for r in sel])
        big = max(sel, key=lambda r: int(r["p"]))
        out["settings"].append({
            "n_threads": k, "cases": len(sel),
            "median_speedup": float(np.median(speed)),
            "speedup_q1_q3": [float(q) for q in np.quantile(speed, [0.25, 0.75])],
            "slower_than_one_thread": int((speed < 1).sum()),
            "p_at_least_10000": {"cases": int(large.sum()), "median_speedup": float(np.median(speed[large]))},
            "largest_case": {"case_id": big["case_id"], "p": int(big["p"]),
                             "one_thread_seconds": float(big["one_thread_median_seconds"]),
                             "threaded_seconds": float(big["threaded_median_seconds"])},
            "alpha_bitwise_equal_cases": sum(r["alpha_bitwise_equal_to_one_thread"] == "True" for r in sel),
            "iteration_or_credible_set_changes": sum(r["iterations_changed"] == "True" or r["credible_sets_changed"] == "True" for r in sel)})
    return out


def main():
    primary = [r for r in read("benchmark_primary.tsv") if r["kind"] == "primary"]
    result = {"windows_per_sweep": windows(), "memory_by_p": memory(primary),
              "per_sweep_exponents": sweep_exponents(primary),
              "exact_column_classes": column_classes(primary), "threads": threads()}
    (OUT / "derived_analyses.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
