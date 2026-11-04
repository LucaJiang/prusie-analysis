"""Summarize frozen paired calls without fitting or dropping any cases."""
import os
import csv
import gzip
import hashlib
import json
import math
import statistics
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parents[1]
DATA = PAPER / "03_results"
OUT = Path(os.environ.get("ANALYSIS_OUTPUT", DATA))
OUT.mkdir(parents=True, exist_ok=True)
EXPECTED = {"GWAS": 200, "eQTL": 14157, "windows": 6}


def describe(values):
    quantiles = np.quantile(values, [0, 0.25, 0.5, 0.75, 1])
    return dict(zip(("min", "q1", "median", "q3", "max"), map(float, quantiles)))


def main():
    records_path = DATA / "benchmark_records.jsonl.gz"
    with gzip.open(records_path, "rt") as handle:
        records = [json.loads(line) for line in handle]
    assert len(records) == sum(EXPECTED.values())
    assert len({record["case_key"] for record in records}) == len(records)
    rows = []
    for record in sorted(records, key=lambda entry: entry["case_key"]):
        row = {name: record[name] for name in ("case_key", "case_id", "panel", "kind", "p", "alpha_pass", "max_abs_alpha_error", "max_abs_pip_error")}
        row["host"] = record["host"]
        row["job_id"] = record["job_id"]
        row["ld_sha256"] = record["ld_sha256"]
        for source_backend, backend in (("R", "susieR"), ("prusie", "prusie")):
            result = record[source_backend]
            for metric, key in (("", "times"), ("fitting_", "fitting_times")):
                times = np.asarray(result[key], dtype=float)
                assert times.shape == (3,) and np.all(np.isfinite(times)) and np.all(times > 0)
                median = float(np.median(times))
                row[f"{backend}_{metric}seconds"] = median
                row[f"{backend}_{metric}minimum"] = float(times.min())
                row[f"{backend}_{metric}maximum"] = float(times.max())
                row[f"{backend}_{metric}q1"] = float(np.quantile(times, 0.25))
                row[f"{backend}_{metric}q3"] = float(np.quantile(times, 0.75))
                row[f"{backend}_{metric}repeat_iqr_over_median"] = float(np.diff(np.quantile(times, [0.25, 0.75]))[0] / median)
            row[f"{backend}_niter"] = result["niter"][1]
            row[f"{backend}_converged"] = result["converged"][1]
            peak = result["peak_rss_diagnostic_bytes"]
            assert isinstance(peak, int) and peak > 0
            row[f"{backend}_peak_rss_bytes"] = peak
        row["speedup"] = row["susieR_seconds"] / row["prusie_seconds"]
        row["fitting_speedup"] = row["susieR_fitting_seconds"] / row["prusie_fitting_seconds"]
        rows.append(row)
    summary = {
        "records_sha256": hashlib.sha256(records_path.read_bytes()).hexdigest(),
        "timing_source_protocol_sha256": json.loads((DATA / "benchmark_manuscript_protocol.json").read_text())["timing_source_protocol_sha256"],
        "retained_calls_per_backend": 3,
        "paired_ratio": "Per case: median of three susieR calls / median of three prusie calls; panel summary: median of these ratios",
        "total_software_seconds": "Sum of per-case median complete-call times, one representative analysis per case; not the sum of repeated calls or a measured parallel makespan",
        "clock_boundary": "Already loaded arrays enter public susie_rss through returned result, including software validation and preparation; excludes disk reads, LD reconstruction, imports, explicit GC, comparisons and serialization",
        "memory_source": "Backend timing-process high-water marks in benchmark_records.jsonl.gz: R.peak_rss_diagnostic_bytes and prusie.peak_rss_diagnostic_bytes",
        "memory_boundary": "One fresh process per case/backend; includes interpreter, libraries, loaded inputs, the warm-up and the three timed calls; captured after the final call before numerical comparison or serialization; not an isolated fitting allocation or an independent single-call memory experiment",
        "panels": {},
    }
    for panel, expected_count in EXPECTED.items():
        selected = [row for row in rows if row["panel"] == panel]
        panel_records = [record for record in records if record["panel"] == panel]
        assert len(selected) == expected_count
        entry = {"cases": expected_count, "metrics": {}}
        for metric, prefix in (("complete", ""), ("fitting", "fitting_")):
            ratios = np.asarray([row[f"{prefix}speedup"] for row in selected])
            entry["metrics"][metric] = {
                "paired_ratio": describe(ratios),
                "wins": int(np.sum(ratios > 1)),
                "ties": int(np.sum(ratios == 1)),
                "regressions": int(np.sum(ratios < 1)),
                "backends": {},
            }
            for source_backend, backend in (("R", "susieR"), ("prusie", "prusie")):
                times = [row[f"{backend}_{prefix}seconds"] for row in selected]
                entry["metrics"][metric]["backends"][backend] = {
                    "total_software_seconds": math.fsum(times),
                    "case_median_seconds": describe(times),
                    "repeat_iqr_over_median": describe([row[f"{backend}_{prefix}repeat_iqr_over_median"] for row in selected]),
                    "all_retained_seconds": math.fsum(value for record in panel_records for value in record[source_backend]["fitting_times" if prefix else "times"]),
                }
        entry["process_peak_rss_bytes"] = {backend: describe([row[f"{backend}_peak_rss_bytes"] for row in selected]) for backend in ("prusie", "susieR")}
        entry["alpha_passes"] = sum(row["alpha_pass"] for row in selected)
        entry["alpha_failures"] = expected_count - entry["alpha_passes"]
        entry["pip_exceedances"] = sum(row["max_abs_pip_error"] > 1e-5 for row in selected)
        entry["max_abs_alpha_error"] = max(row["max_abs_alpha_error"] for row in selected)
        entry["max_abs_pip_error"] = max(row["max_abs_pip_error"] for row in selected)
        entry["nonconverged"] = {backend: sum(not row[f"{backend}_converged"] for row in selected) for backend in ("susieR", "prusie")}
        entry["at_iteration_limit"] = {backend: sum(row[f"{backend}_niter"] == 100 for row in selected) for backend in ("susieR", "prusie")}
        entry["iteration_or_convergence_mismatches"] = [record["case_id"] for record in panel_records if record["R"]["niter"] != record["prusie"]["niter"] or record["R"]["converged"] != record["prusie"]["converged"]]
        entry["exceedances"] = [row for row in selected if not row["alpha_pass"]]
        entry["paired_peak_rss_ratio"] = describe([row["susieR_peak_rss_bytes"] / row["prusie_peak_rss_bytes"] for row in selected])
        summary["panels"][panel] = entry
    # The manuscript window table subsequently uses the dedicated ten-call R replay.
    # Panel summary above retains its original measurement scope.
    replay = json.loads((DATA / "window10000_retime.json").read_text())
    target = next(row for row in rows if row["case_key"] == "window_10000")
    for prefix, field in (("", "retained_complete_seconds"), ("fitting_", "retained_fitting_seconds")):
        vals = replay[field]
        median = statistics.median(vals)
        q1, _, q3 = statistics.quantiles(vals, n=4, method="inclusive")
        for suffix, value in (("seconds", median), ("minimum", min(vals)), ("maximum", max(vals)), ("q1", q1), ("q3", q3), ("repeat_iqr_over_median", (q3-q1)/median)):
            target[f"susieR_{prefix}{suffix}"] = f"{value:.15g}"
        target[f"{prefix}speedup"] = f"{median / float(target[f'prusie_{prefix}seconds']):.15g}"
    target["susieR_peak_rss_bytes"] = replay["peak_rss_bytes"]
    for name, selected in (("benchmark_primary.tsv", [row for row in rows if row["kind"] == "primary"]), ("benchmark_windows.tsv", [row for row in rows if row["panel"] == "windows"])):
        with (OUT / name).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected)
    (OUT / "benchmark_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for panel, entry in summary["panels"].items():
        print(panel, json.dumps({"paired_complete": entry["metrics"]["complete"]["paired_ratio"], "paired_fitting": entry["metrics"]["fitting"]["paired_ratio"], "totals_seconds": {backend: entry["metrics"]["complete"]["backends"][backend]["total_software_seconds"] for backend in ("prusie", "susieR")}, "alpha_failures": entry["alpha_failures"], "pip_exceedances": entry["pip_exceedances"], "nonconverged": entry["nonconverged"], "mismatches": entry["iteration_or_convergence_mismatches"]}))


if __name__ == "__main__":
    main()
