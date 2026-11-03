"""Data summaries for the 0.3.0 benchmark figures.

Reads only the frozen current-release artefacts:
    docs/assets/benchmark/burgundy_records.jsonl.gz
    docs/assets/benchmark/burgundy_summary.json

Pure functions returning plain dicts and NumPy arrays. No fitting, no
re-measurement, no historical-campaign data.
"""

from __future__ import annotations

import gzip
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def load_records(path):
    """Load `burgundy_records.jsonl.gz` as a list of dicts."""
    with gzip.open(Path(path), "rt") as stream:
        return [json.loads(line) for line in stream]


def load_summary(path):
    return json.loads(Path(path).read_text())


def panel_rows(records, panel):
    """Primary (non-window) measurement rows for one panel."""
    return [r for r in records if r["panel"] == panel and r["kind"] == "primary"]


def _call_medians(backend_record, metric):
    key = {"complete": "times", "fitting": "fitting_times"}[metric]
    values = [v for v in backend_record[key] if _finite(v) and v > 0]
    return float(np.median(values)) if values else math.nan


def per_case(records, panel, metric="complete"):
    """One row per case: p, prusie median, susieR median, speedup, convergence.

    `metric` is "complete" or "fitting".
    """
    rows = []
    for r in panel_rows(records, panel):
        prusie_s = _call_medians(r["prusie"], metric)
        susie_s = _call_medians(r["R"], metric)
        p = int(r["p"])
        iters_prusie = [n for n in r["prusie"].get("niter", []) if isinstance(n, int)]
        iters_susie = [n for n in r["R"].get("niter", []) if isinstance(n, int)]
        median_iter_prusie = float(np.median(iters_prusie)) if iters_prusie else math.nan
        median_iter_susie = float(np.median(iters_susie)) if iters_susie else math.nan
        rows.append(
            {
                "case_id": r["case_id"],
                "panel": panel,
                "p": p,
                "prusie_seconds": prusie_s,
                "susieR_seconds": susie_s,
                "speedup": (susie_s / prusie_s) if (prusie_s > 0 and susie_s > 0) else math.nan,
                "max_abs_alpha_error": r.get("max_abs_alpha_error"),
                "alpha_pass": r.get("alpha_pass"),
                "both_nonconverged_at_limit": r.get("both_nonconverged_at_limit", False),
                "numerical_disposition": r.get("numerical_disposition"),
                "median_iter_prusie": median_iter_prusie,
                "median_iter_susie": median_iter_susie,
                "max_median_iter": float(np.nanmax([median_iter_prusie, median_iter_susie])),
            }
        )
    return rows


def alpha_summary(records, panel):
    """Alpha-only agreement summary for a panel: counts and exceedance list."""
    rows = panel_rows(records, panel)
    diffs = [r.get("max_abs_alpha_error") for r in rows]
    exact_zero = sum(1 for d in diffs if d == 0)
    nonzero = [d for d in diffs if _finite(d) and d > 0]
    exceedances = [
        {
            "case_id": r["case_id"],
            "max_abs_alpha_error": r["max_abs_alpha_error"],
            "disposition": r.get("numerical_disposition"),
            "both_nonconverged_at_limit": r.get("both_nonconverged_at_limit", False),
            "p": int(r["p"]),
        }
        for r in rows
        if r.get("alpha_pass") is False
    ]
    both_nc = sum(1 for r in rows if r.get("both_nonconverged_at_limit"))
    exceedances_with_converged = sum(
        1 for e in exceedances if not e["both_nonconverged_at_limit"]
    )
    return {
        "cases": len(rows),
        "max_abs_alpha_error": max((d for d in diffs if _finite(d)), default=math.nan),
        "exact_zero_cases": exact_zero,
        "nonzero_cases": len(nonzero),
        "exceedances": exceedances,
        "exceedance_count": len(exceedances),
        "exceedances_with_at_least_one_converged": exceedances_with_converged,
        "both_nonconverged_at_limit": both_nc,
    }


def _log_bin_edges(values, n_bins):
    values = np.asarray([v for v in values if _finite(v) and v > 0], dtype=float)
    if values.size == 0:
        return np.array([1.0, 10.0])
    lo, hi = float(values.min()), float(values.max())
    if lo == hi:
        hi = lo * 1.01 + 1.0
    return np.geomspace(lo, hi, n_bins + 1)


def binned_runtime(records, panel, metric="complete", n_bins=None):
    """Log-p-binned paired runtime summary.

    Returns a list of bin dicts with p centre, counts, prusie/susieR median
    and IQR, and the aggregate speedup Σ susieR / Σ prusie within the bin.
    """
    if n_bins is None:
        n_bins = 5 if panel == "GWAS" else 12
    rows = [r for r in per_case(records, panel, metric) if r["prusie_seconds"] > 0 and r["susieR_seconds"] > 0]
    if not rows:
        return []
    edges = _log_bin_edges([r["p"] for r in rows], n_bins)
    bins = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        sel = [r for r in rows if (lo <= r["p"] < hi) or (i == n_bins - 1 and r["p"] == hi)]
        if not sel:
            continue
        prusie = np.array([r["prusie_seconds"] for r in sel])
        susie = np.array([r["susieR_seconds"] for r in sel])
        bins.append(
            {
                "p_lo": float(lo),
                "p_hi": float(hi),
                "p_centre": float(math.sqrt(lo * hi)),
                "n": len(sel),
                "prusie_median": float(np.median(prusie)),
                "prusie_q1": float(np.quantile(prusie, 0.25)),
                "prusie_q3": float(np.quantile(prusie, 0.75)),
                "susie_median": float(np.median(susie)),
                "susie_q1": float(np.quantile(susie, 0.25)),
                "susie_q3": float(np.quantile(susie, 0.75)),
                "aggregate_speedup": float(susie.sum() / prusie.sum()),
            }
        )
    return bins


def aggregate_speedup(records, panel, metric="complete"):
    """Σ susieR retained seconds / Σ prusie retained seconds, panel-wide.

    Uses all retained calls (not just medians) to match the headline
    definition in `burgundy_summary.json`.
    """
    rows = panel_rows(records, panel)
    key = {"complete": "times", "fitting": "fitting_times"}[metric]
    prusie_total = 0.0
    susie_total = 0.0
    paired = 0
    for r in rows:
        p_times = [v for v in r["prusie"][key] if _finite(v) and v > 0]
        s_times = [v for v in r["R"][key] if _finite(v) and v > 0]
        if not p_times or not s_times:
            continue
        paired += 1
        prusie_total += sum(p_times)
        susie_total += sum(s_times)
    return {
        "panel": panel,
        "metric": metric,
        "paired_cases": paired,
        "prusie_total_seconds": prusie_total,
        "susie_total_seconds": susie_total,
        "aggregate_speedup": (susie_total / prusie_total) if prusie_total > 0 else math.nan,
    }


def convergence_breakdown(records, panel):
    """Count iteration-limit outcomes by backend for a panel."""
    rows = panel_rows(records, panel)
    counts = Counter()
    for r in rows:
        prusie_hit = any(n is not None and n >= 100 and c is False for n, c in zip(r["prusie"].get("niter", []), r["prusie"].get("converged", []), strict=False))
        susie_hit = any(n is not None and n >= 100 and c is False for n, c in zip(r["R"].get("niter", []), r["R"].get("converged", []), strict=False))
        if prusie_hit and susie_hit:
            counts["both_at_limit"] += 1
        elif prusie_hit:
            counts["prusie_only_at_limit"] += 1
        elif susie_hit:
            counts["susie_only_at_limit"] += 1
        else:
            counts["converged"] += 1
    return dict(counts)
