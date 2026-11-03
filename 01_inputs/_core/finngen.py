"""FinnGen R13 T2D: released-P locus selection and whole-window HLA exclusion."""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import csv, gzip, json, math, time
from .acquire import atomic_json, sha256

FIELDS = {"#chrom", "pos", "ref", "alt", "beta", "sebeta", "pval", "mlogp", "af_alt"}

P_THRESHOLD = 5e-8

SELECTION_POLICY = "published_P_underflow_mlogp_inclusive_finalHLA_v2"

def hla_overlap(chrom: int | str, start: int, end: int) -> bool:
    return str(chrom).removeprefix("chr") == "6" and start <= 34_000_000 and end >= 25_000_000

def significant_logp(p: float, mlogp: float | None) -> float:
    if not math.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("Invalid published P value")
    if p > 0:
        return -math.log10(p)
    if mlogp is not None and math.isfinite(mlogp) and mlogp >= 0:
        return mlogp
    raise ValueError("Underflowed P without finite released mlogp")

def merge_intervals(regions: list[dict]) -> list[dict]:
    merged = []
    for region in sorted(regions, key=lambda r: (int(r["chrom"]), r["start"], r["end"])):
        r = dict(region)
        if merged and merged[-1]["chrom"] == r["chrom"] and r["start"] <= merged[-1]["end"] + 1:
            merged[-1]["end"] = max(merged[-1]["end"], r["end"])
        else:
            merged.append({k: r[k] for k in ("chrom", "start", "end")})
    return merged

def define_regions(signals: list[dict], lengths: dict, half_window: float = 1_500_000,
                   maximum_width: int = 10_000_000, shrink: float = 0.9) -> dict:
    """Greedy strongest leads, overlap merging, recursive oversized-region shrink.

    Signal records contain chrom, pos, mlogp, ref, alt. No reference coverage,
    fine-mapping output, or association-strength subset beyond P<5e-8 is used.
    """
    lengths = {str(k).removeprefix("chr"): int(v) for k, v in lengths.items()}
    ordered = sorted(signals, key=lambda r: (-r["mlogp"], int(r["chrom"]), r["pos"], r["ref"], r["alt"]))
    for s in ordered:
        if not 1 <= s["pos"] <= lengths[str(s["chrom"])]:
            raise ValueError("Signal outside chromosome bounds")
    changes, leads, initial = [], [], []

    def recurse(rows, radius, depth):
        if depth > 250:
            raise RuntimeError("Window shrink failed to terminate")
        remaining = list(rows)
        windows = []
        while remaining:
            lead = remaining[0]
            chrom = str(lead["chrom"])
            start = max(1, int(lead["pos"] - radius))
            end = min(lengths[chrom], int(lead["pos"] + radius))
            window = {"chrom": chrom, "start": start, "end": end}
            windows.append(window)
            leads.append({**window, "lead_pos": lead["pos"], "lead_ref": lead["ref"],
                          "lead_alt": lead["alt"], "lead_mlogp": lead["mlogp"], "depth": depth,
                          "half_window": radius})
            remaining = [s for s in remaining if not (str(s["chrom"]) == chrom and start <= s["pos"] <= end)]
        merged = merge_intervals(windows)
        if depth == 0:
            initial.extend(merged)
        small = [r for r in merged if r["end"] - r["start"] + 1 <= maximum_width]
        large = [r for r in merged if r["end"] - r["start"] + 1 > maximum_width]
        for r in large:
            changes.append({**r, "depth": depth, "old_half_window": radius,
                            "new_half_window": radius * shrink, "reason": "merged_width_exceeds_10Mb"})
        if not large:
            return small
        subset = [s for s in rows if any(str(s["chrom"]) == r["chrom"] and r["start"] <= s["pos"] <= r["end"] for r in large)]
        out = merge_intervals(small + recurse(subset, radius * shrink, depth + 1))
        # A newly selected edge lead can extend outside a previous oversized
        # interval and meet a retained neighbour. Recheck the actual final union.
        if any(r["end"] - r["start"] + 1 > maximum_width for r in out):
            changes.append({"reason": "recursive_union_exceeds_10Mb", "depth": depth,
                            "old_half_window": radius, "new_half_window": radius * shrink})
            return recurse(rows, radius * shrink, depth + 1)
        return out

    final = recurse(ordered, half_window, 0) if ordered else []
    for r in final:
        r.update(case_id=f"T2D_chr{r['chrom']}_{r['start']}_{r['end']}",
                 bed_start=r["start"] - 1, bed_end=r["end"], hla_overlap=hla_overlap(r["chrom"], r["start"], r["end"]))
        r["status"] = "EXCLUDED" if r["hla_overlap"] else "ELIGIBLE"
        r["reason"] = "HLA_overlap" if r["hla_overlap"] else ""
    coverage = Counter()
    for s in ordered:
        matches = [r for r in final if str(s["chrom"]) == r["chrom"] and r["start"] <= s["pos"] <= r["end"]]
        if len(matches) != 1:
            raise RuntimeError("Significant signal lacks exactly one final region")
        r = matches[0]
        if not r["hla_overlap"]:
            coverage["retained"] += 1
        elif hla_overlap(s["chrom"], s["pos"], s["pos"]):
            coverage["inside_HLA"] += 1
        else:
            coverage["non_HLA_flank_lost_by_whole_window_exclusion"] += 1
    return {"initial_merged": initial, "lead_windows": leads, "changes": changes,
            "final": final, "signal_coverage": dict(coverage)}

def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)

def census(summary: Path, lengths: dict, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    receipt_path = output / "gwas_census.json"
    source_hash = sha256(summary)
    if receipt_path.exists():
        old = json.loads(receipt_path.read_text())
        if (old.get("summary_sha256") == source_hash and old.get("chromosome_lengths") == lengths
                and old.get("selection_policy") == SELECTION_POLICY):
            for name, digest in old["ledger_hashes"].items():
                if not (output / name).exists() or sha256(output / name) != digest:
                    raise RuntimeError("Frozen GWAS ledger missing or changed; investigate before reuse")
            return old
        raise RuntimeError("Frozen GWAS roster identity differs; preserve it and use a new output directory")
    chrom_counts, invalid = Counter(), Counter()
    signals = []
    with gzip.open(summary, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not FIELDS.issubset(reader.fieldnames):
            raise ValueError(f"Missing R13 columns: {FIELDS - set(reader.fieldnames)}")
        header = reader.fieldnames
        for row in reader:
            chrom = row["#chrom"].removeprefix("chr")
            chrom_counts[chrom] += 1
            if chrom not in {str(x) for x in range(1, 23)}:
                continue
            try:
                pos = int(row["pos"])
                if not 1 <= pos <= int(lengths[chrom]):
                    raise ValueError("coordinate")
                logp = significant_logp(float(row["pval"]), float(row["mlogp"]))
            except (ValueError, KeyError):
                invalid["invalid_position_or_pvalue"] += 1
                continue
            if logp > -math.log10(P_THRESHOLD):
                signals.append({"chrom": chrom, "pos": pos, "ref": row["ref"], "alt": row["alt"],
                                "pval": row["pval"], "mlogp": logp, "released_mlogp": row["mlogp"],
                                "original_id": row.get("rsids", "")})
    missing = [str(i) for i in range(1, 23) if not chrom_counts[str(i)]]
    if missing or invalid:
        raise ValueError(f"Full-summary census requires investigation: missing={missing}, invalid={dict(invalid)}")
    region_data = define_regions(signals, lengths)
    for name, rows in [("gwas_regions_initial", region_data["initial_merged"]),
                       ("gwas_lead_windows", region_data["lead_windows"]),
                       ("gwas_region_changes", region_data["changes"]),
                       ("gwas_regions", region_data["final"]), ("gwas_signals", signals)]:
        write_tsv(output / f"{name}.tsv", rows)
    frozen = {"release": "R13", "endpoint": "T2D", "summary_path": str(summary), "summary_sha256": source_hash,
              "chromosome_lengths": lengths, "header": header, "chromosome_rows": dict(chrom_counts),
              "significant_signals": len(signals), "initial_regions": len(region_data["initial_merged"]),
              "final_regions_before_hla": len(region_data["final"]),
              "HLA_overlap_regions": sum(r["hla_overlap"] for r in region_data["final"]),
              "eligible_regions": sum(not r["hla_overlap"] for r in region_data["final"]),
              "signal_coverage": region_data["signal_coverage"], "frozen_epoch": time.time(),
              "selection_independent_of_reference_and_fits": True, "selection_policy": SELECTION_POLICY,
              "ledger_hashes": {name: sha256(output / name) for name in
                                ("gwas_regions_initial.tsv", "gwas_lead_windows.tsv", "gwas_region_changes.tsv", "gwas_regions.tsv", "gwas_signals.tsv")}}
    atomic_json(receipt_path, frozen)
    return frozen
