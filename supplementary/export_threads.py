"""Export prusie 0.3.0 n_threads records (exclusive EPYC 7742) to a TSV table."""
import csv, hashlib, json, statistics, sys
from pathlib import Path
from configuration import load
CFG=load()
src=Path(CFG['thread_records'])
out=Path(CFG['thread_table'])
d = json.loads(src.read_text())
rows = []
for r in d["records"]:
    if not r["comparison"].startswith("C4_") or not r["applicable"]:
        continue
    rows.append({"case_id": r["case_id"], "panel": r["data_type"], "p": r["p"],
                 "n_threads": int(r["comparison"].split("_")[1][1:]),
                 "one_thread_median_seconds": statistics.median(r["parent_retained_seconds"]),
                 "threaded_median_seconds": statistics.median(r["candidate_retained_seconds"]),
                 "median_paired_round_ratio": r["ratio"],
                 "alpha_bitwise_equal_to_one_thread": r["alpha_bitwise_equal"],
                 "iterations_changed": r["niter_change"] != 0, "credible_sets_changed": r["cs_changed"]})
rows.sort(key=lambda r: (r["n_threads"], r["panel"], r["p"], r["case_id"]))
with out.open("w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), hashlib.sha256(src.read_bytes()).hexdigest(), d["backend_identities"]["scaling"]["V030_t4"]["wheel_sha256"])
