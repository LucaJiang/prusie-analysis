"""Reference sample selection, source indexing and dosage checks."""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import hashlib, io, json, os, time
import numpy as np
from .acquire import atomic_json, sha256

RELEASE = "IGSR NYGC 20220422_3202_phased_SNV_INDEL_SV"

ALLELES = frozenset("ACGT")

AUTOSOMES = tuple(str(i) for i in range(1, 23))

SHORT_ALLELE_MAX = 50

def _table(path: Path) -> list[dict]:
    with path.open() as handle:
        header = next(handle).lstrip("#").split()
        return [dict(zip(header, line.split(), strict=True)) for line in handle if line.strip()]

def select_samples(source_psam: Path, pedigree: Path, corrected_psam: Path,
                   exclusions: Path, kinship: Path, output: Path) -> dict:
    """Use the official population roster and release-specific KING exclusions.

    All metadata joins are exact sample-ID joins. The precomputed removal list
    is checked against the complete published first/second-degree pair table;
    checking founders is deliberately insufficient.
    """
    output.mkdir(parents=True, exist_ok=True)
    samples = _table(source_psam)
    official = {r["SampleID"]: r for r in _table(pedigree)}
    corrected_rows = _table(corrected_psam)
    corrected = {r["IID"]: r for r in corrected_rows}
    ids = [r["IID"] for r in samples]
    if len(ids) != len(set(ids)) or set(ids) != set(official) or set(ids) != set(corrected):
        raise ValueError("Sample IDs are duplicated or disagree between source and official metadata")
    if ids != [r["IID"] for r in corrected_rows]:
        raise ValueError("PGEN sample order differs from the release's corrected PSAM")
    for row in samples:
        truth = official[row["IID"]]
        if row["Population"] != truth["Population"] or row["SuperPop"] != truth["Superpopulation"]:
            raise ValueError(f"Population mismatch: {row['IID']}")
    removed = {r["IID"] for r in _table(exclusions)}
    if not removed <= set(ids):
        raise ValueError("Kinship exclusion list contains unknown samples")
    pairs = _table(kinship)
    retained = set(ids) - removed
    unresolved = [(r["IID1"], r["IID2"], r["KINSHIP"]) for r in pairs
                  if float(r["KINSHIP"]) > 0.0884 and r["IID1"] in retained and r["IID2"] in retained]
    if unresolved:
        raise ValueError(f"Retained close-relative pairs: {unresolved}")
    # Verify listed pedigree edges too, including shared-parent sibships.
    parent_map = {r["IID"]: {r["PAT"], r["MAT"]} - {"0", "NA"} for r in corrected_rows}
    pedigree_edges = set()
    for iid, parents in parent_map.items():
        for parent in parents & set(ids):
            pedigree_edges.add(tuple(sorted((iid, parent))))
    families: dict[str, list[str]] = {}
    for iid, parents in parent_map.items():
        for parent in parents:
            families.setdefault(parent, []).append(iid)
    for siblings in families.values():
        for i, left in enumerate(siblings):
            for right in siblings[i + 1:]:
                pedigree_edges.add(tuple(sorted((left, right))))
    if any(set(edge) <= retained for edge in pedigree_edges):
        raise ValueError("Precomputed kinship exclusions leave an explicit pedigree relationship")
    report = {
        "release": RELEASE, "source_sample_count": len(ids), "excluded_close_relatives": len(removed),
        "all_retained_count": len(retained), "king_cutoff": 0.0884,
        "rule": "Official IGSR population labels; remove the PLINK resource's release-specific first- and second-degree KING exclusion list; preserve source sample order",
        "reported_related_pairs": len(pairs), "retained_related_pairs": len(unresolved),
        "explicit_pedigree_edges_checked": len(pedigree_edges),
        "metadata_sha256": {str(p): sha256(p) for p in (source_psam, pedigree, corrected_psam, exclusions, kinship)},
        "populations": {},
    }
    with (output / "sample_selection.tsv").open("w") as ledger:
        ledger.write("source_sample_index\tsample_id\tpopulation\tsuperpopulation\tselection\n")
        for idx, row in enumerate(samples):
            ledger.write(f"{idx}\t{row['IID']}\t{row['Population']}\t{row['SuperPop']}\t"
                         f"{'excluded_KING_degree1_or_2' if row['IID'] in removed else 'retained'}\n")
    for pop in ("FIN", "EUR"):
        chosen = [(idx, row) for idx, row in enumerate(samples)
                  if row["IID"] in retained and row["Population" if pop == "FIN" else "SuperPop"] == pop]
        with (output / f"{pop}.samples.tsv").open("w") as handle:
            handle.write("source_sample_index\tsample_id\tpopulation\tsuperpopulation\n")
            for idx, row in chosen:
                handle.write(f"{idx}\t{row['IID']}\t{row['Population']}\t{row['SuperPop']}\n")
        with (output / f"{pop}.psam").open("w") as handle:
            handle.write("#IID\tPAT\tMAT\tSEX\tSuperPop\tPopulation\n")
            for _, row in chosen:
                corr = corrected[row["IID"]]
                handle.write("\t".join(corr[k] for k in ("IID", "PAT", "MAT", "SEX", "SuperPop", "Population")) + "\n")
        report["populations"][pop] = {
            "n": len(chosen), "sample_list": str(output / f"{pop}.samples.tsv"),
            "sample_list_sha256": sha256(output / f"{pop}.samples.tsv"),
            "excluded_ids": [r["IID"] for r in samples if r["IID"] in removed
                             and r["Population" if pop == "FIN" else "SuperPop"] == pop],
        }
    atomic_json(output / "sample_selection.json", report)
    return report

def variant_class(ref: str, alt: str) -> str:
    if "," in alt:
        return "multiallelic"
    if not ref or not alt or not set(ref + alt) <= ALLELES:
        return "symbolic_or_non_ACGT"
    if ref == alt:
        return "identical_alleles"
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt) and max(len(ref), len(alt)) <= SHORT_ALLELE_MAX:
        return "short_indel"
    return "long_or_complex_allele"

def normalize_variant(pos: int, ref: str, alt: str, sequence: str) -> tuple[int, str, str]:
    """Return a left-normalized biallelic variant, checking 1-based REF first.

    The source is already left normalized, but this independently verifies and
    canonicalizes each retained short indel. No strand/swap repair is applied to
    a source REF mismatch. The caller retains original fields in the QC ledger.
    """
    if pos < 1 or pos - 1 + len(ref) > len(sequence):
        raise ValueError("coordinate_out_of_bounds")
    if sequence[pos - 1:pos - 1 + len(ref)].upper() != ref:
        raise ValueError("REF_mismatch")
    if variant_class(ref, alt) not in ("SNV", "short_indel"):
        raise ValueError("unsupported_alleles")
    # Remove common suffixes. If an allele would become empty, extend one base
    # left first; at chromosome start preserve the required VCF anchor.
    while ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
        if min(len(ref), len(alt)) == 1:
            if pos == 1:
                break
            previous = sequence[pos - 2].upper()
            if previous not in ALLELES:
                break
            ref, alt, pos = previous + ref[:-1], previous + alt[:-1], pos - 1
        else:
            ref, alt = ref[:-1], alt[:-1]
    while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    if sequence[pos - 1:pos - 1 + len(ref)].upper() != ref:
        raise ValueError("normalized_REF_mismatch")
    return pos, ref, alt

def census_pvar(source: Path, output: Path, status_path: Path | None = None) -> dict:
    """Stream the complete PVAR once; index every autosomal row without filtering."""
    import pysam

    output.mkdir(parents=True, exist_ok=True)
    receipt = output / "source_census.json"
    stat = source.stat()
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old["source_size"] == stat.st_size and old["source_mtime_ns"] == stat.st_mtime_ns:
            if all((output / f"chr{c}.tsv.gz.tbi").exists() for c in AUTOSOMES):
                return old
    handles = {}
    counts: dict[str, dict] = {}
    nonautosomal = Counter()
    digest = hashlib.sha256()
    idx = 0
    headers = []
    first_rows = {}
    started = time.time()
    try:
        with source.open("rb") as handle:
            for line in handle:
                digest.update(line)
                if line.startswith(b"#"):
                    headers.append(line.decode())
                    continue
                cols = line.split(b"\t", 5)
                if len(cols) < 6:
                    raise ValueError(f"Malformed PVAR row {idx}")
                chrom, pos, vid, ref, alt = (x.decode() for x in cols[:5])
                chrom = chrom.removeprefix("chr")
                if chrom not in AUTOSOMES:
                    nonautosomal[chrom] += 1
                    idx += 1
                    continue
                posi = int(pos)
                if chrom not in handles:
                    handles[chrom] = pysam.BGZFile(str(output / f"chr{chrom}.tsv.gz.part"), "w")
                    handles[chrom].write(b"#chrom\tpos\tid\tref\talt\tsource_index\tvariant_class\n")
                    counts[chrom] = {"count": 0, "min_pos": posi, "max_pos": posi,
                                     "source_start": idx, "source_stop": idx, "classes": Counter()}
                    first_rows[chrom] = line.decode().rstrip().split("\t")
                data = counts[chrom]
                if posi < data["max_pos"]:
                    raise ValueError(f"Unsorted source: chr{chrom}:{pos}")
                kind = variant_class(ref, alt)
                handles[chrom].write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{idx}\t{kind}\n".encode())
                data["count"] += 1
                data["max_pos"], data["source_stop"] = posi, idx + 1
                data["classes"][kind] += 1
                idx += 1
                if idx % 1000000 == 0:
                    state = {"stage": "COMPLETE_PVAR_CENSUS", "source_rows_scanned": idx,
                             "current_chromosome": chrom, "elapsed_seconds": time.time() - started}
                    print(json.dumps(state), flush=True)
                    if status_path:
                        atomic_json(status_path, state)
    finally:
        for handle in handles.values():
            handle.close()
    if set(counts) != set(AUTOSOMES):
        raise ValueError("Full reference is missing an autosome")
    outputs = {}
    for chrom in AUTOSOMES:
        part = output / f"chr{chrom}.tsv.gz.part"
        final = output / f"chr{chrom}.tsv.gz"
        part.replace(final)
        pysam.tabix_index(str(final), seq_col=0, start_col=1, end_col=1, force=True)
        outputs[chrom] = {"path": str(final), "sha256": sha256(final), "size": final.stat().st_size}
    (output / "source_header.txt").write_text("".join(headers))
    atomic_json(output / "source_first_rows.json", first_rows)
    report = {"source_path": str(source), "source_size": stat.st_size, "source_mtime_ns": stat.st_mtime_ns,
              "source_sha256": digest.hexdigest(), "total_variants": idx, "chromosomes": counts,
              "nonautosomal_variants": sum(nonautosomal.values()), "nonautosomal_contigs": len(nonautosomal),
              "outputs": outputs, "elapsed_seconds": time.time() - started, "state": "READY"}
    atomic_json(receipt, report)
    return report

def verify_local_sources(prefix: Path, census_dir: Path) -> dict:
    """Compare zstd decompression identities and PGEN/PVAR row counts."""
    import pgenlib
    import numpy as np
    import zstandard

    census = json.loads((census_dir / "source_census.json").read_text())
    validation_path = census_dir / "source_validation.json"
    if validation_path.exists() and (census_dir / "source_allele_offsets.npy").exists():
        old = json.loads(validation_path.read_text())
        unchanged = old.get("state") == "READY"
        for extension in ("pvar", "pgen"):
            previous = old.get("source", {}).get(extension, {})
            plain = Path(str(prefix) + "." + extension)
            compressed = Path(str(plain) + ".zst")
            unchanged &= (plain.stat().st_size == previous.get("size_bytes")
                          and plain.stat().st_mtime_ns == previous.get("mtime_ns")
                          and compressed.stat().st_size == previous.get("compressed_size_bytes")
                          and compressed.stat().st_mtime_ns == previous.get("compressed_mtime_ns"))
        if unchanged:
            return old
    source_sha = {"pvar": census["source_sha256"]}
    results = {}
    for extension in ("pvar", "pgen"):
        plain = Path(str(prefix) + "." + extension)
        compressed = Path(str(plain) + ".zst")
        if extension not in source_sha:
            source_sha[extension] = sha256(plain)
        digest = hashlib.sha256()
        total = 0
        with compressed.open("rb") as inp, zstandard.ZstdDecompressor().stream_reader(inp) as reader:
            if extension == "pvar":
                offsets_path = census_dir / "source_allele_offsets.npy.part"
                offsets = np.lib.format.open_memmap(offsets_path, mode="w+", dtype=np.uintp,
                                                    shape=(census["total_variants"] + 1,))
                offsets[0] = 0
                row_index = 0
                with io.BufferedReader(reader, buffer_size=4 * 1024**2) as lines:
                    for line in lines:
                        digest.update(line)
                        total += len(line)
                        if line.startswith(b"#"):
                            continue
                        fields = line.split(b"\t", 5)
                        offsets[row_index + 1] = offsets[row_index] + fields[4].count(b",") + 2
                        row_index += 1
                if row_index != census["total_variants"]:
                    raise ValueError("Compressed PVAR record count mismatch")
                offsets.flush()
                del offsets
                offsets_path.replace(census_dir / "source_allele_offsets.npy")
            else:
                for block in iter(lambda: reader.read(4 * 1024**2), b""):
                    digest.update(block)
                    total += len(block)
        if digest.hexdigest() != source_sha[extension] or total != plain.stat().st_size:
            raise ValueError(f"Existing compressed and decompressed {extension} disagree")
        results[extension] = {"path": str(plain), "sha256": source_sha[extension], "size_bytes": total,
                              "mtime_ns": plain.stat().st_mtime_ns,
                              "compressed_path": str(compressed), "compressed_sha256": sha256(compressed),
                              "compressed_size_bytes": compressed.stat().st_size,
                              "compressed_mtime_ns": compressed.stat().st_mtime_ns,
                              "zstd_integrity_and_decompressed_identity": True}
    samples = _table(Path(str(prefix) + ".psam"))
    with pgenlib.PgenReader(os.fsencode(str(prefix) + ".pgen"),
                           allele_idx_offsets=np.load(census_dir / "source_allele_offsets.npy", mmap_mode="r")) as reader:
        if reader.get_raw_sample_ct() != len(samples) or reader.get_variant_ct() != census["total_variants"]:
            raise ValueError("PGEN dimensions do not match the complete PSAM/PVAR")
        first_rows = json.loads((census_dir / "source_first_rows.json").read_text())
        probes = {}
        for chrom, row in first_rows.items():
            index = census["chromosomes"][chrom]["source_start"]
            values = np.empty(len(samples), dtype=np.int8)
            reader.read(index, values)
            info = dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
            ac = int(values[values >= 0].sum())
            an = int((values >= 0).sum() * 2)
            if ac != int(info["AC"]) or an != int(info["AN"]):
                raise ValueError(f"PGEN ALT orientation/order mismatch on chromosome {chrom}")
            probes[chrom] = {"source_index": index, "position": int(row[1]), "ref": row[3], "alt": row[4],
                             "computed_ac": ac, "published_ac": int(info["AC"]), "computed_an": an,
                             "published_an": int(info["AN"]), "match": True}
    report = {"release": RELEASE, "state": "READY", "source": results,
              "sample_count": len(samples), "variant_count": census["total_variants"],
              "orientation_probes": probes,
              "dosage_orientation": "pgenlib allele_idx=1; diploid counts 0/1/2 of PVAR ALT; -9 is missing"}
    atomic_json(validation_path, report)
    return report

def genotype_qc(values, *, maf_min: float = 0.01, missing_max: float = 0.05):
    """Variant-major hard-call arrays; calculate AF exclusively in this subset."""
    import numpy as np

    valid = values >= 0
    if not np.isin(values, [-9, 0, 1, 2]).all():
        raise ValueError("Unexpected diploid genotype encoding")
    observed = valid.sum(axis=1)
    missing = 1.0 - observed / values.shape[1]
    ac = np.where(valid, values, 0).sum(axis=1)
    af = np.divide(ac, 2 * observed, out=np.full(len(values), np.nan), where=observed > 0)
    maf = np.minimum(af, 1 - af)
    # Mean imputation preserves a constant observed column, including the
    # all-heterozygous AF=0.5 case; allele polymorphism alone is insufficient.
    minimum = np.where(valid, values, 3).min(axis=1)
    maximum = np.where(valid, values, -1).max(axis=1)
    zero_variance = (observed < 2) | (minimum == maximum)
    reasons = np.full(len(values), "PASS", dtype=object)
    reasons[maf < maf_min] = "reference_MAF_lt_0.01"
    reasons[zero_variance | ~np.isfinite(maf)] = "reference_zero_variance"
    # Exact integer comparison avoids a binary floating-point boundary accident.
    reasons[(values.shape[1] - observed) * 100 > missing_max * 100 * values.shape[1]] = "reference_missing_rate_gt_0.05"
    return af, maf, missing, reasons
