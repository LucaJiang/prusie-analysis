"""Allele classification and QC ledger utilities."""
from pathlib import Path
from .acquire import sha256

COMPLEMENT = str.maketrans("ACGT", "TGCA")

def classify(a: str, b: str) -> str:
    if "," in a or "," in b:
        return "original_multiallelic"
    if not a or not b or set(a + b) - set("ACGT") or a == b:
        return "invalid_alleles"
    if len(a) != len(b):
        return "indel"
    if len(a) != 1:
        return "MNV"
    if a + b in ("AT", "TA", "CG", "GC"):
        return "palindrome"
    return "PASS"

def allele_detail(other: str, effect: str) -> str:
    """Raw representation categories, independent of the first-failure ledger."""
    if "," in other or "," in effect:return "multiallelic_record"
    if other in ("", ".", "0", "NA") or effect in ("", ".", "0", "NA"):return "missing_allele"
    if any(c in other+effect for c in "<>[]*"):return "symbolic_or_breakend"
    if set(other+effect)-set("ACGT"):return "invalid_or_complex_allele"
    if other==effect:return "identical_alleles"
    if len(effect)>len(other):return "insertion"
    if len(effect)<len(other):return "deletion"
    if len(effect)>1:return "MNV"
    return "palindromic_SNV" if other+effect in ("AT","TA","CG","GC") else "nonpalindromic_SNV"

def orient(other: str, effect: str, ref: str, alt: str) -> tuple[int, str]:
    """Align the source effect to target ALT; a complement alone keeps sign."""
    if classify(other, effect) != "PASS" or classify(ref, alt) != "PASS":
        raise ValueError("ineligible_or_ambiguous_alleles")
    transforms = [(other, effect, 1, "direct"), (effect, other, -1, "swap"),
                  (other.translate(COMPLEMENT), effect.translate(COMPLEMENT), 1, "complement"),
                  (effect.translate(COMPLEMENT), other.translate(COMPLEMENT), -1, "complement_swap")]
    matches = [(sign, label) for a, b, sign, label in transforms if (a, b) == (ref, alt)]
    if len(matches) != 1:
        raise ValueError("allele_mismatch" if not matches else "ambiguous_orientation")
    return matches[0]

def stages(counts: dict, order: list[str]) -> list[dict]:
    remaining = sum(counts.values())
    result = []
    for stage in order:
        if stage == "PASS":
            continue
        removed = int(counts.get(stage, 0))
        result.append(dict(stage=stage, before=remaining, removed=removed, after=remaining-removed))
        remaining -= removed
    if remaining != counts.get("PASS", 0):
        raise ValueError(f"Unaccounted first-failure reasons: {set(counts)-set(order)}")
    return result

def file_record(path: Path, fresh=True) -> dict:
    st = path.stat()
    return dict(path=str(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns,
                sha256=sha256(path) if fresh else None, hash_method="fresh_full_SHA256" if fresh else "stat_only")
