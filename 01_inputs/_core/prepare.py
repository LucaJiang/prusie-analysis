"""Allele orientation, normalized genotype factors and fixed RSS settings."""
from pathlib import Path
import numpy as np
from .reference import normalize_variant

COMPLEMENT = str.maketrans("ACGT", "TGCA")

def orient_variant(pos: int, ref: str, alt: str, sequence: str) -> tuple[int, str, str, int, str]:
    """Validate declared GRCh38 REF first; resolve alternative orientation only uniquely.

    Direct FASTA-supported REF/ALT is unambiguous under both release contracts,
    including palindromes. For a REF mismatch, swaps and SNV strand transforms
    must agree in both canonical key and effect sign, otherwise reject.
    """
    ref, alt = ref.upper(), alt.upper()
    if not ref or not alt or set(ref + alt) - set("ACGT") or ref == alt:
        raise ValueError("unsupported_alleles")
    if max(len(ref), len(alt)) > 50 or (len(ref) == len(alt) and len(ref) != 1):
        raise ValueError("not_SNV_or_short_indel")
    try:
        return (*normalize_variant(pos, ref, alt, sequence), 1, "declared_REF_ALT")
    except ValueError as direct_error:
        if str(direct_error) != "REF_mismatch":
            raise
    alternatives = [(alt, ref, -1, "REF_ALT_swap")]
    if len(ref) == len(alt) == 1:
        alternatives += [(ref.translate(COMPLEMENT), alt.translate(COMPLEMENT), 1, "strand_complement"),
                         (alt.translate(COMPLEMENT), ref.translate(COMPLEMENT), -1, "strand_complement_swap")]
    valid = {}
    for new_ref, new_alt, sign, label in alternatives:
        try:
            key = (*normalize_variant(pos, new_ref, new_alt, sequence), sign)
            valid[key] = label
        except ValueError:
            pass
    if len(valid) != 1:
        raise ValueError("ambiguous_orientation" if valid else "REF_mismatch")
    key, label = next(iter(valid.items()))
    return (*key, label)

def factor_from_dosages(dosages: np.ndarray, *, maf_min: float = .01,
                       missing_max: float = .05) -> tuple[np.ndarray, dict]:
    """Mean imputation, centering and unit column norm within this reference."""
    x = np.asarray(dosages, dtype=np.float64).copy()
    if x.ndim != 2 or x.shape[0] < 2 or np.any(np.isinf(x)):
        raise ValueError("Invalid genotype dimensions/values")
    observed = np.isfinite(x)
    if np.any((x[observed] < 0) | (x[observed] > 2)):
        raise ValueError("ALT genotype outside 0..2")
    missing = (~observed).mean(axis=0)
    if np.any(missing > missing_max + 1e-15):
        raise ValueError("Reference missingness exceeds frozen threshold")
    means = np.nanmean(x, axis=0)
    af = means / 2
    maf = np.minimum(af, 1 - af)
    if np.any(~np.isfinite(af)) or np.any(maf < maf_min - 1e-15):
        raise ValueError("Reference MAF below frozen threshold")
    rows, columns = np.where(~observed)
    x[rows, columns] = means[columns]
    x -= means
    norms = np.sqrt(np.sum(x * x, axis=0))
    if np.any(norms <= 0):
        raise ValueError("Reference zero-variance variant")
    x /= norms
    return x, {"n_ref": len(x), "af": af, "maf": maf, "missing_rate": missing}

def common_parameters(settings: dict, n: int, p: int) -> dict:
    keys = ("L", "max_iter", "tol", "coverage", "min_abs_corr", "scaled_prior_variance",
            "estimate_prior_variance", "estimate_prior_method", "estimate_residual_variance", "residual_variance",
            "null_weight", "prior_tol", "r_tol", "check_null_threshold", "standardize", "check_input", "check_prior", "refine")
    return {**{k: settings[k] for k in keys}, "n": n, "n_purity": p}
