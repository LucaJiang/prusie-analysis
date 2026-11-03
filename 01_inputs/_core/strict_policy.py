"""Current-record biallelic SNV policy."""
from .strict_common import classify as legacy_classify

POLICY = "strict_snp_maf05_current_record_v2"

ALLELE_SCOPE = "current input record; historical unsplit-site ancestry is not a QC dependency"

STUDY_STAGES = ["invalid_coordinate", "current_record_multiallelic", "invalid_alleles",
                "indel", "MNV", "palindrome", "invalid_statistics", "study_MAF_lt_0.05",
                "source_FASTA_mismatch", "mapping_failure", "allele_mismatch",
                "conflicting_duplicate", "identical_duplicate_collapsed",
                "reference_MAF_lt_0.05", "reference_missing_gt_0.02", "reference_HWE_midp_lt_1e-6",
                "reference_conflicting_duplicate", "reference_allele_mismatch", "reference_QC_or_absent", "PASS"]

REFERENCE_STAGES = ["invalid_coordinate", "current_record_multiallelic", "invalid_alleles", "indel", "MNV",
                    "palindrome", "REF_mismatch", "conflicting_duplicate", "identical_duplicate_collapsed",
                    "reference_MAF_lt_0.05", "reference_missing_gt_0.02",
                    "reference_HWE_midp_lt_1e-6", "zero_variance", "PASS"]

def classify(a, b):
    reason = legacy_classify(a, b)
    return "current_record_multiallelic" if reason == "original_multiallelic" else reason
