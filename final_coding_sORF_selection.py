#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd

INPUT_DIR = Path("../")
INPUT_SUFFIX = "_sORF_results.tsv"
OUTPUT_PREFIX = "SNHGs"

# coding potential thresholds
CANONICAL_MIN_PREDICTED_P_SITES = 3
CANONICAL_MIN_FRAME0_RATIO = 50.0
NONCANONICAL_MIN_PREDICTED_P_SITES = 4
NONCANONICAL_MIN_FRAME0_RATIO = 60.0

CANONICAL_START_CODON = "ATG"

OUTPUT_COLUMNS = [
    "chromosome",
    "gene_name",
    "gene_id",
    "transcript_id",
    "strand",
    "genomic_start_loci",
    "genomic_end_loci",
    "orf_length_bp",
    "start_codon",
    "is_canonical",
    "aggregated_predicted_orf_p_sites",
    "aggregated_total_orf_p_sites",
    "aggregated_frame_0_ratio (%)",
    "supported_samples",
    "max_TPM",
    "max_percentage",
    "orf_dna_sequence",
    "amino_acid_sequence",
]


def usage():
    print("Usage: python final_coding_ORF_selection.py <sample1> <sample2> [...]")
    print("Example: python final_coding_ORF_selection.py DLD1_alone DLD1_co")


def natural_chrom_key(chrom):
    chrom = str(chrom)
    c = chrom.replace("chr", "")
    if c.isdigit():
        return (0, int(c))
    if c == "X":
        return (1, 23)
    if c == "Y":
        return (1, 24)
    if c in ["M", "MT"]:
        return (1, 25)
    return (2, c)


def read_sample_file(sample):
    path = INPUT_DIR / f"{sample}{INPUT_SUFFIX}"
    if not path.exists():
        raise FileNotFoundError(f"Input file not found for sample '{sample}': {path}")

    df = pd.read_csv(path, sep="\t")
    df["sample"] = sample

    required = [
        "chromosome",
        "gene_name",
        "gene_id",
        "transcript_id",
        "strand",
        "genomic_start_loci",
        "genomic_end_loci",
        "orf_length_bp",
        "start_codon",
        "predicted_orf_p_sites",
        "total_orf_p_sites",
        "frame_0_ratio (%)",
        "orf_dna_sequence",
        "amino_acid_sequence",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    for col in [
        "genomic_start_loci",
        "genomic_end_loci",
        "orf_length_bp",
        "predicted_orf_p_sites",
        "total_orf_p_sites",
        "frame_0_ratio (%)",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    if "TPM" not in df.columns:
        df["TPM"] = 0
    if "percentage" not in df.columns:
        df["percentage"] = 0

    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce").fillna(0)
    df["percentage"] = pd.to_numeric(df["percentage"], errors="coerce").fillna(0)

    df["is_canonical"] = df["start_codon"].astype(str).str.upper().eq(CANONICAL_START_CODON)

    return df


def make_orf_key(df):
    return (
        df["chromosome"].astype(str)
        + ":"
        + df["genomic_start_loci"].astype(int).astype(str)
        + "-"
        + df["genomic_end_loci"].astype(int).astype(str)
        + ":"
        + df["strand"].astype(str)
        + ":"
        + df["transcript_id"].astype(str)
        + ":"
        + df["start_codon"].astype(str)
    )


def aggregate_identical_orfs(df):

    df = df.copy()
    df["orf_key"] = make_orf_key(df)
    df["weighted_ratio_numerator"] = df["predicted_orf_p_sites"] * df["frame_0_ratio (%)"]

    group_cols = ["orf_key"]

    agg = df.groupby(group_cols, as_index=False).agg(
        chromosome=("chromosome", "first"),
        gene_name=("gene_name", "first"),
        gene_id=("gene_id", "first"),
        transcript_id=("transcript_id", "first"),
        strand=("strand", "first"),
        genomic_start_loci=("genomic_start_loci", "first"),
        genomic_end_loci=("genomic_end_loci", "first"),
        orf_length_bp=("orf_length_bp", "first"),
        start_codon=("start_codon", "first"),
        is_canonical=("is_canonical", "first"),
        aggregated_predicted_orf_p_sites=("predicted_orf_p_sites", "sum"),
        aggregated_total_orf_p_sites=("total_orf_p_sites", "sum"),
        weighted_ratio_numerator=("weighted_ratio_numerator", "sum"),
        supported_samples=("sample", lambda x: ",".join(sorted(set(map(str, x))))),
        n_supported_samples=("sample", lambda x: len(set(map(str, x)))),
        max_TPM=("TPM", "max"),
        max_percentage=("percentage", "max"),
        orf_dna_sequence=("orf_dna_sequence", "first"),
        amino_acid_sequence=("amino_acid_sequence", "first"),
    )

    agg["aggregated_frame_0_ratio (%)"] = agg.apply(
        lambda row: (
            row["weighted_ratio_numerator"] / row["aggregated_predicted_orf_p_sites"]
            if row["aggregated_predicted_orf_p_sites"] > 0
            else 0.0
        ),
        axis=1,
    )

    return agg.drop(columns=["weighted_ratio_numerator"])


def assign_coding_potential(df):
    df = df.copy()

    canonical_pass = (
        (df["is_canonical"])
        & (df["aggregated_predicted_orf_p_sites"] >= CANONICAL_MIN_PREDICTED_P_SITES)
        & (df["aggregated_frame_0_ratio (%)"] >= CANONICAL_MIN_FRAME0_RATIO)
    )

    noncanonical_pass = (
        (~df["is_canonical"])
        & (df["aggregated_predicted_orf_p_sites"] >= NONCANONICAL_MIN_PREDICTED_P_SITES)
        & (df["aggregated_frame_0_ratio (%)"] >= NONCANONICAL_MIN_FRAME0_RATIO)
    )

    df["coding_potential"] = canonical_pass | noncanonical_pass
    return df

# priority: Canonical ORF > Num. of p-site > ORF length
def select_one_orf_per_transcript(df):

    if df.empty:
        return df

    df = df.copy()
    df = df.sort_values(
        by=[
            "transcript_id",
            "is_canonical",
            "aggregated_predicted_orf_p_sites",
            "orf_length_bp",
            "aggregated_frame_0_ratio (%)",
            "genomic_start_loci",
        ],
        ascending=[True, False, False, False, False, True],
    )

    return df.drop_duplicates(subset=["transcript_id"], keep="first").reset_index(drop=True)


def sort_output(df):
    df = df.copy()
    df["_chrom_sort"] = df["chromosome"].apply(natural_chrom_key)
    df = df.sort_values(
        by=["_chrom_sort", "genomic_start_loci", "genomic_end_loci", "gene_name", "transcript_id"],
        ascending=[True, True, True, True, True],
    ).reset_index(drop=True)
    return df.drop(columns=["_chrom_sort"])


def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    samples = sys.argv[1:]
    print(f"Aggregating samples: {', '.join(samples)}")

    dfs = []
    for sample in samples:
        print(f"Loading sample: {sample}")
        dfs.append(read_sample_file(sample))

    combined = pd.concat(dfs, ignore_index=True)
    aggregated = aggregate_identical_orfs(combined)
    aggregated = assign_coding_potential(aggregated)

    all_path = f"{OUTPUT_PREFIX}.all_aggregated_ORFs.tsv"
    aggregated_out = sort_output(aggregated)
    aggregated_out[OUTPUT_COLUMNS].to_csv(all_path, sep="\t", index=False)

    coding = aggregated[aggregated["coding_potential"]].copy()
    final = select_one_orf_per_transcript(coding)
    final = sort_output(final)

    final_path = f"{OUTPUT_PREFIX}.final_coding_ORFs.tsv"
    final[OUTPUT_COLUMNS].to_csv(final_path, sep="\t", index=False)

    print(f"Saved all aggregated ORFs: {all_path}")
    print(f"Number of aggregated ORFs: {len(aggregated_out)}")
    print(f"Saved final coding ORFs: {final_path}")
    print(f"Number of final coding ORFs: {len(final)}")


if __name__ == "__main__":
    main()
