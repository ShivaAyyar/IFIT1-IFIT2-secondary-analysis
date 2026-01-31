"""
Configuration file for IFIT2-IFIT3 eCLIP Analysis Pipeline

This file contains all sample information, directory paths, and parameters
used throughout the analysis pipeline.
"""

from pathlib import Path

# ============================================================================
# Project Directory Structure
# ============================================================================

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw_fastq"
QC_DIR = DATA_DIR / "qc"
TRIMMED_DIR = DATA_DIR / "trimmed"
ALIGNED_DIR = DATA_DIR / "aligned"
PEAKS_DIR = DATA_DIR / "peaks"
REFERENCE_DIR = DATA_DIR / "reference"
RESULTS_DIR = DATA_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
LOGS_DIR = PROJECT_DIR / "logs"

# Create directories
for dir_path in [DATA_DIR, RAW_DIR, QC_DIR, TRIMMED_DIR, ALIGNED_DIR,
                 PEAKS_DIR, REFERENCE_DIR, RESULTS_DIR, FIGURES_DIR, LOGS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Sample Information - GSE284636 Uninfected HEK293T eCLIP
# ============================================================================

# CRITICAL NOTES:
#  - These are SINGLE-END 75nt reads (not paired-end!)
#  - FLAG antibody pulls down IFIT2:FLAG
#  - HA antibody pulls down IFIT3:HA
#  - In co-expression samples, separate IPs distinguish IFIT2 vs IFIT3 binding

SAMPLE_INFO = {
    # ==================== IFIT2+IFIT3 Co-expression ====================
    # When both proteins are expressed, FLAG IP pulls down IFIT2, HA IP pulls down IFIT3

    "IFIT2_IFIT3_FLAG_IP": {
        "description": "IFIT2+IFIT3 co-expression, FLAG IP (IFIT2 pulldown)",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG+IFIT3:HA",
        "ip_antibody": "FLAG",
        "target_protein": "IFIT2",
        "sample_type": "IP",
        "srr_accession": "SRR31773513",
        "gsm_accession": "GSM8689674"
    },
    "IFIT2_IFIT3_HA_IP": {
        "description": "IFIT2+IFIT3 co-expression, HA IP (IFIT3 pulldown)",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG+IFIT3:HA",
        "ip_antibody": "HA",
        "target_protein": "IFIT3",
        "sample_type": "IP",
        "srr_accession": "SRR31773515",
        "gsm_accession": "GSM8689672"
    },
    "IFIT2_IFIT3_FLAG_SMInput": {
        "description": "IFIT2+IFIT3 co-expression, FLAG SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG+IFIT3:HA",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773512",
        "gsm_accession": "GSM8689675"
    },
    "IFIT2_IFIT3_HA_SMInput": {
        "description": "IFIT2+IFIT3 co-expression, HA SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG+IFIT3:HA",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773514",
        "gsm_accession": "GSM8689673"
    },

    # ==================== IFIT2 Alone ====================
    "IFIT2_alone_FLAG_IP": {
        "description": "IFIT2 alone, FLAG IP",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG",
        "ip_antibody": "FLAG",
        "target_protein": "IFIT2",
        "sample_type": "IP",
        "srr_accession": "SRR31773519",
        "gsm_accession": "GSM8689668"
    },
    "IFIT2_alone_SMInput": {
        "description": "IFIT2 alone, SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT2:FLAG",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773518",
        "gsm_accession": "GSM8689669"
    },

    # ==================== IFIT3 Alone ====================
    "IFIT3_alone_HA_IP": {
        "description": "IFIT3 alone, HA IP",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT3:HA",
        "ip_antibody": "HA",
        "target_protein": "IFIT3",
        "sample_type": "IP",
        "srr_accession": "SRR31773517",
        "gsm_accession": "GSM8689670"
    },
    "IFIT3_alone_SMInput": {
        "description": "IFIT3 alone, SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "IFIT3:HA",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773516",
        "gsm_accession": "GSM8689671"
    },

    # ==================== Wild-Type Controls ====================
    "WT_FLAG_IP": {
        "description": "Wild-type, FLAG IP control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "WT",
        "ip_antibody": "FLAG",
        "target_protein": "control",
        "sample_type": "IP",
        "srr_accession": "SRR31773539",
        "gsm_accession": "GSM8689664"
    },
    "WT_HA_IP": {
        "description": "Wild-type, HA IP control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "WT",
        "ip_antibody": "HA",
        "target_protein": "control",
        "sample_type": "IP",
        "srr_accession": "SRR31773541",
        "gsm_accession": "GSM8689662"
    },
    "WT_FLAG_SMInput": {
        "description": "Wild-type, FLAG SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "WT",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773538",
        "gsm_accession": "GSM8689665"
    },
    "WT_HA_SMInput": {
        "description": "Wild-type, HA SMInput control",
        "cell_line": "HEK293T",
        "condition": "uninfected",
        "ifit_expression": "WT",
        "ip_antibody": "none",
        "target_protein": "input",
        "sample_type": "SMInput",
        "srr_accession": "SRR31773540",
        "gsm_accession": "GSM8689663"
    }
}

# ============================================================================
# Reference Genome URLs
# ============================================================================

REFERENCE_URLS = {
    'genome_fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
    'gencode_gtf': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz',
    'chrom_sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
}

# ============================================================================
# Analysis Parameters
# ============================================================================

# UMI and trimming parameters
UMI_LENGTH = 10  # First 10 nt are UMI
MIN_READ_LENGTH = 18  # Minimum length after trimming
QUALITY_CUTOFF = 10  # Quality trimming threshold
ADAPTER_SEQUENCE = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # Illumina TruSeq adapter

# STAR alignment parameters
STAR_THREADS = 8
STAR_SJDB_OVERHANG = 74  # For 75nt single-end reads

# Peak calling parameters
PEAK_MIN_COVERAGE = 5  # Minimum read coverage to call a peak
PEAK_FOLD_ENRICHMENT = 2.0  # Minimum fold enrichment over input

# UTR analysis parameters
UTR_SHORT_THRESHOLD = 50  # Threshold for defining "short" 5' UTRs (nt)

# Required tools
REQUIRED_TOOLS = {
    'fasterq-dump': 'SRA Toolkit (for downloading FASTQ)',
    'fastqc': 'FastQC (quality control)',
    'cutadapt': 'Cutadapt (adapter trimming)',
    'STAR': 'STAR aligner (alignment)',
    'samtools': 'SAMtools (BAM processing)',
    'bedtools': 'BEDTools (genomic intervals)',
    'umi_tools': 'UMI-tools (UMI deduplication)'
}
