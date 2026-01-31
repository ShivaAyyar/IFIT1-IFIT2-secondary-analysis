#!/usr/bin/env python
"""
03_process_samples.py - Process eCLIP samples (QC, trim, align, dedup)

This script processes a single eCLIP sample through the complete pipeline:
1. FastQC quality control
2. UMI extraction
3. Adapter trimming
4. STAR alignment
5. PCR duplicate removal

Usage:
    # Process specific sample
    python scripts/03_process_samples.py --sample IFIT2_IFIT3_FLAG_IP

    # Process with custom threads
    python scripts/03_process_samples.py --sample IFIT2_IFIT3_FLAG_IP --threads 8

    # Process sample by index (for SLURM array jobs)
    python scripts/03_process_samples.py --sample-index 0
"""

import argparse
import logging
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    SAMPLE_INFO, RAW_DIR, QC_DIR, TRIMMED_DIR, ALIGNED_DIR,
    REFERENCE_DIR, UMI_LENGTH, MIN_READ_LENGTH, QUALITY_CUTOFF
)
from utils import (
    run_fastqc, extract_umi_and_trim_single_end,
    run_star_alignment_single_end, deduplicate_bam,
    check_tool_available
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/sample_processing.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def process_sample(sample_name, sample_info, threads=8):
    """Process a single eCLIP sample through the complete pipeline."""

    logger.info("="*70)
    logger.info(f"Processing sample: {sample_name}")
    logger.info("="*70)
    logger.info(f"Description: {sample_info['description']}")
    logger.info(f"SRR: {sample_info['srr_accession']}")
    logger.info("")

    srr = sample_info['srr_accession']

    # Find FASTQ file
    fastq_candidates = list(RAW_DIR.glob(f"{srr}*.fastq.gz"))
    if not fastq_candidates:
        logger.error(f"FASTQ file not found for {srr} in {RAW_DIR}")
        logger.error("Run 02_download_data.py first")
        return False

    fastq = fastq_candidates[0]
    logger.info(f"Input FASTQ: {fastq}")

    # Step 1: FastQC
    logger.info("\n[Step 1/4] Running FastQC...")
    run_fastqc([fastq], QC_DIR, threads=threads)

    # Step 2: UMI extraction and adapter trimming
    logger.info("\n[Step 2/4] Extracting UMI and trimming adapters...")
    trim_result = extract_umi_and_trim_single_end(
        fastq,
        TRIMMED_DIR,
        umi_length=UMI_LENGTH,
        min_length=MIN_READ_LENGTH,
        quality_cutoff=QUALITY_CUTOFF,
        threads=threads
    )

    if not trim_result:
        logger.error("Trimming failed")
        return False

    trimmed_fastq = trim_result['fastq']
    logger.info(f"Trimmed FASTQ: {trimmed_fastq}")

    # Step 3: STAR alignment
    logger.info("\n[Step 3/4] Aligning with STAR...")
    star_index = REFERENCE_DIR / 'STAR_index_hg19'

    if not (star_index / 'SA').exists():
        logger.error(f"STAR index not found at {star_index}")
        logger.error("Run 01_setup_reference.py first")
        return False

    align_result = run_star_alignment_single_end(
        trimmed_fastq,
        star_index,
        ALIGNED_DIR,
        sample_name,
        threads=threads
    )

    if not align_result:
        logger.error("Alignment failed")
        return False

    bam_file = align_result['bam']
    logger.info(f"Aligned BAM: {bam_file}")

    # Step 4: Deduplicate
    logger.info("\n[Step 4/4] Removing PCR duplicates...")
    dedup_result = deduplicate_bam(bam_file, ALIGNED_DIR, use_umi=True)

    if not dedup_result:
        logger.error("Deduplication failed")
        return False

    logger.info(f"Deduplicated BAM: {dedup_result['bam']}")

    logger.info("\n" + "="*70)
    logger.info(f"✓ Sample processing complete: {sample_name}")
    logger.info("="*70)

    return True


def main():
    parser = argparse.ArgumentParser(
        description='Process eCLIP samples (QC, trim, align, dedup)'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--sample',
        help='Sample name to process'
    )
    group.add_argument(
        '--sample-index',
        type=int,
        help='Sample index (for SLURM array jobs)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads (default: 8)'
    )

    args = parser.parse_args()

    # Check required tools
    required_tools = ['fastqc', 'umi_tools', 'cutadapt', 'STAR', 'samtools']
    missing_tools = [t for t in required_tools if not check_tool_available(t)]

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        sys.exit(1)

    # Determine sample to process
    if args.sample:
        sample_name = args.sample
        if sample_name not in SAMPLE_INFO:
            logger.error(f"Sample '{sample_name}' not found")
            logger.error(f"Available: {', '.join(SAMPLE_INFO.keys())}")
            sys.exit(1)
        sample_info = SAMPLE_INFO[sample_name]
    else:
        # Use sample index
        sample_list = list(SAMPLE_INFO.items())
        if args.sample_index >= len(sample_list):
            logger.error(f"Sample index {args.sample_index} out of range (0-{len(sample_list)-1})")
            sys.exit(1)
        sample_name, sample_info = sample_list[args.sample_index]

    # Process sample
    success = process_sample(sample_name, sample_info, threads=args.threads)

    if not success:
        sys.exit(1)


if __name__ == '__main__':
    main()
