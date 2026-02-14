#!/usr/bin/env python
"""
04_call_peaks.py - Call peaks for IP samples using CLIPper

This script calls peaks for IP samples using CLIPper (Yeo lab eCLIP),
then normalizes peaks against size-matched input controls.

Usage:
    # Process all samples
    python scripts/04_call_peaks.py

    # Process single sample by index (for SLURM array jobs)
    python scripts/04_call_peaks.py --sample-index 0

    # List sample indices
    python scripts/04_call_peaks.py --list-samples
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import (
    SAMPLE_INFO, ALIGNED_DIR, PEAKS_DIR, LOGS_DIR,
    CLIPPER_SPECIES, CLIPPER_FDR,
    CLIPPER_MODE, CLIPPER_CONTAINER,
    PEAK_FOLD_ENRICHMENT, PEAK_PVALUE_THRESHOLD
)
from analysis import call_peaks_clipper, normalize_peaks_with_input


def get_ip_sample_pairs():
    """Get list of (IP sample, Input sample) pairs."""
    peak_calls = []

    for sample_name, sample_info in SAMPLE_INFO.items():
        if sample_info['sample_type'] != 'IP':
            continue

        # Find corresponding input
        if 'IFIT2_IFIT3' in sample_name and 'FLAG' in sample_name:
            input_sample = 'IFIT2_IFIT3_FLAG_SMInput'
        elif 'IFIT2_IFIT3' in sample_name and 'HA' in sample_name:
            input_sample = 'IFIT2_IFIT3_HA_SMInput'
        elif 'IFIT2_alone' in sample_name:
            input_sample = 'IFIT2_alone_SMInput'
        elif 'IFIT3_alone' in sample_name:
            input_sample = 'IFIT3_alone_SMInput'
        elif 'WT_FLAG' in sample_name:
            input_sample = 'WT_FLAG_SMInput'
        elif 'WT_HA' in sample_name:
            input_sample = 'WT_HA_SMInput'
        else:
            continue

        peak_calls.append((sample_name, input_sample))

    return peak_calls


def process_sample(ip_sample, input_sample, logger):
    """Process a single IP sample: peak calling + normalization."""
    ip_bam = ALIGNED_DIR / f"{ip_sample}_Aligned.sortedByCoord.dedup.bam"
    input_bam = ALIGNED_DIR / f"{input_sample}_Aligned.sortedByCoord.dedup.bam"

    # Check BAM files exist
    if not ip_bam.exists():
        logger.error(f"IP BAM not found: {ip_bam}")
        return False

    if not input_bam.exists():
        logger.error(f"Input BAM not found: {input_bam}")
        return False

    # Validate BAM file before processing
    logger.info(f"Validating BAM file: {ip_bam}")
    import subprocess
    result = subprocess.run(
        ['samtools', 'quickcheck', str(ip_bam)],
        capture_output=True
    )
    if result.returncode != 0:
        logger.error(f"BAM file failed validation: {ip_bam}")
        logger.error("Run: samtools quickcheck <bam> to diagnose")
        return False

    # Check read count
    result = subprocess.run(
        ['samtools', 'view', '-c', str(ip_bam)],
        capture_output=True, text=True
    )
    read_count = int(result.stdout.strip()) if result.stdout.strip() else 0
    logger.info(f"  Read count: {read_count:,}")

    if read_count == 0:
        logger.error(f"BAM file has no reads: {ip_bam}")
        return False

    logger.info("\n" + "="*70)
    logger.info(f"Processing: {ip_sample}")
    logger.info("="*70)

    # Step 1: Call peaks with CLIPper
    peaks_file = call_peaks_clipper(
        ip_bam,
        PEAKS_DIR,
        ip_sample,
        species=CLIPPER_SPECIES,
        fdr=CLIPPER_FDR,
        mode=CLIPPER_MODE,
        container_path=CLIPPER_CONTAINER
    )

    if not peaks_file:
        logger.error(f"Peak calling failed for {ip_sample}")
        return False

    # Step 2: Normalize peaks against input
    normalized_peaks = normalize_peaks_with_input(
        peaks_file,
        ip_bam,
        input_bam,
        PEAKS_DIR,
        ip_sample,
        fold_change_threshold=PEAK_FOLD_ENRICHMENT,
        pvalue_threshold=PEAK_PVALUE_THRESHOLD
    )

    if normalized_peaks:
        logger.info(f"Peak calling complete for {ip_sample}")
        return True
    else:
        logger.error(f"Normalization failed for {ip_sample}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Call peaks for eCLIP IP samples using CLIPper'
    )
    parser.add_argument(
        '--sample-index', '-i',
        type=int,
        default=None,
        help='Process only the sample at this index (0-5). For SLURM array jobs.'
    )
    parser.add_argument(
        '--list-samples', '-l',
        action='store_true',
        help='List sample indices and exit'
    )

    args = parser.parse_args()

    # Get sample pairs
    peak_calls = get_ip_sample_pairs()

    # List samples mode
    if args.list_samples:
        print("IP Sample Indices:")
        print("-" * 50)
        for i, (ip_sample, input_sample) in enumerate(peak_calls):
            print(f"  {i}: {ip_sample} (input: {input_sample})")
        print("-" * 50)
        print(f"Total: {len(peak_calls)} IP samples")
        print("\nUsage: python scripts/04_call_peaks.py --sample-index <N>")
        return

    # Setup logging
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    if args.sample_index is not None:
        # Single sample mode - log to sample-specific file
        log_file = LOGS_DIR / f'peak_calling_sample_{args.sample_index}.log'
    else:
        log_file = LOGS_DIR / 'peak_calling.log'

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)

    # Header
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Peak Calling with CLIPper")
    logger.info("="*70)
    logger.info("\nUsing Yeo lab CLIPper for peak calling")
    logger.info(f"  Execution mode: {CLIPPER_MODE}")
    if CLIPPER_MODE in ['singularity', 'docker'] and CLIPPER_CONTAINER:
        logger.info(f"  Container: {CLIPPER_CONTAINER}")
    logger.info(f"  Species: {CLIPPER_SPECIES}")
    logger.info(f"  FDR: {CLIPPER_FDR}")
    logger.info(f"  Fold-change threshold: {PEAK_FOLD_ENRICHMENT}")
    logger.info(f"  P-value threshold: {PEAK_PVALUE_THRESHOLD}")

    # Single sample mode
    if args.sample_index is not None:
        if args.sample_index < 0 or args.sample_index >= len(peak_calls):
            logger.error(f"Invalid sample index: {args.sample_index}")
            logger.error(f"Valid range: 0-{len(peak_calls)-1}")
            sys.exit(1)

        ip_sample, input_sample = peak_calls[args.sample_index]
        logger.info(f"\nProcessing single sample (index {args.sample_index}): {ip_sample}")

        success = process_sample(ip_sample, input_sample, logger)

        if success:
            logger.info(f"\nSample {ip_sample} completed successfully!")
            sys.exit(0)
        else:
            logger.error(f"\nSample {ip_sample} failed!")
            sys.exit(1)

    # All samples mode
    logger.info(f"\nCalling peaks for {len(peak_calls)} IP samples\n")

    successful = 0
    failed = 0

    for ip_sample, input_sample in peak_calls:
        if process_sample(ip_sample, input_sample, logger):
            successful += 1
        else:
            failed += 1

    # Summary
    logger.info("\n" + "="*70)
    logger.info("Peak Calling Summary")
    logger.info("="*70)
    logger.info(f"Successful: {successful}/{len(peak_calls)}")
    logger.info(f"Failed: {failed}/{len(peak_calls)}")
    logger.info("\n" + "="*70)
    logger.info("Peak calling complete!")
    logger.info("="*70)
    logger.info("\nOutput files:")
    logger.info(f"  {PEAKS_DIR}/*_clipper_peaks.bed (raw CLIPper peaks)")
    logger.info(f"  {PEAKS_DIR}/*_normalized_peaks.bed (input-normalized peaks)")
    logger.info("\nNext steps:")
    logger.info("  python scripts/05_analyze_utrs.py")

    sys.exit(0 if failed == 0 else 1)


if __name__ == '__main__':
    main()
