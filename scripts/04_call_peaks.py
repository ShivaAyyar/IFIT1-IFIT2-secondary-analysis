#!/usr/bin/env python
"""
04_call_peaks.py - Call peaks for all IP samples using CLIPper

This script calls peaks for all IP samples using CLIPper (Yeo lab eCLIP),
then normalizes peaks against size-matched input controls.

Usage:
    python scripts/04_call_peaks.py
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import (
    SAMPLE_INFO, ALIGNED_DIR, PEAKS_DIR,
    CLIPPER_SPECIES, CLIPPER_FDR,
    PEAK_FOLD_ENRICHMENT, PEAK_PVALUE_THRESHOLD
)
from analysis import call_peaks_clipper, normalize_peaks_with_input

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/peak_calling.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Peak Calling with CLIPper")
    logger.info("="*70)
    logger.info("\nUsing Yeo lab CLIPper for peak calling")
    logger.info(f"  Species: {CLIPPER_SPECIES}")
    logger.info(f"  FDR: {CLIPPER_FDR}")
    logger.info(f"  Fold-change threshold: {PEAK_FOLD_ENRICHMENT}")
    logger.info(f"  P-value threshold: {PEAK_PVALUE_THRESHOLD}")

    # Define IP-input pairs
    peak_calls = []

    # Process each IP sample
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
            logger.warning(f"No input found for {sample_name}, skipping")
            continue

        peak_calls.append((sample_name, input_sample))

    logger.info(f"\nCalling peaks for {len(peak_calls)} IP samples\n")

    successful = 0
    failed = 0

    for ip_sample, input_sample in peak_calls:
        ip_bam = ALIGNED_DIR / f"{ip_sample}_Aligned.sortedByCoord.dedup.bam"
        input_bam = ALIGNED_DIR / f"{input_sample}_Aligned.sortedByCoord.dedup.bam"

        if not ip_bam.exists():
            logger.warning(f"IP BAM not found: {ip_bam}")
            failed += 1
            continue

        if not input_bam.exists():
            logger.warning(f"Input BAM not found: {input_bam}")
            failed += 1
            continue

        logger.info("\n" + "="*70)
        logger.info(f"Processing: {ip_sample}")
        logger.info("="*70)

        # Step 1: Call peaks with CLIPper
        peaks_file = call_peaks_clipper(
            ip_bam,
            PEAKS_DIR,
            ip_sample,
            species=CLIPPER_SPECIES,
            fdr=CLIPPER_FDR
        )

        if not peaks_file:
            logger.error(f"Peak calling failed for {ip_sample}")
            failed += 1
            continue

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
            logger.info(f"✓ Peak calling complete for {ip_sample}")
            successful += 1
        else:
            logger.error(f"✗ Normalization failed for {ip_sample}")
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


if __name__ == '__main__':
    main()
