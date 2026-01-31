#!/usr/bin/env python
"""
04_call_peaks.py - Call peaks for all IP samples

This script calls peaks for all IP samples using their corresponding input controls.

Usage:
    python scripts/04_call_peaks.py
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import SAMPLE_INFO, ALIGNED_DIR, PEAKS_DIR, PEAK_MIN_COVERAGE, PEAK_FOLD_ENRICHMENT
from analysis import call_peaks_simple

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
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Peak Calling")
    logger.info("="*70)

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

    for ip_sample, input_sample in peak_calls:
        ip_bam = ALIGNED_DIR / f"{ip_sample}_Aligned.sortedByCoord.dedup.bam"
        input_bam = ALIGNED_DIR / f"{input_sample}_Aligned.sortedByCoord.dedup.bam"

        if not ip_bam.exists():
            logger.warning(f"IP BAM not found: {ip_bam}")
            continue

        if not input_bam.exists():
            logger.warning(f"Input BAM not found: {input_bam}")
            input_bam = None

        logger.info(f"[{ip_sample}]")
        peaks_file = call_peaks_simple(
            ip_bam,
            input_bam,
            PEAKS_DIR,
            ip_sample,
            min_reads=PEAK_MIN_COVERAGE,
            fold_enrichment=PEAK_FOLD_ENRICHMENT
        )

    logger.info("\n" + "="*70)
    logger.info("Peak calling complete!")
    logger.info("="*70)


if __name__ == '__main__':
    main()
