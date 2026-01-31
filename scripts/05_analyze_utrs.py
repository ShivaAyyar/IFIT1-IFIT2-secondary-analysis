#!/usr/bin/env python
"""
05_analyze_utrs.py - Analyze 5' UTR binding patterns

This script:
1. Parses GENCODE GTF for 5' UTRs
2. Identifies transcripts bound by IFIT2, IFIT3, and complex
3. Analyzes UTR length distributions
4. Saves results to CSV files

Usage:
    python scripts/05_analyze_utrs.py
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import REFERENCE_DIR, PEAKS_DIR, RESULTS_DIR
from analysis import (
    parse_gencode_gtf_for_utrs,
    create_utr5_bed,
    identify_bound_transcripts,
    identify_complex_bound_transcripts,
    analyze_utr_length_distribution
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/utr_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: 5' UTR Analysis")
    logger.info("="*70)

    # Parse GENCODE GTF
    logger.info("\nStep 1: Parsing GENCODE GTF for 5' UTRs...")
    gtf_file = REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz'

    if not gtf_file.exists():
        logger.error(f"GTF file not found: {gtf_file}")
        sys.exit(1)

    utr_df = parse_gencode_gtf_for_utrs(gtf_file)

    # Create BED file
    logger.info("\nStep 2: Creating 5' UTR BED file...")
    utr_bed = REFERENCE_DIR / '5utr_regions.bed'
    create_utr5_bed(utr_df, utr_bed)

    # Analyze IFIT2+IFIT3 complex binding
    logger.info("\nStep 3: Analyzing complex binding...")

    # Use normalized peaks from CLIPper
    flag_peaks = PEAKS_DIR / 'IFIT2_IFIT3_FLAG_IP_normalized_peaks.bed'
    ha_peaks = PEAKS_DIR / 'IFIT2_IFIT3_HA_IP_normalized_peaks.bed'

    if not flag_peaks.exists() or not ha_peaks.exists():
        logger.error("Normalized peak files not found. Run 04_call_peaks.py first")
        logger.error(f"Looking for: {flag_peaks}")
        logger.error(f"Looking for: {ha_peaks}")
        sys.exit(1)

    complex_results = identify_complex_bound_transcripts(
        flag_peaks,
        ha_peaks,
        utr_bed
    )

    # Analyze IFIT2 alone
    logger.info("\nStep 4: Analyzing IFIT2 alone binding...")
    ifit2_alone_peaks = PEAKS_DIR / 'IFIT2_alone_FLAG_IP_normalized_peaks.bed'
    ifit2_alone_bound = set()

    if ifit2_alone_peaks.exists():
        ifit2_alone_bound = identify_bound_transcripts(
            ifit2_alone_peaks,
            utr_bed,
            utr_df
        )

    # Analyze IFIT3 alone
    logger.info("\nStep 5: Analyzing IFIT3 alone binding...")
    ifit3_alone_peaks = PEAKS_DIR / 'IFIT3_alone_HA_IP_normalized_peaks.bed'
    ifit3_alone_bound = set()

    if ifit3_alone_peaks.exists():
        ifit3_alone_bound = identify_bound_transcripts(
            ifit3_alone_peaks,
            utr_bed,
            utr_df
        )

    # Save bound transcript lists
    logger.info("\nStep 6: Saving bound transcript lists...")

    (RESULTS_DIR / 'ifit2_alone_bound_transcripts.txt').write_text(
        '\n'.join(sorted(ifit2_alone_bound))
    )
    (RESULTS_DIR / 'ifit3_alone_bound_transcripts.txt').write_text(
        '\n'.join(sorted(ifit3_alone_bound))
    )
    (RESULTS_DIR / 'complex_bound_transcripts.txt').write_text(
        '\n'.join(sorted(complex_results['complex']))
    )

    # Analyze UTR length distributions
    logger.info("\nStep 7: Analyzing UTR length distributions...")

    # Complex binding
    pc_complex, stats_complex = analyze_utr_length_distribution(
        utr_df,
        complex_results['complex']
    )

    # Save UTR length distribution data
    pc_complex.to_csv(RESULTS_DIR / 'utr_length_distribution.csv', index=False)

    logger.info("\n" + "="*70)
    logger.info("UTR analysis complete!")
    logger.info("="*70)
    logger.info("\nResults saved to:")
    logger.info(f"  {RESULTS_DIR}/ifit2_alone_bound_transcripts.txt")
    logger.info(f"  {RESULTS_DIR}/ifit3_alone_bound_transcripts.txt")
    logger.info(f"  {RESULTS_DIR}/complex_bound_transcripts.txt")
    logger.info(f"  {RESULTS_DIR}/utr_length_distribution.csv")


if __name__ == '__main__':
    main()
