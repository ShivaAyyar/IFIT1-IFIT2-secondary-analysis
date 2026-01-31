#!/usr/bin/env python
"""
06_visualize.py - Generate publication-quality figures

This script creates figures from the analysis results.

Usage:
    python scripts/06_visualize.py
"""

import logging
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from config import RESULTS_DIR, FIGURES_DIR
from visualization import plot_utr_length_analysis, plot_ifit_comparison

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/visualization.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Visualization")
    logger.info("="*70)

    # Load UTR data
    utr_file = RESULTS_DIR / 'utr_length_distribution.csv'
    if not utr_file.exists():
        logger.error(f"UTR data not found: {utr_file}")
        logger.error("Run 05_analyze_utrs.py first")
        sys.exit(1)

    utr_df = pd.read_csv(utr_file)
    logger.info(f"Loaded {len(utr_df)} transcripts")

    # Load bound transcript lists
    ifit2_alone = set((RESULTS_DIR / 'ifit2_alone_bound_transcripts.txt').read_text().strip().split('\n'))
    ifit3_alone = set((RESULTS_DIR / 'ifit3_alone_bound_transcripts.txt').read_text().strip().split('\n'))
    complex_bound = set((RESULTS_DIR / 'complex_bound_transcripts.txt').read_text().strip().split('\n'))

    # Generate IFIT comparison figure
    logger.info("\nGenerating IFIT comparison figure...")
    plot_ifit_comparison(
        utr_df,
        ifit2_alone,
        ifit3_alone,
        complex_bound,
        FIGURES_DIR / 'ifit_comparison.pdf'
    )

    # Generate UTR length analysis for complex
    logger.info("\nGenerating UTR length analysis figure...")

    # Recalculate stats for plotting
    pc = utr_df[utr_df['gene_type'] == 'protein_coding'].copy()
    pc['bound'] = pc['transcript_id'].isin(complex_bound)

    bound = pc[pc['bound']]['utr5_length']
    unbound = pc[~pc['bound']]['utr5_length']

    from scipy.stats import mannwhitneyu
    stat, pval = mannwhitneyu(bound, unbound, alternative='less')

    stats = {
        'bound_median': bound.median(),
        'unbound_median': unbound.median(),
        'mann_whitney_pval': pval
    }

    plot_utr_length_analysis(
        pc,
        stats,
        FIGURES_DIR / 'utr_length_histogram.pdf'
    )

    logger.info("\n" + "="*70)
    logger.info("Visualization complete!")
    logger.info("="*70)
    logger.info("\nFigures saved to:")
    logger.info(f"  {FIGURES_DIR}/ifit_comparison.pdf")
    logger.info(f"  {FIGURES_DIR}/utr_length_histogram.pdf")


if __name__ == '__main__':
    main()
