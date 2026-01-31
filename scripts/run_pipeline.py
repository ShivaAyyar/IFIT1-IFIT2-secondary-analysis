#!/usr/bin/env python
"""
run_pipeline.py - Master script to run the entire eCLIP analysis pipeline

This script orchestrates all steps of the analysis pipeline.

Usage:
    # Run complete pipeline
    python scripts/run_pipeline.py

    # Skip reference setup
    python scripts/run_pipeline.py --skip-reference

    # Skip data download
    python scripts/run_pipeline.py --skip-download
"""

import argparse
import logging
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import SAMPLE_INFO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def run_script(script_path, args=None):
    """Run a Python script and check for errors."""
    cmd = [sys.executable, str(script_path)]
    if args:
        cmd.extend(args)

    logger.info(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd)

    if result.returncode != 0:
        logger.error(f"Script failed: {script_path}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Run complete IFIT2-IFIT3 eCLIP analysis pipeline'
    )
    parser.add_argument(
        '--skip-reference',
        action='store_true',
        help='Skip reference genome setup'
    )
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip data download'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads (default: 8)'
    )

    args = parser.parse_args()

    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Complete Pipeline")
    logger.info("="*70)

    scripts_dir = Path(__file__).parent

    # Step 1: Setup reference genome
    if not args.skip_reference:
        logger.info("\n" + "="*70)
        logger.info("STEP 1: Setup Reference Genome")
        logger.info("="*70)

        if not run_script(
            scripts_dir / '01_setup_reference.py',
            ['--threads', str(args.threads)]
        ):
            sys.exit(1)
    else:
        logger.info("\nSkipping reference setup")

    # Step 2: Download data
    if not args.skip_download:
        logger.info("\n" + "="*70)
        logger.info("STEP 2: Download Data")
        logger.info("="*70)

        if not run_script(
            scripts_dir / '02_download_data.py',
            ['--all', '--threads', str(args.threads)]
        ):
            sys.exit(1)
    else:
        logger.info("\nSkipping data download")

    # Step 3: Process all samples
    logger.info("\n" + "="*70)
    logger.info("STEP 3: Process Samples")
    logger.info("="*70)

    for sample_name in SAMPLE_INFO.keys():
        if not run_script(
            scripts_dir / '03_process_samples.py',
            ['--sample', sample_name, '--threads', str(args.threads)]
        ):
            logger.warning(f"Failed to process {sample_name}, continuing...")

    # Step 4: Call peaks
    logger.info("\n" + "="*70)
    logger.info("STEP 4: Call Peaks")
    logger.info("="*70)

    if not run_script(scripts_dir / '04_call_peaks.py'):
        logger.error("Peak calling failed")
        sys.exit(1)

    # Step 5: Analyze UTRs
    logger.info("\n" + "="*70)
    logger.info("STEP 5: Analyze 5' UTRs")
    logger.info("="*70)

    if not run_script(scripts_dir / '05_analyze_utrs.py'):
        logger.error("UTR analysis failed")
        sys.exit(1)

    # Step 6: Generate figures
    logger.info("\n" + "="*70)
    logger.info("STEP 6: Generate Figures")
    logger.info("="*70)

    if not run_script(scripts_dir / '06_visualize.py'):
        logger.error("Visualization failed")
        sys.exit(1)

    logger.info("\n" + "="*70)
    logger.info("PIPELINE COMPLETE!")
    logger.info("="*70)
    logger.info("\nAll analysis steps completed successfully.")
    logger.info("Check the results/ directory for output files and figures.")


if __name__ == '__main__':
    main()
