#!/usr/bin/env python
"""
01_setup_reference.py - Download and index reference genome

This script downloads the hg19 reference genome and GENCODE v19 annotation,
then builds a STAR genome index for alignment.

Usage:
    python scripts/01_setup_reference.py [--threads 8]

Requirements:
    - wget (for downloading)
    - STAR (for indexing)
    - gunzip (for decompression)
    - ~64 GB RAM for STAR indexing
    - ~100 GB free disk space
"""

import argparse
import logging
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from config import (
    REFERENCE_DIR, REFERENCE_URLS, STAR_THREADS, STAR_SJDB_OVERHANG
)
from utils import download_reference_files, build_star_index, check_tool_available

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/setup_reference.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Download and index reference genome for eCLIP analysis'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=STAR_THREADS,
        help=f'Number of threads for STAR indexing (default: {STAR_THREADS})'
    )
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip downloading reference files (assume already downloaded)'
    )
    parser.add_argument(
        '--skip-index',
        action='store_true',
        help='Skip building STAR index (assume already built)'
    )

    args = parser.parse_args()

    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Reference Genome Setup")
    logger.info("="*70)

    # Check required tools
    logger.info("\nChecking for required tools...")
    required_tools = ['wget', 'STAR', 'gunzip']
    missing_tools = []

    for tool in required_tools:
        if check_tool_available(tool):
            logger.info(f"  {tool}: ✓ Found")
        else:
            logger.error(f"  {tool}: ✗ Missing")
            missing_tools.append(tool)

    if missing_tools:
        logger.error(f"\nMissing required tools: {', '.join(missing_tools)}")
        logger.error("Please install missing tools before proceeding.")
        sys.exit(1)

    # Step 1: Download reference files
    if not args.skip_download:
        logger.info("\n" + "="*70)
        logger.info("Step 1: Downloading reference files")
        logger.info("="*70)
        logger.info("This may take 10-30 minutes depending on your internet speed...")
        logger.info("")

        ref_files = download_reference_files(REFERENCE_URLS, REFERENCE_DIR)

        if not ref_files:
            logger.error("Failed to download reference files")
            sys.exit(1)

        logger.info("✓ Download complete!")
    else:
        logger.info("\nSkipping download (--skip-download specified)")

    # Step 2: Build STAR index
    if not args.skip_index:
        logger.info("\n" + "="*70)
        logger.info("Step 2: Building STAR genome index")
        logger.info("="*70)
        logger.info("⚠️  This requires ~30-35 GB RAM and will take 30-60 minutes...")
        logger.info("If you don't have enough RAM, this will fail.")
        logger.info("Consider using HPC or cloud instance if needed.")
        logger.info("")

        genome_fasta = REFERENCE_DIR / 'hg19.fa.gz'
        gencode_gtf = REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz'

        if not genome_fasta.exists():
            logger.error(f"Genome FASTA not found: {genome_fasta}")
            logger.error("Run without --skip-download first")
            sys.exit(1)

        if not gencode_gtf.exists():
            logger.error(f"GENCODE GTF not found: {gencode_gtf}")
            logger.error("Run without --skip-download first")
            sys.exit(1)

        star_index = build_star_index(
            genome_fasta,
            gencode_gtf,
            REFERENCE_DIR,
            threads=args.threads,
            sjdb_overhang=STAR_SJDB_OVERHANG
        )

        if not star_index:
            logger.error("Failed to build STAR index")
            sys.exit(1)

        logger.info("✓ STAR index built successfully!")
    else:
        logger.info("\nSkipping STAR indexing (--skip-index specified)")

    # Verify setup
    logger.info("\n" + "="*70)
    logger.info("Verification")
    logger.info("="*70)

    genome_fa = REFERENCE_DIR / 'hg19.fa.gz'
    gencode_gtf = REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz'
    star_index_dir = REFERENCE_DIR / 'STAR_index_hg19'

    logger.info(f"hg19.fa.gz:                      {'✓ Found' if genome_fa.exists() else '✗ Missing'}")
    logger.info(f"gencode.v19.annotation.gtf.gz:   {'✓ Found' if gencode_gtf.exists() else '✗ Missing'}")
    logger.info(f"STAR index:                      {'✓ Found' if (star_index_dir / 'SA').exists() else '✗ Missing'}")

    logger.info("\n" + "="*70)
    logger.info("Reference genome setup complete!")
    logger.info("="*70)
    logger.info("\nNext steps:")
    logger.info("  1. Download eCLIP data: python scripts/02_download_data.py")
    logger.info("  2. Process samples: python scripts/03_process_samples.py")
    logger.info("")


if __name__ == '__main__':
    main()
