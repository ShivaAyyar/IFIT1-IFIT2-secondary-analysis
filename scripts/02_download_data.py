#!/usr/bin/env python
"""
02_download_data.py - Download FASTQ files from SRA

This script downloads eCLIP FASTQ files from the SRA database.

Usage:
    # Download all samples
    python scripts/02_download_data.py --all

    # Download specific sample
    python scripts/02_download_data.py --sample IFIT2_IFIT3_FLAG_IP

    # Download with custom threads
    python scripts/02_download_data.py --all --threads 4
"""

import argparse
import logging
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from config import SAMPLE_INFO, RAW_DIR
from utils import download_sra_fastq, check_tool_available

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/download_data.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Download eCLIP FASTQ files from SRA'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--sample',
        help='Sample name to download (e.g., IFIT2_IFIT3_FLAG_IP)'
    )
    group.add_argument(
        '--all',
        action='store_true',
        help='Download all samples'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for download (default: 4)'
    )

    args = parser.parse_args()

    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Data Download")
    logger.info("="*70)

    # Check for fasterq-dump
    if not check_tool_available('fasterq-dump'):
        logger.error("fasterq-dump not found. Please install sra-tools.")
        sys.exit(1)

    # Determine which samples to download
    if args.all:
        samples_to_download = SAMPLE_INFO
        logger.info(f"\nDownloading all {len(samples_to_download)} samples...")
    else:
        if args.sample not in SAMPLE_INFO:
            logger.error(f"Sample '{args.sample}' not found in SAMPLE_INFO")
            logger.error(f"Available samples: {', '.join(SAMPLE_INFO.keys())}")
            sys.exit(1)
        samples_to_download = {args.sample: SAMPLE_INFO[args.sample]}
        logger.info(f"\nDownloading sample: {args.sample}")

    # Download samples
    logger.info("")
    logger.info("IMPORTANT: These are SINGLE-END 75nt reads")
    logger.info("="*70)

    successful = []
    failed = []

    for sample_name, sample_info in samples_to_download.items():
        srr_accession = sample_info['srr_accession']

        logger.info(f"\n[{sample_name}]")
        logger.info(f"  SRR: {srr_accession}")
        logger.info(f"  Description: {sample_info['description']}")

        result = download_sra_fastq(
            srr_accession,
            RAW_DIR,
            threads=args.threads
        )

        if result:
            successful.append(sample_name)
        else:
            failed.append(sample_name)

    # Summary
    logger.info("\n" + "="*70)
    logger.info("Download Summary")
    logger.info("="*70)
    logger.info(f"Successful: {len(successful)}/{len(samples_to_download)}")
    logger.info(f"Failed: {len(failed)}/{len(samples_to_download)}")

    if failed:
        logger.warning(f"\nFailed samples: {', '.join(failed)}")

    logger.info("\n" + "="*70)
    logger.info("Download complete!")
    logger.info("="*70)
    logger.info("\nNext steps:")
    logger.info("  1. Process samples: python scripts/03_process_samples.py --sample <SAMPLE_NAME>")
    logger.info("")


if __name__ == '__main__':
    main()
