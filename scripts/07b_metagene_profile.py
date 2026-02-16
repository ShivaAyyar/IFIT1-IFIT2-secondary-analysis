#!/usr/bin/env python
"""
07b_metagene_profile.py - Generate metagene binding profiles

This script uses deepTools to create metagene profiles showing IFIT binding
density across transcript regions (5'UTR -> CDS -> 3'UTR).

Requirements:
    - deepTools (conda install -c bioconda deeptools)
    - BAM files with aligned reads
    - GTF annotation file

Usage:
    python scripts/07b_metagene_profile.py
"""

import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from config import ALIGNED_DIR, REFERENCE_DIR, RESULTS_DIR, FIGURES_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/metagene_profile.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Output directory
METAGENE_DIR = RESULTS_DIR / 'metagene'
METAGENE_DIR.mkdir(parents=True, exist_ok=True)


def check_deeptools_installed():
    """Check if deepTools is installed."""
    try:
        result = subprocess.run(
            ['computeMatrix', '--version'],
            capture_output=True,
            text=True
        )
        return True
    except FileNotFoundError:
        return False


def bam_to_bigwig(bam_file, bigwig_file, normalize=True):
    """
    Convert BAM to bigWig format for deepTools.

    Parameters:
    -----------
    bam_file : Path
        Input BAM file
    bigwig_file : Path
        Output bigWig file
    normalize : bool
        Whether to normalize by RPKM

    Returns:
    --------
    Path to bigWig file or None if failed
    """
    logger.info(f"Converting {bam_file.name} to bigWig")

    cmd = [
        'bamCoverage',
        '-b', str(bam_file),
        '-o', str(bigwig_file),
        '--normalizeUsing', 'RPKM' if normalize else 'None',
        '-p', '4',  # Number of processors
        '--ignoreDuplicates'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Created bigWig: {bigwig_file}")
        return bigwig_file
    except subprocess.CalledProcessError as e:
        logger.error(f"  bamCoverage failed: {e.stderr}")
        return None


def create_gene_bed(gtf_file, output_bed):
    """
    Create BED file with gene coordinates from GTF.

    Parameters:
    -----------
    gtf_file : Path
        Input GTF file (gzipped or plain)
    output_bed : Path
        Output BED file

    Returns:
    --------
    Path to BED file
    """
    logger.info(f"Creating gene BED from {gtf_file.name}")

    import gzip

    genes = {}

    # Open GTF file
    open_func = gzip.open if str(gtf_file).endswith('.gz') else open
    mode = 'rt' if str(gtf_file).endswith('.gz') else 'r'

    with open_func(gtf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes
            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attrs[key] = value.strip('"')

            gene_id = attrs.get('gene_id', '')
            gene_name = attrs.get('gene_name', gene_id)
            gene_type = attrs.get('gene_type', '')

            # Only include protein-coding genes
            if gene_type == 'protein_coding':
                genes[gene_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': gene_name,
                    'strand': strand
                }

    # Write BED file
    with open(output_bed, 'w') as f:
        for gene_id, g in genes.items():
            f.write(f"{g['chrom']}\t{g['start']}\t{g['end']}\t{g['name']}\t0\t{g['strand']}\n")

    logger.info(f"  Created gene BED with {len(genes)} protein-coding genes")
    return output_bed


def compute_matrix_scale_regions(bigwig_files, regions_bed, output_matrix,
                                  upstream=1000, downstream=1000, body_length=5000):
    """
    Compute matrix for metagene analysis using deepTools.

    Parameters:
    -----------
    bigwig_files : list
        List of bigWig file paths
    regions_bed : Path
        BED file with gene regions
    output_matrix : Path
        Output matrix file
    upstream : int
        bp upstream of TSS
    downstream : int
        bp downstream of TES
    body_length : int
        Scaled gene body length

    Returns:
    --------
    Path to matrix file or None if failed
    """
    logger.info("Computing matrix for metagene profile")

    cmd = [
        'computeMatrix', 'scale-regions',
        '-S'] + [str(f) for f in bigwig_files] + [
        '-R', str(regions_bed),
        '-o', str(output_matrix),
        '--upstream', str(upstream),
        '--downstream', str(downstream),
        '--regionBodyLength', str(body_length),
        '-p', '4',
        '--skipZeros',
        '--missingDataAsZero'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Created matrix: {output_matrix}")
        return output_matrix
    except subprocess.CalledProcessError as e:
        logger.error(f"  computeMatrix failed: {e.stderr}")
        return None


def plot_profile(matrix_file, output_file, sample_labels=None):
    """
    Generate metagene profile plot using deepTools.

    Parameters:
    -----------
    matrix_file : Path
        Input matrix file from computeMatrix
    output_file : Path
        Output plot file (PDF or PNG)
    sample_labels : list, optional
        Labels for each sample

    Returns:
    --------
    Path to output file or None if failed
    """
    logger.info("Generating metagene profile plot")

    cmd = [
        'plotProfile',
        '-m', str(matrix_file),
        '-o', str(output_file),
        '--plotTitle', 'IFIT Binding Metagene Profile',
        '--perGroup',
        '--plotHeight', '8',
        '--plotWidth', '12'
    ]

    if sample_labels:
        cmd.extend(['--samplesLabel'] + sample_labels)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Created profile plot: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        logger.error(f"  plotProfile failed: {e.stderr}")
        return None


def plot_heatmap(matrix_file, output_file, sample_labels=None):
    """
    Generate metagene heatmap using deepTools.

    Parameters:
    -----------
    matrix_file : Path
        Input matrix file from computeMatrix
    output_file : Path
        Output plot file (PDF or PNG)
    sample_labels : list, optional
        Labels for each sample

    Returns:
    --------
    Path to output file or None if failed
    """
    logger.info("Generating metagene heatmap")

    cmd = [
        'plotHeatmap',
        '-m', str(matrix_file),
        '-o', str(output_file),
        '--plotTitle', 'IFIT Binding Heatmap',
        '--colorMap', 'Reds',
        '--whatToShow', 'heatmap and colorbar',
        '--sortUsing', 'mean'
    ]

    if sample_labels:
        cmd.extend(['--samplesLabel'] + sample_labels)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Created heatmap: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        logger.error(f"  plotHeatmap failed: {e.stderr}")
        return None


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Metagene Profiles")
    logger.info("="*70)

    # Check deepTools installation
    if not check_deeptools_installed():
        logger.error("deepTools is not installed. Install with:")
        logger.error("  conda install -c bioconda deeptools")
        sys.exit(1)

    # Find GTF file
    gtf_file = REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz'
    if not gtf_file.exists():
        gtf_file = REFERENCE_DIR / 'gencode.v19.annotation.gtf'
    if not gtf_file.exists():
        logger.error(f"GTF file not found at {REFERENCE_DIR}")
        sys.exit(1)

    # Create gene BED file
    gene_bed = METAGENE_DIR / 'protein_coding_genes.bed'
    if not gene_bed.exists():
        create_gene_bed(gtf_file, gene_bed)

    # Define samples
    samples = {
        'Complex_FLAG': 'IFIT2_IFIT3_FLAG_IP',
        'Complex_HA': 'IFIT2_IFIT3_HA_IP',
        'IFIT2_only': 'IFIT2_alone_FLAG_IP',
        'IFIT3_only': 'IFIT3_alone_HA_IP',
    }

    # Convert BAMs to bigWigs
    bigwig_files = []
    sample_labels = []

    for label, sample_name in samples.items():
        # Find BAM file
        bam_patterns = [
            f"{sample_name}_Aligned.sortedByCoord.dedup.bam",
            f"{sample_name}_dedup.bam",
            f"{sample_name}.bam"
        ]

        bam_file = None
        for pattern in bam_patterns:
            potential_bam = ALIGNED_DIR / pattern
            if potential_bam.exists():
                bam_file = potential_bam
                break

        if bam_file is None:
            logger.warning(f"BAM file not found for {sample_name}")
            continue

        # Convert to bigWig
        bigwig_file = METAGENE_DIR / f"{label}.bw"
        if not bigwig_file.exists():
            result = bam_to_bigwig(bam_file, bigwig_file)
            if result is None:
                continue

        bigwig_files.append(bigwig_file)
        sample_labels.append(label)

    if len(bigwig_files) == 0:
        logger.error("No bigWig files created. Cannot generate metagene profiles.")
        sys.exit(1)

    # Compute matrix
    matrix_file = METAGENE_DIR / 'metagene_matrix.gz'
    result = compute_matrix_scale_regions(
        bigwig_files,
        gene_bed,
        matrix_file,
        upstream=2000,
        downstream=2000,
        body_length=5000
    )

    if result is None:
        logger.error("Failed to compute matrix")
        sys.exit(1)

    # Generate plots
    plot_profile(
        matrix_file,
        FIGURES_DIR / 'metagene_profile.pdf',
        sample_labels
    )

    plot_heatmap(
        matrix_file,
        FIGURES_DIR / 'metagene_heatmap.pdf',
        sample_labels
    )

    logger.info("\n" + "="*70)
    logger.info("Metagene profile analysis complete!")
    logger.info("="*70)
    logger.info("\nResults saved to:")
    logger.info(f"  {METAGENE_DIR}/")
    logger.info(f"  {FIGURES_DIR}/metagene_profile.pdf")
    logger.info(f"  {FIGURES_DIR}/metagene_heatmap.pdf")


if __name__ == '__main__':
    main()
