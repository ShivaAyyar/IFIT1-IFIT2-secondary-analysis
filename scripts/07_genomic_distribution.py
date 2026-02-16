#!/usr/bin/env python
"""
07_genomic_distribution.py - Analyze genomic distribution of IFIT binding sites

This script uses HOMER annotatePeaks.pl to annotate peaks with genomic features
(5'UTR, 3'UTR, CDS, intron, intergenic, etc.) and visualizes the distribution.

Requirements:
    - HOMER (conda install -c bioconda homer)
    - Peaks must be in BED format

Usage:
    python scripts/07_genomic_distribution.py
"""

import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))

from config import PEAKS_DIR, RESULTS_DIR, FIGURES_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/genomic_distribution.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Output directories
GENOMIC_DIST_DIR = RESULTS_DIR / 'genomic_distribution'
GENOMIC_DIST_DIR.mkdir(parents=True, exist_ok=True)


def check_homer_installed():
    """Check if HOMER is installed and accessible."""
    try:
        result = subprocess.run(
            ['annotatePeaks.pl', '-h'],
            capture_output=True,
            text=True
        )
        return True
    except FileNotFoundError:
        return False


def annotate_peaks_with_homer(peaks_bed, output_file, genome='hg19'):
    """
    Annotate peaks using HOMER annotatePeaks.pl.

    Parameters:
    -----------
    peaks_bed : Path
        Input BED file with peaks
    output_file : Path
        Output annotation file
    genome : str
        Genome build (default: hg19)

    Returns:
    --------
    Path to annotation file or None if failed
    """
    logger.info(f"Annotating peaks with HOMER: {peaks_bed.name}")

    cmd = [
        'annotatePeaks.pl',
        str(peaks_bed),
        genome,
        '-annStats', str(output_file.with_suffix('.stats.txt'))
    ]

    try:
        with open(output_file, 'w') as outf:
            result = subprocess.run(
                cmd,
                stdout=outf,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        logger.info(f"  Annotation complete: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        logger.error(f"  HOMER annotation failed: {e.stderr}")
        return None
    except FileNotFoundError:
        logger.error("  HOMER not found. Install with: conda install -c bioconda homer")
        return None


def parse_homer_annotation(annotation_file):
    """
    Parse HOMER annotation output to extract feature distribution.

    Parameters:
    -----------
    annotation_file : Path
        HOMER annotatePeaks.pl output file

    Returns:
    --------
    DataFrame with peak annotations
    """
    # HOMER output is tab-delimited with specific columns
    df = pd.read_csv(annotation_file, sep='\t', comment='#')

    # Rename columns for clarity
    if len(df.columns) >= 8:
        # Standard HOMER columns
        df.columns = ['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak_Score',
                      'Focus_Ratio', 'Annotation'] + list(df.columns[8:])

    return df


def extract_feature_type(annotation_string):
    """
    Extract simplified feature type from HOMER annotation string.

    HOMER annotations look like:
    - "Intergenic"
    - "intron (NM_001234, intron 3 of 10)"
    - "exon (NM_001234, exon 2 of 10)"
    - "promoter-TSS (NM_001234)"
    - "5' UTR (NM_001234)"
    - "3' UTR (NM_001234)"
    - "TTS (NM_001234)"
    """
    ann = str(annotation_string).lower()

    if "5' utr" in ann or "5'utr" in ann:
        return "5' UTR"
    elif "3' utr" in ann or "3'utr" in ann:
        return "3' UTR"
    elif "promoter" in ann or "tss" in ann:
        return "Promoter/TSS"
    elif "exon" in ann:
        return "Exon (CDS)"
    elif "intron" in ann:
        return "Intron"
    elif "tts" in ann:
        return "TTS"
    elif "intergenic" in ann:
        return "Intergenic"
    else:
        return "Other"


def calculate_feature_distribution(annotation_df):
    """
    Calculate the distribution of peaks across genomic features.

    Parameters:
    -----------
    annotation_df : DataFrame
        Parsed HOMER annotation

    Returns:
    --------
    DataFrame with feature counts and percentages
    """
    # Extract simplified feature types
    annotation_df['Feature'] = annotation_df['Annotation'].apply(extract_feature_type)

    # Count features
    feature_counts = annotation_df['Feature'].value_counts()

    # Create summary DataFrame
    summary = pd.DataFrame({
        'Feature': feature_counts.index,
        'Count': feature_counts.values,
        'Percentage': (feature_counts.values / len(annotation_df) * 100).round(2)
    })

    return summary


def plot_genomic_distribution(distributions_dict, output_file):
    """
    Create publication-quality bar plot of genomic feature distributions.

    Parameters:
    -----------
    distributions_dict : dict
        Dictionary mapping sample names to distribution DataFrames
    output_file : Path
        Output file path for the figure
    """
    # Define feature order and colors
    feature_order = ["5' UTR", "Exon (CDS)", "3' UTR", "Intron",
                     "Promoter/TSS", "TTS", "Intergenic", "Other"]
    colors = {
        "5' UTR": "#E41A1C",
        "3' UTR": "#377EB8",
        "Exon (CDS)": "#4DAF4A",
        "Intron": "#984EA3",
        "Promoter/TSS": "#FF7F00",
        "TTS": "#FFFF33",
        "Intergenic": "#A65628",
        "Other": "#999999"
    }

    # Prepare data for plotting
    plot_data = []
    for sample_name, dist_df in distributions_dict.items():
        for _, row in dist_df.iterrows():
            plot_data.append({
                'Sample': sample_name,
                'Feature': row['Feature'],
                'Percentage': row['Percentage']
            })

    plot_df = pd.DataFrame(plot_data)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create grouped bar plot
    samples = list(distributions_dict.keys())
    x = range(len(samples))
    width = 0.1

    for i, feature in enumerate(feature_order):
        feature_data = plot_df[plot_df['Feature'] == feature]
        if len(feature_data) > 0:
            values = []
            for sample in samples:
                sample_data = feature_data[feature_data['Sample'] == sample]
                if len(sample_data) > 0:
                    values.append(sample_data['Percentage'].values[0])
                else:
                    values.append(0)

            offset = (i - len(feature_order)/2) * width
            ax.bar([xi + offset for xi in x], values, width,
                   label=feature, color=colors.get(feature, '#999999'))

    ax.set_xlabel('Sample')
    ax.set_ylabel('Percentage of Peaks (%)')
    ax.set_title('Genomic Distribution of IFIT Binding Sites')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved genomic distribution plot to {output_file}")


def plot_stacked_distribution(distributions_dict, output_file):
    """
    Create stacked bar plot of genomic feature distributions.
    """
    feature_order = ["5' UTR", "Exon (CDS)", "3' UTR", "Intron",
                     "Promoter/TSS", "TTS", "Intergenic", "Other"]
    colors = ["#E41A1C", "#4DAF4A", "#377EB8", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#999999"]

    # Build matrix for stacking
    samples = list(distributions_dict.keys())
    matrix = []

    for feature in feature_order:
        row = []
        for sample in samples:
            dist_df = distributions_dict[sample]
            feature_row = dist_df[dist_df['Feature'] == feature]
            if len(feature_row) > 0:
                row.append(feature_row['Percentage'].values[0])
            else:
                row.append(0)
        matrix.append(row)

    # Create stacked bar plot
    fig, ax = plt.subplots(figsize=(10, 6))

    bottom = [0] * len(samples)
    for i, (feature, row) in enumerate(zip(feature_order, matrix)):
        ax.bar(samples, row, bottom=bottom, label=feature, color=colors[i])
        bottom = [b + r for b, r in zip(bottom, row)]

    ax.set_xlabel('Sample')
    ax.set_ylabel('Percentage of Peaks (%)')
    ax.set_title('Genomic Distribution of IFIT Binding Sites')
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved stacked distribution plot to {output_file}")


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Genomic Distribution")
    logger.info("="*70)

    # Check HOMER installation
    if not check_homer_installed():
        logger.error("HOMER is not installed. Install with:")
        logger.error("  conda install -c bioconda homer")
        logger.error("\nAlternatively, run this script on the HPC where HOMER is available.")
        sys.exit(1)

    # Define samples to analyze
    ip_samples = [
        'IFIT2_IFIT3_FLAG_IP',  # Complex - IFIT2
        'IFIT2_IFIT3_HA_IP',    # Complex - IFIT3
        'IFIT2_alone_FLAG_IP',  # IFIT2 alone
        'IFIT3_alone_HA_IP',    # IFIT3 alone
    ]

    distributions = {}
    all_annotations = []

    for sample in ip_samples:
        peaks_file = PEAKS_DIR / f"{sample}_normalized_peaks.bed"

        if not peaks_file.exists():
            logger.warning(f"Peak file not found: {peaks_file}")
            continue

        # Annotate peaks
        annotation_file = GENOMIC_DIST_DIR / f"{sample}_annotated.txt"
        result = annotate_peaks_with_homer(peaks_file, annotation_file)

        if result is None:
            continue

        # Parse annotations
        ann_df = parse_homer_annotation(annotation_file)
        ann_df['Sample'] = sample
        all_annotations.append(ann_df)

        # Calculate distribution
        dist_df = calculate_feature_distribution(ann_df)
        dist_df.to_csv(GENOMIC_DIST_DIR / f"{sample}_distribution.csv", index=False)
        distributions[sample] = dist_df

        logger.info(f"\n{sample} genomic distribution:")
        for _, row in dist_df.iterrows():
            logger.info(f"  {row['Feature']}: {row['Count']} ({row['Percentage']}%)")

    if len(distributions) == 0:
        logger.error("No samples were successfully annotated")
        sys.exit(1)

    # Combine all annotations
    combined_df = pd.concat(all_annotations, ignore_index=True)
    combined_df.to_csv(GENOMIC_DIST_DIR / 'all_peak_annotations.csv', index=False)

    # Create summary
    summary_data = []
    for sample, dist_df in distributions.items():
        for _, row in dist_df.iterrows():
            summary_data.append({
                'Sample': sample,
                'Feature': row['Feature'],
                'Count': row['Count'],
                'Percentage': row['Percentage']
            })

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(GENOMIC_DIST_DIR / 'distribution_summary.csv', index=False)

    # Create visualizations
    logger.info("\nGenerating visualizations...")

    plot_genomic_distribution(
        distributions,
        FIGURES_DIR / 'genomic_distribution_grouped.pdf'
    )

    plot_stacked_distribution(
        distributions,
        FIGURES_DIR / 'genomic_distribution_stacked.pdf'
    )

    logger.info("\n" + "="*70)
    logger.info("Genomic distribution analysis complete!")
    logger.info("="*70)
    logger.info("\nResults saved to:")
    logger.info(f"  {GENOMIC_DIST_DIR}/")
    logger.info(f"  {FIGURES_DIR}/genomic_distribution_grouped.pdf")
    logger.info(f"  {FIGURES_DIR}/genomic_distribution_stacked.pdf")


if __name__ == '__main__':
    main()
