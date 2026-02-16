#!/usr/bin/env python
"""
07_genomic_distribution.py - Analyze genomic distribution of IFIT binding sites

This script uses HOMER annotatePeaks.pl to annotate peaks with genomic features
(5'UTR, 3'UTR, CDS, intron, intergenic, etc.) and visualizes the distribution.

Requirements:
    - HOMER (conda install -c bioconda homer)
    - HOMER genome data must be installed for hg19
    - Peaks must be in BED format

HOMER Genome Installation:
    After installing HOMER, you MUST also install the genome data:

    # Option 1: Standard installation
    perl ~/.homer/configureHomer.pl -install hg19

    # Option 2: On HPC systems with module
    module load homer
    configureHomer.pl -install hg19

    # Option 3: Specify custom HOMER path
    perl /path/to/homer/configureHomer.pl -install hg19

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

from config import PEAKS_DIR, RESULTS_DIR, FIGURES_DIR, REFERENCE_DIR

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


def check_homer_genome_installed(genome='hg19'):
    """
    Check if HOMER has the specified genome installed.

    Returns:
    --------
    bool indicating if genome is installed
    """
    try:
        # HOMER stores genomes in its data directory
        result = subprocess.run(
            ['findMotifs.pl', '-find', '/dev/null', '-list'],
            capture_output=True,
            text=True
        )
        # Check if hg19 is mentioned in available genomes
        # This is a rough check - HOMER doesn't have a clean way to list installed genomes
        return True  # Can't reliably check, so assume it's there
    except Exception:
        return True  # Assume it's there and let HOMER fail with a clear message


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

    # First check that the input file has content
    if peaks_bed.stat().st_size == 0:
        logger.error(f"  Input peak file is empty: {peaks_bed}")
        return None

    # Count peaks in input
    with open(peaks_bed, 'r') as f:
        peak_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
    logger.info(f"  Input contains {peak_count} peaks")

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

        # Log stderr (HOMER often writes progress/warnings to stderr)
        if result.stderr:
            for line in result.stderr.strip().split('\n'):
                if line.strip():
                    logger.info(f"  HOMER: {line}")

        # Check if output file has content
        if output_file.stat().st_size == 0:
            logger.error(f"  HOMER produced empty output file")
            logger.error(f"  This usually means HOMER genome data is not installed.")
            logger.error(f"  To install hg19 genome for HOMER, run:")
            logger.error(f"    perl ~/.homer/configureHomer.pl -install hg19")
            logger.error(f"  Or on HPC with module system:")
            logger.error(f"    module load homer")
            logger.error(f"    configureHomer.pl -install hg19")
            return None

        # Count output lines (excluding header)
        with open(output_file, 'r') as f:
            output_lines = sum(1 for line in f) - 1  # subtract header
        logger.info(f"  Annotation complete: {output_file} ({output_lines} annotated peaks)")

        return output_file
    except subprocess.CalledProcessError as e:
        logger.error(f"  HOMER annotation failed with exit code {e.returncode}")
        if e.stderr:
            logger.error(f"  HOMER stderr: {e.stderr}")

        # Check for common HOMER errors
        if "Can't find" in str(e.stderr) or "not found" in str(e.stderr).lower():
            logger.error(f"  HOMER genome '{genome}' may not be installed.")
            logger.error(f"  To install, run: perl ~/.homer/configureHomer.pl -install {genome}")
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
    DataFrame with peak annotations, or None if file is empty/invalid
    """
    # Check file exists and has content
    annotation_path = Path(annotation_file)
    if not annotation_path.exists():
        logger.error(f"Annotation file does not exist: {annotation_file}")
        return None

    if annotation_path.stat().st_size == 0:
        logger.error(f"Annotation file is empty: {annotation_file}")
        return None

    # Read first line to check for header
    with open(annotation_file, 'r') as f:
        first_line = f.readline()
        if not first_line.strip():
            logger.error(f"Annotation file has no content: {annotation_file}")
            return None

    try:
        # HOMER output is tab-delimited with specific columns
        df = pd.read_csv(annotation_file, sep='\t', comment='#')

        if len(df) == 0:
            logger.error(f"Annotation file has no data rows: {annotation_file}")
            return None

        # Rename columns for clarity
        if len(df.columns) >= 8:
            # Standard HOMER columns
            df.columns = ['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak_Score',
                          'Focus_Ratio', 'Annotation'] + list(df.columns[8:])

        logger.info(f"  Parsed {len(df)} annotated peaks")
        return df

    except pd.errors.EmptyDataError:
        logger.error(f"Could not parse annotation file (no columns): {annotation_file}")
        logger.error("  This usually means HOMER genome data is not installed.")
        logger.error("  To install hg19 genome for HOMER, run:")
        logger.error("    perl ~/.homer/configureHomer.pl -install hg19")
        return None
    except Exception as e:
        logger.error(f"Error parsing annotation file: {e}")
        return None


def create_feature_beds_from_gtf(gtf_file, output_dir):
    """
    Create BED files for each genomic feature type from GTF.
    This is used as a fallback when HOMER genome is not installed.

    Parameters:
    -----------
    gtf_file : Path
        GENCODE GTF file (gzipped or plain)
    output_dir : Path
        Directory to write feature BED files

    Returns:
    --------
    dict mapping feature names to BED file paths
    """
    import gzip

    logger.info(f"Creating feature BED files from GTF: {gtf_file.name}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Storage for features
    features = {
        'five_prime_utr': [],
        'three_prime_utr': [],
        'CDS': [],
        'exon': [],
        'gene': []
    }

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

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes for gene info
            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attrs[key] = value.strip('"')

            gene_type = attrs.get('gene_type', '')

            # Only include protein-coding genes
            if gene_type != 'protein_coding':
                continue

            gene_id = attrs.get('gene_id', 'unknown')
            transcript_id = attrs.get('transcript_id', gene_id)

            if feature_type in features:
                features[feature_type].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': f"{transcript_id}_{feature_type}",
                    'score': 0,
                    'strand': strand
                })

    # Write BED files
    feature_beds = {}
    for feature_name, entries in features.items():
        if len(entries) == 0:
            continue

        bed_file = output_dir / f'{feature_name}.bed'
        with open(bed_file, 'w') as f:
            for entry in entries:
                f.write(f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t"
                        f"{entry['name']}\t{entry['score']}\t{entry['strand']}\n")

        feature_beds[feature_name] = bed_file
        logger.info(f"  Created {feature_name}.bed with {len(entries)} entries")

    return feature_beds


def annotate_peaks_with_bedtools(peaks_bed, feature_beds, output_file):
    """
    Annotate peaks using bedtools intersect (fallback when HOMER unavailable).

    Parameters:
    -----------
    peaks_bed : Path
        Input BED file with peaks
    feature_beds : dict
        Dictionary mapping feature names to BED file paths
    output_file : Path
        Output annotation file

    Returns:
    --------
    DataFrame with peak annotations
    """
    logger.info(f"Annotating peaks with bedtools: {peaks_bed.name}")

    # Read peaks
    peaks_df = pd.read_csv(peaks_bed, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'name', 'score', 'strand']
                           if pd.read_csv(peaks_bed, sep='\t', header=None, nrows=1).shape[1] >= 6
                           else ['chrom', 'start', 'end', 'name', 'score'])

    # Add peak_id column
    peaks_df['PeakID'] = [f"peak_{i}" for i in range(len(peaks_df))]
    peaks_df['Annotation'] = 'Intergenic'  # Default

    # Priority order: 5'UTR > CDS > 3'UTR > exon > intron > intergenic
    priority = ['five_prime_utr', 'CDS', 'three_prime_utr', 'exon', 'gene']
    feature_labels = {
        'five_prime_utr': "5' UTR",
        'three_prime_utr': "3' UTR",
        'CDS': 'Exon (CDS)',
        'exon': 'Exon (CDS)',
        'gene': 'Intron'
    }

    # Check overlap with each feature type (in reverse priority order)
    for feature_name in reversed(priority):
        if feature_name not in feature_beds:
            continue

        feature_bed = feature_beds[feature_name]

        # Run bedtools intersect
        cmd = [
            'bedtools', 'intersect',
            '-a', str(peaks_bed),
            '-b', str(feature_bed),
            '-wa', '-u'  # Report original peak if it overlaps
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            overlapping_peaks = set()

            for line in result.stdout.strip().split('\n'):
                if line:
                    fields = line.split('\t')
                    peak_key = f"{fields[0]}:{fields[1]}-{fields[2]}"
                    overlapping_peaks.add(peak_key)

            # Update annotation for overlapping peaks
            for idx, row in peaks_df.iterrows():
                peak_key = f"{row['chrom']}:{row['start']}-{row['end']}"
                if peak_key in overlapping_peaks:
                    peaks_df.at[idx, 'Annotation'] = feature_labels.get(feature_name, feature_name)

            logger.info(f"  {feature_name}: {len(overlapping_peaks)} peaks overlap")

        except subprocess.CalledProcessError as e:
            logger.warning(f"  bedtools intersect failed for {feature_name}: {e.stderr}")

    # Identify intergenic (peaks within genes but not in UTR/CDS/exon are intronic)
    # Peaks not overlapping genes at all are intergenic
    if 'gene' in feature_beds:
        cmd = [
            'bedtools', 'intersect',
            '-a', str(peaks_bed),
            '-b', str(feature_beds['gene']),
            '-v'  # Report peaks that don't overlap
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            intergenic_peaks = set()

            for line in result.stdout.strip().split('\n'):
                if line:
                    fields = line.split('\t')
                    peak_key = f"{fields[0]}:{fields[1]}-{fields[2]}"
                    intergenic_peaks.add(peak_key)

            # Mark intergenic peaks
            for idx, row in peaks_df.iterrows():
                peak_key = f"{row['chrom']}:{row['start']}-{row['end']}"
                if peak_key in intergenic_peaks:
                    peaks_df.at[idx, 'Annotation'] = 'Intergenic'

            logger.info(f"  Intergenic: {len(intergenic_peaks)} peaks")

        except subprocess.CalledProcessError as e:
            logger.warning(f"  bedtools intersect failed for intergenic: {e.stderr}")

    # Save output
    peaks_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  Saved annotations to {output_file}")

    return peaks_df


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
    use_homer = check_homer_installed()

    if use_homer:
        logger.info("\nHOMER is installed. Will attempt HOMER annotation first.")
        logger.info("IMPORTANT: HOMER requires genome data to be installed separately.")
        logger.info("If you get empty output files, install hg19 genome with:")
        logger.info("  perl ~/.homer/configureHomer.pl -install hg19")
        logger.info("Or on HPC: configureHomer.pl -install hg19")
        logger.info("\nIf HOMER fails, will fall back to bedtools-based annotation.\n")
    else:
        logger.info("\nHOMER is not installed. Using bedtools-based annotation.")
        logger.info("(Install HOMER with: conda install -c bioconda homer)\n")

    # Prepare bedtools fallback: create feature BED files from GTF
    feature_beds = None
    gtf_file = REFERENCE_DIR / 'gencode.v19.annotation.gtf.gz'
    if not gtf_file.exists():
        gtf_file = REFERENCE_DIR / 'gencode.v19.annotation.gtf'

    if gtf_file.exists():
        feature_beds_dir = GENOMIC_DIST_DIR / 'feature_beds'
        # Only create if not already done
        if not (feature_beds_dir / 'CDS.bed').exists():
            feature_beds = create_feature_beds_from_gtf(gtf_file, feature_beds_dir)
        else:
            logger.info("Using existing feature BED files")
            feature_beds = {
                'five_prime_utr': feature_beds_dir / 'five_prime_utr.bed',
                'three_prime_utr': feature_beds_dir / 'three_prime_utr.bed',
                'CDS': feature_beds_dir / 'CDS.bed',
                'exon': feature_beds_dir / 'exon.bed',
                'gene': feature_beds_dir / 'gene.bed'
            }
            # Filter to existing files
            feature_beds = {k: v for k, v in feature_beds.items() if v.exists()}
    else:
        logger.warning(f"GTF file not found at {REFERENCE_DIR}")
        if not use_homer:
            logger.error("Cannot proceed: neither HOMER nor GTF file available")
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

        annotation_file = GENOMIC_DIST_DIR / f"{sample}_annotated.txt"
        ann_df = None

        # Try HOMER first if available
        if use_homer:
            result = annotate_peaks_with_homer(peaks_file, annotation_file)
            if result is not None:
                ann_df = parse_homer_annotation(annotation_file)

        # Fall back to bedtools if HOMER failed or unavailable
        if ann_df is None and feature_beds:
            logger.info(f"Using bedtools fallback for {sample}")
            ann_df = annotate_peaks_with_bedtools(
                peaks_file,
                feature_beds,
                annotation_file
            )

        if ann_df is None:
            logger.warning(f"  Skipping sample {sample} - annotation failed")
            continue

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
