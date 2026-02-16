#!/usr/bin/env python
"""
09_motif_analysis.py - Motif discovery in IFIT binding sites

This script extracts sequences from peak regions and runs MEME Suite STREME
for de novo motif discovery. Also includes k-mer enrichment analysis as a
Python-native fallback.

Requirements:
    - MEME Suite (conda install -c bioconda meme)
    - bedtools (already installed)
    - biopython (pip install biopython)
    - hg19 reference genome FASTA

Usage:
    python scripts/09_motif_analysis.py
"""

import logging
import subprocess
import sys
from pathlib import Path
from collections import Counter
import itertools

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))

from config import PEAKS_DIR, RESULTS_DIR, FIGURES_DIR, REFERENCE_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/motif_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Output directory
MOTIF_DIR = RESULTS_DIR / 'motifs'
MOTIF_DIR.mkdir(parents=True, exist_ok=True)


def check_tools_installed():
    """Check if required tools are installed."""
    tools = {}

    # Check bedtools
    try:
        subprocess.run(['bedtools', '--version'], capture_output=True, check=True)
        tools['bedtools'] = True
    except (FileNotFoundError, subprocess.CalledProcessError):
        tools['bedtools'] = False

    # Check STREME (MEME Suite)
    try:
        subprocess.run(['streme', '--version'], capture_output=True, check=True)
        tools['streme'] = True
    except (FileNotFoundError, subprocess.CalledProcessError):
        tools['streme'] = False

    return tools


def extract_peak_sequences(peaks_bed, genome_fasta, output_fasta, extend=0):
    """
    Extract sequences from peak regions using bedtools getfasta.

    Parameters:
    -----------
    peaks_bed : Path
        Input BED file with peaks
    genome_fasta : Path
        Reference genome FASTA file
    output_fasta : Path
        Output FASTA file
    extend : int
        Extend peaks by this many bp on each side

    Returns:
    --------
    Path to output FASTA or None if failed
    """
    logger.info(f"Extracting sequences from {peaks_bed.name}")

    # If extending, create temporary extended BED
    if extend > 0:
        extended_bed = output_fasta.with_suffix('.extended.bed')
        cmd_extend = [
            'bedtools', 'slop',
            '-i', str(peaks_bed),
            '-g', str(REFERENCE_DIR / 'hg19.chrom.sizes'),
            '-b', str(extend)
        ]
        with open(extended_bed, 'w') as f:
            subprocess.run(cmd_extend, stdout=f, check=True)
        bed_to_use = extended_bed
    else:
        bed_to_use = peaks_bed

    # Extract sequences with strand information
    cmd = [
        'bedtools', 'getfasta',
        '-fi', str(genome_fasta),
        '-bed', str(bed_to_use),
        '-fo', str(output_fasta),
        '-s',  # Use strand information
        '-name'  # Use BED name field in FASTA header
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Extracted sequences to {output_fasta}")

        # Count sequences
        with open(output_fasta) as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        logger.info(f"  Total sequences: {seq_count}")

        return output_fasta

    except subprocess.CalledProcessError as e:
        logger.error(f"  bedtools getfasta failed: {e.stderr}")
        return None


def run_streme(fasta_file, output_dir, min_width=4, max_width=15):
    """
    Run MEME Suite STREME for de novo motif discovery.

    Parameters:
    -----------
    fasta_file : Path
        Input FASTA file with sequences
    output_dir : Path
        Output directory for STREME results
    min_width : int
        Minimum motif width
    max_width : int
        Maximum motif width

    Returns:
    --------
    Path to STREME output directory or None if failed
    """
    logger.info(f"Running STREME on {fasta_file.name}")

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        'streme',
        '--p', str(fasta_file),
        '--oc', str(output_dir),
        '--rna',  # RNA alphabet (use U instead of T)
        '--minw', str(min_width),
        '--maxw', str(max_width),
        '--thresh', '0.05',  # E-value threshold
        '--nmotifs', '10'  # Maximum number of motifs
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  STREME complete. Results in {output_dir}")
        return output_dir
    except subprocess.CalledProcessError as e:
        logger.error(f"  STREME failed: {e.stderr}")
        return None


def convert_to_rna(fasta_file, output_file):
    """Convert DNA FASTA to RNA (T -> U)."""
    with open(fasta_file) as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(line)
            else:
                f_out.write(line.upper().replace('T', 'U'))
    return output_file


def analyze_kmer_enrichment(fasta_file, k=6):
    """
    Analyze k-mer frequencies in sequences.

    Parameters:
    -----------
    fasta_file : Path
        Input FASTA file
    k : int
        K-mer length

    Returns:
    --------
    DataFrame with k-mer frequencies
    """
    logger.info(f"Analyzing {k}-mer frequencies in {fasta_file.name}")

    # Read sequences
    sequences = []
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_seq:
            sequences.append(''.join(current_seq))

    # Count k-mers
    kmer_counts = Counter()
    total_kmers = 0

    for seq in sequences:
        seq = seq.replace('T', 'U')  # Convert to RNA
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(c in 'ACGU' for c in kmer):  # Only count valid RNA k-mers
                kmer_counts[kmer] += 1
                total_kmers += 1

    # Create DataFrame
    kmer_df = pd.DataFrame([
        {'kmer': kmer, 'count': count, 'frequency': count / total_kmers}
        for kmer, count in kmer_counts.most_common()
    ])

    logger.info(f"  Found {len(kmer_df)} unique {k}-mers")

    return kmer_df


def analyze_gc_content(fasta_file):
    """
    Analyze GC content and nucleotide composition.

    Parameters:
    -----------
    fasta_file : Path
        Input FASTA file

    Returns:
    --------
    dict with composition statistics
    """
    logger.info(f"Analyzing nucleotide composition in {fasta_file.name}")

    # Read sequences
    all_seq = []
    gc_per_seq = []

    with open(fasta_file) as f:
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    seq = ''.join(current_seq).upper()
                    all_seq.append(seq)
                    gc = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
                    gc_per_seq.append(gc)
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            seq = ''.join(current_seq).upper()
            all_seq.append(seq)
            gc = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
            gc_per_seq.append(gc)

    # Overall composition
    combined = ''.join(all_seq)
    total = len(combined)

    composition = {
        'A': combined.count('A') / total,
        'C': combined.count('C') / total,
        'G': combined.count('G') / total,
        'T': combined.count('T') / total,
        'U': combined.count('U') / total,
    }

    stats = {
        'total_sequences': len(all_seq),
        'total_nucleotides': total,
        'gc_content': (composition['G'] + composition['C']),
        'gc_mean': np.mean(gc_per_seq),
        'gc_std': np.std(gc_per_seq),
        'composition': composition,
        'gc_per_sequence': gc_per_seq
    }

    logger.info(f"  Total sequences: {stats['total_sequences']}")
    logger.info(f"  GC content: {stats['gc_content']*100:.1f}%")

    return stats


def plot_kmer_enrichment(kmer_df, output_file, top_n=30, title='K-mer Frequencies'):
    """
    Create bar plot of top k-mer frequencies.
    """
    if kmer_df is None or len(kmer_df) == 0:
        return

    df = kmer_df.head(top_n)

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.bar(range(len(df)), df['frequency'], color='steelblue')
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['kmer'], rotation=90, fontsize=8)
    ax.set_xlabel('K-mer')
    ax.set_ylabel('Frequency')
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved k-mer plot to {output_file}")


def plot_gc_distribution(gc_stats_dict, output_file):
    """
    Plot GC content distribution across samples.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Bar plot of mean GC content
    ax = axes[0]
    samples = list(gc_stats_dict.keys())
    gc_means = [gc_stats_dict[s]['gc_mean'] * 100 for s in samples]
    gc_stds = [gc_stats_dict[s]['gc_std'] * 100 for s in samples]

    ax.bar(samples, gc_means, yerr=gc_stds, capsize=5, color='steelblue')
    ax.set_ylabel('GC Content (%)')
    ax.set_title('Mean GC Content by Sample')
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.axhline(50, color='gray', linestyle='--', alpha=0.5)

    # Histogram of GC distribution for each sample
    ax = axes[1]
    for sample, stats in gc_stats_dict.items():
        gc_values = [x * 100 for x in stats['gc_per_sequence']]
        ax.hist(gc_values, bins=30, alpha=0.5, label=sample, density=True)

    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Density')
    ax.set_title('GC Content Distribution')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved GC distribution plot to {output_file}")


def plot_kmer_comparison(kmer_dfs_dict, output_file, k=6, top_n=20):
    """
    Create heatmap comparing k-mer frequencies across samples.
    """
    # Get union of top k-mers from all samples
    all_top_kmers = set()
    for df in kmer_dfs_dict.values():
        if df is not None and len(df) > 0:
            all_top_kmers.update(df.head(top_n)['kmer'].tolist())

    if len(all_top_kmers) == 0:
        return

    # Build matrix
    samples = list(kmer_dfs_dict.keys())
    matrix_data = []

    for kmer in all_top_kmers:
        row = {'kmer': kmer}
        for sample in samples:
            df = kmer_dfs_dict[sample]
            if df is not None and len(df) > 0:
                kmer_row = df[df['kmer'] == kmer]
                if len(kmer_row) > 0:
                    row[sample] = kmer_row['frequency'].values[0]
                else:
                    row[sample] = 0
            else:
                row[sample] = 0
        matrix_data.append(row)

    matrix_df = pd.DataFrame(matrix_data).set_index('kmer')

    # Sort by variance (most different k-mers first)
    matrix_df['variance'] = matrix_df.var(axis=1)
    matrix_df = matrix_df.sort_values('variance', ascending=False).drop('variance', axis=1)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, max(6, len(matrix_df) * 0.3)))

    sns.heatmap(
        matrix_df.head(30),
        cmap='YlOrRd',
        annot=False,
        linewidths=0.5,
        ax=ax,
        cbar_kws={'label': 'Frequency'}
    )

    ax.set_title(f'{k}-mer Frequency Comparison')
    ax.set_xlabel('Sample')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved k-mer comparison heatmap to {output_file}")


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Motif Discovery")
    logger.info("="*70)

    # Check tools
    tools = check_tools_installed()
    logger.info(f"Tools available: {tools}")

    if not tools['bedtools']:
        logger.error("bedtools is required. Install with: conda install -c bioconda bedtools")
        sys.exit(1)

    # Check for genome FASTA
    genome_fasta = REFERENCE_DIR / 'hg19.fa'
    genome_fasta_gz = REFERENCE_DIR / 'hg19.fa.gz'

    if not genome_fasta.exists() and genome_fasta_gz.exists():
        # bedtools cannot use gzip-compressed FASTA, need to decompress
        logger.info(f"Decompressing {genome_fasta_gz.name} for bedtools...")
        logger.info("  (This may take a few minutes and ~3GB disk space)")

        import gzip
        import shutil

        try:
            with gzip.open(genome_fasta_gz, 'rb') as f_in:
                with open(genome_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"  Decompressed to {genome_fasta}")
        except Exception as e:
            logger.error(f"  Failed to decompress genome: {e}")
            logger.error("  You can manually decompress with: gunzip -k data/reference/hg19.fa.gz")
            sys.exit(1)

    if not genome_fasta.exists():
        logger.error(f"Genome FASTA not found at {REFERENCE_DIR}")
        logger.error("Please download hg19.fa to the reference directory")
        logger.error("  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz")
        logger.error("  gunzip hg19.fa.gz")
        sys.exit(1)

    # Define samples
    samples = {
        'Complex_FLAG': 'IFIT2_IFIT3_FLAG_IP',
        'Complex_HA': 'IFIT2_IFIT3_HA_IP',
        'IFIT2_only': 'IFIT2_alone_FLAG_IP',
        'IFIT3_only': 'IFIT3_alone_HA_IP',
    }

    # Process each sample
    gc_stats = {}
    kmer_results = {}

    for label, sample_name in samples.items():
        logger.info(f"\n{'='*50}")
        logger.info(f"Processing {label} ({sample_name})")
        logger.info("="*50)

        peaks_file = PEAKS_DIR / f"{sample_name}_normalized_peaks.bed"

        if not peaks_file.exists():
            logger.warning(f"Peak file not found: {peaks_file}")
            continue

        # Extract sequences
        fasta_file = MOTIF_DIR / f"{label}_peaks.fa"
        result = extract_peak_sequences(peaks_file, genome_fasta, fasta_file)

        if result is None:
            continue

        # Convert to RNA
        rna_fasta = MOTIF_DIR / f"{label}_peaks_rna.fa"
        convert_to_rna(fasta_file, rna_fasta)

        # Run STREME if available
        if tools['streme']:
            streme_dir = MOTIF_DIR / f"{label}_streme"
            run_streme(rna_fasta, streme_dir)

        # K-mer analysis (Python-native)
        kmer_df = analyze_kmer_enrichment(fasta_file, k=6)
        kmer_df.to_csv(MOTIF_DIR / f"{label}_6mer_frequencies.csv", index=False)
        kmer_results[label] = kmer_df

        # Plot k-mer frequencies
        plot_kmer_enrichment(
            kmer_df,
            FIGURES_DIR / f"{label}_kmer_frequencies.pdf",
            title=f"{label} - Top 6-mer Frequencies"
        )

        # GC content analysis
        gc_stats[label] = analyze_gc_content(fasta_file)

    # Save GC content summary
    gc_summary = []
    for label, stats in gc_stats.items():
        gc_summary.append({
            'Sample': label,
            'Total_Sequences': stats['total_sequences'],
            'Total_Nucleotides': stats['total_nucleotides'],
            'GC_Content': stats['gc_content'],
            'GC_Mean': stats['gc_mean'],
            'GC_Std': stats['gc_std'],
            'A_freq': stats['composition']['A'],
            'C_freq': stats['composition']['C'],
            'G_freq': stats['composition']['G'],
            'T_freq': stats['composition']['T'],
        })

    gc_df = pd.DataFrame(gc_summary)
    gc_df.to_csv(MOTIF_DIR / 'gc_content_summary.csv', index=False)

    # Create comparison plots
    if len(gc_stats) > 0:
        plot_gc_distribution(gc_stats, FIGURES_DIR / 'gc_content_distribution.pdf')

    if len(kmer_results) > 0:
        plot_kmer_comparison(kmer_results, FIGURES_DIR / 'kmer_comparison_heatmap.pdf')

    logger.info("\n" + "="*70)
    logger.info("Motif analysis complete!")
    logger.info("="*70)
    logger.info("\nResults saved to:")
    logger.info(f"  {MOTIF_DIR}/")
    if tools['streme']:
        logger.info(f"  {MOTIF_DIR}/*_streme/ (STREME motif results)")
    logger.info(f"  {FIGURES_DIR}/*_kmer_frequencies.pdf")
    logger.info(f"  {FIGURES_DIR}/gc_content_distribution.pdf")
    logger.info(f"  {FIGURES_DIR}/kmer_comparison_heatmap.pdf")


if __name__ == '__main__':
    main()
