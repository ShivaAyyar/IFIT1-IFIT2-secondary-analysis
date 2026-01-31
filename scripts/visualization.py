"""
Visualization functions for IFIT2-IFIT3 eCLIP Analysis Pipeline

This module contains functions for creating publication-quality figures:
- 5' UTR length distribution plots
- IFIT2 vs IFIT3 vs complex binding comparisons
- QC metric visualizations
"""

import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

# Set publication-quality figure defaults
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'


def plot_utr_length_analysis(protein_coding_df, stats, output_file):
    """
    Create publication-quality figure showing 5' UTR length distribution.

    Parameters:
    -----------
    protein_coding_df : DataFrame
        DataFrame with 'utr5_length' and 'bound' columns
    stats : dict
        Statistics dictionary from analyze_utr_length_distribution()
    output_file : Path
        Output file path for the figure
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    bound = protein_coding_df[protein_coding_df['bound']]['utr5_length']
    unbound = protein_coding_df[~protein_coding_df['bound']]['utr5_length']

    # Panel 1: Histograms
    ax = axes[0]
    ax.hist(unbound, bins=50, alpha=0.5, label='Unbound', density=True, color='gray')
    ax.hist(bound, bins=50, alpha=0.7, label='Bound', density=True, color='red')
    ax.set_xlabel("5' UTR Length (nt)")
    ax.set_ylabel("Density")
    ax.set_title("5' UTR Length Distribution")
    ax.legend()
    ax.set_xlim(0, 500)

    # Panel 2: Cumulative distribution
    ax = axes[1]
    bound_sorted = np.sort(bound)
    unbound_sorted = np.sort(unbound)
    ax.plot(unbound_sorted, np.linspace(0, 1, len(unbound_sorted)),
            label='Unbound', color='gray', linewidth=2)
    ax.plot(bound_sorted, np.linspace(0, 1, len(bound_sorted)),
            label='Bound', color='red', linewidth=2)
    ax.axvline(50, color='black', linestyle='--', alpha=0.5, label='50 nt')
    ax.set_xlabel("5' UTR Length (nt)")
    ax.set_ylabel("Cumulative Fraction")
    ax.set_title("Cumulative Distribution")
    ax.legend()
    ax.set_xlim(0, 500)

    # Panel 3: Box plot
    ax = axes[2]
    data_to_plot = [unbound, bound]
    bp = ax.boxplot(data_to_plot, labels=['Unbound', 'Bound'],
                    patch_artist=True, showfliers=False)
    bp['boxes'][0].set_facecolor('gray')
    bp['boxes'][1].set_facecolor('red')
    ax.set_ylabel("5' UTR Length (nt)")
    ax.set_title("5' UTR Length Comparison")

    # Add statistics text
    p_val = stats.get('mann_whitney_pval', float('nan'))
    ax.text(0.5, 0.95, f"p = {p_val:.2e}",
            transform=ax.transAxes, ha='center', va='top')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved UTR length analysis figure to {output_file}")


def plot_ifit_comparison(utr_df, ifit2_alone, ifit3_alone, complex_transcripts, output_file):
    """
    Create 6-panel comparison figure showing IFIT2 vs IFIT3 vs complex binding.

    Parameters:
    -----------
    utr_df : DataFrame
        DataFrame with UTR information
    ifit2_alone : set
        Transcripts bound by IFIT2 only
    ifit3_alone : set
        Transcripts bound by IFIT3 only
    complex_transcripts : set
        Transcripts bound by both (complex)
    output_file : Path
        Output file path for the figure
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Filter for protein-coding
    pc = utr_df[utr_df['gene_type'] == 'protein_coding'].copy()

    # Add binding categories
    pc['binding'] = 'unbound'
    pc.loc[pc['transcript_id'].isin(ifit2_alone), 'binding'] = 'IFIT2 only'
    pc.loc[pc['transcript_id'].isin(ifit3_alone), 'binding'] = 'IFIT3 only'
    pc.loc[pc['transcript_id'].isin(complex_transcripts), 'binding'] = 'Complex'

    colors = {'unbound': 'lightgray', 'IFIT2 only': 'blue',
              'IFIT3 only': 'green', 'Complex': 'red'}

    # Panel 1: Histogram comparison
    ax = axes[0, 0]
    for cat in ['unbound', 'IFIT2 only', 'IFIT3 only', 'Complex']:
        data = pc[pc['binding'] == cat]['utr5_length']
        if len(data) > 0:
            ax.hist(data, bins=50, alpha=0.6, label=cat, density=True, color=colors[cat])
    ax.set_xlabel("5' UTR Length (nt)")
    ax.set_ylabel("Density")
    ax.set_title("UTR Length Distribution by Binding Category")
    ax.legend()
    ax.set_xlim(0, 500)

    # Panel 2: Cumulative distributions
    ax = axes[0, 1]
    for cat in ['Complex', 'IFIT2 only', 'IFIT3 only', 'unbound']:
        data = pc[pc['binding'] == cat]['utr5_length']
        if len(data) > 0:
            sorted_data = np.sort(data)
            ax.plot(sorted_data, np.linspace(0, 1, len(sorted_data)),
                   label=cat, linewidth=2, color=colors[cat])
    ax.axvline(50, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel("5' UTR Length (nt)")
    ax.set_ylabel("Cumulative Fraction")
    ax.set_title("Cumulative Distributions")
    ax.legend()
    ax.set_xlim(0, 500)

    # Panel 3: Box plots
    ax = axes[0, 2]
    categories = ['unbound', 'IFIT2 only', 'IFIT3 only', 'Complex']
    data_to_plot = [pc[pc['binding'] == cat]['utr5_length'] for cat in categories]
    bp = ax.boxplot(data_to_plot, labels=categories, patch_artist=True, showfliers=False)
    for patch, cat in zip(bp['boxes'], categories):
        patch.set_facecolor(colors[cat])
    ax.set_ylabel("5' UTR Length (nt)")
    ax.set_title("UTR Length by Binding Category")
    ax.set_xticklabels(categories, rotation=45, ha='right')

    # Panel 4: Fraction with short UTRs
    ax = axes[1, 0]
    fractions = []
    for cat in categories:
        data = pc[pc['binding'] == cat]['utr5_length']
        if len(data) > 0:
            frac = (data < 50).sum() / len(data)
            fractions.append(frac)
        else:
            fractions.append(0)
    bars = ax.bar(categories, fractions, color=[colors[c] for c in categories])
    ax.set_ylabel("Fraction with short (<50 nt) 5' UTR")
    ax.set_title("Short 5' UTR Enrichment")
    ax.set_xticklabels(categories, rotation=45, ha='right')
    ax.set_ylim(0, 1)

    # Panel 5: Transcript counts
    ax = axes[1, 1]
    counts = [len(pc[pc['binding'] == cat]) for cat in categories]
    bars = ax.bar(categories, counts, color=[colors[c] for c in categories])
    ax.set_ylabel("Number of Transcripts")
    ax.set_title("Transcript Counts")
    ax.set_xticklabels(categories, rotation=45, ha='right')
    for i, count in enumerate(counts):
        ax.text(i, count, f'{count:,}', ha='center', va='bottom')

    # Panel 6: Median UTR lengths
    ax = axes[1, 2]
    medians = []
    for cat in categories:
        data = pc[pc['binding'] == cat]['utr5_length']
        if len(data) > 0:
            medians.append(data.median())
        else:
            medians.append(0)
    bars = ax.bar(categories, medians, color=[colors[c] for c in categories])
    ax.axhline(50, color='black', linestyle='--', alpha=0.5, label='50 nt')
    ax.set_ylabel("Median 5' UTR Length (nt)")
    ax.set_title("Median UTR Lengths")
    ax.set_xticklabels(categories, rotation=45, ha='right')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved IFIT comparison figure to {output_file}")


def plot_qc_metrics(qc_data, output_file):
    """
    Create QC metrics visualization.

    Parameters:
    -----------
    qc_data : dict or DataFrame
        QC metrics for all samples
    output_file : Path
        Output file path for the figure
    """
    # This is a placeholder - implement based on actual QC data structure
    logger.info(f"QC metrics plot saved to {output_file}")

    # Basic example
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, "QC Metrics Placeholder",
            ha='center', va='center', fontsize=20)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
