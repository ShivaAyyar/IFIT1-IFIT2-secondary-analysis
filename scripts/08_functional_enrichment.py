#!/usr/bin/env python
"""
08_functional_enrichment.py - Functional enrichment analysis of IFIT-bound genes

This script uses g:Profiler to perform GO term, KEGG pathway, and Reactome
enrichment analysis on genes bound by IFIT2, IFIT3, and the IFIT2-IFIT3 complex.

Requirements:
    - gprofiler-official (pip install gprofiler-official)

Usage:
    python scripts/08_functional_enrichment.py
"""

import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))

from config import RESULTS_DIR, FIGURES_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/functional_enrichment.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Output directory
ENRICHMENT_DIR = RESULTS_DIR / 'enrichment'
ENRICHMENT_DIR.mkdir(parents=True, exist_ok=True)


def load_bound_genes(results_dir):
    """
    Load bound transcript lists and convert to gene symbols.

    Returns:
    --------
    dict mapping category names to lists of gene symbols
    """
    gene_sets = {}

    # Load UTR length distribution for gene ID to gene name mapping
    utr_file = results_dir / 'utr_length_distribution.csv'
    if not utr_file.exists():
        logger.error(f"UTR data file not found: {utr_file}")
        return None

    utr_df = pd.read_csv(utr_file)
    transcript_to_gene = dict(zip(utr_df['transcript_id'], utr_df['gene_name']))

    # Load bound transcript lists
    transcript_files = {
        'IFIT2_only': results_dir / 'ifit2_alone_bound_transcripts.txt',
        'IFIT3_only': results_dir / 'ifit3_alone_bound_transcripts.txt',
        'Complex': results_dir / 'complex_bound_transcripts.txt'
    }

    for category, filepath in transcript_files.items():
        if filepath.exists():
            transcripts = filepath.read_text().strip().split('\n')
            transcripts = [t for t in transcripts if t]  # Remove empty strings

            # Convert to gene symbols (unique)
            genes = set()
            for transcript in transcripts:
                if transcript in transcript_to_gene:
                    gene = transcript_to_gene[transcript]
                    if pd.notna(gene) and gene:
                        genes.add(gene)

            gene_sets[category] = list(genes)
            logger.info(f"Loaded {len(transcripts)} transcripts -> {len(genes)} unique genes for {category}")
        else:
            logger.warning(f"Transcript file not found: {filepath}")

    return gene_sets


def get_background_genes(results_dir):
    """
    Get background gene set (all protein-coding genes analyzed).

    Returns:
    --------
    list of gene symbols
    """
    utr_file = results_dir / 'utr_length_distribution.csv'
    utr_df = pd.read_csv(utr_file)

    # Filter to protein-coding genes
    pc_df = utr_df[utr_df['gene_type'] == 'protein_coding']

    # Get unique gene names
    background = pc_df['gene_name'].dropna().unique().tolist()
    logger.info(f"Background gene set: {len(background)} protein-coding genes")

    return background


def run_gprofiler_enrichment(gene_list, background=None, organism='hsapiens',
                             sources=None, significance_threshold=0.05):
    """
    Run g:Profiler enrichment analysis.

    Parameters:
    -----------
    gene_list : list
        List of gene symbols to analyze
    background : list, optional
        Background gene list (if None, uses all genes)
    organism : str
        Organism code (default: hsapiens)
    sources : list, optional
        Data sources to query (default: GO, KEGG, REAC)
    significance_threshold : float
        FDR threshold for significance

    Returns:
    --------
    DataFrame with enrichment results
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        logger.error("gprofiler-official not installed. Install with:")
        logger.error("  pip install gprofiler-official")
        return None

    if sources is None:
        sources = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']

    gp = GProfiler(return_dataframe=True)

    try:
        results = gp.profile(
            organism=organism,
            query=gene_list,
            background=background,
            sources=sources,
            user_threshold=significance_threshold,
            significance_threshold_method='fdr',
            no_evidences=False
        )

        if results is not None and len(results) > 0:
            logger.info(f"  Found {len(results)} significant terms")
        else:
            logger.info("  No significant enrichments found")
            results = pd.DataFrame()

        return results

    except Exception as e:
        logger.error(f"  g:Profiler query failed: {e}")
        return pd.DataFrame()


def create_enrichment_dotplot(enrichment_df, output_file, title='GO Enrichment',
                              top_n=20, source_filter=None):
    """
    Create publication-quality dotplot of enrichment results.

    Parameters:
    -----------
    enrichment_df : DataFrame
        g:Profiler results
    output_file : Path
        Output file path
    title : str
        Plot title
    top_n : int
        Number of top terms to show
    source_filter : str, optional
        Filter to specific source (e.g., 'GO:BP')
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        logger.warning(f"No data to plot for {title}")
        return

    df = enrichment_df.copy()

    # Filter by source if specified
    if source_filter and 'source' in df.columns:
        df = df[df['source'] == source_filter]

    if len(df) == 0:
        logger.warning(f"No data after filtering for {source_filter}")
        return

    # Sort by p-value and take top N
    df = df.sort_values('p_value').head(top_n)

    # Calculate -log10(p-value)
    df['neg_log10_pval'] = -np.log10(df['p_value'])

    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))

    # Create dotplot
    scatter = ax.scatter(
        df['neg_log10_pval'],
        range(len(df)),
        s=df['intersection_size'] * 10,  # Scale by gene count
        c=df['neg_log10_pval'],
        cmap='Reds',
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )

    # Add term names
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['name'].values)
    ax.set_xlabel('-log10(p-value)')
    ax.set_title(title)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='-log10(p-value)')

    # Add size legend
    legend_sizes = [5, 20, 50]
    legend_dots = [ax.scatter([], [], s=s*10, c='gray', alpha=0.7,
                              edgecolors='black', linewidths=0.5)
                   for s in legend_sizes]
    ax.legend(legend_dots, [str(s) for s in legend_sizes],
              title='Gene count', loc='lower right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved enrichment dotplot to {output_file}")


def create_comparison_heatmap(enrichments_dict, output_file, top_n=15):
    """
    Create heatmap comparing enriched terms across binding categories.

    Parameters:
    -----------
    enrichments_dict : dict
        Dictionary mapping category names to enrichment DataFrames
    output_file : Path
        Output file path
    top_n : int
        Number of top terms per category to include
    """
    # Collect top terms from each category
    all_terms = set()
    for category, df in enrichments_dict.items():
        if df is not None and len(df) > 0:
            top_terms = df.sort_values('p_value').head(top_n)['name'].tolist()
            all_terms.update(top_terms)

    if len(all_terms) == 0:
        logger.warning("No terms to compare")
        return

    # Build matrix
    categories = list(enrichments_dict.keys())
    matrix_data = []

    for term in all_terms:
        row = {'Term': term}
        for category in categories:
            df = enrichments_dict[category]
            if df is not None and len(df) > 0:
                term_row = df[df['name'] == term]
                if len(term_row) > 0:
                    row[category] = -np.log10(term_row['p_value'].values[0])
                else:
                    row[category] = 0
            else:
                row[category] = 0
        matrix_data.append(row)

    matrix_df = pd.DataFrame(matrix_data)
    matrix_df = matrix_df.set_index('Term')

    # Sort by maximum value across categories
    matrix_df['max_val'] = matrix_df.max(axis=1)
    matrix_df = matrix_df.sort_values('max_val', ascending=False).drop('max_val', axis=1)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, max(6, len(matrix_df) * 0.3)))

    sns.heatmap(
        matrix_df,
        cmap='YlOrRd',
        annot=True,
        fmt='.1f',
        linewidths=0.5,
        ax=ax,
        cbar_kws={'label': '-log10(p-value)'}
    )

    ax.set_title('Enrichment Comparison Across Binding Categories')
    ax.set_xlabel('Binding Category')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved comparison heatmap to {output_file}")


def main():
    logger.info("="*70)
    logger.info("IFIT2-IFIT3 eCLIP Analysis: Functional Enrichment")
    logger.info("="*70)

    # Load gene sets
    logger.info("\nLoading bound gene sets...")
    gene_sets = load_bound_genes(RESULTS_DIR)

    if gene_sets is None or len(gene_sets) == 0:
        logger.error("Failed to load gene sets")
        sys.exit(1)

    # Get background
    logger.info("\nGetting background gene set...")
    background = get_background_genes(RESULTS_DIR)

    # Run enrichment for each category
    all_enrichments = {}
    go_bp_enrichments = {}
    kegg_enrichments = {}

    for category, genes in gene_sets.items():
        logger.info(f"\n{'='*50}")
        logger.info(f"Running enrichment for {category} ({len(genes)} genes)")
        logger.info("="*50)

        # Run g:Profiler
        results = run_gprofiler_enrichment(
            genes,
            background=background,
            sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']
        )

        if results is not None and len(results) > 0:
            # Save full results
            results.to_csv(ENRICHMENT_DIR / f'{category}_enrichment.csv', index=False)

            all_enrichments[category] = results

            # Separate by source
            if 'source' in results.columns:
                go_bp = results[results['source'] == 'GO:BP']
                go_mf = results[results['source'] == 'GO:MF']
                go_cc = results[results['source'] == 'GO:CC']
                kegg = results[results['source'] == 'KEGG']
                reac = results[results['source'] == 'REAC']

                go_bp_enrichments[category] = go_bp
                kegg_enrichments[category] = kegg

                # Save separated results
                if len(go_bp) > 0:
                    go_bp.to_csv(ENRICHMENT_DIR / f'{category}_GO_BP.csv', index=False)
                if len(go_mf) > 0:
                    go_mf.to_csv(ENRICHMENT_DIR / f'{category}_GO_MF.csv', index=False)
                if len(go_cc) > 0:
                    go_cc.to_csv(ENRICHMENT_DIR / f'{category}_GO_CC.csv', index=False)
                if len(kegg) > 0:
                    kegg.to_csv(ENRICHMENT_DIR / f'{category}_KEGG.csv', index=False)
                if len(reac) > 0:
                    reac.to_csv(ENRICHMENT_DIR / f'{category}_Reactome.csv', index=False)

                # Log top results
                logger.info(f"\nTop GO Biological Process terms for {category}:")
                for _, row in go_bp.head(5).iterrows():
                    logger.info(f"  {row['name']}: p={row['p_value']:.2e}")

                logger.info(f"\nTop KEGG pathways for {category}:")
                for _, row in kegg.head(5).iterrows():
                    logger.info(f"  {row['name']}: p={row['p_value']:.2e}")
        else:
            all_enrichments[category] = pd.DataFrame()
            go_bp_enrichments[category] = pd.DataFrame()
            kegg_enrichments[category] = pd.DataFrame()

    # Create visualizations
    logger.info("\n" + "="*50)
    logger.info("Generating visualizations")
    logger.info("="*50)

    # Individual dotplots for each category
    for category, results in all_enrichments.items():
        if results is not None and len(results) > 0:
            # GO:BP dotplot
            create_enrichment_dotplot(
                results,
                FIGURES_DIR / f'{category}_GO_BP_dotplot.pdf',
                title=f'{category} - GO Biological Process',
                source_filter='GO:BP'
            )

            # KEGG dotplot
            create_enrichment_dotplot(
                results,
                FIGURES_DIR / f'{category}_KEGG_dotplot.pdf',
                title=f'{category} - KEGG Pathways',
                source_filter='KEGG'
            )

    # Comparison heatmaps
    create_comparison_heatmap(
        go_bp_enrichments,
        FIGURES_DIR / 'GO_BP_comparison_heatmap.pdf'
    )

    create_comparison_heatmap(
        kegg_enrichments,
        FIGURES_DIR / 'KEGG_comparison_heatmap.pdf'
    )

    logger.info("\n" + "="*70)
    logger.info("Functional enrichment analysis complete!")
    logger.info("="*70)
    logger.info("\nResults saved to:")
    logger.info(f"  {ENRICHMENT_DIR}/")
    logger.info(f"  {FIGURES_DIR}/*_dotplot.pdf")
    logger.info(f"  {FIGURES_DIR}/*_comparison_heatmap.pdf")


if __name__ == '__main__':
    main()
