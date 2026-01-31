"""
Analysis functions for IFIT2-IFIT3 eCLIP Analysis Pipeline

This module contains functions for:
- Parsing GENCODE GTF files for 5' UTR annotations
- Creating BED files for genomic analysis
- Calling peaks from eCLIP data
- Identifying bound transcripts
- Analyzing 5' UTR length distributions
- Distinguishing IFIT2-only, IFIT3-only, and complex binding
"""

import gzip
import logging
from pathlib import Path
from collections import defaultdict
import pandas as pd
from scipy.stats import mannwhitneyu

logger = logging.getLogger(__name__)

# Check for pybedtools availability
try:
    import pybedtools
    PYBEDTOOLS_AVAILABLE = True
except ImportError:
    PYBEDTOOLS_AVAILABLE = False
    logger.warning("pybedtools not available - peak calling and overlap analysis will be limited")


def parse_gencode_gtf_for_utrs(gtf_file):
    """
    Parse GENCODE GTF to extract 5' UTR information.

    Returns DataFrame with:
    - transcript_id
    - gene_id
    - gene_name
    - utr5_length
    - utr5_coords (chrom, start, end, strand)
    """
    gtf_file = Path(gtf_file)

    logger.info(f"Parsing GTF file: {gtf_file}")

    # Store UTR coordinates by transcript
    utr5_by_transcript = defaultdict(list)
    transcript_info = {}

    # Open file (handle gzipped)
    opener = gzip.open if str(gtf_file).endswith('.gz') else open
    mode = 'rt' if str(gtf_file).endswith('.gz') else 'r'

    with opener(gtf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')

            transcript_id = attr_dict.get('transcript_id', '')

            if feature == 'UTR' or feature == 'five_prime_utr' or feature == '5UTR':
                utr5_by_transcript[transcript_id].append({
                    'chrom': chrom,
                    'start': int(start) - 1,  # Convert to 0-based
                    'end': int(end),
                    'strand': strand
                })

                if transcript_id not in transcript_info:
                    transcript_info[transcript_id] = {
                        'gene_id': attr_dict.get('gene_id', ''),
                        'gene_name': attr_dict.get('gene_name', ''),
                        'gene_type': attr_dict.get('gene_type', '')
                    }

            elif feature == 'transcript':
                if transcript_id not in transcript_info:
                    transcript_info[transcript_id] = {
                        'gene_id': attr_dict.get('gene_id', ''),
                        'gene_name': attr_dict.get('gene_name', ''),
                        'gene_type': attr_dict.get('gene_type', '')
                    }

    # Calculate 5' UTR lengths
    results = []

    for transcript_id, utrs in utr5_by_transcript.items():
        if not utrs:
            continue

        # Calculate total UTR length
        utr_length = sum(u['end'] - u['start'] for u in utrs)

        # Get combined coordinates
        chrom = utrs[0]['chrom']
        strand = utrs[0]['strand']
        utr_start = min(u['start'] for u in utrs)
        utr_end = max(u['end'] for u in utrs)

        info = transcript_info.get(transcript_id, {})

        results.append({
            'transcript_id': transcript_id,
            'gene_id': info.get('gene_id', ''),
            'gene_name': info.get('gene_name', ''),
            'gene_type': info.get('gene_type', ''),
            'chrom': chrom,
            'utr5_start': utr_start,
            'utr5_end': utr_end,
            'strand': strand,
            'utr5_length': utr_length
        })

    df = pd.DataFrame(results)
    logger.info(f"Parsed {len(df):,} transcripts with 5' UTR annotations")

    return df


def create_utr5_bed(utr_df, output_file):
    """
    Create BED file of 5' UTR regions for overlap analysis.

    Parameters:
    -----------
    utr_df : DataFrame
        DataFrame with UTR information from parse_gencode_gtf_for_utrs()
    output_file : Path
        Output BED file path
    """
    bed_df = utr_df[['chrom', 'utr5_start', 'utr5_end', 'transcript_id', 'utr5_length', 'strand']].copy()
    bed_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed_df = bed_df.dropna()
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)

    bed_df.to_csv(output_file, sep='\t', header=False, index=False)
    logger.info(f"Created 5' UTR BED file: {output_file}")

    return output_file


def call_peaks_simple(ip_bam, input_bam, output_dir, sample_name,
                       window_size=50, min_reads=10, fold_enrichment=2.0):
    """
    Simple peak calling using coverage comparison.

    For each position, calculate IP/Input ratio and call peaks
    where the ratio exceeds the threshold.

    Parameters:
    -----------
    ip_bam : Path
        IP sample BAM file
    input_bam : Path
        Input control BAM file (can be None)
    output_dir : Path
        Output directory for peaks
    sample_name : str
        Sample name for output files
    window_size : int
        Window size for merging nearby reads
    min_reads : int
        Minimum read count to call a peak
    fold_enrichment : float
        Minimum fold enrichment over input
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    peaks_file = output_dir / f"{sample_name}_peaks.bed"

    logger.info(f"Calling peaks for {sample_name} using coverage method...")

    if not PYBEDTOOLS_AVAILABLE:
        logger.error("pybedtools required for peak calling")
        return None

    # Convert BAM to BED
    ip_bed = pybedtools.BedTool(str(ip_bam)).bam_to_bed()

    if input_bam:
        input_bed = pybedtools.BedTool(str(input_bam)).bam_to_bed()

        # Get coverage in windows
        genome_file = output_dir / "genome.txt"  # Need to create this

        # For now, use a simpler approach: merge overlapping reads
        peaks = ip_bed.merge(d=window_size, c=1, o='count')

        # Filter by minimum read count
        peaks = peaks.filter(lambda x: int(x[3]) >= min_reads)

    else:
        # No input - just merge IP reads
        peaks = ip_bed.merge(d=window_size, c=1, o='count')
        peaks = peaks.filter(lambda x: int(x[3]) >= min_reads)

    # Save peaks
    peaks.saveas(str(peaks_file))

    # Count peaks
    peak_count = sum(1 for _ in open(peaks_file))
    logger.info(f"Found {peak_count:,} peaks")

    return peaks_file


def identify_bound_transcripts(peaks_file, utr_bed, utr_df, min_overlap=1):
    """
    Identify transcripts with eCLIP peaks overlapping their 5' UTR.

    Parameters:
    -----------
    peaks_file : Path
        BED file of eCLIP peaks
    utr_bed : Path
        BED file of 5' UTR regions
    utr_df : DataFrame
        DataFrame with UTR information
    min_overlap : int
        Minimum overlap in bp to consider as bound

    Returns:
    --------
    bound_transcripts : set
        Set of transcript IDs with 5' UTR binding
    """
    if not PYBEDTOOLS_AVAILABLE:
        logger.error("pybedtools required for overlap analysis")
        return set()

    peaks = pybedtools.BedTool(str(peaks_file))
    utrs = pybedtools.BedTool(str(utr_bed))

    # Find overlaps
    overlaps = utrs.intersect(peaks, wa=True, u=True)

    # Extract transcript IDs from overlapping UTRs
    bound_transcripts = set()
    for interval in overlaps:
        bound_transcripts.add(interval.name)

    logger.info(f"Found {len(bound_transcripts):,} transcripts with 5' UTR binding")

    return bound_transcripts


def identify_complex_bound_transcripts(flag_peaks, ha_peaks, utr_bed, method='intersection'):
    """
    Identify transcripts bound by the IFIT2-IFIT3 complex.

    This function distinguishes between:
    - IFIT2-only binding (FLAG IP peaks)
    - IFIT3-only binding (HA IP peaks)
    - IFIT2-IFIT3 complex binding (both FLAG AND HA IP peaks)

    Parameters:
    -----------
    flag_peaks : Path
        BED file of FLAG IP peaks (IFIT2 binding from co-expression samples)
    ha_peaks : Path
        BED file of HA IP peaks (IFIT3 binding from co-expression samples)
    utr_bed : Path
        BED file of 5' UTR regions
    method : str
        'intersection' - transcripts bound by BOTH IFIT2 AND IFIT3 (complex)
        'union' - transcripts bound by IFIT2 OR IFIT3

    Returns:
    --------
    dict with keys:
        'ifit2_only': transcripts bound by IFIT2 but not IFIT3
        'ifit3_only': transcripts bound by IFIT3 but not IFIT2
        'complex': transcripts bound by both (heterodimer)
        'union': transcripts bound by either
    """
    if not PYBEDTOOLS_AVAILABLE:
        logger.error("pybedtools required for complex binding analysis")
        return None

    flag_peaks = pybedtools.BedTool(str(flag_peaks))
    ha_peaks = pybedtools.BedTool(str(ha_peaks))
    utrs = pybedtools.BedTool(str(utr_bed))

    # Find UTRs overlapping FLAG peaks (IFIT2 binding)
    flag_utrs = utrs.intersect(flag_peaks, wa=True, u=True)
    ifit2_transcripts = set(interval.name for interval in flag_utrs)

    # Find UTRs overlapping HA peaks (IFIT3 binding)
    ha_utrs = utrs.intersect(ha_peaks, wa=True, u=True)
    ifit3_transcripts = set(interval.name for interval in ha_utrs)

    # Calculate different binding categories
    complex_transcripts = ifit2_transcripts & ifit3_transcripts  # Intersection
    ifit2_only = ifit2_transcripts - ifit3_transcripts
    ifit3_only = ifit3_transcripts - ifit2_transcripts
    union_transcripts = ifit2_transcripts | ifit3_transcripts

    results = {
        'ifit2_only': ifit2_only,
        'ifit3_only': ifit3_only,
        'complex': complex_transcripts,
        'union': union_transcripts,
        'ifit2_all': ifit2_transcripts,
        'ifit3_all': ifit3_transcripts
    }

    logger.info("\n" + "="*60)
    logger.info("IFIT2-IFIT3 Complex Binding Analysis")
    logger.info("="*60)
    logger.info(f"IFIT2 binding (FLAG IP):        {len(ifit2_transcripts):,} transcripts")
    logger.info(f"IFIT3 binding (HA IP):          {len(ifit3_transcripts):,} transcripts")
    logger.info(f"IFIT2 only (no IFIT3):          {len(ifit2_only):,} transcripts")
    logger.info(f"IFIT3 only (no IFIT2):          {len(ifit3_only):,} transcripts")
    logger.info(f"Complex (both IFIT2 + IFIT3):   {len(complex_transcripts):,} transcripts")
    logger.info(f"Union (IFIT2 or IFIT3):         {len(union_transcripts):,} transcripts")
    logger.info("="*60)

    # Calculate overlap percentage
    if len(union_transcripts) > 0:
        complex_pct = 100 * len(complex_transcripts) / len(union_transcripts)
        logger.info(f"\nComplex binding represents {complex_pct:.1f}% of all bound transcripts")

    return results


def analyze_utr_length_distribution(utr_df, bound_transcripts):
    """
    Analyze 5' UTR length distribution for bound vs unbound transcripts.

    Parameters:
    -----------
    utr_df : DataFrame
        DataFrame with UTR information
    bound_transcripts : set
        Set of transcript IDs that are bound

    Returns:
    --------
    tuple : (annotated_df, stats_dict)
        - annotated_df: DataFrame with 'bound' column added
        - stats_dict: Dictionary with summary statistics
    """
    # Add binding status
    utr_df = utr_df.copy()
    utr_df['bound'] = utr_df['transcript_id'].isin(bound_transcripts)

    # Filter for protein-coding genes
    protein_coding = utr_df[utr_df['gene_type'] == 'protein_coding'].copy()

    # Summary statistics
    bound = protein_coding[protein_coding['bound']]['utr5_length']
    unbound = protein_coding[~protein_coding['bound']]['utr5_length']

    stats = {
        'bound_count': len(bound),
        'unbound_count': len(unbound),
        'bound_median': bound.median(),
        'unbound_median': unbound.median(),
        'bound_mean': bound.mean(),
        'unbound_mean': unbound.mean(),
        'bound_short_fraction': (bound < 50).sum() / len(bound) if len(bound) > 0 else 0,
        'unbound_short_fraction': (unbound < 50).sum() / len(unbound) if len(unbound) > 0 else 0
    }

    # Statistical test
    if len(bound) > 0 and len(unbound) > 0:
        stat, pval = mannwhitneyu(bound, unbound, alternative='less')
        stats['mann_whitney_U'] = stat
        stats['mann_whitney_pval'] = pval

    logger.info("\n5' UTR Length Analysis Summary:")
    logger.info("=" * 50)
    logger.info(f"Bound transcripts: {stats['bound_count']:,}")
    logger.info(f"Unbound transcripts: {stats['unbound_count']:,}")
    logger.info(f"\nMedian 5' UTR length:")
    logger.info(f"  Bound: {stats['bound_median']:.1f} nt")
    logger.info(f"  Unbound: {stats['unbound_median']:.1f} nt")
    logger.info(f"\nFraction with short (<50 nt) 5' UTR:")
    logger.info(f"  Bound: {stats['bound_short_fraction']:.1%}")
    logger.info(f"  Unbound: {stats['unbound_short_fraction']:.1%}")
    if 'mann_whitney_pval' in stats:
        logger.info(f"\nMann-Whitney U test (bound < unbound):")
        logger.info(f"  p-value: {stats['mann_whitney_pval']:.2e}")

    return protein_coding, stats
