"""
Utility functions for IFIT2-IFIT3 eCLIP Analysis Pipeline

This module contains all utility functions for:
- Checking tool availability
- Downloading SRA data
- Downloading reference genomes
- Building STAR index
- Running FastQC
- UMI extraction and adapter trimming
- STAR alignment
- PCR duplicate removal
"""

import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# eCLIP adapter sequences (Yeo lab protocol)
ECLIP_ADAPTERS = {
    # InvRNA adapters used in eCLIP
    'single_end_3prime': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',  # Illumina TruSeq adapter
    # Additional common adapters
    'illumina_universal': 'AGATCGGAAGAG'
}


def check_tool_available(tool_name):
    """Check if a command-line tool is available."""
    try:
        result = subprocess.run(['which', tool_name], capture_output=True, text=True)
        return result.returncode == 0
    except:
        return False


def download_sra_fastq(srr_accession, output_dir, threads=4):
    """
    Download FASTQ files from SRA using fasterq-dump.

    IMPORTANT: These are single-end 75nt reads, not paired-end!

    Parameters:
    -----------
    srr_accession : str
        SRA run accession (e.g., 'SRR31773513')
    output_dir : Path
        Output directory for FASTQ files
    threads : int
        Number of threads for download
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if files already exist
    existing_files = list(output_dir.glob(f"{srr_accession}*.fastq*"))
    if existing_files:
        logger.info(f"Files already exist for {srr_accession}, skipping download")
        return existing_files

    logger.info(f"Downloading {srr_accession}...")

    # Single-end reads - NO --split-files flag!
    cmd = [
        'fasterq-dump',
        '--threads', str(threads),
        '--outdir', str(output_dir),
        '--progress',
        srr_accession
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully downloaded {srr_accession}")

        # Compress the files
        for fastq_file in output_dir.glob(f"{srr_accession}*.fastq"):
            logger.info(f"Compressing {fastq_file.name}...")
            subprocess.run(['gzip', str(fastq_file)], check=True)

        return list(output_dir.glob(f"{srr_accession}*.fastq.gz"))

    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR downloading {srr_accession}: {e.stderr}")
        return None


def download_reference_files(urls, output_dir):
    """
    Download reference genome files.

    Parameters:
    -----------
    urls : dict
        Dictionary mapping file types to URLs
    output_dir : Path
        Output directory for reference files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded = {}

    for name, url in urls.items():
        filename = url.split('/')[-1]
        output_path = output_dir / filename

        if output_path.exists():
            logger.info(f"{name}: Already exists at {output_path}")
            downloaded[name] = output_path
            continue

        logger.info(f"Downloading {name}...")
        try:
            subprocess.run(['wget', '-O', str(output_path), url], check=True)
            downloaded[name] = output_path
            logger.info(f"Downloaded to {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"ERROR downloading {name}: {e}")

    return downloaded


def build_star_index(genome_fasta, gtf_file, output_dir, threads=8, sjdb_overhang=74):
    """
    Build STAR genome index.

    Parameters:
    -----------
    genome_fasta : Path
        Path to genome FASTA file
    gtf_file : Path
        Path to GTF annotation file
    output_dir : Path
        Output directory for STAR index
    threads : int
        Number of threads
    sjdb_overhang : int
        Read length - 1 (for 75nt eCLIP reads, use 74)
    """
    output_dir = Path(output_dir)
    index_dir = output_dir / "STAR_index_hg19"

    # Check if index already exists
    if (index_dir / "SA").exists():
        logger.info(f"STAR index already exists at {index_dir}")
        return index_dir

    index_dir.mkdir(parents=True, exist_ok=True)

    # Decompress files if needed
    genome_fasta = Path(genome_fasta)
    gtf_file = Path(gtf_file)

    if str(genome_fasta).endswith('.gz'):
        logger.info("Decompressing genome FASTA...")
        subprocess.run(['gunzip', '-k', str(genome_fasta)], check=True)
        genome_fasta = Path(str(genome_fasta)[:-3])

    if str(gtf_file).endswith('.gz'):
        logger.info("Decompressing GTF...")
        subprocess.run(['gunzip', '-k', str(gtf_file)], check=True)
        gtf_file = Path(str(gtf_file)[:-3])

    logger.info(f"Building STAR index (this may take 30+ minutes)...")
    logger.info(f"Using sjdbOverhang={sjdb_overhang} for 75nt reads")

    cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--runThreadN', str(threads),
        '--genomeDir', str(index_dir),
        '--genomeFastaFiles', str(genome_fasta),
        '--sjdbGTFfile', str(gtf_file),
        '--sjdbOverhang', str(sjdb_overhang)
    ]

    try:
        subprocess.run(cmd, check=True)
        logger.info(f"STAR index built successfully at {index_dir}")
        return index_dir
    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR building STAR index: {e}")
        return None


def run_fastqc(fastq_files, output_dir, threads=4):
    """
    Run FastQC on FASTQ files.

    Parameters:
    -----------
    fastq_files : list
        List of FASTQ file paths
    output_dir : Path
        Output directory for FastQC reports
    threads : int
        Number of threads
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running FastQC on {len(fastq_files)} files...")

    cmd = [
        'fastqc',
        '--outdir', str(output_dir),
        '--threads', str(threads),
        '--quiet'
    ] + [str(f) for f in fastq_files]

    try:
        subprocess.run(cmd, check=True)
        logger.info(f"FastQC complete. Reports in {output_dir}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR running FastQC: {e}")
        return False


def parse_cutadapt_output_single_end(log_text):
    """
    Parse cutadapt log output for statistics (single-end version).

    Parameters:
    -----------
    log_text : str
        Standard error output from cutadapt

    Returns:
    --------
    dict : Trimming statistics
    """
    stats = {}

    for line in log_text.split('\n'):
        if 'Total reads processed:' in line:
            stats['reads_input'] = line.split(':')[1].strip().replace(',', '')
        elif 'Reads written (passing filters):' in line:
            stats['reads_written'] = line.split(':')[1].strip().split()[0].replace(',', '')
        elif 'Total basepairs processed:' in line:
            stats['bp_input'] = line.split(':')[1].strip().replace(',', '').split()[0]
        elif 'Total written (filtered):' in line:
            stats['bp_written'] = line.split(':')[1].strip().replace(',', '').split()[0]

    return stats


def extract_umi_and_trim_single_end(fastq, output_dir,
                                     umi_length=10, min_length=18,
                                     quality_cutoff=10, threads=4):
    """
    Extract UMI and trim adapters from SINGLE-END eCLIP reads.

    For single-end eCLIP, the UMI is typically at the 5' end of the read.

    Parameters:
    -----------
    fastq : Path
        Input FASTQ file (single-end)
    output_dir : Path
        Output directory
    umi_length : int
        Length of UMI sequence (default 10)
    min_length : int
        Minimum read length after trimming
    quality_cutoff : int
        Phred quality cutoff for trimming
    threads : int
        Number of threads
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = Path(fastq).stem.replace('.fastq', '').replace('.gz', '')

    # Step 1: Extract UMI using umi_tools (single-end mode)
    logger.info(f"Extracting UMI from {sample_name}...")

    umi_fastq = output_dir / f"{sample_name}_umi.fastq.gz"

    # UMI pattern: first 10 bases (for single-end)
    umi_cmd = [
        'umi_tools', 'extract',
        '--bc-pattern', 'N' * umi_length,  # UMI at 5' end
        '--stdin', str(fastq),
        '--stdout', str(umi_fastq),
        '--log', str(output_dir / f"{sample_name}_umi_extract.log")
    ]

    try:
        subprocess.run(umi_cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"UMI extraction failed, proceeding without UMI: {e}")
        umi_fastq = fastq

    # Step 2: Adapter trimming with cutadapt (single-end mode)
    logger.info(f"Trimming adapters from {sample_name}...")

    trimmed_fastq = output_dir / f"{sample_name}_trimmed.fastq.gz"

    trim_cmd = [
        'cutadapt',
        '-j', str(threads),
        '-a', ECLIP_ADAPTERS['single_end_3prime'],  # 3' adapter
        '-q', str(quality_cutoff),                   # Quality trimming
        '-m', str(min_length),                       # Minimum length
        '--trim-n',                                  # Trim N's
        '-o', str(trimmed_fastq),
        str(umi_fastq)
    ]

    try:
        result = subprocess.run(trim_cmd, check=True, capture_output=True, text=True)

        # Parse cutadapt output for statistics
        stats = parse_cutadapt_output_single_end(result.stderr)
        logger.info(f"Trimming complete: {stats.get('reads_written', 'N/A')} reads passed")

        return {'fastq': trimmed_fastq, 'stats': stats}

    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR trimming: {e.stderr}")
        return None


def parse_star_log(log_file):
    """
    Parse STAR Log.final.out for alignment statistics.

    Parameters:
    -----------
    log_file : Path
        Path to STAR Log.final.out file

    Returns:
    --------
    dict : Alignment statistics
    """
    stats = {}

    if not Path(log_file).exists():
        return stats

    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if '|' in line:
                parts = line.split('|')
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    stats[key] = value

    return stats


def run_star_alignment_single_end(fastq, genome_dir, output_dir,
                                   sample_name, threads=8):
    """
    Run STAR alignment for single-end eCLIP data.

    Parameters optimized for eCLIP:
    - Allow multimapping (up to 10)
    - Soft-clipping enabled
    - Output sorted BAM

    Parameters:
    -----------
    fastq : Path
        Input FASTQ file
    genome_dir : Path
        STAR genome index directory
    output_dir : Path
        Output directory
    sample_name : str
        Sample name for output files
    threads : int
        Number of threads
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / f"{sample_name}_"

    logger.info(f"Aligning {sample_name} with STAR (single-end mode)...")

    cmd = [
        'STAR',
        '--runThreadN', str(threads),
        '--genomeDir', str(genome_dir),
        '--readFilesIn', str(fastq),  # Single FASTQ file only
        '--readFilesCommand', 'zcat',  # For gzipped input
        '--outFileNamePrefix', str(output_prefix),

        # Output format
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outSAMunmapped', 'Within',
        '--outSAMattributes', 'All',

        # Alignment parameters for eCLIP
        '--outFilterMultimapNmax', '10',      # Allow up to 10 multi-mappers
        '--outFilterMismatchNmax', '10',      # Max mismatches
        '--outFilterMismatchNoverLmax', '0.04',
        '--alignIntronMin', '20',
        '--alignIntronMax', '1000000',
        '--alignSJoverhangMin', '8',
        '--alignSJDBoverhangMin', '1',
        '--outFilterType', 'BySJout',

        # Chimeric alignment (for circular RNAs)
        '--chimSegmentMin', '20',

        # Memory and performance
        '--limitBAMsortRAM', '30000000000'
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        bam_file = output_prefix.parent / f"{sample_name}_Aligned.sortedByCoord.out.bam"

        # Index the BAM file
        subprocess.run(['samtools', 'index', str(bam_file)], check=True)

        logger.info(f"Alignment complete: {bam_file}")

        # Parse alignment stats
        stats = parse_star_log(output_prefix.parent / f"{sample_name}_Log.final.out")

        return {'bam': bam_file, 'stats': stats}

    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR in STAR alignment: {e.stderr}")
        return None


def deduplicate_bam(input_bam, output_dir, use_umi=True):
    """
    Remove PCR duplicates using UMI information.

    Parameters:
    -----------
    input_bam : Path
        Input BAM file (coordinate sorted)
    output_dir : Path
        Output directory
    use_umi : bool
        Use UMI for deduplication (requires UMI in read name)
    """
    input_bam = Path(input_bam)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = input_bam.stem.replace('_Aligned.sortedByCoord.out', '')
    dedup_bam = output_dir / f"{sample_name}_dedup.bam"

    logger.info(f"Deduplicating {sample_name}...")

    if use_umi:
        # Use umi_tools dedup
        cmd = [
            'umi_tools', 'dedup',
            '--stdin', str(input_bam),
            '--stdout', str(dedup_bam),
            '--log', str(output_dir / f"{sample_name}_dedup.log"),
            '--output-stats', str(output_dir / f"{sample_name}_dedup_stats")
        ]
    else:
        # Use samtools markdup
        cmd = [
            'samtools', 'markdup',
            '-r',  # Remove duplicates
            '-s',  # Report stats
            str(input_bam),
            str(dedup_bam)
        ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)

        # Index the deduplicated BAM
        subprocess.run(['samtools', 'index', str(dedup_bam)], check=True)

        # Get read counts
        input_count = int(subprocess.run(
            ['samtools', 'view', '-c', str(input_bam)],
            capture_output=True, text=True
        ).stdout.strip())

        output_count = int(subprocess.run(
            ['samtools', 'view', '-c', str(dedup_bam)],
            capture_output=True, text=True
        ).stdout.strip())

        dup_rate = 1 - (output_count / input_count) if input_count > 0 else 0

        logger.info(f"Deduplication complete: {output_count:,} reads ({dup_rate:.1%} duplicates removed)")

        return {
            'bam': dedup_bam,
            'input_reads': input_count,
            'output_reads': output_count,
            'duplicate_rate': dup_rate
        }

    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR in deduplication: {e}")
        return None
