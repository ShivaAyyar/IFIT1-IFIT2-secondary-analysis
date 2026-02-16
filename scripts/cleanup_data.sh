#!/bin/bash
#===============================================================================
# cleanup_data.sh - Remove unnecessary intermediate files after analysis
#
# This script removes raw/intermediate data that is no longer needed after
# the eCLIP analysis pipeline has completed. Results and essential files are kept.
#
# Usage:
#   chmod +x scripts/cleanup_data.sh
#   ./scripts/cleanup_data.sh [--dry-run]
#
# Options:
#   --dry-run    Show what would be deleted without actually deleting
#===============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check for dry-run mode
DRY_RUN=false
if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN=true
    echo -e "${YELLOW}DRY RUN MODE - No files will be deleted${NC}"
    echo ""
fi

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_DIR/data"

echo "=========================================="
echo "eCLIP Data Cleanup Script"
echo "=========================================="
echo "Project directory: $PROJECT_DIR"
echo "Data directory: $DATA_DIR"
echo ""

# Initialize counters
TOTAL_SIZE=0
FILE_COUNT=0

# Function to calculate size and optionally delete
cleanup_files() {
    local pattern="$1"
    local description="$2"

    # Find files matching pattern
    files=$(find "$DATA_DIR" -type f -name "$pattern" 2>/dev/null || true)

    if [[ -n "$files" ]]; then
        size=$(echo "$files" | xargs du -ch 2>/dev/null | tail -1 | cut -f1)
        count=$(echo "$files" | wc -l)

        echo -e "${YELLOW}[$description]${NC}"
        echo "  Pattern: $pattern"
        echo "  Files: $count"
        echo "  Size: $size"

        if [[ "$DRY_RUN" == false ]]; then
            echo "$files" | xargs rm -f
            echo -e "  ${GREEN}DELETED${NC}"
        else
            echo -e "  ${YELLOW}Would delete${NC}"
        fi
        echo ""

        FILE_COUNT=$((FILE_COUNT + count))
    fi
}

# Function to cleanup directories
cleanup_dir() {
    local dir="$1"
    local description="$2"

    if [[ -d "$dir" ]]; then
        size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        count=$(find "$dir" -type f 2>/dev/null | wc -l)

        echo -e "${YELLOW}[$description]${NC}"
        echo "  Directory: $dir"
        echo "  Files: $count"
        echo "  Size: $size"

        if [[ "$DRY_RUN" == false ]]; then
            rm -rf "$dir"
            echo -e "  ${GREEN}DELETED${NC}"
        else
            echo -e "  ${YELLOW}Would delete${NC}"
        fi
        echo ""

        FILE_COUNT=$((FILE_COUNT + count))
    fi
}

echo "=========================================="
echo "FILES TO BE REMOVED"
echo "=========================================="
echo ""

#-------------------------------------------------------------------------------
# 1. RAW FASTQ FILES (can be re-downloaded from SRA)
#-------------------------------------------------------------------------------
echo -e "${RED}=== RAW DATA (Re-downloadable from SRA) ===${NC}"
cleanup_dir "$DATA_DIR/raw_fastq" "Raw FASTQ files"

#-------------------------------------------------------------------------------
# 2. TRIMMED/INTERMEDIATE FASTQ FILES
#-------------------------------------------------------------------------------
echo -e "${RED}=== INTERMEDIATE FASTQ FILES ===${NC}"
cleanup_dir "$DATA_DIR/trimmed" "Trimmed FASTQ files"

#-------------------------------------------------------------------------------
# 3. PRE-DEDUP BAM FILES (keep only dedup BAMs)
#-------------------------------------------------------------------------------
echo -e "${RED}=== PRE-DEDUPLICATION BAM FILES ===${NC}"
cleanup_files "*_Aligned.sortedByCoord.out.bam" "STAR output BAMs (before dedup)"
cleanup_files "*_Aligned.sortedByCoord.out.bam.bai" "STAR output BAM indexes"

#-------------------------------------------------------------------------------
# 4. STAR LOG FILES (optional - keep for QC)
#-------------------------------------------------------------------------------
echo -e "${RED}=== STAR INTERMEDIATE FILES ===${NC}"
cleanup_files "*_Chimeric.out.junction" "Chimeric junction files"
cleanup_files "*_SJ.out.tab" "Splice junction files"
cleanup_files "*_Log.out" "STAR log files"
cleanup_files "*_Log.progress.out" "STAR progress logs"

#-------------------------------------------------------------------------------
# 5. UMI DEDUP STATS (large TSV files)
#-------------------------------------------------------------------------------
echo -e "${RED}=== UMI DEDUP STATISTICS ===${NC}"
cleanup_files "*_dedup_stats_per_umi.tsv" "Per-UMI dedup stats"
cleanup_files "*_dedup_stats_per_umi_per_position.tsv" "Per-position dedup stats"
cleanup_files "*_dedup_stats_edit_distance.tsv" "Edit distance stats"

#-------------------------------------------------------------------------------
# 6. PEAK CALLING INTERMEDIATE FILES
#-------------------------------------------------------------------------------
echo -e "${RED}=== PEAK CALLING INTERMEDIATES ===${NC}"
cleanup_files "*_clipper_peaks_raw.bed" "Raw CLIPper peaks (before filtering)"
cleanup_files "*_overlap.bed" "Overlap intermediate files"
cleanup_files "*_overlap.bed.full" "Full overlap files"
cleanup_files "*_ip_read_count.txt" "IP read count files"
cleanup_files "*_input_read_count.txt" "Input read count files"

#-------------------------------------------------------------------------------
# 7. REFERENCE FILES (can be re-downloaded)
#-------------------------------------------------------------------------------
echo -e "${RED}=== REFERENCE FILES (Re-downloadable) ===${NC}"
# Keep: 5utr_regions.bed (generated), hg19.chrom.sizes (small), STAR_index (needed)
# Remove: Large genome files that can be re-downloaded
if [[ -f "$DATA_DIR/reference/hg19.fa" ]]; then
    size=$(du -h "$DATA_DIR/reference/hg19.fa" | cut -f1)
    echo -e "${YELLOW}[Genome FASTA]${NC}"
    echo "  File: $DATA_DIR/reference/hg19.fa"
    echo "  Size: $size"
    if [[ "$DRY_RUN" == false ]]; then
        rm -f "$DATA_DIR/reference/hg19.fa"
        echo -e "  ${GREEN}DELETED${NC}"
    else
        echo -e "  ${YELLOW}Would delete${NC}"
    fi
    echo ""
fi

if [[ -f "$DATA_DIR/reference/hg19.fa.gz" ]]; then
    size=$(du -h "$DATA_DIR/reference/hg19.fa.gz" | cut -f1)
    echo -e "${YELLOW}[Compressed Genome FASTA]${NC}"
    echo "  File: $DATA_DIR/reference/hg19.fa.gz"
    echo "  Size: $size"
    if [[ "$DRY_RUN" == false ]]; then
        rm -f "$DATA_DIR/reference/hg19.fa.gz"
        echo -e "  ${GREEN}DELETED${NC}"
    else
        echo -e "  ${YELLOW}Would delete${NC}"
    fi
    echo ""
fi

# GTF can be re-downloaded but keep the uncompressed version for analysis
if [[ -f "$DATA_DIR/reference/gencode.v19.annotation.gtf.gz" ]]; then
    size=$(du -h "$DATA_DIR/reference/gencode.v19.annotation.gtf.gz" | cut -f1)
    echo -e "${YELLOW}[Compressed GTF (keeping uncompressed)]${NC}"
    echo "  File: $DATA_DIR/reference/gencode.v19.annotation.gtf.gz"
    echo "  Size: $size"
    if [[ "$DRY_RUN" == false ]]; then
        rm -f "$DATA_DIR/reference/gencode.v19.annotation.gtf.gz"
        echo -e "  ${GREEN}DELETED${NC}"
    else
        echo -e "  ${YELLOW}Would delete${NC}"
    fi
    echo ""
fi

#-------------------------------------------------------------------------------
# 8. QC FILES (optional - useful for troubleshooting)
#-------------------------------------------------------------------------------
echo -e "${RED}=== QC FILES (Optional) ===${NC}"
cleanup_dir "$DATA_DIR/qc" "FastQC reports"

#-------------------------------------------------------------------------------
# 9. PYTHON CACHE
#-------------------------------------------------------------------------------
echo -e "${RED}=== PYTHON CACHE ===${NC}"
cleanup_dir "$PROJECT_DIR/scripts/__pycache__" "Python bytecode cache"

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
echo "=========================================="
echo "CLEANUP SUMMARY"
echo "=========================================="

if [[ "$DRY_RUN" == true ]]; then
    echo -e "${YELLOW}DRY RUN - No files were deleted${NC}"
    echo "Run without --dry-run to actually delete files"
else
    echo -e "${GREEN}Cleanup complete!${NC}"
fi

echo ""
echo "=========================================="
echo "FILES KEPT (Essential for analysis)"
echo "=========================================="
echo ""
echo "data/aligned/*_dedup.bam           - Deduplicated BAM files"
echo "data/aligned/*_dedup.bam.bai       - BAM index files"
echo "data/aligned/*_Log.final.out       - STAR alignment summaries"
echo "data/aligned/*_dedup.log           - Deduplication logs"
echo "data/peaks/*_normalized_peaks.bed  - Final normalized peaks"
echo "data/peaks/*_clipper_peaks.bed     - Filtered CLIPper peaks"
echo "data/reference/5utr_regions.bed    - Generated 5' UTR annotations"
echo "data/reference/hg19.chrom.sizes    - Chromosome sizes"
echo "data/reference/STAR_index_hg19/    - STAR genome index"
echo "data/reference/gencode.v19.annotation.gtf - GTF annotation"
echo "data/results/                      - All analysis results"
echo ""

# Calculate space after cleanup
if [[ "$DRY_RUN" == false ]]; then
    remaining=$(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)
    echo "Remaining data directory size: $remaining"
fi
