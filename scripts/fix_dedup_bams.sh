#!/bin/bash
# Script to re-run deduplication for BAM files that are empty (0 bytes)
# This uses the existing STAR output BAM files (*_Aligned.sortedByCoord.out.bam)
# and re-runs umi_tools dedup to create proper deduplicated BAM files

ALIGNED_DIR="data/aligned"

echo "=============================================="
echo "Fixing Empty Deduplicated BAM Files"
echo "=============================================="
echo ""

# Activate conda environment if not already active
if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" != "eclip" ]]; then
    echo "Activating eclip environment..."
    source ~/.bashrc
    conda activate eclip
fi

# Count files to process
total_empty=0
for dedup_bam in "$ALIGNED_DIR"/*_Aligned.sortedByCoord.dedup.bam; do
    if [ -f "$dedup_bam" ] && [ ! -s "$dedup_bam" ]; then
        ((total_empty++))
    fi
done

echo "Found $total_empty empty deduplicated BAM files to fix"
echo ""

# Process each empty dedup BAM
fixed=0
failed=0

for dedup_bam in "$ALIGNED_DIR"/*_Aligned.sortedByCoord.dedup.bam; do
    if [ ! -f "$dedup_bam" ]; then
        continue
    fi

    # Check if file is empty (0 bytes)
    if [ -s "$dedup_bam" ]; then
        echo "SKIP: $dedup_bam (already has data)"
        continue
    fi

    # Get the source BAM (STAR output)
    sample_name=$(basename "$dedup_bam" _Aligned.sortedByCoord.dedup.bam)
    source_bam="$ALIGNED_DIR/${sample_name}_Aligned.sortedByCoord.out.bam"

    if [ ! -f "$source_bam" ]; then
        echo "ERROR: Source BAM not found: $source_bam"
        ((failed++))
        continue
    fi

    echo "Processing: $sample_name"
    echo "  Source: $source_bam"
    echo "  Output: $dedup_bam"

    # Index source BAM if needed
    if [ ! -f "${source_bam}.bai" ]; then
        echo "  Indexing source BAM..."
        samtools index "$source_bam"
    fi

    # Run umi_tools dedup
    # Note: UMI should be in read name (from umi_tools extract during preprocessing)
    echo "  Running umi_tools dedup..."

    umi_tools dedup \
        --stdin="$source_bam" \
        --stdout="$dedup_bam" \
        --log="${ALIGNED_DIR}/${sample_name}_dedup.log" \
        --extract-umi-method=read_id \
        --method=unique \
        2>&1 | tee -a "${ALIGNED_DIR}/${sample_name}_dedup_stdout.log"

    # Check if output has data
    if [ -s "$dedup_bam" ]; then
        # Index the output
        echo "  Indexing output BAM..."
        samtools index "$dedup_bam"

        # Get read count
        read_count=$(samtools view -c "$dedup_bam")
        echo "  SUCCESS: $read_count reads"
        ((fixed++))
    else
        echo "  ERROR: Deduplication produced empty file"
        ((failed++))
    fi

    echo ""
done

echo "=============================================="
echo "Summary"
echo "=============================================="
echo "Fixed: $fixed"
echo "Failed: $failed"
echo "=============================================="

if [ $failed -gt 0 ]; then
    echo ""
    echo "WARNING: Some files failed. Check the log files for details:"
    echo "  ${ALIGNED_DIR}/*_dedup.log"
    echo "  ${ALIGNED_DIR}/*_dedup_stdout.log"
fi
