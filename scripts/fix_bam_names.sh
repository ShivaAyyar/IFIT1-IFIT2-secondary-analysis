#!/bin/bash
# Script to rename BAM files to match expected naming convention
# and clean up STAR temporary directories

ALIGNED_DIR="data/aligned"

echo "=== Renaming BAM files in $ALIGNED_DIR ==="
echo ""

for bam in "$ALIGNED_DIR"/*_dedup.bam; do
    if [ -f "$bam" ]; then
        # Extract sample name (everything before _dedup.bam)
        sample_name=$(basename "$bam" _dedup.bam)

        # New filename
        new_name="${ALIGNED_DIR}/${sample_name}_Aligned.sortedByCoord.dedup.bam"

        echo "Renaming: $bam -> $new_name"
        mv "$bam" "$new_name"

        # Rename .bai index file if it exists
        if [ -f "${bam}.bai" ]; then
            mv "${bam}.bai" "${new_name}.bai"
            echo "  Also renamed index: ${bam}.bai -> ${new_name}.bai"
        fi
    fi
done

echo ""
echo "=== Cleaning up STAR temporary directories ==="
echo ""

# Remove STAR temporary directories
for tmpdir in "$ALIGNED_DIR"/*_STARtmp; do
    if [ -d "$tmpdir" ]; then
        echo "Removing: $tmpdir"
        rm -rf "$tmpdir"
    fi
done

echo ""
echo "=== Done! Verifying files ==="
echo ""
echo "BAM files:"
ls -lh "$ALIGNED_DIR"/*_Aligned.sortedByCoord.dedup.bam
echo ""
echo "Remaining directories (should be none):"
ls -d "$ALIGNED_DIR"/*_STARtmp 2>/dev/null || echo "  No STAR temporary directories found (good!)"
