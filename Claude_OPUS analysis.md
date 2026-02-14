# CLIPper and eCLIP integration issues with custom pipelines

Your codebase uses CLIPper with the command `clipper --species hg19 --bam IP_BAM --outfile OUTPUT_BED`, which is **syntactically correct**—but there's a critical format mismatch: CLIPper outputs **8 columns (BED8)**, not the standard BED6 your downstream processing expects. This is likely causing silent failures in your normalization scripts. Below is a complete analysis of integration gaps and compatibility issues.

---

## CLI interface is correct but incomplete

Your command syntax using long-form arguments (`--species`, `--bam`, `--outfile`) matches CLIPper's accepted parameters. Both short (`-s`, `-b`, `-o`) and long forms work:

```bash
# Your codebase (valid)
clipper --species hg19 --bam IP_BAM --outfile OUTPUT_BED

# ENCODE official usage (adds recommended flags)
clipper --species hg19 --bam IP_BAM --outfile OUTPUT_BED --bonferroni --superlocal --threshold-method binomial --save-pickle
```

The ENCODE eCLIP pipeline adds `--bonferroni` for multiple testing correction, `--superlocal` for calculating p-values against local context (±500 bp) rather than whole genes, and `--save-pickle` for intermediate files. Without these flags, your peaks may have different statistical properties than ENCODE-standard calls.

---

## Output format is BED8, not BED6

**This is your most likely integration issue.** CLIPper outputs an 8-column format, but your codebase expects BED6:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | chrom | Chromosome |
| 2 | start | Peak start (0-based) |
| 3 | end | Peak end |
| 4 | name | Format: `geneID_peakNumber_readCount` |
| 5 | score | **Raw p-value** (not browser score) |
| 6 | strand | + or - |
| 7 | peak_center_start | **Extra column** |
| 8 | peak_center_end | **Extra column** |

Example CLIPper output:
```
chr1  133723  133804  ENSG00000233750.3_0_4  0.006532397293615632  +  133761  133765
```

**Critical issues for downstream processing:**
- Column 5 contains raw p-values (e.g., `3.452872201838815e-29`), not BED-compatible integer scores
- Columns 7-8 will break any tool expecting exactly 6 columns
- The name field format (`gene_peak_reads`) may not parse correctly if your code expects simple identifiers

**Fix:** Either truncate to 6 columns with `cut -f1-6` or update your downstream parsers to handle BED8.

---

## BAM file requirements are strict

CLIPper enforces specific BAM requirements that may not be documented in your pipeline:

**Required:**
- **Coordinate-sorted** (not name-sorted): `samtools sort -o sorted.bam input.bam`
- **Indexed** (`.bai` file required): CLIPper calls `samtools index` automatically if missing, but this requires samtools in `$PATH`
- For paired-end data, **only R2 reads**: `samtools view -f 128 -b -o r2.bam merged.bam`

**Prohibited:**
- Do NOT run `samtools rmdup`—the wiki explicitly warns "bad things happen"
- Random/alternate contigs (`chr1_KI270708v1_random`) may cause errors

**Supported species codes:** `hg19`, `GRCh38`, `mm9`, `mm10`, `ce10`, `dm3`

---

## Perl normalization scripts have format-sensitive inputs

The `overlap_peakfi_with_bam.pl` and `compress_peaks.pl` scripts have specific requirements your integration must satisfy.

### overlap_peakfi_with_bam.pl usage

```bash
perl overlap_peakfi_with_bam_PE.pl \
  IP.bam \
  INPUT.bam \
  peaks.bed \
  ip_mapped_readnum.txt \
  input_mapped_readnum.txt \
  output.normed.bed
```

**Input requirements:**
- Peak BED file: **Standard BED6 columns** (you must strip columns 7-8 from CLIPper output)
- BAM files: **R2-only reads** with specific SAM flags (147/16 for minus strand, 163/0 for plus)
- Read count files: Plain text with single integer (from `samtools view -cF 4`)

**Required Perl modules:**
```perl
Statistics::Basic (v1.6611)
Statistics::Distributions (v1.02)
Statistics::R (v0.34)
```

**Perl version matters:** Use Perl 5.10.1 for deterministic output. Hash key ordering changed in Perl 5.18+, causing non-reproducible results.

### Common integration errors

| Error | Cause | Fix |
|-------|-------|-----|
| "R1 strand error 83/99" | BAM contains R1 reads | Filter with `samtools view -f 128` |
| "Use of uninitialized value $strand" | Unexpected SAM flags | Ensure eCLIP preprocessing was run |
| Non-deterministic output | Perl 5.18+ hash ordering | Pin Perl to 5.10.1 |

### Correct pipeline order

```
CLIPper output → cut -f1-6 → overlap_peakfi_with_bam.pl → compress_l2foldenrpeakfi.pl
```

---

## Compatibility issues with modern Python environments

CLIPper has significant compatibility problems with recent Python, numpy, and Cython versions that may affect your installation.

### Python version support

| Version | Status |
|---------|--------|
| Python 3.7-3.8 | ✅ Supported (use `environment3.yml`) |
| Python 3.9+ | ⚠️ Untested, may have issues |
| Python 2.7 | ❌ Legacy only |

### Critical numpy incompatibilities

**numpy ≥1.24:** Removed `np.float` alias entirely
```python
AttributeError: module 'numpy' has no attribute 'float'
```
**Fix:** Pin `numpy<1.24` or `numpy==1.23.0`

**numpy 2.0+:** Binary incompatibility with modules compiled against numpy 1.x
```python
A module that was compiled using NumPy 1.x cannot be run in NumPy 2.0.0
```
**Fix:** Pin `numpy<2`

### Cython compilation failures on Python 3

The C extension uses Python 2 API functions:
```
error: 'PyInt_FromLong' was not declared in this scope
error: 'Py_InitModule' was not declared in this scope
```

This affects the `peaks.so` module and may cause:
```
clipper -h  # Shows wrong help message from peaksmodule.cc fallback
```

### pandas deprecations

Older CLIPper versions use deprecated `DataFrame.sort()`:
```python
AttributeError: 'DataFrame' object has no attribute 'sort'
```
**Fix:** Use `pandas==0.19.*` for legacy installations or update to current CLIPper master.

---

## Recommended installation procedure

Avoid `pip install`—it pulls from a potentially outdated package. Use this approach:

```bash
# Create clean environment
unset PYTHONPATH
conda create -n clipper3 python=3.8 numpy=1.23.0 cython pandas scipy

# Clone and install
git clone https://github.com/YeoLab/clipper.git
cd clipper
conda env update -f environment3.yml
conda activate clipper3
python setup.py install  # NOT pip install

# Verify
clipper -h  # Should show full argument list, not legacy peaksmodule message
```

**Docker alternative** (most reliable):
```bash
docker pull brianyee/clipper:5d865bb
```

---

## Integration checklist for your pipeline

Based on this analysis, verify these items in your codebase:

1. **[ ] Strip CLIPper output to BED6:** Add `cut -f1-6` before downstream processing
2. **[ ] Convert p-values to scores:** Column 5 contains raw p-values; use `-log10(pval)*10` for browser compatibility
3. **[ ] Filter BAM to R2-only:** Add `samtools view -f 128` for paired-end data
4. **[ ] Generate read count files:** Add `samtools view -cF 4 $bam > readcount.txt` before perl scripts
5. **[ ] Pin numpy version:** Ensure `numpy<1.24` in requirements
6. **[ ] Add CLIPper flags:** Consider `--bonferroni --superlocal` for ENCODE-compatible calls
7. **[ ] Check Perl modules:** Verify Statistics::Basic, Statistics::Distributions, Statistics::R are installed
8. **[ ] Validate BAM preprocessing:** Coordinate-sorted, indexed, deduplicated (not with `rmdup`)

---

## Conclusion

The primary integration issue is the **BED8 vs BED6 format mismatch**—CLIPper's extra columns and raw p-values in column 5 will break standard BED parsers and normalization scripts expecting integer scores. Secondary issues include missing preprocessing steps (R2-only filtering, read count generation), Python/numpy version incompatibilities, and undocumented Perl module dependencies. Addressing the format conversion and dependency pinning should resolve most pipeline failures.