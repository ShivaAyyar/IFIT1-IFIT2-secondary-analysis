# CLIPper Installation and Configuration Guide

## The Problem: Python 3.10 Compatibility

CLIPper's C extensions were written for Python 2 API and have compilation issues with Python 3.9+. When the C extension fails to compile, CLIPper falls back to showing help from the legacy `peaksmodule.cc` instead of the actual CLIPper tool.

**Symptoms:**
```bash
$ clipper --help
Usage is � [options]
file_in is a file with a list of the length of aligned reads
Options
 -L <int>   Effective Gene Length
```

This is the WRONG tool. The correct CLIPper should show:
```bash
Usage: clipper [options]

Options:
  --bam, -b      BAM file for peak calling
  --species, -s  Species designation (hg19 or mm9)
  --outfile, -o  Output BED file
```

---

## Solution Options

We provide **three solutions** in priority order:

### Option 1: Singularity Container (Recommended for HPC) ⭐

**Pros:**
- Guaranteed to work - uses YeoLab's exact environment
- No compilation issues
- Most reproducible

**Requirements:**
- Singularity installed on HPC (check with `module avail singularity`)

**Setup:**

```bash
# 1. Check if Singularity is available
module avail singularity
module load singularity  # Or whatever command loads it on your HPC

# 2. Pull the official CLIPper container (from eCLIP CWL pipeline)
# The correct image is brianyee/clipper:5d865bb (specified in eCLIP CWL files)
singularity pull docker://brianyee/clipper:5d865bb

# This creates: clipper_5d865bb.sif (approximately 1-2 GB)

# 3. Verify it works
singularity exec clipper_5d865bb.sif clipper --help
# Should show correct CLIPper help with --bam, --species options

# 4. Configure the pipeline
# Edit scripts/config.py:
CLIPPER_MODE = 'singularity'
CLIPPER_CONTAINER = '/path/to/clipper_5d865bb.sif'  # Use absolute path
```

**Usage:**
```bash
# Run peak calling - it will automatically use the container
python scripts/04_call_peaks.py
```

---

### Option 2: Python 3.8 Environment (If no Singularity)

**Pros:**
- Local installation, no containers needed
- CLIPper officially supports Python 3.7-3.8

**Requirements:**
- Conda/Mamba

**Setup:**

```bash
# 1. Remove old environment (if it exists)
conda deactivate
conda env remove -n eclip -y

# 2. Create Python 3.8 environment
conda env create -f environment_py38.yml

# 3. Activate and verify
conda activate eclip_py38
clipper --help  # Should show correct help now

# 4. If it still shows wrong help, reinstall CLIPper:
pip uninstall clipper -y
pip install --no-build-isolation git+https://github.com/YeoLab/clipper.git@5d865bb17b2bc6787b4c382bc857119ae917ad59

# 5. Verify again
clipper --help

# 6. Configure the pipeline (scripts/config.py)
CLIPPER_MODE = 'local'
CLIPPER_CONTAINER = None
```

**Usage:**
```bash
conda activate eclip_py38
python scripts/04_call_peaks.py
```

---

### Option 3: Docker (For Local Machines)

**Pros:**
- Guaranteed compatibility like Singularity
- Works on local machines

**Requirements:**
- Docker installed locally

**Setup:**

```bash
# 1. Pull Docker image (official CLIPper container from eCLIP CWL)
docker pull brianyee/clipper:5d865bb

# 2. Verify
docker run --rm brianyee/clipper:5d865bb clipper --help
# Should show correct help with --bam, --species options

# 3. Configure the pipeline (scripts/config.py)
CLIPPER_MODE = 'docker'
CLIPPER_CONTAINER = 'brianyee/clipper:5d865bb'
```

**Usage:**
```bash
python scripts/04_call_peaks.py
```

---

## Configuration Reference

Edit `scripts/config.py` to set the CLIPper execution mode:

```python
# CLIPper execution mode
# Options: 'local', 'singularity', 'docker'
CLIPPER_MODE = 'local'  # Change this

# Container image path (only for singularity/docker modes)
# For Singularity: '/path/to/eclip_latest.sif'
# For Docker: 'yeolab/eclip:latest'
CLIPPER_CONTAINER = None  # Set this if using containers
```

---

## Verification

After setup, verify CLIPper works correctly:

### Test 1: Help Output
```bash
# For local installation
clipper --help | head -5

# For Singularity
singularity exec /path/to/eclip_latest.sif clipper --help | head -5

# For Docker
docker run --rm yeolab/eclip:latest clipper --help | head -5
```

**Expected output:**
```
Usage: clipper [options]

Options:
  -h, --help            show this help message and exit
  -b BAM, --bam=BAM     BAM file containing aligned CLIP reads
```

**Wrong output (means it's not working):**
```
Usage is � [options]
file_in is a file with a list of the length of aligned reads
```

### Test 2: Argument Parsing
```bash
# Should complain about missing --bam, NOT show read length help
clipper --species hg19
```

**Expected:**
```
error: --bam is required
```

**Wrong (means it's broken):**
```
Usage is � [options]
file_in is a file with a list of the length of aligned reads
```

---

## Troubleshooting

### Problem: Singularity not found

```bash
$ singularity --version
bash: singularity: command not found
```

**Solutions:**
1. Check if it's a module: `module avail singularity` then `module load singularity`
2. Ask your HPC admin if Singularity is available
3. Use Option 2 (Python 3.8) instead

---

### Problem: Python 3.8 still shows wrong help

This means the C extension still failed to compile. Try:

```bash
# 1. Complete cleanup
conda activate eclip_py38
pip uninstall clipper -y
conda remove clipper -y --force 2>/dev/null || true

# 2. Verify it's gone
python -c "import clipper" 2>&1 | grep "No module"  # Should error

# 3. Reinstall with verbose output to see compilation errors
pip install -v --no-build-isolation git+https://github.com/YeoLab/clipper.git@5d865bb17b2bc6787b4c382bc857119ae917ad59

# 4. Check for compilation errors in output
# Look for lines like: "error: 'PyInt_FromLong' was not declared"

# 5. If still broken, use Singularity (Option 1) instead
```

---

### Problem: Container can't access files

```bash
ERROR: BAM file not found: /path/to/file.bam
```

**For Singularity:**
```bash
# The script automatically binds necessary directories
# If issues persist, manually specify bind paths:

# Edit scripts/analysis.py, in call_peaks_clipper(), add to bind_mounts:
bind_mounts = [
    f"{ip_bam.parent}:{ip_bam.parent}",
    f"{output_dir}:{output_dir}",
    "/your/data/dir:/your/data/dir"  # Add your data directory
]
```

**For Docker:**
Similar issue - ensure volumes are mounted. The script handles this automatically.

---

### Problem: Peak calling hangs/runs forever

CLIPper can be slow on large BAM files. Check:

```bash
# Monitor CLIPper process
top -u $USER

# Check if it's actually running
ps aux | grep clipper

# Check log file
tail -f logs/peak_calling.log
```

If it's genuinely stuck:
1. Kill the process
2. Check BAM file isn't corrupted: `samtools quickcheck file.bam`
3. Try with smaller test BAM file first

---

## Performance Comparison

| Mode | Setup Time | Execution Speed | Reliability |
|------|-----------|----------------|-------------|
| Singularity | 5 min (one-time) | ~same as local | ⭐⭐⭐⭐⭐ |
| Python 3.8 | 10 min (one-time) | Baseline | ⭐⭐⭐ |
| Docker | 5 min (one-time) | ~same as local | ⭐⭐⭐⭐⭐ |
| Python 3.10 (local) | Already done | N/A | ❌ (broken) |

---

## Summary of Changes

The following files were modified to support multiple CLIPper execution modes:

1. **scripts/config.py**: Added `CLIPPER_MODE` and `CLIPPER_CONTAINER` configuration
2. **scripts/analysis.py**: Updated `call_peaks_clipper()` to support local/singularity/docker modes
3. **scripts/04_call_peaks.py**: Passes mode and container configuration to peak calling
4. **environment_py38.yml**: New Python 3.8 environment for Option 2

**No user action required** beyond:
1. Choosing which option to use
2. Setting `CLIPPER_MODE` in config.py
3. Setting `CLIPPER_CONTAINER` path if using containers

---

## Recommended Workflow on HPC

```bash
# Step 1: Check if Singularity is available
module avail singularity

# Step 2A: If yes, use Singularity (best option)
module load singularity
singularity pull docker://yeolab/eclip:latest
# Edit config.py: CLIPPER_MODE='singularity', CLIPPER_CONTAINER='/path/to/eclip_latest.sif'

# Step 2B: If no, use Python 3.8
conda env create -f environment_py38.yml
conda activate eclip_py38
# Edit config.py: CLIPPER_MODE='local', CLIPPER_CONTAINER=None

# Step 3: Fix BAM names (if needed)
chmod +x scripts/fix_bam_names.sh
./scripts/fix_bam_names.sh

# Step 4: Run peak calling
python scripts/04_call_peaks.py

# Step 5: Verify peaks were called
ls -lh data/peaks/*_normalized_peaks.bed
```

---

## Additional Resources

- [YeoLab CLIPper GitHub](https://github.com/YeoLab/clipper)
- [YeoLab eCLIP Pipeline](https://github.com/YeoLab/eCLIP)
- [Singularity Documentation](https://sylabs.io/docs/)
- [ENCODE eCLIP Standards](https://www.encodeproject.org/eclip/)

---

## Getting Help

If none of these options work:

1. **Check the diagnostic output:**
   ```bash
   python -c "import clipper.src.main; print(clipper.src.main.__file__)"
   python -c "import inspect, clipper.src.main; print(inspect.getsource(clipper.src.main.call_main)[:500])"
   ```

2. **Open an issue** with:
   - Python version: `python --version`
   - CLIPper installation method used
   - Full error message
   - Output of diagnostic commands above

3. **Contact:**
   - This repository's issue tracker
   - YeoLab CLIPper issues: https://github.com/YeoLab/clipper/issues
