# Gap Filling Workflow Documentation

## Overview

This is a multi-round, iterative gap filling pipeline for genome assemblies using long-read sequencing data (HiFi and ONT reads).

**Core Principle**: The N sequences in an assembly are **placeholders only** - their length does NOT represent the actual biological gap size.

---

## File Structure

```
gapfilling/
├── iterative_gapfiller.py   # Main entry point - orchestrates the workflow
├── gap_filler.py            # Core gap filling engine (OPTIMIZED)
├── gap_validation.py        # Validation module
├── assembly_indexer.py      # Fast indexed FASTA access
├── temp_file_manager.py     # Temporary file management
├── find_gaps.py             # Standalone gap scanning tool
├── read_aligner.py          # Standalone alignment tool
└── WORKFLOW.md              # This documentation
```

---

## Key Optimizations

### 1. Supplementary Alignment Detection

**Problem**: When a read truly crosses a gap, minimap2 splits it into primary + supplementary alignments. The original code filtered these out!

**Solution**: Now detects reads that have alignments on BOTH sides of a gap:
```
Primary alignment:      =====>        (ends at gap left boundary)
Supplementary:                 =====> (starts at gap right boundary)
                              ↑
                        Same read = crosses gap!
```

### 2. Smart Merge for Flanking Assemblies

**Problem**: Original flanking strategy always produced `left_seq + 500N + right_seq`, requiring multiple iterations.

**Solution**: After assembling left and right sides, attempt to find overlap and merge directly:
```
Before: [left_assembly][500N][right_assembly]  → needs more iterations
After:  [left_assembly===overlap===right_assembly] → complete in one step!
```

---

## Workflow Steps

### STEP 0: Gap Normalization

Standardize all N placeholders to 500bp before processing.

**Why this matters**:
- Original assembly may have various N lengths (500bp, 72kb, etc.)
- Normalizing ensures consistent spanning reads detection
- ONT reads can produce continuous alignment across short N regions

### STEP 1: Read Alignment

Align reads to current assembly using minimap2:
- HiFi: `-ax map-hifi`
- ONT: `-ax map-ont`
- Generate sorted, indexed BAM files

### STEP 2: Gap Detection

Scan assembly for N regions:
- Find all `N+` patterns ≥ 100bp
- Create gap records with coordinates

### STEP 3: Gap Filling (Three-Step Strategy)

For each gap:

**Step 3.1: Try Spanning Reads (Two Methods)**

```
Method A - Direct spanning:
  Read alignment:  |=================|
  Gap region:           |--500N--|
  Condition: ref_start <= gap_start AND ref_end >= gap_end

Method B - Supplementary-linked:
  Primary:         |======>        (ends near gap_start)
  Supplementary:          <======| (starts near gap_end)
  Same read name = true gap-crossing read!
```

If ≥3 spanning reads found → assemble with wtdbg2 → complete fill

**Step 3.2: Try Flanking Reads with Merge**

```
1. Collect left flanking reads (clipped at gap_start)
2. Collect right flanking reads (clipped at gap_end)
3. Assemble each side separately
4. TRY TO MERGE:
   - Use minimap2 to find overlap between left and right
   - If overlap ≥100bp with ≥90% identity → merge directly
   - Result: complete fill without 500N!
```

**Step 3.3: Flanking with 500N Placeholder (Fallback)**

If merge fails:
```
Final sequence: [left_assembly] + [500N] + [right_assembly]
```
This 500N will be processed in the next iteration.

### STEP 4: Assembly Update

Apply fills from end to start (reverse order) to prevent coordinate shift.

### STEP 5: Iteration Control

Continue if:
- Progress was made (new complete or partial fills)
- Partially filled gaps remain

Stop when:
- No progress made
- All gaps completely filled
- Max iterations reached

---

## Key Classes

### IterativeGapFiller (iterative_gapfiller.py)

Main orchestrator class.

```python
filler = IterativeGapFiller(
    assembly_file="input.fasta",
    hifi_reads="hifi.fastq.gz",
    ont_reads="ont.fastq.gz",
    output_dir="output",
    threads=8,
    max_iterations=10
)
final_assembly = filler.run()
```

### GapFiller (gap_filler.py) - OPTIMIZED

Core filling engine with three-step strategy.

```python
filler = GapFiller(
    assembly_file="assembly.fasta",
    hifi_bam="hifi.bam",
    ont_bam="ont.bam",
    threads=8,
    min_spanning_reads=3,
    min_overlap=100  # For merge detection
)
result = filler.fill_gap(gap_dict)
```

**Key Methods**:

| Method | Purpose |
|--------|---------|
| `fill_gap()` | Main entry - three-step strategy |
| `_try_spanning_reads_assembly()` | Step 1: spanning (direct + supplementary) |
| `_get_direct_spanning_reads()` | Method A: continuous alignment spans |
| `_get_supplementary_spanning_reads()` | Method B: split alignment detection |
| `_try_flanking_reads_with_merge()` | Step 2: flanking with merge attempt |
| `_try_merge_sequences()` | Overlap detection and merging |
| `_get_flanking_reads()` | Collect clipped reads (includes supplementary!) |

---

## Output Structure

```
output_dir/
├── assembly_normalized.fasta     # Gap-normalized assembly
├── final_assembly.fasta          # Final result
├── final_stats.json              # Statistics
├── validation_results.json       # Remaining gaps
├── summary.txt                   # Human-readable summary
│
└── iteration_N/
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    ├── iteration_stats.json
    └── work/
        └── gap_ChrX_start_end/
            ├── spanning_reads.fasta      # If spanning strategy used
            ├── left_reads.fasta          # If flanking strategy used
            ├── right_reads.fasta
            ├── merge_left.fa             # For merge attempt
            ├── merge_right.fa
            └── *_consensus.fa            # Assembly results
```

---

## Fill Strategies Summary

| Strategy | Condition | Result | Complete? |
|----------|-----------|--------|-----------|
| `spanning` | ≥3 reads span gap (direct or supplementary-linked) | Assembled sequence | Yes |
| `flanking_merged` | Left+right assemblies have overlap | Merged sequence | Yes |
| `flanking` | Left+right without overlap | left + 500N + right | No |
| Failed | No reads found | No change | - |

---

## Command Line Usage

```bash
python iterative_gapfiller.py \
    --assembly input.fasta \
    --hifi-reads hifi.fastq.gz \
    --ont-reads ont.fastq.gz \
    --output output_dir \
    --threads 8 \
    --max-iterations 10 \
    --min-gap-size 100 \
    --min-mapq 20 \
    --verbose
```

---

## Dependencies

### Required
- Python 3.8+
- BioPython, pysam
- minimap2, samtools, wtdbg2 (in PATH)

### Optional
- pyfaidx (faster sequence access)
- seqkit (for hard-clip read extraction)

---

## Expected Improvements

With the optimizations:

1. **More spanning reads detected**: Supplementary alignment method captures reads that were previously missed

2. **Fewer iterations needed**: Merge strategy can complete gaps in one round instead of multiple

3. **Better handling of large gaps**: Split alignments across large gaps are now properly detected

---

## Troubleshooting

### Still no spanning reads found
- Check read length vs gap size
- Very large real gaps may not have any crossing reads
- This is a biological limitation, not a software issue

### Merge not working
- Requires ≥100bp overlap with ≥90% identity
- If gap is truly large, left and right assemblies won't overlap
- Falls back to 500N placeholder correctly

### Low fill rate after multiple iterations
- Some gaps are in repeat regions
- Some gaps lack read coverage
- Consider using ultra-long ONT reads
