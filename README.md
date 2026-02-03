# GapFill

Gap filling tool for haploid and polyploid genome assemblies using long-read sequencing data (HiFi and/or ONT), with optional Hi-C integration.

## Features

- **Multi-round iterative filling** - Progressively fills gaps across multiple iterations
- **Haploid & Polyploid support** - Automatic mode detection based on input assemblies
- **HiFi + ONT integration** - 7-tier strategy prioritizing HiFi accuracy with ONT length advantage
- **Delayed validation** - Validates fills in next iteration using properly aligned reads
- **Automatic flank polishing** - Detects and polishes problematic flanking sequences
- **Hi-C data support** - Gap size estimation, candidate selection, and fill validation
- **Checkpoint/resume** - Resume interrupted runs using .ok file markers
- **Optimized polyploid mode** - Batch alignment reduces computation by 75%+ for polyploid genomes

### Performance Optimizations (v2.0)

- **Reads filtering** - Filter out reads anchored in non-gap regions, reducing alignment data by ~80%
- **Consensus-first** - Direct consensus for highly consistent HiFi reads, bypassing wtdbg2 assembly
- **Parallel gap filling** - Process multiple gaps concurrently using multiprocessing
- **High-confidence validation** - Immediate validation for high-quality fills, reducing iterations

## Installation

### Dependencies

**Python packages:**
```bash
pip install biopython pysam numpy pyfaidx
```

**External tools:**
```bash
# Required
conda install -c bioconda minimap2 samtools wtdbg

# Optional but recommended
conda install -c bioconda racon        # Flank polishing & ONT assembly polish
conda install -c bioconda bwa-mem2     # Hi-C alignment
conda install -c bioconda bcftools     # Polyploid SNP calling
conda install -c bioconda whatshap     # Alternative phasing method
```

### Install GapFill

```bash
git clone https://github.com/unavailable-2374/Gap_filling.git
cd Gap_filling
pip install -e .
```

## Quick Start

### Haploid genome

```bash
# HiFi only
python -m gapfill -a assembly.fa --hifi hifi.fq.gz -o output

# HiFi + ONT
python -m gapfill -a assembly.fa --hifi hifi.fq.gz --ont ont.fq.gz -o output

# With Hi-C data
python -m gapfill -a assembly.fa --hifi hifi.fq.gz --hic hic_R1.fq.gz hic_R2.fq.gz -o output

# Resume interrupted run
python -m gapfill -a assembly.fa --hifi hifi.fq.gz -o output --resume
```

### Polyploid genome (diploid example)

```bash
# Basic
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq.gz -o output

# Optimized mode (75% fewer alignments)
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq.gz --optimized -o output

# With ONT
python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq.gz --ont ont.fq.gz -o output
```

### Tetraploid genome

```bash
python -m gapfill -a hap1.fa hap2.fa hap3.fa hap4.fa --hifi hifi.fq.gz --optimized -o output
```

## Command Line Options

```
Required:
  -a, --assembly FILE(s)    Assembly file(s). 1 file = haploid, 2+ files = polyploid

Reads (at least one required):
  --hifi FILE               HiFi reads (FASTQ/FASTA, gzipped OK)
  --ont FILE                ONT reads (FASTQ/FASTA, gzipped OK)
  --hic R1 R2               Hi-C paired-end reads
  --hic-bam FILE            Pre-aligned Hi-C BAM file

Output:
  -o, --output DIR          Output directory (default: gapfill_output)
  -t, --threads N           Number of threads (default: 8)

Gap filling:
  --max-iterations N        Maximum iterations (default: 10)
  --min-gap-size N          Minimum gap size to process (default: 100)
  --min-mapq N              Minimum mapping quality (default: 20)

Polyploid:
  --phasing METHOD          Phasing method: builtin or whatshap (default: builtin)
  --no-ambiguous-reads      Exclude ambiguous reads in polyploid mode
  --optimized               Use batch alignment optimization (reduces alignments by 75%+)

Checkpoint:
  --resume                  Resume from checkpoint if available
  --clear-checkpoint        Clear existing checkpoint and start fresh

Performance:
  --no-filter-reads         Disable read filtering optimization (all modes)
  --no-parallel             Disable parallel gap filling (haploid only)

Other:
  -v, --verbose             Verbose output
  --version                 Show version
```

## How It Works

### Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  PREPROCESSING PHASE (once, before iterations)                  │
├─────────────────────────────────────────────────────────────────┤
│  1. Normalize gaps → 500N                                       │
│  2. Initial alignment → BAM                                     │
│  3. Pre-assess gap flanks                                       │
│     ├─ OK: ready for filling                                    │
│     └─ Issues: needs polishing                                  │
│  4. Polish problematic flanks                                   │
│  5. [OPT] Filter reads → keep only gap-related reads            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ITERATION LOOP                                                  │
├─────────────────────────────────────────────────────────────────┤
│  1. Align reads → BAM (use filtered reads if enabled)           │
│  2. Validate previous fills (iteration 2+)                      │
│     ├─ Complete fill: check spanning reads                      │
│     │   ├─ Pass → FILLED_COMPLETE                               │
│     │   └─ Fail → revert to 500N gap                            │
│     ├─ Partial fill: independent left/right validation          │
│     │   ├─ Both pass → FILLED_PARTIAL                           │
│     │   ├─ One pass → keep passed side, revert failed side      │
│     │   └─ Both fail → revert to 500N gap                       │
│     └─ 2b. Re-align reads if any reverts occurred               │
│  3. Find remaining gaps (using reverted assembly coordinates)   │
│  4. [OPT] Fill gaps in parallel                                 │
│     ├─ Tier 0: consensus-first (skip wtdbg2)                    │
│     ├─ High-confidence → immediate validation                   │
│     └─ Others → pending validation                              │
│  5. Apply fills to assembly (update downstream coordinates)     │
└─────────────────────────────────────────────────────────────────┘
```

### Core Principles

1. **N-runs are placeholders** - The length of N-runs in assemblies does not represent actual gap size
2. **Gap normalization** - All N-runs are normalized to 500bp to ensure consistent spanning read detection
3. **Delayed validation** - Fills are validated in the next iteration using reads aligned to the filled assembly
4. **Automatic flank polishing** - Problematic flanks (clip accumulation, high mismatches) are polished before filling

### 7-Tier Filling Strategy

GapFill uses a hierarchical approach prioritizing accuracy:

| Tier | Strategy | Description |
|------|----------|-------------|
| **0** | **Direct consensus** | HiFi spanning reads ≥5 with >95% identity → direct consensus (skip wtdbg2) |
| 1 | HiFi spanning | Spanning reads from HiFi assembled with wtdbg2 |
| 2 | ONT spanning | Spanning reads from ONT (length advantage) |
| 3 | Hybrid spanning | Combined HiFi+ONT spanning reads |
| 4 | HiFi flanking + merge | Assemble flanking HiFi reads, merge contigs |
| 5 | ONT flanking + merge | Assemble flanking ONT reads, polish with HiFi |
| 6 | 500N placeholder | Insert placeholder if filling fails |

**Tier 0 advantages:**
- Bypasses wtdbg2 assembly entirely for simple gaps
- Computes consensus from highly consistent reads directly
- Reduces gap filling time by ~50% for suitable gaps

### Gap Status Tracking

| Status | Description | Next Action |
|--------|-------------|-------------|
| PENDING | Not yet attempted | Attempt filling |
| NEEDS_POLISH | Flanks need polishing | Polish in preprocessing |
| FILLED_PENDING | Filled, awaiting validation | Validate next iteration |
| FILLED_COMPLETE | Completely filled and validated | Skip |
| FILLED_PARTIAL | Partially filled and validated | Skip |
| UNFILLABLE | Confirmed unfillable | Skip permanently |
| FAILED | Fill/validation failed | Retry next iteration |

### Validation Logic

**Complete fills** (spanning reads assembled):
- Require spanning reads covering the entire fill
- Short gaps (< 5kb): ≥3 spanning reads
- Long gaps (≥ 5kb): ≥5 spanning reads
- Fail → revert entire fill to 500N gap

**Partial fills** (flanking assembly with 500N placeholder):
- Independent left/right validation
- Check junction coverage (≥5x at connection points)
- Check insert coverage (≥5x average, no breakpoints)
- Left/right validated independently → keep passed side, revert failed side
- Both sides fail → revert to single 500N gap

**Re-alignment after revert:**
When validation fails and fills are reverted, reads are re-aligned to the reverted assembly before the next filling attempt. This ensures gap coordinates are accurate.

### Hi-C Integration

When Hi-C data is provided:

| Stage | Purpose | Method |
|-------|---------|--------|
| Pre-filling | Estimate true gap sizes | Contact frequency analysis |
| During filling | Select best candidate | Hi-C contact consistency scoring |
| Post-filling | Validate fills | Detect contact anomalies |
| Polyploid | Enhance phasing | Long-range haplotype links |

### Polyploid Optimization

**Reads Filtering (default enabled):**

By default, polyploid mode filters reads before phasing:
```
1. Align full reads to reference haplotype (once)
2. Filter out reads anchored in non-gap regions
3. Use filtered reads for phasing and gap filling
```

This reduces alignment data by ~80%, significantly speeding up the phasing step.

**Batch Alignment (--optimized flag):**

The `--optimized` flag enables batch alignment:

```
Standard mode:
  4-ploid × 10 iterations × 2 read types = 80 alignments

Optimized mode:
  10 iterations × 2 read types = 20 alignments (75% reduction)
```

This works because phased reads carry SNP signatures that naturally align to their correct haplotype in a merged reference.

| Ploidy | Iterations | Standard | Optimized | Reduction |
|--------|------------|----------|-----------|-----------|
| 2n | 10 | 40 | 20 | 50% |
| 4n | 10 | 80 | 20 | 75% |
| 6n | 10 | 120 | 20 | 83% |
| 8n | 10 | 160 | 20 | 87.5% |

## Package Structure

```
gapfill/
├── __init__.py
├── __main__.py
├── cli.py                        # Command-line interface
├── core/
│   ├── filler.py                 # GapFiller - 7-tier HiFi/ONT strategy
│   ├── validator.py              # GapValidator + GapStatusTracker
│   ├── polisher.py               # FlankPolisher - flank sequence polishing
│   ├── consensus.py              # ConsensusBuilder - direct consensus (Tier 0)
│   └── parallel.py               # Parallel gap filling support
├── engines/
│   ├── haploid.py                # HaploidEngine - haploid mode + Hi-C
│   ├── polyploid.py              # PolyploidEngine - standard polyploid
│   └── optimized_polyploid.py    # OptimizedPolyploidEngine - batch alignment
└── utils/
    ├── indexer.py                # AssemblyIndexer
    ├── scanner.py                # GapScanner
    ├── tempfiles.py              # TempFileManager
    ├── hic.py                    # HiCAnalyzer - Hi-C data analysis
    ├── checkpoint.py             # CheckpointManager - resume support
    └── reads_cache.py            # ReadsCache - filtered reads caching
```

## Output Structure

### Haploid

```
output/
├── checkpoint.json              # Checkpoint for resume
├── preprocessing/               # Preprocessing phase files
│   ├── hifi.bam, ont.bam
│   └── polish_work/
├── assembly_normalized.fasta    # Gap-normalized assembly
├── assembly_polished.fasta      # Polished assembly (if needed)
├── final_assembly.fasta         # Final filled assembly
├── final_stats.json             # Statistics
└── iteration_N/
    ├── assembly_filled.fasta
    ├── assembly_reverted.fasta  # If validation failures reverted
    ├── hifi.bam, ont.bam
    └── work/gap_*/              # Per-gap working files
```

### Polyploid

```
output/
├── checkpoint.json
├── snp_database.json            # SNP database for phasing
├── phased_*_hifi.fasta          # Phased HiFi reads per haplotype
├── phased_*_ont.fasta           # Phased ONT reads per haplotype
├── hap1/final_assembly.fasta
├── hap2/final_assembly.fasta
└── polyploid_summary.json
```

### Polyploid (Optimized)

```
output/
├── checkpoint.json
├── phase_hifi.bam, phase_ont.bam
├── iteration_N/
│   ├── merged_reference.fasta   # Combined reference
│   ├── merged_hifi.bam          # Single alignment
│   └── hap*/                    # Split by haplotype
├── hap1_filled.fasta
├── hap2_filled.fasta
└── summary.json
```

## Statistics Output

`final_stats.json` includes:

```json
{
  "iterations": 5,
  "total_gaps_initial": 150,
  "gaps_completely_filled": 120,
  "gaps_partially_filled": 20,
  "gaps_failed": 10,
  "gaps_polished": 15,
  "gap_status_summary": {
    "filled_complete": 120,
    "filled_partial": 20,
    "unfillable": 5,
    "failed": 5
  },
  "total_bp_filled": 1500000,
  "hic_validated": 115,
  "hic_failed": 5
}
```

## Performance Tips

### Haploid Mode Optimizations

By default, haploid mode enables all optimizations:

| Optimization | Effect | Flag to Disable |
|--------------|--------|-----------------|
| Reads filtering | Filters out reads anchored in non-gap regions → 80% less alignment data | `--no-filter-reads` |
| Parallel filling | Process multiple gaps concurrently → 3-5x speedup | `--no-parallel` |
| Consensus-first | Direct consensus for simple gaps → 50% fewer wtdbg2 runs | (always enabled) |
| High-confidence validation | Immediate validation for Tier 0/1 fills → 30% fewer iterations | (always enabled) |

**Expected speedup: 4-5x faster** (e.g., 74h → 16h for typical genome)

### General Tips

1. **Use `--optimized` for polyploid** - Significantly reduces runtime (75%+ fewer alignments)
2. **Provide both HiFi and ONT** - Combines accuracy and length advantages
3. **Add Hi-C data** - Improves fill quality and enables validation
4. **Use `--resume`** - Resume interrupted runs without recomputing
5. **Adjust `--max-iterations`** - More iterations may fill more gaps but with diminishing returns
6. **Use SSD storage** - Alignment I/O benefits from fast storage

### When to Disable Optimizations

- `--no-filter-reads`: If memory is very limited (filtering requires loading gap regions). Applies to all modes.
- `--no-parallel`: If running multiple samples simultaneously (avoid CPU competition). Haploid mode only.

## Troubleshooting

### Common Issues

**Out of memory during sorting:**
- The tool uses `-m 2G` per thread for samtools sort
- Reduce `--threads` if memory is limited

**Validation failures:**
- Check `assembly_reverted.fasta` for gaps that failed validation
- These gaps will be retried in subsequent iterations

**Many gaps marked UNFILLABLE:**
- These gaps have correct flanks but no spanning reads available
- May indicate true gaps that cannot be filled with available read lengths

**Resume not working:**
- Ensure `checkpoint.json` exists in output directory
- Use `--clear-checkpoint` to start fresh if checkpoint is corrupted

## Citation

If you use GapFill in your research, please cite:

```
[Citation information to be added]
```

## License

MIT License
