# Gap Filling Pipeline

An optimized iterative gap filling pipeline for genome assemblies using HiFi and ONT long reads.

## Features

- **Iterative multi-round gap filling** - Progressive gap closure across multiple iterations
- **Dual read support** - Works with both PacBio HiFi and Oxford Nanopore reads
- **Gap normalization** - Standardizes N placeholders to ensure consistent processing
- **Optimized spanning detection** - Utilizes supplementary alignments to identify gap-crossing reads
- **Smart merge strategy** - Attempts to merge flanking assemblies to reduce iteration count
- **wtdbg2-based assembly** - Fast and accurate local assembly for gap regions

## Installation

### Dependencies

**Required:**
- Python 3.8+
- BioPython
- pysam
- minimap2
- samtools
- wtdbg2

**Optional:**
- pyfaidx (10-1000x faster sequence access)
- seqkit (for hard-clip read extraction)

### Install Python packages

```bash
pip install biopython pysam pyfaidx
```

### Install external tools

```bash
# Using conda
conda install -c bioconda minimap2 samtools wtdbg

# Or from source
# See respective tool documentation
```

## Usage

### Basic Usage

```bash
python iterative_gapfiller.py \
    --assembly input.fasta \
    --hifi-reads hifi.fastq.gz \
    --ont-reads ont.fastq.gz \
    --output output_dir \
    --threads 8
```

### Full Options

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

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly` | Required | Input assembly FASTA file |
| `--hifi-reads` | Optional | HiFi reads (FASTQ/FASTA, gzipped OK) |
| `--ont-reads` | Optional | ONT reads (FASTQ/FASTA, gzipped OK) |
| `--output` | `output` | Output directory |
| `--threads` | 8 | Number of threads |
| `--max-iterations` | 10 | Maximum iterations |
| `--min-gap-size` | 100 | Minimum gap size to process |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--verbose` | False | Enable debug logging |

## How It Works

### Core Principle

The N sequences in an assembly are **placeholders only** - their length does NOT represent the actual biological gap size.

### Workflow

```
STEP 0: Normalize all gaps to 500N
           ↓
STEP 1: Align reads → BAM
           ↓
STEP 2: Find gaps (N regions)
           ↓
STEP 3: Fill each gap:
        3.1 Try spanning reads → complete fill
        3.2 Try flanking reads with merge → complete fill
        3.3 Flanking with 500N placeholder → partial fill
           ↓
STEP 4: Update assembly
           ↓
STEP 5: Repeat until no progress
```

### Gap Filling Strategies

| Strategy | Description | Result |
|----------|-------------|--------|
| **Spanning** | Reads that cross the entire gap (direct or via supplementary alignments) | Complete fill |
| **Flanking Merged** | Left and right assemblies with detectable overlap | Complete fill |
| **Flanking** | Left and right assemblies without overlap | Partial fill (500N placeholder) |

## Output

```
output_dir/
├── assembly_normalized.fasta   # Gap-normalized assembly
├── final_assembly.fasta        # Final filled assembly
├── final_stats.json            # Statistics
├── validation_results.json     # Remaining unfilled gaps
├── summary.txt                 # Human-readable summary
└── iteration_N/                # Per-iteration results
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    └── work/                   # Per-gap working files
```

## File Structure

```
Gap_filling/
├── iterative_gapfiller.py   # Main entry point
├── gap_filler.py            # Core filling engine
├── gap_validation.py        # Validation module
├── assembly_indexer.py      # Fast FASTA access
├── temp_file_manager.py     # Temp file management
├── find_gaps.py             # Standalone gap scanner
├── read_aligner.py          # Standalone aligner
├── WORKFLOW.md              # Detailed documentation
└── README.md                # This file
```

## Documentation

See [WORKFLOW.md](WORKFLOW.md) for detailed documentation including:
- Complete workflow explanation
- Key optimizations
- Class and method descriptions
- Troubleshooting guide

## License

MIT License

## Citation

If you use this tool in your research, please cite:

```
Gap Filling Pipeline - An optimized iterative gap filling tool for genome assemblies
https://github.com/unavailable-2374/Gap_filling
```
