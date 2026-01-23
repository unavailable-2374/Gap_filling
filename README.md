# Gap Filling Pipeline

An optimized iterative gap filling pipeline for genome assemblies using HiFi and ONT long reads. Supports both haploid and polyploid genomes.

## Features

- **Automatic ploidy detection** - Single file = haploid, multiple files = polyploid
- **Iterative multi-round filling** - Progressive gap closure
- **HiFi + ONT support** - Works with both read types
- **Optimized spanning detection** - Uses supplementary alignments
- **Smart merge strategy** - Reduces iteration count
- **Integrated phasing** - Built-in or WhatsHap for polyploid

## Installation

```bash
# Python packages
pip install biopython pysam pyfaidx

# External tools (via conda)
conda install -c bioconda minimap2 samtools wtdbg
```

## Quick Start

```bash
# Haploid - single assembly file
python gapfill.py -a assembly.fa --hifi reads.fq.gz -o output

# Diploid - two haplotype files
python gapfill.py -a hap1.fa hap2.fa --hifi reads.fq.gz -o output

# Tetraploid - four haplotype files
python gapfill.py -a hap1.fa hap2.fa hap3.fa hap4.fa --ont reads.fq.gz -o output
```

## Usage

```bash
python gapfill.py [OPTIONS]

Required:
  -a, --assembly FILE(s)    Assembly FASTA file(s)
                            - 1 file  → haploid mode
                            - 2+ files → polyploid mode
  --hifi FILE               HiFi reads (at least one of --hifi/--ont required)
  --ont FILE                ONT reads

Options:
  -o, --output DIR          Output directory (default: gapfill_output)
  -t, --threads N           Threads (default: 8)
  --max-iterations N        Max iterations (default: 10)
  --min-gap-size N          Min gap size to process (default: 100)
  --phasing METHOD          Phasing method: builtin|whatshap (default: builtin)
  --no-ambiguous-reads      Exclude ambiguous reads in polyploid mode
  -v, --verbose             Verbose output
```

## Examples

```bash
# Haploid with both read types
python gapfill.py -a genome.fa --hifi hifi.fq.gz --ont ont.fq.gz -o filled -t 16

# Diploid with WhatsHap phasing
python gapfill.py -a hap1.fa hap2.fa --hifi reads.fq.gz --phasing whatshap -o diploid_out

# Hexaploid (6n)
python gapfill.py -a h1.fa h2.fa h3.fa h4.fa h5.fa h6.fa --hifi reads.fq.gz -o hexaploid_out
```

## How It Works

### Core Principle

N sequences are **placeholders only** - their length does NOT represent actual gap size.

### Haploid Workflow

```
Normalize gaps (500N) → Align reads → Find gaps → Fill gaps → Repeat
```

### Polyploid Workflow

```
Detect SNPs → Phase reads → Independent filling per haplotype
```

### Gap Filling Strategy

| Step | Method | Result |
|------|--------|--------|
| 1 | Spanning reads (direct + supplementary) | Complete fill |
| 2 | Flanking reads with merge | Complete fill |
| 3 | Flanking reads with 500N | Partial fill |

## Output

```
output/
├── final_assembly.fasta      # Filled assembly
├── final_stats.json          # Statistics
├── validation_results.json   # Remaining gaps
└── iteration_N/              # Per-iteration files
```

For polyploid:
```
output/
├── hap1/final_assembly.fasta
├── hap2/final_assembly.fasta
├── snp_database.json
└── polyploid_summary.json
```

## File Structure

```
Gap_filling/
├── gapfill.py                # ← UNIFIED ENTRY POINT
├── iterative_gapfiller.py    # Haploid engine
├── polyploid_gap_filler.py   # Polyploid engine
├── gap_filler.py             # Core filling logic
├── gap_validation.py         # Validation
├── assembly_indexer.py       # FASTA indexing
├── temp_file_manager.py      # Temp files
├── find_gaps.py              # Gap scanner
├── read_aligner.py           # Aligner wrapper
└── WORKFLOW.md               # Detailed docs
```

## Documentation

See [WORKFLOW.md](WORKFLOW.md) for detailed documentation.

## License

MIT License
