# GapFill

Gap filling pipeline for genome assemblies using HiFi and ONT long reads.
Supports both haploid and polyploid genomes.

## Installation

```bash
# Clone
git clone https://github.com/unavailable-2374/Gap_filling.git
cd Gap_filling

# Install dependencies
pip install biopython pysam pyfaidx numpy

# External tools (via conda)
conda install -c bioconda minimap2 samtools wtdbg
```

## Quick Start

```bash
# Haploid - single assembly file
python -m gapfill -a assembly.fa --hifi reads.fq.gz -o output

# Diploid - two haplotype files
python -m gapfill -a hap1.fa hap2.fa --hifi reads.fq.gz -o output

# Tetraploid - four haplotype files
python -m gapfill -a hap1.fa hap2.fa hap3.fa hap4.fa --ont reads.fq.gz -o output
```

## Usage

```
python -m gapfill [OPTIONS]

Required:
  -a, --assembly FILE(s)    Assembly file(s): 1 = haploid, 2+ = polyploid

Reads (at least one required):
  --hifi FILE               HiFi reads (FASTQ/FASTA)
  --ont FILE                ONT reads (FASTQ/FASTA)

Options:
  -o, --output DIR          Output directory (default: gapfill_output)
  -t, --threads N           Threads (default: 8)
  --max-iterations N        Max iterations (default: 10)
  --min-gap-size N          Min gap size to process (default: 100)
  --phasing METHOD          Phasing: builtin|whatshap (default: builtin)
  --no-ambiguous-reads      Exclude ambiguous reads (polyploid mode)
  -v, --verbose             Verbose output
```

## Features

- **Automatic mode detection** - Single file = haploid, multiple = polyploid
- **Gap normalization** - Standardizes N placeholders to 500bp
- **Optimized spanning detection** - Uses supplementary alignments
- **Smart merge** - Reduces iterations by merging flanking assemblies
- **Integrated phasing** - SNP-based read phasing for polyploid genomes

## Package Structure

```
gapfill/
├── __init__.py           # Package entry
├── __main__.py           # python -m gapfill
├── cli.py                # Command-line interface
├── core/
│   ├── filler.py         # GapFiller - core filling logic
│   └── validator.py      # GapValidator - validation
├── engines/
│   ├── haploid.py        # HaploidEngine
│   └── polyploid.py      # PolyploidEngine
└── utils/
    ├── indexer.py        # AssemblyIndexer
    ├── scanner.py        # GapScanner
    └── tempfiles.py      # TempFileManager
```

## Workflow

### Haploid Mode
```
Normalize gaps → Align reads → Find gaps → Fill gaps → Repeat
```

### Polyploid Mode
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

**Haploid:**
```
output/
├── final_assembly.fasta
├── final_stats.json
└── iteration_N/
```

**Polyploid:**
```
output/
├── hap1/final_assembly.fasta
├── hap2/final_assembly.fasta
├── snp_database.json
└── polyploid_summary.json
```

## Documentation

See [WORKFLOW.md](WORKFLOW.md) for detailed documentation.

## License

MIT License
