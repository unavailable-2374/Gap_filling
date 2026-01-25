# GapFill

Gap filling tool for haploid and polyploid genome assemblies using long-read sequencing data (HiFi and/or ONT), with optional Hi-C integration.

## Features

- **Multi-round iterative filling** - Progressively fills gaps across multiple iterations
- **Haploid & Polyploid support** - Automatic mode detection based on input assemblies
- **HiFi + ONT integration** - 6-tier strategy prioritizing HiFi accuracy with ONT length advantage
- **Hi-C data support** - Gap size estimation, candidate selection, and fill validation
- **Optimized polyploid mode** - Batch alignment reduces computation by 75%+ for polyploid genomes

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
conda install -c bioconda racon        # HiFi polishing of ONT assemblies
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

Other:
  -v, --verbose             Verbose output
  --version                 Show version
```

## How It Works

### Core Principles

1. **N-runs are placeholders** - The length of N-runs in assemblies does not represent actual gap size
2. **Gap normalization** - All N-runs are normalized to 500bp to ensure consistent spanning read detection
3. **Iterative filling** - Each iteration re-aligns reads to the updated assembly
4. **Supplementary alignment detection** - Uses supplementary alignments to find true spanning reads

### 6-Tier Filling Strategy

GapFill uses a hierarchical approach prioritizing accuracy:

| Tier | Strategy | Description |
|------|----------|-------------|
| 1 | HiFi spanning | Direct spanning reads from HiFi (highest accuracy) |
| 2 | ONT spanning | Spanning reads from ONT (length advantage) |
| 3 | Hybrid spanning | Combined HiFi+ONT spanning reads |
| 4 | HiFi flanking + merge | Assemble flanking HiFi reads, merge contigs |
| 5 | ONT flanking + merge | Assemble flanking ONT reads, polish with HiFi |
| 6 | 500N placeholder | Insert placeholder if filling fails |

### Hi-C Integration

When Hi-C data is provided:

| Stage | Purpose | Method |
|-------|---------|--------|
| Pre-filling | Estimate true gap sizes | Contact frequency analysis |
| During filling | Select best candidate | Hi-C contact consistency scoring |
| Post-filling | Validate fills | Detect contact anomalies |
| Polyploid | Enhance phasing | Long-range haplotype links |

### Polyploid Optimization

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
│   ├── filler.py                 # GapFiller - 6-tier HiFi/ONT strategy + Hi-C
│   └── validator.py              # GapValidator - fill validation
├── engines/
│   ├── haploid.py                # HaploidEngine - haploid mode + Hi-C
│   ├── polyploid.py              # PolyploidEngine - standard polyploid
│   └── optimized_polyploid.py    # OptimizedPolyploidEngine - batch alignment
└── utils/
    ├── indexer.py                # AssemblyIndexer
    ├── scanner.py                # GapScanner
    ├── tempfiles.py              # TempFileManager
    └── hic.py                    # HiCAnalyzer - Hi-C data analysis
```

## Output Structure

### Haploid

```
output/
├── assembly_normalized.fasta    # Gap-normalized assembly
├── final_assembly.fasta         # Final filled assembly
├── final_stats.json             # Statistics
├── hic_aligned.bam              # Hi-C alignments (if provided)
└── iteration_N/
    ├── assembly_filled.fasta
    ├── hifi.bam, ont.bam
    └── work/gap_*/              # Per-gap working files
```

### Polyploid

```
output/
├── snp_database.json            # SNP database for phasing
├── phased_*_hifi.fasta          # Phased HiFi reads per haplotype
├── phased_*_ont.fasta           # Phased ONT reads per haplotype
├── hap1/final_assembly.fasta
├── hap2/final_assembly.fasta
└── summary.json
```

### Polyploid (Optimized)

```
output/
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
  "total_bp_filled": 1500000,
  "hic_validated": 115,
  "hic_failed": 5
}
```

## Performance Tips

1. **Use `--optimized` for polyploid** - Significantly reduces runtime
2. **Provide both HiFi and ONT** - Combines accuracy and length advantages
3. **Add Hi-C data** - Improves fill quality and enables validation
4. **Adjust `--max-iterations`** - More iterations may fill more gaps but with diminishing returns
5. **Use SSD storage** - Alignment I/O benefits from fast storage

## Citation

If you use GapFill in your research, please cite:

```
[Citation information to be added]
```

## License

MIT License
