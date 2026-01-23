# Gap Filling Pipeline - Project Instructions

## Core Principles (MUST FOLLOW)

1. **Gap filling is iterative**: Multi-round, progressive process
2. **N length is meaningless**: The length of N sequences (placeholders) does NOT represent actual gap size. A 72kb N and 500bp N both just mean "there's a gap here"
3. **Gap normalization is critical**: All N placeholders must be normalized to 500bp before processing

## Quick Reference

- **Main entry**: `iterative_gapfiller.py` → `IterativeGapFiller` class
- **Fill engine**: `gap_filler.py` → `GapFiller` class
- **Full documentation**: See `WORKFLOW.md`

## File Structure

```
iterative_gapfiller.py  # Main orchestrator
gap_filler.py           # Core filling (two-step: spanning → flanking)
gap_validation.py       # Validation module
assembly_indexer.py     # Fast FASTA access (pyfaidx)
temp_file_manager.py    # Temp file management
```

## Key Workflow

```
STEP 0: Normalize all gaps to 500N (critical!)
        ↓
STEP 1: Align reads → BAM
        ↓
STEP 2: Find gaps (N regions)
        ↓
STEP 3: Fill each gap:
        3.1 Try spanning reads → complete fill
        3.2 Try flanking reads → partial fill with 500N
        3.3 Mark as failed if both fail
        ↓
STEP 4: Update assembly
        ↓
STEP 5: Repeat until no progress or max iterations
```

## Usage

```bash
python iterative_gapfiller.py \
    --assembly input.fasta \
    --hifi-reads hifi.fastq.gz \
    --ont-reads ont.fastq.gz \
    --output output_dir \
    --threads 8
```

## Dependencies

- minimap2, samtools, wtdbg2 (in PATH)
- BioPython, pysam, pyfaidx (Python)
