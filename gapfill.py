#!/usr/bin/env python3
"""
GapFill - Unified gap filling pipeline for genome assemblies

Automatically detects ploidy based on input:
- Single assembly file → Haploid mode
- Multiple assembly files → Polyploid mode (with phasing)

Usage:
    # Haploid
    python gapfill.py -a assembly.fa --hifi reads.fq.gz -o output

    # Polyploid (diploid, tetraploid, etc.)
    python gapfill.py -a hap1.fa hap2.fa --hifi reads.fq.gz -o output

Author: Gap Filling Pipeline
"""

import argparse
import logging
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="GapFill - Unified gap filling for haploid and polyploid genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Haploid genome
  python gapfill.py -a assembly.fa --hifi hifi.fq.gz -o output

  # Diploid genome (2 haplotypes)
  python gapfill.py -a hap1.fa hap2.fa --hifi hifi.fq.gz -o output

  # Tetraploid genome (4 haplotypes)
  python gapfill.py -a hap1.fa hap2.fa hap3.fa hap4.fa --ont ont.fq.gz -o output

  # With WhatsHap phasing
  python gapfill.py -a hap1.fa hap2.fa --hifi hifi.fq.gz --phasing whatshap -o output
"""
    )

    # Required arguments
    parser.add_argument(
        "--assembly", "-a",
        nargs='+',
        required=True,
        help="Assembly FASTA file(s). Single file = haploid, multiple = polyploid"
    )

    # Reads
    parser.add_argument("--hifi-reads", "--hifi", help="HiFi reads (FASTQ/FASTA)")
    parser.add_argument("--ont-reads", "--ont", help="ONT reads (FASTQ/FASTA)")

    # Output
    parser.add_argument("--output", "-o", default="gapfill_output", help="Output directory")
    parser.add_argument("--threads", "-t", type=int, default=8, help="Threads (default: 8)")

    # Gap filling options
    parser.add_argument("--max-iterations", type=int, default=10, help="Max iterations (default: 10)")
    parser.add_argument("--min-gap-size", type=int, default=100, help="Min gap size (default: 100)")
    parser.add_argument("--min-mapq", type=int, default=20, help="Min MAPQ (default: 20)")

    # Polyploid options
    parser.add_argument(
        "--phasing",
        choices=['builtin', 'whatshap'],
        default='builtin',
        help="Phasing method for polyploid (default: builtin)"
    )
    parser.add_argument(
        "--no-ambiguous-reads",
        action="store_true",
        help="Exclude ambiguous reads in polyploid mode"
    )

    # Other
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    # Validate reads
    if not args.hifi_reads and not args.ont_reads:
        parser.error("At least one of --hifi-reads or --ont-reads is required")

    # Validate assembly files
    for asm in args.assembly:
        if not Path(asm).exists():
            parser.error(f"Assembly file not found: {asm}")

    # Determine mode
    num_assemblies = len(args.assembly)

    if num_assemblies == 1:
        # =========================
        # HAPLOID MODE
        # =========================
        logger.info("=" * 60)
        logger.info("GapFill - HAPLOID MODE")
        logger.info("=" * 60)

        from iterative_gapfiller import IterativeGapFiller

        filler = IterativeGapFiller(
            assembly_file=args.assembly[0],
            hifi_reads=args.hifi_reads,
            ont_reads=args.ont_reads,
            output_dir=args.output,
            threads=args.threads,
            max_iterations=args.max_iterations,
            min_gap_size=args.min_gap_size,
            min_mapq=args.min_mapq
        )

        result = filler.run()
        logger.info(f"\nOutput: {result}")

    else:
        # =========================
        # POLYPLOID MODE
        # =========================
        logger.info("=" * 60)
        logger.info(f"GapFill - POLYPLOID MODE ({num_assemblies} haplotypes)")
        logger.info("=" * 60)

        from polyploid_gap_filler import PolyploidGapFiller

        filler = PolyploidGapFiller(
            haplotype_assemblies=args.assembly,
            hifi_reads=args.hifi_reads,
            ont_reads=args.ont_reads,
            output_dir=args.output,
            threads=args.threads,
            max_iterations=args.max_iterations,
            phasing_method=args.phasing,
            use_ambiguous_reads=not args.no_ambiguous_reads
        )

        results = filler.run()

        logger.info("\nFilled assemblies:")
        for hap_name, path in results.items():
            logger.info(f"  {hap_name}: {path}")


if __name__ == "__main__":
    main()
