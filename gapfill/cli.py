"""
Command-line interface for GapFill

Usage:
    python -m gapfill -a assembly.fa --hifi reads.fq.gz -o output
    python -m gapfill -a hap1.fa hap2.fa --hifi reads.fq.gz -o output
"""

import argparse
import logging
import sys
from pathlib import Path


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        prog="gapfill",
        description="GapFill - Gap filling for haploid and polyploid genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Haploid
  python -m gapfill -a assembly.fa --hifi hifi.fq.gz -o output

  # Diploid (2 haplotypes)
  python -m gapfill -a hap1.fa hap2.fa --hifi hifi.fq.gz -o output

  # Tetraploid (4 haplotypes)
  python -m gapfill -a hap1.fa hap2.fa hap3.fa hap4.fa --ont ont.fq.gz -o output
"""
    )

    # Required
    parser.add_argument(
        "-a", "--assembly",
        nargs='+',
        required=True,
        help="Assembly file(s): 1 file = haploid, 2+ files = polyploid"
    )

    # Reads
    parser.add_argument("--hifi", dest="hifi_reads", help="HiFi reads (FASTQ/FASTA)")
    parser.add_argument("--ont", dest="ont_reads", help="ONT reads (FASTQ/FASTA)")
    parser.add_argument("--hic", dest="hic_reads", nargs=2, metavar=('R1', 'R2'),
                       help="Hi-C reads (R1 and R2 FASTQ files)")
    parser.add_argument("--hic-bam", dest="hic_bam", help="Pre-aligned Hi-C BAM file")

    # Output
    parser.add_argument("-o", "--output", default="gapfill_output", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Threads (default: 8)")

    # Gap filling
    parser.add_argument("--max-iterations", type=int, default=10, help="Max iterations")
    parser.add_argument("--min-gap-size", type=int, default=100, help="Min gap size")
    parser.add_argument("--min-mapq", type=int, default=20, help="Min MAPQ")

    # Polyploid
    parser.add_argument(
        "--phasing",
        choices=['builtin', 'whatshap'],
        default='builtin',
        help="Phasing method (default: builtin)"
    )
    parser.add_argument(
        "--no-ambiguous-reads",
        action="store_true",
        help="Exclude ambiguous reads in polyploid mode"
    )

    # Other
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    # Validate
    if not args.hifi_reads and not args.ont_reads:
        parser.error("At least one of --hifi or --ont is required")

    for asm in args.assembly:
        if not Path(asm).exists():
            parser.error(f"File not found: {asm}")

    # Run
    num_assemblies = len(args.assembly)

    if num_assemblies == 1:
        logger.info("=" * 60)
        logger.info("GapFill - HAPLOID MODE")
        logger.info("=" * 60)

        from gapfill.engines.haploid import HaploidEngine

        engine = HaploidEngine(
            assembly_file=args.assembly[0],
            hifi_reads=args.hifi_reads,
            ont_reads=args.ont_reads,
            hic_reads=args.hic_reads,
            hic_bam=args.hic_bam,
            output_dir=args.output,
            threads=args.threads,
            max_iterations=args.max_iterations,
            min_gap_size=args.min_gap_size,
            min_mapq=args.min_mapq
        )
        result = engine.run()
        logger.info(f"\nOutput: {result}")

    else:
        logger.info("=" * 60)
        logger.info(f"GapFill - POLYPLOID MODE ({num_assemblies} haplotypes)")
        logger.info("=" * 60)

        from gapfill.engines.polyploid import PolyploidEngine

        engine = PolyploidEngine(
            haplotype_assemblies=args.assembly,
            hifi_reads=args.hifi_reads,
            ont_reads=args.ont_reads,
            hic_reads=args.hic_reads,
            hic_bam=args.hic_bam,
            output_dir=args.output,
            threads=args.threads,
            max_iterations=args.max_iterations,
            phasing_method=args.phasing,
            use_ambiguous_reads=not args.no_ambiguous_reads
        )
        results = engine.run()

        logger.info("\nFilled assemblies:")
        for hap, path in results.items():
            logger.info(f"  {hap}: {path}")


if __name__ == "__main__":
    main()
