#!/usr/bin/env python3
"""
Find gaps in assembly and generate BED file
Detects N-runs (stretches of N's) in genome assembly

OPTIMIZED: Uses AssemblyIndexer for efficient sequence access
"""

import argparse
import re
from pathlib import Path

from assembly_indexer import AssemblyIndexer


def find_gaps(assembly_file, output_bed, min_gap_size=1, max_gap_size=None):
    """
    Find gaps (N-runs) in assembly and write to BED file

    Args:
        assembly_file: Input assembly FASTA file
        output_bed: Output BED file
        min_gap_size: Minimum gap size to report (default: 1)
        max_gap_size: Maximum gap size to report (default: None/unlimited)
    """
    gaps = []
    total_gaps = 0
    total_gap_bp = 0

    print(f"Scanning assembly: {assembly_file}")

    # OPTIMIZED: Use AssemblyIndexer
    indexer = AssemblyIndexer(assembly_file)

    for chrom in indexer.get_all_chroms():
        seq = indexer.get_full_sequence(chrom).upper()

        # Find all N-runs using regex
        for match in re.finditer(r'N+', seq):
            gap_start = match.start()
            gap_end = match.end()
            gap_size = gap_end - gap_start

            # Filter by gap size
            if gap_size < min_gap_size:
                continue
            if max_gap_size and gap_size > max_gap_size:
                continue

            gap_name = f"{chrom}_gap{len(gaps)+1}"

            gaps.append({
                'chrom': chrom,
                'start': gap_start,
                'end': gap_end,
                'name': gap_name,
                'size': gap_size
            })

            total_gaps += 1
            total_gap_bp += gap_size

    print(f"Found {total_gaps} gaps totaling {total_gap_bp:,} bp")

    # Write BED file
    with open(output_bed, 'w') as f:
        f.write("# Gaps found in assembly\n")
        f.write("# chrom\tstart\tend\tname\tsize\n")

        for gap in gaps:
            f.write(f"{gap['chrom']}\t{gap['start']}\t{gap['end']}\t"
                   f"{gap['name']}\t{gap['size']}\n")

    print(f"Gaps written to: {output_bed}")

    # Print statistics
    if gaps:
        gap_sizes = [g['size'] for g in gaps]
        print(f"\nGap size statistics:")
        print(f"  Minimum: {min(gap_sizes):,} bp")
        print(f"  Maximum: {max(gap_sizes):,} bp")
        print(f"  Mean:    {sum(gap_sizes)//len(gap_sizes):,} bp")
        print(f"  Median:  {sorted(gap_sizes)[len(gap_sizes)//2]:,} bp")

        # Gap size distribution
        size_bins = {
            '<100': 0,
            '100-500': 0,
            '500-1000': 0,
            '1000-5000': 0,
            '>5000': 0
        }

        for size in gap_sizes:
            if size < 100:
                size_bins['<100'] += 1
            elif size < 500:
                size_bins['100-500'] += 1
            elif size < 1000:
                size_bins['500-1000'] += 1
            elif size < 5000:
                size_bins['1000-5000'] += 1
            else:
                size_bins['>5000'] += 1

        print(f"\nGap size distribution:")
        for bin_name, count in size_bins.items():
            pct = count / len(gap_sizes) * 100
            print(f"  {bin_name:10s}: {count:5d} ({pct:5.1f}%)")

    return gaps


def main():
    parser = argparse.ArgumentParser(
        description="Find gaps (N-runs) in genome assembly and generate BED file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Find all gaps
  %(prog)s assembly.fasta -o gaps.bed

  # Find gaps between 100-10000 bp
  %(prog)s assembly.fasta -o gaps.bed --min-size 100 --max-size 10000

  # Find only large gaps (>1000 bp)
  %(prog)s assembly.fasta -o large_gaps.bed --min-size 1000
        """
    )

    parser.add_argument('assembly', help='Input assembly FASTA file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output BED file')
    parser.add_argument('--min-size', type=int, default=1,
                       help='Minimum gap size to report (default: 1)')
    parser.add_argument('--max-size', type=int, default=None,
                       help='Maximum gap size to report (default: unlimited)')

    args = parser.parse_args()

    find_gaps(args.assembly, args.output, args.min_size, args.max_size)


if __name__ == '__main__':
    main()
