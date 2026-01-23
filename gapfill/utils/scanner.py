#!/usr/bin/env python3
"""
Gap Scanner - Detects N-runs in genome assemblies

Finds gaps (stretches of N's) and generates BED output.
"""

import re
import logging
from pathlib import Path
from typing import List, Dict, Optional

from gapfill.utils.indexer import AssemblyIndexer

logger = logging.getLogger(__name__)


class GapScanner:
    """
    Scans assembly for gaps (N-runs)

    Usage:
        scanner = GapScanner("assembly.fasta")
        gaps = scanner.find_gaps(min_size=100)
    """

    def __init__(self, assembly_file: str):
        self.assembly_file = Path(assembly_file)
        self.indexer = AssemblyIndexer(assembly_file)

    def find_gaps(self,
                  min_size: int = 1,
                  max_size: Optional[int] = None) -> List[Dict]:
        """
        Find all gaps in the assembly

        Args:
            min_size: Minimum gap size to report
            max_size: Maximum gap size to report (None = unlimited)

        Returns:
            List of gap dictionaries with chrom, start, end, name, size
        """
        gaps = []

        for chrom in self.indexer.get_all_chroms():
            seq = self.indexer.get_full_sequence(chrom).upper()

            for match in re.finditer(r'N+', seq):
                gap_start = match.start()
                gap_end = match.end()
                gap_size = gap_end - gap_start

                if gap_size < min_size:
                    continue
                if max_size and gap_size > max_size:
                    continue

                gaps.append({
                    'chrom': chrom,
                    'start': gap_start,
                    'end': gap_end,
                    'name': f"{chrom}_gap{len(gaps)+1}",
                    'size': gap_size
                })

        return gaps

    def write_bed(self, gaps: List[Dict], output_file: str):
        """Write gaps to BED file"""
        with open(output_file, 'w') as f:
            f.write("# Gaps in assembly\n")
            f.write("# chrom\tstart\tend\tname\tsize\n")

            for gap in gaps:
                f.write(f"{gap['chrom']}\t{gap['start']}\t{gap['end']}\t"
                       f"{gap['name']}\t{gap['size']}\n")

    def get_stats(self, gaps: List[Dict]) -> Dict:
        """Get gap statistics"""
        if not gaps:
            return {'count': 0, 'total_bp': 0}

        sizes = [g['size'] for g in gaps]
        return {
            'count': len(gaps),
            'total_bp': sum(sizes),
            'min_size': min(sizes),
            'max_size': max(sizes),
            'mean_size': sum(sizes) // len(sizes),
            'median_size': sorted(sizes)[len(sizes) // 2]
        }

    def close(self):
        self.indexer.close()


def find_gaps(assembly_file: str, output_bed: str,
              min_size: int = 1, max_size: Optional[int] = None) -> List[Dict]:
    """
    Convenience function to find gaps and write BED file

    Args:
        assembly_file: Input assembly FASTA
        output_bed: Output BED file
        min_size: Minimum gap size
        max_size: Maximum gap size

    Returns:
        List of gap dictionaries
    """
    scanner = GapScanner(assembly_file)
    gaps = scanner.find_gaps(min_size, max_size)
    scanner.write_bed(gaps, output_bed)

    stats = scanner.get_stats(gaps)
    logger.info(f"Found {stats['count']} gaps totaling {stats['total_bp']:,} bp")

    scanner.close()
    return gaps
