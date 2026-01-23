#!/usr/bin/env python3
"""
Gap Filler - Core gap filling with three-step strategy

OPTIMIZATIONS:
1. Utilizes supplementary alignments to detect reads spanning gaps
2. Attempts to merge left/right flanking assemblies
3. Falls back to flanking with 500N placeholder

Strategy:
  Step 1: Find spanning reads (direct + supplementary-linked)
  Step 2: Try flanking reads with smart merge
  Step 3: Fall back to flanking with 500N placeholder
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

import pysam
from Bio import SeqIO

from gapfill.utils.indexer import AssemblyIndexer


class GapFiller:
    """
    Optimized gap filler with three-step strategy
    """

    def __init__(self,
                 assembly_file: str,
                 hifi_bam: Optional[str] = None,
                 ont_bam: Optional[str] = None,
                 hifi_reads: Optional[str] = None,
                 ont_reads: Optional[str] = None,
                 threads: int = 8,
                 work_dir: Optional[str] = None,
                 flank_size: int = 500,
                 min_mapq: int = 20,
                 min_spanning_reads: int = 3,
                 min_overlap: int = 100):

        self.assembly_file = Path(assembly_file)
        self.hifi_bam = Path(hifi_bam) if hifi_bam else None
        self.ont_bam = Path(ont_bam) if ont_bam else None
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.threads = threads
        self.work_dir = Path(work_dir) if work_dir else Path(tempfile.mkdtemp())
        self.flank_size = flank_size
        self.min_mapq = min_mapq
        self.min_spanning_reads = min_spanning_reads
        self.min_overlap = min_overlap

        self.logger = logging.getLogger(__name__)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.assembly_indexer = AssemblyIndexer(assembly_file)
        self._read_sequence_cache: Dict[str, str] = {}

    def fill_gap(self, gap: Dict) -> Dict:
        """Fill a single gap using three-step strategy"""
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']
        gap_name = gap.get('name', f"{chrom}_{gap_start}_{gap_end}")

        self.logger.info(f"Filling gap {gap_name}: {chrom}:{gap_start}-{gap_end}")

        gap_work_dir = self.work_dir / gap_name
        gap_work_dir.mkdir(exist_ok=True)

        # Step 1: Spanning reads
        self.logger.info(f"  Step 1: Checking spanning reads...")
        spanning_result = self._try_spanning_reads_assembly(
            chrom, gap_start, gap_end, gap_name, gap_work_dir
        )

        if spanning_result['success']:
            self.logger.info(f"  ✓ Spanning: {len(spanning_result['sequence'])}bp")
            return spanning_result

        # Step 2: Flanking with merge
        self.logger.info(f"  Step 2: Checking flanking reads (with merge)...")
        flanking_result = self._try_flanking_reads_with_merge(
            chrom, gap_start, gap_end, gap_name, gap_work_dir
        )

        if flanking_result['success']:
            strategy = "merged" if flanking_result.get('is_complete') else "500N"
            self.logger.info(f"  ✓ Flanking ({strategy}): {len(flanking_result['sequence'])}bp")
            return flanking_result

        # Step 3: Failed
        self.logger.warning(f"  ✗ All strategies failed for {gap_name}")
        return {
            'success': False,
            'sequence': '',
            'strategy': None,
            'reason': 'No spanning or flanking reads found',
            'gap_name': gap_name
        }

    def _try_spanning_reads_assembly(self, chrom: str, gap_start: int, gap_end: int,
                                      gap_name: str, work_dir: Path) -> Dict:
        """Step 1: Find spanning reads (direct + supplementary-linked)"""
        spanning_reads = []
        seen_reads = set()

        for bam_path, source in [(self.hifi_bam, 'hifi'), (self.ont_bam, 'ont')]:
            if not bam_path or not bam_path.exists():
                continue

            # Method A: Direct spanning
            direct_reads = self._get_direct_spanning_reads(bam_path, chrom, gap_start, gap_end)
            for read_seq, read_name in direct_reads:
                if read_name not in seen_reads:
                    spanning_reads.append((read_seq, read_name, source))
                    seen_reads.add(read_name)

            # Method B: Supplementary-linked spanning
            supp_reads = self._get_supplementary_spanning_reads(bam_path, chrom, gap_start, gap_end)
            for read_seq, read_name in supp_reads:
                if read_name not in seen_reads:
                    spanning_reads.append((read_seq, read_name, source))
                    seen_reads.add(read_name)

        if len(spanning_reads) < self.min_spanning_reads:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'spanning',
                'reason': f'Insufficient spanning reads ({len(spanning_reads)} < {self.min_spanning_reads})'
            }

        reads_fasta = work_dir / "spanning_reads.fasta"
        with open(reads_fasta, 'w') as f:
            for i, (seq, name, source) in enumerate(spanning_reads):
                f.write(f">{name}__{source}__{i}\n{seq}\n")

        assembled_seq = self._run_wtdbg2_assembly(reads_fasta, work_dir, "spanning")

        if assembled_seq and len(assembled_seq) >= 100:
            return {
                'success': True,
                'sequence': assembled_seq,
                'strategy': 'spanning',
                'spanning_reads': len(spanning_reads),
                'gap_name': gap_name,
                'is_complete': True
            }
        else:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'spanning',
                'reason': 'wtdbg2 assembly failed'
            }

    def _get_direct_spanning_reads(self, bam_path: Path, chrom: str,
                                    gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """Method A: Get reads that directly span the gap"""
        spanning_reads = []

        try:
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                if chrom not in bam.references:
                    return []

                for read in bam.fetch(chrom, max(0, gap_start - 100), gap_end + 100):
                    if read.is_unmapped or read.is_secondary:
                        continue
                    if read.mapping_quality < self.min_mapq:
                        continue

                    if read.reference_start <= gap_start and read.reference_end >= gap_end:
                        seq = read.query_sequence
                        if seq and len(seq) >= 500:
                            spanning_reads.append((seq, read.query_name))

        except Exception as e:
            self.logger.warning(f"Error fetching direct spanning reads: {e}")

        return spanning_reads

    def _get_supplementary_spanning_reads(self, bam_path: Path, chrom: str,
                                           gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """Method B: Find reads spanning via supplementary alignments"""
        left_reads = {}
        right_reads = {}
        buffer = 1000

        try:
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                if chrom not in bam.references:
                    return []

                # Find reads ending near gap_start
                for read in bam.fetch(chrom, max(0, gap_start - buffer), gap_start + buffer):
                    if read.is_unmapped or read.mapping_quality < self.min_mapq:
                        continue
                    if read.reference_end and abs(read.reference_end - gap_start) <= buffer:
                        left_reads[read.query_name] = {
                            'ref_end': read.reference_end,
                            'seq': read.query_sequence
                        }

                # Find reads starting near gap_end
                for read in bam.fetch(chrom, max(0, gap_end - buffer), gap_end + buffer):
                    if read.is_unmapped or read.mapping_quality < self.min_mapq:
                        continue
                    if read.reference_start and abs(read.reference_start - gap_end) <= buffer:
                        right_reads[read.query_name] = {
                            'ref_start': read.reference_start,
                            'seq': read.query_sequence
                        }

        except Exception as e:
            self.logger.warning(f"Error fetching supplementary spanning reads: {e}")
            return []

        # Find reads on BOTH sides
        spanning_reads = []
        common_reads = set(left_reads.keys()) & set(right_reads.keys())

        for read_name in common_reads:
            left_info = left_reads[read_name]
            right_info = right_reads[read_name]

            if left_info['ref_end'] <= gap_start + buffer and right_info['ref_start'] >= gap_end - buffer:
                seq = left_info['seq'] or right_info['seq']
                if seq and len(seq) >= 500:
                    spanning_reads.append((seq, read_name))

        return spanning_reads

    def _try_flanking_reads_with_merge(self, chrom: str, gap_start: int, gap_end: int,
                                        gap_name: str, work_dir: Path) -> Dict:
        """Step 2: Flanking reads with merge attempt"""
        left_reads = self._get_flanking_reads(chrom, gap_start, gap_end, 'left')
        right_reads = self._get_flanking_reads(chrom, gap_start, gap_end, 'right')

        self.logger.info(f"    Found {len(left_reads)} left, {len(right_reads)} right flanking reads")

        if not left_reads and not right_reads:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'flanking',
                'reason': 'No flanking reads found'
            }

        # Assemble left side
        left_seq = ''
        if left_reads:
            left_fasta = work_dir / "left_reads.fasta"
            with open(left_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(left_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            left_seq = self._run_wtdbg2_assembly(left_fasta, work_dir, "left") or ''

        # Assemble right side
        right_seq = ''
        if right_reads:
            right_fasta = work_dir / "right_reads.fasta"
            with open(right_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(right_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            right_seq = self._run_wtdbg2_assembly(right_fasta, work_dir, "right") or ''

        if not left_seq and not right_seq:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'flanking',
                'reason': 'Assembly failed for both sides'
            }

        # Try to merge
        if left_seq and right_seq:
            merged_seq = self._try_merge_sequences(left_seq, right_seq, work_dir)

            if merged_seq:
                return {
                    'success': True,
                    'sequence': merged_seq,
                    'strategy': 'flanking_merged',
                    'gap_name': gap_name,
                    'is_complete': True,
                    'has_placeholder': False
                }

        # Fall back to 500N
        final_seq = left_seq + 'N' * 500 + right_seq

        return {
            'success': True,
            'sequence': final_seq,
            'strategy': 'flanking',
            'gap_name': gap_name,
            'is_complete': False,
            'has_placeholder': True
        }

    def _try_merge_sequences(self, left_seq: str, right_seq: str, work_dir: Path) -> Optional[str]:
        """Attempt to merge left and right sequences by finding overlap"""
        if len(left_seq) < self.min_overlap or len(right_seq) < self.min_overlap:
            return None

        try:
            left_fa = work_dir / "merge_left.fa"
            right_fa = work_dir / "merge_right.fa"

            with open(left_fa, 'w') as f:
                f.write(f">left\n{left_seq}\n")
            with open(right_fa, 'w') as f:
                f.write(f">right\n{right_seq}\n")

            result = subprocess.run(
                ['minimap2', '-x', 'asm5', '-c', str(left_fa), str(right_fa)],
                capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0 or not result.stdout.strip():
                return None

            for line in result.stdout.strip().split('\n'):
                fields = line.split('\t')
                if len(fields) < 12:
                    continue

                query_end = int(fields[3])
                strand = fields[4]
                target_len = int(fields[6])
                target_end = int(fields[8])
                matches = int(fields[9])
                block_len = int(fields[10])

                identity = matches / block_len if block_len > 0 else 0

                if identity >= 0.9 and block_len >= self.min_overlap:
                    if target_end >= target_len - 100 and strand == '+':
                        merged = left_seq + right_seq[query_end:]
                        return merged

            return None

        except Exception:
            return None

    def _get_flanking_reads(self, chrom: str, gap_start: int, gap_end: int,
                            side: str) -> List[Tuple[str, str, str]]:
        """Get reads that extend into the gap from one side"""
        flanking_reads = []
        window = 500
        seen_reads = set()

        for bam_path, source in [(self.hifi_bam, 'hifi'), (self.ont_bam, 'ont')]:
            if not bam_path or not bam_path.exists():
                continue

            try:
                with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                    if chrom not in bam.references:
                        continue

                    if side == 'left':
                        search_start = max(0, gap_start - window)
                        search_end = gap_start + window

                        for read in bam.fetch(chrom, search_start, search_end):
                            if read.is_unmapped or read.is_secondary:
                                continue
                            if read.query_name in seen_reads:
                                continue
                            if read.mapping_quality < self.min_mapq:
                                continue

                            should_include = False

                            if read.reference_end and abs(read.reference_end - gap_start) <= window:
                                cigar = read.cigartuples
                                if cigar and cigar[-1][0] in [4, 5]:
                                    should_include = True

                            if (read.reference_start < gap_start < read.reference_end and
                                read.reference_end < gap_end):
                                should_include = True

                            if read.is_supplementary and read.reference_end:
                                if abs(read.reference_end - gap_start) <= window:
                                    should_include = True

                            if should_include:
                                seq = read.query_sequence
                                if seq and len(seq) >= 500:
                                    flanking_reads.append((seq, read.query_name, source))
                                    seen_reads.add(read.query_name)

                    else:  # right
                        search_start = max(0, gap_end - window)
                        search_end = gap_end + window

                        for read in bam.fetch(chrom, search_start, search_end):
                            if read.is_unmapped or read.is_secondary:
                                continue
                            if read.query_name in seen_reads:
                                continue
                            if read.mapping_quality < self.min_mapq:
                                continue

                            should_include = False

                            if read.reference_start and abs(read.reference_start - gap_end) <= window:
                                cigar = read.cigartuples
                                if cigar and cigar[0][0] in [4, 5]:
                                    should_include = True

                            if (read.reference_start < gap_end < read.reference_end and
                                read.reference_start > gap_start):
                                should_include = True

                            if read.is_supplementary and read.reference_start:
                                if abs(read.reference_start - gap_end) <= window:
                                    should_include = True

                            if should_include:
                                seq = read.query_sequence
                                if seq and len(seq) >= 500:
                                    flanking_reads.append((seq, read.query_name, source))
                                    seen_reads.add(read.query_name)

            except Exception as e:
                self.logger.warning(f"Error fetching flanking reads: {e}")

        return flanking_reads

    def _run_wtdbg2_assembly(self, reads_fasta: Path, work_dir: Path,
                             prefix: str) -> Optional[str]:
        """Run wtdbg2 assembly"""
        output_prefix = work_dir / f"{prefix}_wtdbg2"

        try:
            if reads_fasta.stat().st_size < 100:
                return None

            read_count = 0
            total_bases = 0
            with open(reads_fasta) as f:
                for line in f:
                    if line.startswith('>'):
                        read_count += 1
                    else:
                        total_bases += len(line.strip())

            if read_count < 3:
                return None

            estimated_size = max(total_bases // 50, 10000)

            result = subprocess.run(
                ['wtdbg2', '-x', 'ccs', '-g', str(estimated_size),
                 '-t', str(self.threads), '-i', str(reads_fasta),
                 '-o', str(output_prefix), '-L', '500', '-e', '2'],
                capture_output=True, text=True, timeout=300
            )

            if result.returncode != 0:
                return None

            ctg_lay = Path(f"{output_prefix}.ctg.lay.gz")
            if not ctg_lay.exists():
                return None

            consensus_fa = work_dir / f"{prefix}_consensus.fa"
            result = subprocess.run(
                ['wtpoa-cns', '-t', str(self.threads),
                 '-i', str(ctg_lay), '-o', str(consensus_fa)],
                capture_output=True, text=True, timeout=300
            )

            if result.returncode != 0:
                return None

            if consensus_fa.exists():
                sequences = list(SeqIO.parse(consensus_fa, 'fasta'))
                if sequences:
                    best_seq = max(sequences, key=lambda x: len(x.seq))
                    return str(best_seq.seq)

            return None

        except Exception as e:
            self.logger.warning(f"wtdbg2 error: {e}")
            return None

    def close(self):
        if hasattr(self, 'assembly_indexer') and self.assembly_indexer:
            self.assembly_indexer.close()
