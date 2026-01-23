#!/usr/bin/env python3
"""
Gap Filler - Optimized gap filling with three-step strategy

OPTIMIZATIONS:
1. Utilizes supplementary alignments to detect reads spanning gaps
   - A read split across a gap produces primary + supplementary alignments
   - This identifies true gap-crossing reads that were previously missed

2. Attempts to merge left/right flanking assemblies
   - Checks for overlap between left and right sequences
   - If overlap found, merges directly without 500N placeholder
   - Reduces iteration count significantly

Strategy:
  Step 1: Find spanning reads (direct alignment spans + supplementary-linked)
  Step 2: Try flanking reads with smart merge
  Step 3: Fall back to flanking with 500N placeholder

Author: Gap Filling Pipeline
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from assembly_indexer import AssemblyIndexer


class GapFiller:
    """
    Optimized gap filler with three-step strategy:
    1. Spanning reads (including supplementary-linked) -> complete fill
    2. Flanking reads with merge attempt -> complete fill if overlap found
    3. Flanking reads with 500N placeholder -> partial fill
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
        """
        Initialize GapFiller

        Args:
            assembly_file: Path to assembly FASTA (full genome)
            hifi_bam: Path to HiFi BAM aligned to assembly
            ont_bam: Path to ONT BAM aligned to assembly
            hifi_reads: Path to original HiFi reads (for hard-clip extraction)
            ont_reads: Path to original ONT reads (for hard-clip extraction)
            threads: Number of threads for assembly
            work_dir: Working directory for temporary files
            flank_size: Flank size for read collection (default 500)
            min_mapq: Minimum mapping quality for reads (default 20)
            min_spanning_reads: Minimum spanning reads required (default 3)
            min_overlap: Minimum overlap for merging left/right sequences (default 100)
        """
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

        # Create work directory
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Assembly indexer for sequence access
        self.assembly_indexer = AssemblyIndexer(assembly_file)

        # Cache for read sequences
        self._read_sequence_cache: Dict[str, str] = {}

        self.logger.info(f"GapFiller initialized (optimized)")
        self.logger.info(f"  Assembly: {self.assembly_file}")
        self.logger.info(f"  HiFi BAM: {self.hifi_bam}")
        self.logger.info(f"  ONT BAM: {self.ont_bam}")

    def fill_gap(self, gap: Dict) -> Dict:
        """
        Fill a single gap using the optimized three-step strategy

        Args:
            gap: Gap dictionary with 'chrom', 'start', 'end', 'name'

        Returns:
            Fill result dictionary
        """
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']
        gap_name = gap.get('name', f"{chrom}_{gap_start}_{gap_end}")

        self.logger.info(f"Filling gap {gap_name}: {chrom}:{gap_start}-{gap_end}")

        # Create gap-specific work directory
        gap_work_dir = self.work_dir / gap_name
        gap_work_dir.mkdir(exist_ok=True)

        # =================================================================
        # STEP 1: Check for spanning reads (including supplementary-linked)
        # =================================================================
        self.logger.info(f"  Step 1: Checking for spanning reads...")
        spanning_result = self._try_spanning_reads_assembly(
            chrom, gap_start, gap_end, gap_name, gap_work_dir
        )

        if spanning_result['success']:
            self.logger.info(f"  ✓ Spanning reads assembly succeeded: {len(spanning_result['sequence'])}bp")
            return spanning_result

        self.logger.info(f"  Step 1 failed: {spanning_result.get('reason', 'No spanning reads')}")

        # =================================================================
        # STEP 2: Check for flanking reads with merge attempt
        # =================================================================
        self.logger.info(f"  Step 2: Checking for flanking reads (with merge)...")
        flanking_result = self._try_flanking_reads_with_merge(
            chrom, gap_start, gap_end, gap_name, gap_work_dir
        )

        if flanking_result['success']:
            strategy = "merged" if flanking_result.get('is_complete') else "500N placeholder"
            self.logger.info(f"  ✓ Flanking assembly succeeded ({strategy}): {len(flanking_result['sequence'])}bp")
            return flanking_result

        self.logger.info(f"  Step 2 failed: {flanking_result.get('reason', 'No flanking reads')}")

        # =================================================================
        # STEP 3: All strategies failed
        # =================================================================
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
        """
        Step 1: Find spanning reads using TWO methods:

        Method A: Direct spanning (read.ref_start <= gap_start AND read.ref_end >= gap_end)
        Method B: Supplementary-linked spanning (same read has alignments on both sides of gap)

        The supplementary method is crucial because:
        - When a read crosses a gap, minimap2 often splits it into primary + supplementary
        - Primary alignment ends at gap left boundary
        - Supplementary alignment starts at gap right boundary
        - These are the SAME READ and truly span the gap
        """
        self.logger.debug(f"    Collecting spanning reads for {gap_name}")

        spanning_reads = []
        seen_reads = set()

        for bam_path, source in [(self.hifi_bam, 'hifi'), (self.ont_bam, 'ont')]:
            if not bam_path or not bam_path.exists():
                continue

            # Method A: Direct spanning reads
            direct_reads = self._get_direct_spanning_reads(bam_path, chrom, gap_start, gap_end)
            for read_seq, read_name in direct_reads:
                if read_name not in seen_reads:
                    spanning_reads.append((read_seq, read_name, source))
                    seen_reads.add(read_name)

            # Method B: Supplementary-linked spanning reads
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
                'reason': f'Insufficient spanning reads ({len(spanning_reads)} < {self.min_spanning_reads})',
                'spanning_reads_count': len(spanning_reads)
            }

        self.logger.info(f"    Found {len(spanning_reads)} spanning reads")

        # Write reads to FASTA for wtdbg2
        reads_fasta = work_dir / "spanning_reads.fasta"
        with open(reads_fasta, 'w') as f:
            for i, (seq, name, source) in enumerate(spanning_reads):
                f.write(f">{name}__{source}__{i}\n{seq}\n")

        # Run wtdbg2 assembly
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
                'reason': 'wtdbg2 assembly failed or produced too short sequence',
                'spanning_reads_count': len(spanning_reads)
            }

    def _get_direct_spanning_reads(self, bam_path: Path, chrom: str,
                                    gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """
        Method A: Get reads that directly span the entire gap in alignment coordinates

        This works when:
        - The N region is short enough that minimap2 produces continuous alignment
        - Typically works for gaps < 1kb with long ONT reads
        """
        spanning_reads = []

        try:
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                if chrom not in bam.references:
                    return []

                for read in bam.fetch(chrom, max(0, gap_start - 100), gap_end + 100):
                    # Skip unmapped and secondary, but NOT supplementary (we handle those separately)
                    if read.is_unmapped or read.is_secondary:
                        continue

                    if read.mapping_quality < self.min_mapq:
                        continue

                    # Check if read spans the entire gap
                    if read.reference_start <= gap_start and read.reference_end >= gap_end:
                        seq = read.query_sequence
                        if seq and len(seq) >= 500:
                            spanning_reads.append((seq, read.query_name))

        except Exception as e:
            self.logger.warning(f"Error fetching direct spanning reads: {e}")

        return spanning_reads

    def _get_supplementary_spanning_reads(self, bam_path: Path, chrom: str,
                                           gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """
        Method B: Find reads that span the gap via supplementary alignments

        When a read truly crosses a gap:
        - minimap2 produces a primary alignment ending near gap_start
        - minimap2 produces a supplementary alignment starting near gap_end
        - Both alignments are from the SAME READ

        This is the KEY insight: supplementary alignments indicate gap-crossing reads!
        """
        # Track reads that have alignments on left and right of gap
        left_reads = {}   # read_name -> alignment info
        right_reads = {}  # read_name -> alignment info

        buffer = 1000  # Search buffer around gap boundaries

        try:
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                if chrom not in bam.references:
                    return []

                # Find reads ending near gap_start (left side)
                for read in bam.fetch(chrom, max(0, gap_start - buffer), gap_start + buffer):
                    if read.is_unmapped:
                        continue
                    if read.mapping_quality < self.min_mapq:
                        continue

                    # Check if alignment ends near gap_start
                    if read.reference_end and abs(read.reference_end - gap_start) <= buffer:
                        left_reads[read.query_name] = {
                            'ref_end': read.reference_end,
                            'seq': read.query_sequence,
                            'is_supplementary': read.is_supplementary
                        }

                # Find reads starting near gap_end (right side)
                for read in bam.fetch(chrom, max(0, gap_end - buffer), gap_end + buffer):
                    if read.is_unmapped:
                        continue
                    if read.mapping_quality < self.min_mapq:
                        continue

                    # Check if alignment starts near gap_end
                    if read.reference_start and abs(read.reference_start - gap_end) <= buffer:
                        right_reads[read.query_name] = {
                            'ref_start': read.reference_start,
                            'seq': read.query_sequence,
                            'is_supplementary': read.is_supplementary
                        }

        except Exception as e:
            self.logger.warning(f"Error fetching supplementary spanning reads: {e}")
            return []

        # Find reads that appear on BOTH sides of the gap
        # These are true gap-crossing reads!
        spanning_reads = []
        common_reads = set(left_reads.keys()) & set(right_reads.keys())

        self.logger.debug(f"    Found {len(common_reads)} reads with alignments on both sides of gap")

        for read_name in common_reads:
            left_info = left_reads[read_name]
            right_info = right_reads[read_name]

            # Verify the alignments are on correct sides
            if left_info['ref_end'] <= gap_start + buffer and right_info['ref_start'] >= gap_end - buffer:
                # Use the full sequence (either from primary or supplementary)
                seq = left_info['seq'] or right_info['seq']
                if seq and len(seq) >= 500:
                    spanning_reads.append((seq, read_name))

        return spanning_reads

    def _try_flanking_reads_with_merge(self, chrom: str, gap_start: int, gap_end: int,
                                        gap_name: str, work_dir: Path) -> Dict:
        """
        Step 2: Collect and assemble flanking reads, then attempt to merge

        OPTIMIZATION: After assembling left and right sides separately,
        check if they have overlapping sequences. If yes, merge them
        directly without the 500N placeholder.

        This can complete a gap in ONE iteration instead of multiple.
        """
        self.logger.debug(f"    Collecting flanking reads for {gap_name}")

        # Collect flanking reads
        left_reads = self._get_flanking_reads(chrom, gap_start, gap_end, 'left')
        right_reads = self._get_flanking_reads(chrom, gap_start, gap_end, 'right')

        self.logger.info(f"    Found {len(left_reads)} left, {len(right_reads)} right flanking reads")

        if not left_reads and not right_reads:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'flanking',
                'reason': 'No flanking reads found on either side'
            }

        # Assemble left side
        left_seq = ''
        if left_reads:
            left_fasta = work_dir / "left_reads.fasta"
            with open(left_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(left_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            left_seq = self._run_wtdbg2_assembly(left_fasta, work_dir, "left") or ''
            self.logger.info(f"    Left assembly: {len(left_seq)}bp" if left_seq else "    Left assembly failed")

        # Assemble right side
        right_seq = ''
        if right_reads:
            right_fasta = work_dir / "right_reads.fasta"
            with open(right_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(right_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            right_seq = self._run_wtdbg2_assembly(right_fasta, work_dir, "right") or ''
            self.logger.info(f"    Right assembly: {len(right_seq)}bp" if right_seq else "    Right assembly failed")

        if not left_seq and not right_seq:
            return {
                'success': False,
                'sequence': '',
                'strategy': 'flanking',
                'reason': 'Assembly failed for both sides',
                'left_reads_count': len(left_reads),
                'right_reads_count': len(right_reads)
            }

        # =================================================================
        # OPTIMIZATION: Try to merge left and right sequences
        # =================================================================
        if left_seq and right_seq:
            merged_seq = self._try_merge_sequences(left_seq, right_seq, work_dir)

            if merged_seq:
                self.logger.info(f"    ✓ Successfully merged left+right: {len(merged_seq)}bp")
                return {
                    'success': True,
                    'sequence': merged_seq,
                    'strategy': 'flanking_merged',
                    'left_length': len(left_seq),
                    'right_length': len(right_seq),
                    'merged_length': len(merged_seq),
                    'left_reads': len(left_reads),
                    'right_reads': len(right_reads),
                    'gap_name': gap_name,
                    'is_complete': True,  # Merged = complete fill!
                    'has_placeholder': False
                }

        # =================================================================
        # Fall back: Use 500N placeholder
        # =================================================================
        connector = 'N' * 500
        final_seq = left_seq + connector + right_seq

        self.logger.info(f"    Final flanking fill: {len(left_seq)}bp + 500N + {len(right_seq)}bp = {len(final_seq)}bp")

        return {
            'success': True,
            'sequence': final_seq,
            'strategy': 'flanking',
            'left_length': len(left_seq),
            'right_length': len(right_seq),
            'left_reads': len(left_reads),
            'right_reads': len(right_reads),
            'gap_name': gap_name,
            'is_complete': False,
            'has_placeholder': True
        }

    def _try_merge_sequences(self, left_seq: str, right_seq: str, work_dir: Path) -> Optional[str]:
        """
        Attempt to merge left and right sequences by finding overlap

        Strategy:
        1. Use minimap2 to align right_seq to left_seq
        2. If there's significant overlap, merge them
        3. Return merged sequence or None if no good overlap
        """
        if len(left_seq) < self.min_overlap or len(right_seq) < self.min_overlap:
            return None

        try:
            # Write sequences to temp files
            left_fa = work_dir / "merge_left.fa"
            right_fa = work_dir / "merge_right.fa"

            with open(left_fa, 'w') as f:
                f.write(f">left\n{left_seq}\n")
            with open(right_fa, 'w') as f:
                f.write(f">right\n{right_seq}\n")

            # Use minimap2 to find overlap
            result = subprocess.run(
                ['minimap2', '-x', 'asm5', '-c', str(left_fa), str(right_fa)],
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0 or not result.stdout.strip():
                return None

            # Parse PAF output
            for line in result.stdout.strip().split('\n'):
                fields = line.split('\t')
                if len(fields) < 12:
                    continue

                # PAF format: query_name, query_len, query_start, query_end, strand,
                #             target_name, target_len, target_start, target_end, matches, block_len, mapq
                query_len = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_len = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                block_len = int(fields[10])

                # Check for good overlap
                # Right sequence should align to the END of left sequence
                # with the START of right sequence
                identity = matches / block_len if block_len > 0 else 0

                if identity >= 0.9 and block_len >= self.min_overlap:
                    # Check if this is an end-to-start overlap
                    # target_end should be near left_seq end
                    # query_start should be near right_seq start
                    if target_end >= target_len - 100 and query_start <= 100:
                        if strand == '+':
                            # Merge: left_seq[:target_start] + right_seq
                            # Or: left_seq + right_seq[query_end:]
                            # Use the one that gives longer sequence
                            merged = left_seq + right_seq[query_end:]
                            self.logger.debug(f"    Merge found: {len(left_seq)}bp + {len(right_seq)}bp -> {len(merged)}bp (overlap: {block_len}bp)")
                            return merged

            return None

        except subprocess.TimeoutExpired:
            self.logger.debug("    Merge alignment timed out")
            return None
        except Exception as e:
            self.logger.debug(f"    Merge error: {e}")
            return None

    def _get_flanking_reads(self, chrom: str, gap_start: int, gap_end: int,
                            side: str) -> List[Tuple[str, str, str]]:
        """
        Get reads that extend into the gap from one side

        INCLUDES both primary and supplementary alignments - supplementary
        alignments often contain gap-crossing reads!
        """
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
                            # NOTE: We now INCLUDE supplementary alignments!
                            if read.query_name in seen_reads:
                                continue
                            if read.mapping_quality < self.min_mapq:
                                continue

                            should_include = False

                            # Case 1: Read ends near gap_start with right-side clip
                            if read.reference_end and abs(read.reference_end - gap_start) <= window:
                                cigar = read.cigartuples
                                if cigar and cigar[-1][0] in [4, 5]:
                                    should_include = True

                            # Case 2: Read crosses gap_start
                            if (read.reference_start < gap_start < read.reference_end and
                                read.reference_end < gap_end):
                                should_include = True

                            # Case 3: Supplementary alignment near gap boundary
                            if read.is_supplementary and read.reference_end:
                                if abs(read.reference_end - gap_start) <= window:
                                    should_include = True

                            if should_include:
                                seq = self._get_full_read_sequence(read)
                                if seq and len(seq) >= 500:
                                    flanking_reads.append((seq, read.query_name, source))
                                    seen_reads.add(read.query_name)

                    else:  # side == 'right'
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

                            # Case 1: Read starts near gap_end with left-side clip
                            if read.reference_start and abs(read.reference_start - gap_end) <= window:
                                cigar = read.cigartuples
                                if cigar and cigar[0][0] in [4, 5]:
                                    should_include = True

                            # Case 2: Read crosses gap_end
                            if (read.reference_start < gap_end < read.reference_end and
                                read.reference_start > gap_start):
                                should_include = True

                            # Case 3: Supplementary alignment near gap boundary
                            if read.is_supplementary and read.reference_start:
                                if abs(read.reference_start - gap_end) <= window:
                                    should_include = True

                            if should_include:
                                seq = self._get_full_read_sequence(read)
                                if seq and len(seq) >= 500:
                                    flanking_reads.append((seq, read.query_name, source))
                                    seen_reads.add(read.query_name)

            except Exception as e:
                self.logger.warning(f"Error fetching flanking reads from {bam_path}: {e}")

        return flanking_reads

    def _get_full_read_sequence(self, read: pysam.AlignedSegment) -> Optional[str]:
        """
        Get full read sequence, including hard-clipped portions if possible
        """
        cigar = read.cigartuples
        if not cigar:
            return read.query_sequence

        has_hard_clip = any(op == 5 for op, _ in cigar)

        if has_hard_clip:
            read_name = read.query_name.split()[0]

            if read_name in self._read_sequence_cache:
                return self._read_sequence_cache[read_name]

            if self.hifi_reads or self.ont_reads:
                extracted = self._extract_read_sequence(read_name)
                if extracted:
                    return extracted

            return read.query_sequence
        else:
            return read.query_sequence

    def _extract_read_sequence(self, read_name: str) -> Optional[str]:
        """
        Extract read sequence from original reads file
        """
        for reads_file in [self.hifi_reads, self.ont_reads]:
            if not reads_file or not reads_file.exists():
                continue

            try:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
                    f.write(f"{read_name}\n")
                    ids_file = f.name

                result = subprocess.run(
                    ['seqkit', 'grep', '-f', ids_file, str(reads_file)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                Path(ids_file).unlink()

                if result.returncode == 0 and result.stdout:
                    from io import StringIO
                    fmt = 'fastq' if str(reads_file).endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) else 'fasta'
                    for record in SeqIO.parse(StringIO(result.stdout), fmt):
                        seq = str(record.seq)
                        self._read_sequence_cache[read_name] = seq
                        return seq

            except FileNotFoundError:
                pass
            except Exception as e:
                self.logger.debug(f"Error extracting read {read_name}: {e}")

        return None

    def _run_wtdbg2_assembly(self, reads_fasta: Path, work_dir: Path,
                             prefix: str) -> Optional[str]:
        """
        Run wtdbg2 assembly on reads
        """
        output_prefix = work_dir / f"{prefix}_wtdbg2"

        try:
            reads_size = reads_fasta.stat().st_size
            if reads_size < 100:
                self.logger.debug(f"    Skipping wtdbg2: input too small ({reads_size} bytes)")
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
                self.logger.debug(f"    Skipping wtdbg2: too few reads ({read_count})")
                return None

            self.logger.debug(f"    wtdbg2 input: {read_count} reads, {total_bases:,} bp")

            estimated_size = max(total_bases // 50, 10000)

            wtdbg2_cmd = [
                'wtdbg2',
                '-x', 'ccs',
                '-g', str(estimated_size),
                '-t', str(self.threads),
                '-i', str(reads_fasta),
                '-o', str(output_prefix),
                '-L', '500',
                '-e', '2',
            ]

            result = subprocess.run(
                wtdbg2_cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                self.logger.debug(f"    wtdbg2 failed: {result.stderr[:200] if result.stderr else 'unknown'}")
                return None

            ctg_lay = Path(f"{output_prefix}.ctg.lay.gz")
            if not ctg_lay.exists():
                self.logger.debug(f"    wtdbg2 produced no contigs")
                return None

            consensus_fa = work_dir / f"{prefix}_consensus.fa"
            wtpoa_cmd = [
                'wtpoa-cns',
                '-t', str(self.threads),
                '-i', str(ctg_lay),
                '-o', str(consensus_fa)
            ]

            result = subprocess.run(
                wtpoa_cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                self.logger.debug(f"    wtpoa-cns failed")
                return None

            if consensus_fa.exists():
                sequences = list(SeqIO.parse(consensus_fa, 'fasta'))
                if sequences:
                    best_seq = max(sequences, key=lambda x: len(x.seq))
                    return str(best_seq.seq)

            return None

        except subprocess.TimeoutExpired:
            self.logger.warning(f"    wtdbg2 timed out")
            return None
        except Exception as e:
            self.logger.warning(f"    wtdbg2 error: {e}")
            return None

    def close(self):
        """Clean up resources"""
        if hasattr(self, 'assembly_indexer') and self.assembly_indexer:
            self.assembly_indexer.close()
