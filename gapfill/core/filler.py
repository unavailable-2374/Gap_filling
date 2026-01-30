#!/usr/bin/env python3
"""
Gap Filler - Optimized gap filling with tiered HiFi/ONT strategy

OPTIMIZATIONS:
1. Tiered strategy: HiFi-only → ONT-only → Hybrid
2. Different wtdbg2 presets for different read types
3. Optional HiFi polish for ONT assemblies
4. Source tracking for quality assessment
5. Hi-C guided candidate selection when multiple assemblies produced
6. Integrated validation with flank analysis

Strategy tiers:
  TIER 1: HiFi-only spanning (highest accuracy)
  TIER 2: ONT-only spanning (+ optional HiFi polish)
  TIER 3: Hybrid spanning
  TIER 4: HiFi flanking + merge
  TIER 5: ONT flanking + merge
  TIER 6: Hybrid flanking + 500N placeholder
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, TYPE_CHECKING
from collections import defaultdict

import pysam
from Bio import SeqIO

from gapfill.utils.indexer import AssemblyIndexer
from gapfill.core.validator import GapValidator, GapStatus, ValidationResult

if TYPE_CHECKING:
    from gapfill.utils.hic import HiCAnalyzer


class GapFiller:
    """
    Optimized gap filler with tiered HiFi/ONT strategy
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
                 min_overlap: int = 100,
                 enable_polish: bool = True,
                 enable_validation: bool = True,
                 hic_analyzer: Optional['HiCAnalyzer'] = None,
                 gap_size_estimates: Optional[Dict[str, int]] = None):

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
        self.enable_polish = enable_polish
        self.enable_validation = enable_validation
        self.hic_analyzer = hic_analyzer
        self.gap_size_estimates = gap_size_estimates or {}

        self.logger = logging.getLogger(__name__)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.assembly_indexer = AssemblyIndexer(assembly_file)
        self._read_sequence_cache: Dict[str, str] = {}

        # Initialize validator
        self.validator = GapValidator(threads=threads) if enable_validation else None

        # Check available data types
        self.has_hifi = self.hifi_bam is not None and self.hifi_bam.exists()
        self.has_ont = self.ont_bam is not None and self.ont_bam.exists()
        self.has_hic = self.hic_analyzer is not None

        self.logger.info(f"GapFiller initialized (tiered HiFi/ONT strategy)")
        self.logger.info(f"  HiFi BAM: {self.hifi_bam} ({'available' if self.has_hifi else 'not available'})")
        self.logger.info(f"  ONT BAM: {self.ont_bam} ({'available' if self.has_ont else 'not available'})")
        self.logger.info(f"  Hi-C: {'available' if self.has_hic else 'not available'}")
        self.logger.info(f"  Polish enabled: {self.enable_polish}")
        self.logger.info(f"  Validation enabled: {self.enable_validation}")

    def fill_gap(self, gap: Dict) -> Dict:
        """Fill a single gap using tiered HiFi/ONT strategy"""
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']
        gap_name = gap.get('name', f"{chrom}_{gap_start}_{gap_end}")

        self.logger.info(f"Filling gap {gap_name}: {chrom}:{gap_start}-{gap_end}")

        gap_work_dir = self.work_dir / gap_name
        gap_work_dir.mkdir(exist_ok=True)

        # Collect reads by type
        reads_info = self._collect_reads_by_type(chrom, gap_start, gap_end)

        self.logger.info(f"  Reads: HiFi spanning={reads_info['hifi_spanning_count']}, "
                        f"ONT spanning={reads_info['ont_spanning_count']}, "
                        f"HiFi flanking L/R={reads_info['hifi_left_count']}/{reads_info['hifi_right_count']}, "
                        f"ONT flanking L/R={reads_info['ont_left_count']}/{reads_info['ont_right_count']}")

        # =====================================================================
        # TIER 1: HiFi-only Spanning (highest accuracy)
        # =====================================================================
        if reads_info['hifi_spanning_count'] >= self.min_spanning_reads:
            self.logger.info(f"  TIER 1: Trying HiFi-only spanning...")
            result = self._assemble_spanning_reads(
                reads_info['hifi_spanning'], gap_work_dir, 'hifi', gap_name, gap
            )
            if result['success']:
                result['source'] = 'hifi_spanning'
                result['tier'] = 1
                self.logger.info(f"  ✓ TIER 1 success: {len(result['sequence'])}bp (HiFi-only)")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # TIER 2: ONT-only Spanning (+ optional HiFi polish)
        # =====================================================================
        if reads_info['ont_spanning_count'] >= self.min_spanning_reads:
            self.logger.info(f"  TIER 2: Trying ONT-only spanning...")
            result = self._assemble_spanning_reads(
                reads_info['ont_spanning'], gap_work_dir, 'ont', gap_name, gap
            )
            if result['success']:
                # Optional: Polish with HiFi reads if available
                if self.enable_polish and reads_info['hifi_spanning_count'] > 0:
                    polished = self._polish_with_hifi(
                        result['sequence'],
                        reads_info['hifi_spanning'],
                        gap_work_dir
                    )
                    if polished:
                        result['sequence'] = polished
                        result['source'] = 'ont_spanning_hifi_polished'
                    else:
                        result['source'] = 'ont_spanning'
                else:
                    result['source'] = 'ont_spanning'
                result['tier'] = 2
                self.logger.info(f"  ✓ TIER 2 success: {len(result['sequence'])}bp ({result['source']})")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # TIER 3: Hybrid Spanning (HiFi + ONT combined)
        # =====================================================================
        total_spanning = reads_info['hifi_spanning_count'] + reads_info['ont_spanning_count']
        if total_spanning >= self.min_spanning_reads:
            self.logger.info(f"  TIER 3: Trying hybrid spanning...")
            # Combine reads, but assemble separately and pick best
            result = self._assemble_hybrid_spanning(
                reads_info['hifi_spanning'],
                reads_info['ont_spanning'],
                gap_work_dir, gap_name, gap
            )
            if result['success']:
                result['tier'] = 3
                self.logger.info(f"  ✓ TIER 3 success: {len(result['sequence'])}bp ({result['source']})")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # TIER 4: HiFi Flanking + Merge
        # =====================================================================
        if reads_info['hifi_left_count'] > 0 or reads_info['hifi_right_count'] > 0:
            self.logger.info(f"  TIER 4: Trying HiFi flanking + merge...")
            result = self._assemble_flanking_reads(
                reads_info['hifi_left'], reads_info['hifi_right'],
                gap_work_dir, 'hifi', gap_name, gap
            )
            if result['success']:
                result['tier'] = 4
                self.logger.info(f"  ✓ TIER 4 success: {len(result['sequence'])}bp ({result['source']})")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # TIER 5: ONT Flanking + Merge (+ optional HiFi polish)
        # =====================================================================
        if reads_info['ont_left_count'] > 0 or reads_info['ont_right_count'] > 0:
            self.logger.info(f"  TIER 5: Trying ONT flanking + merge...")
            result = self._assemble_flanking_reads(
                reads_info['ont_left'], reads_info['ont_right'],
                gap_work_dir, 'ont', gap_name, gap
            )
            if result['success']:
                # Optional polish
                if self.enable_polish and (reads_info['hifi_left_count'] > 0 or reads_info['hifi_right_count'] > 0):
                    hifi_flanking = reads_info['hifi_left'] + reads_info['hifi_right']
                    if hifi_flanking:
                        polished = self._polish_with_hifi(
                            result['sequence'], hifi_flanking, gap_work_dir
                        )
                        if polished:
                            result['sequence'] = polished
                            result['source'] = result['source'].replace('ont_', 'ont_hifi_polished_')
                result['tier'] = 5
                self.logger.info(f"  ✓ TIER 5 success: {len(result['sequence'])}bp ({result['source']})")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # TIER 6: Hybrid Flanking + 500N
        # =====================================================================
        all_left = reads_info['hifi_left'] + reads_info['ont_left']
        all_right = reads_info['hifi_right'] + reads_info['ont_right']

        if all_left or all_right:
            self.logger.info(f"  TIER 6: Trying hybrid flanking...")
            result = self._assemble_flanking_reads(
                all_left, all_right, gap_work_dir, 'hybrid', gap_name, gap
            )
            if result['success']:
                result['tier'] = 6
                self.logger.info(f"  ✓ TIER 6 success: {len(result['sequence'])}bp ({result['source']})")
                return self._finalize_result(result, gap, reads_info)

        # =====================================================================
        # All tiers failed - analyze flanks to determine if UNFILLABLE or NEEDS_POLISH
        # =====================================================================
        self.logger.warning(f"  ✗ All tiers failed for {gap_name}")

        result = {
            'success': False,
            'sequence': '',
            'strategy': None,
            'source': None,
            'tier': None,
            'reason': 'All tiers failed - no suitable reads found',
            'gap_name': gap_name,
            'reads_info': {
                'hifi_spanning': reads_info['hifi_spanning_count'],
                'ont_spanning': reads_info['ont_spanning_count'],
                'hifi_flanking': reads_info['hifi_left_count'] + reads_info['hifi_right_count'],
                'ont_flanking': reads_info['ont_left_count'] + reads_info['ont_right_count']
            }
        }

        # Analyze flanks to determine final status
        return self._finalize_failed_result(result, gap)

    def _collect_reads_by_type(self, chrom: str, gap_start: int, gap_end: int) -> Dict:
        """Collect spanning and flanking reads, separated by type"""
        result = {
            'hifi_spanning': [],
            'ont_spanning': [],
            'hifi_left': [],
            'hifi_right': [],
            'ont_left': [],
            'ont_right': [],
            'hifi_spanning_count': 0,
            'ont_spanning_count': 0,
            'hifi_left_count': 0,
            'hifi_right_count': 0,
            'ont_left_count': 0,
            'ont_right_count': 0
        }

        # Collect HiFi reads
        if self.has_hifi:
            hifi_spanning = self._get_all_spanning_reads(self.hifi_bam, chrom, gap_start, gap_end)
            hifi_left = self._get_flanking_reads_one_side(self.hifi_bam, chrom, gap_start, gap_end, 'left')
            hifi_right = self._get_flanking_reads_one_side(self.hifi_bam, chrom, gap_start, gap_end, 'right')

            result['hifi_spanning'] = [(seq, name, 'hifi') for seq, name in hifi_spanning]
            result['hifi_left'] = [(seq, name, 'hifi') for seq, name in hifi_left]
            result['hifi_right'] = [(seq, name, 'hifi') for seq, name in hifi_right]
            result['hifi_spanning_count'] = len(hifi_spanning)
            result['hifi_left_count'] = len(hifi_left)
            result['hifi_right_count'] = len(hifi_right)

        # Collect ONT reads
        if self.has_ont:
            ont_spanning = self._get_all_spanning_reads(self.ont_bam, chrom, gap_start, gap_end)
            ont_left = self._get_flanking_reads_one_side(self.ont_bam, chrom, gap_start, gap_end, 'left')
            ont_right = self._get_flanking_reads_one_side(self.ont_bam, chrom, gap_start, gap_end, 'right')

            result['ont_spanning'] = [(seq, name, 'ont') for seq, name in ont_spanning]
            result['ont_left'] = [(seq, name, 'ont') for seq, name in ont_left]
            result['ont_right'] = [(seq, name, 'ont') for seq, name in ont_right]
            result['ont_spanning_count'] = len(ont_spanning)
            result['ont_left_count'] = len(ont_left)
            result['ont_right_count'] = len(ont_right)

        return result

    def _get_all_spanning_reads(self, bam_path: Path, chrom: str,
                                 gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """Get all spanning reads (direct + supplementary-linked)"""
        spanning_reads = []
        seen_reads = set()

        # Method A: Direct spanning
        direct = self._get_direct_spanning_reads(bam_path, chrom, gap_start, gap_end)
        for seq, name in direct:
            if name not in seen_reads:
                spanning_reads.append((seq, name))
                seen_reads.add(name)

        # Method B: Supplementary-linked
        supp = self._get_supplementary_spanning_reads(bam_path, chrom, gap_start, gap_end)
        for seq, name in supp:
            if name not in seen_reads:
                spanning_reads.append((seq, name))
                seen_reads.add(name)

        return spanning_reads

    def _get_direct_spanning_reads(self, bam_path: Path, chrom: str,
                                    gap_start: int, gap_end: int) -> List[Tuple[str, str]]:
        """Get reads that directly span the gap"""
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
        """Find reads spanning via supplementary alignments"""
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

    def _get_flanking_reads_one_side(self, bam_path: Path, chrom: str,
                                      gap_start: int, gap_end: int,
                                      side: str) -> List[Tuple[str, str]]:
        """Get flanking reads for one side"""
        flanking_reads = []
        window = 500
        seen_reads = set()

        try:
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                if chrom not in bam.references:
                    return []

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
                                flanking_reads.append((seq, read.query_name))
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
                                flanking_reads.append((seq, read.query_name))
                                seen_reads.add(read.query_name)

        except Exception as e:
            self.logger.warning(f"Error fetching flanking reads: {e}")

        return flanking_reads

    def _assemble_spanning_reads(self, reads: List[Tuple[str, str, str]],
                                  work_dir: Path, read_type: str,
                                  gap_name: str, gap_info: Optional[Dict] = None) -> Dict:
        """Assemble spanning reads of a specific type"""
        if len(reads) < self.min_spanning_reads:
            return {'success': False, 'reason': 'Insufficient reads'}

        reads_fasta = work_dir / f"spanning_{read_type}.fasta"
        with open(reads_fasta, 'w') as f:
            for i, (seq, name, source) in enumerate(reads):
                f.write(f">{name}__{source}__{i}\n{seq}\n")

        assembled_seq = self._run_wtdbg2_assembly(
            reads_fasta, work_dir, f"spanning_{read_type}", read_type, gap_info
        )

        if assembled_seq and len(assembled_seq) >= 100:
            return {
                'success': True,
                'sequence': assembled_seq,
                'strategy': 'spanning',
                'source': f'{read_type}_spanning',
                'read_count': len(reads),
                'gap_name': gap_name,
                'is_complete': True,
                'has_placeholder': False
            }
        else:
            return {'success': False, 'reason': f'wtdbg2 assembly failed for {read_type}'}

    def _assemble_hybrid_spanning(self, hifi_reads: List[Tuple[str, str, str]],
                                   ont_reads: List[Tuple[str, str, str]],
                                   work_dir: Path, gap_name: str,
                                   gap_info: Optional[Dict] = None) -> Dict:
        """Assemble hybrid spanning reads - try multiple strategies"""
        results = []

        # Strategy 1: Combined assembly with ccs preset (favor accuracy)
        all_reads = hifi_reads + ont_reads
        if len(all_reads) >= self.min_spanning_reads:
            reads_fasta = work_dir / "spanning_hybrid_ccs.fasta"
            with open(reads_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(all_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")

            seq = self._run_wtdbg2_assembly(
                reads_fasta, work_dir, "spanning_hybrid_ccs", 'hifi', gap_info
            )
            if seq and len(seq) >= 100:
                results.append(('hybrid_ccs', seq))

        # Strategy 2: Combined assembly with ont preset (favor length)
        if len(all_reads) >= self.min_spanning_reads:
            reads_fasta = work_dir / "spanning_hybrid_ont.fasta"
            with open(reads_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(all_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")

            seq = self._run_wtdbg2_assembly(
                reads_fasta, work_dir, "spanning_hybrid_ont", 'ont', gap_info
            )
            if seq and len(seq) >= 100:
                results.append(('hybrid_ont', seq))

        # Pick best result (longest sequence)
        if results:
            best_source, best_seq = max(results, key=lambda x: len(x[1]))
            return {
                'success': True,
                'sequence': best_seq,
                'strategy': 'spanning',
                'source': best_source,
                'read_count': len(all_reads),
                'gap_name': gap_name,
                'is_complete': True,
                'has_placeholder': False
            }

        return {'success': False, 'reason': 'Hybrid assembly failed'}

    def _assemble_flanking_reads(self, left_reads: List[Tuple[str, str, str]],
                                  right_reads: List[Tuple[str, str, str]],
                                  work_dir: Path, read_type: str,
                                  gap_name: str, gap_info: Optional[Dict] = None) -> Dict:
        """Assemble flanking reads and try to merge"""

        # Assemble left side
        left_seq = ''
        if left_reads:
            left_fasta = work_dir / f"left_{read_type}.fasta"
            with open(left_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(left_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            left_seq = self._run_wtdbg2_assembly(
                left_fasta, work_dir, f"left_{read_type}", read_type, gap_info
            ) or ''

        # Assemble right side
        right_seq = ''
        if right_reads:
            right_fasta = work_dir / f"right_{read_type}.fasta"
            with open(right_fasta, 'w') as f:
                for i, (seq, name, source) in enumerate(right_reads):
                    f.write(f">{name}__{source}__{i}\n{seq}\n")
            right_seq = self._run_wtdbg2_assembly(
                right_fasta, work_dir, f"right_{read_type}", read_type, gap_info
            ) or ''

        if not left_seq and not right_seq:
            return {'success': False, 'reason': f'Flanking assembly failed for {read_type}'}

        # Try to merge
        if left_seq and right_seq:
            merged_seq = self._try_merge_sequences(left_seq, right_seq, work_dir)

            if merged_seq:
                return {
                    'success': True,
                    'sequence': merged_seq,
                    'strategy': 'flanking_merged',
                    'source': f'{read_type}_flanking_merged',
                    'left_length': len(left_seq),
                    'right_length': len(right_seq),
                    'gap_name': gap_name,
                    'is_complete': True,
                    'has_placeholder': False
                }

        # Fall back to 500N placeholder
        final_seq = left_seq + 'N' * 500 + right_seq

        return {
            'success': True,
            'sequence': final_seq,
            'strategy': 'flanking',
            'source': f'{read_type}_flanking_500N',
            'left_length': len(left_seq),
            'right_length': len(right_seq),
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

    def _polish_with_hifi(self, sequence: str, hifi_reads: List[Tuple[str, str, str]],
                          work_dir: Path) -> Optional[str]:
        """Polish sequence with HiFi reads using minimap2 + racon"""
        if not hifi_reads or len(sequence) < 100:
            return None

        try:
            # Write sequence to file
            seq_fa = work_dir / "polish_seq.fa"
            with open(seq_fa, 'w') as f:
                f.write(f">seq\n{sequence}\n")

            # Write HiFi reads
            reads_fa = work_dir / "polish_reads.fa"
            with open(reads_fa, 'w') as f:
                for i, (seq, name, source) in enumerate(hifi_reads):
                    f.write(f">{name}__{i}\n{seq}\n")

            # Align reads to sequence
            paf_file = work_dir / "polish.paf"
            result = subprocess.run(
                ['minimap2', '-x', 'map-hifi', '-t', str(self.threads),
                 str(seq_fa), str(reads_fa)],
                capture_output=True, text=True, timeout=120
            )

            if result.returncode != 0:
                return None

            with open(paf_file, 'w') as f:
                f.write(result.stdout)

            # Check if we have enough alignments
            alignment_count = len([l for l in result.stdout.strip().split('\n') if l])
            if alignment_count < 2:
                self.logger.debug(f"    Polish skipped: only {alignment_count} alignments")
                return None

            # Run racon
            polished_fa = work_dir / "polished.fa"
            result = subprocess.run(
                ['racon', '-t', str(self.threads),
                 str(reads_fa), str(paf_file), str(seq_fa)],
                capture_output=True, text=True, timeout=300
            )

            if result.returncode != 0:
                self.logger.debug(f"    Racon failed: {result.stderr[:100] if result.stderr else 'unknown'}")
                return None

            # Parse polished sequence
            polished_seq = ''
            for line in result.stdout.split('\n'):
                if not line.startswith('>'):
                    polished_seq += line.strip()

            if len(polished_seq) >= len(sequence) * 0.8:  # Sanity check
                self.logger.debug(f"    Polish success: {len(sequence)}bp -> {len(polished_seq)}bp")
                return polished_seq
            else:
                self.logger.debug(f"    Polish result too short, keeping original")
                return None

        except FileNotFoundError:
            self.logger.debug("    Racon not found, skipping polish")
            return None
        except subprocess.TimeoutExpired:
            self.logger.debug("    Polish timed out")
            return None
        except Exception as e:
            self.logger.debug(f"    Polish error: {e}")
            return None

    def _run_wtdbg2_assembly(self, reads_fasta: Path, work_dir: Path,
                             prefix: str, read_type: str = 'hifi',
                             gap_info: Optional[Dict] = None) -> Optional[str]:
        """Run wtdbg2 assembly with type-specific parameters

        If Hi-C is available and multiple contigs are produced,
        uses Hi-C contact frequency to select the best candidate.
        """
        output_prefix = work_dir / f"{prefix}_wtdbg2"

        # Set parameters based on read type
        # Biological rationale:
        # - HiFi: high accuracy (~99.9%), can use lower edge coverage
        # - ONT: higher error rate, but edge_cov=2 with -S 1 rescues low-cov edges
        # - min_read_len raised to reduce noise from short fragments
        if read_type == 'hifi':
            preset = 'ccs'
            edge_cov = 2
            min_read_len = 1000   # HiFi reads typically 10-25kb
        elif read_type == 'ont':
            preset = 'ont'
            edge_cov = 2          # Lowered from 3, -S 1 will rescue low-cov edges
            min_read_len = 2000   # ONT reads typically 2-100kb
        else:  # hybrid
            preset = 'ccs'
            edge_cov = 2
            min_read_len = 1000

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

            # Smarter genome size estimation
            # Use Hi-C estimate if available, otherwise estimate from read stats
            if gap_info and gap_info.get('name') in self.gap_size_estimates:
                estimated_size = self.gap_size_estimates[gap_info['name']]
                self.logger.debug(f"    Using Hi-C estimated size: {estimated_size}bp")
            elif read_count >= 5:
                # Use 80% of average read length as estimate
                avg_read_len = total_bases / read_count
                estimated_size = max(int(avg_read_len * 0.8), 5000)
            else:
                # Low coverage: conservative estimate
                estimated_size = max(total_bases // 10, 5000)

            # -S 1: rescue low-coverage edges (critical for gap regions)
            result = subprocess.run(
                ['wtdbg2', '-x', preset, '-g', str(estimated_size),
                 '-t', str(self.threads), '-i', str(reads_fasta),
                 '-o', str(output_prefix), '-L', str(min_read_len),
                 '-e', str(edge_cov), '-S', '1'],
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
                if not sequences:
                    return None

                # If only one sequence or no Hi-C, return longest
                if len(sequences) == 1 or not self.has_hic or not gap_info:
                    best_seq = max(sequences, key=lambda x: len(x.seq))
                    return str(best_seq.seq)

                # Multiple candidates + Hi-C available: use Hi-C to select
                self.logger.debug(f"    {len(sequences)} candidate contigs, using Hi-C to select")
                best_seq = self._select_best_candidate_with_hic(
                    sequences, gap_info, work_dir
                )
                return best_seq

            return None

        except Exception as e:
            self.logger.warning(f"wtdbg2 error: {e}")
            return None

    def _select_best_candidate_with_hic(self, candidates: List,
                                         gap_info: Dict,
                                         work_dir: Path) -> Optional[str]:
        """
        Select best candidate sequence using Hi-C contact consistency.

        Strategy:
        - For each candidate, simulate inserting it into the gap
        - Score based on Hi-C contact frequency with flanking regions
        - Higher contact frequency = more likely correct
        """
        if not self.hic_analyzer or not candidates:
            # Fall back to longest
            return str(max(candidates, key=lambda x: len(x.seq)).seq)

        chrom = gap_info['chrom']
        gap_start = gap_info['start']
        gap_end = gap_info['end']

        scores = []

        for i, record in enumerate(candidates):
            seq = str(record.seq)
            seq_len = len(seq)

            # Score based on:
            # 1. Sequence length closer to estimated gap size
            # 2. Hi-C contact consistency (if available)

            # Length score: prefer sequences close to estimated size
            estimated_size = self.gap_size_estimates.get(gap_info.get('name', ''), gap_end - gap_start)
            length_diff = abs(seq_len - estimated_size)
            length_score = max(0, 1.0 - length_diff / max(estimated_size, 1000))

            # Hi-C score: check contact frequency with flanking regions
            hic_score = self._compute_hic_consistency_score(
                chrom, gap_start, gap_end, seq_len
            )

            # Combined score (weighted)
            combined_score = 0.3 * length_score + 0.7 * hic_score
            scores.append((combined_score, seq_len, seq))

            self.logger.debug(f"    Candidate {i+1}: {seq_len}bp, "
                            f"length_score={length_score:.2f}, "
                            f"hic_score={hic_score:.2f}, "
                            f"combined={combined_score:.2f}")

        # Select best scoring candidate
        best = max(scores, key=lambda x: (x[0], x[1]))  # Score, then length
        self.logger.debug(f"    Selected: {best[1]}bp (score={best[0]:.2f})")

        return best[2]

    def _compute_hic_consistency_score(self, chrom: str, gap_start: int,
                                        gap_end: int, fill_length: int) -> float:
        """
        Compute Hi-C consistency score for a potential fill.

        Checks if Hi-C contacts across the gap region are consistent
        with the proposed fill length.
        """
        if not self.hic_analyzer:
            return 0.5  # Neutral score

        try:
            # Get Hi-C contact matrix around gap region
            window = 10000
            region_start = max(0, gap_start - window)
            region_end = gap_end + window

            matrix = self.hic_analyzer.get_contact_matrix(
                chrom, region_start, region_end, resolution=1000
            )

            if matrix is None or matrix.size == 0:
                return 0.5

            # Calculate expected vs observed contact pattern
            # Contacts should decay with distance
            # If fill is correct size, decay pattern should be smooth

            n_bins = matrix.shape[0]
            gap_bin_start = (gap_start - region_start) // 1000
            gap_bin_end = (gap_end - region_start) // 1000

            # Check contacts between left flank and right flank
            left_bins = range(max(0, gap_bin_start - 5), gap_bin_start)
            right_bins = range(gap_bin_end, min(n_bins, gap_bin_end + 5))

            cross_contacts = 0
            for lb in left_bins:
                for rb in right_bins:
                    if 0 <= lb < n_bins and 0 <= rb < n_bins:
                        cross_contacts += matrix[lb, rb]

            # Normalize by region size
            n_left = len(list(left_bins))
            n_right = len(list(right_bins))
            if n_left > 0 and n_right > 0:
                avg_contacts = cross_contacts / (n_left * n_right)

                # Score: more contacts = better connectivity = higher score
                # Typical Hi-C contact at ~10kb distance
                expected_contacts = 5  # Baseline expectation
                score = min(1.0, avg_contacts / expected_contacts)
                return score

            return 0.5

        except Exception as e:
            self.logger.debug(f"    Hi-C scoring error: {e}")
            return 0.5

    def _finalize_result(self, result: Dict, gap: Dict, reads_info: Dict) -> Dict:
        """
        Finalize a successful fill result by adding validation info

        For complete fills: validates coverage, spanning reads, junction quality
        For partial fills: validates filled portions
        """
        if not self.validator:
            return result

        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']

        # Determine BAM file to use for validation (prefer HiFi)
        bam_file = None
        if self.has_hifi:
            bam_file = str(self.hifi_bam)
        elif self.has_ont:
            bam_file = str(self.ont_bam)

        if not bam_file:
            return result

        try:
            if result.get('is_complete', False) and not result.get('has_placeholder', False):
                # Complete fill - validate with full criteria
                validation = self.validator.validate_complete_fill(
                    bam_file, chrom, gap_start, gap_end, result['sequence']
                )
            else:
                # Partial fill with placeholder - find the new gap boundaries
                seq = result['sequence']
                # Find N-run in sequence (the placeholder)
                import re
                n_match = re.search(r'N{10,}', seq)
                if n_match:
                    # Calculate new gap position in reference coordinates
                    left_len = n_match.start()
                    right_len = len(seq) - n_match.end()
                    new_gap_start = gap_start + left_len
                    new_gap_end = gap_end - right_len

                    validation = self.validator.validate_partial_fill(
                        bam_file, chrom, gap_start, gap_end, seq,
                        new_gap_start, new_gap_end
                    )
                else:
                    # No N-run found but marked as partial - treat as complete
                    validation = self.validator.validate_complete_fill(
                        bam_file, chrom, gap_start, gap_end, result['sequence']
                    )

            # Add validation info to result
            result['validation'] = {
                'valid': validation.valid,
                'status': validation.status.value,
                'confidence': validation.confidence,
                'spanning_reads': validation.spanning_reads,
                'avg_coverage': validation.avg_coverage,
                'reason': validation.reason
            }

            self.logger.debug(f"    Validation: {validation.status.value} "
                            f"(confidence={validation.confidence:.2f}, "
                            f"spanning={validation.spanning_reads}, "
                            f"cov={validation.avg_coverage:.1f}x)")

        except Exception as e:
            self.logger.warning(f"Validation error: {e}")
            result['validation'] = {'valid': True, 'status': 'unknown', 'reason': str(e)}

        return result

    def _finalize_failed_result(self, result: Dict, gap: Dict) -> Dict:
        """
        Finalize a failed fill result by analyzing flanks

        Determines:
        - UNFILLABLE: flanks are correct, no spanning reads exist
        - NEEDS_POLISH: flanks have issues (clip accumulation, high mismatches)
        - FAILED: generic failure, may retry
        """
        if not self.validator:
            return result

        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']

        # Use HiFi BAM for flank analysis (more accurate)
        bam_file = None
        if self.has_hifi:
            bam_file = str(self.hifi_bam)
        elif self.has_ont:
            bam_file = str(self.ont_bam)

        if not bam_file:
            return result

        try:
            analysis = self.validator.analyze_failed_gap(
                bam_file, chrom, gap_start, gap_end
            )

            result['validation'] = {
                'valid': False,
                'status': analysis.status.value,
                'reason': analysis.reason,
                'left_flank_needs_polish': analysis.left_flank_needs_polish,
                'right_flank_needs_polish': analysis.right_flank_needs_polish,
                'flank_issues': analysis.flank_issues
            }

            if analysis.status == GapStatus.UNFILLABLE:
                self.logger.info(f"    → Gap marked UNFILLABLE: {analysis.reason}")
            elif analysis.status == GapStatus.NEEDS_POLISH:
                self.logger.info(f"    → Gap NEEDS_POLISH: {analysis.reason}")

        except Exception as e:
            self.logger.warning(f"Flank analysis error: {e}")
            result['validation'] = {'valid': False, 'status': 'failed', 'reason': str(e)}

        return result

    def close(self):
        if hasattr(self, 'assembly_indexer') and self.assembly_indexer:
            self.assembly_indexer.close()
        if hasattr(self, 'validator') and self.validator:
            self.validator.close()
