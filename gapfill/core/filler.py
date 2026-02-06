#!/usr/bin/env python3
"""
Gap Filler - Optimized gap filling with tiered HiFi/ONT strategy

OPTIMIZATIONS:
1. Consensus-first: Skip wtdbg2 for high-consistency spanning reads
2. Tiered strategy: HiFi-only → ONT-only → Hybrid
3. Different wtdbg2 presets for different read types
4. Optional HiFi polish for ONT assemblies
5. Source tracking for quality assessment
6. Integrated validation with flank analysis

Strategy tiers:
  TIER 0: Direct consensus (for highly consistent spanning reads)
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
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

import pysam
from Bio import SeqIO

from gapfill.utils.indexer import AssemblyIndexer
from gapfill.core.validator import (GapValidator, GapStatus, ValidationResult,
                                     PartialFillValidationResult)
from gapfill.core.consensus import ConsensusBuilder, try_consensus_fill


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
                 enable_consensus_first: bool = True):

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
        self.enable_consensus_first = enable_consensus_first

        self.logger = logging.getLogger(__name__)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.assembly_indexer = AssemblyIndexer(assembly_file)
        self._read_sequence_cache: Dict[str, str] = {}

        # Initialize validator
        self.validator = GapValidator(threads=threads) if enable_validation else None

        # Check available data types
        self.has_hifi = self.hifi_bam is not None and self.hifi_bam.exists()
        self.has_ont = self.ont_bam is not None and self.ont_bam.exists()

        self.logger.info(f"GapFiller initialized (tiered HiFi/ONT strategy)")
        self.logger.info(f"  HiFi BAM: {self.hifi_bam} ({'available' if self.has_hifi else 'not available'})")
        self.logger.info(f"  ONT BAM: {self.ont_bam} ({'available' if self.has_ont else 'not available'})")
        self.logger.info(f"  Polish enabled: {self.enable_polish}")
        self.logger.info(f"  Validation enabled: {self.enable_validation}")

    def _validate_fill_locally(self, result, gap, work_dir):
        """
        Validate a fill locally before accepting it.

        Constructs a small local reference (flank + fill + flank), aligns nearby
        reads to it, and checks spanning/coverage metrics.

        Returns ValidationResult/PartialFillValidationResult, or None on error.
        """
        if not self.validator:
            return None

        is_complete = result.get('is_complete', False) and not result.get('has_placeholder', False)

        return self.validator.validate_fill_locally(
            hifi_bam=str(self.hifi_bam) if self.has_hifi else None,
            ont_bam=str(self.ont_bam) if self.has_ont else None,
            assembly_file=str(self.assembly_file),
            chrom=gap['chrom'],
            gap_start=gap['start'],
            gap_end=gap['end'],
            fill_sequence=result['sequence'],
            is_complete=is_complete,
            work_dir=str(work_dir / "local_validation"),
            threads=min(2, self.threads)
        )

    @staticmethod
    def _validation_to_dict(validation):
        """Convert a validation result to a serializable dict."""
        if validation is None:
            return {'valid': True, 'status': 'not_validated', 'reason': 'validation skipped'}

        from gapfill.core.validator import PartialFillValidationResult
        if isinstance(validation, PartialFillValidationResult):
            return {
                'valid': validation.valid,
                'status': validation.status.value,
                'left_valid': validation.left_valid,
                'right_valid': validation.right_valid,
                'left_fill_length': validation.left_fill_length,
                'right_fill_length': validation.right_fill_length,
                'avg_coverage': validation.avg_coverage,
                'reason': validation.reason,
            }
        else:
            return {
                'valid': validation.valid,
                'status': validation.status.value,
                'confidence': getattr(validation, 'confidence', 0.5),
                'spanning_reads': getattr(validation, 'spanning_reads', 0),
                'avg_coverage': getattr(validation, 'avg_coverage', 0),
                'reason': getattr(validation, 'reason', ''),
            }

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
        # TIER 0: Direct Consensus (skip wtdbg2 for highly consistent reads)
        # =====================================================================
        if self.enable_consensus_first and reads_info['hifi_spanning_count'] >= 3:
            self.logger.info(f"  TIER 0: Trying direct consensus...")
            consensus_result = try_consensus_fill(
                reads_info['hifi_spanning'], gap, self.threads
            )
            if consensus_result:
                consensus_result['tier'] = 0
                validation = self._validate_fill_locally(consensus_result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    consensus_result['validated'] = True
                    consensus_result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 0 success: {len(consensus_result['sequence'])}bp "
                                   f"({consensus_result['source']})")
                    return self._finalize_result(consensus_result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 0: validation failed: {validation.reason}, trying next...")

        # ONT consensus attempt
        if self.enable_consensus_first and reads_info['ont_spanning_count'] >= 5:
            self.logger.info(f"  TIER 0b: Trying ONT consensus...")
            consensus_result = try_consensus_fill(
                reads_info['ont_spanning'], gap, self.threads, read_type='ont'
            )
            if consensus_result:
                consensus_result['tier'] = 0
                validation = self._validate_fill_locally(consensus_result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    consensus_result['validated'] = True
                    consensus_result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 0b success: {len(consensus_result['sequence'])}bp "
                                   f"({consensus_result['source']})")
                    return self._finalize_result(consensus_result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 0b: validation failed: {validation.reason}, trying next...")

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
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 1 success: {len(result['sequence'])}bp (HiFi-only)")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 1: validation failed: {validation.reason}, trying next...")

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
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 2 success: {len(result['sequence'])}bp ({result['source']})")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 2: validation failed: {validation.reason}, trying next...")

        # =====================================================================
        # TIER 3: Hybrid Spanning (HiFi + ONT combined)
        # =====================================================================
        total_spanning = reads_info['hifi_spanning_count'] + reads_info['ont_spanning_count']
        if total_spanning >= self.min_spanning_reads:
            self.logger.info(f"  TIER 3: Trying hybrid spanning...")
            result = self._assemble_hybrid_spanning(
                reads_info['hifi_spanning'],
                reads_info['ont_spanning'],
                gap_work_dir, gap_name, gap
            )
            if result['success']:
                result['tier'] = 3
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 3 success: {len(result['sequence'])}bp ({result['source']})")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 3: validation failed: {validation.reason}, trying next...")

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
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 4 success: {len(result['sequence'])}bp ({result['source']})")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 4: validation failed: {validation.reason}, trying next...")

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
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 5 success: {len(result['sequence'])}bp ({result['source']})")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 5: validation failed: {validation.reason}, trying next...")

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
                validation = self._validate_fill_locally(result, gap, gap_work_dir)
                if validation is None or validation.valid:
                    result['validated'] = True
                    result['validation'] = self._validation_to_dict(validation)
                    self.logger.info(f"  ✓ TIER 6 success: {len(result['sequence'])}bp ({result['source']})")
                    return self._finalize_result(result, gap, reads_info)
                else:
                    self.logger.info(f"  TIER 6: validation failed: {validation.reason}, trying next...")

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
            'reason': 'All tiers failed - no suitable reads or validation failed',
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

    def _cluster_reads(self, reads: List[Tuple[str, str, str]],
                        work_dir: Path, prefix: str,
                        min_identity: float = 0.9) -> List[List[Tuple[str, str, str]]]:
        """
        Cluster reads by sequence similarity to separate reads from different sources.

        This prevents mixing reads from different haplotypes/paralogs during assembly,
        which can cause wtdbg2 to produce fragmented assemblies with internal Ns.

        Args:
            reads: List of (sequence, name, type) tuples
            work_dir: Working directory
            prefix: Prefix for temp files
            min_identity: Minimum identity threshold for clustering (default 0.9)

        Returns:
            List of read clusters, each cluster is a list of reads
        """
        if len(reads) <= 3:
            # Too few reads to cluster meaningfully
            return [reads]

        # Write reads to fasta
        reads_fasta = work_dir / f"{prefix}_for_cluster.fasta"
        with open(reads_fasta, 'w') as f:
            for i, (seq, name, source) in enumerate(reads):
                f.write(f">{i}\n{seq}\n")

        # Use minimap2 all-vs-all alignment to compute similarity
        try:
            result = subprocess.run(
                ['minimap2', '-x', 'ava-pb', '-t', str(min(4, self.threads)),
                 str(reads_fasta), str(reads_fasta)],
                capture_output=True, text=True, timeout=120
            )

            if result.returncode != 0:
                return [reads]

            # Parse PAF output to build similarity graph
            # PAF format: qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq ...
            similarities = defaultdict(dict)
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 11:
                    continue

                qname = int(fields[0])
                tname = int(fields[5])
                if qname == tname:
                    continue

                matches = int(fields[9])
                alnlen = int(fields[10])
                identity = matches / alnlen if alnlen > 0 else 0

                # Store max identity between pairs
                if tname not in similarities[qname] or identity > similarities[qname][tname]:
                    similarities[qname][tname] = identity
                    similarities[tname][qname] = identity

            # Simple greedy clustering
            n_reads = len(reads)
            assigned = [False] * n_reads
            clusters = []

            for i in range(n_reads):
                if assigned[i]:
                    continue

                # Start new cluster with read i
                cluster_indices = [i]
                assigned[i] = True

                # Add all reads similar to any read in the cluster
                changed = True
                while changed:
                    changed = False
                    for j in range(n_reads):
                        if assigned[j]:
                            continue
                        # Check if j is similar to any read in cluster
                        for k in cluster_indices:
                            if similarities[j].get(k, 0) >= min_identity:
                                cluster_indices.append(j)
                                assigned[j] = True
                                changed = True
                                break

                clusters.append([reads[idx] for idx in cluster_indices])

            # Sort clusters by size (largest first)
            clusters.sort(key=len, reverse=True)

            if len(clusters) > 1:
                self.logger.info(f"    Clustered {len(reads)} reads into {len(clusters)} groups: "
                               f"{[len(c) for c in clusters]}")

            return clusters

        except Exception as e:
            self.logger.debug(f"Clustering failed: {e}")
            return [reads]

    def _assemble_flanking_reads(self, left_reads: List[Tuple[str, str, str]],
                                  right_reads: List[Tuple[str, str, str]],
                                  work_dir: Path, read_type: str,
                                  gap_name: str, gap_info: Optional[Dict] = None) -> Dict:
        """Assemble flanking reads and try to merge"""

        # Assemble left side - with clustering to avoid mixed-source assembly
        left_seq = ''
        if left_reads:
            # Cluster reads to separate different sources (haplotypes, paralogs)
            left_clusters = self._cluster_reads(left_reads, work_dir, f"left_{read_type}")

            # Try assembling each cluster, use the best (longest clean) result
            best_left_seq = ''
            for ci, cluster in enumerate(left_clusters):
                if len(cluster) < 3:
                    continue

                left_fasta = work_dir / f"left_{read_type}_c{ci}.fasta"
                with open(left_fasta, 'w') as f:
                    for i, (seq, name, source) in enumerate(cluster):
                        f.write(f">{name}__{source}__{i}\n{seq}\n")

                assembled = self._run_wtdbg2_assembly(
                    left_fasta, work_dir, f"left_{read_type}_c{ci}", read_type, gap_info
                )

                if assembled and len(assembled) > len(best_left_seq):
                    best_left_seq = assembled

            left_seq = best_left_seq

        # Assemble right side - with clustering
        right_seq = ''
        if right_reads:
            right_clusters = self._cluster_reads(right_reads, work_dir, f"right_{read_type}")

            best_right_seq = ''
            for ci, cluster in enumerate(right_clusters):
                if len(cluster) < 3:
                    continue

                right_fasta = work_dir / f"right_{read_type}_c{ci}.fasta"
                with open(right_fasta, 'w') as f:
                    for i, (seq, name, source) in enumerate(cluster):
                        f.write(f">{name}__{source}__{i}\n{seq}\n")

                assembled = self._run_wtdbg2_assembly(
                    right_fasta, work_dir, f"right_{read_type}_c{ci}", read_type, gap_info
                )

                if assembled and len(assembled) > len(best_right_seq):
                    best_right_seq = assembled

            right_seq = best_right_seq

        if not left_seq and not right_seq:
            return {'success': False, 'reason': f'Flanking assembly failed for {read_type}'}

        # Trim flank overlap from assembled sequences
        # This removes the portion that duplicates the existing flank region
        if gap_info and (left_seq or right_seq):
            left_seq, right_seq = self._trim_flank_overlap(
                left_seq, right_seq, gap_info, work_dir
            )

        if not left_seq and not right_seq:
            return {'success': False, 'reason': f'No gap-extending sequence after trimming'}

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

    def _trim_flank_overlap(self, left_seq: str, right_seq: str,
                             gap_info: Dict, work_dir: Path) -> Tuple[str, str]:
        """
        Remove flank overlap from assembled sequences.

        The assembled sequences may contain portions that overlap with the
        existing flank regions. This method aligns the assembled sequences
        to the flanks and trims the overlapping portions.

        Args:
            left_seq: Assembled sequence from left flanking reads
            right_seq: Assembled sequence from right flanking reads
            gap_info: Gap information dict with 'chrom', 'start', 'end'
            work_dir: Working directory for temporary files

        Returns:
            Tuple of (trimmed_left_seq, trimmed_right_seq)
        """
        chrom = gap_info['chrom']
        gap_start = gap_info['start']
        gap_end = gap_info['end']

        # Get flank sequences from assembly
        flank_size = 2000  # Check 2kb of flank for overlap
        left_flank = self.assembly_indexer.get_sequence(
            chrom, max(0, gap_start - flank_size), gap_start
        )
        right_flank = self.assembly_indexer.get_sequence(
            chrom, gap_end, gap_end + flank_size
        )

        trimmed_left = left_seq
        trimmed_right = right_seq

        # Trim left_seq: remove the part that overlaps with left flank
        if left_seq and left_flank:
            trimmed_left = self._trim_left_overlap(left_seq, left_flank, work_dir)

        # Trim right_seq: remove the part that overlaps with right flank
        if right_seq and right_flank:
            trimmed_right = self._trim_right_overlap(right_seq, right_flank, work_dir)

        return trimmed_left, trimmed_right

    def _trim_left_overlap(self, assembled_seq: str, flank_seq: str,
                           work_dir: Path) -> str:
        """
        Trim the portion of assembled_seq that overlaps with left flank.

        For left side assembly, the overlap is at the START of assembled_seq.
        We want to keep only the part that extends INTO the gap.

        Example:
            flank_seq:     ...ACGTACGT (ends at gap_start)
            assembled_seq: ACGTACGT XXXXXXX (starts with flank, extends into gap)
                          |overlap| |keep |
            result:                 XXXXXXX
        """
        if len(assembled_seq) < 100 or len(flank_seq) < 100:
            return assembled_seq

        try:
            # Write sequences to files
            assembled_fa = work_dir / "trim_left_assembled.fa"
            flank_fa = work_dir / "trim_left_flank.fa"

            with open(assembled_fa, 'w') as f:
                f.write(f">assembled\n{assembled_seq}\n")
            with open(flank_fa, 'w') as f:
                f.write(f">flank\n{flank_seq}\n")

            # Align assembled to flank
            result = subprocess.run(
                ['minimap2', '-x', 'asm5', '-c', str(flank_fa), str(assembled_fa)],
                capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0 or not result.stdout.strip():
                return assembled_seq

            # Parse PAF output to find overlap
            # PAF format: query_name, query_len, query_start, query_end, strand,
            #             target_name, target_len, target_start, target_end, ...
            best_query_end = 0

            for line in result.stdout.strip().split('\n'):
                fields = line.split('\t')
                if len(fields) < 12:
                    continue

                query_start = int(fields[2])  # assembled_seq start
                query_end = int(fields[3])    # assembled_seq end
                strand = fields[4]
                target_end = int(fields[8])   # flank end
                target_len = int(fields[6])   # flank length
                matches = int(fields[9])
                block_len = int(fields[10])

                identity = matches / block_len if block_len > 0 else 0

                # We want alignments where:
                # 1. The flank alignment reaches near the end of flank (near gap_start)
                # 2. The assembled_seq alignment starts near the beginning
                # 3. Good identity
                if (identity >= 0.9 and
                    target_end >= target_len - 100 and  # Flank alignment near gap
                    query_start < 200 and               # Assembled starts from beginning
                    strand == '+'):
                    if query_end > best_query_end:
                        best_query_end = query_end

            if best_query_end > 0 and best_query_end < len(assembled_seq) - 100:
                # Trim the overlapping portion
                trimmed = assembled_seq[best_query_end:]
                self.logger.debug(f"    Trimmed left overlap: {len(assembled_seq)}bp -> {len(trimmed)}bp")
                return trimmed

            return assembled_seq

        except Exception as e:
            self.logger.debug(f"    Left trim error: {e}")
            return assembled_seq

    def _trim_right_overlap(self, assembled_seq: str, flank_seq: str,
                            work_dir: Path) -> str:
        """
        Trim the portion of assembled_seq that overlaps with right flank.

        For right side assembly, the overlap is at the END of assembled_seq.
        We want to keep only the part that extends INTO the gap.

        Example:
            assembled_seq: XXXXXXX GCTAGCTA (ends with flank, starts in gap)
                          | keep | |overlap|
            flank_seq:             GCTAGCTA... (starts at gap_end)
            result:        XXXXXXX
        """
        if len(assembled_seq) < 100 or len(flank_seq) < 100:
            return assembled_seq

        try:
            # Write sequences to files
            assembled_fa = work_dir / "trim_right_assembled.fa"
            flank_fa = work_dir / "trim_right_flank.fa"

            with open(assembled_fa, 'w') as f:
                f.write(f">assembled\n{assembled_seq}\n")
            with open(flank_fa, 'w') as f:
                f.write(f">flank\n{flank_seq}\n")

            # Align assembled to flank
            result = subprocess.run(
                ['minimap2', '-x', 'asm5', '-c', str(flank_fa), str(assembled_fa)],
                capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0 or not result.stdout.strip():
                return assembled_seq

            # Parse PAF output to find overlap
            best_query_start = len(assembled_seq)
            assembled_len = len(assembled_seq)

            for line in result.stdout.strip().split('\n'):
                fields = line.split('\t')
                if len(fields) < 12:
                    continue

                query_start = int(fields[2])  # assembled_seq start
                query_end = int(fields[3])    # assembled_seq end
                query_len = int(fields[1])    # assembled_seq length
                strand = fields[4]
                target_start = int(fields[7]) # flank start
                matches = int(fields[9])
                block_len = int(fields[10])

                identity = matches / block_len if block_len > 0 else 0

                # We want alignments where:
                # 1. The flank alignment starts near the beginning of flank (near gap_end)
                # 2. The assembled_seq alignment ends near the end
                # 3. Good identity
                if (identity >= 0.9 and
                    target_start < 100 and              # Flank alignment near gap
                    query_end > query_len - 200 and     # Assembled ends near end
                    strand == '+'):
                    if query_start < best_query_start:
                        best_query_start = query_start

            if best_query_start > 100 and best_query_start < assembled_len:
                # Trim the overlapping portion
                trimmed = assembled_seq[:best_query_start]
                self.logger.debug(f"    Trimmed right overlap: {len(assembled_seq)}bp -> {len(trimmed)}bp")
                return trimmed

            return assembled_seq

        except Exception as e:
            self.logger.debug(f"    Right trim error: {e}")
            return assembled_seq

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

        Returns the longest contig sequence if multiple are produced.
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

            # Genome size estimation from read stats
            if read_count >= 5:
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

                # Return longest sequence
                best_seq = max(sequences, key=lambda x: len(x.seq))
                raw_seq = str(best_seq.seq)

                # Check for internal N-runs in wtdbg2 output
                # If the assembly contains Ns even after clustering, the region is too complex
                # Reject rather than insert potentially incorrect sequence
                if 'N' in raw_seq.upper():
                    import re
                    n_count = len(re.findall(r'N+', raw_seq, re.IGNORECASE))
                    n_total = raw_seq.upper().count('N')
                    self.logger.info(f"    Assembly contains {n_count} N-runs ({n_total} Ns total), "
                                   f"region too complex, skipping")
                    return None

                return raw_seq

            return None

        except Exception as e:
            self.logger.warning(f"wtdbg2 error: {e}")
            return None

    def _finalize_result(self, result: Dict, gap: Dict, reads_info: Dict) -> Dict:
        """
        Finalize a successful fill result.

        Since validation already happened during the tier loop (validate-before-apply),
        this just adds reads info metadata.
        """
        if 'validation' not in result:
            result['validation'] = {'valid': True, 'status': 'validated'}

        result['reads_info'] = {
            'hifi_spanning': reads_info.get('hifi_spanning_count', 0),
            'ont_spanning': reads_info.get('ont_spanning_count', 0),
            'hifi_flanking': reads_info.get('hifi_left_count', 0) + reads_info.get('hifi_right_count', 0),
            'ont_flanking': reads_info.get('ont_left_count', 0) + reads_info.get('ont_right_count', 0),
        }

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
