#!/usr/bin/env python3
"""
Gap Validator - Validates filled gaps and determines gap fillability

Key responsibilities:
1. Validate filled sequences (Complete/Partial)
2. Detect flank issues that may cause filling failures
3. Determine if a gap is truly unfillable
4. Guide flank polishing decisions
"""

import logging
import re
import subprocess
from enum import Enum
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union
from collections import defaultdict

import pysam
import numpy as np
from Bio import SeqIO


class GapStatus(Enum):
    """Gap filling status"""
    PENDING = "pending"                    # Not yet attempted
    FILLED_PENDING = "filled_pending"      # Filled, awaiting validation in next iteration
    FILLED_COMPLETE = "filled_complete"    # Completely filled, validated
    FILLED_PARTIAL = "filled_partial"      # Partially filled, validated
    UNFILLABLE = "unfillable"              # Confirmed unfillable (correct flanks, no spanning reads)
    FAILED = "failed"                      # Failed, may retry
    NEEDS_POLISH = "needs_polish"          # Flanks need polishing before retry


@dataclass
class ValidationResult:
    """Result of gap fill validation"""
    valid: bool
    status: GapStatus
    confidence: float = 0.0  # 0-1

    # Coverage metrics
    avg_coverage: float = 0.0
    min_coverage: float = 0.0
    zero_coverage_ratio: float = 0.0

    # Spanning evidence
    spanning_reads: int = 0
    left_boundary_reads: int = 0
    right_boundary_reads: int = 0

    # Junction quality
    left_clip_ratio: float = 0.0
    right_clip_ratio: float = 0.0

    # Flank analysis (for failed gaps)
    left_flank_needs_polish: bool = False
    right_flank_needs_polish: bool = False
    flank_issues: List[str] = field(default_factory=list)

    # Reason
    reason: str = ""


@dataclass
class PartialFillValidationResult:
    """Result of partial fill validation with independent left/right results"""
    # Overall result
    valid: bool  # True if at least one side is valid
    status: GapStatus

    # Left side validation
    left_valid: bool = True
    left_fill_length: int = 0
    left_avg_coverage: float = 0.0
    left_junction_coverage: float = 0.0
    left_zero_coverage_ratio: float = 0.0
    left_issues: List[str] = field(default_factory=list)

    # Right side validation
    right_valid: bool = True
    right_fill_length: int = 0
    right_avg_coverage: float = 0.0
    right_junction_coverage: float = 0.0
    right_zero_coverage_ratio: float = 0.0
    right_issues: List[str] = field(default_factory=list)

    # Combined metrics
    avg_coverage: float = 0.0
    zero_coverage_ratio: float = 0.0
    reason: str = ""


@dataclass
class FlankAnalysis:
    """Analysis of flanking sequence quality"""
    needs_polish: bool = False
    clip_positions: List[int] = field(default_factory=list)  # Positions with clip accumulation
    mismatch_density: float = 0.0
    consensus_differs: bool = False
    issues: List[str] = field(default_factory=list)


class GapValidator:
    """
    Validates gap filling results and determines fillability

    Validation modes:
    - Complete fill: Check spanning reads, coverage continuity, junction quality
    - Partial fill: Check filled portion coverage, boundary evidence
    - Failed fill: Analyze flanks to determine if polish needed or truly unfillable
    """

    # Thresholds for validation
    MIN_SPANNING_READS_SHORT = 1    # For gaps < 1kb
    MIN_SPANNING_READS_LONG = 2     # For gaps >= 1kb
    MIN_AVG_COVERAGE = 2.0
    MAX_ZERO_COVERAGE_RATIO = 0.05  # 5% for complete
    MAX_ZERO_COVERAGE_RATIO_PARTIAL = 0.20  # 20% for partial
    MAX_JUNCTION_CLIP_RATIO = 0.25  # 25% clips at junction = suspicious

    # Thresholds for flank analysis
    FLANK_ANALYSIS_SIZE = 1000  # bp to analyze on each side
    CLIP_ACCUMULATION_THRESHOLD = 5  # reads with clip at same position
    HIGH_MISMATCH_DENSITY = 0.05  # 5% mismatch = suspicious

    def __init__(self, threads: int = 8):
        self.threads = threads
        self.logger = logging.getLogger(__name__)
        self._bam_handles: Dict[str, pysam.AlignmentFile] = {}

    # Thresholds for junction and coverage validation
    MIN_JUNCTION_COVERAGE = 5    # Minimum coverage at junction points
    MIN_INSERT_COVERAGE = 5      # Minimum average coverage for inserted sequence
    MAX_COVERAGE_GAP_RATIO = 0.05  # Maximum ratio of zero-coverage positions (5%)

    def validate_complete_fill(
        self,
        bam_file: str,
        chrom: str,
        gap_start: int,
        gap_end: int,
        filled_sequence: str
    ) -> ValidationResult:
        """
        Validate a completely filled gap.

        Complete fills require spanning reads to confirm the fill bridges
        the original gap correctly.

        Criteria:
        1. Has spanning reads (reads that span the entire filled region)
        2. Coverage level is adequate
        """
        fill_length = gap_end - gap_start
        bam = self._get_bam(bam_file)
        if bam is None:
            return ValidationResult(
                valid=False,
                status=GapStatus.FAILED,
                reason="Cannot open BAM file"
            )

        # Collect alignment statistics
        stats = self._analyze_region(bam, chrom, gap_start, gap_end)

        issues = []

        # Complete fills: REQUIRE spanning reads
        min_spanning = (self.MIN_SPANNING_READS_SHORT
                       if fill_length < 1000
                       else self.MIN_SPANNING_READS_LONG)
        if stats['spanning_reads'] < min_spanning:
            issues.append(f"Insufficient spanning reads: {stats['spanning_reads']} < {min_spanning}")

        # Also check basic coverage
        if stats['avg_coverage'] < self.MIN_AVG_COVERAGE:
            issues.append(f"Low coverage: {stats['avg_coverage']:.1f}x < {self.MIN_AVG_COVERAGE}x")

        # Determine result
        valid = len(issues) == 0
        confidence = self._calculate_confidence(stats, fill_length, is_complete=True)

        return ValidationResult(
            valid=valid,
            status=GapStatus.FILLED_COMPLETE if valid else GapStatus.FAILED,
            confidence=confidence,
            avg_coverage=stats['avg_coverage'],
            min_coverage=stats['min_coverage'],
            zero_coverage_ratio=stats['zero_coverage_ratio'],
            spanning_reads=stats['spanning_reads'],
            left_boundary_reads=stats['left_boundary_reads'],
            right_boundary_reads=stats['right_boundary_reads'],
            left_clip_ratio=stats['left_clip_ratio'],
            right_clip_ratio=stats['right_clip_ratio'],
            reason="; ".join(issues) if issues else "Validation passed"
        )

    def validate_partial_fill(
        self,
        bam_file: str,
        chrom: str,
        original_gap_start: int,
        original_gap_end: int,
        filled_sequence: str,
        new_gap_start: int,
        new_gap_end: int
    ) -> PartialFillValidationResult:
        """
        Validate a partially filled gap with independent left/right validation.

        Each side (left and right) is validated independently:
        1. Junction coverage >= 5 (where fill meets original sequence)
        2. Average coverage of inserted sequence >= 5
        3. No coverage breakpoints in inserted sequence (zero_cov_ratio < 5%)

        Returns:
            PartialFillValidationResult with independent left/right results.
            - If both sides pass: keep entire fill
            - If one side passes: keep that side, revert the failed side
            - If both fail: revert entire fill
        """
        bam = self._get_bam(bam_file)
        if bam is None:
            return PartialFillValidationResult(
                valid=False,
                status=GapStatus.FAILED,
                left_valid=False,
                right_valid=False,
                reason="Cannot open BAM file"
            )

        # Initialize result
        result = PartialFillValidationResult(
            valid=False,
            status=GapStatus.FAILED,
            left_valid=True,
            right_valid=True
        )

        # Determine filled portions
        # new_gap_start/end are where the 500N placeholder is in the filled sequence
        left_fill_length = new_gap_start - original_gap_start  # Left inserted portion
        right_fill_length = original_gap_end - new_gap_end     # Right inserted portion

        result.left_fill_length = left_fill_length
        result.right_fill_length = right_fill_length

        # =====================================================================
        # Validate LEFT side independently
        # =====================================================================
        if left_fill_length > 0:
            left_issues = []

            # Left junction coverage
            left_junction_cov = self._get_junction_coverage(
                bam, chrom, original_gap_start, window=50
            )
            result.left_junction_coverage = left_junction_cov

            if left_junction_cov < self.MIN_JUNCTION_COVERAGE:
                left_issues.append(f"junction coverage {left_junction_cov:.1f}x < {self.MIN_JUNCTION_COVERAGE}")

            # Left inserted sequence coverage
            left_stats = self._analyze_region(bam, chrom, original_gap_start, new_gap_start)
            result.left_avg_coverage = left_stats['avg_coverage']
            result.left_zero_coverage_ratio = left_stats['zero_coverage_ratio']

            if left_stats['avg_coverage'] < self.MIN_INSERT_COVERAGE:
                left_issues.append(f"avg coverage {left_stats['avg_coverage']:.1f}x < {self.MIN_INSERT_COVERAGE}")

            if left_stats['zero_coverage_ratio'] > self.MAX_COVERAGE_GAP_RATIO:
                left_issues.append(f"coverage gaps {left_stats['zero_coverage_ratio']*100:.1f}%")

            result.left_valid = len(left_issues) == 0
            result.left_issues = left_issues
        else:
            # No left fill - mark as valid (nothing to validate)
            result.left_valid = True
            result.left_fill_length = 0

        # =====================================================================
        # Validate RIGHT side independently
        # =====================================================================
        if right_fill_length > 0:
            right_issues = []

            # Right junction coverage (at the end of the right fill)
            right_junction_cov = self._get_junction_coverage(
                bam, chrom, new_gap_end + right_fill_length, window=50
            )
            result.right_junction_coverage = right_junction_cov

            if right_junction_cov < self.MIN_JUNCTION_COVERAGE:
                right_issues.append(f"junction coverage {right_junction_cov:.1f}x < {self.MIN_JUNCTION_COVERAGE}")

            # Right inserted sequence coverage
            right_stats = self._analyze_region(bam, chrom, new_gap_end, new_gap_end + right_fill_length)
            result.right_avg_coverage = right_stats['avg_coverage']
            result.right_zero_coverage_ratio = right_stats['zero_coverage_ratio']

            if right_stats['avg_coverage'] < self.MIN_INSERT_COVERAGE:
                right_issues.append(f"avg coverage {right_stats['avg_coverage']:.1f}x < {self.MIN_INSERT_COVERAGE}")

            if right_stats['zero_coverage_ratio'] > self.MAX_COVERAGE_GAP_RATIO:
                right_issues.append(f"coverage gaps {right_stats['zero_coverage_ratio']*100:.1f}%")

            result.right_valid = len(right_issues) == 0
            result.right_issues = right_issues
        else:
            # No right fill - mark as valid (nothing to validate)
            result.right_valid = True
            result.right_fill_length = 0

        # =====================================================================
        # Determine overall result
        # =====================================================================
        # Valid if at least one side with actual fill is valid
        has_left = left_fill_length > 0
        has_right = right_fill_length > 0

        if has_left and has_right:
            # Both sides have fills
            result.valid = result.left_valid or result.right_valid
        elif has_left:
            result.valid = result.left_valid
        elif has_right:
            result.valid = result.right_valid
        else:
            # No fill on either side (shouldn't happen)
            result.valid = False

        # Set status
        if result.left_valid and result.right_valid:
            result.status = GapStatus.FILLED_PARTIAL
        elif result.valid:
            result.status = GapStatus.FILLED_PARTIAL  # Partial success
        else:
            result.status = GapStatus.FAILED

        # Combined metrics
        coverages = []
        zero_ratios = []
        if has_left:
            coverages.append(result.left_avg_coverage)
            zero_ratios.append(result.left_zero_coverage_ratio)
        if has_right:
            coverages.append(result.right_avg_coverage)
            zero_ratios.append(result.right_zero_coverage_ratio)

        result.avg_coverage = sum(coverages) / len(coverages) if coverages else 0
        result.zero_coverage_ratio = sum(zero_ratios) / len(zero_ratios) if zero_ratios else 0

        # Build reason string
        reasons = []
        if has_left:
            if result.left_valid:
                reasons.append(f"Left OK ({left_fill_length}bp, cov={result.left_avg_coverage:.1f}x)")
            else:
                reasons.append(f"Left FAIL ({left_fill_length}bp): {'; '.join(result.left_issues)}")
        if has_right:
            if result.right_valid:
                reasons.append(f"Right OK ({right_fill_length}bp, cov={result.right_avg_coverage:.1f}x)")
            else:
                reasons.append(f"Right FAIL ({right_fill_length}bp): {'; '.join(result.right_issues)}")

        result.reason = " | ".join(reasons)

        return result

    def validate_fill_locally(self, hifi_bam, ont_bam, assembly_file, chrom,
                              gap_start, gap_end, fill_sequence, is_complete,
                              work_dir, threads=1, flank_size=5000):
        """
        Validate a fill by constructing a local reference and aligning nearby reads.

        Instead of relying on the next iteration's BAM (delayed validation), this
        builds a small reference (flank + fill + flank), extracts nearby reads from
        the current BAM(s), aligns them locally, and validates.

        Args:
            hifi_bam: Path to HiFi BAM (or None)
            ont_bam: Path to ONT BAM (or None)
            assembly_file: Path to current assembly FASTA
            chrom: Chromosome name
            gap_start: Gap start position in assembly
            gap_end: Gap end position in assembly
            fill_sequence: The fill sequence to validate
            is_complete: True if complete fill (no N placeholder)
            work_dir: Working directory for temp files
            threads: Threads for minimap2
            flank_size: Flank size to extract (default 5000)

        Returns:
            ValidationResult or PartialFillValidationResult, or None on error
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Extract flanks from assembly
            assembly_seqs = {}
            for record in SeqIO.parse(assembly_file, 'fasta'):
                if record.id == chrom:
                    assembly_seqs[record.id] = str(record.seq)
                    break

            if chrom not in assembly_seqs:
                self.logger.warning(f"Chromosome {chrom} not found in assembly")
                return None

            seq = assembly_seqs[chrom]
            left_start = max(0, gap_start - flank_size)
            left_flank = seq[left_start:gap_start]
            right_flank = seq[gap_end:gap_end + flank_size]

            # Build local reference: left_flank + fill_sequence + right_flank
            local_ref_seq = left_flank + fill_sequence + right_flank
            local_ref_file = work_dir / "local_ref.fa"
            with open(local_ref_file, 'w') as f:
                f.write(f">local_ref\n{local_ref_seq}\n")

            # Extract reads from BAM(s) in region around the gap
            local_reads_file = work_dir / "local_reads.fa"
            read_count = self._extract_reads_from_bams(
                hifi_bam, ont_bam, chrom,
                max(0, gap_start - flank_size),
                gap_end + flank_size,
                local_reads_file
            )

            if read_count < 1:
                self.logger.debug(f"  No reads extracted for local validation")
                return None

            # Align reads to local reference
            local_bam_file = work_dir / "local.bam"
            align_cmd = (
                f"minimap2 -ax map-hifi -t {threads} "
                f"{local_ref_file} {local_reads_file} | "
                f"samtools sort -@ {threads} -o {local_bam_file} - && "
                f"samtools index {local_bam_file}"
            )
            result = subprocess.run(
                align_cmd, shell=True, capture_output=True, text=True, timeout=120
            )
            if result.returncode != 0:
                self.logger.debug(f"  Local alignment failed: {result.stderr[:200]}")
                return None

            # Compute local coordinates
            fill_start_local = len(left_flank)
            fill_end_local = len(left_flank) + len(fill_sequence)

            local_bam_str = str(local_bam_file)

            if is_complete:
                # Complete fill validation
                validation = self.validate_complete_fill(
                    local_bam_str, "local_ref",
                    fill_start_local, fill_end_local,
                    fill_sequence
                )
                return validation
            else:
                # Partial fill — find N-run position in fill_sequence
                n_match = re.search(r'N{10,}', fill_sequence)
                if n_match:
                    # N-run positions in local ref coordinates
                    new_gap_start_local = fill_start_local + n_match.start()
                    new_gap_end_local = fill_start_local + n_match.end()

                    validation = self.validate_partial_fill(
                        local_bam_str, "local_ref",
                        fill_start_local, fill_end_local,
                        fill_sequence,
                        new_gap_start_local, new_gap_end_local
                    )
                    return validation
                else:
                    # No N-run but marked partial — validate as complete
                    validation = self.validate_complete_fill(
                        local_bam_str, "local_ref",
                        fill_start_local, fill_end_local,
                        fill_sequence
                    )
                    return validation

        except Exception as e:
            self.logger.warning(f"Local validation error: {e}")
            return None
        finally:
            # Close any BAM handles opened for local validation
            local_bam_str = str(work_dir / "local.bam")
            if local_bam_str in self._bam_handles:
                try:
                    self._bam_handles[local_bam_str].close()
                except:
                    pass
                del self._bam_handles[local_bam_str]

    def _extract_reads_from_bams(self, hifi_bam, ont_bam, chrom,
                                  region_start, region_end,
                                  output_file) -> int:
        """
        Extract read sequences from BAM file(s) in a region.

        Args:
            hifi_bam: Path to HiFi BAM (or None)
            ont_bam: Path to ONT BAM (or None)
            chrom: Chromosome name
            region_start: Region start
            region_end: Region end
            output_file: Output FASTA file

        Returns:
            Number of reads extracted
        """
        read_count = 0
        seen_names = set()

        with open(output_file, 'w') as f:
            for bam_path in [hifi_bam, ont_bam]:
                if not bam_path:
                    continue
                bam_path = str(bam_path)

                try:
                    with pysam.AlignmentFile(bam_path, 'rb') as bam:
                        if chrom not in bam.references:
                            continue

                        for read in bam.fetch(chrom, region_start, region_end):
                            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                                continue
                            if read.query_name in seen_names:
                                continue

                            seq = read.query_sequence
                            if seq and len(seq) >= 100:
                                f.write(f">{read.query_name}\n{seq}\n")
                                read_count += 1
                                seen_names.add(read.query_name)
                except Exception as e:
                    self.logger.debug(f"  Error extracting reads from {bam_path}: {e}")

        return read_count

    def _get_junction_coverage(self, bam: pysam.AlignmentFile, chrom: str,
                                position: int, window: int = 50) -> float:
        """
        Get average coverage at a junction point.

        Args:
            bam: BAM file handle
            chrom: Chromosome name
            position: Junction position
            window: Window size around junction (default 50bp each side)

        Returns:
            Average coverage at the junction
        """
        start = max(0, position - window)
        end = position + window

        try:
            coverage = np.zeros(end - start, dtype=np.uint16)

            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                read_start = read.reference_start
                read_end = read.reference_end

                if read_start is None or read_end is None:
                    continue

                cov_start = max(0, read_start - start)
                cov_end = min(len(coverage), read_end - start)
                if cov_start < cov_end:
                    coverage[cov_start:cov_end] += 1

            return float(np.mean(coverage))

        except Exception as e:
            self.logger.warning(f"Error getting junction coverage at {chrom}:{position}: {e}")
            return 0.0

    def analyze_failed_gap(
        self,
        bam_file: str,
        chrom: str,
        gap_start: int,
        gap_end: int
    ) -> ValidationResult:
        """
        Analyze a gap that failed to fill

        Determines:
        1. Are the flanking sequences correct? (no polish needed)
        2. Do flanks need polishing? (retry after polish)
        3. Is this gap truly unfillable? (skip in future)
        """
        bam = self._get_bam(bam_file)
        if bam is None:
            return ValidationResult(
                valid=False,
                status=GapStatus.FAILED,
                reason="Cannot open BAM file"
            )

        # Analyze left flank
        left_flank = self._analyze_flank(
            bam, chrom,
            max(0, gap_start - self.FLANK_ANALYSIS_SIZE),
            gap_start,
            is_left=True
        )

        # Analyze right flank
        right_flank = self._analyze_flank(
            bam, chrom,
            gap_end,
            gap_end + self.FLANK_ANALYSIS_SIZE,
            is_left=False
        )

        # Combine analysis
        needs_polish = left_flank.needs_polish or right_flank.needs_polish
        flank_issues = left_flank.issues + right_flank.issues

        if needs_polish:
            status = GapStatus.NEEDS_POLISH
            reason = f"Flanks need polishing: {'; '.join(flank_issues)}"
        else:
            # Flanks are correct but still can't fill → truly unfillable
            status = GapStatus.UNFILLABLE
            reason = "Flanks verified correct, no spanning reads available - gap is unfillable"

        return ValidationResult(
            valid=False,
            status=status,
            left_flank_needs_polish=left_flank.needs_polish,
            right_flank_needs_polish=right_flank.needs_polish,
            flank_issues=flank_issues,
            reason=reason
        )

    def pre_assess_gap(
        self,
        bam_file: str,
        chrom: str,
        gap_start: int,
        gap_end: int
    ) -> ValidationResult:
        """
        Pre-assess a gap BEFORE attempting to fill it.

        This should be called before the first filling iteration to identify
        gaps that need flank polishing or reassembly before filling can succeed.

        Returns:
        - PENDING: Flanks look good, can attempt filling
        - NEEDS_POLISH: Flanks have issues (clip accumulation, high mismatches)
        - FAILED: Cannot assess (BAM unavailable)
        """
        bam = self._get_bam(bam_file)
        if bam is None:
            return ValidationResult(
                valid=True,
                status=GapStatus.PENDING,
                reason="Cannot open BAM file for pre-assessment, will attempt filling"
            )

        # Analyze left flank
        left_flank = self._analyze_flank(
            bam, chrom,
            max(0, gap_start - self.FLANK_ANALYSIS_SIZE),
            gap_start,
            is_left=True
        )

        # Analyze right flank
        right_flank = self._analyze_flank(
            bam, chrom,
            gap_end,
            gap_end + self.FLANK_ANALYSIS_SIZE,
            is_left=False
        )

        # Combine analysis
        needs_polish = left_flank.needs_polish or right_flank.needs_polish
        flank_issues = left_flank.issues + right_flank.issues

        if needs_polish:
            return ValidationResult(
                valid=False,
                status=GapStatus.NEEDS_POLISH,
                left_flank_needs_polish=left_flank.needs_polish,
                right_flank_needs_polish=right_flank.needs_polish,
                flank_issues=flank_issues,
                reason=f"Pre-assessment: flanks need polishing - {'; '.join(flank_issues)}"
            )
        else:
            return ValidationResult(
                valid=True,
                status=GapStatus.PENDING,
                left_flank_needs_polish=False,
                right_flank_needs_polish=False,
                flank_issues=[],
                reason="Pre-assessment: flanks look good, ready for filling"
            )

    def pre_assess_gaps(
        self,
        bam_file: str,
        gaps: List[Dict]
    ) -> Dict[str, ValidationResult]:
        """
        Pre-assess multiple gaps before filling.

        Args:
            bam_file: BAM file for analysis
            gaps: List of gap dictionaries with 'chrom', 'start', 'end', 'name'

        Returns:
            Dict mapping gap_name to ValidationResult
        """
        results = {}
        needs_polish_count = 0
        ready_count = 0

        for gap in gaps:
            gap_name = gap.get('name', f"{gap['chrom']}_{gap['start']}_{gap['end']}")
            result = self.pre_assess_gap(
                bam_file, gap['chrom'], gap['start'], gap['end']
            )
            results[gap_name] = result

            if result.status == GapStatus.NEEDS_POLISH:
                needs_polish_count += 1
            else:
                ready_count += 1

        self.logger.info(f"Pre-assessment complete: {ready_count} ready, {needs_polish_count} need polish")
        return results

    def _analyze_region(
        self,
        bam: pysam.AlignmentFile,
        chrom: str,
        start: int,
        end: int
    ) -> Dict:
        """Analyze read alignments in a region"""
        region_length = end - start
        if region_length <= 0:
            return self._empty_stats()

        coverage = np.zeros(region_length, dtype=np.uint16)

        spanning_reads = 0
        left_boundary_reads = 0
        right_boundary_reads = 0
        left_clips = 0
        right_clips = 0
        total_left_reads = 0
        total_right_reads = 0

        junction_zone = 100  # bp around junction to check for clips

        try:
            for read in bam.fetch(chrom, max(0, start - 500), end + 500):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                read_start = read.reference_start
                read_end = read.reference_end

                if read_start is None or read_end is None:
                    continue

                # Spanning check
                if read_start <= start and read_end >= end:
                    spanning_reads += 1

                # Boundary crossing
                if read_start < start < read_end:
                    left_boundary_reads += 1
                if read_start < end < read_end:
                    right_boundary_reads += 1

                # Coverage within region
                cov_start = max(0, read_start - start)
                cov_end = min(region_length, read_end - start)
                if cov_start < cov_end:
                    coverage[cov_start:cov_end] += 1

                # Junction clip analysis
                cigar = read.cigartuples
                if cigar:
                    # Left junction (reads ending near start with right clip)
                    if start - junction_zone <= read_end <= start + junction_zone:
                        total_left_reads += 1
                        if cigar[-1][0] in (4, 5):  # S or H at end
                            left_clips += 1

                    # Right junction (reads starting near end with left clip)
                    if end - junction_zone <= read_start <= end + junction_zone:
                        total_right_reads += 1
                        if cigar[0][0] in (4, 5):  # S or H at start
                            right_clips += 1

        except Exception as e:
            self.logger.warning(f"Error analyzing region {chrom}:{start}-{end}: {e}")
            return self._empty_stats()

        return {
            'spanning_reads': spanning_reads,
            'left_boundary_reads': left_boundary_reads,
            'right_boundary_reads': right_boundary_reads,
            'avg_coverage': float(np.mean(coverage)),
            'min_coverage': float(np.min(coverage)),
            'zero_coverage_ratio': float(np.sum(coverage == 0) / region_length),
            'left_clip_ratio': left_clips / total_left_reads if total_left_reads > 0 else 0,
            'right_clip_ratio': right_clips / total_right_reads if total_right_reads > 0 else 0,
        }

    def _analyze_flank(
        self,
        bam: pysam.AlignmentFile,
        chrom: str,
        start: int,
        end: int,
        is_left: bool
    ) -> FlankAnalysis:
        """
        Analyze a flanking region to determine if it needs polishing

        Signs that flank needs polish:
        1. Clip accumulation: many reads clip at the same position
        2. High mismatch density: reads disagree with reference
        3. Coverage anomaly: sudden drop not at the gap boundary
        """
        result = FlankAnalysis()

        if end <= start:
            return result

        region_length = end - start

        # Track clip positions
        clip_positions = defaultdict(int)
        mismatch_count = 0
        total_aligned_bases = 0

        try:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                read_start = read.reference_start
                read_end = read.reference_end
                cigar = read.cigartuples

                if not cigar or read_start is None:
                    continue

                # Track soft clips
                if cigar[0][0] == 4:  # Left soft clip
                    clip_pos = read_start
                    if start <= clip_pos < end:
                        clip_positions[clip_pos] += 1

                if cigar[-1][0] == 4:  # Right soft clip
                    clip_pos = read_end
                    if start <= clip_pos < end:
                        clip_positions[clip_pos] += 1

                # Count mismatches using MD tag if available
                if read.has_tag('MD'):
                    md = read.get_tag('MD')
                    # Count non-digit characters (mismatches/deletions)
                    mismatches = len(re.findall(r'[ACGT]', str(md)))
                    mismatch_count += mismatches
                    total_aligned_bases += read.query_alignment_length

        except Exception as e:
            self.logger.warning(f"Error analyzing flank {chrom}:{start}-{end}: {e}")
            return result

        # Check for clip accumulation
        for pos, count in clip_positions.items():
            if count >= self.CLIP_ACCUMULATION_THRESHOLD:
                result.clip_positions.append(pos)
                result.needs_polish = True
                result.issues.append(
                    f"Clip accumulation at {chrom}:{pos} ({count} reads)"
                )

        # Check mismatch density
        if total_aligned_bases > 0:
            result.mismatch_density = mismatch_count / total_aligned_bases
            if result.mismatch_density > self.HIGH_MISMATCH_DENSITY:
                result.needs_polish = True
                result.issues.append(
                    f"High mismatch density: {result.mismatch_density*100:.1f}%"
                )

        return result

    def _calculate_confidence(
        self,
        stats: Dict,
        gap_length: int,
        is_complete: bool
    ) -> float:
        """Calculate confidence score (0-1) based on validation metrics"""
        score = 0.0

        # Spanning reads contribution (0-0.4)
        if is_complete:
            if stats['spanning_reads'] >= 5:
                score += 0.4
            elif stats['spanning_reads'] >= 2:
                score += 0.3
            elif stats['spanning_reads'] >= 1:
                score += 0.2

        # Coverage contribution (0-0.3)
        if stats['avg_coverage'] >= 10:
            score += 0.3
        elif stats['avg_coverage'] >= 5:
            score += 0.2
        elif stats['avg_coverage'] >= 2:
            score += 0.1

        # Coverage continuity (0-0.2)
        if stats['zero_coverage_ratio'] == 0:
            score += 0.2
        elif stats['zero_coverage_ratio'] < 0.02:
            score += 0.1

        # Junction quality (0-0.1)
        avg_clip = (stats['left_clip_ratio'] + stats['right_clip_ratio']) / 2
        if avg_clip < 0.1:
            score += 0.1
        elif avg_clip < 0.2:
            score += 0.05

        return min(1.0, score)

    def _empty_stats(self) -> Dict:
        """Return empty statistics dict"""
        return {
            'spanning_reads': 0,
            'left_boundary_reads': 0,
            'right_boundary_reads': 0,
            'avg_coverage': 0,
            'min_coverage': 0,
            'zero_coverage_ratio': 1.0,
            'left_clip_ratio': 0,
            'right_clip_ratio': 0,
        }

    def _get_bam(self, bam_file: str) -> Optional[pysam.AlignmentFile]:
        """Get or open a BAM file handle"""
        if bam_file not in self._bam_handles:
            try:
                self._bam_handles[bam_file] = pysam.AlignmentFile(bam_file, 'rb')
            except Exception as e:
                self.logger.error(f"Cannot open BAM {bam_file}: {e}")
                return None
        return self._bam_handles[bam_file]

    def close(self):
        """Close all BAM handles"""
        for handle in self._bam_handles.values():
            try:
                handle.close()
            except:
                pass
        self._bam_handles.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


@dataclass
class PendingFill:
    """Information about a fill awaiting validation.

    DEPRECATED: Retained for backwards compatibility with old checkpoint files.
    The validate-before-apply architecture no longer uses pending fills.
    """
    gap_id: str
    chrom: str
    original_start: int
    original_end: int
    filled_start: int  # New coordinates after filling
    filled_end: int
    sequence: str
    is_complete: bool  # True if no placeholder N's
    source: str  # e.g., "hifi_spanning", "ont_flanking_merged"
    tier: int


class GapStatusTracker:
    """
    Tracks gap status across iterations.

    With the validate-before-apply architecture, fills are validated locally
    before being applied. No more pending fills or delayed validation.
    """

    def __init__(self):
        self._status: Dict[str, GapStatus] = {}
        self._history: Dict[str, List[str]] = defaultdict(list)
        self.logger = logging.getLogger(__name__)

    def get_status(self, gap_id: str) -> GapStatus:
        """Get current status of a gap"""
        return self._status.get(gap_id, GapStatus.PENDING)

    def set_status(self, gap_id: str, status: GapStatus, reason: str = ""):
        """Set gap status with optional reason"""
        self._status[gap_id] = status
        self._history[gap_id].append(f"{status.value}: {reason}")

        if status == GapStatus.UNFILLABLE:
            self.logger.info(f"Gap {gap_id} marked as UNFILLABLE - will skip in future iterations")

    def should_attempt(self, gap_id: str) -> bool:
        """Check if gap should be attempted in current iteration"""
        status = self.get_status(gap_id)

        # Skip: already filled (complete or partial), unfillable
        if status in (GapStatus.UNFILLABLE, GapStatus.FILLED_COMPLETE,
                      GapStatus.FILLED_PARTIAL):
            return False

        return True

    def get_summary(self) -> Dict[str, int]:
        """Get summary counts by status"""
        summary = defaultdict(int)
        for status in self._status.values():
            summary[status.value] += 1
        return dict(summary)

    def get_gaps_by_status(self, status: GapStatus) -> List[str]:
        """Get list of gap IDs with given status"""
        return [gid for gid, s in self._status.items() if s == status]

    def get_history(self, gap_id: str) -> List[str]:
        """Get history of status changes for a gap"""
        return self._history.get(gap_id, [])

    def to_dict(self) -> Dict:
        """Serialize to dictionary"""
        return {
            'status': {k: v.value for k, v in self._status.items()},
            'history': dict(self._history),
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'GapStatusTracker':
        """Deserialize from dictionary, handling old format gracefully"""
        tracker = cls()
        for gap_id, status_str in data.get('status', {}).items():
            try:
                status = GapStatus(status_str)
            except ValueError:
                status = GapStatus.FAILED
            # Convert old FILLED_PENDING status to FAILED (needs re-fill)
            if status == GapStatus.FILLED_PENDING:
                tracker.logger.info(f"Converting {gap_id} from FILLED_PENDING to FAILED (old format)")
                status = GapStatus.FAILED
            tracker._status[gap_id] = status
        for gap_id, history in data.get('history', {}).items():
            tracker._history[gap_id] = history
        # Gracefully ignore old pending_fills data
        if 'pending_fills' in data and data['pending_fills']:
            tracker.logger.info(f"Ignoring {len(data['pending_fills'])} old pending_fills entries")
        return tracker
