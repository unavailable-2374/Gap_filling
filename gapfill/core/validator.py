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
from enum import Enum
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import pysam
import numpy as np


class GapStatus(Enum):
    """Gap filling status"""
    PENDING = "pending"                    # Not yet attempted
    FILLED_COMPLETE = "filled_complete"    # Completely filled, validated
    FILLED_PARTIAL = "filled_partial"      # Partially filled, validated
    UNFILLABLE = "unfillable"              # Confirmed unfillable (correct flanks, no spanning reads)
    FAILED = "failed"                      # Failed, may retry after flank polish
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

    def validate_complete_fill(
        self,
        bam_file: str,
        chrom: str,
        gap_start: int,
        gap_end: int,
        filled_sequence: str
    ) -> ValidationResult:
        """
        Validate a completely filled gap

        Criteria:
        1. Has spanning reads
        2. Coverage is continuous (no zero-coverage regions)
        3. Junctions are clean (low clip ratio)
        4. Coverage level is adequate
        """
        gap_length = gap_end - gap_start
        bam = self._get_bam(bam_file)
        if bam is None:
            return ValidationResult(
                valid=False,
                status=GapStatus.FAILED,
                reason="Cannot open BAM file"
            )

        # Collect alignment statistics
        stats = self._analyze_region(bam, chrom, gap_start, gap_end)

        # Determine thresholds based on gap size
        min_spanning = (self.MIN_SPANNING_READS_SHORT
                       if gap_length < 1000
                       else self.MIN_SPANNING_READS_LONG)

        # Validation checks
        issues = []

        # Check 1: Spanning reads
        if stats['spanning_reads'] < min_spanning:
            issues.append(f"Insufficient spanning reads: {stats['spanning_reads']} < {min_spanning}")

        # Check 2: Coverage continuity
        if stats['zero_coverage_ratio'] > self.MAX_ZERO_COVERAGE_RATIO:
            issues.append(f"Coverage gaps: {stats['zero_coverage_ratio']*100:.1f}% zero coverage")

        # Check 3: Average coverage
        if stats['avg_coverage'] < self.MIN_AVG_COVERAGE:
            issues.append(f"Low coverage: {stats['avg_coverage']:.1f}x < {self.MIN_AVG_COVERAGE}x")

        # Check 4: Junction quality
        if stats['left_clip_ratio'] > self.MAX_JUNCTION_CLIP_RATIO:
            issues.append(f"Left junction clips: {stats['left_clip_ratio']*100:.1f}%")
        if stats['right_clip_ratio'] > self.MAX_JUNCTION_CLIP_RATIO:
            issues.append(f"Right junction clips: {stats['right_clip_ratio']*100:.1f}%")

        # Determine result
        valid = len(issues) == 0
        confidence = self._calculate_confidence(stats, gap_length, is_complete=True)

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
    ) -> ValidationResult:
        """
        Validate a partially filled gap

        Criteria:
        1. Filled portion has coverage
        2. New boundary has read evidence (reads terminating there)
        3. Filled portion coverage is reasonable
        """
        bam = self._get_bam(bam_file)
        if bam is None:
            return ValidationResult(
                valid=False,
                status=GapStatus.FAILED,
                reason="Cannot open BAM file"
            )

        # Determine which side was filled
        left_filled = new_gap_start > original_gap_start
        right_filled = new_gap_end < original_gap_end

        issues = []

        # Analyze filled portions
        if left_filled:
            left_stats = self._analyze_region(
                bam, chrom, original_gap_start, new_gap_start
            )
            if left_stats['avg_coverage'] < self.MIN_AVG_COVERAGE:
                issues.append(f"Left fill low coverage: {left_stats['avg_coverage']:.1f}x")
            if left_stats['zero_coverage_ratio'] > self.MAX_ZERO_COVERAGE_RATIO_PARTIAL:
                issues.append(f"Left fill has coverage gaps")
            # Check for reads terminating at new boundary (expected for partial)
            if left_stats['right_boundary_reads'] < 1:
                issues.append("No reads at new left boundary")

        if right_filled:
            right_stats = self._analyze_region(
                bam, chrom, new_gap_end, original_gap_end
            )
            if right_stats['avg_coverage'] < self.MIN_AVG_COVERAGE:
                issues.append(f"Right fill low coverage: {right_stats['avg_coverage']:.1f}x")
            if right_stats['zero_coverage_ratio'] > self.MAX_ZERO_COVERAGE_RATIO_PARTIAL:
                issues.append("Right fill has coverage gaps")
            if right_stats['left_boundary_reads'] < 1:
                issues.append("No reads at new right boundary")

        valid = len(issues) == 0

        # Combined stats
        avg_cov = 0
        if left_filled and right_filled:
            avg_cov = (left_stats['avg_coverage'] + right_stats['avg_coverage']) / 2
        elif left_filled:
            avg_cov = left_stats['avg_coverage']
        elif right_filled:
            avg_cov = right_stats['avg_coverage']

        return ValidationResult(
            valid=valid,
            status=GapStatus.FILLED_PARTIAL if valid else GapStatus.FAILED,
            confidence=0.7 if valid else 0.3,  # Partial fills have lower base confidence
            avg_coverage=avg_cov,
            reason="; ".join(issues) if issues else "Partial fill validated"
        )

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


class GapStatusTracker:
    """
    Tracks gap status across iterations

    Ensures unfillable gaps are skipped in future iterations
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

        # Skip unfillable and already filled gaps
        if status in (GapStatus.UNFILLABLE, GapStatus.FILLED_COMPLETE):
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
            'history': dict(self._history)
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'GapStatusTracker':
        """Deserialize from dictionary"""
        tracker = cls()
        for gap_id, status_str in data.get('status', {}).items():
            tracker._status[gap_id] = GapStatus(status_str)
        for gap_id, history in data.get('history', {}).items():
            tracker._history[gap_id] = history
        return tracker
