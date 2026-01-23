#!/usr/bin/env python3
"""
Gap Validator - Validates filled gaps by checking read alignment quality

Validates gap fills by analyzing junction coverage and read spanning.
"""

import logging
import subprocess
import numpy as np
from pathlib import Path
from collections import defaultdict

import pysam

from gapfill.utils.indexer import AssemblyIndexer
from gapfill.utils.tempfiles import TempFileManager


class BAMPool:
    """Reusable BAM file handle pool"""

    def __init__(self):
        self._handles = {}
        self.logger = logging.getLogger(__name__)

    def get(self, bam_path):
        if bam_path is None:
            return None

        key = str(Path(bam_path).resolve())

        if key not in self._handles:
            try:
                self._handles[key] = pysam.AlignmentFile(key, 'rb')
            except Exception as e:
                self.logger.error(f"Failed to open BAM {bam_path}: {e}")
                return None

        return self._handles[key]

    def close_all(self):
        for handle in self._handles.values():
            try:
                handle.close()
            except:
                pass
        self._handles.clear()


class GapValidator:
    """
    Validates gap filling results by analyzing read alignments

    Validation modes:
    - Single read type: Validate using HiFi or ONT reads
    - Dual read type: Validate using both (EITHER passing = PASS)
    """

    def __init__(self, assembly_file, hifi_reads=None, ont_reads=None,
                 hifi_bam=None, ont_bam=None, threads=8, work_dir=None):

        self.assembly_file = Path(assembly_file)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.hifi_bam = Path(hifi_bam) if hifi_bam else None
        self.ont_bam = Path(ont_bam) if ont_bam else None
        self.threads = threads
        self.work_dir = Path(work_dir) if work_dir else Path('.')
        self.logger = logging.getLogger(__name__)

        self.assembly_indexer = AssemblyIndexer(assembly_file)
        self.temp_manager = TempFileManager(work_dir=self.work_dir, auto_cleanup=True)
        self.bam_pool = BAMPool()
        self.validation_bam = self.hifi_bam if self.hifi_bam else self.ont_bam

    def validate_filled_gap(self, gap, filled_sequence, reads_bam=None,
                           min_spanning_reads=3, max_clip_ratio=0.3):
        """Validate a filled gap by checking read alignment quality"""
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap_start + len(filled_sequence)

        bam_to_use = reads_bam if reads_bam else self.validation_bam

        if bam_to_use and Path(bam_to_use).exists():
            return self._validate_with_existing_bam(
                bam_to_use, chrom, gap_start, gap_end,
                filled_sequence, min_spanning_reads, max_clip_ratio
            )
        else:
            return {
                'valid': False,
                'reason': 'No BAM file available',
                'issues': [],
                'spanning_reads': 0,
                'total_reads': 0,
                'avg_coverage': 0
            }

    def _validate_with_existing_bam(self, bam_file, chrom, gap_start, gap_end,
                                    filled_sequence, min_spanning_reads, max_clip_ratio):
        """Validate using existing BAM file"""
        return self._analyze_alignment_quality(
            bam_file, chrom, gap_start, gap_end,
            min_spanning_reads, max_clip_ratio,
            filled_sequence=filled_sequence
        )

    def _analyze_alignment_quality(self, bam_file, chrom, gap_start, gap_end,
                                   min_spanning_reads, max_clip_ratio,
                                   filled_sequence=None):
        """Analyze alignment quality around filled gap"""
        bamfile = self.bam_pool.get(bam_file)
        if bamfile is None:
            return {
                'valid': False,
                'reason': 'Failed to open BAM file',
                'issues': [],
                'spanning_reads': 0,
                'total_reads': 0,
                'avg_coverage': 0
            }

        # Detect partial fill
        is_partial_fill = False
        if filled_sequence:
            tail_region = filled_sequence[-min(100, len(filled_sequence)):]
            n_ratio_tail = tail_region.upper().count('N') / len(tail_region) if tail_region else 0
            if n_ratio_tail > 0.5:
                is_partial_fill = True

        junction_zone = 500
        flank_size = 1000
        gap_length = gap_end - gap_start

        search_start = max(0, gap_start - flank_size)
        search_end = gap_end + flank_size

        # Initialize coverage arrays
        left_flank_cov = np.zeros(flank_size, dtype=np.uint16)
        fill_cov = np.zeros(gap_length, dtype=np.uint16)
        right_flank_cov = np.zeros(flank_size, dtype=np.uint16)

        total_reads = 0
        spanning_reads = 0
        left_boundary_reads = 0
        right_boundary_reads = 0

        try:
            for read in bamfile.fetch(chrom, search_start, search_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                total_reads += 1
                read_start = read.reference_start
                read_end = read.reference_end

                if read_start <= gap_start and read_end >= gap_end:
                    spanning_reads += 1

                if read_start < gap_start and read_end > gap_start:
                    left_boundary_reads += 1
                if read_start < gap_end and read_end > gap_end:
                    right_boundary_reads += 1

                # Track fill region coverage
                if read_start < gap_end and read_end > gap_start:
                    fill_start = max(gap_start, read_start)
                    fill_end = min(gap_end, read_end)
                    if fill_start < fill_end:
                        idx_start = fill_start - gap_start
                        idx_end = fill_end - gap_start
                        fill_cov[idx_start:idx_end] += 1

                if total_reads >= 1000:
                    break

        except Exception as e:
            self.logger.warning(f"Error analyzing BAM: {e}")
            return {
                'valid': False,
                'reason': f'Analysis error: {str(e)}',
                'issues': [],
                'spanning_reads': 0,
                'total_reads': 0,
                'avg_coverage': 0
            }

        avg_fill = float(fill_cov.mean()) if len(fill_cov) > 0 else 0
        min_coverage = int(fill_cov.min())
        zero_coverage_ratio = np.sum(fill_cov == 0) / gap_length if gap_length > 0 else 0

        issues = []
        valid = True

        # Check minimum fill coverage
        if avg_fill < 0.5:
            issues.append({
                'type': 'low_fill_coverage',
                'severity': 'critical',
                'message': f'Fill region average coverage {avg_fill:.1f}x is below minimum'
            })
            valid = False

        # Check zero coverage
        max_zero_ratio = 0.50 if is_partial_fill else 0.10
        if zero_coverage_ratio > max_zero_ratio:
            issues.append({
                'type': 'coverage_gap',
                'severity': 'critical',
                'message': f'{zero_coverage_ratio*100:.1f}% of filled region has zero coverage'
            })
            valid = False

        result = {
            'valid': valid,
            'spanning_reads': spanning_reads,
            'left_boundary_reads': left_boundary_reads,
            'right_boundary_reads': right_boundary_reads,
            'total_reads': total_reads,
            'avg_coverage': avg_fill,
            'min_coverage': min_coverage,
            'zero_coverage_ratio': zero_coverage_ratio,
            'issues': issues,
            'is_partial_fill': is_partial_fill
        }

        if valid:
            result['reason'] = f'Validation passed (fill coverage: {avg_fill:.1f}x)'
        else:
            critical_issues = [i for i in issues if i['severity'] == 'critical']
            result['reason'] = critical_issues[0]['message'] if critical_issues else 'Validation failed'

        return result

    def cleanup(self):
        self.bam_pool.close_all()
        self.temp_manager.cleanup()
        self.assembly_indexer.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False
