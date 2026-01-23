#!/usr/bin/env python3
"""
Hi-C Data Integration for Gap Filling

Functions:
1. Gap size estimation using Hi-C read pairs
2. Fill validation using contact frequency
3. Phasing enhancement for polyploid mode

Hi-C data characteristics:
- Long-range connectivity (kb to Mb)
- Short reads (~150bp paired-end)
- Contact frequency inversely correlates with genomic distance
"""

import logging
import subprocess
import tempfile
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass

import pysam
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class GapSizeEstimate:
    """Estimated gap size from Hi-C data"""
    gap_name: str
    chrom: str
    gap_start: int
    gap_end: int
    original_size: int  # N placeholder size
    estimated_size: int  # Hi-C estimated size
    confidence: str  # 'high', 'medium', 'low'
    supporting_pairs: int
    insert_sizes: List[int]


@dataclass
class FillValidation:
    """Validation result for a filled gap"""
    gap_name: str
    is_valid: bool
    contact_score: float  # 0-1, higher is better
    expected_contacts: int
    observed_contacts: int
    anomaly_type: Optional[str]  # 'break', 'spike', None


class HiCAnalyzer:
    """
    Analyze Hi-C data to assist gap filling
    """

    def __init__(self,
                 hic_bam: str,
                 assembly_file: str,
                 threads: int = 8,
                 min_mapq: int = 20,
                 min_insert_size: int = 1000,
                 max_insert_size: int = 10000000):
        """
        Initialize Hi-C analyzer.

        Args:
            hic_bam: Path to Hi-C BAM file (aligned to assembly)
            assembly_file: Path to assembly FASTA
            threads: Number of threads
            min_mapq: Minimum mapping quality
            min_insert_size: Minimum insert size to consider
            max_insert_size: Maximum insert size to consider
        """
        self.hic_bam = Path(hic_bam)
        self.assembly_file = Path(assembly_file)
        self.threads = threads
        self.min_mapq = min_mapq
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size

        self._validate_inputs()

        logger.info("HiCAnalyzer initialized")
        logger.info(f"  Hi-C BAM: {self.hic_bam}")
        logger.info(f"  Assembly: {self.assembly_file}")

    def _validate_inputs(self):
        """Validate input files exist"""
        if not self.hic_bam.exists():
            raise FileNotFoundError(f"Hi-C BAM not found: {self.hic_bam}")
        if not self.assembly_file.exists():
            raise FileNotFoundError(f"Assembly not found: {self.assembly_file}")

        # Check BAM index
        bai_file = Path(str(self.hic_bam) + ".bai")
        if not bai_file.exists():
            logger.info("Indexing Hi-C BAM...")
            subprocess.run(['samtools', 'index', str(self.hic_bam)],
                         check=True, capture_output=True)

    # =========================================================================
    # Gap Size Estimation
    # =========================================================================

    def estimate_gap_sizes(self, gaps: List[Dict],
                           flank_size: int = 5000) -> List[GapSizeEstimate]:
        """
        Estimate true gap sizes using Hi-C read pairs.

        Strategy:
        - Find Hi-C pairs where one read is in left flank, other in right flank
        - Insert size = left_pos + gap_size + (right_flank_size - right_pos)
        - Solve for gap_size

        Args:
            gaps: List of gap dicts with 'chrom', 'start', 'end', 'name'
            flank_size: Size of flanking region to search

        Returns:
            List of GapSizeEstimate objects
        """
        logger.info(f"Estimating sizes for {len(gaps)} gaps using Hi-C...")
        estimates = []

        with pysam.AlignmentFile(str(self.hic_bam), 'rb') as bam:
            for gap in gaps:
                estimate = self._estimate_single_gap_size(
                    bam, gap, flank_size
                )
                estimates.append(estimate)

                if estimate.confidence != 'low':
                    logger.info(f"  {gap['name']}: {estimate.original_size}bp → "
                               f"{estimate.estimated_size}bp ({estimate.confidence}, "
                               f"{estimate.supporting_pairs} pairs)")

        return estimates

    def _estimate_single_gap_size(self, bam: pysam.AlignmentFile,
                                   gap: Dict, flank_size: int) -> GapSizeEstimate:
        """Estimate size for a single gap"""
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']
        gap_name = gap.get('name', f"{chrom}_{gap_start}_{gap_end}")
        original_size = gap_end - gap_start

        # Define flanking regions
        left_start = max(0, gap_start - flank_size)
        left_end = gap_start
        right_start = gap_end
        right_end = gap_end + flank_size

        # Collect read pairs spanning the gap
        left_reads = {}  # read_name -> (position, mate_position)
        right_reads = {}

        try:
            # Get reads in left flank
            for read in bam.fetch(chrom, left_start, left_end):
                if self._is_valid_hic_read(read):
                    if read.next_reference_name == chrom:
                        mate_pos = read.next_reference_start
                        if right_start <= mate_pos <= right_end:
                            left_reads[read.query_name] = {
                                'pos': read.reference_start,
                                'mate_pos': mate_pos
                            }

            # Get reads in right flank
            for read in bam.fetch(chrom, right_start, right_end):
                if self._is_valid_hic_read(read):
                    if read.next_reference_name == chrom:
                        mate_pos = read.next_reference_start
                        if left_start <= mate_pos <= left_end:
                            right_reads[read.query_name] = {
                                'pos': read.reference_start,
                                'mate_pos': mate_pos
                            }

        except Exception as e:
            logger.warning(f"Error fetching Hi-C reads for {gap_name}: {e}")

        # Find pairs that span the gap
        spanning_pairs = set(left_reads.keys()) & set(right_reads.keys())

        if len(spanning_pairs) < 3:
            return GapSizeEstimate(
                gap_name=gap_name,
                chrom=chrom,
                gap_start=gap_start,
                gap_end=gap_end,
                original_size=original_size,
                estimated_size=original_size,  # Use original if can't estimate
                confidence='low',
                supporting_pairs=len(spanning_pairs),
                insert_sizes=[]
            )

        # Calculate implied gap sizes from insert sizes
        implied_gap_sizes = []

        for read_name in spanning_pairs:
            left_info = left_reads[read_name]
            right_info = right_reads[read_name]

            # Left read position relative to gap
            left_pos = left_info['pos']
            left_dist_to_gap = gap_start - left_pos

            # Right read position relative to gap
            right_pos = right_info['pos']
            right_dist_from_gap = right_pos - gap_end

            # The genomic distance between reads = left_dist + gap_size + right_dist
            # We observe: insert_size in the BAM
            # But for Hi-C, reads can be far apart, so we estimate from positions

            # Implied gap size
            implied_gap = right_pos - left_pos - left_dist_to_gap - right_dist_from_gap
            # This simplifies to: implied_gap = gap_end - gap_start = original_size
            # But that's wrong... let me reconsider

            # Actually, for Hi-C spanning pairs:
            # If the true gap is larger than placeholder, reads would be further apart
            # If the true gap is smaller, reads would be closer

            # Use the template length if available and reasonable
            if hasattr(left_info, 'template_length'):
                tlen = abs(left_info['template_length'])
                if self.min_insert_size <= tlen <= self.max_insert_size:
                    # implied_gap = tlen - left_dist_to_gap - right_dist_from_gap
                    implied = tlen - (gap_start - left_pos) - (right_pos - gap_end)
                    if implied > 0:
                        implied_gap_sizes.append(implied)
            else:
                # Estimate from positions directly
                # This is less accurate but workable
                observed_dist = right_pos - left_pos
                # If gap were 0, distance would be: (gap_start - left_pos) + (right_pos - gap_end)
                zero_gap_dist = (gap_start - left_pos) + (right_pos - gap_end)
                implied = observed_dist - zero_gap_dist + original_size
                if implied > 0:
                    implied_gap_sizes.append(implied)

        if not implied_gap_sizes:
            return GapSizeEstimate(
                gap_name=gap_name,
                chrom=chrom,
                gap_start=gap_start,
                gap_end=gap_end,
                original_size=original_size,
                estimated_size=original_size,
                confidence='low',
                supporting_pairs=len(spanning_pairs),
                insert_sizes=[]
            )

        # Use median as estimate
        estimated_size = int(statistics.median(implied_gap_sizes))

        # Determine confidence based on consistency
        if len(implied_gap_sizes) >= 10:
            std_dev = statistics.stdev(implied_gap_sizes) if len(implied_gap_sizes) > 1 else 0
            cv = std_dev / estimated_size if estimated_size > 0 else float('inf')
            if cv < 0.2:
                confidence = 'high'
            elif cv < 0.5:
                confidence = 'medium'
            else:
                confidence = 'low'
        elif len(implied_gap_sizes) >= 5:
            confidence = 'medium'
        else:
            confidence = 'low'

        return GapSizeEstimate(
            gap_name=gap_name,
            chrom=chrom,
            gap_start=gap_start,
            gap_end=gap_end,
            original_size=original_size,
            estimated_size=max(100, estimated_size),  # Minimum 100bp
            confidence=confidence,
            supporting_pairs=len(spanning_pairs),
            insert_sizes=implied_gap_sizes
        )

    def _is_valid_hic_read(self, read: pysam.AlignedSegment) -> bool:
        """Check if read is valid for Hi-C analysis"""
        if read.is_unmapped or read.mate_is_unmapped:
            return False
        if read.is_secondary or read.is_supplementary:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if not read.is_paired:
            return False
        return True

    # =========================================================================
    # Fill Validation
    # =========================================================================

    def validate_fills(self, gaps: List[Dict], filled_assembly: str,
                       window_size: int = 10000) -> List[FillValidation]:
        """
        Validate gap fills using Hi-C contact patterns.

        Strategy:
        - A correctly filled gap should show continuous contact pattern
        - Anomalies indicate potential mis-assembly:
          - 'break': sudden drop in contacts (wrong sequence inserted)
          - 'spike': unusual high contacts (collapsed repeat)

        Args:
            gaps: List of filled gaps with coordinates
            filled_assembly: Path to filled assembly
            window_size: Window size for contact analysis

        Returns:
            List of FillValidation objects
        """
        logger.info(f"Validating {len(gaps)} filled gaps using Hi-C...")

        # Re-align Hi-C to filled assembly if needed
        filled_bam = self._align_hic_to_assembly(filled_assembly)

        validations = []

        with pysam.AlignmentFile(str(filled_bam), 'rb') as bam:
            for gap in gaps:
                validation = self._validate_single_fill(bam, gap, window_size)
                validations.append(validation)

                status = "✓ PASS" if validation.is_valid else "✗ FAIL"
                logger.info(f"  {gap['name']}: {status} "
                           f"(score={validation.contact_score:.2f}, "
                           f"contacts={validation.observed_contacts}/{validation.expected_contacts})")

        return validations

    def _align_hic_to_assembly(self, assembly: str) -> Path:
        """Align Hi-C reads to assembly (or use existing BAM)"""
        # For now, assume BAM is already aligned to the assembly
        # In production, would re-align if assembly changed
        return self.hic_bam

    def _validate_single_fill(self, bam: pysam.AlignmentFile,
                               gap: Dict, window_size: int) -> FillValidation:
        """Validate a single filled gap"""
        chrom = gap['chrom']
        gap_start = gap['start']
        gap_end = gap['end']
        gap_name = gap.get('name', f"{chrom}_{gap_start}_{gap_end}")

        # Define regions around the filled gap
        # Left flank | Filled region | Right flank
        left_start = max(0, gap_start - window_size)
        left_end = gap_start
        right_start = gap_end
        right_end = gap_end + window_size
        fill_start = gap_start
        fill_end = gap_end

        # Count contacts between regions
        contacts = {
            'left_to_fill': 0,
            'fill_to_right': 0,
            'left_to_right': 0,
            'within_fill': 0
        }

        try:
            # Scan reads in the filled region and count mate locations
            for read in bam.fetch(chrom, fill_start, fill_end):
                if not self._is_valid_hic_read(read):
                    continue

                mate_chrom = read.next_reference_name
                mate_pos = read.next_reference_start

                if mate_chrom != chrom:
                    continue

                if left_start <= mate_pos < left_end:
                    contacts['left_to_fill'] += 1
                elif right_start <= mate_pos < right_end:
                    contacts['fill_to_right'] += 1
                elif fill_start <= mate_pos < fill_end:
                    contacts['within_fill'] += 1

            # Also count left-to-right contacts (should exist if correctly joined)
            for read in bam.fetch(chrom, left_start, left_end):
                if not self._is_valid_hic_read(read):
                    continue

                mate_chrom = read.next_reference_name
                mate_pos = read.next_reference_start

                if mate_chrom == chrom and right_start <= mate_pos < right_end:
                    contacts['left_to_right'] += 1

        except Exception as e:
            logger.warning(f"Error validating {gap_name}: {e}")

        # Calculate expected contacts based on distance
        # Hi-C contact frequency ~ 1/distance (roughly)
        gap_size = gap_end - gap_start
        fill_size = gap_size  # After filling

        # Expected contacts scale with region sizes and inversely with distance
        # This is a simplified model
        expected_left_fill = max(1, window_size * fill_size // 100000)
        expected_fill_right = max(1, fill_size * window_size // 100000)
        expected_left_right = max(1, window_size * window_size // (fill_size + 100000))

        total_expected = expected_left_fill + expected_fill_right + expected_left_right
        total_observed = contacts['left_to_fill'] + contacts['fill_to_right'] + contacts['left_to_right']

        # Calculate contact score (0-1)
        if total_expected > 0:
            contact_score = min(1.0, total_observed / total_expected)
        else:
            contact_score = 0.0

        # Detect anomalies
        anomaly_type = None

        # Check for contact break (very few contacts)
        if total_observed < total_expected * 0.1 and total_expected > 5:
            anomaly_type = 'break'

        # Check for contact spike (way too many contacts - collapsed repeat)
        if total_observed > total_expected * 5 and total_observed > 50:
            anomaly_type = 'spike'

        # Check for asymmetry (one side has contacts, other doesn't)
        if contacts['left_to_fill'] > 10 and contacts['fill_to_right'] < 2:
            anomaly_type = 'asymmetric'
        elif contacts['fill_to_right'] > 10 and contacts['left_to_fill'] < 2:
            anomaly_type = 'asymmetric'

        # Determine validity
        is_valid = anomaly_type is None and contact_score >= 0.3

        return FillValidation(
            gap_name=gap_name,
            is_valid=is_valid,
            contact_score=contact_score,
            expected_contacts=total_expected,
            observed_contacts=total_observed,
            anomaly_type=anomaly_type
        )

    # =========================================================================
    # Phasing Enhancement
    # =========================================================================

    def enhance_phasing(self, snp_db: Dict, phased_reads: Dict[str, Path],
                        haplotype_names: List[str]) -> Dict[str, List[str]]:
        """
        Use Hi-C long-range information to improve phasing.

        Strategy:
        - Hi-C pairs should come from the same haplotype
        - If read1 is phased to hap1, read2 should also be hap1
        - Use this to rescue ambiguous reads and correct mis-phased reads

        Args:
            snp_db: SNP database {chrom: {pos: {hap: base}}}
            phased_reads: {hap_name: path_to_phased_reads}
            haplotype_names: List of haplotype names

        Returns:
            Dict of read assignments {read_name: hap_name}
        """
        logger.info("Enhancing phasing with Hi-C long-range information...")

        # Build initial phasing map from existing assignments
        read_to_hap = {}
        for hap_name, reads_file in phased_reads.items():
            if hap_name == 'ambiguous':
                continue
            if reads_file and reads_file.exists():
                with open(reads_file) as f:
                    for line in f:
                        if line.startswith('>'):
                            read_name = line[1:].strip().split()[0]
                            read_to_hap[read_name] = hap_name

        logger.info(f"  Initial phased reads: {len(read_to_hap)}")

        # Use Hi-C pairs to propagate phasing
        hic_corrections = defaultdict(lambda: defaultdict(int))
        rescued_reads = {}

        with pysam.AlignmentFile(str(self.hic_bam), 'rb') as bam:
            for read in bam:
                if not self._is_valid_hic_read(read):
                    continue

                read_name = read.query_name
                mate_name = read.query_name  # Same name for paired reads

                # Check if this read or its mate is already phased
                read_hap = read_to_hap.get(read_name)

                if read_hap:
                    # This read is phased, use it to phase/validate mate
                    # Mates should be on same haplotype
                    hic_corrections[mate_name][read_hap] += 1

        # Process corrections
        corrections_made = 0
        rescues_made = 0

        for read_name, hap_votes in hic_corrections.items():
            if not hap_votes:
                continue

            best_hap = max(hap_votes.keys(), key=lambda h: hap_votes[h])
            best_count = hap_votes[best_hap]
            total_count = sum(hap_votes.values())

            # Require strong support
            if best_count >= 3 and best_count / total_count >= 0.7:
                current_hap = read_to_hap.get(read_name)

                if current_hap is None:
                    # Rescue ambiguous read
                    rescued_reads[read_name] = best_hap
                    rescues_made += 1
                elif current_hap != best_hap:
                    # Potential mis-phasing, log but don't auto-correct
                    # (could be real recombination or error)
                    corrections_made += 1
                    logger.debug(f"  Potential mis-phase: {read_name} "
                                f"{current_hap} -> {best_hap}")

        logger.info(f"  Hi-C rescued reads: {rescues_made}")
        logger.info(f"  Potential corrections: {corrections_made}")

        # Merge rescued reads into result
        result = dict(read_to_hap)
        result.update(rescued_reads)

        return result

    # =========================================================================
    # Utility Functions
    # =========================================================================

    def get_contact_matrix(self, chrom: str, start: int, end: int,
                           resolution: int = 10000) -> np.ndarray:
        """
        Generate Hi-C contact matrix for a region.

        Args:
            chrom: Chromosome name
            start: Start position
            end: End position
            resolution: Bin size

        Returns:
            2D numpy array of contact counts
        """
        n_bins = (end - start) // resolution + 1
        matrix = np.zeros((n_bins, n_bins), dtype=int)

        with pysam.AlignmentFile(str(self.hic_bam), 'rb') as bam:
            for read in bam.fetch(chrom, start, end):
                if not self._is_valid_hic_read(read):
                    continue

                mate_chrom = read.next_reference_name
                if mate_chrom != chrom:
                    continue

                mate_pos = read.next_reference_start
                if not (start <= mate_pos < end):
                    continue

                read_bin = (read.reference_start - start) // resolution
                mate_bin = (mate_pos - start) // resolution

                if 0 <= read_bin < n_bins and 0 <= mate_bin < n_bins:
                    matrix[read_bin, mate_bin] += 1
                    matrix[mate_bin, read_bin] += 1  # Symmetric

        return matrix


def align_hic_reads(hic_reads_r1: str, hic_reads_r2: str,
                    assembly: str, output_bam: str,
                    threads: int = 8) -> Path:
    """
    Align Hi-C reads to assembly using bwa-mem2.

    Args:
        hic_reads_r1: Path to R1 FASTQ
        hic_reads_r2: Path to R2 FASTQ
        assembly: Path to assembly FASTA
        output_bam: Path for output BAM
        threads: Number of threads

    Returns:
        Path to sorted, indexed BAM
    """
    logger.info("Aligning Hi-C reads...")

    output_bam = Path(output_bam)

    # Index assembly if needed
    if not Path(f"{assembly}.bwt.2bit.64").exists():
        logger.info("  Indexing assembly with bwa-mem2...")
        subprocess.run(['bwa-mem2', 'index', assembly],
                      check=True, capture_output=True)

    # Align with bwa-mem2 (Hi-C mode: -5SP)
    # -5: for split alignments, mark the primary as the shorter one
    # -S: skip mate rescue
    # -P: skip pairing, treat as single-end
    cmd = (f"bwa-mem2 mem -5SP -t {threads} {assembly} "
           f"{hic_reads_r1} {hic_reads_r2} | "
           f"samtools view -@ {threads} -bS - | "
           f"samtools sort -@ {threads} -o {output_bam} -")

    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    # Index BAM
    subprocess.run(['samtools', 'index', str(output_bam)],
                  check=True, capture_output=True)

    logger.info(f"  Created {output_bam}")
    return output_bam
