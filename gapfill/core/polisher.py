#!/usr/bin/env python3
"""
Flank Polisher - Polish gap flanking sequences to improve fillability

When gap flanks have issues (clip accumulation, high mismatch density),
polishing them can improve the chances of successful gap filling.

Polishing strategies:
1. Racon polish: Use reads to correct flank sequence errors
2. Reassembly: Re-assemble the flank region from reads (for severe issues)
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import pysam
from Bio import SeqIO


@dataclass
class PolishResult:
    """Result of flank polishing"""
    success: bool
    left_polished: bool = False
    right_polished: bool = False
    left_sequence: str = ""
    right_sequence: str = ""
    original_left: str = ""
    original_right: str = ""
    reason: str = ""


class FlankPolisher:
    """
    Polish gap flanking sequences using reads.

    Uses racon for consensus polishing. If racon is not available,
    falls back to majority-vote consensus from pileup.
    """

    # Default parameters
    FLANK_SIZE = 1000  # bp to polish on each side
    MIN_COVERAGE = 3   # Minimum coverage for polish
    MIN_READS = 5      # Minimum reads for polish attempt

    def __init__(self,
                 threads: int = 8,
                 flank_size: int = 1000,
                 work_dir: Optional[Path] = None):
        self.threads = threads
        self.flank_size = flank_size
        self.work_dir = work_dir or Path(tempfile.mkdtemp())
        self.logger = logging.getLogger(__name__)

        # Check for racon availability
        self.has_racon = self._check_racon()
        if not self.has_racon:
            self.logger.warning("racon not found, will use simple consensus")

    def _check_racon(self) -> bool:
        """Check if racon is available"""
        try:
            subprocess.run(['racon', '--version'],
                          capture_output=True, timeout=5)
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def polish_gap_flanks(self,
                          assembly_file: str,
                          bam_file: str,
                          chrom: str,
                          gap_start: int,
                          gap_end: int,
                          polish_left: bool = True,
                          polish_right: bool = True) -> PolishResult:
        """
        Polish the flanking sequences of a gap.

        Args:
            assembly_file: Path to assembly FASTA
            bam_file: Path to aligned reads BAM
            chrom: Chromosome name
            gap_start: Gap start position
            gap_end: Gap end position
            polish_left: Whether to polish left flank
            polish_right: Whether to polish right flank

        Returns:
            PolishResult with polished sequences
        """
        result = PolishResult(success=False)

        # Load assembly sequence
        assembly_seqs = {}
        for record in SeqIO.parse(assembly_file, 'fasta'):
            assembly_seqs[record.id] = str(record.seq)

        if chrom not in assembly_seqs:
            result.reason = f"Chromosome {chrom} not found in assembly"
            return result

        seq = assembly_seqs[chrom]

        # Define flank regions
        left_start = max(0, gap_start - self.flank_size)
        left_end = gap_start
        right_start = gap_end
        right_end = min(len(seq), gap_end + self.flank_size)

        # Extract original flanks
        result.original_left = seq[left_start:left_end]
        result.original_right = seq[right_start:right_end]

        try:
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                # Polish left flank
                if polish_left and left_end > left_start:
                    left_reads = self._collect_flank_reads(
                        bam, chrom, left_start, left_end
                    )
                    if len(left_reads) >= self.MIN_READS:
                        polished = self._polish_region(
                            result.original_left, left_reads, 'left'
                        )
                        if polished and polished != result.original_left:
                            result.left_sequence = polished
                            result.left_polished = True
                            self.logger.info(f"    Left flank polished: "
                                           f"{len(result.original_left)}bp -> {len(polished)}bp")
                    else:
                        self.logger.debug(f"    Left flank: insufficient reads ({len(left_reads)})")

                # Polish right flank
                if polish_right and right_end > right_start:
                    right_reads = self._collect_flank_reads(
                        bam, chrom, right_start, right_end
                    )
                    if len(right_reads) >= self.MIN_READS:
                        polished = self._polish_region(
                            result.original_right, right_reads, 'right'
                        )
                        if polished and polished != result.original_right:
                            result.right_sequence = polished
                            result.right_polished = True
                            self.logger.info(f"    Right flank polished: "
                                           f"{len(result.original_right)}bp -> {len(polished)}bp")
                    else:
                        self.logger.debug(f"    Right flank: insufficient reads ({len(right_reads)})")

            result.success = result.left_polished or result.right_polished
            if result.success:
                result.reason = "Polishing successful"
            else:
                result.reason = "No changes made (insufficient reads or no improvement)"

        except Exception as e:
            result.reason = f"Polishing error: {e}"
            self.logger.warning(result.reason)

        return result

    def _collect_flank_reads(self,
                             bam: pysam.AlignmentFile,
                             chrom: str,
                             start: int,
                             end: int) -> List[Tuple[str, str]]:
        """Collect reads covering the flank region"""
        reads = []
        seen_names = set()

        try:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.query_name in seen_names:
                    continue

                seq = read.query_sequence
                if seq and len(seq) >= 500:
                    reads.append((seq, read.query_name))
                    seen_names.add(read.query_name)
        except Exception as e:
            self.logger.debug(f"Error collecting reads: {e}")

        return reads

    def _polish_region(self,
                       reference: str,
                       reads: List[Tuple[str, str]],
                       name: str) -> Optional[str]:
        """
        Polish a region using reads.

        Uses racon if available, otherwise falls back to simple consensus.
        """
        if len(reference) < 100 or len(reads) < self.MIN_READS:
            return None

        work_dir = Path(self.work_dir) / f"polish_{name}"
        work_dir.mkdir(parents=True, exist_ok=True)

        # Write reference
        ref_file = work_dir / "ref.fa"
        with open(ref_file, 'w') as f:
            f.write(f">ref\n{reference}\n")

        # Write reads
        reads_file = work_dir / "reads.fa"
        with open(reads_file, 'w') as f:
            for seq, read_name in reads:
                f.write(f">{read_name}\n{seq}\n")

        if self.has_racon:
            return self._polish_with_racon(ref_file, reads_file, work_dir)
        else:
            return self._polish_with_pileup(ref_file, reads_file, work_dir)

    def _polish_with_racon(self,
                           ref_file: Path,
                           reads_file: Path,
                           work_dir: Path) -> Optional[str]:
        """Polish using minimap2 + racon"""
        try:
            # Align reads to reference
            paf_file = work_dir / "align.paf"
            result = subprocess.run(
                ['minimap2', '-x', 'map-hifi', '-t', str(self.threads),
                 str(ref_file), str(reads_file)],
                capture_output=True, text=True, timeout=120
            )

            if result.returncode != 0:
                return None

            with open(paf_file, 'w') as f:
                f.write(result.stdout)

            # Check alignment count
            alignment_count = len([l for l in result.stdout.strip().split('\n') if l])
            if alignment_count < 3:
                self.logger.debug(f"    Too few alignments ({alignment_count})")
                return None

            # Run racon
            result = subprocess.run(
                ['racon', '-t', str(self.threads),
                 str(reads_file), str(paf_file), str(ref_file)],
                capture_output=True, text=True, timeout=300
            )

            if result.returncode != 0:
                return None

            # Parse polished sequence
            polished = ''
            for line in result.stdout.split('\n'):
                if not line.startswith('>'):
                    polished += line.strip()

            # Sanity check: polished should be similar length
            if polished and 0.8 <= len(polished) / len(open(ref_file).read().split('\n')[1]) <= 1.2:
                return polished

            return None

        except subprocess.TimeoutExpired:
            self.logger.debug("    Racon timed out")
            return None
        except Exception as e:
            self.logger.debug(f"    Racon error: {e}")
            return None

    def _polish_with_pileup(self,
                            ref_file: Path,
                            reads_file: Path,
                            work_dir: Path) -> Optional[str]:
        """
        Simple pileup-based consensus as fallback.

        Uses minimap2 alignment and majority vote at each position.
        """
        try:
            # Read reference
            with open(ref_file) as f:
                f.readline()  # Skip header
                ref_seq = f.read().replace('\n', '')

            # Align and create BAM
            bam_file = work_dir / "align.bam"
            cmd = (
                f"minimap2 -ax map-hifi -t {self.threads} "
                f"{ref_file} {reads_file} | "
                f"samtools sort -@ {self.threads} -o {bam_file} - && "
                f"samtools index {bam_file}"
            )

            result = subprocess.run(cmd, shell=True, capture_output=True,
                                   text=True, timeout=120)
            if result.returncode != 0:
                return None

            # Pileup consensus
            consensus = list(ref_seq)

            with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
                for pileup_col in bam.pileup('ref', 0, len(ref_seq)):
                    pos = pileup_col.reference_pos
                    if pos >= len(consensus):
                        continue

                    # Count bases
                    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                    for read in pileup_col.pileups:
                        if read.is_del or read.is_refskip:
                            continue
                        base = read.alignment.query_sequence[read.query_position].upper()
                        if base in base_counts:
                            base_counts[base] += 1

                    # Majority vote
                    total = sum(base_counts.values())
                    if total >= self.MIN_COVERAGE:
                        max_base = max(base_counts, key=base_counts.get)
                        if base_counts[max_base] / total >= 0.6:
                            consensus[pos] = max_base

            return ''.join(consensus)

        except Exception as e:
            self.logger.debug(f"    Pileup consensus error: {e}")
            return None

    def polish_assembly_flanks(self,
                               assembly_file: str,
                               bam_file: str,
                               gaps_to_polish: List[Dict],
                               output_file: str) -> Dict[str, PolishResult]:
        """
        Polish flanks for multiple gaps and write updated assembly.

        Args:
            assembly_file: Input assembly FASTA
            bam_file: Aligned reads BAM
            gaps_to_polish: List of gap dicts with 'chrom', 'start', 'end', 'name',
                           and optionally 'polish_left', 'polish_right'
            output_file: Output polished assembly FASTA

        Returns:
            Dict mapping gap_name to PolishResult
        """
        results = {}

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly_file, 'fasta'):
            sequences[record.id] = list(str(record.seq))

        # Process each gap
        for gap in gaps_to_polish:
            gap_name = gap.get('name', f"{gap['chrom']}_{gap['start']}_{gap['end']}")
            chrom = gap['chrom']
            gap_start = gap['start']
            gap_end = gap['end']
            polish_left = gap.get('polish_left', True)
            polish_right = gap.get('polish_right', True)

            self.logger.info(f"  Polishing flanks for {gap_name}")

            result = self.polish_gap_flanks(
                assembly_file, bam_file, chrom, gap_start, gap_end,
                polish_left, polish_right
            )
            results[gap_name] = result

            if not result.success:
                continue

            # Apply polished sequences
            if chrom not in sequences:
                continue

            seq = sequences[chrom]

            # Apply right first (to avoid coordinate shift)
            if result.right_polished:
                right_start = gap_end
                right_end = min(len(seq), gap_end + self.flank_size)
                # Replace the right flank region
                seq[right_start:right_end] = list(result.right_sequence)

            if result.left_polished:
                left_start = max(0, gap_start - self.flank_size)
                left_end = gap_start
                # Replace the left flank region
                seq[left_start:left_end] = list(result.left_sequence)

            sequences[chrom] = seq

        # Write polished assembly
        with open(output_file, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                seq_str = ''.join(seq)
                for i in range(0, len(seq_str), 80):
                    f.write(seq_str[i:i+80] + '\n')

        # Summary
        polished_count = sum(1 for r in results.values() if r.success)
        self.logger.info(f"  Polished {polished_count}/{len(gaps_to_polish)} gaps")

        return results
