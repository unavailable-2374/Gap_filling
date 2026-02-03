#!/usr/bin/env python3
"""
Polyploid Gap Filling Engine

Gap filling for polyploid genomes (diploid, tetraploid, hexaploid, etc.)
with integrated phasing support.

Key features:
1. Automatic ploidy detection from number of haplotype assemblies
2. Gap normalization BEFORE SNP detection (all gaps -> 500N)
3. Alignment-based SNP detection (handles different gap positions)
4. Support for haplotype-specific gaps
5. SNP-based read phasing (builtin or WhatsHap)
6. Independent gap filling per haplotype
7. Full HiFi + ONT utilization (both data types phased independently)
"""

import logging
import subprocess
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
from datetime import datetime
from collections import defaultdict

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from gapfill.engines.haploid import HaploidEngine, _mark_done, _is_done, _clear_marker
from gapfill.core.validator import GapStatusTracker, GapStatus
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads
from gapfill.utils.checkpoint import PolyploidCheckpointManager, CheckpointState
from gapfill.utils.reads_cache import ReadsCache


@dataclass
class GapRegion:
    """Represents a gap region in an assembly"""
    chrom: str
    start: int
    end: int
    size: int

    @property
    def name(self) -> str:
        return f"{self.chrom}_{self.start}_{self.end}"


class ReadPhaser:
    """Phase reads to haplotypes based on haplotype-specific SNPs"""

    def __init__(self, haplotype_assemblies: List[Path],
                 threads: int = 8,
                 min_snp_qual: int = 20,
                 min_read_snps: int = 1,  # Reduced from 2 to 1 for better phasing
                 work_dir: Optional[Path] = None):

        self.haplotypes = haplotype_assemblies
        self.num_haplotypes = len(haplotype_assemblies)
        self.threads = threads
        self.min_snp_qual = min_snp_qual
        self.min_read_snps = min_read_snps
        self.work_dir = work_dir or Path('.')

        self.logger = logging.getLogger(__name__)
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        # Store normalized assemblies and gap info
        self.normalized_assemblies: Dict[str, Path] = {}
        self.gap_regions: Dict[str, List[GapRegion]] = {}  # hap_name -> gaps
        self.haplotype_specific_gaps: Dict[str, Set[str]] = {}  # hap_name -> gap names unique to this hap

        # Phasing statistics
        self.phasing_stats = {
            'total_snp_positions': 0,
            'reads_with_0_snps': 0,
            'reads_with_1_snp': 0,
            'reads_with_2plus_snps': 0,
            'reads_with_tie': 0
        }

    def normalize_all_assemblies(self, min_gap_size: int = 100) -> Dict[str, Path]:
        """
        Normalize gaps in all haplotype assemblies to 500N.
        Must be called BEFORE SNP detection.

        Returns:
            Dict mapping hap_name to normalized assembly path
        """
        self.logger.info("Normalizing gaps in all haplotype assemblies...")

        for i, hap_file in enumerate(self.haplotypes):
            hap_name = self.hap_names[i]

            # Find gaps
            gaps = self._find_gaps(hap_file, min_gap_size)
            self.gap_regions[hap_name] = gaps
            self.logger.info(f"  {hap_name}: found {len(gaps)} gaps")

            # Normalize
            normalized_file = self.work_dir / f"{hap_name}_normalized.fasta"
            self._normalize_gaps(hap_file, gaps, normalized_file)
            self.normalized_assemblies[hap_name] = normalized_file

        # Identify haplotype-specific gaps
        self._identify_haplotype_specific_gaps()

        return self.normalized_assemblies

    def _find_gaps(self, assembly_file: Path, min_size: int = 100) -> List[GapRegion]:
        """Find all gaps (N-runs) in an assembly"""
        gaps = []

        for record in SeqIO.parse(assembly_file, 'fasta'):
            seq = str(record.seq).upper()

            for match in re.finditer(r'N+', seq):
                gap_size = match.end() - match.start()
                if gap_size >= min_size:
                    gaps.append(GapRegion(
                        chrom=record.id,
                        start=match.start(),
                        end=match.end(),
                        size=gap_size
                    ))

        return gaps

    def _normalize_gaps(self, input_file: Path, gaps: List[GapRegion],
                        output_file: Path, normalized_size: int = 500):
        """Normalize all gaps to a fixed size (default 500N)"""
        sequences = {}
        for record in SeqIO.parse(input_file, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group gaps by chromosome
        gaps_by_chrom = defaultdict(list)
        for gap in gaps:
            gaps_by_chrom[gap.chrom].append(gap)

        # Normalize each chromosome (process from end to start)
        for chrom, chrom_gaps in gaps_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]
            chrom_gaps_sorted = sorted(chrom_gaps, key=lambda x: x.start, reverse=True)

            for gap in chrom_gaps_sorted:
                seq = seq[:gap.start] + 'N' * normalized_size + seq[gap.end:]

            sequences[chrom] = seq

        # Write normalized assembly
        with open(output_file, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

    def _identify_haplotype_specific_gaps(self):
        """
        Identify gaps that exist in only one haplotype.
        Uses chromosome name and approximate position to match gaps.
        """
        # Build gap index: chrom -> list of (hap_name, gap)
        gap_index = defaultdict(list)
        for hap_name, gaps in self.gap_regions.items():
            for gap in gaps:
                gap_index[gap.chrom].append((hap_name, gap))

        # For each haplotype, find gaps that don't have a corresponding gap in other haplotypes
        for hap_name in self.hap_names:
            self.haplotype_specific_gaps[hap_name] = set()

        position_tolerance = 1000  # Allow 1kb tolerance for matching gaps

        for chrom, hap_gaps in gap_index.items():
            for hap_name, gap in hap_gaps:
                # Check if any other haplotype has a gap at similar position
                has_match = False
                for other_hap, other_gap in hap_gaps:
                    if other_hap == hap_name:
                        continue
                    # Check if gaps overlap or are close enough
                    if abs(gap.start - other_gap.start) <= position_tolerance:
                        has_match = True
                        break

                if not has_match:
                    self.haplotype_specific_gaps[hap_name].add(gap.name)
                    self.logger.info(f"  Haplotype-specific gap: {hap_name} {gap.name}")

    def detect_haplotype_snps(self, method: str = 'builtin') -> Dict:
        """
        Detect haplotype-specific SNPs using alignment.
        Must call normalize_all_assemblies() first.
        """
        self.logger.info("Detecting haplotype-specific SNPs...")

        if not self.normalized_assemblies:
            raise RuntimeError("Must call normalize_all_assemblies() before SNP detection")

        if method == 'whatshap':
            return self._detect_snps_whatshap()
        else:
            return self._detect_snps_alignment()

    def _detect_snps_alignment(self) -> Dict:
        """
        Alignment-based SNP detection.
        Aligns each haplotype to hap1 (reference) using minimap2.
        Correctly handles different gap positions between haplotypes.

        Two-pass approach:
        1. First pass: Detect all SNP positions from all haplotype comparisons
        2. Second pass: Re-scan alignments to fill in matching bases at SNP positions

        This correctly handles coordinate differences caused by gap normalization.
        """
        snp_db = {}

        ref_hap = self.hap_names[0]
        ref_file = self.normalized_assemblies[ref_hap]

        # Load ALL haplotype sequences upfront
        all_hap_seqs = {}
        for hap_name in self.hap_names:
            hap_file = self.normalized_assemblies[hap_name]
            all_hap_seqs[hap_name] = {}
            for record in SeqIO.parse(hap_file, 'fasta'):
                all_hap_seqs[hap_name][record.id] = str(record.seq).upper()

        ref_seqs = all_hap_seqs[ref_hap]

        # Build gap regions set for quick lookup
        ref_gap_regions = self._build_gap_position_set(ref_seqs)

        # Initialize SNP database
        for chrom in ref_seqs:
            snp_db[chrom] = {}

        # Store BAM files for second pass
        bam_files = {}

        # PASS 1: Detect all SNP positions from all haplotype comparisons
        self.logger.info("  Pass 1: Detecting SNP positions...")
        for other_hap in self.hap_names[1:]:
            other_file = self.normalized_assemblies[other_hap]

            self.logger.info(f"    Aligning {other_hap} to {ref_hap}...")

            other_seqs = all_hap_seqs[other_hap]
            other_gap_regions = self._build_gap_position_set(other_seqs)

            # Run minimap2 alignment
            bam_file = self.work_dir / f"align_{other_hap}_to_{ref_hap}.bam"
            self._align_assemblies(other_file, ref_file, bam_file)
            bam_files[other_hap] = bam_file

            # Parse alignment to find SNPs (only detect differences)
            snps_found = self._extract_snps_from_alignment(
                bam_file, ref_seqs, other_seqs,
                ref_gap_regions, other_gap_regions,
                ref_hap, other_hap, snp_db,
                record_matches=False  # First pass: only record differences
            )

            self.logger.info(f"      Found {snps_found} SNPs between {ref_hap} and {other_hap}")

        # PASS 2: Fill in matching bases at known SNP positions
        # This uses alignment to get correct coordinates for each haplotype
        self.logger.info("  Pass 2: Filling all haplotype bases at SNP positions...")
        positions_completed = 0

        for other_hap in self.hap_names[1:]:
            bam_file = bam_files.get(other_hap)
            if not bam_file or not bam_file.exists():
                continue

            other_seqs = all_hap_seqs[other_hap]
            other_gap_regions = self._build_gap_position_set(other_seqs)

            # Re-scan alignment to fill matching bases at known SNP positions
            filled = self._extract_snps_from_alignment(
                bam_file, ref_seqs, other_seqs,
                ref_gap_regions, other_gap_regions,
                ref_hap, other_hap, snp_db,
                record_matches=True  # Second pass: record matches at known SNP positions
            )
            positions_completed += filled

            # Clean up BAM file
            bam_file.unlink(missing_ok=True)
            Path(str(bam_file) + ".bai").unlink(missing_ok=True)

        self.logger.info(f"  Filled {positions_completed} matching haplotype entries")

        total_snps = sum(len(positions) for positions in snp_db.values())
        self.logger.info(f"  Total: {total_snps} haplotype-specific SNP positions")

        return snp_db

    def _build_gap_position_set(self, sequences: Dict[str, str]) -> Dict[str, Set[int]]:
        """Build a set of positions that are within gap regions (N's)"""
        gap_positions = {}
        for chrom, seq in sequences.items():
            gap_positions[chrom] = set()
            for match in re.finditer(r'N+', seq):
                for pos in range(match.start(), match.end()):
                    gap_positions[chrom].add(pos)
        return gap_positions

    def _align_assemblies(self, query_file: Path, ref_file: Path, output_bam: Path):
        """Align two assemblies using minimap2 asm5 preset"""
        try:
            # asm5: for ~0.1% divergence (same-species haplotypes)
            # -m 2G: memory limit per thread for samtools sort
            cmd = (
                f"minimap2 -ax asm5 -t {self.threads} "
                f"{ref_file} {query_file} | "
                f"samtools sort -@ {self.threads} -m 2G -o {output_bam} - && "
                f"samtools index {output_bam}"
            )

            subprocess.run(cmd, shell=True, check=True,
                          capture_output=True, text=True)

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed: {e.stderr}")
            raise

    def _extract_snps_from_alignment(self, bam_file: Path,
                                      ref_seqs: Dict[str, str],
                                      query_seqs: Dict[str, str],
                                      ref_gap_positions: Dict[str, Set[int]],
                                      query_gap_positions: Dict[str, Set[int]],
                                      ref_hap: str, query_hap: str,
                                      snp_db: Dict,
                                      record_matches: bool = False) -> int:
        """
        Extract SNPs from alignment BAM file.
        Skips positions within gap regions.

        Args:
            record_matches: If False (pass 1), only record positions where bases differ.
                          If True (pass 2), also record matching bases at known SNP positions.

        Returns:
            Number of entries added to snp_db
        """
        entries_added = 0

        with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                ref_name = read.reference_name
                if ref_name not in snp_db:
                    continue

                query_name = read.query_name
                if query_name not in query_seqs:
                    continue

                query_seq = query_seqs[query_name]
                ref_gap_set = ref_gap_positions.get(ref_name, set())
                query_gap_set = query_gap_positions.get(query_name, set())

                # Use aligned pairs to find mismatches
                try:
                    for query_pos, ref_pos in read.get_aligned_pairs():
                        if query_pos is None or ref_pos is None:
                            continue

                        # Skip gap regions
                        if ref_pos in ref_gap_set:
                            continue
                        if query_pos in query_gap_set:
                            continue

                        ref_base = ref_seqs[ref_name][ref_pos].upper()
                        query_base = read.query_sequence[query_pos].upper()

                        # Skip if either is N
                        if ref_base == 'N' or query_base == 'N':
                            continue

                        if not record_matches:
                            # PASS 1: Only record positions where bases differ
                            if ref_base != query_base:
                                if ref_pos not in snp_db[ref_name]:
                                    snp_db[ref_name][ref_pos] = {ref_hap: ref_base}
                                snp_db[ref_name][ref_pos][query_hap] = query_base
                                entries_added += 1
                        else:
                            # PASS 2: Record this haplotype's base at known SNP positions
                            # This handles cases where this haplotype matches ref at a SNP position
                            if ref_pos in snp_db[ref_name]:
                                if query_hap not in snp_db[ref_name][ref_pos]:
                                    snp_db[ref_name][ref_pos][query_hap] = query_base
                                    entries_added += 1

                except Exception as e:
                    self.logger.debug(f"Error processing alignment: {e}")
                    continue

        return entries_added

    def _detect_snps_builtin(self) -> Dict:
        """
        Legacy builtin SNP detection (position-by-position comparison).
        Only works correctly if gaps are at identical positions.
        Kept for backwards compatibility.
        """
        self.logger.warning("Using legacy position-based SNP detection. "
                           "Consider using alignment-based detection for robustness.")

        snp_db = {}

        # Load all haplotype sequences from NORMALIZED assemblies
        hap_seqs = {}
        for hap_name, hap_file in self.normalized_assemblies.items():
            hap_seqs[hap_name] = {}
            for record in SeqIO.parse(hap_file, 'fasta'):
                hap_seqs[hap_name][record.id] = str(record.seq).upper()

        ref_hap = self.hap_names[0]

        for chrom in hap_seqs[ref_hap]:
            snp_db[chrom] = {}
            ref_seq = hap_seqs[ref_hap][chrom]

            for other_hap in self.hap_names[1:]:
                if chrom not in hap_seqs[other_hap]:
                    continue

                other_seq = hap_seqs[other_hap][chrom]
                min_len = min(len(ref_seq), len(other_seq))

                for pos in range(min_len):
                    ref_base = ref_seq[pos]
                    other_base = other_seq[pos]

                    if ref_base != other_base and ref_base != 'N' and other_base != 'N':
                        if pos not in snp_db[chrom]:
                            snp_db[chrom][pos] = {ref_hap: ref_base}
                        snp_db[chrom][pos][other_hap] = other_base

        total_snps = sum(len(positions) for positions in snp_db.values())
        self.logger.info(f"  Found {total_snps} haplotype-specific SNP positions")

        return snp_db

    def _detect_snps_whatshap(self) -> Dict:
        """Use WhatsHap for SNP detection"""
        self.logger.info("  Using WhatsHap for SNP detection...")
        # Simplified - would need bcftools and whatshap installed
        return self._detect_snps_alignment()

    def phase_reads_from_bam(self, bam_file: Path, snp_db: Dict,
                             output_prefix: Path, read_type: str = 'reads',
                             filtered_read_names: Optional[Set[str]] = None) -> Dict[str, Path]:
        """
        Phase reads based on SNP profiles.

        Args:
            bam_file: BAM file with aligned reads
            snp_db: SNP database
            output_prefix: Output path prefix
            read_type: 'hifi' or 'ont'
            filtered_read_names: If provided, only output reads in this set
                                (but still use all reads for phasing decisions)
        """
        self.logger.info(f"Phasing {read_type} reads to haplotypes...")
        self.logger.info(f"  min_read_snps threshold: {self.min_read_snps}")
        if filtered_read_names:
            self.logger.info(f"  Will output only {len(filtered_read_names):,} filtered reads")

        # Count total SNP positions for diagnostics
        total_snp_positions = sum(len(positions) for positions in snp_db.values())
        self.logger.info(f"  Total SNP positions in database: {total_snp_positions}")

        # Diagnostic: Check SNP database coverage for each haplotype
        hap_coverage = {hap: 0 for hap in self.hap_names}
        for chrom, positions in snp_db.items():
            for ref_pos, hap_bases in positions.items():
                for hap_name in hap_bases:
                    if hap_name in hap_coverage:
                        hap_coverage[hap_name] += 1
        self.logger.info(f"  SNP entries per haplotype: {hap_coverage}")

        phased_reads = {hap: [] for hap in self.hap_names}
        ambiguous_reads = []

        stats = {
            'total': 0,
            'phased': {hap: 0 for hap in self.hap_names},
            'ambiguous': 0,
            'no_snp_overlap': 0,  # Reads that don't overlap any SNP
            'tie_scores': 0,      # Reads with tied scores
            'below_threshold': 0,  # Reads below min_read_snps
            'filtered_out': 0     # Reads not in filtered set
        }

        try:
            with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
                for read in bam:
                    if read.is_unmapped or read.is_secondary:
                        continue

                    stats['total'] += 1

                    hap_assignment, reason = self._assign_read_to_haplotype_detailed(read, snp_db)

                    # Check if read passes filter (if filtering enabled)
                    passes_filter = (filtered_read_names is None or
                                    read.query_name in filtered_read_names)

                    if hap_assignment:
                        if passes_filter:
                            phased_reads[hap_assignment].append(read)
                            stats['phased'][hap_assignment] += 1
                        else:
                            stats['filtered_out'] += 1
                    else:
                        if passes_filter:
                            ambiguous_reads.append(read)
                            stats['ambiguous'] += 1
                        else:
                            stats['filtered_out'] += 1

                        if reason == 'no_snp':
                            stats['no_snp_overlap'] += 1
                        elif reason == 'tie':
                            stats['tie_scores'] += 1
                        elif reason == 'below_threshold':
                            stats['below_threshold'] += 1

        except Exception as e:
            self.logger.error(f"Error phasing reads: {e}")
            raise

        # Write phased reads
        output_files = {}

        for hap_name in self.hap_names:
            output_file = Path(f"{output_prefix}_{hap_name}_{read_type}.fasta")

            with open(output_file, 'w') as f:
                for read in phased_reads[hap_name]:
                    seq = read.query_sequence
                    if seq:
                        f.write(f">{read.query_name}\n{seq}\n")

            output_files[hap_name] = output_file
            pct = stats['phased'][hap_name] / stats['total'] * 100 if stats['total'] > 0 else 0
            self.logger.info(f"  {hap_name}: {stats['phased'][hap_name]} reads ({pct:.1f}%)")

        # Write ambiguous reads
        ambiguous_file = Path(f"{output_prefix}_ambiguous_{read_type}.fasta")
        with open(ambiguous_file, 'w') as f:
            for read in ambiguous_reads:
                seq = read.query_sequence
                if seq:
                    f.write(f">{read.query_name}\n{seq}\n")

        output_files['ambiguous'] = ambiguous_file

        # Detailed statistics
        total = stats['total']
        phased_total = sum(stats['phased'].values())
        self.logger.info(f"  --- Phasing Summary ---")
        self.logger.info(f"  Total reads: {total}")
        self.logger.info(f"  Phased (output): {phased_total} ({phased_total/total*100:.1f}%)" if total > 0 else "  Phased: 0")
        self.logger.info(f"  Ambiguous (output): {stats['ambiguous']} ({stats['ambiguous']/total*100:.1f}%)" if total > 0 else "  Ambiguous: 0")
        if stats['filtered_out'] > 0:
            self.logger.info(f"  Filtered out (far from gaps): {stats['filtered_out']}")
        self.logger.info(f"    - No SNP overlap: {stats['no_snp_overlap']}")
        self.logger.info(f"    - Below threshold ({self.min_read_snps}): {stats['below_threshold']}")
        self.logger.info(f"    - Tied scores: {stats['tie_scores']}")

        return output_files

    def _assign_read_to_haplotype_detailed(self, read: pysam.AlignedSegment, snp_db: Dict) -> tuple:
        """
        Assign a read to a haplotype based on SNP profile.
        Returns (haplotype_name, None) or (None, reason).
        """
        chrom = read.reference_name

        if chrom not in snp_db:
            return None, 'no_snp'

        hap_scores = {hap: 0 for hap in self.hap_names}
        snps_checked = 0

        try:
            for query_pos, ref_pos in read.get_aligned_pairs():
                if ref_pos is None or query_pos is None:
                    continue

                if ref_pos in snp_db[chrom]:
                    snps_checked += 1
                    read_base = read.query_sequence[query_pos].upper()

                    # Get reference (hap1) base for fallback
                    ref_hap = self.hap_names[0]
                    ref_base_at_pos = snp_db[chrom][ref_pos].get(ref_hap)

                    # Score each haplotype
                    for hap_name in self.hap_names:
                        if hap_name in snp_db[chrom][ref_pos]:
                            hap_base = snp_db[chrom][ref_pos][hap_name]
                        else:
                            # CRITICAL FIX: If haplotype not in database at this position,
                            # it means it has the SAME base as reference (hap1).
                            # Only differences are recorded during SNP detection.
                            hap_base = ref_base_at_pos

                        if hap_base and read_base == hap_base:
                            hap_scores[hap_name] += 1

        except Exception:
            return None, 'error'

        if snps_checked == 0:
            return None, 'no_snp'

        if snps_checked < self.min_read_snps:
            return None, 'below_threshold'

        max_score = max(hap_scores.values())
        if max_score == 0:
            return None, 'no_match'

        best_haps = [h for h, s in hap_scores.items() if s == max_score]

        if len(best_haps) == 1:
            return best_haps[0], None
        else:
            return None, 'tie'

class PolyploidEngine:
    """
    Gap filler for polyploid genomes

    Coordinates phasing and independent gap filling for each haplotype.
    Now fully utilizes both HiFi and ONT data when available.

    Key workflow changes:
    1. Normalize gaps in ALL haplotypes BEFORE SNP detection
    2. Use alignment-based SNP detection (handles gap position differences)
    3. Support haplotype-specific gaps
    """

    def __init__(self,
                 haplotype_assemblies: List[str],
                 hifi_reads: Optional[str] = None,
                 ont_reads: Optional[str] = None,
                 hic_reads: Optional[List[str]] = None,
                 hic_bam: Optional[str] = None,
                 output_dir: str = "polyploid_output",
                 threads: int = 8,
                 max_iterations: int = 10,
                 phasing_method: str = 'builtin',
                 use_ambiguous_reads: bool = True,
                 min_gap_size: int = 100,
                 min_mapq: int = 20,
                 min_read_snps: int = 1,
                 resume: bool = False,
                 clear_checkpoint: bool = False,
                 optimized_mode: bool = True):

        self.haplotypes = [Path(h) for h in haplotype_assemblies]
        self.num_haplotypes = len(self.haplotypes)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.hic_reads = hic_reads  # [R1, R2] or None
        self.hic_bam = Path(hic_bam) if hic_bam else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.max_iterations = max_iterations
        self.phasing_method = phasing_method
        self.use_ambiguous_reads = use_ambiguous_reads
        self.min_gap_size = min_gap_size
        self.min_mapq = min_mapq
        self.min_read_snps = min_read_snps
        self.resume = resume
        self.optimized_mode = optimized_mode

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        # Checkpoint manager
        self.checkpoint = PolyploidCheckpointManager(str(self.output_dir), self.hap_names)
        if clear_checkpoint:
            self.checkpoint.clear()

        # Hi-C analyzer (initialized when needed)
        self.hic_analyzer: Optional[HiCAnalyzer] = None

        # Store normalized assemblies
        self.normalized_assemblies: Dict[str, Path] = {}

        # Reads cache for optimized mode (filter once, use for phasing)
        self.reads_cache: Optional[ReadsCache] = None
        if self.optimized_mode:
            self.reads_cache = ReadsCache(self.output_dir, threads=threads)

        self._validate_inputs()

        self.logger.info("=" * 60)
        self.logger.info("PolyploidEngine initialized")
        self.logger.info(f"  Ploidy: {self.num_haplotypes}n")
        self.logger.info(f"  Haplotypes: {[h.name for h in self.haplotypes]}")
        self.logger.info(f"  HiFi reads: {self.hifi_reads}")
        self.logger.info(f"  ONT reads: {self.ont_reads}")
        self.logger.info(f"  Hi-C reads: {self.hic_reads}")
        self.logger.info(f"  Hi-C BAM: {self.hic_bam}")
        self.logger.info(f"  Phasing method: {self.phasing_method}")
        self.logger.info(f"  Optimized mode: {self.optimized_mode}")
        self.logger.info(f"  Resume: {self.resume}")
        self.logger.info("=" * 60)

    def _validate_inputs(self):
        for hap in self.haplotypes:
            if not hap.exists():
                raise FileNotFoundError(f"Haplotype file not found: {hap}")

        if not self.hifi_reads and not self.ont_reads:
            raise ValueError("At least one reads file required")

    def run(self) -> Dict[str, Path]:
        """Run polyploid gap filling with full HiFi + ONT utilization"""
        self.logger.info("Starting polyploid gap filling...")

        # Check for existing checkpoint
        checkpoint_state = None
        if self.resume:
            if self.checkpoint.exists():
                checkpoint_state = self.checkpoint.load()
                if checkpoint_state:
                    self.logger.info("=" * 60)
                    self.logger.info("RESUMING FROM CHECKPOINT")
                    self.logger.info(f"  Phase: {checkpoint_state.phase}")
                    self.logger.info("=" * 60)
            else:
                # No checkpoint.json but --resume specified
                # Scan for existing files from previous run
                self.logger.info("=" * 60)
                self.logger.info("No checkpoint.json found, scanning existing files...")
                checkpoint_state = self.checkpoint.scan_existing_files()
                if checkpoint_state.phase != "init":
                    self.logger.info(f"  Detected phase: {checkpoint_state.phase}")
                    self.checkpoint.save(checkpoint_state)
                else:
                    self.logger.info("  No usable files found, starting fresh")
                    checkpoint_state = None
                self.logger.info("=" * 60)

        # Initialize checkpoint state if not resuming
        if not checkpoint_state:
            checkpoint_state = CheckpointState(
                engine="polyploid",
                phase="init",
                max_iterations=self.max_iterations
            )
            self.checkpoint.save(checkpoint_state)
            self.checkpoint.init_haplotype_states()

        # Initialize phaser
        phaser = ReadPhaser(
            self.haplotypes,
            threads=self.threads,
            min_read_snps=self.min_read_snps,
            work_dir=self.output_dir
        )

        # STEP 0: Normalize gaps in ALL haplotypes FIRST
        if _is_done(self.output_dir, "normalized"):
            self.logger.info("STEP 0: Reusing normalized assemblies")
            for hap_name in self.hap_names:
                normalized_file = self.output_dir / f"{hap_name}_normalized.fasta"
                self.normalized_assemblies[hap_name] = normalized_file
                phaser.normalized_assemblies[hap_name] = normalized_file
                phaser.gap_regions[hap_name] = phaser._find_gaps(
                    normalized_file, self.min_gap_size
                )
        else:
            self.logger.info("STEP 0: Normalizing gaps in all haplotypes")
            self.normalized_assemblies = phaser.normalize_all_assemblies(
                min_gap_size=self.min_gap_size
            )
            _mark_done(self.output_dir, "normalized")

        # Report haplotype-specific gaps
        for hap_name, specific_gaps in phaser.haplotype_specific_gaps.items():
            if specific_gaps:
                self.logger.info(f"  {hap_name} has {len(specific_gaps)} haplotype-specific gaps")

        # Use first normalized haplotype as reference for read alignment
        ref_hap = self.hap_names[0]
        ref_assembly = self.normalized_assemblies[ref_hap]

        # Step 0b: Prepare Hi-C data if available (use normalized assembly)
        if self.hic_reads or self.hic_bam:
            self.logger.info("STEP 0b: Preparing Hi-C data")
            self._prepare_hic_data(ref_assembly)

        # NOTE: For polyploid mode, reads filtering is done AFTER phasing, not before.
        # Phasing requires full reads to ensure adequate SNP coverage.
        # The HaploidEngine in STEP 3 will filter phased reads during gap filling.

        # Step 1: Detect haplotype-specific SNPs using alignment
        snp_file = self.output_dir / "snp_database.json"
        if _is_done(self.output_dir, "snp_detection") and snp_file.exists():
            self.logger.info("STEP 1: Reusing SNP database")
            snp_db = self._load_snp_db(snp_file)
        else:
            self.logger.info("STEP 1: Detecting haplotype-specific SNPs (alignment-based)")
            snp_db = phaser.detect_haplotype_snps(method=self.phasing_method)
            self._save_snp_db(snp_db, snp_file)
            _mark_done(self.output_dir, "snp_detection")

        # Step 2: Phase reads (check for existing phased reads)
        phased_hifi = {}
        phased_ont = {}

        if _is_done(self.output_dir, "phasing"):
            self.logger.info("STEP 2: Reusing phased reads")
            # Reconstruct phased_hifi and phased_ont from files
            for hap_name in self.hap_names:
                hifi_path = self.output_dir / f"phased_{hap_name}_hifi.fasta"
                ont_path = self.output_dir / f"phased_{hap_name}_ont.fasta"
                if hifi_path.exists():
                    phased_hifi[hap_name] = hifi_path
                if ont_path.exists():
                    phased_ont[hap_name] = ont_path
        else:
            # Collect gap regions for filtering
            gap_regions_for_filter = {}
            for hap_name in self.hap_names:
                gaps = phaser.gap_regions.get(hap_name, [])
                gap_regions_for_filter[hap_name] = [
                    {'chrom': g.chrom, 'start': g.start, 'end': g.end}
                    for g in gaps
                ]

            # Step 2a: Align HiFi reads
            hifi_bam = None
            filtered_hifi_names = None
            if self.hifi_reads:
                self.logger.info("STEP 2a: Aligning HiFi reads")

                hifi_bam = self.output_dir / "hifi_aligned.bam"
                if not (self.resume and hifi_bam.exists() and hifi_bam.stat().st_size > 0):
                    self._align_reads_to_ref(self.hifi_reads, ref_assembly, hifi_bam, 'map-hifi')
                else:
                    self.logger.info(f"  Reusing existing BAM: {hifi_bam}")

                # Filter reads based on gaps (only if optimized mode)
                if self.optimized_mode:
                    self.logger.info("  Filtering HiFi reads based on gap proximity...")
                    filtered_hifi_names = self._filter_reads_by_gaps(
                        hifi_bam, gap_regions_for_filter, self.min_mapq
                    )

            # Step 2b: Align ONT reads
            ont_bam = None
            filtered_ont_names = None
            if self.ont_reads:
                self.logger.info("STEP 2b: Aligning ONT reads")

                ont_bam = self.output_dir / "ont_aligned.bam"
                if not (self.resume and ont_bam.exists() and ont_bam.stat().st_size > 0):
                    self._align_reads_to_ref(self.ont_reads, ref_assembly, ont_bam, 'map-ont')
                else:
                    self.logger.info(f"  Reusing existing BAM: {ont_bam}")

                # Filter reads based on gaps (only if optimized mode)
                if self.optimized_mode:
                    self.logger.info("  Filtering ONT reads based on gap proximity...")
                    filtered_ont_names = self._filter_reads_by_gaps(
                        ont_bam, gap_regions_for_filter, self.min_mapq
                    )

            # Step 2c: Phase HiFi reads
            if hifi_bam:
                self.logger.info("STEP 2c: Phasing HiFi reads")
                phased_hifi = phaser.phase_reads_from_bam(
                    hifi_bam, snp_db,
                    self.output_dir / "phased",
                    read_type='hifi',
                    filtered_read_names=filtered_hifi_names  # Only output filtered reads
                )

            # Step 2d: Phase ONT reads
            if ont_bam:
                self.logger.info("STEP 2d: Phasing ONT reads")
                phased_ont = phaser.phase_reads_from_bam(
                    ont_bam, snp_db,
                    self.output_dir / "phased",
                    read_type='ont',
                    filtered_read_names=filtered_ont_names  # Only output filtered reads
                )

            # Save phased reads to checkpoint
            phased_paths = {}
            for hap in self.hap_names:
                if hap in phased_hifi:
                    phased_paths[f"{hap}_hifi"] = str(phased_hifi[hap])
                if hap in phased_ont:
                    phased_paths[f"{hap}_ont"] = str(phased_ont[hap])
            if 'ambiguous' in phased_hifi:
                phased_paths['ambiguous_hifi'] = str(phased_hifi['ambiguous'])
            if 'ambiguous' in phased_ont:
                phased_paths['ambiguous_ont'] = str(phased_ont['ambiguous'])

            _mark_done(self.output_dir, "phasing")

        # Step 2c: Enhance phasing with Hi-C (if available)
        if self.hic_analyzer:
            self.logger.info("STEP 2c: Enhancing phasing with Hi-C long-range information")

            # Combine phased reads info
            all_phased = {}
            if phased_hifi:
                all_phased.update({k: v for k, v in phased_hifi.items() if k != 'ambiguous'})
            if phased_ont:
                for k, v in phased_ont.items():
                    if k != 'ambiguous' and k not in all_phased:
                        all_phased[k] = v

            if all_phased:
                enhanced_assignments = self.hic_analyzer.enhance_phasing(
                    snp_db, all_phased, self.hap_names
                )

                # Write enhanced phased reads
                rescued_counts = {hap: 0 for hap in self.hap_names}
                for read_name, hap in enhanced_assignments.items():
                    if hap in rescued_counts:
                        rescued_counts[hap] += 1

                self.logger.info(f"  Hi-C rescued reads per haplotype: {rescued_counts}")

        # Step 3: Run gap filling for each haplotype with both read types
        # Use NORMALIZED assemblies (skip_normalization=True in HaploidEngine)
        self.logger.info("STEP 3: Running gap filling for each haplotype")

        filled_assemblies = {}

        for i, hap_name in enumerate(self.hap_names):
            hap_file = self.normalized_assemblies[hap_name]

            self.logger.info(f"\n{'='*40}")
            self.logger.info(f"Processing {hap_name}: {hap_file.name}")
            self.logger.info(f"{'='*40}")

            # Prepare HiFi reads for this haplotype
            hap_hifi_reads = None
            if phased_hifi:
                hap_hifi = phased_hifi.get(hap_name)
                if hap_hifi and hap_hifi.exists() and hap_hifi.stat().st_size > 0:
                    if self.use_ambiguous_reads and 'ambiguous' in phased_hifi:
                        # Combine phased + ambiguous HiFi reads
                        combined_hifi = self.output_dir / f"{hap_name}_combined_hifi.fasta"
                        self._combine_fasta_files(
                            [hap_hifi, phased_hifi['ambiguous']],
                            combined_hifi
                        )
                        hap_hifi_reads = combined_hifi
                    else:
                        hap_hifi_reads = hap_hifi

                    self.logger.info(f"  HiFi reads: {hap_hifi_reads}")

            # Prepare ONT reads for this haplotype
            hap_ont_reads = None
            if phased_ont:
                hap_ont = phased_ont.get(hap_name)
                if hap_ont and hap_ont.exists() and hap_ont.stat().st_size > 0:
                    if self.use_ambiguous_reads and 'ambiguous' in phased_ont:
                        # Combine phased + ambiguous ONT reads
                        combined_ont = self.output_dir / f"{hap_name}_combined_ont.fasta"
                        self._combine_fasta_files(
                            [hap_ont, phased_ont['ambiguous']],
                            combined_ont
                        )
                        hap_ont_reads = combined_ont
                    else:
                        hap_ont_reads = hap_ont

                    self.logger.info(f"  ONT reads: {hap_ont_reads}")

            # Check we have at least one read type
            if not hap_hifi_reads and not hap_ont_reads:
                self.logger.warning(f"  No phased reads for {hap_name}, skipping")
                filled_assemblies[hap_name] = hap_file
                continue

            hap_output = self.output_dir / hap_name

            filled_assembly = self._run_gap_filling_for_haplotype(
                hap_file,
                hap_hifi_reads,
                hap_ont_reads,
                hap_output,
                hap_name=hap_name,
                skip_normalization=True  # Already normalized
            )

            filled_assemblies[hap_name] = filled_assembly

        # Step 4: Generate summary
        self.logger.info("\nSTEP 4: Generating summary report")
        self._generate_summary(filled_assemblies, phaser)

        # Mark as complete
        _mark_done(self.output_dir, "complete")
        self.logger.info("Polyploid gap filling complete")

        return filled_assemblies

    def _load_snp_db(self, snp_file: Path) -> Dict:
        """Load SNP database from JSON file"""
        with open(snp_file) as f:
            data = json.load(f)

        # Convert string keys back to integers
        snp_db = {}
        for chrom, positions in data.items():
            snp_db[chrom] = {int(pos): bases for pos, bases in positions.items()}
        return snp_db

    def _prepare_hic_data(self, ref_assembly: Path):
        """
        Prepare Hi-C BAM files for all haplotypes using merged reference alignment.

        Strategy:
        1. Create merged reference with all haplotypes (hap1__Chr1, hap2__Chr1, ...)
        2. Align Hi-C reads once to merged reference
        3. Split BAM by haplotype prefix
        4. Create per-haplotype HiCAnalyzer
        """
        # Store per-haplotype analyzers
        self.hic_analyzers: Dict[str, HiCAnalyzer] = {}

        # Check if user provided pre-aligned BAM
        if self.hic_bam:
            # User provided BAM - assume it's aligned to first haplotype (backward compatible)
            self.logger.info("  Using provided Hi-C BAM (aligned to single haplotype)")
            self.hic_analyzer = HiCAnalyzer(
                hic_bam=str(self.hic_bam),
                assembly_file=str(ref_assembly),
                threads=self.threads
            )
            # Use same analyzer for all haplotypes
            for hap_name in self.hap_names:
                self.hic_analyzers[hap_name] = self.hic_analyzer
            return

        if not self.hic_reads:
            self.logger.warning("  No Hi-C data provided")
            return

        # Step 1: Create merged reference
        self.logger.info("  Creating merged reference for Hi-C alignment...")
        merged_ref = self._create_merged_hic_reference()

        # Step 2: Align Hi-C to merged reference
        merged_hic_bam = self.output_dir / "hic_merged_aligned.bam"
        if not merged_hic_bam.exists():
            self.logger.info("  Aligning Hi-C reads to merged reference...")
            align_hic_reads(
                self.hic_reads[0],
                self.hic_reads[1],
                str(merged_ref),
                str(merged_hic_bam),
                threads=self.threads
            )
        else:
            self.logger.info(f"  Reusing existing Hi-C BAM: {merged_hic_bam}")

        # Step 3: Split BAM by haplotype
        self.logger.info("  Splitting Hi-C BAM by haplotype...")
        for hap_name in self.hap_names:
            hap_hic_bam = self.output_dir / f"hic_{hap_name}.bam"
            hap_assembly = self.normalized_assemblies[hap_name]

            if not hap_hic_bam.exists():
                self._split_hic_bam_by_haplotype(merged_hic_bam, hap_hic_bam, hap_name)

            # Step 4: Create analyzer for this haplotype
            if hap_hic_bam.exists() and hap_hic_bam.stat().st_size > 0:
                self.hic_analyzers[hap_name] = HiCAnalyzer(
                    hic_bam=str(hap_hic_bam),
                    assembly_file=str(hap_assembly),
                    threads=self.threads
                )
                self.logger.info(f"    {hap_name}: Hi-C analyzer ready")

        # Set default analyzer to hap1 for backward compatibility
        if self.hap_names[0] in self.hic_analyzers:
            self.hic_analyzer = self.hic_analyzers[self.hap_names[0]]

        self.logger.info(f"  Hi-C analyzers ready for {len(self.hic_analyzers)} haplotypes")

    def _create_merged_hic_reference(self) -> Path:
        """Create merged reference with haplotype prefixes for Hi-C alignment"""
        merged_ref = self.output_dir / "hic_merged_reference.fasta"

        if merged_ref.exists():
            return merged_ref

        with open(merged_ref, 'w') as out:
            for hap_name in self.hap_names:
                hap_file = self.normalized_assemblies[hap_name]
                for record in SeqIO.parse(hap_file, 'fasta'):
                    # Add haplotype prefix: hap1__Chr1
                    new_id = f"{hap_name}__{record.id}"
                    out.write(f">{new_id}\n")
                    seq = str(record.seq)
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i+80] + '\n')

        self.logger.info(f"    Created merged reference: {merged_ref}")
        return merged_ref

    def _split_hic_bam_by_haplotype(self, input_bam: Path, output_bam: Path, hap_name: str):
        """
        Split Hi-C BAM to only include alignments to specific haplotype.
        Also removes haplotype prefix from contig names.
        """
        prefix = f"{hap_name}__"

        with pysam.AlignmentFile(str(input_bam), 'rb') as inp:
            # Find contigs for this haplotype
            contigs = [ref for ref in inp.references if ref.startswith(prefix)]

            if not contigs:
                self.logger.warning(f"    No contigs found for {hap_name}")
                # Create empty BAM
                with pysam.AlignmentFile(str(output_bam), 'wb', template=inp) as out:
                    pass
                return

            # Create new header without prefix
            new_header = inp.header.to_dict()
            new_refs = []
            for sq in new_header.get('SQ', []):
                if sq['SN'].startswith(prefix):
                    sq['SN'] = sq['SN'][len(prefix):]
                    new_refs.append(sq)
            new_header['SQ'] = new_refs

            with pysam.AlignmentFile(str(output_bam), 'wb', header=new_header) as out:
                for contig in contigs:
                    for read in inp.fetch(contig):
                        # Create new alignment with renamed reference
                        a = pysam.AlignedSegment()
                        a.query_name = read.query_name
                        a.query_sequence = read.query_sequence
                        a.flag = read.flag
                        a.reference_id = out.get_tid(contig[len(prefix):])
                        a.reference_start = read.reference_start
                        a.mapping_quality = read.mapping_quality
                        a.cigar = read.cigar
                        a.query_qualities = read.query_qualities

                        # Handle mate information
                        if read.is_paired and not read.mate_is_unmapped:
                            mate_ref = read.next_reference_name
                            if mate_ref and mate_ref.startswith(prefix):
                                a.next_reference_id = out.get_tid(mate_ref[len(prefix):])
                                a.next_reference_start = read.next_reference_start
                            else:
                                # Mate on different haplotype - mark as unmapped mate
                                a.flag |= 0x8  # mate unmapped
                                a.next_reference_id = -1
                                a.next_reference_start = 0

                        if a.reference_id >= 0:
                            out.write(a)

        # Index the output BAM
        pysam.index(str(output_bam))
        self.logger.info(f"    Split BAM for {hap_name}: {output_bam}")

    def _align_reads_to_ref(self, reads_file: Path, ref_file: Path,
                            output_bam: Path, preset: str) -> bool:
        """Align reads to reference for phasing"""
        try:
            # --secondary=no: ensure unique assignment for phasing
            # -m 2G: memory limit per thread for samtools sort
            cmd = f"minimap2 -ax {preset} -t {self.threads} --secondary=no " \
                  f"{ref_file} {reads_file} | " \
                  f"samtools sort -@ {self.threads} -m 2G -o {output_bam} - && " \
                  f"samtools index {output_bam}"

            subprocess.run(cmd, shell=True, check=True,
                         capture_output=True, text=True)

            self.logger.info(f"  Created {output_bam}")
            return True

        except Exception as e:
            self.logger.error(f"Alignment error: {e}")
            return False

    def _save_snp_db(self, snp_db: Dict, output_file: Path):
        """Save SNP database to JSON"""
        serializable = {}
        for chrom, positions in snp_db.items():
            serializable[chrom] = {str(pos): bases for pos, bases in positions.items()}

        with open(output_file, 'w') as f:
            json.dump(serializable, f, indent=2)

    def _filter_reads_by_gaps(self, bam_file: Path, gap_regions: Dict[str, List],
                               min_mapq: int = 20, gap_proximity: int = 1000) -> Set[str]:
        """
        Filter reads based on proximity to gaps from ALL haplotypes.

        Only filter out reads that are:
        1. Anchored (high mapq, no significant soft-clips)
        2. Far from ALL gaps across all haplotypes

        Args:
            bam_file: BAM file with aligned reads
            gap_regions: Dict mapping hap_name -> list of gap dicts
            min_mapq: Minimum mapping quality for "anchored" reads
            gap_proximity: Distance threshold for "near gap"

        Returns:
            Set of read names that should be kept
        """
        # Collect all gap intervals from all haplotypes
        # Use reference (hap1) coordinates since reads are aligned to hap1
        all_gap_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

        # Only use hap1's gaps since reads are aligned to hap1
        ref_hap = self.hap_names[0]
        if ref_hap in gap_regions:
            for gap in gap_regions[ref_hap]:
                chrom = gap['chrom']
                all_gap_intervals[chrom].append((gap['start'], gap['end']))

        # Sort intervals for each chromosome
        for chrom in all_gap_intervals:
            all_gap_intervals[chrom].sort()

        self.logger.info(f"  Filtering based on {sum(len(v) for v in all_gap_intervals.values())} gap regions")

        kept_reads = set()
        stats = {'total': 0, 'kept': 0, 'filtered': 0}

        with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary:
                    continue

                stats['total'] += 1

                # Unmapped reads: keep (not anchored)
                if read.is_unmapped:
                    kept_reads.add(read.query_name)
                    stats['kept'] += 1
                    continue

                # Check if read is anchored (high mapq + no significant clips)
                is_anchored = self._is_read_anchored(read, min_mapq)

                if not is_anchored:
                    # Not anchored: keep
                    kept_reads.add(read.query_name)
                    stats['kept'] += 1
                    continue

                # Anchored read: check if near any gap
                chrom = read.reference_name
                read_start = read.reference_start
                read_end = read.reference_end

                if read_start is None or read_end is None:
                    kept_reads.add(read.query_name)
                    stats['kept'] += 1
                    continue

                # Check distance to gaps
                near_gap = False
                if chrom in all_gap_intervals:
                    for gap_start, gap_end in all_gap_intervals[chrom]:
                        # Calculate distance
                        if read_end <= gap_start:
                            distance = gap_start - read_end
                        elif read_start >= gap_end:
                            distance = read_start - gap_end
                        else:
                            distance = 0  # Overlaps

                        if distance <= gap_proximity:
                            near_gap = True
                            break

                if near_gap:
                    kept_reads.add(read.query_name)
                    stats['kept'] += 1
                else:
                    stats['filtered'] += 1

        self.logger.info(f"  Reads filtering: {stats['total']:,} total → "
                        f"{stats['kept']:,} kept ({100*stats['kept']/max(1,stats['total']):.1f}%), "
                        f"{stats['filtered']:,} filtered")

        return kept_reads

    def _is_read_anchored(self, read: pysam.AlignedSegment, min_mapq: int = 20) -> bool:
        """Check if a read is anchored (high mapq + no significant soft-clips)"""
        if read.mapping_quality < min_mapq:
            return False

        cigar = read.cigartuples
        if not cigar:
            return False

        CLIP_THRESHOLD = 100

        # Check for significant clips
        if cigar[0][0] in (4, 5) and cigar[0][1] >= CLIP_THRESHOLD:
            return False
        if cigar[-1][0] in (4, 5) and cigar[-1][1] >= CLIP_THRESHOLD:
            return False

        return True

    def _combine_fasta_files(self, input_files: List[Path], output_file: Path):
        """Combine multiple FASTA files"""
        with open(output_file, 'w') as out:
            for input_file in input_files:
                if input_file and input_file.exists():
                    with open(input_file) as inp:
                        out.write(inp.read())

    def _run_gap_filling_for_haplotype(self, assembly: Path,
                                        hifi_reads: Optional[Path],
                                        ont_reads: Optional[Path],
                                        output_dir: Path,
                                        hap_name: str,
                                        skip_normalization: bool = False) -> Path:
        """Run gap filling for a single haplotype with both read types"""

        # Get Hi-C BAM for this haplotype (if available)
        hic_bam = None
        if hasattr(self, 'hic_analyzers') and hap_name in self.hic_analyzers:
            hic_bam = self.output_dir / f"hic_{hap_name}.bam"
            if not hic_bam.exists():
                hic_bam = None

        self.logger.info(f"  Starting HaploidEngine with:")
        self.logger.info(f"    Assembly: {assembly}")
        self.logger.info(f"    HiFi: {hifi_reads}")
        self.logger.info(f"    ONT: {ont_reads}")
        self.logger.info(f"    Hi-C BAM: {hic_bam}")
        self.logger.info(f"    Skip normalization: {skip_normalization}")

        # Disable reads filtering in HaploidEngine because:
        # - Reads are already filtered at polyploid STEP 2 (based on BAM alignment)
        # - No need to filter again, saves time
        # - Keep parallel filling enabled for speedup
        engine = HaploidEngine(
            assembly_file=str(assembly),
            hifi_reads=str(hifi_reads) if hifi_reads else None,
            ont_reads=str(ont_reads) if ont_reads else None,
            hic_bam=str(hic_bam) if hic_bam else None,
            output_dir=str(output_dir),
            threads=self.threads,
            max_iterations=self.max_iterations,
            skip_normalization=skip_normalization,
            optimized_mode=False,  # Already filtered at polyploid level
            parallel_filling=True  # Parallel gap filling
        )

        return engine.run()

    def _generate_summary(self, filled_assemblies: Dict[str, Path], phaser: ReadPhaser):
        """Generate summary report"""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'num_haplotypes': self.num_haplotypes,
            'phasing_method': self.phasing_method,
            'data_types': {
                'hifi': self.hifi_reads is not None,
                'ont': self.ont_reads is not None
            },
            'gap_info': {
                hap_name: {
                    'total_gaps': len(phaser.gap_regions.get(hap_name, [])),
                    'haplotype_specific_gaps': len(phaser.haplotype_specific_gaps.get(hap_name, set()))
                }
                for hap_name in self.hap_names
            },
            'haplotypes': {}
        }

        for hap_name, assembly_path in filled_assemblies.items():
            stats_file = assembly_path.parent / "final_stats.json"
            if stats_file.exists():
                with open(stats_file) as f:
                    stats = json.load(f)
                summary['haplotypes'][hap_name] = stats
            else:
                summary['haplotypes'][hap_name] = {'assembly': str(assembly_path)}

        summary_file = self.output_dir / "polyploid_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        self.logger.info("\n" + "=" * 60)
        self.logger.info("POLYPLOID GAP FILLING COMPLETE")
        self.logger.info("=" * 60)
        self.logger.info(f"Data types used: HiFi={self.hifi_reads is not None}, ONT={self.ont_reads is not None}")
        for hap_name, assembly_path in filled_assemblies.items():
            self.logger.info(f"  {hap_name}: {assembly_path}")


# Backwards compatibility
PolyploidGapFiller = PolyploidEngine
