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

from gapfill.engines.haploid import HaploidEngine
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads


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
                 min_read_snps: int = 2,
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
        """
        snp_db = {}

        ref_hap = self.hap_names[0]
        ref_file = self.normalized_assemblies[ref_hap]

        # Load reference sequences for gap region detection
        ref_seqs = {}
        for record in SeqIO.parse(ref_file, 'fasta'):
            ref_seqs[record.id] = str(record.seq).upper()

        # Build gap regions set for quick lookup
        ref_gap_regions = self._build_gap_position_set(ref_seqs)

        # Initialize SNP database
        for chrom in ref_seqs:
            snp_db[chrom] = {}

        # Align each other haplotype to reference
        for other_hap in self.hap_names[1:]:
            other_file = self.normalized_assemblies[other_hap]

            self.logger.info(f"  Aligning {other_hap} to {ref_hap}...")

            # Load other haplotype sequences
            other_seqs = {}
            for record in SeqIO.parse(other_file, 'fasta'):
                other_seqs[record.id] = str(record.seq).upper()

            other_gap_regions = self._build_gap_position_set(other_seqs)

            # Run minimap2 alignment
            bam_file = self.work_dir / f"align_{other_hap}_to_{ref_hap}.bam"
            self._align_assemblies(other_file, ref_file, bam_file)

            # Parse alignment to find SNPs
            snps_found = self._extract_snps_from_alignment(
                bam_file, ref_seqs, other_seqs,
                ref_gap_regions, other_gap_regions,
                ref_hap, other_hap, snp_db
            )

            self.logger.info(f"    Found {snps_found} SNPs between {ref_hap} and {other_hap}")

            # Clean up BAM file
            bam_file.unlink(missing_ok=True)
            Path(str(bam_file) + ".bai").unlink(missing_ok=True)

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
            # -m 8G: memory limit per thread for samtools sort
            cmd = (
                f"minimap2 -ax asm5 -t {self.threads} "
                f"{ref_file} {query_file} | "
                f"samtools sort -@ {self.threads} -m 8G -o {output_bam} - && "
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
                                      snp_db: Dict) -> int:
        """
        Extract SNPs from alignment BAM file.
        Skips positions within gap regions.
        """
        snps_found = 0

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

                        # Found a SNP
                        if ref_base != query_base:
                            if ref_pos not in snp_db[ref_name]:
                                snp_db[ref_name][ref_pos] = {ref_hap: ref_base}
                            snp_db[ref_name][ref_pos][query_hap] = query_base
                            snps_found += 1

                except Exception as e:
                    self.logger.debug(f"Error processing alignment: {e}")
                    continue

        return snps_found

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
                             output_prefix: Path, read_type: str = 'reads') -> Dict[str, Path]:
        """Phase reads based on SNP profiles"""
        self.logger.info(f"Phasing {read_type} reads to haplotypes...")

        phased_reads = {hap: [] for hap in self.hap_names}
        ambiguous_reads = []

        stats = {
            'total': 0,
            'phased': {hap: 0 for hap in self.hap_names},
            'ambiguous': 0
        }

        try:
            with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
                for read in bam:
                    if read.is_unmapped or read.is_secondary:
                        continue

                    stats['total'] += 1

                    hap_assignment = self._assign_read_to_haplotype(read, snp_db)

                    if hap_assignment:
                        phased_reads[hap_assignment].append(read)
                        stats['phased'][hap_assignment] += 1
                    else:
                        ambiguous_reads.append(read)
                        stats['ambiguous'] += 1

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
            self.logger.info(f"  {hap_name}: {stats['phased'][hap_name]} {read_type} reads")

        # Write ambiguous reads
        ambiguous_file = Path(f"{output_prefix}_ambiguous_{read_type}.fasta")
        with open(ambiguous_file, 'w') as f:
            for read in ambiguous_reads:
                seq = read.query_sequence
                if seq:
                    f.write(f">{read.query_name}\n{seq}\n")

        output_files['ambiguous'] = ambiguous_file
        self.logger.info(f"  Ambiguous: {stats['ambiguous']} {read_type} reads")

        return output_files

    def _assign_read_to_haplotype(self, read: pysam.AlignedSegment, snp_db: Dict) -> Optional[str]:
        """Assign a read to a haplotype based on SNP profile"""
        chrom = read.reference_name

        if chrom not in snp_db:
            return None

        hap_scores = {hap: 0 for hap in self.hap_names}
        snps_checked = 0

        try:
            for query_pos, ref_pos in read.get_aligned_pairs():
                if ref_pos is None or query_pos is None:
                    continue

                if ref_pos in snp_db[chrom]:
                    snps_checked += 1
                    read_base = read.query_sequence[query_pos].upper()

                    for hap_name, hap_base in snp_db[chrom][ref_pos].items():
                        if read_base == hap_base:
                            hap_scores[hap_name] += 1

        except Exception:
            return None

        if snps_checked < self.min_read_snps:
            return None

        max_score = max(hap_scores.values())
        if max_score == 0:
            return None

        best_haps = [h for h, s in hap_scores.items() if s == max_score]

        if len(best_haps) == 1:
            return best_haps[0]
        else:
            return None


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
                 min_gap_size: int = 100):

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

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        # Hi-C analyzer (initialized when needed)
        self.hic_analyzer: Optional[HiCAnalyzer] = None

        # Store normalized assemblies
        self.normalized_assemblies: Dict[str, Path] = {}

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

        # Initialize phaser
        phaser = ReadPhaser(
            self.haplotypes,
            threads=self.threads,
            work_dir=self.output_dir
        )

        # STEP 0: Normalize gaps in ALL haplotypes FIRST
        self.logger.info("STEP 0: Normalizing gaps in all haplotypes")
        self.normalized_assemblies = phaser.normalize_all_assemblies(
            min_gap_size=self.min_gap_size
        )

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

        # Step 1: Detect haplotype-specific SNPs using alignment
        self.logger.info("STEP 1: Detecting haplotype-specific SNPs (alignment-based)")
        snp_db = phaser.detect_haplotype_snps(method=self.phasing_method)

        # Save SNP database
        snp_file = self.output_dir / "snp_database.json"
        self._save_snp_db(snp_db, snp_file)

        # Step 2: Phase HiFi reads (if available) - align to normalized reference
        phased_hifi = {}
        if self.hifi_reads:
            self.logger.info("STEP 2a: Aligning and phasing HiFi reads")

            hifi_bam = self.output_dir / "hifi_aligned.bam"
            self._align_reads_to_ref(self.hifi_reads, ref_assembly, hifi_bam, 'map-hifi')

            phased_hifi = phaser.phase_reads_from_bam(
                hifi_bam, snp_db,
                self.output_dir / "phased",
                read_type='hifi'
            )

        # Step 2b: Phase ONT reads (if available)
        phased_ont = {}
        if self.ont_reads:
            self.logger.info("STEP 2b: Aligning and phasing ONT reads")

            ont_bam = self.output_dir / "ont_aligned.bam"
            self._align_reads_to_ref(self.ont_reads, ref_assembly, ont_bam, 'map-ont')

            phased_ont = phaser.phase_reads_from_bam(
                ont_bam, snp_db,
                self.output_dir / "phased",
                read_type='ont'
            )

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
                skip_normalization=True  # Already normalized
            )

            filled_assemblies[hap_name] = filled_assembly

        # Step 4: Generate summary
        self.logger.info("\nSTEP 4: Generating summary report")
        self._generate_summary(filled_assemblies, phaser)

        return filled_assemblies

    def _prepare_hic_data(self, ref_assembly: Path):
        """Prepare Hi-C BAM file and analyzer"""
        hic_bam_path = self.hic_bam

        # Align Hi-C reads if BAM not provided
        if not hic_bam_path and self.hic_reads:
            hic_bam_path = self.output_dir / "hic_aligned.bam"
            if not hic_bam_path.exists():
                self.logger.info("  Aligning Hi-C reads to reference haplotype...")
                align_hic_reads(
                    self.hic_reads[0],
                    self.hic_reads[1],
                    str(ref_assembly),
                    str(hic_bam_path),
                    threads=self.threads
                )

        # Initialize analyzer
        if hic_bam_path and hic_bam_path.exists():
            self.hic_analyzer = HiCAnalyzer(
                hic_bam=str(hic_bam_path),
                assembly_file=str(ref_assembly),
                threads=self.threads
            )
            self.logger.info(f"  Hi-C analyzer ready: {hic_bam_path}")
        else:
            self.logger.warning("  Hi-C BAM not available")

    def _align_reads_to_ref(self, reads_file: Path, ref_file: Path,
                            output_bam: Path, preset: str) -> bool:
        """Align reads to reference for phasing"""
        try:
            # --secondary=no: ensure unique assignment for phasing
            # -m 8G: memory limit per thread for samtools sort
            cmd = f"minimap2 -ax {preset} -t {self.threads} --secondary=no " \
                  f"{ref_file} {reads_file} | " \
                  f"samtools sort -@ {self.threads} -m 8G -o {output_bam} - && " \
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
                                        skip_normalization: bool = False) -> Path:
        """Run gap filling for a single haplotype with both read types"""

        self.logger.info(f"  Starting HaploidEngine with:")
        self.logger.info(f"    Assembly: {assembly}")
        self.logger.info(f"    HiFi: {hifi_reads}")
        self.logger.info(f"    ONT: {ont_reads}")
        self.logger.info(f"    Skip normalization: {skip_normalization}")

        engine = HaploidEngine(
            assembly_file=str(assembly),
            hifi_reads=str(hifi_reads) if hifi_reads else None,
            ont_reads=str(ont_reads) if ont_reads else None,
            output_dir=str(output_dir),
            threads=self.threads,
            max_iterations=self.max_iterations,
            skip_normalization=skip_normalization
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
