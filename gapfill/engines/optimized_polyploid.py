#!/usr/bin/env python3
"""
Optimized Polyploid Gap Filling Engine

Key optimization: Batch alignment per iteration
- Merge all haplotype assemblies into one reference
- Align once per iteration (not once per haplotype)
- Split alignments by contig prefix to each haplotype

Comparison (4-ploid, 10 iterations):
  Original:  4 haps × 10 iters × 2 types = 80 alignments
  Optimized: 10 iters × 2 types = 20 alignments (75% reduction)

Why this works:
- Phased reads have SNP signatures of their assigned haplotype
- They align best to their correct haplotype in merged reference
- Coverage is NOT diluted because reads self-select correct target

Updates:
- Now uses alignment-based SNP detection from ReadPhaser
- Normalizes gaps before SNP detection
- Supports haplotype-specific gaps
"""

import logging
import subprocess
import json
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from collections import defaultdict

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gapfill.core.filler import GapFiller
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads
# Import ReadPhaser from polyploid module (with alignment-based SNP detection)
from gapfill.engines.polyploid import ReadPhaser
from gapfill.engines.haploid import _mark_done, _is_done, _clear_marker
from gapfill.utils.checkpoint import PolyploidCheckpointManager, CheckpointState
from gapfill.utils.reads_cache import ReadsCache


class OptimizedPolyploidEngine:
    """
    Optimized polyploid gap filler with batch alignment.

    Key optimization: Merge all haplotypes and align once per iteration.

    Workflow:
    1. Normalize gaps in ALL haplotypes FIRST
    2. Detect SNPs using alignment (handles different gap positions)
    3. Phase reads once
    4. Iterate: merge references → batch align → split by haplotype → fill gaps
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
                 min_gap_size: int = 100,
                 use_ambiguous_reads: bool = True,
                 min_read_snps: int = 1,
                 resume: bool = False,
                 clear_checkpoint: bool = False,
                 filter_reads: bool = True):

        self.haplotypes = [Path(h) for h in haplotype_assemblies]
        self.num_haplotypes = len(self.haplotypes)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.hic_reads = hic_reads
        self.hic_bam = Path(hic_bam) if hic_bam else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.max_iterations = max_iterations
        self.min_gap_size = min_gap_size
        self.use_ambiguous_reads = use_ambiguous_reads
        self.min_read_snps = min_read_snps
        self.resume = resume
        self.filter_reads = filter_reads

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        # Checkpoint manager
        self.checkpoint = PolyploidCheckpointManager(str(self.output_dir), self.hap_names)
        if clear_checkpoint:
            self.checkpoint.clear()

        self.hic_analyzer: Optional[HiCAnalyzer] = None
        self.hic_analyzers: Dict[str, HiCAnalyzer] = {}  # Per-haplotype analyzers

        # Store normalized assemblies
        self.normalized_assemblies: Dict[str, Path] = {}

        # Reads cache for filtering (filter once, use for phasing)
        self.reads_cache: Optional[ReadsCache] = None
        if self.filter_reads:
            self.reads_cache = ReadsCache(self.output_dir, threads=threads)

        self._validate_inputs()

        self.logger.info("=" * 60)
        self.logger.info("OptimizedPolyploidEngine initialized")
        self.logger.info(f"  Ploidy: {self.num_haplotypes}n")
        self.logger.info(f"  OPTIMIZATION: Batch alignment (1 align per iteration)")
        self.logger.info(f"  OPTIMIZATION: Reads filtering: {self.filter_reads}")
        self.logger.info(f"  Expected alignments: {self.max_iterations * 2} "
                        f"(vs {self.num_haplotypes * self.max_iterations * 2} original)")
        self.logger.info(f"  Resume: {self.resume}")
        self.logger.info("=" * 60)

    def _validate_inputs(self):
        for hap in self.haplotypes:
            if not hap.exists():
                raise FileNotFoundError(f"Haplotype not found: {hap}")
        if not self.hifi_reads and not self.ont_reads:
            raise ValueError("At least one reads file required")

    def run(self) -> Dict[str, Path]:
        """Run optimized polyploid gap filling with checkpoint support"""
        self.logger.info("Starting optimized polyploid gap filling...")

        # Check for existing checkpoint
        checkpoint_state = None
        if self.resume:
            if self.checkpoint.exists():
                checkpoint_state = self.checkpoint.load()
                if checkpoint_state:
                    self.logger.info("=" * 60)
                    self.logger.info("RESUMING FROM CHECKPOINT")
                    self.logger.info(f"  Phase: {checkpoint_state.phase}")
                    self.logger.info(f"  Iteration: {checkpoint_state.iteration}")
                    self.logger.info("=" * 60)
            else:
                # No checkpoint.json but --resume specified
                # Scan for existing files from previous run
                self.logger.info("=" * 60)
                self.logger.info("No checkpoint.json found, scanning existing files...")
                checkpoint_state = self.checkpoint.scan_existing_files()
                if checkpoint_state.phase != "init":
                    self.logger.info(f"  Detected phase: {checkpoint_state.phase}")
                    self.logger.info(f"  Detected iteration: {checkpoint_state.iteration}")
                    self.checkpoint.save(checkpoint_state)
                else:
                    self.logger.info("  No usable files found, starting fresh")
                    checkpoint_state = None
                self.logger.info("=" * 60)

        # Initialize checkpoint state if not resuming
        if not checkpoint_state:
            checkpoint_state = CheckpointState(
                engine="optimized_polyploid",
                phase="init",
                max_iterations=self.max_iterations
            )
            self.checkpoint.save(checkpoint_state)
            self.checkpoint.init_haplotype_states()

        # =====================================================================
        # PHASE 1: One-time setup (normalization + phasing)
        # =====================================================================
        self.logger.info("\n" + "=" * 60)
        self.logger.info("PHASE 1: Initial Setup (one-time)")
        self.logger.info("=" * 60)

        # Initialize phaser
        phaser = ReadPhaser(self.haplotypes, self.threads,
                            min_read_snps=self.min_read_snps, work_dir=self.output_dir)

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

        # Prepare Hi-C data if available (using merged reference)
        if self.hic_reads or self.hic_bam:
            self.logger.info("STEP 0b: Preparing Hi-C data (merged reference)")
            self._prepare_hic_data()

        # Step 0c: Filter reads (OPTIMIZATION - keep only gap-related reads)
        if self.filter_reads and self.reads_cache:
            self.logger.info("STEP 0c: Filtering reads (keeping gap-related reads only)...")

            # Collect all gaps from all haplotypes
            all_gaps = []
            for hap_name in self.hap_names:
                hap_gaps = phaser.gap_regions.get(hap_name, [])
                for gap in hap_gaps:
                    all_gaps.append({
                        'chrom': gap.chrom,
                        'start': gap.start,
                        'end': gap.end,
                        'size': gap.size,
                        'name': gap.name
                    })

            self.logger.info(f"  Total gaps across all haplotypes: {len(all_gaps)}")
            self.reads_cache.set_gap_regions(all_gaps)

            # Check for existing BAM files first, then filter
            # Priority: existing BAMs > preprocessing/full BAMs > new alignment
            preprocessing_dir = self.output_dir / "preprocessing"
            preprocessing_dir.mkdir(exist_ok=True)

            if self.hifi_reads and not self.reads_cache.is_cached('hifi'):
                self.logger.info("  Filtering HiFi reads...")
                # Check for existing BAM files (various naming conventions)
                existing_bams = [
                    self.output_dir / "hifi_aligned.bam",
                    self.output_dir / "phase_hifi.bam",
                    preprocessing_dir / "hifi_full.bam"
                ]
                source_bam = None
                for bam in existing_bams:
                    if bam.exists() and bam.stat().st_size > 0:
                        self.logger.info(f"    Reusing existing BAM: {bam}")
                        source_bam = bam
                        break

                if source_bam is None:
                    self.logger.info("    Aligning HiFi reads...")
                    source_bam = preprocessing_dir / "hifi_full.bam"
                    self._align_reads(self.hifi_reads, ref_assembly, source_bam, 'map-hifi')

                self.reads_cache.filter_bam(source_bam, 'hifi', min_mapq=20)
                cache_summary = self.reads_cache.get_summary()
                self.logger.info(f"  HiFi: kept {cache_summary['hifi']['stats']['kept']:,} / "
                               f"{cache_summary['hifi']['stats']['total']:,} reads "
                               f"({cache_summary['hifi']['stats']['kept_ratio']*100:.1f}%)")

            if self.ont_reads and not self.reads_cache.is_cached('ont'):
                self.logger.info("  Filtering ONT reads...")
                # Check for existing BAM files (various naming conventions)
                existing_bams = [
                    self.output_dir / "ont_aligned.bam",
                    self.output_dir / "phase_ont.bam",
                    preprocessing_dir / "ont_full.bam"
                ]
                source_bam = None
                for bam in existing_bams:
                    if bam.exists() and bam.stat().st_size > 0:
                        self.logger.info(f"    Reusing existing BAM: {bam}")
                        source_bam = bam
                        break

                if source_bam is None:
                    self.logger.info("    Aligning ONT reads...")
                    source_bam = preprocessing_dir / "ont_full.bam"
                    self._align_reads(self.ont_reads, ref_assembly, source_bam, 'map-ont')

                self.reads_cache.filter_bam(source_bam, 'ont', min_mapq=20)
                cache_summary = self.reads_cache.get_summary()
                self.logger.info(f"  ONT: kept {cache_summary['ont']['stats']['kept']:,} / "
                               f"{cache_summary['ont']['stats']['total']:,} reads "
                               f"({cache_summary['ont']['stats']['kept_ratio']*100:.1f}%)")

        # Detect SNPs using alignment-based method
        snp_file = self.output_dir / "snp_database.json"
        if _is_done(self.output_dir, "snp_detection") and snp_file.exists():
            self.logger.info("STEP 1: Reusing SNP database")
            snp_db = self._load_snp_db(snp_file)
        else:
            self.logger.info("STEP 1: Detecting haplotype-specific SNPs (alignment-based)")
            snp_db = phaser.detect_haplotype_snps()
            self._save_snp_db(snp_db, snp_file)
            _mark_done(self.output_dir, "snp_detection")

        # Phase reads (one-time, using normalized hap1 as reference)
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
            # Also load ambiguous reads
            ambiguous_hifi = self.output_dir / "phased_ambiguous_hifi.fasta"
            ambiguous_ont = self.output_dir / "phased_ambiguous_ont.fasta"
            if ambiguous_hifi.exists():
                phased_hifi['ambiguous'] = ambiguous_hifi
            if ambiguous_ont.exists():
                phased_ont['ambiguous'] = ambiguous_ont
        else:
            if self.hifi_reads:
                self.logger.info("STEP 2a: Phasing HiFi reads...")
                hifi_bam = self.output_dir / "phase_hifi.bam"
                # Check if BAM exists
                if not (self.resume and hifi_bam.exists() and hifi_bam.stat().st_size > 0):
                    # Use filtered reads if available
                    if self.filter_reads and self.reads_cache:
                        filtered_hifi = self.reads_cache.get_filtered_reads_path('hifi')
                        if filtered_hifi and filtered_hifi.exists():
                            self.logger.info("  Using filtered HiFi reads for phasing...")
                            self._align_reads(filtered_hifi, ref_assembly, hifi_bam, 'map-hifi')
                        else:
                            self._align_reads(self.hifi_reads, ref_assembly, hifi_bam, 'map-hifi')
                    else:
                        self._align_reads(self.hifi_reads, ref_assembly, hifi_bam, 'map-hifi')
                else:
                    self.logger.info(f"  Reusing existing BAM: {hifi_bam}")
                phased_hifi = phaser.phase_reads_from_bam(
                    hifi_bam, snp_db, self.output_dir / "phased", 'hifi'
                )

            if self.ont_reads:
                self.logger.info("STEP 2b: Phasing ONT reads...")
                ont_bam = self.output_dir / "phase_ont.bam"
                # Check if BAM exists
                if not (self.resume and ont_bam.exists() and ont_bam.stat().st_size > 0):
                    # Use filtered reads if available
                    if self.filter_reads and self.reads_cache:
                        filtered_ont = self.reads_cache.get_filtered_reads_path('ont')
                        if filtered_ont and filtered_ont.exists():
                            self.logger.info("  Using filtered ONT reads for phasing...")
                            self._align_reads(filtered_ont, ref_assembly, ont_bam, 'map-ont')
                        else:
                            self._align_reads(self.ont_reads, ref_assembly, ont_bam, 'map-ont')
                    else:
                        self._align_reads(self.ont_reads, ref_assembly, ont_bam, 'map-ont')
                else:
                    self.logger.info(f"  Reusing existing BAM: {ont_bam}")
                phased_ont = phaser.phase_reads_from_bam(
                    ont_bam, snp_db, self.output_dir / "phased", 'ont'
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

        # Prepare combined phased reads per haplotype
        combined_reads = self._prepare_combined_reads(phased_hifi, phased_ont)

        # =====================================================================
        # PHASE 2: Iterative gap filling with batch alignment
        # =====================================================================
        self.logger.info("\n" + "=" * 60)
        self.logger.info("PHASE 2: Iterative Gap Filling (batch alignment)")
        self.logger.info("=" * 60)

        # Current assemblies for each haplotype (start with normalized)
        current_assemblies = {
            hap: self.normalized_assemblies[hap]
            for hap in self.hap_names
        }

        # Track filled gaps (load from checkpoint if resuming)
        filled_gaps = {hap: self.checkpoint.get_hap_completed_gaps(hap) for hap in self.hap_names}
        failed_gaps = {hap: self.checkpoint.get_hap_failed_gaps(hap) for hap in self.hap_names}

        # Determine starting iteration
        start_iteration = 0
        if self.resume and checkpoint_state.iteration > 0:
            start_iteration = max(0, checkpoint_state.iteration - 1)
            self.logger.info(f"Resuming from iteration {start_iteration + 1}")

            # Load current assemblies from checkpoint if available
            for hap_name in self.hap_names:
                hap_state = checkpoint_state.haplotype_states.get(hap_name, {})
                if hap_state.get('current_assembly'):
                    hap_assembly = Path(hap_state['current_assembly'])
                    if hap_assembly.exists():
                        current_assemblies[hap_name] = hap_assembly

        for iteration in range(start_iteration + 1, self.max_iterations + 1):
            self.logger.info(f"\n--- Iteration {iteration} ---")

            # Update checkpoint
            checkpoint_state.iteration = iteration
            self.checkpoint.save(checkpoint_state)

            iter_dir = self.output_dir / f"iteration_{iteration}"
            iter_dir.mkdir(exist_ok=True)

            # Step 1: Create merged reference (all haplotypes with prefixes)
            merged_ref = self._create_merged_reference(current_assemblies, iter_dir)

            # Step 2: Merge all phased reads
            merged_hifi = self._merge_phased_reads(combined_reads, 'hifi', iter_dir)
            merged_ont = self._merge_phased_reads(combined_reads, 'ont', iter_dir)

            # Step 3: SINGLE alignment per read type (the key optimization!)
            self.logger.info("Batch alignment (all haplotypes at once)...")
            hifi_bam = None
            ont_bam = None

            if merged_hifi and merged_hifi.stat().st_size > 0:
                hifi_bam = iter_dir / "merged_hifi.bam"
                hifi_bai = iter_dir / "merged_hifi.bam.bai"
                # Check if BAM exists and is valid
                if self.resume and hifi_bam.exists() and hifi_bai.exists() and hifi_bam.stat().st_size > 0:
                    self.logger.info(f"  Reusing existing HiFi BAM: {hifi_bam}")
                else:
                    self._align_reads(merged_hifi, merged_ref, hifi_bam, 'map-hifi')

            if merged_ont and merged_ont.stat().st_size > 0:
                ont_bam = iter_dir / "merged_ont.bam"
                ont_bai = iter_dir / "merged_ont.bam.bai"
                # Check if BAM exists and is valid
                if self.resume and ont_bam.exists() and ont_bai.exists() and ont_bam.stat().st_size > 0:
                    self.logger.info(f"  Reusing existing ONT BAM: {ont_bam}")
                else:
                    self._align_reads(merged_ont, merged_ref, ont_bam, 'map-ont')

            # Step 4: Split BAM by haplotype and fill gaps
            new_fills = 0

            for hap_name in self.hap_names:
                hap_dir = iter_dir / hap_name
                hap_dir.mkdir(exist_ok=True)

                # Split alignments for this haplotype
                hap_hifi_bam = None
                hap_ont_bam = None

                if hifi_bam:
                    hap_hifi_bam = hap_dir / "hifi.bam"
                    self._split_bam_by_prefix(hifi_bam, hap_hifi_bam, hap_name)

                if ont_bam:
                    hap_ont_bam = hap_dir / "ont.bam"
                    self._split_bam_by_prefix(ont_bam, hap_ont_bam, hap_name)

                # Find gaps
                gaps = self._find_gaps(current_assemblies[hap_name])
                remaining_gaps = [
                    g for g in gaps
                    if g['name'] not in filled_gaps[hap_name]
                    and g['name'] not in failed_gaps[hap_name]
                ]

                if not remaining_gaps:
                    continue

                self.logger.info(f"  {hap_name}: {len(remaining_gaps)} gaps to fill")

                # Create GapFiller for this haplotype
                work_dir = hap_dir / "work"
                work_dir.mkdir(exist_ok=True)

                # Need to create a haplotype-specific assembly (without prefix)
                hap_assembly = self._extract_haplotype_assembly(
                    current_assemblies[hap_name], hap_name, hap_dir
                )

                # Get Hi-C analyzer for this haplotype (if available)
                hap_hic_analyzer = self.hic_analyzers.get(hap_name)

                filler = GapFiller(
                    assembly_file=str(hap_assembly),
                    hifi_bam=str(hap_hifi_bam) if hap_hifi_bam else None,
                    ont_bam=str(hap_ont_bam) if hap_ont_bam else None,
                    threads=self.threads,
                    work_dir=str(work_dir),
                    hic_analyzer=hap_hic_analyzer
                )

                # Fill gaps
                fill_results = {}
                for gap in remaining_gaps:
                    result = filler.fill_gap(gap)
                    fill_results[gap['name']] = {'gap': gap, 'result': result}

                    if result['success']:
                        if result.get('is_complete') and not result.get('has_placeholder'):
                            filled_gaps[hap_name].add(gap['name'])
                            new_fills += 1
                            # Save to checkpoint
                            self.checkpoint.add_completed_gap_for_hap(
                                hap_name, gap['name'], result.get('sequence', '')
                            )
                    else:
                        self.checkpoint.add_failed_gap_for_hap(hap_name, gap['name'])

                filler.close()

                # Apply fills
                if fill_results:
                    current_assemblies[hap_name] = self._apply_fills(
                        current_assemblies[hap_name], fill_results, hap_dir
                    )
                    # Update checkpoint with current assembly
                    self.checkpoint.set_hap_assembly(hap_name, str(current_assemblies[hap_name]))

            self.logger.info(f"  Iteration {iteration}: {new_fills} new complete fills")

            if new_fills == 0:
                self.logger.info("No progress, stopping early")
                break

        # =====================================================================
        # PHASE 3: Finalize
        # =====================================================================
        self.logger.info("\n" + "=" * 60)
        self.logger.info("PHASE 3: Finalizing")
        self.logger.info("=" * 60)

        final_assemblies = {}
        for hap_name, assembly in current_assemblies.items():
            final_path = self.output_dir / f"{hap_name}_filled.fasta"
            shutil.copy(assembly, final_path)
            final_assemblies[hap_name] = final_path
            self.logger.info(f"  {hap_name}: {final_path}")

        self._save_summary(final_assemblies, filled_gaps, failed_gaps, phaser)

        # Mark as complete
        _mark_done(self.output_dir, "complete")
        self.logger.info("Optimized polyploid gap filling complete")

        return final_assemblies

    def _load_snp_db(self, snp_file: Path) -> Dict:
        """Load SNP database from JSON file"""
        with open(snp_file) as f:
            data = json.load(f)

        # Convert string keys back to integers
        snp_db = {}
        for chrom, positions in data.items():
            snp_db[chrom] = {int(pos): bases for pos, bases in positions.items()}
        return snp_db

    def _save_snp_db(self, snp_db: Dict, output_file: Path):
        """Save SNP database to JSON"""
        serializable = {}
        for chrom, positions in snp_db.items():
            serializable[chrom] = {str(pos): bases for pos, bases in positions.items()}

        with open(output_file, 'w') as f:
            json.dump(serializable, f, indent=2)

    def _create_merged_reference(self, assemblies: Dict[str, Path],
                                  output_dir: Path) -> Path:
        """Create merged reference with haplotype prefixes"""
        merged_ref = output_dir / "merged_reference.fasta"

        with open(merged_ref, 'w') as out:
            for hap_name, assembly in assemblies.items():
                for record in SeqIO.parse(assembly, 'fasta'):
                    # Add haplotype prefix to contig name
                    new_id = f"{hap_name}__{record.id}"
                    out.write(f">{new_id}\n")
                    seq = str(record.seq)
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i+80] + '\n')

        self.logger.info(f"  Merged reference: {merged_ref}")
        return merged_ref

    def _merge_phased_reads(self, combined_reads: Dict[str, Dict[str, Path]],
                            read_type: str, output_dir: Path) -> Optional[Path]:
        """Merge all phased reads of a type into one file"""
        merged_file = output_dir / f"merged_{read_type}.fasta"

        total_reads = 0
        with open(merged_file, 'w') as out:
            for hap_name in self.hap_names:
                if hap_name in combined_reads and read_type in combined_reads[hap_name]:
                    reads_file = combined_reads[hap_name][read_type]
                    if reads_file and reads_file.exists():
                        with open(reads_file) as inp:
                            for line in inp:
                                out.write(line)
                                if line.startswith('>'):
                                    total_reads += 1

        if total_reads > 0:
            self.logger.info(f"  Merged {read_type}: {total_reads} reads")
            return merged_file
        return None

    def _split_bam_by_prefix(self, input_bam: Path, output_bam: Path, hap_name: str):
        """Split BAM to only include alignments to specific haplotype"""
        prefix = f"{hap_name}__"

        # Use samtools view with region filtering
        # First, get list of contigs for this haplotype
        contigs = []
        with pysam.AlignmentFile(str(input_bam), 'rb') as bam:
            for ref in bam.references:
                if ref.startswith(prefix):
                    contigs.append(ref)

        if not contigs:
            # Create empty BAM
            with pysam.AlignmentFile(str(input_bam), 'rb') as inp:
                with pysam.AlignmentFile(str(output_bam), 'wb', template=inp) as out:
                    pass
            return

        # Extract alignments and rename contigs (remove prefix)
        with pysam.AlignmentFile(str(input_bam), 'rb') as inp:
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
                        # Update reference name (remove prefix)
                        a = pysam.AlignedSegment()
                        a.query_name = read.query_name
                        a.query_sequence = read.query_sequence
                        a.flag = read.flag
                        a.reference_id = out.get_tid(contig[len(prefix):])
                        a.reference_start = read.reference_start
                        a.mapping_quality = read.mapping_quality
                        a.cigar = read.cigar
                        a.query_qualities = read.query_qualities
                        if a.reference_id >= 0:
                            out.write(a)

        # Index
        pysam.index(str(output_bam))

    def _prepare_combined_reads(self, phased_hifi: Dict, phased_ont: Dict) -> Dict:
        """Combine phased reads, optionally including ambiguous"""
        combined = {hap: {} for hap in self.hap_names}

        for hap in self.hap_names:
            # HiFi
            if hap in phased_hifi:
                hifi_file = phased_hifi[hap]
                if self.use_ambiguous_reads and 'ambiguous' in phased_hifi:
                    combined_file = self.output_dir / f"{hap}_combined_hifi.fasta"
                    self._concat_files([hifi_file, phased_hifi['ambiguous']], combined_file)
                    combined[hap]['hifi'] = combined_file
                else:
                    combined[hap]['hifi'] = hifi_file

            # ONT
            if hap in phased_ont:
                ont_file = phased_ont[hap]
                if self.use_ambiguous_reads and 'ambiguous' in phased_ont:
                    combined_file = self.output_dir / f"{hap}_combined_ont.fasta"
                    self._concat_files([ont_file, phased_ont['ambiguous']], combined_file)
                    combined[hap]['ont'] = combined_file
                else:
                    combined[hap]['ont'] = ont_file

        return combined

    def _concat_files(self, files: List[Path], output: Path):
        """Concatenate files"""
        with open(output, 'w') as out:
            for f in files:
                if f and f.exists():
                    with open(f) as inp:
                        out.write(inp.read())

    def _align_reads(self, reads: Path, ref: Path, output_bam: Path, preset: str):
        """Align reads to reference"""
        # -m 2G: memory limit per thread for samtools sort
        cmd = (f"minimap2 -ax {preset} -t {self.threads} {ref} {reads} | "
               f"samtools sort -@ {self.threads} -m 2G -o {output_bam} - && "
               f"samtools index {output_bam}")
        subprocess.run(cmd, shell=True, check=True, capture_output=True)

    def _extract_haplotype_assembly(self, assembly: Path, hap_name: str,
                                     output_dir: Path) -> Path:
        """Extract assembly for gap filling (same file, just copy)"""
        output = output_dir / f"{hap_name}_assembly.fasta"
        shutil.copy(assembly, output)
        return output

    def _find_gaps(self, assembly: Path) -> List[Dict]:
        """Find gaps in assembly"""
        gaps = []
        for record in SeqIO.parse(assembly, 'fasta'):
            seq = str(record.seq).upper()
            for match in re.finditer(r'N+', seq):
                if match.end() - match.start() >= self.min_gap_size:
                    gaps.append({
                        'chrom': record.id,
                        'start': match.start(),
                        'end': match.end(),
                        'name': f"{record.id}_{match.start()}_{match.end()}",
                        'size': match.end() - match.start()
                    })
        return gaps

    def _apply_fills(self, assembly: Path, fill_results: Dict,
                      output_dir: Path) -> Path:
        """Apply fills to assembly"""
        output = output_dir / "filled_assembly.fasta"

        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        fills_by_chrom = defaultdict(list)
        for gap_name, data in fill_results.items():
            if not data['result']['success']:
                continue
            gap = data['gap']
            fills_by_chrom[gap['chrom']].append({
                'start': gap['start'],
                'end': gap['end'],
                'sequence': data['result']['sequence']
            })

        for chrom, fills in fills_by_chrom.items():
            if chrom not in sequences:
                continue
            seq = sequences[chrom]
            for fill in sorted(fills, key=lambda x: x['start'], reverse=True):
                seq = seq[:fill['start']] + fill['sequence'] + seq[fill['end']:]
            sequences[chrom] = seq

        with open(output, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        return output

    def _save_summary(self, assemblies: Dict, filled: Dict, failed: Dict, phaser: ReadPhaser):
        """Save summary statistics"""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'num_haplotypes': self.num_haplotypes,
            'optimization': 'batch_alignment',
            'gap_info': {
                hap_name: {
                    'total_gaps': len(phaser.gap_regions.get(hap_name, [])),
                    'haplotype_specific_gaps': len(phaser.haplotype_specific_gaps.get(hap_name, set()))
                }
                for hap_name in self.hap_names
            },
            'haplotypes': {}
        }

        for hap in self.hap_names:
            summary['haplotypes'][hap] = {
                'assembly': str(assemblies.get(hap, '')),
                'filled_gaps': len(filled.get(hap, set())),
                'failed_gaps': len(failed.get(hap, set()))
            }

        with open(self.output_dir / "summary.json", 'w') as f:
            json.dump(summary, f, indent=2)

    def _prepare_hic_data(self):
        """
        Prepare Hi-C BAM files for all haplotypes using merged reference alignment.

        Strategy:
        1. Create merged reference with all haplotypes (hap1__Chr1, hap2__Chr1, ...)
        2. Align Hi-C reads once to merged reference
        3. Split BAM by haplotype prefix
        4. Create per-haplotype HiCAnalyzer
        """
        # Check if user provided pre-aligned BAM
        if self.hic_bam:
            # User provided BAM - assume it's aligned to first haplotype (backward compatible)
            self.logger.info("  Using provided Hi-C BAM (aligned to single haplotype)")
            ref_assembly = self.normalized_assemblies[self.hap_names[0]]
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
