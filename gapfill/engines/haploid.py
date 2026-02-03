#!/usr/bin/env python3
"""
Haploid Gap Filling Engine

Iterative multi-round gap filling for haploid genomes.

Key features:
1. Gap normalization (all N placeholders -> 500N)
2. OPTIMIZED: One-time alignment + filtered reads cache
3. OPTIMIZED: Parallel gap filling
4. OPTIMIZED: Consensus-first strategy (skip wtdbg2 for simple gaps)
5. OPTIMIZED: High-confidence immediate validation
6. Correct coordinate update after each fill
7. Stop when all fillable gaps are completely filled
"""

import logging
import subprocess
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Set
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from gapfill.core.filler import GapFiller
from gapfill.core.validator import (GapValidator, GapStatusTracker, GapStatus, PendingFill,
                                     PartialFillValidationResult)
from gapfill.core.parallel import fill_gaps_parallel, fill_gaps_sequential, GapBatcher
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads
from gapfill.utils.checkpoint import CheckpointManager, CheckpointState
from gapfill.utils.reads_cache import ReadsCache, GapReadsExtractor


# =============================================================================
# Progress Marker (.ok file) utilities
# =============================================================================

def _mark_done(path: Path, step: str) -> Path:
    """Create a .ok marker file to indicate step completion"""
    ok_file = path / f"{step}.ok"
    ok_file.touch()
    return ok_file


def _is_done(path: Path, step: str) -> bool:
    """Check if a step has been completed (has .ok file)"""
    ok_file = path / f"{step}.ok"
    return ok_file.exists()


def _clear_marker(path: Path, step: str):
    """Remove a .ok marker file"""
    ok_file = path / f"{step}.ok"
    if ok_file.exists():
        ok_file.unlink()


class HaploidEngine:
    """
    Iterative gap filler for haploid genomes
    """

    def __init__(self,
                 assembly_file: str,
                 hifi_reads: Optional[str] = None,
                 ont_reads: Optional[str] = None,
                 hic_reads: Optional[List[str]] = None,
                 hic_bam: Optional[str] = None,
                 output_dir: str = "output",
                 threads: int = 8,
                 max_iterations: int = 10,
                 min_gap_size: int = 100,
                 min_mapq: int = 20,
                 skip_normalization: bool = False,
                 resume: bool = False,
                 clear_checkpoint: bool = False,
                 optimized_mode: bool = True,
                 parallel_filling: bool = True):

        self.initial_assembly = Path(assembly_file)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.hic_reads = hic_reads  # [R1, R2] or None
        self.hic_bam = Path(hic_bam) if hic_bam else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.max_iterations = max_iterations
        self.min_gap_size = min_gap_size
        self.min_mapq = min_mapq
        self.skip_normalization = skip_normalization
        self.resume = resume
        self.optimized_mode = optimized_mode
        self.parallel_filling = parallel_filling

        # Cache for gap scan results
        self._gap_cache: Dict[str, List[Dict]] = {}

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Checkpoint manager
        self.checkpoint = CheckpointManager(str(self.output_dir))
        if clear_checkpoint:
            self.checkpoint.clear()

        # Hi-C analyzer (initialized when needed)
        self.hic_analyzer: Optional[HiCAnalyzer] = None
        self.gap_size_estimates: Dict[str, int] = {}

        # Gap status tracker (tracks UNFILLABLE, NEEDS_POLISH, etc.)
        self.gap_tracker = GapStatusTracker()
        self.gap_tracker_file = self.output_dir / "gap_tracker.json"

        # Reads cache for optimized mode (filter once, use many times)
        self.reads_cache: Optional[ReadsCache] = None
        if self.optimized_mode:
            self.reads_cache = ReadsCache(self.output_dir, threads=threads)

        # Preprocessing BAM tracking (for potential reuse in iteration 1)
        self._preprocessing_hifi_bam: Optional[Path] = None
        self._preprocessing_ont_bam: Optional[Path] = None
        self._preprocessing_assembly: Optional[Path] = None
        self._can_reuse_preprocessing_bam: bool = False

        # Store initial gaps for later Hi-C validation
        self.initial_gaps: List[Dict] = []

        self.stats = {
            'iterations': 0,
            'total_gaps_initial': 0,
            'gaps_completely_filled': 0,
            'gaps_partially_filled': 0,
            'gaps_failed': 0,
            'total_bp_filled': 0,
            'hic_validated': 0,
            'hic_failed': 0
        }

        self.logger.info("=" * 60)
        self.logger.info("HaploidEngine initialized")
        self.logger.info(f"  Assembly: {self.initial_assembly}")
        self.logger.info(f"  HiFi reads: {self.hifi_reads}")
        self.logger.info(f"  ONT reads: {self.ont_reads}")
        self.logger.info(f"  Hi-C reads: {self.hic_reads}")
        self.logger.info(f"  Hi-C BAM: {self.hic_bam}")
        self.logger.info(f"  Output dir: {self.output_dir}")
        self.logger.info(f"  Resume: {self.resume}")
        self.logger.info(f"  Optimized mode: {self.optimized_mode}")
        self.logger.info(f"  Parallel filling: {self.parallel_filling}")
        self.logger.info("=" * 60)

    def run(self) -> Path:
        """Run iterative gap filling with delayed validation"""
        current_assembly = self.initial_assembly
        start_iteration = 0

        # Track gap states
        completely_filled_gaps: Set[str] = set()
        partially_filled_gaps: Set[str] = set()
        failed_gaps: Set[str] = set()

        # =====================================================================
        # RESUME DETECTION (using .ok files)
        # =====================================================================
        if self.resume:
            self.logger.info("=" * 60)
            self.logger.info("CHECKING FOR RESUME POINTS...")

            # Load gap tracker state if exists
            if self._load_gap_tracker():
                self.logger.info("  Loaded gap tracker state")

            # Check preprocessing
            if _is_done(self.output_dir, "preprocessing"):
                self.logger.info("  Preprocessing: DONE")
                # Find the preprocessed assembly
                polished = self.output_dir / "assembly_polished.fasta"
                normalized = self.output_dir / "assembly_normalized.fasta"
                if polished.exists():
                    current_assembly = polished
                elif normalized.exists():
                    current_assembly = normalized
            else:
                self.logger.info("  Preprocessing: NOT DONE")

            # Find latest completed iteration
            for i in range(self.max_iterations, 0, -1):
                iter_dir = self.output_dir / f"iteration_{i}"
                if _is_done(iter_dir, "iteration"):
                    start_iteration = i
                    # Use the filled assembly from this iteration
                    filled_asm = iter_dir / "assembly_filled.fasta"
                    reverted_asm = iter_dir / "assembly_reverted.fasta"
                    if reverted_asm.exists():
                        current_assembly = reverted_asm
                    elif filled_asm.exists():
                        current_assembly = filled_asm
                    self.logger.info(f"  Last completed iteration: {i}")
                    break

            self.logger.info(f"  Current assembly: {current_assembly}")
            self.logger.info("=" * 60)

        # =====================================================================
        # PREPROCESSING PHASE (before iteration loop)
        # =====================================================================
        if not _is_done(self.output_dir, "preprocessing"):
            current_assembly = self._run_preprocessing(current_assembly)
            _mark_done(self.output_dir, "preprocessing")
            self._save_gap_tracker()
        else:
            self.logger.info("Skipping preprocessing (already done)")
            # Load initial gaps for later Hi-C validation
            self.initial_gaps = self._find_gaps(current_assembly)
            self.stats['total_gaps_initial'] = len(self.initial_gaps)

        # =====================================================================
        # ITERATION LOOP
        # =====================================================================
        iteration = start_iteration
        while iteration < self.max_iterations:
            iteration += 1
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"ITERATION {iteration}")
            self.logger.info(f"{'='*60}")

            iter_dir = self.output_dir / f"iteration_{iteration}"
            iter_dir.mkdir(exist_ok=True)

            # Check if this iteration is already complete
            if _is_done(iter_dir, "iteration"):
                self.logger.info("  Iteration already complete, skipping...")
                # Load assembly from this iteration
                reverted_asm = iter_dir / "assembly_reverted.fasta"
                filled_asm = iter_dir / "assembly_filled.fasta"
                if reverted_asm.exists():
                    current_assembly = reverted_asm
                elif filled_asm.exists():
                    current_assembly = filled_asm
                continue

            # Step 1: Align reads (use filtered reads in optimized mode)
            hifi_bam = None
            ont_bam = None

            if not _is_done(iter_dir, "alignment"):
                self.logger.info("Step 1: Aligning reads...")

                # OPTIMIZATION: Reuse preprocessing BAM for iteration 1 if no polishing occurred
                if (iteration == 1 and
                    self._can_reuse_preprocessing_bam and
                    self._preprocessing_assembly == current_assembly):

                    self.logger.info("  Reusing preprocessing BAM (no polishing occurred)...")
                    hifi_bam = self._preprocessing_hifi_bam
                    ont_bam = self._preprocessing_ont_bam

                    # Create symlinks in iteration directory for consistency
                    if hifi_bam and hifi_bam.exists():
                        iter_hifi_bam = iter_dir / "hifi.bam"
                        iter_hifi_bai = iter_dir / "hifi.bam.bai"
                        if not iter_hifi_bam.exists():
                            iter_hifi_bam.symlink_to(hifi_bam)
                        if hifi_bam.with_suffix('.bam.bai').exists() and not iter_hifi_bai.exists():
                            iter_hifi_bai.symlink_to(hifi_bam.with_suffix('.bam.bai'))
                        hifi_bam = iter_hifi_bam

                    if ont_bam and ont_bam.exists():
                        iter_ont_bam = iter_dir / "ont.bam"
                        iter_ont_bai = iter_dir / "ont.bam.bai"
                        if not iter_ont_bam.exists():
                            iter_ont_bam.symlink_to(ont_bam)
                        if ont_bam.with_suffix('.bam.bai').exists() and not iter_ont_bai.exists():
                            iter_ont_bai.symlink_to(ont_bam.with_suffix('.bam.bai'))
                        ont_bam = iter_ont_bam

                elif self.optimized_mode and self.reads_cache:
                    # Use filtered reads for faster alignment
                    hifi_bam, ont_bam = self._align_filtered_reads(current_assembly, iter_dir)
                else:
                    hifi_bam, ont_bam = self._align_reads_with_cache(current_assembly, iter_dir)

                _mark_done(iter_dir, "alignment")
            else:
                self.logger.info("Step 1: Alignment already done, reusing...")
                hifi_bam = iter_dir / "hifi.bam" if (iter_dir / "hifi.bam").exists() else None
                ont_bam = iter_dir / "ont.bam" if (iter_dir / "ont.bam").exists() else None

            if not hifi_bam and not ont_bam:
                self.logger.error("No BAM files generated")
                break

            # Step 2: Validate previous iteration's fills (if any)
            if not _is_done(iter_dir, "validation"):
                assembly_before_validation = current_assembly
                if iteration > 1:
                    pending_fills = self.gap_tracker.get_pending_fills()
                    if pending_fills:
                        self.logger.info(f"Step 2: Validating {len(pending_fills)} fills from previous iteration...")
                        current_assembly = self._validate_previous_fills(
                            current_assembly, pending_fills, hifi_bam, ont_bam, iter_dir
                        )

                        # If assembly changed (reverts occurred), re-align reads
                        if current_assembly != assembly_before_validation:
                            self.logger.info("Step 2b: Re-aligning reads to reverted assembly...")
                            # Clear alignment.ok since we need to re-align
                            _clear_marker(iter_dir, "alignment")
                            if self.optimized_mode and self.reads_cache:
                                hifi_bam, ont_bam = self._align_filtered_reads(current_assembly, iter_dir)
                            else:
                                hifi_bam, ont_bam = self._align_reads_with_cache(current_assembly, iter_dir)
                            _mark_done(iter_dir, "alignment")

                _mark_done(iter_dir, "validation")
                self._save_gap_tracker()
            else:
                self.logger.info("Step 2: Validation already done, skipping...")
                # Check if there's a reverted assembly from previous validation
                reverted_asm = iter_dir / "assembly_reverted.fasta"
                if reverted_asm.exists():
                    current_assembly = reverted_asm

            # Step 3: Find gaps
            self.logger.info("Step 3: Finding gaps...")
            gaps = self._find_gaps(current_assembly)

            # Use gap tracker to determine which gaps to attempt
            remaining_gaps = [g for g in gaps
                            if self.gap_tracker.should_attempt(g['name'])]

            self.logger.info(f"  Remaining to process: {len(remaining_gaps)}")

            # Log skip summary
            tracker_summary = self.gap_tracker.get_summary()
            if tracker_summary:
                self.logger.info(f"  Status summary: {tracker_summary}")

            if not remaining_gaps:
                self.logger.info("No remaining gaps, stopping")
                break

            # Step 4: Fill gaps
            self.logger.info("Step 4: Filling gaps...")
            work_dir = iter_dir / "work"
            work_dir.mkdir(exist_ok=True)

            if self.parallel_filling and len(remaining_gaps) > 1:
                # Parallel filling for multiple gaps
                self.logger.info(f"  Using parallel filling ({self.threads} workers)...")
                parallel_results = fill_gaps_parallel(
                    gaps=remaining_gaps,
                    assembly_file=str(current_assembly),
                    hifi_bam=str(hifi_bam) if hifi_bam else None,
                    ont_bam=str(ont_bam) if ont_bam else None,
                    hifi_reads=str(self.hifi_reads) if self.hifi_reads else None,
                    ont_reads=str(self.ont_reads) if self.ont_reads else None,
                    work_dir=str(work_dir),
                    threads=self.threads,
                    min_mapq=self.min_mapq,
                    use_consensus_first=True
                )
                # Convert parallel results to expected format
                fill_results = {}
                gap_lookup = {g['name']: g for g in remaining_gaps}
                for gap_name, pr in parallel_results.items():
                    if pr.get('result'):
                        fill_results[gap_name] = {
                            'gap': gap_lookup.get(gap_name, {'name': gap_name}),
                            'result': pr['result']
                        }
                    else:
                        fill_results[gap_name] = {
                            'gap': gap_lookup.get(gap_name, {'name': gap_name}),
                            'result': {'success': False, 'reason': pr.get('error', 'Unknown error')}
                        }
            else:
                # Sequential filling
                gap_filler = GapFiller(
                    assembly_file=str(current_assembly),
                    hifi_bam=str(hifi_bam) if hifi_bam else None,
                    ont_bam=str(ont_bam) if ont_bam else None,
                    hifi_reads=str(self.hifi_reads) if self.hifi_reads else None,
                    ont_reads=str(self.ont_reads) if self.ont_reads else None,
                    threads=self.threads,
                    work_dir=str(work_dir),
                    min_mapq=self.min_mapq,
                    hic_analyzer=self.hic_analyzer,
                    gap_size_estimates=self.gap_size_estimates,
                    enable_consensus_first=True
                )

                fill_results = {}
                for gap in remaining_gaps:
                    result = gap_filler.fill_gap(gap)
                    fill_results[gap['name']] = {'gap': gap, 'result': result}

                gap_filler.close()

            # Step 5: Process fill results (with HIGH-CONFIDENCE IMMEDIATE VALIDATION)
            new_pending = 0
            new_complete = 0  # High-confidence fills that skip delayed validation
            new_unfillable = 0
            new_failed = 0

            # Track successful fills to apply to assembly
            successful_fills = {}

            for gap_name, data in fill_results.items():
                gap = data['gap']
                result = data['result']

                if result.get('success'):
                    # Fill succeeded
                    successful_fills[gap_name] = data
                    is_complete = result.get('is_complete', False) and not result.get('has_placeholder', False)

                    # Check if this fill can skip delayed validation (high confidence)
                    validation = result.get('validation', {})
                    skip_delayed = validation.get('skip_delayed_validation', False)

                    fill_seq = result.get('sequence', '')
                    fill_length = len(fill_seq)
                    original_gap_size = gap['end'] - gap['start']

                    if skip_delayed and is_complete:
                        # HIGH CONFIDENCE: Mark as FILLED_COMPLETE immediately
                        new_complete += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FILLED_COMPLETE,
                            f"High-confidence fill: {result.get('source', 'unknown')} "
                            f"(tier={result.get('tier', 0)}, cov={validation.get('avg_coverage', 0):.1f}x)"
                        )
                        self.checkpoint.add_completed_gap(gap_name, fill_seq)
                        self.logger.info(f"    {gap_name}: ✓ filled & validated "
                                       f"({result.get('source', 'unknown')}, "
                                       f"{original_gap_size}bp -> {fill_length}bp)")
                    else:
                        # Standard: Mark as PENDING validation
                        new_pending += 1
                        self.gap_tracker.add_pending_fill(
                            gap_id=gap_name,
                            chrom=gap['chrom'],
                            original_start=gap['start'],
                            original_end=gap['end'],
                            filled_start=gap['start'],
                            filled_end=gap['start'] + fill_length,
                            sequence=fill_seq,
                            is_complete=is_complete,
                            source=result.get('source', 'unknown'),
                            tier=result.get('tier', 0)
                        )
                        self.logger.info(f"    {gap_name}: filled ({result.get('source', 'unknown')}, "
                                       f"{original_gap_size}bp -> {fill_length}bp), pending validation")
                else:
                    # Fill failed → analyze why
                    validation = result.get('validation', {})
                    validation_status = validation.get('status', '')

                    if validation_status == 'unfillable':
                        new_unfillable += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.UNFILLABLE,
                            validation.get('reason', 'No spanning reads available')
                        )
                    else:
                        new_failed += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FAILED,
                            validation.get('reason', result.get('reason', 'Unknown failure'))
                        )

                    failed_gaps.add(gap_name)
                    self.checkpoint.add_failed_gap(gap_name)

            self.logger.info(f"  Filled (confirmed): {new_complete}, "
                           f"Filled (pending validation): {new_pending}, "
                           f"Unfillable: {new_unfillable}, Failed: {new_failed}")

            # Step 6: Apply fills to assembly
            if successful_fills:
                current_assembly = self._apply_fills(
                    current_assembly, successful_fills, iter_dir
                )

            # Save iteration stats
            self._save_iteration_stats(iter_dir, iteration, fill_results)

            # Log status summary
            tracker_summary = self.gap_tracker.get_summary()
            self.logger.info(f"  Status summary: {tracker_summary}")

            # Save gap tracker state and mark iteration complete
            self._save_gap_tracker()
            _mark_done(iter_dir, "iteration")

            # Check for progress (no new fills and no pending validation)
            pending_count = len(self.gap_tracker.get_pending_fills())
            if new_pending == 0 and pending_count == 0:
                self.logger.info("No progress in this iteration, stopping")
                break

        # Save final results
        final_assembly = self.output_dir / "final_assembly.fasta"
        import shutil
        shutil.copy(current_assembly, final_assembly)

        self.stats['iterations'] = iteration
        self.stats['gaps_completely_filled'] = len(completely_filled_gaps)
        self.stats['gaps_partially_filled'] = len(partially_filled_gaps)
        self.stats['gaps_failed'] = len(failed_gaps)
        self.stats['gap_status_summary'] = self.gap_tracker.get_summary()

        # Hi-C validation of filled gaps
        if self.hic_analyzer and completely_filled_gaps:
            self.logger.info("\nValidating fills with Hi-C...")
            filled_gaps_list = [
                gap for gap in self.initial_gaps
                if gap['name'] in completely_filled_gaps
            ]
            if filled_gaps_list:
                validations = self.hic_analyzer.validate_fills(
                    filled_gaps_list, str(final_assembly)
                )
                for v in validations:
                    if v.is_valid:
                        self.stats['hic_validated'] += 1
                    else:
                        self.stats['hic_failed'] += 1
                        self.logger.warning(f"  Hi-C validation failed: {v.gap_name} "
                                          f"(anomaly: {v.anomaly_type})")

        self._save_final_stats()
        self._save_gap_tracker()

        # Mark as complete
        _mark_done(self.output_dir, "complete")
        self.logger.info("Gap filling complete")

        self.logger.info(f"\nFinal assembly: {final_assembly}")
        return final_assembly

    def _align_reads_with_cache(self, assembly: Path, output_dir: Path) -> tuple:
        """Align reads to assembly, reusing existing BAMs if valid"""
        hifi_bam = None
        ont_bam = None

        if self.hifi_reads:
            hifi_bam = output_dir / "hifi.bam"
            hifi_bai = output_dir / "hifi.bam.bai"
            # Check if BAM exists and is valid
            if self.resume and hifi_bam.exists() and hifi_bai.exists() and hifi_bam.stat().st_size > 0:
                self.logger.info(f"  Reusing existing HiFi BAM: {hifi_bam}")
            else:
                self._run_minimap2(assembly, self.hifi_reads, hifi_bam, 'map-hifi')

        if self.ont_reads:
            ont_bam = output_dir / "ont.bam"
            ont_bai = output_dir / "ont.bam.bai"
            # Check if BAM exists and is valid
            if self.resume and ont_bam.exists() and ont_bai.exists() and ont_bam.stat().st_size > 0:
                self.logger.info(f"  Reusing existing ONT BAM: {ont_bam}")
            else:
                self._run_minimap2(assembly, self.ont_reads, ont_bam, 'map-ont')

        return hifi_bam, ont_bam

    def _align_filtered_reads(self, assembly: Path, output_dir: Path) -> tuple:
        """
        Align filtered reads to assembly (OPTIMIZED).

        Uses pre-filtered reads from ReadsCache instead of original reads.
        This significantly reduces alignment time as we only align reads
        that are potentially useful for gap filling.
        """
        hifi_bam = None
        ont_bam = None

        if not self.reads_cache:
            # Fall back to regular alignment
            return self._align_reads_with_cache(assembly, output_dir)

        if self.hifi_reads:
            hifi_bam = output_dir / "hifi.bam"
            hifi_bai = output_dir / "hifi.bam.bai"

            if self.resume and hifi_bam.exists() and hifi_bai.exists() and hifi_bam.stat().st_size > 0:
                self.logger.info(f"  Reusing existing HiFi BAM: {hifi_bam}")
            else:
                filtered_fasta = self.reads_cache.get_filtered_reads_path('hifi')
                if filtered_fasta and filtered_fasta.exists():
                    self.logger.info(f"  Aligning filtered HiFi reads...")
                    self._run_minimap2(assembly, filtered_fasta, hifi_bam, 'map-hifi')
                else:
                    # No filtered reads, use original
                    self.logger.info(f"  Aligning original HiFi reads (no cache)...")
                    self._run_minimap2(assembly, self.hifi_reads, hifi_bam, 'map-hifi')

        if self.ont_reads:
            ont_bam = output_dir / "ont.bam"
            ont_bai = output_dir / "ont.bam.bai"

            if self.resume and ont_bam.exists() and ont_bai.exists() and ont_bam.stat().st_size > 0:
                self.logger.info(f"  Reusing existing ONT BAM: {ont_bam}")
            else:
                filtered_fasta = self.reads_cache.get_filtered_reads_path('ont')
                if filtered_fasta and filtered_fasta.exists():
                    self.logger.info(f"  Aligning filtered ONT reads...")
                    self._run_minimap2(assembly, filtered_fasta, ont_bam, 'map-ont')
                else:
                    # No filtered reads, use original
                    self.logger.info(f"  Aligning original ONT reads (no cache)...")
                    self._run_minimap2(assembly, self.ont_reads, ont_bam, 'map-ont')

        return hifi_bam, ont_bam

    def _run_preprocessing(self, assembly: Path) -> Path:
        """
        Run preprocessing phase before iteration loop.

        Steps:
        1. Hi-C preparation (if available)
        2. Gap normalization
        3. Initial alignment
        4. Pre-assess gap flanks
        5. Polish problematic flanks
        6. Re-align after polish (if needed)

        Returns:
            Path to preprocessed assembly (normalized and polished)
        """
        current_assembly = assembly
        preprocess_dir = self.output_dir / "preprocessing"
        preprocess_dir.mkdir(exist_ok=True)

        # Step 0a: Find initial gaps
        initial_gaps = self._find_gaps(current_assembly)
        self.initial_gaps = initial_gaps  # Store for later Hi-C validation
        self.stats['total_gaps_initial'] = len(initial_gaps)
        self.logger.info(f"Initial gaps found: {len(initial_gaps)}")

        # Step 0b: Prepare Hi-C data if available
        if self.hic_reads or self.hic_bam:
            hic_bam_file = preprocess_dir / "hic.bam"
            if _is_done(preprocess_dir, "hic") and hic_bam_file.exists():
                self.logger.info(f"STEP 0b: Reusing existing Hi-C BAM")
                self.hic_bam = hic_bam_file
            else:
                self.logger.info("STEP 0b: Preparing Hi-C data")
                self._prepare_hic_data(current_assembly)
                _mark_done(preprocess_dir, "hic")

            # Estimate gap sizes using Hi-C
            if self.hic_analyzer and initial_gaps:
                self.logger.info("  Estimating gap sizes with Hi-C...")
                estimates = self.hic_analyzer.estimate_gap_sizes(initial_gaps)
                for est in estimates:
                    if est.confidence in ('high', 'medium'):
                        self.gap_size_estimates[est.gap_name] = est.estimated_size
                        self.logger.info(f"    {est.gap_name}: estimated {est.estimated_size}bp")

        # Step 0c: Normalize gaps
        normalized_file = self.output_dir / "assembly_normalized.fasta"
        if _is_done(self.output_dir, "normalized") and normalized_file.exists():
            self.logger.info(f"STEP 0c: Reusing normalized assembly")
            current_assembly = normalized_file
        elif initial_gaps and not self.skip_normalization:
            self.logger.info("STEP 0c: Normalizing gaps to 500N...")
            current_assembly = self._normalize_gaps(current_assembly, initial_gaps)
            _mark_done(self.output_dir, "normalized")
        elif self.skip_normalization:
            self.logger.info("STEP 0c: Skipping normalization (already done)")
            import shutil
            shutil.copy(current_assembly, normalized_file)
            current_assembly = normalized_file
            _mark_done(self.output_dir, "normalized")

        # Step 0d: Initial alignment for pre-assessment
        if not _is_done(preprocess_dir, "alignment"):
            self.logger.info("STEP 0d: Initial alignment for pre-assessment...")
            hifi_bam, ont_bam = self._align_reads_with_cache(current_assembly, preprocess_dir)
            _mark_done(preprocess_dir, "alignment")
        else:
            self.logger.info("STEP 0d: Reusing existing alignment...")
            hifi_bam = preprocess_dir / "hifi.bam" if (preprocess_dir / "hifi.bam").exists() else None
            ont_bam = preprocess_dir / "ont.bam" if (preprocess_dir / "ont.bam").exists() else None

        # Track preprocessing BAM paths for potential reuse in iteration 1
        self._preprocessing_hifi_bam = hifi_bam
        self._preprocessing_ont_bam = ont_bam
        self._preprocessing_assembly = current_assembly

        if not hifi_bam and not ont_bam:
            self.logger.warning("  No BAM files for pre-assessment, skipping polish")
            return current_assembly

        # Step 0e: Filter reads (OPTIMIZATION - filter once, use many times)
        if self.optimized_mode and self.reads_cache:
            if not _is_done(preprocess_dir, "filter"):
                self.logger.info("STEP 0e: Filtering reads (keeping gap-related reads only)...")
                gaps = self._find_gaps(current_assembly)
                self.reads_cache.set_gap_regions(gaps)

                if hifi_bam and not self.reads_cache.is_cached('hifi'):
                    self.reads_cache.filter_bam(hifi_bam, 'hifi', self.min_mapq)
                    cache_summary = self.reads_cache.get_summary()
                    self.logger.info(f"  HiFi: kept {cache_summary['hifi']['stats']['kept']:,} / "
                                   f"{cache_summary['hifi']['stats']['total']:,} reads "
                                   f"({cache_summary['hifi']['stats']['kept_ratio']*100:.1f}%)")

                if ont_bam and not self.reads_cache.is_cached('ont'):
                    self.reads_cache.filter_bam(ont_bam, 'ont', self.min_mapq)
                    cache_summary = self.reads_cache.get_summary()
                    self.logger.info(f"  ONT: kept {cache_summary['ont']['stats']['kept']:,} / "
                                   f"{cache_summary['ont']['stats']['total']:,} reads "
                                   f"({cache_summary['ont']['stats']['kept_ratio']*100:.1f}%)")

                _mark_done(preprocess_dir, "filter")
            else:
                self.logger.info("STEP 0e: Reads already filtered")
                # Re-initialize reads_cache with existing gaps
                gaps = self._find_gaps(current_assembly)
                self.reads_cache.set_gap_regions(gaps)

        # Step 0f: Pre-assess gap flanks
        polished_file = self.output_dir / "assembly_polished.fasta"
        polishing_occurred = False

        if not _is_done(preprocess_dir, "polish"):
            self.logger.info("STEP 0f: Pre-assessing gap flanks...")
            gaps = self._find_gaps(current_assembly)
            gaps_needing_polish = self._pre_assess_gaps(gaps, hifi_bam, ont_bam)

            # Step 0g: Polish problematic flanks
            if gaps_needing_polish:
                self.logger.info(f"STEP 0g: Polishing {len(gaps_needing_polish)} gap flanks...")
                polished_assembly = self._polish_gap_flanks(
                    current_assembly, gaps_needing_polish, hifi_bam, ont_bam, preprocess_dir
                )

                if polished_assembly != current_assembly:
                    current_assembly = polished_assembly
                    self.logger.info("  Polishing complete, will re-align in first iteration")
                    polishing_occurred = True
            else:
                self.logger.info("  No flanks need polishing")

            _mark_done(preprocess_dir, "polish")
        else:
            self.logger.info("STEP 0f-0g: Polish already done")
            if polished_file.exists():
                current_assembly = polished_file
                polishing_occurred = True

        # Mark whether preprocessing BAM can be reused for iteration 1
        self._can_reuse_preprocessing_bam = not polishing_occurred

        return current_assembly

    def _validate_previous_fills(self,
                                  assembly: Path,
                                  pending_fills: Dict,
                                  hifi_bam: Optional[Path],
                                  ont_bam: Optional[Path],
                                  output_dir: Path) -> Path:
        """
        Validate fills from previous iteration using current BAM.

        For each pending fill:
        - If validation passes → mark as FILLED_COMPLETE/PARTIAL
        - If validation fails → revert the fill (restore gap)

        Returns:
            Path to assembly (may be modified if fills reverted)
        """
        validator = GapValidator(threads=self.threads)

        # Prefer HiFi BAM for validation
        bam_file = None
        if hifi_bam and hifi_bam.exists():
            bam_file = str(hifi_bam)
        elif ont_bam and ont_bam.exists():
            bam_file = str(ont_bam)

        if not bam_file:
            self.logger.warning("  No BAM for validation, accepting all fills")
            for gap_id, pf in pending_fills.items():
                if pf.is_complete:
                    self.gap_tracker.set_status(gap_id, GapStatus.FILLED_COMPLETE,
                                               "Accepted without validation (no BAM)")
                else:
                    self.gap_tracker.set_status(gap_id, GapStatus.FILLED_PARTIAL,
                                               "Accepted without validation (no BAM)")
            validator.close()
            return assembly

        validated_count = 0
        failed_count = 0
        fills_to_revert = []

        # Track different types of reverts
        full_reverts = []      # Completely failed fills
        partial_reverts = []   # Partial fills where one side failed

        try:
            for gap_id, pf in pending_fills.items():
                self.logger.info(f"    Validating {gap_id}...")

                fill_length = pf.filled_end - pf.filled_start

                # Validate based on fill type
                if pf.is_complete:
                    # Complete fills: require spanning reads
                    result = validator.validate_complete_fill(
                        bam_file, pf.chrom, pf.filled_start, pf.filled_end, pf.sequence
                    )

                    if result.valid:
                        validated_count += 1
                        self.gap_tracker.set_status(gap_id, GapStatus.FILLED_COMPLETE,
                                                   f"Validated: {result.reason}")
                        self.checkpoint.add_completed_gap(gap_id, pf.sequence)
                        self.logger.info(f"      ✓ Validated {fill_length}bp complete fill "
                                       f"(cov={result.avg_coverage:.1f}x, spanning={result.spanning_reads})")
                    else:
                        failed_count += 1
                        full_reverts.append(pf)
                        self.gap_tracker.set_status(gap_id, GapStatus.FAILED,
                                                   f"Validation failed: {result.reason}")
                        self.logger.warning(f"      ✗ Failed {fill_length}bp complete fill: {result.reason}")
                else:
                    # Partial fills: check junctions and coverage independently for each side
                    n_match = re.search(r'N{100,}', pf.sequence)
                    if n_match:
                        new_gap_start = pf.filled_start + n_match.start()
                        new_gap_end = pf.filled_start + n_match.end()

                        result = validator.validate_partial_fill(
                            bam_file, pf.chrom,
                            pf.original_start, pf.original_end,
                            pf.sequence,
                            new_gap_start, new_gap_end
                        )

                        # Handle independent left/right validation
                        if isinstance(result, PartialFillValidationResult):
                            if result.left_valid and result.right_valid:
                                # Both sides passed
                                validated_count += 1
                                self.gap_tracker.set_status(gap_id, GapStatus.FILLED_PARTIAL,
                                                           f"Validated: {result.reason}")
                                self.checkpoint.add_partial_gap(gap_id, pf.sequence)
                                self.logger.info(f"      ✓ Validated {fill_length}bp partial fill: {result.reason}")

                            elif result.left_valid or result.right_valid:
                                # One side passed, one failed - partial revert
                                validated_count += 1  # Count as partial success
                                partial_reverts.append({
                                    'pf': pf,
                                    'result': result,
                                    'n_match': n_match,
                                    'new_gap_start': new_gap_start,
                                    'new_gap_end': new_gap_end
                                })
                                side_kept = "left" if result.left_valid else "right"
                                side_failed = "right" if result.left_valid else "left"
                                self.gap_tracker.set_status(gap_id, GapStatus.FILLED_PARTIAL,
                                                           f"Partial: {side_kept} OK, {side_failed} reverted")
                                self.logger.info(f"      ⚠ Partial validation: {result.reason}")
                                self.logger.info(f"        Keeping {side_kept} side, reverting {side_failed} side")

                            else:
                                # Both sides failed
                                failed_count += 1
                                full_reverts.append(pf)
                                self.gap_tracker.set_status(gap_id, GapStatus.FAILED,
                                                           f"Validation failed: {result.reason}")
                                self.logger.warning(f"      ✗ Failed {fill_length}bp partial fill: {result.reason}")
                        else:
                            # Fallback for non-PartialFillValidationResult
                            if result.valid:
                                validated_count += 1
                                self.gap_tracker.set_status(gap_id, GapStatus.FILLED_PARTIAL,
                                                           f"Validated: {result.reason}")
                                self.checkpoint.add_partial_gap(gap_id, pf.sequence)
                                self.logger.info(f"      ✓ Validated {fill_length}bp partial fill")
                            else:
                                failed_count += 1
                                full_reverts.append(pf)
                                self.gap_tracker.set_status(gap_id, GapStatus.FAILED,
                                                           f"Validation failed: {result.reason}")
                                self.logger.warning(f"      ✗ Failed {fill_length}bp partial fill: {result.reason}")
                    else:
                        # No N-run found but marked as partial - treat as complete
                        result = validator.validate_complete_fill(
                            bam_file, pf.chrom, pf.filled_start, pf.filled_end, pf.sequence
                        )
                        if result.valid:
                            validated_count += 1
                            self.gap_tracker.set_status(gap_id, GapStatus.FILLED_COMPLETE,
                                                       f"Validated: {result.reason}")
                            self.checkpoint.add_completed_gap(gap_id, pf.sequence)
                            self.logger.info(f"      ✓ Validated {fill_length}bp fill")
                        else:
                            failed_count += 1
                            full_reverts.append(pf)
                            self.gap_tracker.set_status(gap_id, GapStatus.FAILED,
                                                       f"Validation failed: {result.reason}")
                            self.logger.warning(f"      ✗ Failed {fill_length}bp fill: {result.reason}")

        finally:
            validator.close()

        self.logger.info(f"  Validation complete: {validated_count} passed, {failed_count} failed")

        # Apply reverts
        if full_reverts:
            self.logger.info(f"  Reverting {len(full_reverts)} fully failed fills...")
            assembly = self._revert_fills(assembly, full_reverts, output_dir)

        if partial_reverts:
            self.logger.info(f"  Applying {len(partial_reverts)} partial reverts...")
            assembly = self._partial_revert_fills(assembly, partial_reverts, output_dir)

        return assembly

    def _revert_fills(self, assembly: Path, fills_to_revert: List, output_dir: Path) -> Path:
        """
        Revert failed fills by restoring the original gap (500N).

        This is necessary when a fill fails validation - we need to
        restore the gap so it can be retried in future iterations.
        """
        if not fills_to_revert:
            return assembly

        reverted_assembly = output_dir / "assembly_reverted.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Revert each failed fill (process from end to avoid coordinate shift)
        fills_by_chrom = {}
        for pf in fills_to_revert:
            if pf.chrom not in fills_by_chrom:
                fills_by_chrom[pf.chrom] = []
            fills_by_chrom[pf.chrom].append(pf)

        for chrom, fills in fills_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]
            # Sort by position, descending
            fills_sorted = sorted(fills, key=lambda x: x.filled_start, reverse=True)

            for pf in fills_sorted:
                # Replace filled region with 500N gap
                seq = seq[:pf.filled_start] + 'N' * 500 + seq[pf.filled_end:]
                self.logger.info(f"    Reverted {pf.gap_id}: restored 500N gap")

            sequences[chrom] = seq

        # Write reverted assembly
        with open(reverted_assembly, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        return reverted_assembly

    def _partial_revert_fills(self, assembly: Path, partial_reverts: List[Dict],
                               output_dir: Path) -> Path:
        """
        Partially revert fills where one side passed and one side failed.

        For each partial revert:
        - If left failed, right passed: keep right_fill, replace left_fill+gap with 500N
        - If right failed, left passed: keep left_fill, replace gap+right_fill with 500N

        Args:
            assembly: Current assembly file
            partial_reverts: List of dicts with 'pf', 'result', 'n_match', etc.
            output_dir: Output directory

        Returns:
            Path to modified assembly
        """
        if not partial_reverts:
            return assembly

        reverted_assembly = output_dir / "assembly_partial_reverted.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group by chromosome
        reverts_by_chrom = {}
        for revert_info in partial_reverts:
            pf = revert_info['pf']
            if pf.chrom not in reverts_by_chrom:
                reverts_by_chrom[pf.chrom] = []
            reverts_by_chrom[pf.chrom].append(revert_info)

        # Process each chromosome (from end to start to avoid coordinate issues)
        for chrom, reverts in reverts_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]
            # Sort by position, descending
            reverts_sorted = sorted(reverts, key=lambda x: x['pf'].filled_start, reverse=True)

            for revert_info in reverts_sorted:
                pf = revert_info['pf']
                result = revert_info['result']
                n_match = revert_info['n_match']

                # Calculate positions within the fill sequence
                # pf.sequence structure: left_fill + 500N + right_fill
                left_fill_end = n_match.start()      # End of left_fill in pf.sequence
                right_fill_start = n_match.end()     # Start of right_fill in pf.sequence

                # Absolute positions in assembly
                abs_left_fill_start = pf.filled_start
                abs_left_fill_end = pf.filled_start + left_fill_end
                abs_gap_start = abs_left_fill_end
                abs_gap_end = pf.filled_start + right_fill_start
                abs_right_fill_start = abs_gap_end
                abs_right_fill_end = pf.filled_end

                if not result.left_valid and result.right_valid:
                    # Left failed, right passed
                    # Keep: right_fill (from abs_right_fill_start to abs_right_fill_end)
                    # Replace: left_fill + gap with 500N
                    # New sequence: 500N + right_fill

                    right_fill = seq[abs_right_fill_start:abs_right_fill_end]
                    new_fill = 'N' * 500 + right_fill

                    seq = seq[:abs_left_fill_start] + new_fill + seq[abs_right_fill_end:]

                    self.logger.info(f"    Partial revert {pf.gap_id}: kept right ({result.right_fill_length}bp), "
                                   f"reverted left to 500N")

                elif result.left_valid and not result.right_valid:
                    # Left passed, right failed
                    # Keep: left_fill (from abs_left_fill_start to abs_left_fill_end)
                    # Replace: gap + right_fill with 500N
                    # New sequence: left_fill + 500N

                    left_fill = seq[abs_left_fill_start:abs_left_fill_end]
                    new_fill = left_fill + 'N' * 500

                    seq = seq[:abs_left_fill_start] + new_fill + seq[abs_right_fill_end:]

                    self.logger.info(f"    Partial revert {pf.gap_id}: kept left ({result.left_fill_length}bp), "
                                   f"reverted right to 500N")

                # Update checkpoint with the new partial sequence
                new_sequence = seq[pf.filled_start:pf.filled_start + len(new_fill)]
                self.checkpoint.add_partial_gap(pf.gap_id, new_sequence)

            sequences[chrom] = seq

        # Write modified assembly
        with open(reverted_assembly, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        return reverted_assembly

    def _prepare_hic_data(self, assembly: Path):
        """Prepare Hi-C BAM file and analyzer"""
        hic_bam_path = self.hic_bam

        # Align Hi-C reads if BAM not provided
        if not hic_bam_path and self.hic_reads:
            hic_bam_path = self.output_dir / "hic_aligned.bam"
            if not hic_bam_path.exists():
                self.logger.info("  Aligning Hi-C reads...")
                align_hic_reads(
                    self.hic_reads[0],
                    self.hic_reads[1],
                    str(assembly),
                    str(hic_bam_path),
                    threads=self.threads
                )

        # Initialize analyzer
        if hic_bam_path and hic_bam_path.exists():
            self.hic_analyzer = HiCAnalyzer(
                hic_bam=str(hic_bam_path),
                assembly_file=str(assembly),
                threads=self.threads
            )
            self.logger.info(f"  Hi-C analyzer ready: {hic_bam_path}")
        else:
            self.logger.warning("  Hi-C BAM not available")

    def _find_gaps(self, assembly_file: Path) -> List[Dict]:
        """Find gaps (N-runs) in assembly, with caching"""
        cache_key = str(assembly_file)

        # Check cache first
        if cache_key in self._gap_cache:
            return self._gap_cache[cache_key]

        gaps = []

        for record in SeqIO.parse(assembly_file, 'fasta'):
            seq = str(record.seq).upper()

            for match in re.finditer(r'N+', seq):
                gap_size = match.end() - match.start()
                if gap_size >= self.min_gap_size:
                    gaps.append({
                        'chrom': record.id,
                        'start': match.start(),
                        'end': match.end(),
                        'name': f"{record.id}_{match.start()}_{match.end()}",
                        'size': gap_size
                    })

        # Cache the result
        self._gap_cache[cache_key] = gaps
        return gaps

    def _normalize_gaps(self, assembly_file: Path, gaps: List[Dict]) -> Path:
        """Normalize all gaps to 500N"""
        normalized_file = self.output_dir / "assembly_normalized.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly_file, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group gaps by chromosome
        gaps_by_chrom = {}
        for gap in gaps:
            chrom = gap['chrom']
            if chrom not in gaps_by_chrom:
                gaps_by_chrom[chrom] = []
            gaps_by_chrom[chrom].append(gap)

        # Normalize each chromosome
        for chrom, chrom_gaps in gaps_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]

            # Process gaps from end to start to avoid coordinate shift
            chrom_gaps_sorted = sorted(chrom_gaps, key=lambda x: x['start'], reverse=True)

            for gap in chrom_gaps_sorted:
                seq = seq[:gap['start']] + 'N' * 500 + seq[gap['end']:]

            sequences[chrom] = seq

        # Write normalized assembly
        with open(normalized_file, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        self.logger.info(f"Normalized {len(gaps)} gaps to 500N")
        return normalized_file

    def _align_reads(self, assembly: Path, output_dir: Path) -> tuple:
        """Align reads to assembly"""
        hifi_bam = None
        ont_bam = None

        if self.hifi_reads:
            hifi_bam = output_dir / "hifi.bam"
            self._run_minimap2(assembly, self.hifi_reads, hifi_bam, 'map-hifi')

        if self.ont_reads:
            ont_bam = output_dir / "ont.bam"
            self._run_minimap2(assembly, self.ont_reads, ont_bam, 'map-ont')

        return hifi_bam, ont_bam

    def _run_minimap2(self, ref: Path, reads: Path, output_bam: Path, preset: str):
        """Run minimap2 alignment"""
        try:
            # -m 2G: memory limit per thread for samtools sort
            cmd = f"minimap2 -ax {preset} -t {self.threads} {ref} {reads} | " \
                  f"samtools sort -@ {self.threads} -m 2G -o {output_bam} - && " \
                  f"samtools index {output_bam}"

            subprocess.run(cmd, shell=True, check=True,
                         capture_output=True, text=True)
            self.logger.info(f"  Created {output_bam.name}")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed: {e}")

    def _apply_fills(self, assembly: Path, fill_results: Dict, output_dir: Path) -> Path:
        """
        Apply gap fills to assembly and update pending fill coordinates.

        Fills are applied from end to start to avoid coordinate shifts affecting
        fills yet to be applied. After applying, we update the pending fills'
        coordinates to reflect their actual positions in the new assembly.
        """
        filled_assembly = output_dir / "assembly_filled.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group fills by chromosome with gap names for tracking
        fills_by_chrom = {}
        for gap_name, data in fill_results.items():
            if not data['result']['success']:
                continue
            gap = data['gap']
            chrom = gap['chrom']
            if chrom not in fills_by_chrom:
                fills_by_chrom[chrom] = []
            fills_by_chrom[chrom].append({
                'gap_name': gap_name,
                'start': gap['start'],
                'end': gap['end'],
                'sequence': data['result']['sequence'],
                'fill_length': len(data['result']['sequence'])
            })

        # Track coordinate shifts per chromosome for updating pending fills
        # Key: (chrom, original_start) -> new_start, new_end
        coordinate_updates = {}

        # Apply fills (from end to start) and track coordinate shifts
        for chrom, fills in fills_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]
            fills_sorted = sorted(fills, key=lambda x: x['start'], reverse=True)

            # Track cumulative offset for this chromosome
            cumulative_offset = 0

            for fill in fills_sorted:
                original_start = fill['start']
                original_end = fill['end']
                fill_seq = fill['sequence']
                fill_length = len(fill_seq)
                gap_length = original_end - original_start

                # Apply fill
                seq = seq[:original_start] + fill_seq + seq[original_end:]

                # Calculate actual position in new assembly
                # Since we're going from end to start, this fill's position is not affected
                # by previous fills in this loop (they were all after this one)
                actual_start = original_start
                actual_end = original_start + fill_length

                # Store the coordinate update
                coordinate_updates[(chrom, original_start, original_end)] = (actual_start, actual_end)

            sequences[chrom] = seq

        # Update pending fills with correct coordinates
        pending_fills = self.gap_tracker.get_pending_fills()
        for gap_id, pf in pending_fills.items():
            key = (pf.chrom, pf.original_start, pf.original_end)
            if key in coordinate_updates:
                new_start, new_end = coordinate_updates[key]
                # Update the pending fill's coordinates
                pf.filled_start = new_start
                pf.filled_end = new_end
                self.logger.debug(f"    Updated {gap_id} coordinates: {pf.original_start}-{pf.original_end} -> {new_start}-{new_end}")

        # Write filled assembly
        with open(filled_assembly, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

        return filled_assembly

    def _save_iteration_stats(self, output_dir: Path, iteration: int, fill_results: Dict):
        """Save iteration statistics including fill lengths"""
        successful_fills = []
        failed_fills = []
        total_bp_filled = 0

        for gap_name, data in fill_results.items():
            gap = data['gap']
            result = data['result']
            original_size = gap['end'] - gap['start']

            if result.get('success'):
                fill_length = len(result.get('sequence', ''))
                total_bp_filled += fill_length
                successful_fills.append({
                    'gap_name': gap_name,
                    'chrom': gap['chrom'],
                    'original_start': gap['start'],
                    'original_end': gap['end'],
                    'original_size': original_size,
                    'fill_length': fill_length,
                    'source': result.get('source', 'unknown'),
                    'tier': result.get('tier', 0),
                    'is_complete': result.get('is_complete', False),
                    'has_placeholder': result.get('has_placeholder', False)
                })
            else:
                failed_fills.append({
                    'gap_name': gap_name,
                    'chrom': gap['chrom'],
                    'original_size': original_size,
                    'reason': result.get('reason', 'unknown')
                })

        stats = {
            'iteration': iteration,
            'gaps_processed': len(fill_results),
            'gaps_filled': len(successful_fills),
            'gaps_failed': len(failed_fills),
            'total_bp_filled': total_bp_filled,
            'successful_fills': successful_fills,
            'failed_fills': failed_fills,
            'timestamp': datetime.now().isoformat()
        }

        with open(output_dir / "iteration_stats.json", 'w') as f:
            json.dump(stats, f, indent=2)

        # Log summary
        self.logger.info(f"  Iteration {iteration} summary: "
                        f"{len(successful_fills)} filled ({total_bp_filled:,} bp total), "
                        f"{len(failed_fills)} failed")

    def _save_final_stats(self):
        """Save final statistics"""
        self.stats['timestamp'] = datetime.now().isoformat()

        with open(self.output_dir / "final_stats.json", 'w') as f:
            json.dump(self.stats, f, indent=2)

        self.logger.info("\nFinal Statistics:")
        self.logger.info(f"  Iterations: {self.stats['iterations']}")
        self.logger.info(f"  Initial gaps: {self.stats['total_gaps_initial']}")
        self.logger.info(f"  Completely filled: {self.stats['gaps_completely_filled']}")
        self.logger.info(f"  Partially filled: {self.stats['gaps_partially_filled']}")
        self.logger.info(f"  Failed: {self.stats['gaps_failed']}")

    def _save_gap_tracker(self):
        """Save gap tracker state to file for resume support"""
        try:
            tracker_data = self.gap_tracker.to_dict()
            with open(self.gap_tracker_file, 'w') as f:
                json.dump(tracker_data, f, indent=2)
            self.logger.debug(f"Saved gap tracker state to {self.gap_tracker_file}")
        except Exception as e:
            self.logger.warning(f"Failed to save gap tracker: {e}")

    def _load_gap_tracker(self) -> bool:
        """Load gap tracker state from file. Returns True if loaded successfully."""
        if not self.gap_tracker_file.exists():
            return False
        try:
            with open(self.gap_tracker_file) as f:
                tracker_data = json.load(f)
            self.gap_tracker = GapStatusTracker.from_dict(tracker_data)
            pending_count = len(self.gap_tracker.get_pending_fills())
            self.logger.info(f"Loaded gap tracker state: {pending_count} pending fills")
            return True
        except Exception as e:
            self.logger.warning(f"Failed to load gap tracker: {e}")
            return False

    def _pre_assess_gaps(self, gaps: List[Dict], hifi_bam: Optional[Path], ont_bam: Optional[Path]) -> List[Dict]:
        """
        Pre-assess all gaps before first filling iteration.

        Analyzes flank quality for each gap and identifies those needing polish.

        Returns:
            List of gap dicts that need polishing, with 'polish_left' and 'polish_right' flags
        """
        gaps_needing_polish = []

        # Prefer HiFi BAM for assessment (more accurate)
        bam_file = None
        if hifi_bam and hifi_bam.exists():
            bam_file = str(hifi_bam)
        elif ont_bam and ont_bam.exists():
            bam_file = str(ont_bam)

        if not bam_file:
            self.logger.warning("  No BAM file available for pre-assessment, skipping")
            return gaps_needing_polish

        validator = GapValidator(threads=self.threads)

        try:
            results = validator.pre_assess_gaps(bam_file, gaps)

            needs_polish_count = 0
            ready_count = 0

            # Build gap lookup for adding polish flags
            gap_lookup = {g.get('name', f"{g['chrom']}_{g['start']}_{g['end']}"): g for g in gaps}

            for gap_name, result in results.items():
                if result.status == GapStatus.NEEDS_POLISH:
                    needs_polish_count += 1

                    # Create gap dict with polish flags
                    if gap_name in gap_lookup:
                        gap_info = gap_lookup[gap_name].copy()
                        gap_info['polish_left'] = result.left_flank_needs_polish
                        gap_info['polish_right'] = result.right_flank_needs_polish
                        gaps_needing_polish.append(gap_info)

                    self.logger.info(f"    {gap_name}: NEEDS_POLISH "
                                   f"(left={result.left_flank_needs_polish}, "
                                   f"right={result.right_flank_needs_polish})")
                else:
                    ready_count += 1

            self.logger.info(f"  Pre-assessment: {ready_count} ready, {needs_polish_count} need polish")

            # Update stats
            self.stats['gaps_needs_polish'] = needs_polish_count

        except Exception as e:
            self.logger.warning(f"  Pre-assessment error: {e}")
        finally:
            validator.close()

        return gaps_needing_polish

    def _polish_gap_flanks(self,
                           assembly: Path,
                           gaps_to_polish: List[Dict],
                           hifi_bam: Optional[Path],
                           ont_bam: Optional[Path],
                           output_dir: Path) -> Path:
        """
        Polish flanks for gaps that need it.

        Args:
            assembly: Current assembly file
            gaps_to_polish: List of gaps needing polish (with polish_left/polish_right flags)
            hifi_bam: HiFi BAM file (preferred)
            ont_bam: ONT BAM file (fallback)
            output_dir: Output directory

        Returns:
            Path to polished assembly (or original if no polish performed)
        """
        from gapfill.core.polisher import FlankPolisher

        # Prefer HiFi BAM for polishing
        bam_file = None
        if hifi_bam and hifi_bam.exists():
            bam_file = str(hifi_bam)
        elif ont_bam and ont_bam.exists():
            bam_file = str(ont_bam)

        if not bam_file:
            self.logger.warning("  No BAM file available for polishing, skipping")
            return assembly

        polished_assembly = output_dir / "assembly_polished.fasta"
        work_dir = output_dir / "polish_work"
        work_dir.mkdir(exist_ok=True)

        polisher = FlankPolisher(
            threads=self.threads,
            flank_size=1000,
            work_dir=work_dir
        )

        try:
            results = polisher.polish_assembly_flanks(
                str(assembly),
                bam_file,
                gaps_to_polish,
                str(polished_assembly)
            )

            # Count successful polishes
            polished_count = sum(1 for r in results.values() if r.success)

            if polished_count > 0:
                self.logger.info(f"  Polished {polished_count}/{len(gaps_to_polish)} gap flanks")

                # Reset status for polished gaps to PENDING (so they can be filled)
                for gap_name, result in results.items():
                    if result.success:
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.PENDING,
                            "Flanks polished, ready for filling"
                        )

                # Update stats
                self.stats['gaps_polished'] = polished_count

                return polished_assembly
            else:
                self.logger.info("  No gaps successfully polished, using original assembly")
                return assembly

        except Exception as e:
            self.logger.warning(f"  Polishing error: {e}")
            return assembly


# Backwards compatibility
IterativeGapFiller = HaploidEngine
