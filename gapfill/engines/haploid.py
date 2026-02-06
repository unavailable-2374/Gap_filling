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
from gapfill.core.validator import GapValidator, GapStatusTracker, GapStatus
from gapfill.core.parallel import fill_gaps_parallel, fill_gaps_sequential, GapBatcher
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

        self.stats = {
            'iterations': 0,
            'total_gaps_initial': 0,
            'gaps_completely_filled': 0,
            'gaps_partially_filled': 0,
            'gaps_failed': 0,
            'total_bp_filled': 0
        }

        self.logger.info("=" * 60)
        self.logger.info("HaploidEngine initialized")
        self.logger.info(f"  Assembly: {self.initial_assembly}")
        self.logger.info(f"  HiFi reads: {self.hifi_reads}")
        self.logger.info(f"  ONT reads: {self.ont_reads}")
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
                    if filled_asm.exists():
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
            gaps = self._find_gaps(current_assembly)
            self.stats['total_gaps_initial'] = len(gaps)

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
                filled_asm = iter_dir / "assembly_filled.fasta"
                if filled_asm.exists():
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

            # Step 2: Find gaps
            self.logger.info("Step 2: Finding gaps...")
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

            # Step 3: Fill gaps (with local validation inside each fill)
            if not _is_done(iter_dir, "filling"):
                self.logger.info("Step 3: Filling gaps (with local validation)...")
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
                        enable_consensus_first=True
                    )

                    fill_results = {}
                    for gap in remaining_gaps:
                        result = gap_filler.fill_gap(gap)
                        fill_results[gap['name']] = {'gap': gap, 'result': result}

                    gap_filler.close()

                # Step 4: Process fill results (all locally validated)
                new_complete = 0
                new_partial = 0
                new_unfillable = 0
                new_failed = 0

                # Track successful fills to apply to assembly
                successful_fills = {}

                for gap_name, data in fill_results.items():
                    gap = data['gap']
                    result = data['result']

                    if result.get('success') and result.get('validated', False):
                        # Fill succeeded and passed local validation
                        successful_fills[gap_name] = data
                        is_complete = result.get('is_complete', False) and not result.get('has_placeholder', False)

                        fill_seq = result.get('sequence', '')
                        fill_length = len(fill_seq)
                        original_gap_size = gap['end'] - gap['start']

                        if is_complete:
                            new_complete += 1
                            self.gap_tracker.set_status(
                                gap_name, GapStatus.FILLED_COMPLETE,
                                f"Locally validated: {result.get('source', 'unknown')} "
                                f"(tier={result.get('tier', 0)})"
                            )
                            self.logger.info(f"    {gap_name}: ✓ filled & validated "
                                           f"({result.get('source', 'unknown')}, "
                                           f"{original_gap_size}bp -> {fill_length}bp)")
                        else:
                            # Handle partial fill — truncate failed sides
                            partial_val = result.get('validation', {})
                            if not partial_val.get('left_valid', True) or not partial_val.get('right_valid', True):
                                fill_seq = self._truncate_partial_fill(fill_seq, partial_val)
                                result['sequence'] = fill_seq
                                data['result'] = result
                            new_partial += 1
                            self.gap_tracker.set_status(
                                gap_name, GapStatus.FILLED_PARTIAL,
                                f"Locally validated (partial): {result.get('source', 'unknown')} "
                                f"(tier={result.get('tier', 0)})"
                            )
                            self.logger.info(f"    {gap_name}: ✓ partially filled & validated "
                                           f"({result.get('source', 'unknown')}, "
                                           f"{original_gap_size}bp -> {fill_length}bp)")

                        self.checkpoint.add_completed_gap(gap_name, result.get('sequence', ''))

                    elif result.get('success'):
                        # Assembled but all tiers failed validation
                        new_failed += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FAILED,
                            "All tiers failed local validation"
                        )
                        failed_gaps.add(gap_name)
                        self.checkpoint.add_failed_gap(gap_name)

                    else:
                        # Fill failed entirely → analyze why
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

                self.logger.info(f"  Complete fills: {new_complete}, "
                               f"Partial fills: {new_partial}, "
                               f"Unfillable: {new_unfillable}, Failed: {new_failed}")

                # Save fill results for potential resume
                self._save_iteration_stats(iter_dir, iteration, fill_results)
                self._save_gap_tracker()
                _mark_done(iter_dir, "filling")
            else:
                self.logger.info("Step 3-4: Filling already done")
                # If filling is done but apply is not done, this is an inconsistent state.
                if not _is_done(iter_dir, "apply"):
                    self.logger.warning("  Apply step not done but filling done - inconsistent state")
                    self.logger.warning("  Re-running filling to get actual sequences...")
                    _clear_marker(iter_dir, "filling")
                    continue
                successful_fills = {}
                new_complete = 0
                new_partial = 0

            # Step 5: Apply fills to assembly
            if not _is_done(iter_dir, "apply"):
                if successful_fills:
                    self.logger.info("Step 5: Applying fills to assembly...")
                    current_assembly = self._apply_fills(
                        current_assembly, successful_fills, iter_dir
                    )
                _mark_done(iter_dir, "apply")
            else:
                self.logger.info("Step 5: Apply already done, loading assembly...")
                filled_asm = iter_dir / "assembly_filled.fasta"
                if filled_asm.exists():
                    current_assembly = filled_asm

            # Log status summary
            tracker_summary = self.gap_tracker.get_summary()
            self.logger.info(f"  Status summary: {tracker_summary}")

            # Mark iteration complete
            self._save_gap_tracker()
            _mark_done(iter_dir, "iteration")

            # Check for progress
            if new_complete + new_partial == 0:
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

        self._save_final_stats()
        self._save_gap_tracker()

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
        1. Gap normalization
        2. Initial alignment
        3. Filter reads (if optimized mode)
        4. Pre-assess gap flanks
        5. Polish problematic flanks

        Returns:
            Path to preprocessed assembly (normalized and polished)
        """
        current_assembly = assembly
        preprocess_dir = self.output_dir / "preprocessing"
        preprocess_dir.mkdir(exist_ok=True)

        # Step 0a: Find initial gaps
        initial_gaps = self._find_gaps(current_assembly)
        self.stats['total_gaps_initial'] = len(initial_gaps)
        self.logger.info(f"Initial gaps found: {len(initial_gaps)}")

        # Step 0b: Normalize gaps
        normalized_file = self.output_dir / "assembly_normalized.fasta"
        if _is_done(self.output_dir, "normalized") and normalized_file.exists():
            self.logger.info(f"STEP 0b: Reusing normalized assembly")
            current_assembly = normalized_file
        elif initial_gaps and not self.skip_normalization:
            self.logger.info("STEP 0b: Normalizing gaps to 500N...")
            current_assembly = self._normalize_gaps(current_assembly, initial_gaps)
            _mark_done(self.output_dir, "normalized")
        elif self.skip_normalization:
            self.logger.info("STEP 0b: Skipping normalization (already done)")
            import shutil
            shutil.copy(current_assembly, normalized_file)
            current_assembly = normalized_file
            _mark_done(self.output_dir, "normalized")

        # Step 0c: Initial alignment for pre-assessment
        if not _is_done(preprocess_dir, "alignment"):
            self.logger.info("STEP 0c: Initial alignment for pre-assessment...")
            hifi_bam, ont_bam = self._align_reads_with_cache(current_assembly, preprocess_dir)
            _mark_done(preprocess_dir, "alignment")
        else:
            self.logger.info("STEP 0c: Reusing existing alignment...")
            hifi_bam = preprocess_dir / "hifi.bam" if (preprocess_dir / "hifi.bam").exists() else None
            ont_bam = preprocess_dir / "ont.bam" if (preprocess_dir / "ont.bam").exists() else None

        # Track preprocessing BAM paths for potential reuse in iteration 1
        self._preprocessing_hifi_bam = hifi_bam
        self._preprocessing_ont_bam = ont_bam
        self._preprocessing_assembly = current_assembly

        if not hifi_bam and not ont_bam:
            self.logger.warning("  No BAM files for pre-assessment, skipping polish")
            return current_assembly

        # Step 0d: Filter reads (OPTIMIZATION - filter once, use many times)
        if self.optimized_mode and self.reads_cache:
            if not _is_done(preprocess_dir, "filter"):
                self.logger.info("STEP 0d: Filtering reads (keeping gap-related reads only)...")
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
                self.logger.info("STEP 0d: Reads already filtered")
                # Re-initialize reads_cache with existing gaps
                gaps = self._find_gaps(current_assembly)
                self.reads_cache.set_gap_regions(gaps)

        # Step 0e: Pre-assess gap flanks
        polished_file = self.output_dir / "assembly_polished.fasta"
        polishing_occurred = False

        if not _is_done(preprocess_dir, "polish"):
            self.logger.info("STEP 0e: Pre-assessing gap flanks...")
            gaps = self._find_gaps(current_assembly)
            gaps_needing_polish = self._pre_assess_gaps(gaps, hifi_bam, ont_bam)

            # Step 0f: Polish problematic flanks
            if gaps_needing_polish:
                self.logger.info(f"STEP 0f: Polishing {len(gaps_needing_polish)} gap flanks...")
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
            self.logger.info("STEP 0e-0f: Polish already done")
            if polished_file.exists():
                current_assembly = polished_file
                polishing_occurred = True

        # Mark whether preprocessing BAM can be reused for iteration 1
        self._can_reuse_preprocessing_bam = not polishing_occurred

        return current_assembly

    def _truncate_partial_fill(self, fill_seq: str, validation: Dict) -> str:
        """
        Truncate a partial fill string by removing the side(s) that failed validation.

        Operates on the fill string BEFORE it is applied to the assembly.
        """
        n_match = re.search(r'N{10,}', fill_seq)
        if not n_match:
            return fill_seq

        left = fill_seq[:n_match.start()]
        gap = fill_seq[n_match.start():n_match.end()]
        right = fill_seq[n_match.end():]

        if not validation.get('left_valid', True):
            left = ''
        if not validation.get('right_valid', True):
            right = ''

        return left + gap + right

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
        Apply locally-validated gap fills to assembly.

        Fills are applied from end to start to avoid coordinate shifts affecting
        fills yet to be applied.
        """
        filled_assembly = output_dir / "assembly_filled.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group fills by chromosome
        fills_by_chrom = {}
        for gap_name, data in fill_results.items():
            if not data['result'].get('success'):
                continue
            gap = data['gap']
            chrom = gap['chrom']
            if chrom not in fills_by_chrom:
                fills_by_chrom[chrom] = []
            fills_by_chrom[chrom].append({
                'start': gap['start'],
                'end': gap['end'],
                'sequence': data['result']['sequence'],
            })

        # Apply fills (from end to start)
        for chrom, fills in fills_by_chrom.items():
            if chrom not in sequences:
                continue

            seq = sequences[chrom]
            fills_sorted = sorted(fills, key=lambda x: x['start'], reverse=True)

            for fill in fills_sorted:
                seq = seq[:fill['start']] + fill['sequence'] + seq[fill['end']:]

            sequences[chrom] = seq

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
            summary = self.gap_tracker.get_summary()
            self.logger.info(f"Loaded gap tracker state: {summary}")
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
