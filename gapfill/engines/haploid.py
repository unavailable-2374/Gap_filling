#!/usr/bin/env python3
"""
Haploid Gap Filling Engine

Iterative multi-round gap filling for haploid genomes.

Key features:
1. Gap normalization (all N placeholders -> 500N)
2. Re-generate BAM files each iteration
3. Correct coordinate update after each fill
4. Stop when all fillable gaps are completely filled
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
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads
from gapfill.utils.checkpoint import CheckpointManager, CheckpointState


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
                 clear_checkpoint: bool = False):

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
        self.logger.info("=" * 60)

    def run(self) -> Path:
        """Run iterative gap filling with checkpoint support"""
        current_assembly = self.initial_assembly
        start_iteration = 0

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
                    self.logger.info(f"  Completed gaps: {len(checkpoint_state.completed_gaps)}")
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
                    self.checkpoint.save(checkpoint_state)  # Save scanned state
                else:
                    self.logger.info("  No usable files found, starting fresh")
                    checkpoint_state = None
                self.logger.info("=" * 60)

        # Initialize checkpoint state if not resuming
        if not checkpoint_state:
            checkpoint_state = CheckpointState(
                engine="haploid",
                phase="init",
                max_iterations=self.max_iterations
            )
            self.checkpoint.save(checkpoint_state)

        # Find initial gaps
        initial_gaps = self._find_gaps(current_assembly)
        self.stats['total_gaps_initial'] = len(initial_gaps)
        self.logger.info(f"Initial gaps found: {len(initial_gaps)}")

        # Step 0a: Prepare Hi-C data if available
        if self.hic_reads or self.hic_bam:
            # Check if Hi-C BAM already exists
            existing_hic_bam = self.checkpoint.get_intermediate_file('hic_bam')
            if existing_hic_bam and self.resume:
                self.logger.info(f"STEP 0a: Reusing existing Hi-C BAM: {existing_hic_bam}")
                self.hic_bam = existing_hic_bam
            else:
                self.logger.info("STEP 0a: Preparing Hi-C data")
            self._prepare_hic_data(current_assembly)

            # Estimate gap sizes using Hi-C
            if self.hic_analyzer and initial_gaps:
                self.logger.info("STEP 0b: Estimating gap sizes with Hi-C")
                estimates = self.hic_analyzer.estimate_gap_sizes(initial_gaps)
                for est in estimates:
                    if est.confidence in ('high', 'medium'):
                        self.gap_size_estimates[est.gap_name] = est.estimated_size
                        self.logger.info(f"  {est.gap_name}: estimated {est.estimated_size}bp "
                                        f"(was {est.original_size}bp, {est.confidence})")

        # Step 0c: Normalize gaps (skip if already normalized by PolyploidEngine)
        if initial_gaps and not self.skip_normalization:
            # Check if normalized assembly already exists
            existing_normalized = self.checkpoint.get_intermediate_file('normalized_assembly')
            if existing_normalized and self.resume:
                self.logger.info(f"STEP 0c: Reusing normalized assembly: {existing_normalized}")
                current_assembly = existing_normalized
            else:
                self.logger.info("STEP 0c: Gap Normalization")
                current_assembly = self._normalize_gaps(current_assembly, initial_gaps)
                self.checkpoint.set_intermediate_file('normalized_assembly', str(current_assembly))
        elif self.skip_normalization:
            self.logger.info("STEP 0c: Skipping gap normalization (already normalized)")
            # Copy the assembly to output dir for consistency
            import shutil
            normalized_file = self.output_dir / "assembly_normalized.fasta"
            shutil.copy(current_assembly, normalized_file)
            current_assembly = normalized_file
            self.checkpoint.set_intermediate_file('normalized_assembly', str(current_assembly))

        # Update checkpoint phase
        checkpoint_state.phase = "filling"
        checkpoint_state.current_assembly = str(current_assembly)
        self.checkpoint.save(checkpoint_state)

        # Load gap states from checkpoint
        completely_filled_gaps = self.checkpoint.get_completed_gap_names()
        partially_filled_gaps = set(checkpoint_state.partially_filled_gaps.keys())
        failed_gaps = self.checkpoint.get_failed_gap_names()

        # Determine starting iteration
        if self.resume and checkpoint_state.iteration > 0:
            start_iteration = self.checkpoint.get_resume_iteration()
            self.logger.info(f"Resuming from iteration {start_iteration + 1}")

            # Load current assembly from checkpoint
            if checkpoint_state.current_assembly:
                checkpoint_assembly = Path(checkpoint_state.current_assembly)
                if checkpoint_assembly.exists():
                    current_assembly = checkpoint_assembly

        iteration = start_iteration
        while iteration < self.max_iterations:
            iteration += 1
            self.logger.info(f"\nITERATION {iteration}")

            # Update checkpoint
            checkpoint_state.iteration = iteration
            self.checkpoint.save(checkpoint_state)

            iter_dir = self.output_dir / f"iteration_{iteration}"
            iter_dir.mkdir(exist_ok=True)

            # Step 1: Align reads (check for existing BAMs)
            self.logger.info("Step 1: Aligning reads...")
            hifi_bam, ont_bam = self._align_reads_with_cache(current_assembly, iter_dir)

            if not hifi_bam and not ont_bam:
                self.logger.error("No BAM files generated")
                break

            # Step 2: Find gaps
            self.logger.info("Step 2: Finding gaps...")
            gaps = self._find_gaps(current_assembly)

            # Step 2b: Pre-assess gaps and polish if needed (only on first iteration)
            if iteration == 1 and not self.resume:
                self.logger.info("Step 2b: Pre-assessing gap flanks...")
                gaps_needing_polish = self._pre_assess_gaps(gaps, hifi_bam, ont_bam)

                # Step 2c: Polish flanks for gaps that need it
                if gaps_needing_polish:
                    self.logger.info("Step 2c: Polishing gap flanks...")
                    current_assembly = self._polish_gap_flanks(
                        current_assembly, gaps_needing_polish, hifi_bam, ont_bam, iter_dir
                    )
                    # Re-find gaps after polish (coordinates may have shifted slightly)
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

            # Step 3: Fill gaps
            self.logger.info("Step 3: Filling gaps...")
            work_dir = iter_dir / "work"
            work_dir.mkdir(exist_ok=True)

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
                gap_size_estimates=self.gap_size_estimates
            )

            fill_results = {}
            for gap in remaining_gaps:
                result = gap_filler.fill_gap(gap)
                fill_results[gap['name']] = {'gap': gap, 'result': result}

            gap_filler.close()

            # Step 4: Update assembly, checkpoint, and gap tracker
            new_completely_filled = 0
            new_partially_filled = 0
            new_validation_failed = 0
            new_unfillable = 0
            new_needs_polish = 0
            new_failed = 0

            # Track which fills passed validation (for applying to assembly)
            validated_fills = {}

            for gap_name, data in fill_results.items():
                result = data['result']
                validation = result.get('validation', {})
                validation_status = validation.get('status', '')
                validation_valid = validation.get('valid', True)  # Default True if no validation

                if result['success']:
                    # CONSERVATIVE: Check validation result before accepting fill
                    if not validation_valid:
                        # Assembly succeeded but validation failed → reject fill
                        new_validation_failed += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FAILED,
                            f"Validation failed: {validation.get('reason', 'Unknown')}"
                        )
                        self.logger.warning(f"  Gap {gap_name}: fill rejected (validation failed)")
                        # Don't add to validated_fills - will retry next iteration
                        continue

                    # Validation passed - accept the fill
                    validated_fills[gap_name] = data

                    if result.get('is_complete', False) and not result.get('has_placeholder', False):
                        completely_filled_gaps.add(gap_name)
                        new_completely_filled += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FILLED_COMPLETE,
                            validation.get('reason', 'Fill successful')
                        )
                        # Save to checkpoint
                        self.checkpoint.add_completed_gap(gap_name, result.get('sequence', ''))
                    else:
                        partially_filled_gaps.add(gap_name)
                        new_partially_filled += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FILLED_PARTIAL,
                            validation.get('reason', 'Partial fill')
                        )
                        self.checkpoint.add_partial_gap(gap_name, result.get('sequence', ''))
                else:
                    # Check validation status for failed gaps
                    if validation_status == 'unfillable':
                        new_unfillable += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.UNFILLABLE,
                            validation.get('reason', 'No spanning reads available')
                        )
                    elif validation_status == 'needs_polish':
                        new_needs_polish += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.NEEDS_POLISH,
                            validation.get('reason', 'Flanks need polishing')
                        )
                    else:
                        new_failed += 1
                        self.gap_tracker.set_status(
                            gap_name, GapStatus.FAILED,
                            validation.get('reason', result.get('reason', 'Unknown failure'))
                        )

                    failed_gaps.add(gap_name)
                    self.checkpoint.add_failed_gap(gap_name)

            self.logger.info(f"  Complete: {new_completely_filled}, Partial: {new_partially_filled}, "
                           f"ValidationFailed: {new_validation_failed}, "
                           f"Unfillable: {new_unfillable}, NeedsPolish: {new_needs_polish}, Failed: {new_failed}")

            # Apply only validated fills to assembly
            current_assembly = self._apply_fills(
                current_assembly, validated_fills, iter_dir
            )

            # Update checkpoint with current assembly
            self.checkpoint.set_current_assembly(str(current_assembly))

            # Save iteration stats
            self._save_iteration_stats(iter_dir, iteration, fill_results)

            # Check for progress
            if new_completely_filled == 0 and new_partially_filled == 0:
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
                gap for gap in initial_gaps
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

        # Mark checkpoint as complete
        self.checkpoint.mark_complete()
        self.logger.info("Checkpoint marked as complete")

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
        """Find gaps (N-runs) in assembly"""
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
        """Apply gap fills to assembly"""
        filled_assembly = output_dir / "assembly_filled.fasta"

        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly, 'fasta'):
            sequences[record.id] = str(record.seq)

        # Group fills by chromosome
        fills_by_chrom = {}
        for gap_name, data in fill_results.items():
            if not data['result']['success']:
                continue
            gap = data['gap']
            chrom = gap['chrom']
            if chrom not in fills_by_chrom:
                fills_by_chrom[chrom] = []
            fills_by_chrom[chrom].append({
                'start': gap['start'],
                'end': gap['end'],
                'sequence': data['result']['sequence']
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
        """Save iteration statistics"""
        stats = {
            'iteration': iteration,
            'gaps_processed': len(fill_results),
            'timestamp': datetime.now().isoformat()
        }

        with open(output_dir / "iteration_stats.json", 'w') as f:
            json.dump(stats, f, indent=2)

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
