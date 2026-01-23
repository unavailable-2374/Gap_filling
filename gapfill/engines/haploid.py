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
from gapfill.core.validator import GapValidator
from gapfill.utils.hic import HiCAnalyzer, align_hic_reads


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
                 min_mapq: int = 20):

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

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Hi-C analyzer (initialized when needed)
        self.hic_analyzer: Optional[HiCAnalyzer] = None
        self.gap_size_estimates: Dict[str, int] = {}

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
        self.logger.info("=" * 60)

    def run(self) -> Path:
        """Run iterative gap filling"""
        current_assembly = self.initial_assembly
        iteration = 0

        # Find initial gaps
        initial_gaps = self._find_gaps(current_assembly)
        self.stats['total_gaps_initial'] = len(initial_gaps)
        self.logger.info(f"Initial gaps found: {len(initial_gaps)}")

        # Step 0a: Prepare Hi-C data if available
        if self.hic_reads or self.hic_bam:
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

        # Step 0c: Normalize gaps
        if initial_gaps:
            self.logger.info("STEP 0c: Gap Normalization")
            current_assembly = self._normalize_gaps(current_assembly, initial_gaps)

        completely_filled_gaps = set()
        partially_filled_gaps = set()
        failed_gaps = set()

        while iteration < self.max_iterations:
            iteration += 1
            self.logger.info(f"\nITERATION {iteration}")

            iter_dir = self.output_dir / f"iteration_{iteration}"
            iter_dir.mkdir(exist_ok=True)

            # Step 1: Align reads
            self.logger.info("Step 1: Aligning reads...")
            hifi_bam, ont_bam = self._align_reads(current_assembly, iter_dir)

            if not hifi_bam and not ont_bam:
                self.logger.error("No BAM files generated")
                break

            # Step 2: Find gaps
            self.logger.info("Step 2: Finding gaps...")
            gaps = self._find_gaps(current_assembly)

            remaining_gaps = [g for g in gaps
                            if g['name'] not in completely_filled_gaps
                            and g['name'] not in failed_gaps]

            self.logger.info(f"  Remaining to process: {len(remaining_gaps)}")

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
                min_mapq=self.min_mapq
            )

            fill_results = {}
            for gap in remaining_gaps:
                result = gap_filler.fill_gap(gap)
                fill_results[gap['name']] = {'gap': gap, 'result': result}

            gap_filler.close()

            # Step 4: Update assembly
            new_completely_filled = 0
            new_partially_filled = 0
            new_failed = 0

            for gap_name, data in fill_results.items():
                result = data['result']
                if result['success']:
                    if result.get('is_complete', False) and not result.get('has_placeholder', False):
                        completely_filled_gaps.add(gap_name)
                        new_completely_filled += 1
                    else:
                        partially_filled_gaps.add(gap_name)
                        new_partially_filled += 1
                else:
                    failed_gaps.add(gap_name)
                    new_failed += 1

            self.logger.info(f"  Complete: {new_completely_filled}, Partial: {new_partially_filled}, Failed: {new_failed}")

            # Apply fills to assembly
            current_assembly = self._apply_fills(
                current_assembly, fill_results, iter_dir
            )

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

        self.logger.info(f"\nFinal assembly: {final_assembly}")
        return final_assembly

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
            cmd = f"minimap2 -ax {preset} -t {self.threads} {ref} {reads} | " \
                  f"samtools sort -@ {self.threads} -o {output_bam} - && " \
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


# Backwards compatibility
IterativeGapFiller = HaploidEngine
