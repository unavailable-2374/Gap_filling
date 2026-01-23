#!/usr/bin/env python3
"""
Iterative Gap Filler V2 - Simplified multi-round gap filling

Key features:
1. No mini-reference support (always use full genome)
2. Re-generate BAM files each iteration (no BAM reuse)
3. Correct coordinate update after each fill
4. Stop when all fillable gaps are completely filled by spanning reads

Workflow per iteration:
1. Align reads to current assembly -> generate new BAM
2. Scan assembly for remaining gaps (N regions)
3. For each gap:
   - Try spanning reads -> complete fill
   - Try flanking reads -> partial fill with 500N placeholder
   - No reads -> mark as failed
4. Update assembly and coordinates
5. Continue until no more progress or max iterations reached

Author: Gap Filling Pipeline
"""

import logging
import subprocess
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from gap_filler import GapFiller
from gap_validation import GapValidator


class IterativeGapFiller:
    """
    Iterative gap filler with simplified workflow
    """

    def __init__(self,
                 assembly_file: str,
                 hifi_reads: Optional[str] = None,
                 ont_reads: Optional[str] = None,
                 output_dir: str = "output",
                 threads: int = 8,
                 max_iterations: int = 10,
                 min_gap_size: int = 100,
                 min_mapq: int = 20):
        """
        Initialize IterativeGapFiller

        Args:
            assembly_file: Path to initial assembly FASTA
            hifi_reads: Path to HiFi reads (FASTQ/FASTA, can be gzipped)
            ont_reads: Path to ONT reads (FASTQ/FASTA, can be gzipped)
            output_dir: Output directory
            threads: Number of threads
            max_iterations: Maximum number of iterations
            min_gap_size: Minimum gap size to process (default 100)
            min_mapq: Minimum mapping quality for reads (default 20)
        """
        self.initial_assembly = Path(assembly_file)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.max_iterations = max_iterations
        self.min_gap_size = min_gap_size
        self.min_mapq = min_mapq

        self.logger = logging.getLogger(__name__)

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Track statistics
        self.stats = {
            'iterations': 0,
            'total_gaps_initial': 0,
            'gaps_completely_filled': 0,
            'gaps_partially_filled': 0,
            'gaps_failed': 0,
            'total_bp_filled': 0
        }

        self.logger.info("=" * 60)
        self.logger.info("IterativeGapFiller initialized")
        self.logger.info(f"  Assembly: {self.initial_assembly}")
        self.logger.info(f"  HiFi reads: {self.hifi_reads}")
        self.logger.info(f"  ONT reads: {self.ont_reads}")
        self.logger.info(f"  Output dir: {self.output_dir}")
        self.logger.info(f"  Max iterations: {self.max_iterations}")
        self.logger.info("=" * 60)

    def run(self) -> Path:
        """
        Run iterative gap filling

        Returns:
            Path to final filled assembly
        """
        current_assembly = self.initial_assembly
        iteration = 0

        # Find initial gaps
        initial_gaps = self._find_gaps(current_assembly)
        self.stats['total_gaps_initial'] = len(initial_gaps)
        self.logger.info(f"Initial gaps found: {len(initial_gaps)}")

        # =================================================================
        # STEP 0: Normalize all gaps to 500N
        # This is critical because N length is just a placeholder and does
        # not represent the actual gap size. Normalizing to 500N ensures
        # spanning reads detection works correctly.
        # =================================================================
        if initial_gaps:
            self.logger.info("=" * 60)
            self.logger.info("STEP 0: Gap Normalization")
            self.logger.info("=" * 60)
            current_assembly = self._normalize_gaps(current_assembly, initial_gaps)
            self.logger.info(f"Normalized assembly: {current_assembly}")
            self.logger.info("")

        # Track filled gaps across iterations
        completely_filled_gaps = set()  # Gaps that were completely filled (no placeholder)
        partially_filled_gaps = set()   # Gaps that still have 500N placeholder
        failed_gaps = set()             # Gaps that couldn't be filled

        while iteration < self.max_iterations:
            iteration += 1
            self.logger.info("=" * 60)
            self.logger.info(f"ITERATION {iteration}")
            self.logger.info("=" * 60)

            # Create iteration directory
            iter_dir = self.output_dir / f"iteration_{iteration}"
            iter_dir.mkdir(exist_ok=True)

            # =================================================================
            # Step 1: Align reads to current assembly (ALWAYS re-generate BAM)
            # =================================================================
            self.logger.info("Step 1: Aligning reads to assembly...")
            hifi_bam, ont_bam = self._align_reads(current_assembly, iter_dir)

            if not hifi_bam and not ont_bam:
                self.logger.error("No BAM files generated, cannot proceed")
                break

            # =================================================================
            # Step 2: Find remaining gaps in current assembly
            # =================================================================
            self.logger.info("Step 2: Finding gaps in current assembly...")
            gaps = self._find_gaps(current_assembly)

            # Filter out gaps that are already completely filled or failed
            remaining_gaps = [g for g in gaps
                            if g['name'] not in completely_filled_gaps
                            and g['name'] not in failed_gaps]

            self.logger.info(f"  Total gaps: {len(gaps)}")
            self.logger.info(f"  Remaining to process: {len(remaining_gaps)}")
            self.logger.info(f"  Already complete: {len(completely_filled_gaps)}")
            self.logger.info(f"  Failed: {len(failed_gaps)}")

            if not remaining_gaps:
                self.logger.info("No remaining gaps to process, stopping")
                break

            # =================================================================
            # Step 3: Fill gaps
            # =================================================================
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

            # Process each gap and collect results
            fill_results = {}
            for gap in remaining_gaps:
                result = gap_filler.fill_gap(gap)
                fill_results[gap['name']] = {
                    'gap': gap,
                    'result': result
                }

            gap_filler.close()

            # =================================================================
            # Step 4: Apply fills and update assembly
            # =================================================================
            self.logger.info("Step 4: Applying fills to assembly...")

            new_assembly_path = iter_dir / "assembly_filled.fasta"
            new_completely_filled = set()
            new_partially_filled = set()
            new_failed = set()
            bp_filled_this_iter = 0

            # Load current assembly
            assembly_seqs = {}
            for record in SeqIO.parse(current_assembly, 'fasta'):
                assembly_seqs[record.id] = str(record.seq)

            # Sort gaps by position (descending) to apply from end to start
            # This prevents coordinate shift issues
            sorted_gaps = sorted(
                remaining_gaps,
                key=lambda g: (g['chrom'], g['start']),
                reverse=True
            )

            for gap in sorted_gaps:
                gap_name = gap['name']
                chrom = gap['chrom']
                gap_start = gap['start']
                gap_end = gap['end']

                if gap_name not in fill_results:
                    continue

                result = fill_results[gap_name]['result']

                if result['success']:
                    seq = assembly_seqs[chrom]
                    new_seq = seq[:gap_start] + result['sequence'] + seq[gap_end:]
                    assembly_seqs[chrom] = new_seq

                    bp_filled = len(result['sequence']) - (gap_end - gap_start)
                    bp_filled_this_iter += abs(bp_filled)

                    if result.get('is_complete', False):
                        # Completely filled by spanning reads
                        new_completely_filled.add(gap_name)
                        self.logger.info(f"  ✓ {gap_name}: COMPLETE ({len(result['sequence'])}bp)")
                    else:
                        # Partially filled with 500N placeholder
                        new_partially_filled.add(gap_name)
                        self.logger.info(f"  ~ {gap_name}: PARTIAL ({len(result['sequence'])}bp, has 500N)")
                else:
                    new_failed.add(gap_name)
                    self.logger.info(f"  ✗ {gap_name}: FAILED ({result.get('reason', 'unknown')})")

            # Write new assembly
            records = []
            for chrom, seq in assembly_seqs.items():
                records.append(SeqRecord(Seq(seq), id=chrom, description=''))
            SeqIO.write(records, new_assembly_path, 'fasta')

            # Update tracking sets
            completely_filled_gaps.update(new_completely_filled)
            partially_filled_gaps.update(new_partially_filled)
            partially_filled_gaps -= new_completely_filled  # Remove if now complete
            failed_gaps.update(new_failed)

            # Update stats
            self.stats['gaps_completely_filled'] = len(completely_filled_gaps)
            self.stats['gaps_partially_filled'] = len(partially_filled_gaps)
            self.stats['gaps_failed'] = len(failed_gaps)
            self.stats['total_bp_filled'] += bp_filled_this_iter

            # =================================================================
            # Step 5: Check if we should continue
            # =================================================================
            self.logger.info(f"Iteration {iteration} summary:")
            self.logger.info(f"  Newly complete: {len(new_completely_filled)}")
            self.logger.info(f"  Newly partial: {len(new_partially_filled)}")
            self.logger.info(f"  Newly failed: {len(new_failed)}")
            self.logger.info(f"  Total complete: {len(completely_filled_gaps)}")
            self.logger.info(f"  Total partial (need more work): {len(partially_filled_gaps)}")
            self.logger.info(f"  Total failed: {len(failed_gaps)}")

            # Check termination conditions
            # Continue if there are partially filled gaps (they might become complete next round)
            if not new_completely_filled and not new_partially_filled:
                self.logger.info("No progress made this iteration, stopping")
                break

            if not partially_filled_gaps:
                self.logger.info("All fillable gaps are completely filled, stopping")
                break

            # Update current assembly for next iteration
            current_assembly = new_assembly_path

            # Save iteration statistics
            iter_stats = {
                'iteration': iteration,
                'new_complete': list(new_completely_filled),
                'new_partial': list(new_partially_filled),
                'new_failed': list(new_failed),
                'total_complete': len(completely_filled_gaps),
                'total_partial': len(partially_filled_gaps),
                'total_failed': len(failed_gaps)
            }
            with open(iter_dir / "iteration_stats.json", 'w') as f:
                json.dump(iter_stats, f, indent=2)

        # =================================================================
        # Final output
        # =================================================================
        self.stats['iterations'] = iteration

        # Copy final assembly to output directory
        final_assembly = self.output_dir / "final_assembly.fasta"
        import shutil
        shutil.copy(current_assembly, final_assembly)

        # Run validation on final assembly
        self.logger.info("Running validation on final assembly...")
        self._run_validation(final_assembly)

        # Save final statistics
        self._save_final_stats()

        self.logger.info("=" * 60)
        self.logger.info("ITERATIVE GAP FILLING COMPLETE")
        self.logger.info(f"  Total iterations: {self.stats['iterations']}")
        self.logger.info(f"  Initial gaps: {self.stats['total_gaps_initial']}")
        self.logger.info(f"  Completely filled: {self.stats['gaps_completely_filled']}")
        self.logger.info(f"  Partially filled: {self.stats['gaps_partially_filled']}")
        self.logger.info(f"  Failed: {self.stats['gaps_failed']}")
        self.logger.info(f"  Total bp filled: {self.stats['total_bp_filled']:,}")
        self.logger.info(f"  Final assembly: {final_assembly}")
        self.logger.info("=" * 60)

        return final_assembly

    def _align_reads(self, assembly: Path, output_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Align reads to assembly using minimap2

        IMPORTANT: Always generates new BAM files, never reuses

        Args:
            assembly: Path to assembly FASTA
            output_dir: Output directory for BAM files

        Returns:
            Tuple of (hifi_bam_path, ont_bam_path)
        """
        hifi_bam = None
        ont_bam = None

        # Align HiFi reads
        if self.hifi_reads and self.hifi_reads.exists():
            hifi_bam = output_dir / "hifi.bam"
            self.logger.info(f"  Aligning HiFi reads to {assembly.name}...")

            success = self._run_minimap2(
                assembly, self.hifi_reads, hifi_bam, 'hifi'
            )
            if not success:
                self.logger.warning("HiFi alignment failed")
                hifi_bam = None

        # Align ONT reads
        if self.ont_reads and self.ont_reads.exists():
            ont_bam = output_dir / "ont.bam"
            self.logger.info(f"  Aligning ONT reads to {assembly.name}...")

            success = self._run_minimap2(
                assembly, self.ont_reads, ont_bam, 'ont'
            )
            if not success:
                self.logger.warning("ONT alignment failed")
                ont_bam = None

        return hifi_bam, ont_bam

    def _run_minimap2(self, assembly: Path, reads: Path, output_bam: Path,
                      read_type: str) -> bool:
        """
        Run minimap2 alignment

        Args:
            assembly: Path to assembly FASTA
            reads: Path to reads file
            output_bam: Output BAM path
            read_type: 'hifi' or 'ont'

        Returns:
            True if successful
        """
        preset = 'map-hifi' if read_type == 'hifi' else 'map-ont'

        try:
            # minimap2 | samtools sort
            minimap2_cmd = [
                'minimap2',
                '-ax', preset,
                '-t', str(self.threads),
                '--secondary=no',
                str(assembly),
                str(reads)
            ]

            samtools_cmd = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(output_bam),
                '-'
            ]

            minimap2_proc = subprocess.Popen(
                minimap2_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            samtools_proc = subprocess.Popen(
                samtools_cmd,
                stdin=minimap2_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            minimap2_proc.stdout.close()
            _, stderr = samtools_proc.communicate()

            if samtools_proc.returncode != 0:
                self.logger.error(f"Alignment failed: {stderr.decode()[:500]}")
                return False

            # Index BAM
            subprocess.run(
                ['samtools', 'index', str(output_bam)],
                check=True
            )

            return True

        except Exception as e:
            self.logger.error(f"Alignment error: {e}")
            return False

    def _find_gaps(self, assembly: Path) -> List[Dict]:
        """
        Find gaps (N regions) in assembly

        Args:
            assembly: Path to assembly FASTA

        Returns:
            List of gap dictionaries
        """
        gaps = []
        gap_id = 0

        for record in SeqIO.parse(assembly, 'fasta'):
            chrom = record.id
            seq = str(record.seq).upper()

            # Find all N runs
            for match in re.finditer(r'N+', seq):
                start = match.start()
                end = match.end()
                size = end - start

                if size >= self.min_gap_size:
                    gap_id += 1
                    gaps.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'size': size,
                        'name': f"gap_{chrom}_{start}_{end}"
                    })

        self.logger.info(f"  Found {len(gaps)} gaps >= {self.min_gap_size}bp")
        return gaps

    def _normalize_gaps(self, assembly: Path, gaps: List[Dict], normalized_size: int = 500) -> Path:
        """
        Normalize all gaps to fixed size (500N)

        IMPORTANT: The N's in the assembly are just placeholders - their length
        does NOT represent the actual biological gap size. A 72kb N region and
        a 500bp N region both just mean "there's a gap here". Normalizing to 500N
        ensures that spanning reads detection works correctly regardless of the
        original placeholder size.

        Args:
            assembly: Input assembly path
            gaps: List of gap dictionaries from _find_gaps()
            normalized_size: Target gap size (default: 500)

        Returns:
            Path to normalized assembly
        """
        self.logger.info(f"Normalizing {len(gaps)} gaps to {normalized_size}N markers...")
        self.logger.info("(N's are placeholders - real gap sizes are determined by reads)")

        # Group gaps by chromosome
        gaps_by_chrom = {}
        for gap in gaps:
            chrom = gap['chrom']
            if chrom not in gaps_by_chrom:
                gaps_by_chrom[chrom] = []
            gaps_by_chrom[chrom].append(gap)

        # Read assembly and normalize gaps
        normalized_records = []
        total_normalized = 0

        for record in SeqIO.parse(assembly, 'fasta'):
            chrom = record.id
            seq = str(record.seq)

            # Get gaps for this chromosome
            chrom_gaps = gaps_by_chrom.get(chrom, [])

            if not chrom_gaps:
                # No gaps in this chromosome, keep as is
                normalized_records.append(record)
                continue

            # Sort gaps by position (reverse order to avoid coordinate shift)
            chrom_gaps_sorted = sorted(chrom_gaps, key=lambda x: x['start'], reverse=True)

            # Replace each gap with normalized size
            for gap in chrom_gaps_sorted:
                gap_start = gap['start']
                gap_end = gap['end']
                original_size = gap_end - gap_start

                if original_size != normalized_size:
                    # Replace with normalized N's
                    seq = seq[:gap_start] + 'N' * normalized_size + seq[gap_end:]
                    self.logger.debug(
                        f"  {gap['name']}: {original_size}bp -> {normalized_size}bp"
                    )
                    total_normalized += 1

            # Create new record
            normalized_record = SeqRecord(
                Seq(seq),
                id=record.id,
                description=record.description
            )
            normalized_records.append(normalized_record)

        # Write normalized assembly
        normalized_assembly = self.output_dir / "assembly_normalized.fasta"
        with open(normalized_assembly, 'w') as f:
            SeqIO.write(normalized_records, f, 'fasta')

        self.logger.info(f"  Normalized {total_normalized} gaps")
        self.logger.info(f"  Wrote normalized assembly to {normalized_assembly}")

        return normalized_assembly

    def _run_validation(self, assembly: Path):
        """
        Run validation on filled assembly

        Args:
            assembly: Path to assembly FASTA
        """
        try:
            # Find remaining gaps
            remaining_gaps = self._find_gaps(assembly)

            validation_results = {
                'total_remaining_gaps': len(remaining_gaps),
                'gaps': []
            }

            for gap in remaining_gaps:
                validation_results['gaps'].append({
                    'name': gap['name'],
                    'chrom': gap['chrom'],
                    'start': gap['start'],
                    'end': gap['end'],
                    'size': gap['size']
                })

            # Save validation results
            validation_file = self.output_dir / "validation_results.json"
            with open(validation_file, 'w') as f:
                json.dump(validation_results, f, indent=2)

            self.logger.info(f"  Validation complete: {len(remaining_gaps)} gaps remaining")

        except Exception as e:
            self.logger.warning(f"Validation error: {e}")

    def _save_final_stats(self):
        """Save final statistics to file"""
        stats_file = self.output_dir / "final_stats.json"

        output = {
            'timestamp': datetime.now().isoformat(),
            'initial_assembly': str(self.initial_assembly),
            'statistics': self.stats
        }

        with open(stats_file, 'w') as f:
            json.dump(output, f, indent=2)

        # Also write summary text file
        summary_file = self.output_dir / "summary.txt"
        with open(summary_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("GAP FILLING SUMMARY\n")
            f.write("=" * 60 + "\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Initial assembly: {self.initial_assembly}\n")
            f.write(f"\n")
            f.write(f"Iterations: {self.stats['iterations']}\n")
            f.write(f"Initial gaps: {self.stats['total_gaps_initial']}\n")
            f.write(f"Completely filled: {self.stats['gaps_completely_filled']}\n")
            f.write(f"Partially filled: {self.stats['gaps_partially_filled']}\n")
            f.write(f"Failed: {self.stats['gaps_failed']}\n")
            f.write(f"Total bp filled: {self.stats['total_bp_filled']:,}\n")
            f.write("=" * 60 + "\n")


def main():
    """Command line interface"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Iterative Gap Filler V2 - Multi-round gap filling"
    )
    parser.add_argument(
        "--assembly", "-a",
        required=True,
        help="Input assembly FASTA file"
    )
    parser.add_argument(
        "--hifi-reads",
        help="HiFi reads file (FASTQ/FASTA)"
    )
    parser.add_argument(
        "--ont-reads",
        help="ONT reads file (FASTQ/FASTA)"
    )
    parser.add_argument(
        "--output", "-o",
        default="output",
        help="Output directory (default: output)"
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=8,
        help="Number of threads (default: 8)"
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=10,
        help="Maximum number of iterations (default: 10)"
    )
    parser.add_argument(
        "--min-gap-size",
        type=int,
        default=100,
        help="Minimum gap size to process (default: 100)"
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum mapping quality (default: 20)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Run gap filling
    filler = IterativeGapFiller(
        assembly_file=args.assembly,
        hifi_reads=args.hifi_reads,
        ont_reads=args.ont_reads,
        output_dir=args.output,
        threads=args.threads,
        max_iterations=args.max_iterations,
        min_gap_size=args.min_gap_size,
        min_mapq=args.min_mapq
    )

    final_assembly = filler.run()
    print(f"\nFinal assembly: {final_assembly}")


if __name__ == "__main__":
    main()
