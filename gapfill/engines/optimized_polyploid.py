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
                 use_ambiguous_reads: bool = True):

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

        self.logger = logging.getLogger(__name__)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        self.hic_analyzer: Optional[HiCAnalyzer] = None

        # Store normalized assemblies
        self.normalized_assemblies: Dict[str, Path] = {}

        self._validate_inputs()

        self.logger.info("=" * 60)
        self.logger.info("OptimizedPolyploidEngine initialized")
        self.logger.info(f"  Ploidy: {self.num_haplotypes}n")
        self.logger.info(f"  OPTIMIZATION: Batch alignment (1 align per iteration)")
        self.logger.info(f"  Expected alignments: {self.max_iterations * 2} "
                        f"(vs {self.num_haplotypes * self.max_iterations * 2} original)")
        self.logger.info("=" * 60)

    def _validate_inputs(self):
        for hap in self.haplotypes:
            if not hap.exists():
                raise FileNotFoundError(f"Haplotype not found: {hap}")
        if not self.hifi_reads and not self.ont_reads:
            raise ValueError("At least one reads file required")

    def run(self) -> Dict[str, Path]:
        """Run optimized polyploid gap filling"""
        self.logger.info("Starting optimized polyploid gap filling...")

        # =====================================================================
        # PHASE 1: One-time setup (normalization + phasing)
        # =====================================================================
        self.logger.info("\n" + "=" * 60)
        self.logger.info("PHASE 1: Initial Setup (one-time)")
        self.logger.info("=" * 60)

        # Initialize phaser
        phaser = ReadPhaser(self.haplotypes, self.threads, work_dir=self.output_dir)

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

        # Detect SNPs using alignment-based method
        self.logger.info("STEP 1: Detecting haplotype-specific SNPs (alignment-based)")
        snp_db = phaser.detect_haplotype_snps()

        # Save SNP database
        snp_file = self.output_dir / "snp_database.json"
        self._save_snp_db(snp_db, snp_file)

        # Phase reads (one-time, using normalized hap1 as reference)
        phased_hifi = {}
        if self.hifi_reads:
            self.logger.info("STEP 2a: Phasing HiFi reads...")
            hifi_bam = self.output_dir / "phase_hifi.bam"
            self._align_reads(self.hifi_reads, ref_assembly, hifi_bam, 'map-hifi')
            phased_hifi = phaser.phase_reads_from_bam(
                hifi_bam, snp_db, self.output_dir / "phased", 'hifi'
            )

        phased_ont = {}
        if self.ont_reads:
            self.logger.info("STEP 2b: Phasing ONT reads...")
            ont_bam = self.output_dir / "phase_ont.bam"
            self._align_reads(self.ont_reads, ref_assembly, ont_bam, 'map-ont')
            phased_ont = phaser.phase_reads_from_bam(
                ont_bam, snp_db, self.output_dir / "phased", 'ont'
            )

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

        # Track filled gaps
        filled_gaps = {hap: set() for hap in self.hap_names}
        failed_gaps = {hap: set() for hap in self.hap_names}

        for iteration in range(1, self.max_iterations + 1):
            self.logger.info(f"\n--- Iteration {iteration} ---")

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
                self._align_reads(merged_hifi, merged_ref, hifi_bam, 'map-hifi')

            if merged_ont and merged_ont.stat().st_size > 0:
                ont_bam = iter_dir / "merged_ont.bam"
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

                filler = GapFiller(
                    assembly_file=str(hap_assembly),
                    hifi_bam=str(hap_hifi_bam) if hap_hifi_bam else None,
                    ont_bam=str(hap_ont_bam) if hap_ont_bam else None,
                    threads=self.threads,
                    work_dir=str(work_dir)
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

                filler.close()

                # Apply fills
                if fill_results:
                    current_assemblies[hap_name] = self._apply_fills(
                        current_assemblies[hap_name], fill_results, hap_dir
                    )

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

        return final_assemblies

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
        cmd = (f"minimap2 -ax {preset} -t {self.threads} {ref} {reads} | "
               f"samtools sort -@ {self.threads} -o {output_bam} - && "
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
