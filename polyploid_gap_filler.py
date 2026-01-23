#!/usr/bin/env python3
"""
Polyploid Gap Filler - Gap filling for polyploid genomes with phasing support

Supports diploid (2n), tetraploid (4n), hexaploid (6n), and higher ploidy levels.

Workflow:
1. Align reads to all haplotype assemblies
2. Detect haplotype-specific variants (SNPs)
3. Phase reads to haplotypes based on SNP profiles
4. Run independent gap filling for each haplotype
5. Output filled multi-haplotype assembly

Phasing methods:
- WhatsHap: External tool for sophisticated phasing
- Built-in SNP-based: Lightweight integrated solution

Author: Gap Filling Pipeline
"""

import logging
import subprocess
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from iterative_gapfiller import IterativeGapFiller
from assembly_indexer import AssemblyIndexer


class ReadPhaser:
    """
    Phase reads to haplotypes based on haplotype-specific SNPs

    Supports two methods:
    1. WhatsHap-based phasing (external tool)
    2. Built-in SNP-based phasing (integrated)
    """

    def __init__(self, haplotype_assemblies: List[Path],
                 threads: int = 8,
                 min_snp_qual: int = 20,
                 min_read_snps: int = 2,
                 work_dir: Optional[Path] = None):
        """
        Initialize ReadPhaser

        Args:
            haplotype_assemblies: List of haplotype FASTA files
            threads: Number of threads
            min_snp_qual: Minimum SNP quality score
            min_read_snps: Minimum SNPs required to assign a read
            work_dir: Working directory
        """
        self.haplotypes = haplotype_assemblies
        self.num_haplotypes = len(haplotype_assemblies)
        self.threads = threads
        self.min_snp_qual = min_snp_qual
        self.min_read_snps = min_read_snps
        self.work_dir = work_dir or Path('.')

        self.logger = logging.getLogger(__name__)

        # Haplotype names
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        self.logger.info(f"ReadPhaser initialized with {self.num_haplotypes} haplotypes")

    def detect_haplotype_snps(self, reads_file: Path,
                              method: str = 'builtin') -> Dict[str, Dict]:
        """
        Detect haplotype-specific SNPs

        Args:
            reads_file: Path to reads file
            method: 'builtin' or 'whatshap'

        Returns:
            Dictionary of SNP information per haplotype
        """
        self.logger.info("Detecting haplotype-specific SNPs...")

        if method == 'whatshap':
            return self._detect_snps_whatshap(reads_file)
        else:
            return self._detect_snps_builtin(reads_file)

    def _detect_snps_builtin(self, reads_file: Path) -> Dict[str, Dict]:
        """
        Built-in SNP detection by comparing haplotype assemblies

        Strategy:
        1. Align haplotypes to each other
        2. Find positions where haplotypes differ
        3. These are haplotype-specific markers
        """
        self.logger.info("  Using built-in SNP detection...")

        snp_db = {}  # {chrom: {pos: {hap1: base, hap2: base, ...}}}

        # Load all haplotype sequences
        hap_seqs = {}
        for i, hap_file in enumerate(self.haplotypes):
            hap_name = self.hap_names[i]
            hap_seqs[hap_name] = {}
            for record in SeqIO.parse(hap_file, 'fasta'):
                hap_seqs[hap_name][record.id] = str(record.seq).upper()

        # Use first haplotype as reference for chromosome names
        ref_hap = self.hap_names[0]

        for chrom in hap_seqs[ref_hap]:
            snp_db[chrom] = {}
            ref_seq = hap_seqs[ref_hap][chrom]

            # Compare with other haplotypes
            for other_hap in self.hap_names[1:]:
                if chrom not in hap_seqs[other_hap]:
                    continue

                other_seq = hap_seqs[other_hap][chrom]

                # Find SNPs (simple position-by-position comparison)
                # Note: This assumes haplotypes are well-aligned
                min_len = min(len(ref_seq), len(other_seq))

                for pos in range(min_len):
                    ref_base = ref_seq[pos]
                    other_base = other_seq[pos]

                    if ref_base != other_base and ref_base != 'N' and other_base != 'N':
                        if pos not in snp_db[chrom]:
                            snp_db[chrom][pos] = {ref_hap: ref_base}
                        snp_db[chrom][pos][other_hap] = other_base

        # Count SNPs
        total_snps = sum(len(positions) for positions in snp_db.values())
        self.logger.info(f"  Found {total_snps} haplotype-specific SNP positions")

        return snp_db

    def _detect_snps_whatshap(self, reads_file: Path) -> Dict[str, Dict]:
        """
        Use WhatsHap for SNP detection and phasing

        Requires: bcftools, whatshap
        """
        self.logger.info("  Using WhatsHap for SNP detection...")

        snp_db = {}

        # For each haplotype, call variants
        vcf_files = []

        for i, hap_file in enumerate(self.haplotypes):
            hap_name = self.hap_names[i]
            bam_file = self.work_dir / f"{hap_name}_aligned.bam"
            vcf_file = self.work_dir / f"{hap_name}_variants.vcf.gz"

            # Align reads to this haplotype
            self._align_reads(reads_file, hap_file, bam_file)

            # Call variants
            self._call_variants(hap_file, bam_file, vcf_file)

            vcf_files.append(vcf_file)

        # Merge VCFs and identify haplotype-specific variants
        snp_db = self._merge_vcfs_to_snp_db(vcf_files)

        return snp_db

    def _align_reads(self, reads_file: Path, ref_file: Path,
                     output_bam: Path) -> bool:
        """Align reads to reference"""
        try:
            # Determine read type from filename
            preset = 'map-hifi' if 'hifi' in str(reads_file).lower() else 'map-ont'

            cmd = f"minimap2 -ax {preset} -t {self.threads} {ref_file} {reads_file} | " \
                  f"samtools sort -@ {self.threads} -o {output_bam} - && " \
                  f"samtools index {output_bam}"

            subprocess.run(cmd, shell=True, check=True,
                         capture_output=True, text=True)
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed: {e}")
            return False

    def _call_variants(self, ref_file: Path, bam_file: Path,
                       vcf_file: Path) -> bool:
        """Call variants using bcftools"""
        try:
            cmd = f"bcftools mpileup -Ou -f {ref_file} {bam_file} | " \
                  f"bcftools call -mv -Oz -o {vcf_file} && " \
                  f"bcftools index {vcf_file}"

            subprocess.run(cmd, shell=True, check=True,
                         capture_output=True, text=True)
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Variant calling failed: {e}")
            return False

    def _merge_vcfs_to_snp_db(self, vcf_files: List[Path]) -> Dict[str, Dict]:
        """Merge multiple VCFs into SNP database"""
        snp_db = {}

        for i, vcf_file in enumerate(vcf_files):
            hap_name = self.hap_names[i]

            if not vcf_file.exists():
                continue

            try:
                vcf = pysam.VariantFile(str(vcf_file))

                for record in vcf:
                    chrom = record.chrom
                    pos = record.pos - 1  # Convert to 0-based

                    if chrom not in snp_db:
                        snp_db[chrom] = {}

                    if record.qual and record.qual >= self.min_snp_qual:
                        ref = record.ref
                        alt = record.alts[0] if record.alts else ref

                        if pos not in snp_db[chrom]:
                            snp_db[chrom][pos] = {}

                        snp_db[chrom][pos][hap_name] = alt

            except Exception as e:
                self.logger.warning(f"Error reading VCF {vcf_file}: {e}")

        return snp_db

    def phase_reads(self, bam_file: Path, snp_db: Dict[str, Dict],
                    output_prefix: Path) -> Dict[str, Path]:
        """
        Phase reads based on SNP profiles

        Args:
            bam_file: Aligned reads BAM
            snp_db: SNP database from detect_haplotype_snps()
            output_prefix: Output prefix for phased read files

        Returns:
            Dictionary mapping haplotype name to phased reads file
        """
        self.logger.info("Phasing reads to haplotypes...")

        # Initialize output files
        phased_reads = {hap: [] for hap in self.hap_names}
        ambiguous_reads = []

        # Track statistics
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

                    # Get haplotype assignment
                    hap_assignment = self._assign_read_to_haplotype(
                        read, snp_db
                    )

                    if hap_assignment:
                        phased_reads[hap_assignment].append(read)
                        stats['phased'][hap_assignment] += 1
                    else:
                        ambiguous_reads.append(read)
                        stats['ambiguous'] += 1

        except Exception as e:
            self.logger.error(f"Error phasing reads: {e}")
            raise

        # Write phased reads to files
        output_files = {}

        for hap_name in self.hap_names:
            output_file = Path(f"{output_prefix}_{hap_name}_reads.fasta")

            with open(output_file, 'w') as f:
                for read in phased_reads[hap_name]:
                    seq = read.query_sequence
                    if seq:
                        f.write(f">{read.query_name}\n{seq}\n")

            output_files[hap_name] = output_file
            self.logger.info(f"  {hap_name}: {stats['phased'][hap_name]} reads")

        # Write ambiguous reads (can be used by all haplotypes)
        ambiguous_file = Path(f"{output_prefix}_ambiguous_reads.fasta")
        with open(ambiguous_file, 'w') as f:
            for read in ambiguous_reads:
                seq = read.query_sequence
                if seq:
                    f.write(f">{read.query_name}\n{seq}\n")

        self.logger.info(f"  Ambiguous: {stats['ambiguous']} reads")
        self.logger.info(f"  Total processed: {stats['total']} reads")

        return output_files

    def _assign_read_to_haplotype(self, read: pysam.AlignedSegment,
                                   snp_db: Dict[str, Dict]) -> Optional[str]:
        """
        Assign a read to a haplotype based on SNP profile

        Returns haplotype name or None if ambiguous
        """
        chrom = read.reference_name

        if chrom not in snp_db:
            return None

        # Count matches to each haplotype
        hap_scores = {hap: 0 for hap in self.hap_names}
        snps_checked = 0

        # Get aligned pairs (query_pos, ref_pos)
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

        # Need minimum SNPs to make assignment
        if snps_checked < self.min_read_snps:
            return None

        # Find best haplotype
        max_score = max(hap_scores.values())
        if max_score == 0:
            return None

        # Check if there's a clear winner
        best_haps = [h for h, s in hap_scores.items() if s == max_score]

        if len(best_haps) == 1:
            # Clear assignment
            return best_haps[0]
        else:
            # Ambiguous - multiple haplotypes tie
            return None


class PolyploidGapFiller:
    """
    Gap filler for polyploid genomes

    Coordinates phasing and independent gap filling for each haplotype
    """

    def __init__(self,
                 haplotype_assemblies: List[str],
                 hifi_reads: Optional[str] = None,
                 ont_reads: Optional[str] = None,
                 output_dir: str = "polyploid_output",
                 threads: int = 8,
                 max_iterations: int = 10,
                 phasing_method: str = 'builtin',
                 use_ambiguous_reads: bool = True):
        """
        Initialize PolyploidGapFiller

        Args:
            haplotype_assemblies: List of haplotype FASTA files
            hifi_reads: Path to HiFi reads
            ont_reads: Path to ONT reads
            output_dir: Output directory
            threads: Number of threads
            max_iterations: Max gap filling iterations per haplotype
            phasing_method: 'builtin' or 'whatshap'
            use_ambiguous_reads: Whether to use ambiguous reads for all haplotypes
        """
        self.haplotypes = [Path(h) for h in haplotype_assemblies]
        self.num_haplotypes = len(self.haplotypes)
        self.hifi_reads = Path(hifi_reads) if hifi_reads else None
        self.ont_reads = Path(ont_reads) if ont_reads else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.max_iterations = max_iterations
        self.phasing_method = phasing_method
        self.use_ambiguous_reads = use_ambiguous_reads

        self.logger = logging.getLogger(__name__)

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Haplotype names
        self.hap_names = [f"hap{i+1}" for i in range(self.num_haplotypes)]

        # Validate inputs
        self._validate_inputs()

        self.logger.info("=" * 60)
        self.logger.info("PolyploidGapFiller initialized")
        self.logger.info(f"  Ploidy: {self.num_haplotypes}n")
        self.logger.info(f"  Haplotypes: {[h.name for h in self.haplotypes]}")
        self.logger.info(f"  Phasing method: {self.phasing_method}")
        self.logger.info(f"  Output: {self.output_dir}")
        self.logger.info("=" * 60)

    def _validate_inputs(self):
        """Validate input files exist"""
        for hap in self.haplotypes:
            if not hap.exists():
                raise FileNotFoundError(f"Haplotype file not found: {hap}")

        if not self.hifi_reads and not self.ont_reads:
            raise ValueError("At least one reads file (HiFi or ONT) required")

    def run(self) -> Dict[str, Path]:
        """
        Run polyploid gap filling

        Returns:
            Dictionary mapping haplotype names to filled assembly paths
        """
        self.logger.info("Starting polyploid gap filling...")

        # =================================================================
        # STEP 1: Align reads to reference haplotype for phasing
        # =================================================================
        self.logger.info("=" * 60)
        self.logger.info("STEP 1: Aligning reads for phasing")
        self.logger.info("=" * 60)

        # Use first haplotype as reference for alignment
        ref_hap = self.haplotypes[0]
        reads_file = self.hifi_reads or self.ont_reads

        aligned_bam = self.output_dir / "reads_aligned.bam"
        self._align_reads_to_ref(reads_file, ref_hap, aligned_bam)

        # =================================================================
        # STEP 2: Detect haplotype-specific SNPs
        # =================================================================
        self.logger.info("=" * 60)
        self.logger.info("STEP 2: Detecting haplotype-specific SNPs")
        self.logger.info("=" * 60)

        phaser = ReadPhaser(
            self.haplotypes,
            threads=self.threads,
            work_dir=self.output_dir
        )

        snp_db = phaser.detect_haplotype_snps(reads_file, method=self.phasing_method)

        # Save SNP database
        snp_file = self.output_dir / "snp_database.json"
        self._save_snp_db(snp_db, snp_file)

        # =================================================================
        # STEP 3: Phase reads to haplotypes
        # =================================================================
        self.logger.info("=" * 60)
        self.logger.info("STEP 3: Phasing reads to haplotypes")
        self.logger.info("=" * 60)

        phased_reads = phaser.phase_reads(
            aligned_bam,
            snp_db,
            self.output_dir / "phased"
        )

        # =================================================================
        # STEP 4: Run gap filling for each haplotype
        # =================================================================
        self.logger.info("=" * 60)
        self.logger.info("STEP 4: Running gap filling for each haplotype")
        self.logger.info("=" * 60)

        filled_assemblies = {}

        for i, (hap_name, hap_file) in enumerate(zip(self.hap_names, self.haplotypes)):
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Processing {hap_name}: {hap_file.name}")
            self.logger.info(f"{'='*60}")

            # Get phased reads for this haplotype
            hap_reads = phased_reads.get(hap_name)

            # Optionally combine with ambiguous reads
            if self.use_ambiguous_reads:
                ambiguous_file = self.output_dir / "phased_ambiguous_reads.fasta"
                if ambiguous_file.exists():
                    combined_reads = self.output_dir / f"{hap_name}_combined_reads.fasta"
                    self._combine_fasta_files([hap_reads, ambiguous_file], combined_reads)
                    hap_reads = combined_reads

            # Create haplotype-specific output directory
            hap_output = self.output_dir / hap_name

            # Run gap filling
            filled_assembly = self._run_gap_filling_for_haplotype(
                hap_file, hap_reads, hap_output
            )

            filled_assemblies[hap_name] = filled_assembly

        # =================================================================
        # STEP 5: Generate summary report
        # =================================================================
        self.logger.info("=" * 60)
        self.logger.info("STEP 5: Generating summary report")
        self.logger.info("=" * 60)

        self._generate_summary(filled_assemblies)

        return filled_assemblies

    def _align_reads_to_ref(self, reads_file: Path, ref_file: Path,
                            output_bam: Path) -> bool:
        """Align reads to reference for phasing"""
        self.logger.info(f"  Aligning {reads_file.name} to {ref_file.name}...")

        try:
            preset = 'map-hifi' if self.hifi_reads else 'map-ont'

            minimap2_cmd = [
                'minimap2', '-ax', preset,
                '-t', str(self.threads),
                '--secondary=no',
                str(ref_file), str(reads_file)
            ]

            samtools_sort = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(output_bam), '-'
            ]

            # Run pipeline
            p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            p2 = subprocess.Popen(samtools_sort, stdin=p1.stdout,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p1.stdout.close()
            _, stderr = p2.communicate()

            if p2.returncode != 0:
                self.logger.error(f"Alignment failed: {stderr.decode()[:500]}")
                return False

            # Index
            subprocess.run(['samtools', 'index', str(output_bam)], check=True)

            self.logger.info(f"  Created {output_bam}")
            return True

        except Exception as e:
            self.logger.error(f"Alignment error: {e}")
            return False

    def _save_snp_db(self, snp_db: Dict, output_file: Path):
        """Save SNP database to JSON"""
        # Convert to serializable format
        serializable = {}
        for chrom, positions in snp_db.items():
            serializable[chrom] = {str(pos): bases for pos, bases in positions.items()}

        with open(output_file, 'w') as f:
            json.dump(serializable, f, indent=2)

        self.logger.info(f"  Saved SNP database to {output_file}")

    def _combine_fasta_files(self, input_files: List[Path], output_file: Path):
        """Combine multiple FASTA files"""
        with open(output_file, 'w') as out:
            for input_file in input_files:
                if input_file and input_file.exists():
                    with open(input_file) as inp:
                        out.write(inp.read())

    def _run_gap_filling_for_haplotype(self, assembly: Path,
                                        reads: Path,
                                        output_dir: Path) -> Path:
        """Run gap filling for a single haplotype"""

        # Determine read type
        hifi_reads = None
        ont_reads = None

        if self.hifi_reads:
            hifi_reads = str(reads)
        else:
            ont_reads = str(reads)

        # Create and run gap filler
        filler = IterativeGapFiller(
            assembly_file=str(assembly),
            hifi_reads=hifi_reads,
            ont_reads=ont_reads,
            output_dir=str(output_dir),
            threads=self.threads,
            max_iterations=self.max_iterations
        )

        return filler.run()

    def _generate_summary(self, filled_assemblies: Dict[str, Path]):
        """Generate summary report"""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'num_haplotypes': self.num_haplotypes,
            'phasing_method': self.phasing_method,
            'haplotypes': {}
        }

        for hap_name, assembly_path in filled_assemblies.items():
            # Read stats if available
            stats_file = assembly_path.parent / "final_stats.json"
            if stats_file.exists():
                with open(stats_file) as f:
                    stats = json.load(f)
                summary['haplotypes'][hap_name] = stats
            else:
                summary['haplotypes'][hap_name] = {
                    'assembly': str(assembly_path)
                }

        # Save summary
        summary_file = self.output_dir / "polyploid_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        # Print summary
        self.logger.info("\n" + "=" * 60)
        self.logger.info("POLYPLOID GAP FILLING COMPLETE")
        self.logger.info("=" * 60)

        for hap_name, assembly_path in filled_assemblies.items():
            self.logger.info(f"  {hap_name}: {assembly_path}")

        self.logger.info(f"\nSummary: {summary_file}")
        self.logger.info("=" * 60)


def main():
    """Command line interface"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Polyploid Gap Filler - Gap filling for polyploid genomes"
    )
    parser.add_argument(
        "--haplotypes", "-H",
        nargs='+',
        required=True,
        help="Haplotype assembly FASTA files (e.g., -H hap1.fa hap2.fa)"
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
        default="polyploid_output",
        help="Output directory (default: polyploid_output)"
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
        help="Maximum iterations per haplotype (default: 10)"
    )
    parser.add_argument(
        "--phasing-method",
        choices=['builtin', 'whatshap'],
        default='builtin',
        help="Phasing method (default: builtin)"
    )
    parser.add_argument(
        "--no-ambiguous-reads",
        action="store_true",
        help="Do not use ambiguous reads for gap filling"
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

    # Run
    filler = PolyploidGapFiller(
        haplotype_assemblies=args.haplotypes,
        hifi_reads=args.hifi_reads,
        ont_reads=args.ont_reads,
        output_dir=args.output,
        threads=args.threads,
        max_iterations=args.max_iterations,
        phasing_method=args.phasing_method,
        use_ambiguous_reads=not args.no_ambiguous_reads
    )

    results = filler.run()

    print(f"\nFilled assemblies:")
    for hap_name, path in results.items():
        print(f"  {hap_name}: {path}")


if __name__ == "__main__":
    main()
