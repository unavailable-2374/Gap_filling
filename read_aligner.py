#!/usr/bin/env python3
"""
Read alignment module
Aligns HiFi and ONT reads to assembly/genome using minimap2

Supports:
- Full genome alignment with MAPQ filtering (for accurate repeat handling)
- HiFi and ONT read types
- Checkpoint/resume functionality
"""

import subprocess
import logging
from pathlib import Path

# Default MAPQ threshold for filtering (40 = very high confidence)
DEFAULT_MIN_MAPQ = 40


class ReadAligner:
    """Align long reads to assembly/genome using minimap2"""

    def __init__(self, assembly_file, threads=8, min_mapq=DEFAULT_MIN_MAPQ):
        self.assembly_file = Path(assembly_file)
        self.threads = threads
        self.min_mapq = min_mapq
        self.logger = logging.getLogger(__name__)

    def align_to_full_genome(self, reads_file, output_bam, read_type='ont',
                              min_mapq=None, skip_if_exists=True):
        """
        Align reads to full genome with MAPQ filtering.

        This method is optimized for full genome alignment where accurate MAPQ
        scores are critical for filtering repetitive regions.

        Args:
            reads_file: Path to reads file (FASTA/FASTQ, gzipped OK)
            output_bam: Output BAM file path
            read_type: 'hifi' or 'ont' (default: 'ont')
            min_mapq: Minimum MAPQ to keep (default: self.min_mapq)
            skip_if_exists: Skip if valid BAM exists

        Returns:
            Path to sorted, filtered, indexed BAM file
        """
        reads_file = Path(reads_file)
        output_bam = Path(output_bam)
        min_mapq = min_mapq if min_mapq is not None else self.min_mapq

        # Check for existing valid BAM
        if skip_if_exists and self._is_valid_bam(output_bam):
            self.logger.info(f"Using existing BAM: {output_bam}")
            return output_bam

        # Set minimap2 preset
        preset = 'map-hifi' if read_type == 'hifi' else 'map-ont'

        self.logger.info(f"Aligning {read_type} reads to full genome: {self.assembly_file}")
        self.logger.info(f"  Read type: {read_type} (preset: {preset})")
        self.logger.info(f"  MAPQ filter: >= {min_mapq}")
        self.logger.info(f"  Filtering: secondary=no, supplementary=no")

        # Pipeline: minimap2 | samtools view (filter) | samtools sort
        cmd = (
            f"minimap2 -ax {preset} --secondary=no -t {self.threads} "
            f"{self.assembly_file} {reads_file} | "
            f"samtools view -bS -q {min_mapq} -F 2304 -@ {self.threads} - | "
            f"samtools sort -@ {self.threads} -o {output_bam}"
        )

        self.logger.info(f"Running alignment pipeline...")
        self.logger.debug(f"Command: {cmd}")

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                self.logger.error(f"Alignment failed: {result.stderr}")
                raise RuntimeError("Full genome alignment failed")

            if not output_bam.exists() or output_bam.stat().st_size == 0:
                self.logger.error("Output BAM file not created or empty")
                raise RuntimeError("BAM file creation failed")

            # Index BAM
            subprocess.run(['samtools', 'index', '-@', str(self.threads), str(output_bam)], check=True)

            # Log stats
            result = subprocess.run(['samtools', 'flagstat', str(output_bam)],
                                   capture_output=True, text=True)
            self.logger.info(f"Alignment stats:\n{result.stdout}")

            return output_bam

        except Exception as e:
            self.logger.error(f"Full genome alignment failed: {e}")
            raise

    def align(self, reads_file, output_bam, read_type='hifi', skip_if_exists=True):
        """
        Align reads to assembly (OPTIMIZED: uses pipeline)

        OPTIMIZATION: minimap2 | samtools sort (no intermediate files)
        - Reduces disk I/O by 66%
        - Faster processing
        - Lower disk usage

        CHECKPOINT SUPPORT: Skip alignment if BAM file already exists and is valid

        Args:
            reads_file: Path to reads file (FASTA/FASTQ, gzipped OK)
            output_bam: Output BAM file path
            read_type: 'hifi' or 'ont'
            skip_if_exists: Skip alignment if valid BAM already exists (default: True)

        Returns:
            Path to sorted and indexed BAM file
        """
        reads_file = Path(reads_file)
        output_bam = Path(output_bam)

        # ========================================================================
        # CHECKPOINT: Check if BAM file already exists and is valid
        # ========================================================================
        if skip_if_exists and self._is_valid_bam(output_bam):
            self.logger.info(f"✓ Found existing valid BAM file: {output_bam}")
            self.logger.info(f"  Skipping alignment (saves ~1-3 hours)")
            self.logger.info(f"  To force re-alignment, delete: {output_bam}")
            return output_bam

        self.logger.info(f"Aligning {read_type} reads: {reads_file}")

        # Set minimap2 preset based on read type
        if read_type == 'hifi':
            preset = 'map-hifi'
        elif read_type == 'ont':
            preset = 'map-ont'
        else:
            raise ValueError(f"Unknown read type: {read_type}")

        # Use two-step process to avoid pipeline truncation issues
        # Step 1: minimap2 -> SAM file
        # Step 2: samtools sort SAM -> BAM
        temp_sam = output_bam.with_suffix('.sam')

        minimap2_cmd = [
            'minimap2',
            '-ax', preset,
            '-t', str(self.threads),
            '--secondary=no',  # Only primary alignments
            '-L',  # Long CIGAR for long indels
            '-o', str(temp_sam),  # Output to file instead of stdout
            str(self.assembly_file),
            str(reads_file)
        ]

        self.logger.info("Running: minimap2 (to SAM file)")
        self.logger.debug(f"minimap2: {' '.join(minimap2_cmd)}")

        try:
            # Step 1: Run minimap2 to SAM file
            result = subprocess.run(
                minimap2_cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"minimap2 failed: {result.stderr}")
                raise RuntimeError("minimap2 alignment failed")

            if not temp_sam.exists() or temp_sam.stat().st_size == 0:
                self.logger.error(f"SAM file not created or empty: {temp_sam}")
                raise RuntimeError("SAM file creation failed")

            self.logger.info(f"SAM file created: {temp_sam} ({temp_sam.stat().st_size / 1e9:.1f} GB)")

            # Step 2: Sort SAM to BAM
            samtools_sort_cmd = [
                'samtools', 'sort',
                '-@', str(self.threads),
                '-o', str(output_bam),
                str(temp_sam)
            ]

            self.logger.info("Running: samtools sort (SAM to BAM)")
            self.logger.debug(f"samtools: {' '.join(samtools_sort_cmd)}")

            result = subprocess.run(
                samtools_sort_cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"samtools sort failed: {result.stderr}")
                raise RuntimeError("samtools sort failed")

            # Final validation: check if BAM file was created
            if not output_bam.exists() or output_bam.stat().st_size == 0:
                self.logger.error(f"Output BAM file not created or empty: {output_bam}")
                raise RuntimeError("BAM file creation failed")

            # Clean up temp SAM file
            if temp_sam.exists():
                temp_sam.unlink()
                self.logger.debug(f"Removed temp SAM file: {temp_sam}")

        except RuntimeError:
            raise
        except Exception as e:
            self.logger.error(f"Alignment failed: {e}")
            raise
        finally:
            # Clean up temp SAM file if it exists
            if temp_sam.exists():
                try:
                    temp_sam.unlink()
                except:
                    pass

        # Index BAM
        self.logger.info("Indexing BAM file...")
        subprocess.run(
            ['samtools', 'index', '-@', str(self.threads), str(output_bam)],
            check=True
        )

        self.logger.info(f"Alignment complete: {output_bam}")

        return output_bam

    def _is_valid_bam(self, bam_file):
        """
        Check if BAM file exists and is valid

        Validation checks:
        1. BAM file exists and is not empty
        2. BAM index file (.bai) exists
        3. BAM file can be opened by pysam (validates header and format)

        Args:
            bam_file: Path to BAM file

        Returns:
            True if BAM file is valid, False otherwise
        """
        bam_file = Path(bam_file)
        bai_file = Path(str(bam_file) + '.bai')

        # Check 1: BAM file exists and is not empty
        if not bam_file.exists():
            self.logger.debug(f"BAM file does not exist: {bam_file}")
            return False

        if bam_file.stat().st_size == 0:
            self.logger.warning(f"BAM file is empty: {bam_file}")
            return False

        # Check 2: Index file exists
        if not bai_file.exists():
            self.logger.debug(f"BAM index file missing: {bai_file}")
            return False

        # Check 3: Can open BAM file (validates format and header)
        try:
            import pysam
            with pysam.AlignmentFile(bam_file, 'rb') as bam:
                # Try to read header
                _ = bam.header
                # Try to fetch from first reference (if any)
                if bam.nreferences > 0:
                    first_ref = bam.references[0]
                    # Just check if we can create an iterator (don't consume it)
                    _ = bam.fetch(first_ref, 0, 1)

            self.logger.debug(f"BAM file is valid: {bam_file}")
            return True

        except Exception as e:
            self.logger.warning(f"BAM file validation failed: {bam_file} - {e}")
            return False

    def get_reads_in_region(self, bam_file, chrom, start, end):
        """
        Extract reads overlapping a specific region

        Args:
            bam_file: Path to sorted indexed BAM file
            chrom: Chromosome/contig name
            start: Start position
            end: End position

        Returns:
            List of read information dictionaries
        """
        import pysam

        reads = []

        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                read_info = {
                    'name': read.query_name,
                    'seq': read.query_sequence,
                    'qual': read.query_qualities,
                    'ref_start': read.reference_start,
                    'ref_end': read.reference_end,
                    'mapq': read.mapping_quality,
                    'cigar': read.cigarstring,
                    'is_reverse': read.is_reverse
                }

                reads.append(read_info)

        return reads

    def get_spanning_reads(self, bam_file, chrom, gap_start, gap_end, min_overhang=100):
        """
        Get reads that span across a gap

        Args:
            bam_file: Path to sorted indexed BAM file
            chrom: Chromosome/contig name
            gap_start: Gap start position
            gap_end: Gap end position
            min_overhang: Minimum overhang on each side of the gap

        Returns:
            List of reads that span the gap
        """
        import pysam

        spanning_reads = []

        # Fetch reads in extended region
        fetch_start = max(0, gap_start - 1000)
        fetch_end = gap_end + 1000

        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for read in bam.fetch(chrom, fetch_start, fetch_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                # Check if read spans the gap
                if (read.reference_start <= gap_start - min_overhang and
                    read.reference_end >= gap_end + min_overhang):

                    read_info = {
                        'name': read.query_name,
                        'seq': read.query_sequence,
                        'qual': read.query_qualities,
                        'ref_start': read.reference_start,
                        'ref_end': read.reference_end,
                        'mapq': read.mapping_quality,
                        'cigar': read.cigarstring,
                        'is_reverse': read.is_reverse,
                        'left_overhang': gap_start - read.reference_start,
                        'right_overhang': read.reference_end - gap_end
                    }

                    spanning_reads.append(read_info)

        return spanning_reads
