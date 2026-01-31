#!/usr/bin/env python3
"""
Reads Cache - Filter and cache reads for gap filling

Core optimization: Instead of re-aligning reads every iteration,
we filter out reads that are anchored in non-gap regions (definitely
not useful for gap filling) and keep reads that may be useful.

Key insight:
- Reads fully anchored in non-gap regions cannot help fill gaps
- Reads that overlap with gap regions or have soft-clips near gaps
  may be useful for filling (spanning or flanking)
- Gap regions may shrink across iterations, so we keep all potentially
  useful reads rather than trying to predict which gaps they belong to
"""

import logging
import gzip
import os
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from collections import defaultdict

import pysam


@dataclass
class GapRegion:
    """Represents a gap region in the assembly"""
    chrom: str
    start: int
    end: int
    name: str


@dataclass
class CacheStats:
    """Statistics about reads filtering"""
    total_reads: int = 0
    anchored_reads: int = 0  # Far from all gaps (filtered out)
    kept_reads: int = 0  # Near at least one gap (kept)


class ReadsCache:
    """
    Manages filtered reads for gap filling.

    Strategy (conservative - reverse filtering):
    - ONLY filter out reads that have high-quality alignment FAR from all gaps
    - Keep ALL other reads

    This is more conservative than positive filtering (looking for useful reads)
    and avoids accidentally discarding potentially useful reads.
    """

    # Distance threshold: reads within this distance of a gap are kept
    GAP_PROXIMITY = 1000  # bp

    def __init__(self, output_dir: Path, threads: int = 8):
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.logger = logging.getLogger(__name__)

        # Create cache directory
        self.cache_dir = self.output_dir / "reads_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Gap regions (set during filtering)
        self.gap_regions: List[GapRegion] = []
        self.gap_intervals: Dict[str, List[Tuple[int, int, str]]] = {}  # chrom -> [(start, end, name)]

        # Cached file paths
        self.filtered_hifi_path: Optional[Path] = None
        self.filtered_ont_path: Optional[Path] = None

        # Statistics
        self.hifi_stats = CacheStats()
        self.ont_stats = CacheStats()

    def set_gap_regions(self, gaps: List[Dict]):
        """
        Set gap regions for filtering.

        Args:
            gaps: List of gap dicts with 'chrom', 'start', 'end', 'name'
        """
        self.gap_regions = []
        self.gap_intervals = defaultdict(list)

        for gap in gaps:
            region = GapRegion(
                chrom=gap['chrom'],
                start=gap['start'],
                end=gap['end'],
                name=gap.get('name', f"{gap['chrom']}_{gap['start']}_{gap['end']}")
            )
            self.gap_regions.append(region)
            self.gap_intervals[gap['chrom']].append((gap['start'], gap['end'], region.name))

        # Sort intervals for efficient lookup
        for chrom in self.gap_intervals:
            self.gap_intervals[chrom].sort(key=lambda x: x[0])

        self.logger.info(f"Set {len(self.gap_regions)} gap regions for filtering")

    def filter_bam(self, bam_file: Path, read_type: str, min_mapq: int = 20) -> Path:
        """
        Filter BAM file to keep only potentially useful reads.

        Args:
            bam_file: Input BAM file
            read_type: 'hifi' or 'ont'
            min_mapq: Minimum mapping quality

        Returns:
            Path to filtered FASTA file
        """
        if not self.gap_regions:
            raise ValueError("Gap regions not set. Call set_gap_regions() first.")

        output_fasta = self.cache_dir / f"filtered_{read_type}.fasta.gz"
        stats = CacheStats()

        self.logger.info(f"Filtering {read_type} reads from {bam_file}...")

        # Track seen read names to avoid duplicates
        seen_reads: Set[str] = set()
        kept_reads: Dict[str, str] = {}  # name -> sequence

        try:
            with pysam.AlignmentFile(str(bam_file), 'rb') as bam:
                for read in bam.fetch():
                    stats.total_reads += 1

                    # Skip unmapped, secondary, supplementary
                    if read.is_unmapped:
                        continue

                    # Skip if already seen (keep primary alignment info)
                    if read.query_name in seen_reads:
                        continue
                    seen_reads.add(read.query_name)

                    # Check if read is anchored far from gaps
                    # Only filter out reads that are: high-mapq + no-softclip + far from gaps
                    keep, reason = self._is_useful_read(read, min_mapq)

                    if keep:
                        stats.kept_reads += 1
                        # Store read sequence
                        seq = read.query_sequence
                        if seq and len(seq) >= 500:
                            kept_reads[read.query_name] = seq
                    else:
                        stats.anchored_reads += 1

                    # Progress logging
                    if stats.total_reads % 1000000 == 0:
                        self.logger.info(f"  Processed {stats.total_reads:,} reads, "
                                       f"kept {stats.kept_reads:,} ({100*stats.kept_reads/stats.total_reads:.1f}%)")

        except Exception as e:
            self.logger.error(f"Error filtering BAM: {e}")
            raise

        # Write filtered reads to FASTA
        self.logger.info(f"Writing {len(kept_reads):,} filtered reads to {output_fasta}...")

        with gzip.open(output_fasta, 'wt') as f:
            for name, seq in kept_reads.items():
                f.write(f">{name}\n{seq}\n")

        # Update stats
        if read_type == 'hifi':
            self.hifi_stats = stats
            self.filtered_hifi_path = output_fasta
        else:
            self.ont_stats = stats
            self.filtered_ont_path = output_fasta

        self.logger.info(f"Filtering complete: {stats.total_reads:,} total → "
                        f"{stats.kept_reads:,} kept ({100*stats.kept_reads/max(1,stats.total_reads):.1f}%)")
        self.logger.info(f"  Kept (near gap): {stats.kept_reads:,}, "
                        f"Filtered (far from gaps): {stats.anchored_reads:,}")

        return output_fasta

    def _is_useful_read(self, read: pysam.AlignedSegment, min_mapq: int = 20) -> Tuple[bool, str]:
        """
        Determine if a read should be kept for gap filling.

        Strategy: Only filter out reads that are ANCHORED far from gaps.

        "Anchored" means:
        - High mapping quality (confident alignment)
        - No significant soft-clips (fully aligned)

        If a read is anchored AND far from all gaps → filter out
        Otherwise → keep

        Returns:
            (keep, reason) where reason describes why kept/filtered
        """
        chrom = read.reference_name
        read_start = read.reference_start
        read_end = read.reference_end

        # Keep reads with missing coordinates
        if read_start is None or read_end is None:
            return True, 'no_coords'

        # Filter out reads on chromosomes without gaps
        if chrom not in self.gap_intervals:
            return False, 'no_gaps'

        # Check if read is "anchored" (high-quality, full alignment)
        is_anchored = self._is_read_anchored(read, min_mapq)

        # If not anchored, keep it (might be useful)
        if not is_anchored:
            return True, 'not_anchored'

        # Read is anchored. Check if it's far from all gaps.
        for gap_start, gap_end, gap_name in self.gap_intervals[chrom]:
            # Calculate distance from read to this gap
            if read_end <= gap_start:
                distance = gap_start - read_end
            elif read_start >= gap_end:
                distance = read_start - gap_end
            else:
                # Read overlaps with gap - definitely keep
                return True, 'overlaps_gap'

            # If near any gap, keep it
            if distance <= self.GAP_PROXIMITY:
                return True, 'near_gap'

        # Anchored AND far from all gaps → filter out
        return False, 'anchored_far'

    def _is_read_anchored(self, read: pysam.AlignedSegment, min_mapq: int = 20) -> bool:
        """
        Check if a read is "anchored" - has confident, full alignment.

        A read is NOT anchored if:
        - Low mapping quality (alignment is uncertain)
        - Has significant soft-clips (part of read didn't align)

        Returns:
            True if read is anchored (confident full alignment)
        """
        # Low mapq = uncertain alignment = not anchored
        if read.mapping_quality < min_mapq:
            return False

        # Check for significant soft-clips
        cigar = read.cigartuples
        if not cigar:
            return False

        # Soft-clip (4) or hard-clip (5) at either end
        CLIP_THRESHOLD = 100  # bp

        # Left clip
        if cigar[0][0] in (4, 5) and cigar[0][1] >= CLIP_THRESHOLD:
            return False

        # Right clip
        if cigar[-1][0] in (4, 5) and cigar[-1][1] >= CLIP_THRESHOLD:
            return False

        # High quality, no significant clips = anchored
        return True

    def get_filtered_reads_path(self, read_type: str) -> Optional[Path]:
        """Get path to filtered reads FASTA"""
        if read_type == 'hifi':
            return self.filtered_hifi_path
        else:
            return self.filtered_ont_path

    def get_stats(self, read_type: str) -> CacheStats:
        """Get filtering statistics"""
        if read_type == 'hifi':
            return self.hifi_stats
        else:
            return self.ont_stats

    def align_filtered_reads(self, assembly_file: Path, read_type: str) -> Path:
        """
        Align filtered reads to assembly.

        This uses the filtered FASTA instead of the original reads,
        significantly reducing alignment time.

        Args:
            assembly_file: Assembly FASTA file
            read_type: 'hifi' or 'ont'

        Returns:
            Path to aligned BAM file
        """
        import subprocess

        filtered_fasta = self.get_filtered_reads_path(read_type)
        if not filtered_fasta or not filtered_fasta.exists():
            raise ValueError(f"Filtered {read_type} reads not found. Run filter_bam() first.")

        output_bam = self.cache_dir / f"filtered_{read_type}.bam"

        preset = 'map-hifi' if read_type == 'hifi' else 'map-ont'

        self.logger.info(f"Aligning filtered {read_type} reads...")

        cmd = (
            f"minimap2 -ax {preset} -t {self.threads} {assembly_file} {filtered_fasta} | "
            f"samtools sort -@ {self.threads} -m 2G -o {output_bam} - && "
            f"samtools index {output_bam}"
        )

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            self.logger.info(f"  Created {output_bam}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Alignment failed: {e.stderr}")
            raise

        return output_bam

    def is_cached(self, read_type: str) -> bool:
        """Check if filtered reads exist"""
        path = self.get_filtered_reads_path(read_type)
        return path is not None and path.exists()

    def get_summary(self) -> Dict:
        """Get summary of cache status"""
        return {
            'gap_count': len(self.gap_regions),
            'hifi': {
                'cached': self.is_cached('hifi'),
                'path': str(self.filtered_hifi_path) if self.filtered_hifi_path else None,
                'stats': {
                    'total': self.hifi_stats.total_reads,
                    'kept': self.hifi_stats.kept_reads,
                    'filtered': self.hifi_stats.anchored_reads,
                    'kept_ratio': self.hifi_stats.kept_reads / max(1, self.hifi_stats.total_reads)
                }
            },
            'ont': {
                'cached': self.is_cached('ont'),
                'path': str(self.filtered_ont_path) if self.filtered_ont_path else None,
                'stats': {
                    'total': self.ont_stats.total_reads,
                    'kept': self.ont_stats.kept_reads,
                    'filtered': self.ont_stats.anchored_reads,
                    'kept_ratio': self.ont_stats.kept_reads / max(1, self.ont_stats.total_reads)
                }
            }
        }


class GapReadsExtractor:
    """
    Extract reads for specific gaps from filtered BAM.

    This is used during gap filling to get reads for a specific gap
    from the pre-filtered BAM file.
    """

    def __init__(self, threads: int = 8, min_mapq: int = 20):
        self.threads = threads
        self.min_mapq = min_mapq
        self.logger = logging.getLogger(__name__)

        # BAM handle cache
        self._bam_handles: Dict[str, pysam.AlignmentFile] = {}

    def extract_gap_reads(self,
                          bam_file: Path,
                          chrom: str,
                          gap_start: int,
                          gap_end: int,
                          read_type: str = 'hifi') -> Dict:
        """
        Extract reads for a specific gap.

        Returns dict with:
        - spanning: reads that span the gap
        - left: reads flanking left side
        - right: reads flanking right side
        """
        result = {
            'spanning': [],
            'left': [],
            'right': [],
            'spanning_count': 0,
            'left_count': 0,
            'right_count': 0
        }

        bam = self._get_bam(str(bam_file))
        if bam is None:
            return result

        seen_names = set()
        window = 500

        try:
            # Fetch reads around gap region
            for read in bam.fetch(chrom, max(0, gap_start - window), gap_end + window):
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.mapping_quality < self.min_mapq:
                    continue
                if read.query_name in seen_names:
                    continue

                seq = read.query_sequence
                if not seq or len(seq) < 500:
                    continue

                seen_names.add(read.query_name)
                read_start = read.reference_start
                read_end = read.reference_end

                if read_start is None or read_end is None:
                    continue

                # Classify read
                if read_start <= gap_start and read_end >= gap_end:
                    # Spanning read
                    result['spanning'].append((seq, read.query_name, read_type))
                    result['spanning_count'] += 1
                elif read_end and read_end <= gap_start + window and read_end >= gap_start - window:
                    # Left flanking
                    cigar = read.cigartuples
                    if cigar and cigar[-1][0] in (4, 5):
                        result['left'].append((seq, read.query_name, read_type))
                        result['left_count'] += 1
                elif read_start and read_start >= gap_end - window and read_start <= gap_end + window:
                    # Right flanking
                    cigar = read.cigartuples
                    if cigar and cigar[0][0] in (4, 5):
                        result['right'].append((seq, read.query_name, read_type))
                        result['right_count'] += 1

        except Exception as e:
            self.logger.warning(f"Error extracting reads: {e}")

        return result

    def _get_bam(self, bam_file: str) -> Optional[pysam.AlignmentFile]:
        """Get or open BAM handle"""
        if bam_file not in self._bam_handles:
            try:
                self._bam_handles[bam_file] = pysam.AlignmentFile(bam_file, 'rb')
            except Exception as e:
                self.logger.error(f"Cannot open BAM {bam_file}: {e}")
                return None
        return self._bam_handles[bam_file]

    def close(self):
        """Close all BAM handles"""
        for handle in self._bam_handles.values():
            try:
                handle.close()
            except:
                pass
        self._bam_handles.clear()
