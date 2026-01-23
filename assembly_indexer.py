#!/usr/bin/env python3
"""
AssemblyIndexer - Unified assembly sequence accessor with optimization
Provides fast indexed access to assembly sequences, reducing memory usage
"""

import logging
from pathlib import Path
from Bio import SeqIO

try:
    from pyfaidx import Fasta
    HAVE_PYFAIDX = True
except ImportError:
    HAVE_PYFAIDX = False


class AssemblyIndexer:
    """
    Unified assembly indexer with pyfaidx optimization

    Features:
    - Fast indexed access using pyfaidx
    - Automatic fallback to BioPython SeqIO
    - Memory-efficient sequence access
    - Global caching with automatic cleanup

    Usage:
        indexer = AssemblyIndexer("assembly.fasta")
        seq = indexer.get_sequence("chr1", 1000, 2000)
        length = indexer.get_length("chr1")
    """

    # Class-level cache for shared indexers across modules
    _global_cache = {}
    _max_cache_size = 3  # Keep max 3 different assemblies indexed

    def __init__(self, assembly_file, use_cache=True):
        """
        Initialize assembly indexer

        Args:
            assembly_file: Path to assembly FASTA file
            use_cache: Use global cache to share indices (default: True)
        """
        self.assembly_file = Path(assembly_file)
        self.use_cache = use_cache
        self.logger = logging.getLogger(__name__)

        # Validate file exists
        if not self.assembly_file.exists():
            raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

        # Get canonical path for caching
        self.canonical_path = str(self.assembly_file.resolve())

        # Initialize index
        self.fasta_index = None
        self._seqio_cache = {}  # Fallback cache for SeqIO

        # Try to use pyfaidx for fast access
        if HAVE_PYFAIDX:
            self._init_pyfaidx()
        else:
            self.logger.warning(
                "pyfaidx not available. Install with 'pip install pyfaidx' for "
                "10-1000x faster sequence access and lower memory usage."
            )
            self.logger.info("Falling back to BioPython SeqIO (slower)")

    def _init_pyfaidx(self):
        """Initialize pyfaidx index"""
        try:
            # Check if already in global cache
            if self.use_cache and self.canonical_path in self._global_cache:
                self.fasta_index = self._global_cache[self.canonical_path]
                self.logger.debug(f"Using cached index for {self.assembly_file.name}")
                return

            # Create new index
            self.logger.debug(f"Creating pyfaidx index for {self.assembly_file.name}")
            self.fasta_index = Fasta(str(self.assembly_file))

            # Add to global cache
            if self.use_cache:
                # Clean up cache if too large
                if len(self._global_cache) >= self._max_cache_size:
                    # Remove oldest entry
                    oldest_key = next(iter(self._global_cache))
                    old_index = self._global_cache.pop(oldest_key)
                    try:
                        old_index.close()
                    except:
                        pass
                    self.logger.debug(f"Removed old index from cache: {oldest_key}")

                self._global_cache[self.canonical_path] = self.fasta_index

        except Exception as e:
            self.logger.warning(f"Failed to create pyfaidx index: {e}")
            self.logger.info("Falling back to BioPython SeqIO")
            self.fasta_index = None

    def get_sequence(self, chrom, start, end):
        """
        Get sequence from assembly (optimized)

        Args:
            chrom: Chromosome/contig name
            start: Start position (0-based)
            end: End position (exclusive)

        Returns:
            str: Sequence string
        """
        # Try pyfaidx first (fast)
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return str(self.fasta_index[chrom][start:end])
            except Exception as e:
                self.logger.debug(f"pyfaidx access failed: {e}, falling back to SeqIO")

        # Fallback to SeqIO (slower but reliable)
        return self._get_sequence_seqio(chrom, start, end)

    def _get_sequence_seqio(self, chrom, start, end):
        """Fallback sequence access using SeqIO"""
        # Check cache first
        if chrom in self._seqio_cache:
            seq = self._seqio_cache[chrom]
            return seq[start:end]

        # Load from file
        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            if record.id == chrom:
                seq_str = str(record.seq)
                # Cache for future use (but don't cache if very large)
                if len(seq_str) < 100_000_000:  # 100MB limit
                    self._seqio_cache[chrom] = seq_str
                return seq_str[start:end]

        raise ValueError(f"Chromosome {chrom} not found in assembly")

    def get_full_sequence(self, chrom):
        """
        Get full chromosome/contig sequence

        Args:
            chrom: Chromosome/contig name

        Returns:
            str: Full sequence
        """
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return str(self.fasta_index[chrom][:])
            except Exception as e:
                self.logger.debug(f"pyfaidx access failed: {e}")

        # Fallback
        return self._get_sequence_seqio(chrom, 0, None)

    def _get_sequence_seqio(self, chrom, start, end):
        """Fallback to get sequence (full or partial)"""
        # Check cache
        if chrom in self._seqio_cache:
            seq = self._seqio_cache[chrom]
            if end is None:
                return seq
            return seq[start:end]

        # Load from file
        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            if record.id == chrom:
                seq_str = str(record.seq)
                # Cache if not too large
                if len(seq_str) < 100_000_000:
                    self._seqio_cache[chrom] = seq_str

                if end is None:
                    return seq_str
                return seq_str[start:end]

        raise ValueError(f"Chromosome {chrom} not found")

    def get_length(self, chrom):
        """
        Get chromosome/contig length

        Args:
            chrom: Chromosome/contig name

        Returns:
            int: Length in bp
        """
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return len(self.fasta_index[chrom])
            except:
                pass

        # Fallback
        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            if record.id == chrom:
                return len(record.seq)

        raise ValueError(f"Chromosome {chrom} not found")

    def get_all_chroms(self):
        """
        Get list of all chromosome/contig names

        Returns:
            list: Chromosome names
        """
        if self.fasta_index:
            try:
                return list(self.fasta_index.keys())
            except:
                pass

        # Fallback
        chroms = []
        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            chroms.append(record.id)
        return chroms

    def count_n_bases(self, chrom=None):
        """
        Count N bases in assembly or specific chromosome

        Args:
            chrom: Chromosome name (None for all)

        Returns:
            int: Number of N bases
        """
        if chrom:
            seq = self.get_full_sequence(chrom)
            return seq.upper().count('N')

        # Count in all chromosomes
        total = 0
        if self.fasta_index:
            try:
                for chrom_name in self.fasta_index.keys():
                    seq = str(self.fasta_index[chrom_name][:])
                    total += seq.upper().count('N')
                return total
            except:
                pass

        # Fallback
        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            total += str(record.seq).upper().count('N')
        return total

    def get_stats(self):
        """
        Get assembly statistics

        Returns:
            dict: Statistics including total size, N count, etc.
        """
        chroms = self.get_all_chroms()
        total_size = 0
        total_n = 0

        for chrom in chroms:
            length = self.get_length(chrom)
            n_count = self.count_n_bases(chrom)
            total_size += length
            total_n += n_count

        return {
            'num_sequences': len(chroms),
            'total_size': total_size,
            'total_n_bases': total_n,
            'n_percentage': (total_n / total_size * 100) if total_size > 0 else 0
        }

    def close(self):
        """Close the index and free resources"""
        if self.fasta_index and not self.use_cache:
            try:
                self.fasta_index.close()
            except:
                pass
            self.fasta_index = None

        # Clear local cache
        self._seqio_cache.clear()

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()
        return False

    @classmethod
    def clear_global_cache(cls):
        """Clear the global index cache"""
        for index in cls._global_cache.values():
            try:
                index.close()
            except:
                pass
        cls._global_cache.clear()

    @classmethod
    def set_cache_size(cls, size):
        """Set maximum cache size"""
        cls._max_cache_size = max(1, size)


def get_indexer(assembly_file, use_cache=True):
    """
    Convenience function to get an assembly indexer

    Args:
        assembly_file: Path to assembly FASTA
        use_cache: Use global cache (default: True)

    Returns:
        AssemblyIndexer instance
    """
    return AssemblyIndexer(assembly_file, use_cache=use_cache)
