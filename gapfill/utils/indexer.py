#!/usr/bin/env python3
"""
Assembly Indexer - Fast indexed access to assembly sequences

Provides memory-efficient sequence access using pyfaidx with automatic fallback.
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

    _global_cache = {}
    _max_cache_size = 3

    def __init__(self, assembly_file, use_cache=True):
        self.assembly_file = Path(assembly_file)
        self.use_cache = use_cache
        self.logger = logging.getLogger(__name__)

        if not self.assembly_file.exists():
            raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

        self.canonical_path = str(self.assembly_file.resolve())
        self.fasta_index = None
        self._seqio_cache = {}

        if HAVE_PYFAIDX:
            self._init_pyfaidx()
        else:
            self.logger.warning(
                "pyfaidx not available. Install with 'pip install pyfaidx' for "
                "10-1000x faster sequence access."
            )

    def _init_pyfaidx(self):
        try:
            if self.use_cache and self.canonical_path in self._global_cache:
                self.fasta_index = self._global_cache[self.canonical_path]
                return

            self.fasta_index = Fasta(str(self.assembly_file))

            if self.use_cache:
                if len(self._global_cache) >= self._max_cache_size:
                    oldest_key = next(iter(self._global_cache))
                    old_index = self._global_cache.pop(oldest_key)
                    try:
                        old_index.close()
                    except:
                        pass
                self._global_cache[self.canonical_path] = self.fasta_index

        except Exception as e:
            self.logger.warning(f"Failed to create pyfaidx index: {e}")
            self.fasta_index = None

    def get_sequence(self, chrom, start, end):
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return str(self.fasta_index[chrom][start:end])
            except Exception:
                pass
        return self._get_sequence_seqio(chrom, start, end)

    def _get_sequence_seqio(self, chrom, start, end):
        if chrom in self._seqio_cache:
            seq = self._seqio_cache[chrom]
            if end is None:
                return seq
            return seq[start:end]

        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            if record.id == chrom:
                seq_str = str(record.seq)
                if len(seq_str) < 100_000_000:
                    self._seqio_cache[chrom] = seq_str
                if end is None:
                    return seq_str
                return seq_str[start:end]

        raise ValueError(f"Chromosome {chrom} not found")

    def get_full_sequence(self, chrom):
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return str(self.fasta_index[chrom][:])
            except Exception:
                pass
        return self._get_sequence_seqio(chrom, 0, None)

    def get_length(self, chrom):
        if self.fasta_index:
            try:
                if chrom in self.fasta_index:
                    return len(self.fasta_index[chrom])
            except:
                pass

        for record in SeqIO.parse(self.assembly_file, 'fasta'):
            if record.id == chrom:
                return len(record.seq)

        raise ValueError(f"Chromosome {chrom} not found")

    def get_all_chroms(self):
        if self.fasta_index:
            try:
                return list(self.fasta_index.keys())
            except:
                pass

        return [record.id for record in SeqIO.parse(self.assembly_file, 'fasta')]

    def close(self):
        if self.fasta_index and not self.use_cache:
            try:
                self.fasta_index.close()
            except:
                pass
            self.fasta_index = None
        self._seqio_cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    @classmethod
    def clear_global_cache(cls):
        for index in cls._global_cache.values():
            try:
                index.close()
            except:
                pass
        cls._global_cache.clear()
