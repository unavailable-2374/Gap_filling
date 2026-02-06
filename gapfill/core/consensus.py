#!/usr/bin/env python3
"""
Consensus-based gap filling - Skip wtdbg2 for high-quality cases

When spanning reads are highly consistent (>95% identity), we can
directly compute a consensus sequence instead of running a full
assembler. This is much faster and often more accurate for simple gaps.

Strategies:
1. Direct consensus: For highly consistent spanning reads
2. POA (Partial Order Alignment): For moderately consistent reads
3. Fall back to wtdbg2: For complex cases
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from collections import Counter

import numpy as np


class ConsensusBuilder:
    """
    Build consensus sequences from reads.

    Optimized for gap filling where we have a small number of
    high-quality spanning reads.
    """

    # Thresholds for consensus strategies
    MIN_READS_DIRECT = 3
    MIN_READS_POA = 3
    MIN_IDENTITY_DIRECT = 0.95
    MIN_IDENTITY_POA = 0.85

    def __init__(self, threads: int = 8):
        self.threads = threads
        self.logger = logging.getLogger(__name__)

        # Check for external tools
        self.has_spoa = self._check_tool('spoa')
        self.has_abpoa = self._check_tool('abpoa')

    def _check_tool(self, tool: str) -> bool:
        """Check if tool is available"""
        try:
            subprocess.run([tool, '--version'], capture_output=True, timeout=5)
            return True
        except:
            return False

    def build_consensus(self, reads: List[Tuple[str, str, str]],
                        gap_info: Optional[Dict] = None,
                        read_type: str = 'hifi') -> Optional[Dict]:
        """
        Attempt to build consensus from reads.

        Args:
            reads: List of (sequence, name, type) tuples
            gap_info: Optional gap information for size estimation
            read_type: 'hifi' or 'ont' — affects strategy selection

        Returns:
            Dict with 'sequence', 'method', 'confidence' or None if failed
        """
        if len(reads) < self.MIN_READS_DIRECT:
            return None

        sequences = [r[0] for r in reads]

        # Calculate pairwise identity
        identity = self._estimate_identity(sequences)
        self.logger.debug(f"  Consensus: {len(reads)} {read_type} reads, identity={identity:.2%}")

        if read_type == 'ont':
            # ONT reads: higher error rate, only use POA (no direct consensus)
            if identity >= self.MIN_IDENTITY_POA and len(reads) >= 5:
                consensus = self._poa_consensus(sequences)
                if consensus:
                    return {
                        'sequence': consensus,
                        'method': 'ont_poa_consensus',
                        'confidence': identity,
                        'read_count': len(reads),
                        'identity': identity
                    }
            return None

        # HiFi reads: use both strategies
        # Strategy 1: Direct consensus for highly consistent reads
        if identity >= self.MIN_IDENTITY_DIRECT and len(reads) >= self.MIN_READS_DIRECT:
            consensus = self._direct_consensus(sequences)
            if consensus:
                return {
                    'sequence': consensus,
                    'method': 'direct_consensus',
                    'confidence': min(1.0, identity + 0.1),
                    'read_count': len(reads),
                    'identity': identity
                }

        # Strategy 2: POA for moderately consistent reads
        if identity >= self.MIN_IDENTITY_POA and len(reads) >= self.MIN_READS_POA:
            consensus = self._poa_consensus(sequences)
            if consensus:
                return {
                    'sequence': consensus,
                    'method': 'poa_consensus',
                    'confidence': identity,
                    'read_count': len(reads),
                    'identity': identity
                }

        # Cannot build consensus - need assembly
        return None

    def _estimate_identity(self, sequences: List[str]) -> float:
        """
        Estimate average pairwise identity of sequences.

        Uses a sampling approach for efficiency.
        """
        if len(sequences) < 2:
            return 1.0

        # Sample pairs for efficiency
        n_samples = min(10, len(sequences) * (len(sequences) - 1) // 2)
        identities = []

        import random
        pairs = []
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                pairs.append((i, j))

        if len(pairs) > n_samples:
            pairs = random.sample(pairs, n_samples)

        for i, j in pairs:
            identity = self._pairwise_identity(sequences[i], sequences[j])
            identities.append(identity)

        return np.mean(identities) if identities else 0.0

    def _pairwise_identity(self, seq1: str, seq2: str) -> float:
        """
        Calculate pairwise identity between two sequences.

        Uses a fast k-mer based approximation.
        """
        k = 15
        if len(seq1) < k or len(seq2) < k:
            return 0.0

        # Get k-mers
        kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
        kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))

        # Jaccard similarity as identity proxy
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)

        if union == 0:
            return 0.0

        return intersection / union

    def _direct_consensus(self, sequences: List[str]) -> Optional[str]:
        """
        Build consensus by direct alignment and majority vote.

        For highly consistent sequences (>95% identity).
        """
        if not sequences:
            return None

        # Use the median-length sequence as reference
        sorted_seqs = sorted(sequences, key=len)
        ref_seq = sorted_seqs[len(sorted_seqs) // 2]

        # For very high consistency, just return the reference
        if len(sequences) <= 3:
            return ref_seq

        # Simple majority vote at each position
        # This works well when sequences are very similar
        try:
            # Align all sequences to reference
            aligned = self._simple_align_to_ref(ref_seq, sequences)
            if not aligned:
                return ref_seq

            # Build consensus from aligned sequences
            consensus = self._majority_vote_consensus(aligned)
            return consensus if consensus else ref_seq

        except Exception as e:
            self.logger.debug(f"Direct consensus failed: {e}")
            return ref_seq

    def _simple_align_to_ref(self, ref: str, sequences: List[str]) -> Optional[List[str]]:
        """
        Simple alignment of sequences to reference.

        Uses minimap2 for fast alignment.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write reference
            ref_fa = tmpdir / "ref.fa"
            with open(ref_fa, 'w') as f:
                f.write(f">ref\n{ref}\n")

            # Write sequences
            seqs_fa = tmpdir / "seqs.fa"
            with open(seqs_fa, 'w') as f:
                for i, seq in enumerate(sequences):
                    f.write(f">seq{i}\n{seq}\n")

            # Align with minimap2
            try:
                result = subprocess.run(
                    ['minimap2', '-ax', 'map-hifi', '-t', '1', str(ref_fa), str(seqs_fa)],
                    capture_output=True, text=True, timeout=60
                )
                if result.returncode != 0:
                    return None

                # Parse SAM output to get aligned sequences
                aligned = []
                for line in result.stdout.split('\n'):
                    if line.startswith('@') or not line.strip():
                        continue
                    fields = line.split('\t')
                    if len(fields) >= 10:
                        seq = fields[9]
                        if seq and seq != '*':
                            aligned.append(seq)

                return aligned if len(aligned) >= 2 else None

            except:
                return None

    def _majority_vote_consensus(self, sequences: List[str]) -> Optional[str]:
        """
        Build consensus by majority vote at each position.
        """
        if not sequences:
            return None

        max_len = max(len(s) for s in sequences)
        consensus = []

        for i in range(max_len):
            bases = []
            for seq in sequences:
                if i < len(seq):
                    bases.append(seq[i])

            if bases:
                # Majority vote
                counter = Counter(bases)
                most_common = counter.most_common(1)[0][0]
                consensus.append(most_common)

        return ''.join(consensus)

    def _poa_consensus(self, sequences: List[str]) -> Optional[str]:
        """
        Build consensus using Partial Order Alignment (POA).

        Uses spoa or abpoa if available.
        """
        if not sequences:
            return None

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Write sequences
            seqs_fa = tmpdir / "seqs.fa"
            with open(seqs_fa, 'w') as f:
                for i, seq in enumerate(sequences):
                    f.write(f">seq{i}\n{seq}\n")

            try:
                if self.has_abpoa:
                    # Use abpoa (faster)
                    result = subprocess.run(
                        ['abpoa', '-r', '1', str(seqs_fa)],
                        capture_output=True, text=True, timeout=120
                    )
                elif self.has_spoa:
                    # Use spoa
                    result = subprocess.run(
                        ['spoa', '-r', '1', str(seqs_fa)],
                        capture_output=True, text=True, timeout=120
                    )
                else:
                    # Fall back to direct consensus
                    return self._direct_consensus(sequences)

                if result.returncode != 0:
                    return self._direct_consensus(sequences)

                # Parse FASTA output
                consensus = ''
                for line in result.stdout.split('\n'):
                    if not line.startswith('>'):
                        consensus += line.strip()

                # Sanity check
                avg_len = np.mean([len(s) for s in sequences])
                if consensus and 0.8 <= len(consensus) / avg_len <= 1.2:
                    return consensus

                return self._direct_consensus(sequences)

            except Exception as e:
                self.logger.debug(f"POA consensus failed: {e}")
                return self._direct_consensus(sequences)


def try_consensus_fill(reads: List[Tuple[str, str, str]],
                       gap_info: Dict,
                       threads: int = 8,
                       read_type: str = 'hifi') -> Optional[Dict]:
    """
    Try to fill a gap using consensus-based approach.

    This is faster than wtdbg2 for simple gaps with consistent reads.

    Args:
        reads: List of (sequence, name, type) tuples
        gap_info: Gap information dict
        threads: Number of threads
        read_type: 'hifi' or 'ont'

    Returns:
        Dict with fill result or None if consensus approach failed
    """
    builder = ConsensusBuilder(threads=threads)
    result = builder.build_consensus(reads, gap_info, read_type=read_type)

    if result:
        return {
            'success': True,
            'sequence': result['sequence'],
            'strategy': 'consensus',
            'source': result['method'],
            'read_count': result['read_count'],
            'confidence': result['confidence'],
            'is_complete': True,
            'has_placeholder': False
        }

    return None
