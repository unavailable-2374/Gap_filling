"""
GapFill - Gap filling pipeline for genome assemblies

Supports both haploid and polyploid genomes with HiFi/ONT long reads.

Usage:
    python -m gapfill -a assembly.fa --hifi reads.fq.gz -o output

    # Or as library
    from gapfill import HaploidEngine, PolyploidEngine
"""

__version__ = "1.0.0"
__author__ = "Gap Filling Pipeline"

from gapfill.engines.haploid import HaploidEngine
from gapfill.engines.polyploid import PolyploidEngine
from gapfill.core.filler import GapFiller
from gapfill.core.validator import GapValidator

__all__ = [
    "HaploidEngine",
    "PolyploidEngine",
    "GapFiller",
    "GapValidator",
]
