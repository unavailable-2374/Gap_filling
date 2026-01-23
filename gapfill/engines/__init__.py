"""
Gap filling engines for different ploidy levels
"""
from gapfill.engines.haploid import HaploidEngine
from gapfill.engines.polyploid import PolyploidEngine

__all__ = ["HaploidEngine", "PolyploidEngine"]
