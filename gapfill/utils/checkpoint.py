#!/usr/bin/env python3
"""
Checkpoint module for gap filling state persistence.

Enables resume functionality after interruption.

Checkpoint file structure:
{
    "version": "1.0",
    "engine": "haploid|polyploid|optimized_polyploid",
    "phase": "normalization|phasing|filling",
    "iteration": 3,
    "completed_gaps": {"gap_name": "sequence", ...},
    "failed_gaps": ["gap_name", ...],
    "current_assembly": "path/to/assembly.fasta",
    "intermediate_files": {
        "normalized_assembly": "path",
        "hifi_bam": "path",
        "ont_bam": "path",
        "snp_database": "path",
        "phased_reads": {"hap1_hifi": "path", ...}
    },
    "timestamp": "2024-01-01T00:00:00"
}
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass, field, asdict


@dataclass
class CheckpointState:
    """Represents the state of a gap filling run"""
    version: str = "1.0"
    engine: str = ""
    phase: str = "init"  # init, normalization, phasing, filling, complete
    iteration: int = 0
    max_iterations: int = 10

    # Gap tracking
    completed_gaps: Dict[str, str] = field(default_factory=dict)  # gap_name -> sequence
    partially_filled_gaps: Dict[str, str] = field(default_factory=dict)
    failed_gaps: List[str] = field(default_factory=list)

    # Current state
    current_assembly: str = ""

    # Intermediate files (to check for reuse)
    intermediate_files: Dict[str, str] = field(default_factory=dict)

    # Polyploid specific
    haplotype_states: Dict[str, Dict] = field(default_factory=dict)
    snp_database_path: str = ""
    phased_reads: Dict[str, str] = field(default_factory=dict)

    timestamp: str = ""

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization"""
        d = asdict(self)
        # Convert sets to lists if any
        return d

    @classmethod
    def from_dict(cls, data: Dict) -> 'CheckpointState':
        """Create from dictionary"""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


class CheckpointManager:
    """Manages checkpoint saving and loading for gap filling"""

    CHECKPOINT_FILE = "checkpoint.json"

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.checkpoint_file = self.output_dir / self.CHECKPOINT_FILE
        self.logger = logging.getLogger(__name__)
        self.state: Optional[CheckpointState] = None

    def exists(self) -> bool:
        """Check if a checkpoint exists"""
        return self.checkpoint_file.exists()

    def load(self) -> Optional[CheckpointState]:
        """Load checkpoint from file"""
        if not self.exists():
            self.logger.info("No checkpoint found")
            return None

        try:
            with open(self.checkpoint_file) as f:
                data = json.load(f)

            self.state = CheckpointState.from_dict(data)
            self.logger.info(f"Loaded checkpoint: phase={self.state.phase}, "
                           f"iteration={self.state.iteration}")
            return self.state

        except Exception as e:
            self.logger.error(f"Failed to load checkpoint: {e}")
            return None

    def save(self, state: CheckpointState):
        """Save checkpoint to file"""
        state.timestamp = datetime.now().isoformat()
        self.state = state

        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)

            with open(self.checkpoint_file, 'w') as f:
                json.dump(state.to_dict(), f, indent=2)

            self.logger.debug(f"Checkpoint saved: phase={state.phase}, "
                            f"iteration={state.iteration}")

        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")

    def update_phase(self, phase: str):
        """Update current phase"""
        if self.state:
            self.state.phase = phase
            self.save(self.state)

    def update_iteration(self, iteration: int):
        """Update current iteration"""
        if self.state:
            self.state.iteration = iteration
            self.save(self.state)

    def add_completed_gap(self, gap_name: str, sequence: str):
        """Record a completely filled gap"""
        if self.state:
            self.state.completed_gaps[gap_name] = sequence
            # Remove from partial/failed if present
            self.state.partially_filled_gaps.pop(gap_name, None)
            if gap_name in self.state.failed_gaps:
                self.state.failed_gaps.remove(gap_name)
            self.save(self.state)

    def add_partial_gap(self, gap_name: str, sequence: str):
        """Record a partially filled gap"""
        if self.state:
            if gap_name not in self.state.completed_gaps:
                self.state.partially_filled_gaps[gap_name] = sequence
            self.save(self.state)

    def add_failed_gap(self, gap_name: str):
        """Record a failed gap"""
        if self.state:
            if gap_name not in self.state.failed_gaps:
                self.state.failed_gaps.append(gap_name)
            self.save(self.state)

    def set_intermediate_file(self, key: str, path: str):
        """Record an intermediate file path"""
        if self.state:
            self.state.intermediate_files[key] = path
            self.save(self.state)

    def get_intermediate_file(self, key: str) -> Optional[Path]:
        """Get an intermediate file if it exists and is valid"""
        if not self.state:
            return None

        path_str = self.state.intermediate_files.get(key)
        if not path_str:
            return None

        path = Path(path_str)
        if path.exists() and path.stat().st_size > 0:
            return path
        return None

    def set_current_assembly(self, path: str):
        """Update current assembly path"""
        if self.state:
            self.state.current_assembly = path
            self.save(self.state)

    def mark_complete(self):
        """Mark the run as complete"""
        if self.state:
            self.state.phase = "complete"
            self.save(self.state)

    def clear(self):
        """Remove checkpoint file"""
        if self.checkpoint_file.exists():
            self.checkpoint_file.unlink()
            self.logger.info("Checkpoint cleared")
        self.state = None

    def can_skip_normalization(self) -> bool:
        """Check if normalization can be skipped"""
        if not self.state:
            return False

        # Check if we're past normalization phase
        if self.state.phase in ('phasing', 'filling', 'complete'):
            normalized = self.get_intermediate_file('normalized_assembly')
            if normalized:
                self.logger.info(f"Reusing normalized assembly: {normalized}")
                return True
        return False

    def can_skip_phasing(self) -> bool:
        """Check if phasing can be skipped (polyploid)"""
        if not self.state:
            return False

        if self.state.phase in ('filling', 'complete'):
            # Check if phased reads exist
            if self.state.phased_reads:
                all_exist = all(
                    Path(p).exists()
                    for p in self.state.phased_reads.values()
                    if p
                )
                if all_exist:
                    self.logger.info("Reusing phased reads from checkpoint")
                    return True
        return False

    def get_resume_iteration(self) -> int:
        """Get the iteration to resume from"""
        if not self.state:
            return 0

        # Resume from current iteration (re-run it to be safe)
        return max(0, self.state.iteration - 1)

    def get_completed_gap_names(self) -> Set[str]:
        """Get set of completed gap names"""
        if not self.state:
            return set()
        return set(self.state.completed_gaps.keys())

    def get_failed_gap_names(self) -> Set[str]:
        """Get set of failed gap names"""
        if not self.state:
            return set()
        return set(self.state.failed_gaps)


class PolyploidCheckpointManager(CheckpointManager):
    """Extended checkpoint manager for polyploid engines"""

    def __init__(self, output_dir: str, hap_names: List[str]):
        super().__init__(output_dir)
        self.hap_names = hap_names

    def init_haplotype_states(self):
        """Initialize per-haplotype states"""
        if self.state:
            for hap in self.hap_names:
                if hap not in self.state.haplotype_states:
                    self.state.haplotype_states[hap] = {
                        'completed_gaps': {},
                        'failed_gaps': [],
                        'current_assembly': ''
                    }
            self.save(self.state)

    def add_completed_gap_for_hap(self, hap_name: str, gap_name: str, sequence: str):
        """Record a completed gap for a specific haplotype"""
        if self.state and hap_name in self.state.haplotype_states:
            self.state.haplotype_states[hap_name]['completed_gaps'][gap_name] = sequence
            self.save(self.state)

    def add_failed_gap_for_hap(self, hap_name: str, gap_name: str):
        """Record a failed gap for a specific haplotype"""
        if self.state and hap_name in self.state.haplotype_states:
            failed = self.state.haplotype_states[hap_name]['failed_gaps']
            if gap_name not in failed:
                failed.append(gap_name)
            self.save(self.state)

    def set_hap_assembly(self, hap_name: str, path: str):
        """Set current assembly for a haplotype"""
        if self.state and hap_name in self.state.haplotype_states:
            self.state.haplotype_states[hap_name]['current_assembly'] = path
            self.save(self.state)

    def get_hap_completed_gaps(self, hap_name: str) -> Set[str]:
        """Get completed gaps for a haplotype"""
        if not self.state or hap_name not in self.state.haplotype_states:
            return set()
        return set(self.state.haplotype_states[hap_name].get('completed_gaps', {}).keys())

    def get_hap_failed_gaps(self, hap_name: str) -> Set[str]:
        """Get failed gaps for a haplotype"""
        if not self.state or hap_name not in self.state.haplotype_states:
            return set()
        return set(self.state.haplotype_states[hap_name].get('failed_gaps', []))

    def set_snp_database(self, path: str):
        """Set SNP database path"""
        if self.state:
            self.state.snp_database_path = path
            self.save(self.state)

    def get_snp_database(self) -> Optional[Path]:
        """Get SNP database if exists"""
        if not self.state or not self.state.snp_database_path:
            return None
        path = Path(self.state.snp_database_path)
        if path.exists():
            return path
        return None

    def set_phased_reads(self, phased_reads: Dict[str, str]):
        """Set phased reads paths"""
        if self.state:
            self.state.phased_reads = phased_reads
            self.save(self.state)

    def get_phased_reads(self) -> Dict[str, Path]:
        """Get phased reads that exist"""
        if not self.state:
            return {}

        result = {}
        for key, path_str in self.state.phased_reads.items():
            if path_str:
                path = Path(path_str)
                if path.exists() and path.stat().st_size > 0:
                    result[key] = path
        return result
