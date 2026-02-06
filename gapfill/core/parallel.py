#!/usr/bin/env python3
"""
Parallel Gap Filling - Process multiple gaps concurrently

Uses multiprocessing to fill gaps in parallel, with proper
resource management and error handling.
"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Callable
from dataclasses import dataclass
import traceback
import os


@dataclass
class ParallelConfig:
    """Configuration for parallel gap filling"""
    max_workers: int = 8
    batch_size: int = 50  # Process gaps in batches
    timeout_per_gap: int = 300  # 5 minutes per gap
    use_consensus_first: bool = True


def _fill_single_gap(args: tuple) -> Dict:
    """
    Worker function to fill a single gap.

    This runs in a separate process, so we need to recreate
    the GapFiller instance with the provided configuration.
    """
    gap, config = args

    try:
        # Import here to avoid multiprocessing issues
        from gapfill.core.filler import GapFiller
        from gapfill.core.consensus import try_consensus_fill

        # Create filler instance
        filler = GapFiller(
            assembly_file=config['assembly_file'],
            hifi_bam=config.get('hifi_bam'),
            ont_bam=config.get('ont_bam'),
            hifi_reads=config.get('hifi_reads'),
            ont_reads=config.get('ont_reads'),
            threads=1,  # Single thread per worker
            work_dir=config['work_dir'],
            min_mapq=config.get('min_mapq', 20),
            enable_validation=True,  # Local validation runs per worker
        )

        # Fill the gap
        result = filler.fill_gap(gap)
        filler.close()

        return {
            'gap_name': gap['name'],
            'result': result,
            'error': None
        }

    except Exception as e:
        return {
            'gap_name': gap.get('name', 'unknown'),
            'result': None,
            'error': str(e),
            'traceback': traceback.format_exc()
        }


def fill_gaps_parallel(gaps: List[Dict],
                       assembly_file: str,
                       hifi_bam: Optional[str] = None,
                       ont_bam: Optional[str] = None,
                       hifi_reads: Optional[str] = None,
                       ont_reads: Optional[str] = None,
                       work_dir: str = "work",
                       threads: int = 8,
                       min_mapq: int = 20,
                       use_consensus_first: bool = True,
                       progress_callback: Optional[Callable] = None) -> Dict[str, Dict]:
    """
    Fill multiple gaps in parallel.

    Args:
        gaps: List of gap dictionaries
        assembly_file: Assembly FASTA path
        hifi_bam: HiFi BAM path
        ont_bam: ONT BAM path
        hifi_reads: HiFi reads path (for assembly)
        ont_reads: ONT reads path (for assembly)
        work_dir: Working directory
        threads: Total threads to use
        min_mapq: Minimum mapping quality
        use_consensus_first: Try consensus before wtdbg2
        progress_callback: Optional callback for progress updates

    Returns:
        Dict mapping gap_name to fill result
    """
    logger = logging.getLogger(__name__)

    if not gaps:
        return {}

    # Determine number of workers
    # Each worker gets 1 thread, so max workers = total threads
    # But also limit to avoid too many processes
    max_workers = min(threads, len(gaps), 16)

    logger.info(f"Filling {len(gaps)} gaps with {max_workers} parallel workers")

    # Prepare configuration for workers
    config = {
        'assembly_file': assembly_file,
        'hifi_bam': hifi_bam,
        'ont_bam': ont_bam,
        'hifi_reads': hifi_reads,
        'ont_reads': ont_reads,
        'work_dir': work_dir,
        'min_mapq': min_mapq,
        'use_consensus_first': use_consensus_first,
    }

    # Prepare work items
    work_items = [(gap, config) for gap in gaps]

    # Results storage
    results = {}
    completed = 0

    # Use ProcessPoolExecutor for parallel execution
    # Note: We use 'spawn' context to avoid fork issues with pysam
    try:
        ctx = mp.get_context('spawn')
        with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as executor:
            # Submit all jobs
            future_to_gap = {
                executor.submit(_fill_single_gap, item): item[0]['name']
                for item in work_items
            }

            # Collect results as they complete
            for future in as_completed(future_to_gap):
                gap_name = future_to_gap[future]

                try:
                    result = future.result(timeout=300)
                    results[gap_name] = result

                    if result['error']:
                        logger.warning(f"  Gap {gap_name} failed: {result['error']}")
                    else:
                        fill_result = result['result']
                        if fill_result and fill_result.get('success'):
                            logger.debug(f"  Gap {gap_name}: filled ({fill_result.get('source', 'unknown')})")
                        else:
                            logger.debug(f"  Gap {gap_name}: not filled")

                except Exception as e:
                    logger.error(f"  Gap {gap_name} exception: {e}")
                    results[gap_name] = {
                        'gap_name': gap_name,
                        'result': None,
                        'error': str(e)
                    }

                completed += 1

                if progress_callback:
                    progress_callback(completed, len(gaps))

                if completed % 100 == 0:
                    logger.info(f"  Progress: {completed}/{len(gaps)} gaps processed")

    except Exception as e:
        logger.error(f"Parallel execution failed: {e}")
        # Fall back to sequential processing
        logger.info("Falling back to sequential processing...")
        return fill_gaps_sequential(gaps, assembly_file, hifi_bam, ont_bam,
                                   hifi_reads, ont_reads, work_dir, threads,
                                   min_mapq, use_consensus_first)

    logger.info(f"Parallel filling complete: {len(results)} gaps processed")

    # Count successes
    successes = sum(1 for r in results.values()
                   if r.get('result') and r['result'].get('success'))
    logger.info(f"  Successful fills: {successes}/{len(results)}")

    return results


def fill_gaps_sequential(gaps: List[Dict],
                         assembly_file: str,
                         hifi_bam: Optional[str] = None,
                         ont_bam: Optional[str] = None,
                         hifi_reads: Optional[str] = None,
                         ont_reads: Optional[str] = None,
                         work_dir: str = "work",
                         threads: int = 8,
                         min_mapq: int = 20,
                         use_consensus_first: bool = True) -> Dict[str, Dict]:
    """
    Fill gaps sequentially (fallback for parallel).
    """
    from gapfill.core.filler import GapFiller
    from gapfill.core.consensus import try_consensus_fill

    logger = logging.getLogger(__name__)
    logger.info(f"Filling {len(gaps)} gaps sequentially...")

    # Create single filler instance
    filler = GapFiller(
        assembly_file=assembly_file,
        hifi_bam=hifi_bam,
        ont_bam=ont_bam,
        hifi_reads=hifi_reads,
        ont_reads=ont_reads,
        threads=threads,
        work_dir=work_dir,
        min_mapq=min_mapq,
        enable_validation=True,  # Local validation runs per fill
    )

    results = {}

    for i, gap in enumerate(gaps):
        gap_name = gap.get('name', f"{gap['chrom']}_{gap['start']}_{gap['end']}")

        try:
            result = filler.fill_gap(gap)
            results[gap_name] = {
                'gap_name': gap_name,
                'result': result,
                'error': None
            }

        except Exception as e:
            logger.warning(f"  Gap {gap_name} failed: {e}")
            results[gap_name] = {
                'gap_name': gap_name,
                'result': None,
                'error': str(e)
            }

        if (i + 1) % 100 == 0:
            logger.info(f"  Progress: {i + 1}/{len(gaps)} gaps processed")

    filler.close()

    # Count successes
    successes = sum(1 for r in results.values()
                   if r.get('result') and r['result'].get('success'))
    logger.info(f"  Successful fills: {successes}/{len(results)}")

    return results


class GapBatcher:
    """
    Batch gaps for efficient parallel processing.

    Groups gaps by chromosome and estimated complexity.
    """

    def __init__(self, batch_size: int = 50):
        self.batch_size = batch_size
        self.logger = logging.getLogger(__name__)

    def create_batches(self, gaps: List[Dict],
                       priority_scores: Optional[Dict[str, float]] = None) -> List[List[Dict]]:
        """
        Create batches of gaps for processing.

        Args:
            gaps: List of gaps
            priority_scores: Optional dict of gap_name -> priority score
                            (higher = process first)

        Returns:
            List of gap batches
        """
        if not gaps:
            return []

        # Sort by priority if provided
        if priority_scores:
            sorted_gaps = sorted(
                gaps,
                key=lambda g: priority_scores.get(g.get('name', ''), 0),
                reverse=True
            )
        else:
            # Default: sort by gap size (smaller first = easier)
            sorted_gaps = sorted(gaps, key=lambda g: g.get('size', 0))

        # Create batches
        batches = []
        for i in range(0, len(sorted_gaps), self.batch_size):
            batch = sorted_gaps[i:i + self.batch_size]
            batches.append(batch)

        self.logger.info(f"Created {len(batches)} batches from {len(gaps)} gaps")

        return batches

    def estimate_priority(self, gap: Dict, reads_info: Optional[Dict] = None) -> float:
        """
        Estimate gap filling priority.

        Higher priority = more likely to fill successfully.
        """
        score = 0.0

        # Small gaps are easier
        size = gap.get('size', 500)
        if size < 1000:
            score += 50
        elif size < 5000:
            score += 30
        elif size < 10000:
            score += 10

        # Reads info if available
        if reads_info:
            # HiFi spanning is best
            if reads_info.get('hifi_spanning', 0) >= 5:
                score += 100
            elif reads_info.get('hifi_spanning', 0) >= 3:
                score += 50

            # ONT spanning
            if reads_info.get('ont_spanning', 0) >= 3:
                score += 30

            # Flanking reads
            if (reads_info.get('hifi_left', 0) >= 3 and
                reads_info.get('hifi_right', 0) >= 3):
                score += 20

        return score
