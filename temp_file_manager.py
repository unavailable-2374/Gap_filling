#!/usr/bin/env python3
"""
Temporary File Manager

Provides reliable management and cleanup of temporary files used during gap filling.
Ensures temp files are properly cleaned up even when errors occur.
"""

import os
import tempfile
import logging
import atexit
from pathlib import Path
from contextlib import contextmanager
from typing import Optional, Set, Union

logger = logging.getLogger(__name__)


class TempFileManager:
    """
    Manages temporary files with automatic cleanup.

    Features:
    - Context manager support for automatic cleanup
    - Tracks all created temp files
    - Cleanup on exit (atexit handler)
    - Optional prefix/suffix for temp files
    - Configurable temp directory
    """

    def __init__(self,
                 work_dir: Optional[Union[str, Path]] = None,
                 base_dir: Optional[Union[str, Path]] = None,
                 prefix: str = "gapfill_",
                 auto_cleanup: bool = True):
        """
        Initialize TempFileManager.

        Args:
            work_dir: Working directory for temp files (alias for base_dir)
            base_dir: Base directory for temp files (default: system temp)
            prefix: Prefix for temp file names
            auto_cleanup: Register atexit handler for cleanup
        """
        # work_dir is an alias for base_dir
        effective_dir = work_dir if work_dir else base_dir
        self.base_dir = Path(effective_dir) if effective_dir else None
        self.prefix = prefix
        self.tracked_files: Set[Path] = set()
        self.tracked_dirs: Set[Path] = set()

        # Create base directory if specified
        if self.base_dir:
            self.base_dir.mkdir(parents=True, exist_ok=True)

        # Register cleanup on exit
        if auto_cleanup:
            atexit.register(self.cleanup)

    @contextmanager
    def temp_file(self,
                  suffix: str = "",
                  prefix: Optional[str] = None,
                  delete: bool = False):
        """
        Context manager for creating a temporary file.

        Args:
            suffix: File suffix (e.g., ".fasta", ".bam")
            prefix: File prefix (overrides default)
            delete: Delete file when context exits (default: False, cleanup on manager cleanup)

        Yields:
            Path to temporary file

        Example:
            with temp_manager.temp_file(suffix=".fasta") as temp_path:
                # Use temp_path
                with open(temp_path, 'w') as f:
                    f.write(">seq1\nACGT\n")
        """
        file_prefix = prefix if prefix else self.prefix

        # Create temp file
        fd, temp_path = tempfile.mkstemp(
            suffix=suffix,
            prefix=file_prefix,
            dir=self.base_dir
        )
        os.close(fd)  # Close file descriptor, we just need the path

        temp_path = Path(temp_path)
        self.tracked_files.add(temp_path)

        try:
            yield temp_path
        finally:
            if delete:
                self._safe_delete_file(temp_path)
                self.tracked_files.discard(temp_path)

    @contextmanager
    def temp_dir(self,
                 suffix: str = "",
                 prefix: Optional[str] = None,
                 delete: bool = False):
        """
        Context manager for creating a temporary directory.

        Args:
            suffix: Directory suffix
            prefix: Directory prefix (overrides default)
            delete: Delete directory when context exits

        Yields:
            Path to temporary directory
        """
        dir_prefix = prefix if prefix else self.prefix

        temp_dir = tempfile.mkdtemp(
            suffix=suffix,
            prefix=dir_prefix,
            dir=self.base_dir
        )
        temp_dir = Path(temp_dir)
        self.tracked_dirs.add(temp_dir)

        try:
            yield temp_dir
        finally:
            if delete:
                self._safe_delete_dir(temp_dir)
                self.tracked_dirs.discard(temp_dir)

    def create_temp_file(self,
                         suffix: str = "",
                         prefix: Optional[str] = None) -> Path:
        """
        Create a temporary file without context manager.
        File will be tracked and cleaned up on manager cleanup.

        Args:
            suffix: File suffix
            prefix: File prefix

        Returns:
            Path to temporary file
        """
        file_prefix = prefix if prefix else self.prefix

        fd, temp_path = tempfile.mkstemp(
            suffix=suffix,
            prefix=file_prefix,
            dir=self.base_dir
        )
        os.close(fd)

        temp_path = Path(temp_path)
        self.tracked_files.add(temp_path)

        return temp_path

    def create_temp_dir(self,
                        suffix: str = "",
                        prefix: Optional[str] = None) -> Path:
        """
        Create a temporary directory without context manager.
        Directory will be tracked and cleaned up on manager cleanup.

        Args:
            suffix: Directory suffix
            prefix: Directory prefix

        Returns:
            Path to temporary directory
        """
        dir_prefix = prefix if prefix else self.prefix

        temp_dir = tempfile.mkdtemp(
            suffix=suffix,
            prefix=dir_prefix,
            dir=self.base_dir
        )
        temp_dir = Path(temp_dir)
        self.tracked_dirs.add(temp_dir)

        return temp_dir

    def track_file(self, file_path: Union[str, Path]):
        """
        Add an existing file to be tracked for cleanup.

        Args:
            file_path: Path to file to track
        """
        self.tracked_files.add(Path(file_path))

    def track_dir(self, dir_path: Union[str, Path]):
        """
        Add an existing directory to be tracked for cleanup.

        Args:
            dir_path: Path to directory to track
        """
        self.tracked_dirs.add(Path(dir_path))

    def untrack_file(self, file_path: Union[str, Path]):
        """
        Remove a file from tracking (will not be deleted on cleanup).

        Args:
            file_path: Path to file to untrack
        """
        self.tracked_files.discard(Path(file_path))

    def untrack_dir(self, dir_path: Union[str, Path]):
        """
        Remove a directory from tracking.

        Args:
            dir_path: Path to directory to untrack
        """
        self.tracked_dirs.discard(Path(dir_path))

    def delete_file(self, file_path: Union[str, Path]):
        """
        Delete a specific tracked file immediately.

        Args:
            file_path: Path to file to delete
        """
        file_path = Path(file_path)
        self._safe_delete_file(file_path)
        self.tracked_files.discard(file_path)

    def delete_dir(self, dir_path: Union[str, Path]):
        """
        Delete a specific tracked directory immediately.

        Args:
            dir_path: Path to directory to delete
        """
        dir_path = Path(dir_path)
        self._safe_delete_dir(dir_path)
        self.tracked_dirs.discard(dir_path)

    def cleanup(self):
        """
        Clean up all tracked temporary files and directories.
        Called automatically on exit if auto_cleanup=True.
        """
        # Delete files first
        deleted_files = 0
        for file_path in list(self.tracked_files):
            if self._safe_delete_file(file_path):
                deleted_files += 1
            self.tracked_files.discard(file_path)

        # Then delete directories
        deleted_dirs = 0
        for dir_path in list(self.tracked_dirs):
            if self._safe_delete_dir(dir_path):
                deleted_dirs += 1
            self.tracked_dirs.discard(dir_path)

        if deleted_files > 0 or deleted_dirs > 0:
            logger.debug(f"TempFileManager cleanup: {deleted_files} files, {deleted_dirs} directories")

    def _safe_delete_file(self, file_path: Path) -> bool:
        """
        Safely delete a file, handling errors gracefully.

        Returns:
            True if file was deleted, False otherwise
        """
        try:
            if file_path.exists():
                file_path.unlink()
                return True
        except Exception as e:
            logger.warning(f"Failed to delete temp file {file_path}: {e}")
        return False

    def _safe_delete_dir(self, dir_path: Path) -> bool:
        """
        Safely delete a directory and its contents.

        Returns:
            True if directory was deleted, False otherwise
        """
        try:
            if dir_path.exists():
                import shutil
                shutil.rmtree(dir_path, ignore_errors=True)
                return True
        except Exception as e:
            logger.warning(f"Failed to delete temp directory {dir_path}: {e}")
        return False

    def get_stats(self) -> dict:
        """
        Get statistics about tracked temp files.

        Returns:
            Dictionary with file/dir counts and total size
        """
        total_size = 0
        existing_files = 0

        for file_path in self.tracked_files:
            if file_path.exists():
                existing_files += 1
                try:
                    total_size += file_path.stat().st_size
                except:
                    pass

        return {
            'tracked_files': len(self.tracked_files),
            'existing_files': existing_files,
            'tracked_dirs': len(self.tracked_dirs),
            'total_size_bytes': total_size,
            'total_size_mb': total_size / (1024 * 1024)
        }

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - cleanup all temp files."""
        self.cleanup()
        return False

    def __del__(self):
        """Destructor - attempt cleanup."""
        try:
            self.cleanup()
        except:
            pass  # Ignore errors during destruction
