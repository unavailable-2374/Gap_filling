#!/usr/bin/env python3
"""
Temporary File Manager

Provides reliable management and cleanup of temporary files during gap filling.
"""

import os
import tempfile
import logging
import atexit
import shutil
from pathlib import Path
from contextlib import contextmanager
from typing import Optional, Set, Union

logger = logging.getLogger(__name__)


class TempFileManager:
    """
    Manages temporary files with automatic cleanup.

    Features:
    - Context manager support
    - Tracks all created temp files
    - Cleanup on exit (atexit handler)
    - Configurable temp directory
    """

    def __init__(self,
                 work_dir: Optional[Union[str, Path]] = None,
                 base_dir: Optional[Union[str, Path]] = None,
                 prefix: str = "gapfill_",
                 auto_cleanup: bool = True):
        effective_dir = work_dir if work_dir else base_dir
        self.base_dir = Path(effective_dir) if effective_dir else None
        self.prefix = prefix
        self.tracked_files: Set[Path] = set()
        self.tracked_dirs: Set[Path] = set()

        if self.base_dir:
            self.base_dir.mkdir(parents=True, exist_ok=True)

        if auto_cleanup:
            atexit.register(self.cleanup)

    @contextmanager
    def temp_file(self, suffix: str = "", prefix: Optional[str] = None, delete: bool = False):
        file_prefix = prefix if prefix else self.prefix

        fd, temp_path = tempfile.mkstemp(
            suffix=suffix,
            prefix=file_prefix,
            dir=self.base_dir
        )
        os.close(fd)

        temp_path = Path(temp_path)
        self.tracked_files.add(temp_path)

        try:
            yield temp_path
        finally:
            if delete:
                self._safe_delete_file(temp_path)
                self.tracked_files.discard(temp_path)

    @contextmanager
    def temp_dir(self, suffix: str = "", prefix: Optional[str] = None, delete: bool = False):
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

    def create_temp_file(self, suffix: str = "", prefix: Optional[str] = None) -> Path:
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

    def cleanup(self):
        for file_path in list(self.tracked_files):
            self._safe_delete_file(file_path)
            self.tracked_files.discard(file_path)

        for dir_path in list(self.tracked_dirs):
            self._safe_delete_dir(dir_path)
            self.tracked_dirs.discard(dir_path)

    def _safe_delete_file(self, file_path: Path) -> bool:
        try:
            if file_path.exists():
                file_path.unlink()
                return True
        except Exception as e:
            logger.warning(f"Failed to delete temp file {file_path}: {e}")
        return False

    def _safe_delete_dir(self, dir_path: Path) -> bool:
        try:
            if dir_path.exists():
                shutil.rmtree(dir_path, ignore_errors=True)
                return True
        except Exception as e:
            logger.warning(f"Failed to delete temp directory {dir_path}: {e}")
        return False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False
