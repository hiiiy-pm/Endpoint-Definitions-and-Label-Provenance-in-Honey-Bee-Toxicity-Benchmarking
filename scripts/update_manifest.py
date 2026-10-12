"""Refresh the archive checksums after intentional edits and verification."""
from pathlib import Path
import runpy

ROOT=Path(__file__).resolve().parents[1]
runpy.run_path(str(ROOT/'scripts/make_file_manifest.py'), run_name='__main__')
