# mspprojintegritykeeper.py
"""
Project integrity helper (stubs).
Later we will:
- compute & store file hashes per referenced path
- re-check current files vs stored hashes
- help relink / bulk change project roots
"""

from __future__ import annotations
import os, json, hashlib
from typing import Dict, Any, Iterable, Optional

# ---- core planned APIs (stubs) ----

def writehashtojson(json_path: str, debug: bool=False) -> Optional[str]:
    """
    Read the current project JSON (e.g., .exp.json), compute a hash for each
    referenced file, and write them back into the JSON under a new section (e.g., 'hashes').
    Returns the path written, or None on failure.
    """
    return None  # TODO

def file_hashcheck(hashfromjson: str, currentfile: str, debug: bool=False) -> bool:
    """
    Compare the stored hash (hashfromjson) and the current file's hash.
    Returns True if match, False otherwise.
    """
    return False  # TODO

def relocatemissingfile(lost_ref_entry: Dict[str, Any], debug: bool=False) -> Optional[str]:
    """
    Given a missing reference entry (e.g., {'key': 'csv', 'old_path': '...'}),
    let the user pick a replacement path (tk file dialog) and return it.
    """
    return None  # TODO

# ---- utility ideas (also stubs) ----

def compute_hash(path: str, algo: str="sha256", chunk: int=1024*1024) -> Optional[str]:
    """Stream a file and return its hex digest (or None if path is invalid)."""
    return None  # TODO

def extract_paths_from_exp(exp_json: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    """
    Return {sample_name: {key: absolute_path, ...}, ...} from an exp payload.
    Useful for integrity scans.
    """
    return {}  # TODO

def check_integrity(exp_json_path: str, debug: bool=False) -> Dict[str, Any]:
    """
    High-level: load exp, compute/compare hashes, detect missing/modified files.
    Returns a report dict for the GUI to display.
    """
    return {}  # TODO

def bulk_rewrite_root(exp_json_path: str, old_root: str, new_root: str, dry_run: bool=True) -> Dict[str, Any]:
    """
    Replace a common prefix (old_root) with (new_root) for all relative/absolute paths in the exp.
    If dry_run, return a preview of changes; else write back to disk.
    """
    return {}  # TODO