# pathcanon.py
# GPTmade for solving path inconsistent issue
from __future__ import annotations
from pathlib import Path, PureWindowsPath, PurePosixPath
import os

def to_posix_str(p: str | os.PathLike) -> str:
    """
    Canonicalize for STORAGE: forward slashes only, collapse // and mixed slashes.
    Works even if the file doesn't exist.
    """
    s = str(p)
    # Use Windows parser to handle drive letters/UNC, then express as POSIX
    try:
        return PureWindowsPath(s).as_posix()
    except Exception:
        # Fallback: manual replace
        return s.replace("\\", "/").replace("//", "/").replace(":/", ":/")

def to_native_path(p: str | os.PathLike) -> Path:
    """
    For RUNTIME: convert JSON string to the current-OS Path.
    """
    return Path(str(p))

def norm_join(root: str | os.PathLike, *parts: str) -> Path:
    """
    Join using pathlib (never string-concat); returns Path on current OS.
    """
    cur = Path(root)
    for part in parts:
        cur = cur / part
    return cur

def ensure_dir(d: str | os.PathLike) -> Path:
    p = Path(d)
    p.mkdir(parents=True, exist_ok=True)
    return p

def strip_trailing_sep(s: str) -> str:
    return s.rstrip("\\/")

def is_mixed_sep(s: str) -> bool:
    return ("\\" in s) and ("/" in s)