#!/usr/bin/env python3
# epitome_utils.py

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import re
import statistics
from collections import defaultdict
from datetime import datetime
from typing import Any, Dict, Iterable, List, Mapping, Tuple
import logging
import sys

def logging_config(log_level='INFO'):
    # Get script name
    script_name = os.path.basename(sys.argv[0]).replace('.py', '')
    LOGGER = logging.getLogger(script_name)

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    return LOGGER


def sanitize_filename(s: str) -> str:
    """Remove/replace characters unsafe for filesystem paths."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', s)


def normalize_keys(d: Mapping[str, Any]) -> Dict[str, Any]:
    """Lowercase snake format"""

    return {(re.sub(r'([a-z])([A-Z])', r'\1_\2', str(k)).lower()): v for k, v in d.items()}


# ---------- Lightweight Symmetric Cache ---------- #

class DistanceCache:
    """Symmetric key cache: cache[(scope, a, b)] == cache[(scope, b, a)]."""
    def __init__(self):
        self._d: Dict[Tuple[str, Any, Any], Any] = {}

    def _key(self, scope: str, a: Any, b: Any) -> Tuple[str, Any, Any]:
        return (scope, a, b) if a <= b else (scope, b, a)

    def get(self, scope: str, a: Any, b: Any) -> Any:
        return self._d.get(self._key(scope, a, b))

    def set(self, scope: str, a: Any, b: Any, val: Any) -> None:
        self._d[self._key(scope, a, b)] = val


# ---------- JSON / CSV Readers ---------- #

def _open_maybe_gzip(path: str, mode: str):
    """Open normally or via gzip based on file suffix."""
    if path.endswith(".gz"):
        if "t" in mode:
            return gzip.open(path, mode, encoding="utf-8")
        return gzip.open(path, mode)
    # Plain file
    if "t" in mode:
        return open(path, mode, encoding="utf-8", newline="")
    return open(path, mode)


def read_json(path: str) -> List[Dict[str, Any]]:
    """Load a JSON or JSON.gz file whose root is a list of objects."""
    with _open_maybe_gzip(path, "rt") as f:
        obj = json.load(f)

    if not isinstance(obj, list):
        raise ValueError(f"{path}: JSON root must be a list of objects")

    for i, rec in enumerate(obj):
        if not isinstance(rec, dict):
            raise ValueError(f"{path}: element {i} is not an object")

    return obj


def read_csv(path: str) -> List[Dict[str, Any]]:
    """Load a CSV or CSV.gz into list of dicts; trims strings & converts empty to None."""
    recs: List[Dict[str, Any]] = []

    with _open_maybe_gzip(path, "rt") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header row")

        for row in reader:
            clean: Dict[str, Any] = {}
            for k, v in row.items():
                k = k.strip() if isinstance(k, str) else k
                if isinstance(v, str):
                    v = v.strip()
                    v = v if v != "" else None
                clean[k] = v
            recs.append(clean)

    return recs


def detect_and_read(path: str) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Read JSON/CSV (optionally gzipped), with JSON→CSV fallback.
    Extensions handled:
      - .json, .json.gz
      - .csv, .csv.gz
    """
    base, ext = os.path.splitext(path)
    ext = ext.lower()

    # Handle gz suffix: ext = ".gz" and base endswith .json/.csv
    if ext == ".gz":
        base2, ext2 = os.path.splitext(base)
        ext2 = ext2.lower()
        if ext2 == ".json":
            return path, read_json(path)
        if ext2 == ".csv":
            return path, read_csv(path)
        # Otherwise fallback below

    # Normal non-gzip cases
    if ext == ".json":
        return path, read_json(path)
    if ext == ".csv":
        return path, read_csv(path)

    # Fallback: try JSON then CSV
    try:
        return path, read_json(path)
    except Exception:
        return path, read_csv(path)

# ---------- JSONL Helpers ---------- #

def read_jsonl(path: str) -> List[Dict[str, Any]]:
    """Load JSON Lines (one JSON object per line)."""
    out: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    logging.info(f"Loaded JSONL file: {path}")
    return out


def load_json(path: str) -> List[Dict[str, Any]]:
    """
    For compatibility with prior usage: wrap a single JSON object into a list.
    If the file already contains a list, it becomes [that_list].
    """
    with open(path, "r", encoding="utf-8") as f:
        out = [json.load(f)]
    logging.info(f"Loaded JSON file: {path}")
    return out


def load_json_or_jsonl(path: str) -> List[Dict[str, Any]]:
    """Try JSON first, then JSONL."""
    try:
        return load_json(path)
    except json.JSONDecodeError:
        return read_jsonl(path)


def load_jsonl_gz_by_key(filepath: str, key: str) -> Dict[str, Dict[str, Any]]:
    """
    Load a .jsonl.gz into a dict keyed by `key`, removing `key` from each value dict.
    """
    logging.info(f"Loading file: {filepath} using key: {key}")
    data: Dict[str, Dict[str, Any]] = {}
    count = 0
    with gzip.open(filepath, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if isinstance(entry, dict) and key in entry:
                data[str(entry[key])] = {k: v for k, v in entry.items() if k != key}
                count += 1
    logging.info(f"Loaded {count} records from {filepath}")
    return data

def is_date_string(v: Any) -> bool:
    """True if string parses via datetime.fromisoformat()."""
    try:
        if isinstance(v, str):
            datetime.fromisoformat(v)
            return True
    except Exception:
        return False
    return False


__all__ = [
    "sanitize_filename",
    "normalize_keys",
    "DistanceCache",
    "read_json",
    "read_csv",
    "detect_and_read",
    "read_jsonl",
    "load_json",
    "load_json_or_jsonl",
    "load_jsonl_gz_by_key",
    "is_date_string"
]
