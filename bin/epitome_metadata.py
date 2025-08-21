#!/usr/bin/env python3

# epitome_metadata.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import json
import logging
import gzip
import os
import sys
import re
from datetime import datetime
from collections import OrderedDict
from typing import List, Dict, Any, Tuple, Optional

from epitome_utils import sanitize_filename, detect_and_read, normalize_keys, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

STRING_KEYS = ['accession', 'segment', 'tax_id']
YEAR_KEYS    = ['collection_date']
MISSING_DATA = ["na", "n/a", "null", "none", "nan"]

# -------------------------------
#  FUNCTIONS
# -------------------------------

def _extract_year(rec: Dict[str, Any]) -> None:
    """Extract a 4-digit year from YEAR_KEYS fields and replace the value with an int if found.

    Args:
        rec: Metadata record (dict) to modify in-place.

    Returns:
        None
    """
    year_re = re.compile(r"^\d{4}")
    for k in YEAR_KEYS:
        v = rec.get(k, None)
        if v:
            m = year_re.search(str(v))
            if m:
                rec[k] = int(m.group())


def _coerce_to_str(rec: Dict[str, Any]) -> None:
    """Convert certain keys in a record to strings.

    Args:
        rec: Metadata record (dict) to modify in-place.

    Returns:
        None
    """
    for k in STRING_KEYS:
        if k in rec and rec[k] is not None and rec[k] != "":
            rec[k] = str(rec[k])


def _update_missing(rec: Dict[str, Any]) -> None:
    """Replace certain missing value markers with None.

    Args:
        rec: Metadata record (dict) to modify in-place.

    Returns:
        None
    """
    for k, v in rec.items():
        if isinstance(v, str) and v.lower() in MISSING_DATA:
            rec[k] = None


def merge_records(
    inputs: List[Tuple[str, List[Dict[str, Any]]]]
) -> List[Dict[str, Any]]:
    """Merge multiple metadata sources into a single unified list.

    Args:
        inputs: List of tuples (source_path, list_of_records).

    Returns:
        List of merged metadata records.
    """
    acc_key = "accession"
    tax_key = "taxon"
    seg_key = "segment"

    index: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    all_keys = set([acc_key, tax_key])

    for src_path, recs in inputs:
        for rec in recs:

            rec = normalize_keys(rec)

            _coerce_to_str(rec)
            _extract_year(rec)
            _update_missing(rec)

            acc = rec.get(acc_key)
            tax = rec.get(tax_key)
            seg = rec.get(seg_key)

            if acc in (None, ""):
                LOGGER.error("Error: record in %s missing required '%s'.", src_path, acc_key)
                sys.exit(2)
            if tax in (None, ""):
                LOGGER.error("Error: record for accession '%s' in %s missing required '%s'.", acc, src_path, tax_key)
                sys.exit(2)

            if seg in (None, ""):
                rec[seg_key] = "wg"

            if acc not in index:
                index[acc] = {acc_key: acc, tax_key: tax, seg_key: rec[seg_key]}
            target = index[acc]

            for k, v in rec.items():
                if k == acc_key:
                    continue
                is_empty = (v is None) or (isinstance(v, str) and v == "") or (v == [])
                if not is_empty:
                    if (k not in target) or (target[k] in (None, "", [])):
                        target[k] = v
                all_keys.add(k)

    for rec in index.values():
        for k in all_keys:
            if k not in rec:
                rec[k] = None

    return list(index.values())


def drop_na(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove keys that are None for all records.

    Args:
        data: List of metadata records.

    Returns:
        List of records with all-None keys removed.
    """
    n_rec = len(data)
    na_counts: Dict[str, int] = {}
    rm_set = set()

    for rec in data:
        for k, v in rec.items():
            if v is None:
                na_counts[k] = na_counts.get(k, 0) + 1
                if na_counts[k] == n_rec:
                    rm_set.add(k)
    for rec in data:
        for k in rm_set:
            if k in rec:
                rec.pop(k)
    return data


def save_output(
    records: List[Dict[str, Any]], taxon: str, segment: str, outdir: str
) -> str:
    """Write merged metadata to gzipped JSONL.

    Args:
        records: List of merged metadata records.
        taxon: Taxon name.
        segment: Segment name.
        outdir: Output directory.

    Returns:
        Path to the written file.
    """
    safe_taxon = sanitize_filename(taxon)
    safe_segment = sanitize_filename(segment)
    prefix = os.path.join(outdir, f"{safe_taxon}-{safe_segment}")
    filename = prefix + ".metadata.jsonl.gz"

    with gzip.open(filename, "wt", encoding="utf-8") as f:
        for row in records:
            if isinstance(row, dict):
                f.write(json.dumps(row, sort_keys=True) + "\n")

    LOGGER.info("Wrote %d records -> %s", len(records), filename)
    return filename

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """
    CLI entry point for merging JSON/CSV metadata into JSONL.GZ.
    """
    version = "1.0"

    parser = argparse.ArgumentParser(
        description="Merge JSON/CSV by 'accession' into flat JSONL.GZ metadata."
    )
    parser.add_argument("inputs", nargs="+", help="Input files (.json or .csv).")
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--segment", default="null", help="Segment name.")
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    os.makedirs(args.outdir, exist_ok=True)

    # Load all inputs via detector
    loaded: List[Tuple[str, List[Dict[str, Any]]]] = []
    for p in args.inputs:
        try:
            loaded.append(detect_and_read(p))
        except Exception as e:
            LOGGER.exception("Error reading %s: %s", p, e)
            sys.exit(2)

    # Merge and write
    merged = merge_records(loaded)
    no_na = drop_na(merged)
    save_output(no_na, args.taxon, args.segment, args.outdir)


if __name__ == "__main__":
    main()
