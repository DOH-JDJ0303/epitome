#!/usr/bin/env python3

# epitome_qc.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import os
import re
import statistics
import screed
from collections import defaultdict
import json
import logging
import gzip
import textwrap
from typing import List, Dict, Any, Tuple

from epitome_utils import sanitize_filename, read_csv, normalize_keys, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
#  FUNCTIONS
# -------------------------------

LEGAL_BASE_RE = re.compile(r"[-ATCGRYSWKMBDHVN]")

def load_seqs(fasta: str) -> List[Dict[str, Any]]:
    """Load sequences from FASTA and return list of dicts.

    Args:
        fasta: Path to input FASTA file.

    Returns:
        List of records, each with keys: accession (str), length (int), sequence (str).
    """
    sequences: List[Dict[str, Any]] = []
    for record in screed.open(fasta):
        accession = record.name.split()[0]
        seq = record.sequence.upper()
        sequences.append({
            "accession": accession,
            "length": len(seq),
            "sequence": seq
        })
    LOGGER.info(f'Loaded {len(sequences)} sequences from {fasta}')
    return sequences


def exclude_seqs(sequences: List[Dict[str, Any]], exclusions_path: str) -> List[Dict[str, Any]]:
    """Exclude sequences by accession from a CSV file.

    Args:
        sequences: Input list of sequence records.
        exclusions_path: CSV file path containing an 'accession' column.

    Returns:
        Filtered list of sequence records.
    """
    initial_count = len(sequences)

    rows = read_csv(exclusions_path)
    # Be tolerant of header case
    rows = [normalize_keys(r, lower_keys=True) for r in rows]
    excluded = {r["accession"] for r in rows if "accession" in r and r["accession"]}

    if not excluded:
        LOGGER.info(f'No valid "accession" values found in {exclusions_path}; skipping exclusion.')
        return sequences

    keep = [s for s in sequences if s["accession"] not in excluded]
    LOGGER.info(f'Excluded {initial_count - len(keep)} sequences using {exclusions_path}')
    return keep


def consolidate_seqs(sequences: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Group identical sequences, assigning shared seq_id and accessions list.

    Args:
        sequences: List of sequence records with 'sequence' and 'length'.

    Returns:
        List of consolidated records with keys: seq_id (int), sequence (str),
        length (int), accessions (List[str]).
    """
    initial_count = len(sequences)
    grouped = defaultdict(list)
    for s in sequences:
        key = (s["sequence"], s["length"])
        grouped[key].append(s["accession"])
    data = [
        {
            "seq_id": i,
            "sequence": key[0],
            "length": key[1],
            "accessions": accessions
        }
        for i, (key, accessions) in enumerate(grouped.items(), start=1)
    ]
    LOGGER.info(f'Consolidated {initial_count} to {len(data)} unique sequences')
    return data


def filter_seqs(data: List[Dict[str, Any]], amb_thresh: float, len_thresh: float) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Apply illegal base, ambiguous base (N), and length filters.

    Args:
        data: Consolidated sequence records.
        amb_thresh: Maximum allowed fraction of 'N' bases.
        len_thresh: Fractional tolerance around median length for pass.

    Returns:
        (Updated records with QC fields, list of FASTA records that passed).
    """
    lengths = [d["length"] for d in data]
    median = statistics.median(lengths) if lengths else 0
    lower = median * (1 - len_thresh)
    upper = median * (1 + len_thresh)

    def test(flag: bool) -> str:
        return "fail" if flag else "pass"

    passing_fasta_records: List[str] = []
    for d in data:
        seq = d["sequence"]
        length = d["length"]
        illegal = re.sub(LEGAL_BASE_RE, "", seq)
        amb_ratio = (seq.count("N") / length) if length else 1.0

        fail_illegal = len(illegal) > 0
        fail_amb = amb_ratio > amb_thresh
        fail_len = (length < lower) or (length > upper)

        d["illegal_bases"] = {"filter": illegal, "value": illegal, "status": test(fail_illegal)}
        d["amb_ratio"] = {"filter": amb_thresh, "value": amb_ratio, "status": test(fail_amb)}
        d["length"] = {"filter": len_thresh, "value": length, "status": test(fail_len)}

        if not (fail_illegal or fail_amb or fail_len):
            passing_fasta_records.append(f">{d['seq_id']}\n{seq}")

    LOGGER.info(f'{len(passing_fasta_records)} sequences passed all filters')
    return data, passing_fasta_records


def save_output(json_data: List[Dict[str, Any]], fasta_data: List[str], taxon: str, segment: str, outdir: str) -> None:
    """Save JSONL and FASTA outputs to compressed files.

    Args:
        json_data: QC result records to write (JSONL.GZ).
        fasta_data: FASTA entries (header+sequence lines) that passed QC.
        taxon: Taxon label for file prefix.
        segment: Segment label for file prefix.
        outdir: Output directory.

    Returns:
        None
    """
    prefix = os.path.join(outdir, f"{sanitize_filename(taxon)}-{sanitize_filename(segment)}")

    # Write JSONL
    json_file = f"{prefix}.qc.jsonl.gz"
    with gzip.open(json_file, 'wt', encoding='utf-8') as f:
        for row in json_data:
            out = {
                k: row.get(k)
                for k in ["seq_id", "illegal_bases", "amb_ratio", "length", "accessions"]
            }
            out.update({"taxon": taxon, "segment": segment})
            f.write(json.dumps(out, sort_keys=True) + "\n")
    LOGGER.info(f"Saved QC JSON to {json_file}")

    # Write FASTA
    fasta_file = f"{prefix}.qc.fa.gz"
    with gzip.open(fasta_file, 'wt', encoding='utf-8') as f:
        f.write("\n".join(fasta_data) + "\n")
    LOGGER.info(f"Saved passing sequences to {fasta_file}")


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point for QC and deduplication of FASTA sequences.

    Parses arguments, loads FASTA, optionally excludes accessions,
    consolidates identical sequences, applies QC filters, and writes
    compressed JSONL and FASTA outputs.

    """
    version = "2.0"

    parser = argparse.ArgumentParser(description="QC and deduplication for FASTA sequences")
    parser.add_argument("--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--taxon", default='null', help="Taxon name.")
    parser.add_argument("--segment", default='null', help="Segment name.")
    parser.add_argument("--amb_threshold", type=float, default=0.02, help="Max N base ratio.")
    parser.add_argument("--len_threshold", type=float, default=0.20, help="Length tolerance fraction.")
    parser.add_argument("--exclusions", help="CSV file with accessions to exclude.")
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    os.makedirs(args.outdir, exist_ok=True)

    seqs = load_seqs(args.fasta)
    if args.exclusions:
        seqs = exclude_seqs(seqs, args.exclusions)

    consolidated = consolidate_seqs(seqs)
    qc_results, fasta_pass = filter_seqs(consolidated, args.amb_threshold, args.len_threshold)
    save_output(qc_results, fasta_pass, args.taxon, args.segment, args.outdir)

    print(f"Raw Count: {len(seqs)}\nFinal Count: {len(fasta_pass)}")


if __name__ == "__main__":
    main()
