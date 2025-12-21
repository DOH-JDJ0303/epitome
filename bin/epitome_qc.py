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
from typing import List, Dict, Any, Tuple, Optional

from epitome_utils import sanitize_filename, read_csv, normalize_keys, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
#  FUNCTIONS
# -------------------------------

LEGAL_BASE_RE = re.compile(r"[-ATCGRYSWKMBDHVN]")


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
    rows = [normalize_keys(r) for r in rows]
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


def filter_seqs(
    data: List[Dict[str, Any]],
    amb_thresh: float,
    z_threshold: float,
    metadata_map: Optional[Dict[str, Dict[str, Any]]] = None
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Apply illegal base, ambiguous base (N), and z-score length filters.

    Args:
        data: Consolidated sequence records.
        amb_thresh: Maximum allowed fraction of 'N' bases.
        z_threshold: Maximum absolute z-score for length (standard deviations from mean).
        metadata_map: Optional accession -> metadata mapping (JSONL).

    Returns:
        (Updated records with QC fields, list of FASTA records that passed).
    """
    # Determine which sequences to use for the length statistics calculation.
    lengths: List[int] = []

    if metadata_map:
        for d in data:
            for acc in d.get("accessions", []):
                m = metadata_map.get(acc)
                if isinstance(m, dict) and m.get("canonical") is True:
                    lengths.append(d["length"])

    if not lengths:
        LOGGER.info("Using all sequences for length filtering (No canonical sequences supplied)")
        lengths = [d["length"] for d in data]
    else:
        LOGGER.info(f"Using canonical sequences for length filtering (n={len(lengths)})")

    if len(lengths) >= 2:
        min_len = min(lengths)
        max_len = max(lengths)
        mean_len = statistics.mean(lengths)
        stdev_len = statistics.stdev(lengths)
        lower = mean_len - (z_threshold * stdev_len)
        upper = mean_len + (z_threshold * stdev_len)
    elif len(lengths) == 1:
        min_len = max_len = mean_len = lengths[0]
        stdev_len = 0.0
        lower = upper = mean_len
    else:
        min_len = max_len = mean_len = stdev_len = lower = upper = 0

    LOGGER.info(
        f"Lengths: min={min_len} max={max_len} mean={mean_len:.2f} stdev={stdev_len:.2f} | "
        f"range=[{lower:.2f}, {upper:.2f}] (±{z_threshold} SD)"
    )

    def test(flag: bool) -> str:
        return "fail" if flag else "pass"

    passing_fasta_records: List[str] = []
    for d in data:
        seq = d["sequence"]
        length = d["length"]
        illegal = re.sub(LEGAL_BASE_RE, "", seq)
        amb_ratio = (seq.count("N") / length) if length else 1.0

        # Calculate z-score for this sequence
        if stdev_len > 0:
            z_score = (length - mean_len) / stdev_len
        else:
            z_score = 0.0

        fail_illegal = len(illegal) > 0
        fail_amb = amb_ratio > amb_thresh
        fail_len = abs(z_score) > z_threshold

        d["illegal_bases"] = {"filter": illegal, "value": illegal, "status": test(fail_illegal)}
        d["amb_ratio"] = {"filter": amb_thresh, "value": amb_ratio, "status": test(fail_amb)}
        d["length"] = {"filter": z_threshold, "value": length, "z_score": z_score, "status": test(fail_len)}

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

    if fasta_data:
        # Write FASTA
        fasta_file = f"{prefix}.qc.fa.gz"
        with gzip.open(fasta_file, 'wt', encoding='utf-8') as f:
            f.write("\n".join(fasta_data) + "\n")
        LOGGER.info(f"Saved passing sequences to {fasta_file}")
    else:
        LOGGER.info(f"All sequences failed QC. Nothing to save.")

def load_input(filepath):
    LOGGER.info(f"Loading file: {filepath}")
    metadata: Dict[str, Dict[str, Any]] = {}
    seqs: list = []
    count = 0
    taxon_segment = set()

    with gzip.open(filepath, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if isinstance(entry, dict) and 'accession' in entry:
                metadata[str(entry['accession'])] = {k: v for k, v in entry.items() if k not in ['accession', 'sequence']}
                seqs.append(
                    {'sequence': entry['sequence'], 
                    'length': len(entry['sequence']), 
                    'accession': entry['accession']
                    })
                count += 1
                if entry.get('taxon') and entry.get('segment'):
                    taxon_segment.add((entry['taxon'], entry['segment']))
    
    if len(taxon_segment) == 1:
        taxon_segment = list(taxon_segment)[0]
    else:
        raise ValueError("One taxon / segment required")

    LOGGER.info(f"Loaded {count} records from {filepath}")
    return metadata, seqs, taxon_segment


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point for QC and deduplication of FASTA sequences.

    Parses arguments, loads FASTA, optionally excludes accessions,
    consolidates identical sequences, applies QC filters, and writes
    compressed JSONL and FASTA outputs.

    """
    version = "2.1"

    parser = argparse.ArgumentParser(description="QC and deduplication for FASTA sequences")
    parser.add_argument("input", help="Optional JSONL with canonical sequence metadata")
    parser.add_argument("--amb_threshold", type=float, default=0.02, help="Max N base ratio.")
    parser.add_argument("--z_threshold", type=float, default=2.58, help="Maximum z-score (standard deviations from mean) for length filtering.")
    parser.add_argument("--exclusions", help="CSV file with accessions to exclude.")
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    os.makedirs(args.outdir, exist_ok=True)

    metadata_map, seqs, (taxon, segment) = load_input(args.input)

    if args.exclusions:
        seqs = exclude_seqs(seqs, args.exclusions)

    consolidated = consolidate_seqs(seqs)
    qc_results, fasta_pass = filter_seqs(consolidated, args.amb_threshold, args.z_threshold, metadata_map=metadata_map)
    save_output(qc_results, fasta_pass, taxon, segment, args.outdir)

    print(f"Raw Count: {len(seqs)}\nFinal Count: {len(fasta_pass)}")


if __name__ == "__main__":
    main()