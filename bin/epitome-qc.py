#!/usr/bin/env python3

# epitome-qc.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import os
import re
import csv
import statistics
import screed
from collections import defaultdict
import json
import datetime
import textwrap

# ---- FUNCTIONS ----
def loadSeqs(fasta):
    """
    Load sequences from a FASTA file using screed and return as a list of dictionaries.
    Each entry contains accession, sequence length, and sequence string.
    """
    sequences = []
    for record in screed.open(fasta):
        accession = record.name.split()[0]  # First token in header (e.g., accession ID)
        seq = record.sequence.upper()       # Convert sequence to uppercase
        sequences.append({
            "accession": accession,
            "length": len(seq),
            "seqString": seq
        })
    print(f'{datetime.datetime.now()}: Loaded {len(sequences)} sequences from {fasta}')
    return sequences


def excludeSeqs(sequences, exclusions):
    """
    Exclude sequences whose accession appears in the exclusions CSV.
    """
    start_len = len(sequences)  # Corrected typo from 'star_len'
    excluded_accessions = set()
    with open(exclusions, newline='') as excl_file:
        reader = csv.DictReader(excl_file)
        if "accession" in reader.fieldnames:
            excluded_accessions = {row["accession"] for row in reader}  # Load exclusions

    # Filter out excluded sequences
    sequences = [s for s in sequences if s["accession"] not in excluded_accessions]
    print(f'{datetime.datetime.now()}: Excluded {start_len - len(sequences)} sequences')
    return sequences


def consolidateSeqs(sequences):
    """
    Collapse duplicate sequences (identical string and length) into one entry,
    aggregating accessions for each unique sequence.
    """
    start_len = len(sequences)
    seq_dict = {}
    for s in sequences:
        key = (s["seqString"], s["length"])  # Use sequence and length as unique key
        if key not in seq_dict:
            seq_dict[key] = []
        seq_dict[key].append(s["accession"])

    # Convert to a list of dictionaries with grouped accessions
    data = []
    for i, ((seqString, length), accessions) in enumerate(seq_dict.items(), start=1):
        data.append({
            "seq": i,  # Unique integer ID
            "seqString": seqString,
            "length": length,
            "accessions": accessions
        })

    print(f'{datetime.datetime.now()}: Consolidated from {start_len} to {len(data)} sequences')
    return data

def filterSeqs(data, amb_thresh, len_thresh):
    """
    Apply sequence quality filters:
      - Illegal bases (non-IUPAC characters)
      - Ambiguous base ratio (N content)
      - Length deviation from the median
    Returns:
      - Updated data with filter metadata
      - List of sequence IDs that passed all filters
    """
    print(f'{datetime.datetime.now()}: Applying filters')
    legal_base_regex = re.compile(r"[-ATCGRYSWKMBDHVN]")  # Allowed base characters
    lengths = [entry["length"] for entry in data]
    median_len = statistics.median(lengths)
    upper_len = median_len * (1 + len_thresh)
    lower_len = median_len * (1 - len_thresh)

    def filter_test(flag):
        return "fail" if flag else "pass"

    passing = []
    for entry in data:
        seq = entry["seqString"]
        length = entry["length"]
        illegal_bases = re.sub(legal_base_regex, "", seq)  # Remove allowed bases
        amb_ratio = seq.count("N") / length
        filter_illegal = len(illegal_bases) > 0
        filter_amb = amb_ratio > amb_thresh
        filter_len = length >= upper_len or length <= lower_len

        # Annotate entry with detailed filter metadata
        entry.update({
            "illegalBases": {
                'filter': illegal_bases,
                'value': illegal_bases,
                'status': filter_test(filter_illegal)
            },
            "ambRatio": {
                'filter': amb_thresh,
                'value': amb_ratio,
                'status': filter_test(filter_amb)
            },
            "length": {
                'filter': len_thresh,
                'value': length,
                'status': filter_test(filter_len)
            }
        })

        if not filter_illegal and not filter_amb and not filter_len:
            passing.append(entry['seq'])

    return data, passing


def dict2Fasta(data, inclusions, filename):
    """
    Write only the sequences that passed all filters to a FASTA file.
    Uses the unique numeric `seq` ID as the FASTA header.
    """
    print(f'{datetime.datetime.now()}: Writing passing sequences to {filename}')
    with open(filename, "w") as fasta_out:
        for entry in data:
            if 'seq' in entry and 'seqString' in entry:
                if entry['seq'] in inclusions:
                    fasta_out.write(f">{entry['seq']}\n{entry['seqString']}\n")


def dict2Json(data, taxon, segment, fieldnames, filename):
    """
    Restructure flat list of data entries into nested JSON:
      taxon -> segment -> sequence_id -> metadata fields
    """
    print(f'{datetime.datetime.now()}: Writing metrics to {filename}')
    nested = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for entry in data:
        filtered_entry = {key: value for key, value in entry.items()
                          if key in fieldnames and key != 'seq'}
        nested[taxon][segment][entry['seq']] = filtered_entry

    with open(filename, 'w') as f:
        json.dump(nested, f, indent=4)


def dict2Csv(data, taxon, segment, fieldnames, filename):
    """
    Write filtered metadata to a flat CSV using only selected fieldnames.
    """
    with open(filename, "w", newline='') as csv_out:
        writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
        writer.writeheader()
        for entry in data:
            filtered_entry = {key: value for key, value in entry.items() if key in fieldnames}
            writer.writerow(filtered_entry)

def main():
    version = "1.0"

    parser = argparse.ArgumentParser(description="Input QC script for DNA sequences (no pandas)")
    parser.add_argument("--fasta", help="Path to the FASTA file")
    parser.add_argument("--taxon", default='null', help="Taxon name or 'version'")
    parser.add_argument("--segment", default='null', help="Segment name")
    parser.add_argument("--amb_threshold", type=float, default=0.02, help="Threshold for ambiguous base (N) ratio")
    parser.add_argument("--len_threshold", type=float, default=0.20, help="Threshold for sequence length deviation")
    parser.add_argument("--exclusions", help="Path to CSV file with accessions to exclude")
    parser.add_argument("--outdir", default='./', type=str, help="Output directory")
    parser.add_argument('--version', action='version', version=version)

    args = parser.parse_args()

    start = f"""
    epitome-qc.py v{version}
    Written by Jared Johnson
    Jared.johnson@doh.wa.gov
    """
    print(textwrap.dedent(start))

    os.makedirs(f'{args.outdir}/', exist_ok=True)
    file_basename = f"{args.taxon.replace(' ', '_')}-{args.segment}"

    seqs = loadSeqs(args.fasta)
    if args.exclusions:
        seqs = excludeSeqs(seqs, args.exclusions)
    res = consolidateSeqs(seqs)
    res, passing = filterSeqs(res, args.amb_threshold, args.len_threshold)

    # ---- SAVE TO FILE ----
    fields = ['seq','illegalBases','ambRatio', 'length', 'accessions']
    dict2Fasta(res, passing, f'{args.outdir}/{file_basename}.qc.fa')
    dict2Json(res, args.taxon, args.segment, fields, f'{args.outdir}/{file_basename}.qc.json')
    # dict2Csv(res, args.taxon, args.segment, fields, f'{args.outdir}/{file_basename}.qc.csv')

    end = f"""
    Raw Count: {len(seqs)}
    Final Count: {len(passing)}
    """
    print(textwrap.dedent(end))

if __name__ == "__main__":
    main()