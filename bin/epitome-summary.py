#!/usr/bin/env python3

# epitome-summary.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import numpy as np
import csv
import datetime
import argparse
import os
import json
from collections import defaultdict
import statistics
import logging

# ----------------------------
# File Loading Utilities
# ----------------------------

def load_json_files(files):
    """Load and merge multiple JSON files."""
    data = {}
    for file in files:
        with open(file) as f:
            data |= json.load(f)
    return {
        taxon: {
            segment: {str(id): d3 for id, d3 in segments.items()}
            for segment, segments in taxon_data.items()
        } for taxon, taxon_data in data.items()
    }

def load_csv_files(files):
    """Load metadata from CSV files into a nested dictionary."""
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for file in files:
        with open(file, newline='') as f:
            reader = csv.reader(f)
            header = next(reader)
            for row in reader:
                row_data = {header[i]: row[i] for i in range(len(header) - 1)}
                taxon = row_data['taxon']
                segment = row_data.get('segment', 'wg')
                accession = row_data['accession']
                data[taxon][segment][accession] = {
                    k: v for k, v in row_data.items()
                    if k not in ['taxon', 'segment', 'accession']
                }
    return data

# ----------------------------
# Data Merging and Conversion
# ----------------------------

def merge_data(target, source):
    """Merge nested data structures into a target dictionary."""
    for key, value in source.items():
        if isinstance(value, dict):
            for key2, value2 in value.items():
                target.setdefault(key, {}).setdefault(key2, set()).add(value2)
        elif isinstance(value, list):
            for item in value:
                target.setdefault(key, set()).add(item)
        else:
            target.setdefault(key, set()).add(value)

def convert_sets(obj):
    """Recursively convert sets to lists for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_sets(v) for k, v in obj.items()}
    elif isinstance(obj, set):
        return list(obj)
    elif isinstance(obj, list):
        return [convert_sets(item) for item in obj]
    else:
        return obj

def get_stats(data):
    sample = list(data)[0]
    if isinstance(sample, float) and not sample.is_integer():
        digits = 3
    else:
        digits = 0
    mean = round(statistics.mean(data), digits)
    stdev = round(statistics.stdev(data), digits) if len(data) > 2 else 'null'

    return {
        'min': round(min(data), digits),
        'max': round(max(data), digits),
        'mean': mean,
        'stdev': stdev if isinstance(stdev, (int, float)) else stdev
    }


# ----------------------------
# Data Processing
# ----------------------------

def summarize_clusters(data_condensed, cluster_lookup, data_qc, data_meta):
    """Summarize and merge QC, cluster, and metadata information."""
    res = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for taxon, segments in data_condensed.items():
        for segment, clusters in segments.items():
            for clust, info in clusters.items():
                clust_key = info['condensed'] or clust
                for seq in cluster_lookup[taxon][segment][clust]:
                    cluster_data = res[taxon][segment].setdefault(clust_key, {})
                    cluster_data.setdefault('seq', set()).add(seq)
                    merge_data(cluster_data, data_qc[taxon][segment][seq])
                    if clust_key == clust:
                        merge_data(cluster_data, data_condensed[taxon][segment][clust])

    for taxon, segments in res.items():
        for segment, clusters in segments.items():
            for clust_key, cluster_data in clusters.items():
                for acc in cluster_data.get('accessions', []):
                    merge_data(cluster_data, data_meta[taxon][segment][acc])
                cluster_data['n_raw'] = len(cluster_data.get('accessions', []))
                cluster_data['n_qc'] = len(cluster_data.get('seq', []))
                finalize_cluster_stats(cluster_data)
    
    return res

def finalize_cluster_stats(cluster_data):
    """Calculate statistics or simplify values for each field in a cluster."""
    for key, value in list(cluster_data.items()):
        if isinstance(value, set):
            cluster_data[key] = try_stats(value)
        elif isinstance(value, dict):
            for subkey, subval in list(value.items()):
                if isinstance(subval, set) and subkey not in ['filter']:
                    value[subkey] = try_stats(subval)

def try_stats(value):
    """Attempt to return statistics for a set, otherwise return as list or single value."""
    try:
        return get_stats(value)
    except Exception:
        if len(value) > 1:
            return list(value)
        else:
            return list(value)[0]

def remove_keys(data, keys_to_remove):
    """Recursively remove specific keys from nested dictionaries."""
    if isinstance(data, dict):
        for key in list(data.keys()):
            if key in keys_to_remove:
                del data[key]
            else:
                remove_keys(data[key], keys_to_remove)
    elif isinstance(data, list):
        for item in data:
            remove_keys(item, keys_to_remove)

def summary_to_csv_rows(summary_data):
    """
    Converts nested summary data to a flat list of CSV rows.
    Skips the 'accessions' key and formats 'length' and 'ambRatio' as min-max (mean±stdev).
    """
    rows = []
    for taxon, segments in summary_data.items():
        for segment, clusters in segments.items():
            for cluster_id, data in clusters.items():
                row = {
                    "taxon": taxon,
                    "segment": segment,
                    "cluster": cluster_id,
                    "ref": f'{taxon}-{segment}-{cluster_id}'
                }

                for key, value in data.items():
                    if key == 'accessions':
                        continue
                    if key in ('length', 'ambRatio') and isinstance(value, dict):
                        row[key] = format_stats(value)
                    elif isinstance(value, dict):
                        # Flatten nested dicts (e.g., quality fields)
                        for subkey, subval in value.items():
                            colname = f"{key}.{subkey}"
                            row[colname] = format_stats(subval) if isinstance(subval, dict) else subval
                    else:
                        row[key] = value
                rows.append(row)
    return rows

def format_stats(stats_dict):
    """
    Format statistics dict as: min-max (mean±stdev)
    """
    try:
        return f"{stats_dict['value']['min']}-{stats_dict['value']['max']}; {stats_dict['value']['mean']}±{stats_dict['value']['stdev']}"
    except (KeyError, TypeError):
        return str(stats_dict['value'])

def write_summary_csv(rows, filename):
    """Write the flattened summary data to a CSV file with fixed first columns."""
    if not rows:
        print("No data to write.")
        return

    logging.info(f'Writing simple metrics to {filename}')
    # Collect all field names from data
    all_fields = set(k for row in rows for k in row.keys())

    # Define fixed order for first columns
    fixed_order = ['taxon', 'segment', 'cluster','ref','length','ambRatio','n_raw','n_qc','neighbor_ani','neighbor']

    # Remaining columns, sorted alphabetically and excluding the fixed ones
    remaining = sorted(f for f in all_fields if f not in fixed_order)

    # Final header order
    fieldnames = fixed_order + remaining

    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ----------------------------
# Output
# ----------------------------

def write_json(data, filename):
    logging.info(f'Writing full metrics to {filename}')
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

# ----------------------------
# Main Function
# ----------------------------

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main():
    version = "1.0"
    parser = argparse.ArgumentParser(description="Input QC script for DNA sequences (no pandas)")
    parser.add_argument("--qc", nargs='+', help="Path to JSON file(s) from eptiome-qc.py")
    parser.add_argument("--clusters", nargs='+', help="Path to JSON file(s) from eptiome-cluster.py")
    parser.add_argument("--condensed", nargs='+', help="Path to JSON file(s) from eptiome-condense.py")
    parser.add_argument("--meta", nargs='+', help="CSV files with sample metadata")
    parser.add_argument("--outdir", default='.', help="Output directory")
    parser.add_argument('--version', action='version', version=version)
    args = parser.parse_args()

    print(f"""
    epitome-summary.py v{version}
    Written by Jared Johnson
    Jared.johnson@doh.wa.gov
    """)

    data_qc = load_json_files(args.qc)
    data_cluster = load_json_files(args.clusters)
    data_condensed = load_json_files(args.condensed)
    data_meta = load_csv_files(args.meta)

    summarized_data = summarize_clusters(data_condensed, data_cluster, data_qc, data_meta)
    remove_keys(summarized_data, keys_to_remove={'seq', 'illegalBases', 'condensed'})
    serializable_res = convert_sets(summarized_data)
    write_json(serializable_res, f'{args.outdir}/summary.json')
    csv_rows = summary_to_csv_rows(serializable_res)
    write_summary_csv(csv_rows, f'{args.outdir}/summary.csv')

if __name__ == "__main__":
    main()
