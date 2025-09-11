#!/usr/bin/env python3

# epitome_summary.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import os
import json
import gzip
import logging
import statistics
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Mapping
from datetime import datetime

from epitome_utils import sanitize_filename, load_jsonl_gz_by_key, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()
now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

# -------------------------------
#  FUNCTIONS
# -------------------------------

def summarize_metadata(metadata: Mapping[str, Iterable[Any]]) -> Dict[str, Any]:
    """Summarize a mapping of metadata key -> iterable of values.

    For each field:
      - Numeric floats: min, max, mean, sd
      - Integers: min, max
      - Other types: unique values preserving first-seen order

    Args:
        metadata: Mapping from field name to iterable of values.

    Returns:
        Dictionary mapping field name to summary object.
    """

    def summarize_value_list(values_in: Iterable[Any]) -> Any:
        values = [v for v in values_in if v not in (None, "")]
        if not values:
            return None

        # Numeric float summary
        if all(isinstance(v, float) for v in values):
            return {
                "min": min(values),
                "max": max(values),
                "mean": statistics.fmean(values),
                "sd": statistics.stdev(values) if len(values) >= 2 else None,
            }

        # Integer summary
        if all(isinstance(v, int) for v in values):
            return {
                "min": min(values),
                "max": max(values),
            }

        # Categorical / mixed: unique values preserving type+value uniqueness
        seen = set()
        uniq = []
        for v in values:
            key = (type(v), v)
            if key not in seen:
                seen.add(key)
                uniq.append(v)
        return uniq

    return {
        k: summary
        for k, vs in metadata.items()
        if (summary := summarize_value_list(vs)) is not None
    }


def save_output(
    json_data: List[Dict[str, Any]],
    taxon: str,
    segment: str,
    method: str,
    outdir: str
) -> None:
    """Save summary JSONL.GZ output.

    Args:
        json_data: List of summary row dictionaries.
        taxon: Taxon name for file prefix.
        segment: Segment name for file prefix.
        method: Method label for file prefix.
        outdir: Output directory.

    Returns:
        None
    """
    safe_taxon = sanitize_filename(taxon)
    safe_segment = sanitize_filename(segment)
    prefix = os.path.join(outdir, f"{safe_taxon}-{safe_segment}.{method}")
    filename = prefix + '.jsonl.gz'
    with gzip.open(filename, 'wt', encoding='utf-8') as f:
        for row in json_data:
            if isinstance(row, dict):
                json.dump(row, f)
                f.write('\n')
    LOGGER.info(f"Saved summary to {filename}")

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point for summarizing QC, cluster, and metadata JSONL files.

    Loads per-sequence QC results, cluster assignments, and metadata records;
    aggregates per-cluster accession lists; merges metadata values across
    accessions; produces field-level summaries; and writes JSONL.GZ output.

    """
    version = "2.0"

    parser = argparse.ArgumentParser(description="Summarize QC + cluster + metadata.")
    parser.add_argument("--qc", required=True, help="Path to JSON file from epitome-qc.py.")
    parser.add_argument("--clusters", required=True, help="Path to JSON file from epitome-cluster.py or epitome-condense.py.")
    parser.add_argument("--meta", required=True, help="Path to JSON file with metadata.")
    parser.add_argument("--taxon", required=True, help="Taxon name.")
    parser.add_argument("--segment", required=True, help="Segment name.")
    parser.add_argument("--method", required=True, help="Reference creation method (e.g. centroid, consensus).")
    parser.add_argument("--outdir", default='.', help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    os.makedirs(args.outdir, exist_ok=True)

    # Load input data using shared helpers
    qc_data = load_jsonl_gz_by_key(args.qc, "seq_id")
    cluster_data = load_jsonl_gz_by_key(args.clusters, "cluster")
    meta_data = load_jsonl_gz_by_key(args.meta, "accession")

    summary_list: List[Dict[str, Any]] = []

    for cid, cdata in cluster_data.items():
        seq_string = cdata.get("sequence")
        sids = cdata.get("members") or []

        # Get unique accessions
        accessions: List[str] = []
        for sid in sids:
            accs = qc_data.get(sid, {}).get("accessions", [])
            if isinstance(accs, str):
                accs = [accs]
            accessions.extend([a for a in accs if a])

        unique_accessions = list(dict.fromkeys(accessions))

        # Merge metadata values for all accessions
        per_field_values: Dict[str, List[Any]] = defaultdict(list)
        for aid in unique_accessions:
            fields = meta_data.get(aid)
            if fields:
                for k, v in fields.items():
                    per_field_values[k].append(v)

        metadata_summary = summarize_metadata(per_field_values)
        metadata_summary['accessions'] = unique_accessions

        summary_row = {
            "taxon": args.taxon,
            "segment": args.segment,
            "variant": cid,
            "timestamp": timestamp,
            "sequence": seq_string,
            "method": args.method,
            "metadata": metadata_summary
        }

        summary_list.append(summary_row)

    save_output(summary_list, args.taxon, args.segment, args.method, args.outdir)


if __name__ == "__main__":
    main()
