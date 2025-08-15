#!/usr/bin/env python3

# epitome_cluster.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import json
import gzip
import argparse
import random
from collections import defaultdict
from pathlib import Path
import os
from typing import Dict, Tuple, List, Any, Optional

import screed
import sourmash
import numpy as np
from sklearn.cluster import DBSCAN

from epitome_utils import sanitize_filename, DistanceCache, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------
LOGGER = logging_config()


# -------------------------------
#  FUNCTIONS
# -------------------------------

def build_full_minhash_map(seqs: Dict[str, str], ksize: int, scaled: int) -> Dict[str, sourmash.MinHash]:
    """
    Build a MinHash object for each sequence.

    Args:
        seqs: Mapping of sequence IDs to nucleotide sequences.
        ksize: K-mer size for MinHash.
        scaled: Scaling factor for MinHash.

    Returns:
        Dict mapping sequence ID to sourmash.MinHash object.
    """
    mh_map = {}
    for sid, s in seqs.items():
        mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        mh.add_sequence(s, force=True)
        mh_map[sid] = mh
    return mh_map


def load_seqs(fasta: str, ksize: int, scaled: int, window_size: int) -> Tuple[Dict[str, Dict[str, Dict[str, Any]]], Dict[str, str]]:
    """
    Load sequences from FASTA and split into MinHash windows.

    Args:
        fasta: Path to FASTA file.
        ksize: K-mer size for MinHash.
        scaled: Scaling factor for MinHash.
        window_size: Window length to split sequences.

    Returns:
        Tuple of:
            windows: Dict[window_name -> {seq_id -> {"mh": MinHash}}]
            sequences: Dict[seq_id -> sequence string]
    """
    sequences = {rec.name: rec.sequence for rec in screed.open(fasta)}
    LOGGER.info(f'Loaded {len(sequences)} sequences from {fasta}')

    windows = defaultdict(dict)
    seq_len = min(len(v) for v in sequences.values())
    for seq_id, seq in sequences.items():
        win_size = min(window_size, seq_len)
        start = 0
        for end in range(win_size, seq_len + win_size, win_size):
            mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(seq[start:end], force=True)
            window_name = f'{start}-{end}'
            windows[window_name][seq_id] = {'mh': mh}
            start = end

    LOGGER.info(f'Split sequences into {len(windows)} windows')

    shared_windows = {
        w: seqs for w, seqs in windows.items()
        if len(seqs) == len(sequences)
    }
    if len(windows) != len(shared_windows):
        LOGGER.info(f'Reduced to {len(shared_windows)} shared windows')

    return shared_windows, sequences


def main_and_remainder(data: Dict[str, dict], threshold: int, round_num: int) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """
    Split a dataset into main and remainder subsets based on threshold.

    Args:
        data: Mapping of IDs to sequence metadata dicts.
        threshold: Max number of items in main set.
        round_num: Current clustering round number.

    Returns:
        Tuple of main subset and remainder subset.
    """
    if len(data) > threshold:
        random_keys = random.sample(list(data.keys()), threshold)
        subset = {k: data[k] for k in random_keys}
        for v in subset.values():
            v['source'] = f'main_{round_num}'
        remainder = {k: data[k] for k in data if k not in subset}
        for v in remainder.values():
            v['source'] = f'assigned_{round_num}'
    else:
        subset = data
        for v in subset.values():
            v['source'] = f'main_{round_num}'
        remainder = {}
    return subset, remainder


def select_reps(data: Dict[str, dict]) -> Dict[str, dict]:
    """
    Select representative sequences for each cluster.

    Args:
        data: Mapping of IDs to sequence metadata dicts (with 'cluster' key).

    Returns:
        Dict mapping representative ID to metadata.
    """
    result, seen = {}, set()
    for k, v in data.items():
        c = v['cluster']
        if c not in seen:
            result[k] = v
            seen.add(c)
    return result


def compute_matrix(data: Dict[str, dict], window_scope: str, dist_cache: DistanceCache) -> Tuple[np.ndarray, List[str]]:
    """
    Compute pairwise distance matrix for sequences in a window.

    Args:
        data: Mapping of IDs to dict containing 'mh' key.
        window_scope: Identifier for caching distances.
        dist_cache: DistanceCache object.

    Returns:
        Tuple of (distance matrix, list of IDs in order).
    """
    ids = sorted(data)
    n = len(ids)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = ids[i], ids[j]
            cached = dist_cache.get(window_scope, a, b)
            if cached is None:
                d = data[a]['mh'].containment_ani(data[b]['mh']).dist
                dist_cache.set(window_scope, a, b, d)
            else:
                d = cached
            mat[i, j] = mat[j, i] = d
    np.fill_diagonal(mat, 0)
    return mat, ids


def create_clusters(data: Dict[str, dict], threshold: float, start: int, stage: str, outdir: str, window_scope: str, dist_cache: DistanceCache) -> Tuple[Dict[str, dict], int]:
    """
    Cluster sequences using DBSCAN.

    Args:
        data: Mapping of IDs to metadata dicts.
        threshold: Distance threshold for clustering.
        start: Starting cluster index.
        stage: Stage name for logging.
        outdir: Output directory.
        window_scope: Distance cache scope name.
        dist_cache: DistanceCache object.

    Returns:
        Tuple of updated data and last cluster number used.
    """
    mat, ids = compute_matrix(data, window_scope, dist_cache)
    db = DBSCAN(eps=threshold, min_samples=1, metric='precomputed').fit(mat)
    clusters = [int(i) + 1 for i in db.labels_]
    max_cluster = max(clusters) + start + 1
    for id_, cluster in zip(ids, clusters):
        if cluster == 0:
            max_cluster += 1
            cluster = max_cluster
        data[id_]['cluster'] = cluster
    return data, max_cluster


def assign_clusters(remainder: Dict[str, dict], reps: Dict[str, dict], data: Dict[str, dict], threshold: float, window_scope: str, dist_cache: DistanceCache) -> Dict[str, dict]:
    """
    Assign remaining sequences to nearest cluster representative.

    Args:
        remainder: Unassigned sequences.
        reps: Cluster representatives.
        data: Existing clustered data.
        threshold: Distance threshold.
        window_scope: Distance cache scope.
        dist_cache: DistanceCache object.

    Returns:
        Updated data dict with assigned clusters.
    """
    for k1, v1 in remainder.items():
        min_dist = 1.0
        assigned_cluster = None
        for k2, v2 in reps.items():
            cached = dist_cache.get(window_scope, k1, k2)
            if cached is None:
                d = v1['mh'].containment_ani(v2['mh']).dist
                dist_cache.set(window_scope, k1, k2, d)
            else:
                d = cached
            if d < threshold and d < min_dist:
                min_dist = d
                assigned_cluster = v2['cluster']
        if assigned_cluster is not None:
            v1['cluster'] = assigned_cluster
            data[k1] = v1
    return data


def cluster_seqs(data: Dict[str, dict], max_cluster: int, threshold: float, round_num: str, start: int, outdir: str, window: str, dist_cache: DistanceCache) -> Tuple[Dict[str, dict], Dict[str, dict], int]:
    """
    Perform clustering and assignment in multiple rounds if necessary.

    Args:
        data: Mapping of IDs to sequence metadata.
        max_cluster: Max cluster size per round.
        threshold: Distance threshold for DBSCAN.
        round_num: Current round identifier.
        start: Starting cluster index.
        outdir: Output directory.
        window: Window identifier.
        dist_cache: DistanceCache object.

    Returns:
        Tuple of (clusters dict, loose/unassigned dict, last cluster number).
    """
    if start != 0:
        LOGGER.info(f'Round {round_num} - Continuing from cluster {start}')
    if len(data) > 1:
        main, rem = main_and_remainder(data, max_cluster, round_num)
        LOGGER.info(f'Round {round_num} - Clustering {len(main)} sequences')
        clusters, last_cluster = create_clusters(
            main, threshold, start, f'main_{round_num}', outdir,
            window_scope=f"win:{window}", dist_cache=dist_cache
        )
        LOGGER.info(f'Round {round_num} - Clustered into {last_cluster - start} clusters')
        if rem:
            reps = select_reps(clusters)
            LOGGER.info(f'Round {round_num} - Assigning {len(rem)} sequences')
            clusters = assign_clusters(rem, reps, clusters, threshold, f"win:{window}", dist_cache)
        loose = {k: rem[k] for k in rem if k not in clusters}
    else:
        loose, last_cluster = {}, start + 1
        for v in data.values():
            v['cluster'] = last_cluster
        clusters = data
    return clusters, loose, last_cluster


def compute_centroids(seq_ids: List[str], full_mh_map: Dict[str, sourmash.MinHash], dist_cache: DistanceCache, seqs: Dict[str, str]) -> Tuple[Optional[str], str]:
    """
    Compute centroid sequence ID for a cluster based on average distance.

    Args:
        seq_ids: List of sequence IDs in cluster.
        full_mh_map: Mapping of IDs to MinHash objects.
        dist_cache: DistanceCache object.
        seqs: Mapping of IDs to sequences.

    Returns:
        Tuple of (centroid sequence ID or None, sequence string).
    """
    if not seq_ids:
        return None, ""
    if len(seq_ids) == 1:
        sid = seq_ids[0]
        return sid, seqs.get(sid, "")

    min_avg = float('inf')
    cen = None
    for sid in seq_ids:
        dsum, cnt = 0.0, 0
        for oid in seq_ids:
            if sid == oid:
                continue
            cached = dist_cache.get("full", sid, oid)
            if cached is None:
                d = full_mh_map[sid].containment_ani(full_mh_map[oid]).dist
                dist_cache.set("full", sid, oid, d)
            else:
                d = cached
            dsum += d
            cnt += 1
        avg = dsum / cnt if cnt else 0.0
        if avg < min_avg:
            min_avg = avg
            cen = sid
    return cen, seqs.get(cen, "")


def save_output(json_data: List[dict], seq_data: Dict[str, str], taxon: str, segment: str, outdir: str) -> None:
    """
    Save clustering results and sequences to disk.

    Args:
        json_data: List of cluster result dicts.
        seq_data: Mapping of IDs to sequences.
        taxon: Taxon name.
        segment: Segment name.
        outdir: Output directory.

    Returns:
        None
    """
    safe_taxon = sanitize_filename(taxon)
    safe_segment = sanitize_filename(segment)
    prefix = Path(outdir) / f'{safe_taxon}-{safe_segment}'

    filename = prefix.with_suffix('.clusters.jsonl.gz')
    with gzip.open(filename, 'wt', encoding='utf-8') as f:
        for row in json_data:
            if isinstance(row, dict):
                f.write(json.dumps(row, sort_keys=True) + '\n')
    LOGGER.info(f"Saved cluster metrics to {filename}")

    for row in json_data:
        ids = row.get("members", [])
        cluster = row.get('cluster')
        var = 'multi' if len(ids) > 1 else 'single'
        fname = prefix.with_name(f"{prefix.name}-{cluster}.{var}.fa.gz")
        if var == 'single':
            out = f'>{prefix.name}-{cluster}\n{seq_data[ids[0]]}\n'
        else:
            out = '\n'.join([f'>{id_}\n{seq_data[id_]}' for id_ in ids]) + '\n'
        with gzip.open(fname, 'wt', encoding='utf-8') as f:
            f.write(out)
        LOGGER.info(f"Saved cluster sequence(s) to {fname}")


# -------------------------------
#  MAIN
# -------------------------------

def main():
    """
    Command-line entry point for EPITOME sequence clustering.

    Parses arguments, loads sequences, performs windowed MinHash-based
    clustering with DBSCAN, optionally computes cluster centroids, and
    writes output files.
    """
    version = "2.0"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Path to mutli-FASTA file.")
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--segment", default="null", help="Segment name.")
    parser.add_argument("--dist", type=float, default=0.02, help="Distance threshold (1 - ANI/100).")
    parser.add_argument("--scaled", type=int, default=1000, help="MinHash scaled factor.")
    parser.add_argument("--ksize", type=int, default=31, help="MinHash k-mer size.")
    parser.add_argument("--window_size", type=int, default=8000, help="Window size for MinHash splitting.")
    parser.add_argument("--max_cluster", type=int, default=1000, help="Maximum number of sequences to cluster per round.")
    parser.add_argument("--centroid", action="store_true", help="Calculate the centroid for each cluster. Centroid sequences are returned in the 'sequence' field.")
    parser.add_argument("--outdir", default='.', help="Output directory.")
    parser.add_argument("--seed", type=int, default=11, help="Random seed.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    random.seed(args.seed)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    windows, seqs = load_seqs(args.fasta, args.ksize, args.scaled, args.window_size)

    dist_cache = DistanceCache()
    full_mh_map = build_full_minhash_map(seqs, args.ksize, args.scaled)

    results = defaultdict(list)
    for window, data in windows.items():
        start, round_num, win_res = 0, 1, {}
        while data:
            clusters, data, start = cluster_seqs(
                data, args.max_cluster, args.dist, f'{round_num}_{window}',
                start, args.outdir, window, dist_cache
            )
            win_res.update(clusters)
            round_num += 1
        for sample, info in win_res.items():
            results[sample].append(str(info['cluster']))

    unique, new_id = {}, 1
    cluster_map = {}
    for key, cluster_ids in results.items():
        value = ":".join(cluster_ids)
        if value not in unique:
            unique[value] = new_id
            new_id += 1
        cid = unique[value]
        cluster_map.setdefault(cid, []).append(key)

    cluster_list = []
    for key, value in cluster_map.items():
        row = {
            'taxon': args.taxon,
            'segment': args.segment,
            'cluster': key,
            'members': value
        }
        cluster_list.append(row)

    if args.centroid:
        for i, row in enumerate(cluster_list):
            members = row.get("members", [])
            cid, ref = compute_centroids(members, full_mh_map, dist_cache, seqs)
            row["sequence"] = ref
            cluster_list[i] = row

    save_output(cluster_list, seqs, args.taxon, args.segment, args.outdir)

    # Summary
    LOGGER.info(f"Total input sequences: {len(seqs)}")
    LOGGER.info(f"Total clusters: {new_id - 1}")
    LOGGER.info(f"Results saved to {args.outdir}")


if __name__ == "__main__":
    main()
