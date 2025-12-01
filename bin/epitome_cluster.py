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
#  HELPERS (de-duplicated)
# -------------------------------

def _distance(
    a: str,
    b: str,
    *,
    scope: str,
    mh_map: Dict[str, sourmash.MinHash],
    cache: DistanceCache
) -> float:
    """
    Cached pairwise distance between two IDs for a given scope ('full' or 'win:XYZ').
    """
    cached = cache.get(scope, a, b)
    if cached is not None:
        return float(cached)
    d = mh_map[a].containment_ani(mh_map[b]).dist
    cache.set(scope, a, b, float(d))
    return float(d)


# -------------------------------
#  FUNCTIONS
# -------------------------------

def build_full_minhash_map(seqs: Dict[str, str], ksize: int, scaled: int) -> Dict[str, sourmash.MinHash]:
    """
    Build a MinHash object for each sequence.
    """
    mh_map: Dict[str, sourmash.MinHash] = {}
    for sid, s in seqs.items():
        mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        mh.add_sequence(s, force=True)
        mh_map[sid] = mh
    return mh_map


def load_seqs(fasta: str, ksize: int, scaled: int, window_size: int) -> Tuple[Dict[str, Dict[str, Dict[str, Any]]], Dict[str, str]]:
    """
    Load sequences from FASTA and split into MinHash windows.
    """
    sequences = {rec.name: rec.sequence for rec in screed.open(fasta)}
    LOGGER.info(f'Loaded {len(sequences)} sequences from {fasta}')

    windows: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
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

    shared_windows = {w: seqs for w, seqs in windows.items() if len(seqs) == len(sequences)}
    if len(windows) != len(shared_windows):
        LOGGER.info(f'Reduced to {len(shared_windows)} shared windows')

    return shared_windows, sequences


def main_and_remainder(data: Dict[str, dict], threshold: int, round_num: int) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """
    Split a dataset into main and remainder subsets based on threshold.
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
    Compute pairwise distance matrix for sequences in a window using cached distances.
    """
    ids = sorted(data)
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    mh_map = {sid: data[sid]['mh'] for sid in ids}
    for i in range(n):
        ai = ids[i]
        for j in range(i + 1, n):
            aj = ids[j]
            d = _distance(ai, aj, scope=window_scope, mh_map=mh_map, cache=dist_cache)
            mat[i, j] = mat[j, i] = d
    np.fill_diagonal(mat, 0.0)
    return mat, ids


def create_clusters(
    data: Dict[str, dict],
    threshold: float,
    start: int,
    window_scope: str,
    dist_cache: DistanceCache,
    *,
    min_samples_primary: int = 10,
    min_samples_secondary: int = 1
) -> Tuple[Dict[str, dict], int]:
    """
    Cluster sequences using a two-stage DBSCAN:
      1) Run DBSCAN over all points with a larger min_samples (min_samples_primary).
      2) Run DBSCAN over stage-1 noise only with a smaller min_samples (min_samples_secondary).
    Any remaining noise after stage 2 becomes singleton clusters.

    Cluster numbering continues from `start`.
    Returns (updated_data, last_cluster_id).
    """
    mat, ids = compute_matrix(data, window_scope, dist_cache)
    n = len(ids)
    if n == 0:
        return data, start

    db1 = DBSCAN(eps=threshold, min_samples=min_samples_primary, metric='precomputed').fit(mat)
    labels1 = db1.labels_  # -1 for noise, >=0 for clusters

    # Map non-noise labels to consecutive cluster IDs continuing from `start`
    current_max = start
    label_map_1: Dict[int, int] = {}
    for lab in sorted(set(labels1)):
        if lab == -1:
            continue
        current_max += 1
        label_map_1[lab] = current_max

    # Assign stage-1 non-noise clusters; collect noise indices for stage 2
    noise_indices: List[int] = []
    for idx, lab in enumerate(labels1):
        sid = ids[idx]
        if lab == -1:
            noise_indices.append(idx)
        else:
            data[sid]['cluster'] = label_map_1[lab]

    # ---- Stage 2: Re-cluster noise only (if any) ----
    if noise_indices:
        submat = mat[np.ix_(noise_indices, noise_indices)]
        subids = [ids[i] for i in noise_indices]

        db2 = DBSCAN(eps=threshold, min_samples=min_samples_secondary, metric='precomputed').fit(submat)
        labels2 = db2.labels_

        # Map stage-2 non-noise labels to new cluster IDs (continue counting)
        label_map_2: Dict[int, int] = {}
        for lab in sorted(set(labels2)):
            if lab == -1:
                continue
            current_max += 1
            label_map_2[lab] = current_max

        # Assign stage-2 clusters; leftover noise -> singleton clusters
        for idx2, lab2 in enumerate(labels2):
            sid = subids[idx2]
            if lab2 == -1:
                current_max += 1
                data[sid]['cluster'] = current_max
            else:
                data[sid]['cluster'] = label_map_2[lab2]

    return data, current_max


def assign_clusters(
    remainder: Dict[str, dict],
    reps: Dict[str, dict],
    data: Dict[str, dict],
    threshold: float,
    window_scope: str,
    dist_cache: DistanceCache
) -> Dict[str, dict]:
    """
    Assign remaining sequences to nearest cluster representative.
    """
    if not remainder:
        return data

    mh_map: Dict[str, sourmash.MinHash] = {}
    mh_map.update({k: v['mh'] for k, v in remainder.items()})
    mh_map.update({k: v['mh'] for k, v in reps.items()})

    for k1, v1 in remainder.items():
        min_dist = 1.0
        assigned_cluster = None
        for k2, v2 in reps.items():
            d = _distance(k1, k2, scope=window_scope, mh_map=mh_map, cache=dist_cache)
            if d < threshold and d < min_dist:
                min_dist = d
                assigned_cluster = v2['cluster']
        if assigned_cluster is not None:
            v1['cluster'] = assigned_cluster
            data[k1] = v1
    return data


def cluster_seqs(
    data: Dict[str, dict],
    max_per_round: int,
    threshold: float,
    round_num: str,
    start: int,
    outdir: str,
    window: str,
    dist_cache: DistanceCache,
    min_samples_primary: int,
    min_samples_secondary: int
) -> Tuple[Dict[str, dict], Dict[str, dict], int]:
    """
    Perform clustering and assignment in multiple rounds if necessary.
    """
    if start != 0:
        LOGGER.info(f'Round {round_num} - Continuing from cluster {start}')
    if len(data) > 1:
        main, rem = main_and_remainder(data, max_per_round, round_num)
        LOGGER.info(f'Round {round_num} - Clustering {len(main)} sequences')
        clusters, last_cluster = create_clusters(
            main, threshold, start,
            window_scope=f"win:{window}",
            dist_cache=dist_cache,
            min_samples_primary=min_samples_primary,
            min_samples_secondary=min_samples_secondary
        )
        LOGGER.info(f'Round {round_num} - Clustered into {last_cluster - start} clusters')
        if rem:
            reps = select_reps(clusters)
            LOGGER.info(f'Round {round_num} - Assigning {len(rem)} sequences')
            clusters = assign_clusters(
                rem, reps, clusters, threshold, f"win:{window}", dist_cache
            )
        loose = {k: rem[k] for k in rem if k not in clusters}
    else:
        loose, last_cluster = {}, start + 1
        for v in data.values():
            v['cluster'] = last_cluster
        clusters = data
    return clusters, loose, last_cluster


def compute_centroids(
    seq_ids: List[str],
    full_mh_map: Dict[str, sourmash.MinHash],
    dist_cache: DistanceCache,
    seqs: Dict[str, str]
) -> Tuple[Optional[str], str]:
    """
    Compute centroid sequence ID for a cluster based on average distance.
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
            d = _distance(sid, oid, scope="full", mh_map=full_mh_map, cache=dist_cache)
            dsum += d
            cnt += 1
        avg = dsum / cnt if cnt else 0.0
        if avg < min_avg:
            min_avg = avg
            cen = sid
    return cen, seqs.get(cen, "")


def compute_cluster_min_max(
    members: List[str],
    full_mh_map: Dict[str, sourmash.MinHash],
    dist_cache: DistanceCache
) -> Tuple[float, float]:
    """
    Compute the minimum and maximum pairwise distances among cluster members
    using the full-length MinHashes (scope 'full').

    Returns (min_dist, max_dist). For singleton clusters, returns (0.0, 0.0).
    """
    if not members or len(members) == 1:
        return 0.0, 0.0

    min_d = float('inf')
    max_d = 0.0
    n = len(members)
    for i in range(n):
        a = members[i]
        for j in range(i + 1, n):
            b = members[j]
            d = _distance(a, b, scope="full", mh_map=full_mh_map, cache=dist_cache)
            if d < min_d:
                min_d = d
            if d > max_d:
                max_d = d

    if min_d == float('inf'):
        min_d = 0.0
    return float(min_d), float(max_d)


def save_output(json_data: List[dict], seq_data: Dict[str, str], taxon: str, segment: str, outdir: str) -> None:
    """
    Save clustering results and sequences to disk.
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
        ids = row.get("subset", [])
        cluster = row.get('cluster')
        var = 'multi' if len(ids) > 1 else 'single'
        fname = prefix.with_name(f"{prefix.name}-{cluster}.{var}.fa.gz")
        if var == 'single':
            out = f'>{prefix.name}-{cluster}\n{seq_data[ids[0]]}\n'
        else:
            out = '\n'.join([f'>{id_}\n{seq_data[id_]}' for id_ in ids]) + '\n'
        with gzip.open(fname, 'wt', encoding='utf-8') as f:
            f.write(out)
    LOGGER.info("Saved cluster sequence(s)")


# -------------------------------
#  TROUBLESHOOTING
# -------------------------------

def _compute_full_distance_matrix_for_members(
    members: List[str],
    full_mh_map: Dict[str, sourmash.MinHash],
    dist_cache: DistanceCache
) -> Tuple[np.ndarray, List[str]]:
    """
    Build a pairwise distance matrix for the given members using the 'full' sketches.
    Returns (matrix, ids) where ids is the row/col order.
    """
    ids = sorted(members)
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        ai = ids[i]
        for j in range(i + 1, n):
            aj = ids[j]
            d = _distance(ai, aj, scope="full", mh_map=full_mh_map, cache=dist_cache)
            mat[i, j] = mat[j, i] = float(d)
    np.fill_diagonal(mat, 0.0)
    return mat, ids


def _save_distance_matrix_tsv_gz(
    ids: List[str],
    mat: np.ndarray,
    out_path: Path
) -> None:
    """
    Save a distance matrix as gzipped TSV with header:
      first row:  ID <tab> id1 <tab> id2 ...
      rows:       idX <tab> d11 <tab> d12 ...
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8") as f:
        f.write("ID\t" + "\t".join(ids) + "\n")
        for i, row_id in enumerate(ids):
            row_vals = [f"{mat[i, j]:.6f}" for j in range(len(ids))]
            f.write(row_id + "\t" + "\t".join(row_vals) + "\n")


def export_problematic_clusters(
    *,
    taxon: str,
    segment: str,
    cluster_id: int,
    members: List[str],
    threshold: float,
    full_mh_map: Dict[str, sourmash.MinHash],
    dist_cache: DistanceCache,
    outdir: str
) -> Optional[Path]:
    """
    If the max intra-cluster distance of 'members' exceeds 'threshold',
    compute and save a gzipped TSV distance matrix. Returns output path or None.
    """
    if len(members) <= 1:
        return None

    _, max_d = compute_cluster_min_max(members, full_mh_map, dist_cache)
    if max_d <= threshold:
        return None

    mat, ids = _compute_full_distance_matrix_for_members(members, full_mh_map, dist_cache)

    safe_taxon = sanitize_filename(taxon)
    safe_segment = sanitize_filename(segment)
    out_path = Path(outdir) / f"{safe_taxon}-{safe_segment}-cluster{cluster_id}.dists.tsv.gz"
    _save_distance_matrix_tsv_gz(ids, mat, out_path)
    return out_path


# -------------------------------
#  MAIN
# -------------------------------

def main():
    """
    Command-line entry point for EPITOME sequence clustering.
    """
    version = "2.0"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Path to mutli-FASTA file.")
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--segment", default="null", help="Segment name.")
    parser.add_argument("--dist", type=float, default=0.02, help="Distance threshold (1 - ANI/100).")
    parser.add_argument("--scaled", type=int, default=100, help="MinHash scaled factor.")
    parser.add_argument("--ksize", type=int, default=31, help="MinHash k-mer size.")
    parser.add_argument("--window_size", type=int, default=8000, help="Window size for MinHash splitting.")
    parser.add_argument("--max_per_round", type=int, default=1000, help="Maximum number of sequences to cluster per round.")
    parser.add_argument("--max_per_cluster", type=int, default=3000, help="Maximum number of sequences to include in a cluster.")
    parser.add_argument("--min_samples_primary", type=int, default=5, help="Minimum numbers of samples for to create a cluster on primary DBSCAN stage.")
    parser.add_argument("--min_samples_secondary", type=int, default=1, help="Minimum numbers of samples for to create a cluster on secondary DBSCAN stage (noise from primary stage).")
    parser.add_argument("--centroid", action="store_true", help="Calculate the centroid for each cluster. Centroid sequences are returned in the 'sequence' field.")
    parser.add_argument("--outdir", default='.', help="Output directory.")
    parser.add_argument("--seed", type=int, default=11, help="Random seed.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    LOGGER.info(vars(args))

    random.seed(args.seed)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    windows, seqs = load_seqs(args.fasta, args.ksize, args.scaled, args.window_size)

    dist_cache = DistanceCache()
    full_mh_map = build_full_minhash_map(seqs, args.ksize, args.scaled)

    results: Dict[str, List[str]] = defaultdict(list)
    for window, data in windows.items():
        start, round_num, win_res = 0, 1, {}
        while data:
            clusters, data, start = cluster_seqs(
                data, args.max_per_round, args.dist, f'{round_num}_{window}',
                start, args.outdir, window, dist_cache,
                args.min_samples_primary, args.min_samples_secondary
            )
            win_res.update(clusters)
            round_num += 1
        for sample, info in win_res.items():
            results[sample].append(str(info['cluster']))

    unique, new_id = {}, 1
    cluster_map: Dict[int, List[str]] = {}
    for key, cluster_ids in results.items():
        value = ":".join(cluster_ids)
        if value not in unique:
            unique[value] = new_id
            new_id += 1
        cid = unique[value]
        cluster_map.setdefault(cid, []).append(key)

    cluster_list: List[dict] = []
    for key, members in cluster_map.items():
        # create subset based on specified max cluster size
        subset = random.sample(members, args.max_per_cluster) if len(members) > args.max_per_cluster else members

        # compute per-cluster min/max distances on full sketches
        min_d, max_d = compute_cluster_min_max(subset, full_mh_map, dist_cache)
        row = {
            'taxon': args.taxon,
            'segment': args.segment,
            'cluster': key,
            'members': members,
            'subset': subset,
            'dist_range': [min_d, max_d]
        }
        cluster_list.append(row)

    max_dists = [c['dist_range'][1] for c in cluster_list]
    avg_max_dist = sum(max_dists) / len(cluster_list) if cluster_list else 0.0
    max_max_dist = max(max_dists) if cluster_list else 0.0
    min_max_dist = min(max_dists) if cluster_list else 0.0

    # Emit distance matrices for problematic clusters
    for row in cluster_list:
        export_problematic_clusters(
            taxon=args.taxon,
            segment=args.segment,
            cluster_id=row["cluster"],
            members=row.get("subset", []),
            threshold=args.dist,
            full_mh_map=full_mh_map,
            dist_cache=dist_cache,
            outdir=args.outdir
        )

    if args.centroid:
        for i, row in enumerate(cluster_list):
            subset = row.get("subset", [])
            cid, ref = compute_centroids(subset, full_mh_map, dist_cache, seqs)
            row["sequence"] = ref
            cluster_list[i] = row

    save_output(cluster_list, seqs, args.taxon, args.segment, args.outdir)

    # Summary
    LOGGER.info(f"Total input sequences: {len(seqs)}")
    LOGGER.info(f"Total clusters: {new_id - 1}")
    LOGGER.info(f"Max. intra-cluster distance: Range = {min_max_dist} - {max_max_dist}; Avg. = {avg_max_dist}")
    LOGGER.info(f"Results saved to {args.outdir}")


if __name__ == "__main__":
    main()
