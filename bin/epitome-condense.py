#!/usr/bin/env python3

# epitome-condense.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import sourmash
import screed
import numpy as np
from scipy.spatial.distance import squareform
from sklearn.cluster import DBSCAN
import random
import csv
import datetime
import argparse
import matplotlib.pyplot as plt
import itertools
import json
import os
from collections import defaultdict
import textwrap
import logging

#----- FUNCTIONS -----#
def computeDistance(pair, data):
    """Compute the ANI distance between two sequences identified by a pair of IDs."""
    id1, id2 = pair
    return data[id1]['mh'].containment_ani(data[id2]['mh']).dist


def computeMatrix(data):
    """
    Compute a condensed distance matrix and list of distances for pairwise ANI computations.
    
    Returns:
        - condensed distance matrix (numpy array)
        - sorted list of sequence IDs
        - list of dicts with distances for each pair
    """
    ids = sorted(data.keys())
    pairs = list(itertools.combinations(ids, 2))
    mat = np.zeros((len(ids), len(ids)))
    dist_long = []

    for id1, id2 in pairs:
        dist = computeDistance((id1, id2), data)
        dist_long.append({'id1': id1, 'id2': id2, 'dist': dist})
        i, j = ids.index(id1), ids.index(id2)
        mat[i, j] = dist
        mat[j, i] = dist

    return mat, ids, dist_long


def createClusters(data, threshold, start, stage, outdir):
    mat, ids, _ = computeMatrix(data)
    db = DBSCAN(eps=threshold, min_samples=1, metric='precomputed').fit(mat)
    for id, c in zip(ids, db.labels_):
        data[id]['cluster'] = c
    return data

def selectBest(data):
    """
    Select the best representative sequence per cluster based on 'n' and 'length'.
    
    Returns:
        - list of selected reference sequence IDs
        - updated data with 'condensed' field set for non-reference sequences
    """
    logging.info(f'Selecting best sequences for each cluster.')

    clusters = {}
    for seq_id, info in data.items():
        cluster = info.get('cluster')
        if cluster is not None:
            clusters.setdefault(cluster, []).append(seq_id)

    refs = []
    for cluster, members in clusters.items():
        # Select member with max 'n', then max 'length'
        max_n = max(int(data[m]['n']) for m in members)
        candidates = [m for m in members if int(data[m]['n']) == max_n]

        if len(candidates) > 1:
            max_len = max(int(data[m]['length']) for m in candidates)
            candidates = [m for m in candidates if int(data[m]['length']) == max_len]

        ref = candidates[0]
        refs.append(ref)

        # Mark others as condensed to ref
        for m in members:
            if m != ref:
                data[m]['condensed'] = ref

    return refs, data


def findNeighbor(queries, exclusions, data):
    """
    For each sequence in queries, find its closest neighbor (based on ANI) excluding those in exclusions.
    
    Returns:
        - updated data with 'neighbor_ani' and 'neighbor' fields
        - highest neighbor ANI found
    """
    _, _, dist_long = computeMatrix(data)
    closest_neighbor_ani = 0

    for seq_id in queries:
        min_dist = 1
        neighbors = []

        for pair in dist_long:
            if seq_id == pair['id1']:
                other_id, dist = pair['id2'], pair['dist']
            elif seq_id == pair['id2']:
                other_id, dist = pair['id1'], pair['dist']
            else:
                continue

            if dist < min_dist:
                neighbors = [other_id]
                min_dist = dist
            elif dist == min_dist:
                neighbors.append(other_id)

        neighbors = list(set(neighbors))
        ani = round(100 * (1 - min_dist))
        data[seq_id]['neighbor_ani'] = ani
        data[seq_id]['neighbor'] = neighbors if ani > 0 else []

        if ani > closest_neighbor_ani:
            closest_neighbor_ani = ani

    return data, closest_neighbor_ani


def dictToCsv(dict_list, filename):
    """
    Write a list of dictionaries to a CSV file with specified headers.
    """
    headers = ["taxon", "segment", "cluster", "ref", "condensed", "n", "length", "neighbor_ani", "neighbor"]

    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(dict_list)


def condenseSeqs(windows, data, threshold, stage, outdir):
    """
    Perform one stage of condensing sequences by clustering windows and selecting best representatives.
    
    Returns:
        - list of references (representative sequence IDs)
        - updated data dict
        - status string ('done' or 'not done')
    """
    print(f'\n{datetime.datetime.now()}: Stage {stage}')
    result = data.copy()
    seq_count = len(data)

    if seq_count == 1:
        logging.info(f'Only one sequence provided.')
        return list(data.keys()), data, 'done'

    leftover = [k for k, v in result.items() if v.get('condensed', '') == '']

    for window_key, window_data in windows.items():
        minhash = {k: window_data[k] for k in leftover if k in window_data}
        if len(minhash) == seq_count:
            logging.info(f'Clustering window {window_key}')
            window_clusters = createClusters(minhash, threshold, 0, f'{window_key}_{stage}', outdir)

            for seq_id, cluster_id in window_clusters.items():
                if 'cluster' in result[seq_id]:
                    result[seq_id]['cluster'] = str(result[seq_id]['cluster']) + str(cluster_id)
                else:
                    result[seq_id]['cluster'] = str(cluster_id)

    clusters = [v['cluster'] for v in result.values() if 'cluster' in v]
    logging.info(f'Created {len(set(clusters))} merged clusters!')

    if len(leftover) == len(set(clusters)):
        logging.info(f'Nothing to condense!')
        return list(result.keys()), result, 'done'

    refs, updated_data = selectBest(result)

    # Update condensed info
    for seq_id, info in updated_data.items():
        result[seq_id]['condensed'] = info.get('condensed', '')

    # Clean up cluster keys
    for info in result.values():
        info.pop('cluster', None)

    return refs, result, 'not done'

def dict2Json(data, taxon, segment, fieldnames, filename):
    """
    Restructure flat list of data entries into nested JSON:
      taxon -> segment -> sequence_id -> metadata fields
    """
    logging.info(f'Writing metrics to {filename}')
    nested = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for key, value in data.items():
        nested[taxon][segment][key] = { k: value[k] for k in fieldnames }

    with open(filename, 'w') as f:
        json.dump(nested, f, indent=4)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main():
    version = "1.1"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, nargs='+', help="Path to FASTA file.")
    parser.add_argument("--taxon", default='null', help="Taxon name or 'version'")
    parser.add_argument("--segment", default='null', help="Segment name")
    parser.add_argument("--clusters", type=str, help="Path to the cluster info.")
    parser.add_argument("--dist", default=0.02, type=float,
                        help="Distance threshold (1-%ANI/100) (default: 0.02)")
    parser.add_argument("--scaled", default=1000, type=int, help="Scaled value used by Sourmash")
    parser.add_argument("--ksize", default=31, type=int, help="K-mer size used by Sourmash")
    parser.add_argument("--window", default=8000, type=int, 
                        help="Window size to create clusters. Sequences exceeding this will be split.")
    parser.add_argument("--outdir", default='./', type=str, help="Output directory")
    parser.add_argument('--version', action='version', version=version)
    args = parser.parse_args()

    start = f"""
    epitome-condense.py v{version}
    Written by Jared Johnson
    jared.johnson@doh.wa.gov
    """
    print(textwrap.dedent(start))
    os.makedirs(args.outdir, exist_ok=True)

    # Load sequence cluster sizes
    seqs_per_cluster = {}
    with open(args.clusters) as f:
        clusters = json.load(f)
        for cluster, seqs in clusters[args.taxon][args.segment].items():
            seqs_per_cluster[cluster] = len(seqs)

    # Read FASTA and generate MinHash sketches
    seqs, data = {}, {}
    window_size = args.window
    for fasta_file in args.fasta:
        for record in screed.open(fasta_file):
            name, sequence = record.name, record.sequence
            seqs[name] = sequence
            length = len(sequence)
            window_size = min(window_size, length)
            cluster = name.split('-')[-1]
            data[name] = {'condensed': '', 'length': length, 'n': seqs_per_cluster[cluster]}

    # Generate MinHash for windows
    windows, minhashes = {}, {}
    for name, sequence in seqs.items():
        seq_len = len(sequence)
        j = 0
        mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
        mh.add_sequence(sequence, True)
        data[name]['mh'] = mh
        for i in range(window_size - 1, seq_len, window_size):
            win_mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
            win_mh.add_sequence(sequence[j:(i - 1)], True)
            minhashes[name] = {'mh': win_mh}
            window = f'{j}-{i}'
            windows.setdefault(window, {}).update(minhashes)
            j = i

    # Perform iterative condensing
    stage = 1
    refs, data, status = condenseSeqs(windows, data, args.dist, stage, args.outdir)
    while status != 'done':
        stage += 1
        subset = {k: data[k] for k in refs}
        refs, subset, status = condenseSeqs(windows, subset, args.dist, stage, args.outdir)
        for key, value in subset.items():
            data[key] = value
            if value['condensed'] != '':
                for key2, value2 in data.items():
                    if value2['condensed'] == key:
                        print(f'Update ({key2}): {value2["condensed"]} --> {value["condensed"]}')
                        value2['condensed'] = value['condensed']

    # Final processing
    refs = list(set(refs))
    condensed = [k for k in data if k not in refs]
    print(f'Finding next closest neighbors')
    data, closest_neighbor = findNeighbor(refs, condensed, data)
    if condensed:
        data, _ = findNeighbor(condensed, [], data)
    
    prefix = f'{args.taxon}-{args.segment}'
    dict2Json(data, args.taxon, args.segment, ["condensed", "neighbor_ani", "neighbor"], f'{args.outdir}/{prefix}.condensed.json' )

    end = f"""
    Raw: {len(minhashes)} seqs
    Condensed: {len(refs)} seqs
    Closest Neighbor: {closest_neighbor}% ANI
    
    Results saved to {args.outdir}
    """
    print(textwrap.dedent(end))

if __name__ == "__main__":
    main()