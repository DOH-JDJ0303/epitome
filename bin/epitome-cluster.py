#!/usr/bin/env python3

# epitome-cluster.py
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
import os
import textwrap
from collections import defaultdict
import json
import logging
import gzip




# ----- FUNCTIONS ----- #
def loadSeqs(fasta, ksize, scaled, window_size):
    # Load sequences
    sequences = {rec.name: rec.sequence for rec in screed.open(fasta)}
    logging.info(f'Loaded {len(sequences)} sequences from {fasta}')
    
    # Initialize windows
    windows = defaultdict(dict)
    seq_len = min([ len(v) for v in sequences.values() ])
    for seq_id, seq in sequences.items():
        win_size = min(window_size, seq_len)
        start = 0
        for end in range(win_size, seq_len + win_size, win_size):
            mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(seq[start:end], force=True)
            window_name = f'{start}-{end}'
            windows[window_name][seq_id] = {'mh': mh}
            start = end
    
    logging.info(f'Split sequences into {len(windows)} windows')

    # Keep only shared windows
    shared_windows = {
        w: seqs for w, seqs in windows.items()
        if len(seqs) == len(sequences)
    }
    if len(windows) != len(shared_windows): 
        logging.info(f'Reduced to {len(shared_windows)} shared windows')

    return shared_windows, sequences

def mainAndRemainder(data, threshold, round):
    if len(data) > threshold:
        random_keys = random.sample(list(data.keys()), threshold)
        subset = {k: data[k] for k in random_keys}
        for v in subset.values(): v['source'] = f'main_{round}'
        remainder = {k: data[k] for k in data if k not in subset}
        for v in remainder.values(): v['source'] = f'assigned_{round}'
    else:
        subset = data
        for v in subset.values(): v['source'] = f'main_{round}'
        remainder = {}
    return subset, remainder

def selectReps(data):
    result, seen = {}, set()
    for k, v in data.items():
        c = v['cluster']
        if c not in seen:
            result[k] = v
            seen.add(c)
    return result

def computeDistance(pair, data):
    id1, id2 = pair
    x = data[id1]['mh'].containment_ani(data[id2]['mh'])
    return {'id1': id1, 'id2': id2, 'dist': x.dist}

def computeMatrix(data):
    ids = sorted(data)
    pairs = list(itertools.combinations(ids, 2))
    dist_long = [computeDistance(p, data) for p in pairs]
    mat = np.zeros((len(ids), len(ids)))
    for d in dist_long:
        i, j = ids.index(d['id1']), ids.index(d['id2'])
        mat[i][j] = mat[j][i] = d['dist']
    np.fill_diagonal(mat, 0)
    return mat, ids

def createClusters(data, threshold, start, stage, outdir):
    mat, ids = computeMatrix(data)
    db = DBSCAN(eps=threshold, min_samples=1, metric='precomputed').fit(mat)
    clusters = [ int(i) + 1 for i in db.labels_ ]
    max_cluster = max(clusters) + start
    for id, cluster in zip(ids, clusters):
        if cluster == 0:
            c = max_cluster = max_cluster
        else:
            c = int(cluster) + start
        data[id]['cluster'] = c
    return data, max_cluster

def assignClusters(remainder, reps, data, threshold, margin):
    for k1, v1 in remainder.items():
        min_dist = 1
        matches = []

        for k2, v2 in reps.items():
            d1 = v1['mh'].containment_ani(v2['mh']).dist
            if d1 < threshold:
                v1['cluster'] = v2['cluster']
                data[k1] = v1
                break
            elif d1 < margin*threshold and d1 <= min_dist:
                min_dist = d1
                matches.append( v2['cluster'] )

        if k1 in data:
            continue

        for cluster in matches:
            for k3, v3 in data.items():
                if v3['cluster'] == cluster:
                    d2 = v1['mh'].containment_ani(v3['mh']).dist
                    if d2 < threshold:
                        v1['cluster'] = v3['cluster']
                        data[k1] = v1
                        break
            if k1 in data:
                break
    return data

def assignRefs(seq_windows, ref_windows, threshold):
    _, first_seq = next(iter(seq_windows.items()))
    _, first_ref = next(iter(ref_windows.items()))
    window_set = list(set(seq_windows.keys()) & set(ref_windows.keys()))
    assigned = {}
    for seq in first_seq.keys():
        for ref in first_ref.keys():
            assigned_win = {'reference': ref}
            for w in window_set:
                window_dist = seq_windows[w][seq]['mh'].containment_ani(ref_windows[w][ref]['mh']).dist
                if window_dist >= threshold:
                    break
                else:
                    assigned_win[w] = window_dist
            else:
                assigned[seq] = assigned_win
                break

    logging.info(f'Assigned {len(assigned.keys())} sequences to an existing reference.')
    with open('ref_assigned.json', 'w') as f:
        json.dump(assigned, f, indent=2)

    remainder = defaultdict(lambda: defaultdict(dict))
    n_seqs = 0
    for window, data in seq_windows.items():
        for key, value in data.items():
            if key not in assigned.keys():
                n_seqs += 1
                remainder[window][key] = value

    return remainder

def clusterSeqs(data, max_cluster, threshold, assign_margin, round, start, outdir):
    if start != 0:
        logging.info(f'Round {round} - Continuing from cluster {start}')
    if len(data) > 1:
        main, rem = mainAndRemainder(data, max_cluster, round)
        logging.info(f'Round {round} - Clustering {len(main)} sequences')
        clusters, last_cluster = createClusters(main, threshold, start, f'main_{round}', outdir)
        logging.info(f'Round {round} - Clustered into {last_cluster - start} clusters')
        if rem:
            reps = selectReps(clusters)
            logging.info(f'Round {round} - Assigning {len(rem)} sequences')
            n_pre = len(clusters)
            clusters = assignClusters(rem, reps, clusters, threshold, assign_margin)
            logging.info(f'Round {round} - Assiged {len(clusters) - n_pre} sequences')
        loose = {k: rem[k] for k in rem if k not in clusters}
    else:
        loose, last_cluster = {}, start + 1
        for v in data.values(): v['cluster'] = last_cluster
        clusters = data
    return clusters, loose, last_cluster

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main():
    version = "1.1"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, required=True, help="Path to FASTA file.")
    parser.add_argument("--refs", type=str, required=False, help="Path to existing references.")
    parser.add_argument("--taxon", default='null', help="Taxon name or 'version'")
    parser.add_argument("--segment", default='null', help="Segment name")
    parser.add_argument("--dist", default=0.02, type=float, help="Distance threshold (1 - ANI/100)")
    parser.add_argument("--assign_margin", default=2, type=float, help="Distance threshold multiplier when comparing cluster representatives during the assignment phase.")
    parser.add_argument("--max_cluster", default=1000, type=int, help="Max sequences in initial clustering round")
    parser.add_argument("--ksize", default=31, type=int, help="K-mer size for Sourmash")
    parser.add_argument("--scaled", default=10, type=int, help="Scaled value for Sourmash")
    parser.add_argument("--window_size", default=10000, type=int, help="Sequence window size")
    parser.add_argument("--outdir", default='./', type=str, help="Output directory")
    parser.add_argument('--version', action='version', version=version)
    
    args = parser.parse_args()

    start = f"""

    epitome-cluster.py v{version}
    Written by Jared Johnson
    jared.johnson@doh.wa.gov
    """
    print(textwrap.dedent(start))

    # Load sequences
    windows, seqs = loadSeqs(args.fasta, args.ksize, args.scaled, args.window_size)
    
    # Assign sequences to existing references
    if args.refs:
        logging.info(f'Assinging to references')
        # Load references
        windows_refs, _ = loadSeqs(args.refs, args.ksize, args.scaled, args.window_size)
        windows = assignRefs(windows, windows_refs, args.dist)

    # Create clusters
    results = {}
    for window, data in windows.items():
        logging.info(f'Clustering window {window}')
        start, round, win_res = 0, 1, {}
        while data:
            clusters, data, start = clusterSeqs(data, args.max_cluster, args.dist, args.assign_margin, f'{round}_{window}', start, args.outdir)
            win_res.update(clusters)
            round += 1
        for sample, info in win_res.items():
            results[sample] = results.get(sample, '') + ':' + str(info['cluster'])

    unique, new_id, final = {}, 1, defaultdict(lambda: defaultdict(dict))
    for key, value in results.items():
        if value not in unique:
            unique[value] = new_id
            new_id += 1
        final[args.taxon][args.segment].setdefault(unique[value], []).append(key)

    # Write clusters to JSON and FASTA
    prefix = f'{args.taxon}-{args.segment}'
    with open(f'{args.outdir}/{prefix}.clusters.json', 'w') as f:
        json.dump(final, f, indent=4)

    for cluster, ids in final[args.taxon][args.segment].items():
        if len(ids) > 1:
            var = 'multi'
        else:
            var = 'single'
        prefix_fa = f'{prefix}-{cluster}.{var}'
        fas = []
        for id in ids:
            fas.append(f'>{id}\n{seqs[id]}')
        with gzip.open(f'{args.outdir}/{prefix_fa}.fa.gz', 'wt') as f:
            f.write('\n'.join(fas) + '\n')

    end = f"""
    Total clusters {new_id - 1}
    Results saved to {args.outdir}
    """
    print(textwrap.dedent(end))

if __name__ == "__main__":
    main()