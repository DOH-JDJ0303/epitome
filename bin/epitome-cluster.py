#!/usr/bin/env python3

# epitome-cluster.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import sourmash
import screed
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
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



# ----- FUNCTIONS ----- #
def loadSeqs(fasta, ksize, scaled, window_size):
    seqs, windows = {}, {}
    for rec in screed.open(fasta):
        seqs[rec.name] = rec.sequence
    print(f'{datetime.datetime.now()}: Loaded {len(seqs)} sequences from {fasta}')
    for k, seq in seqs.items():
        j, seq_size = 0, len(seq)
        win_size = min(window_size, seq_size)
        for i in range(win_size - 1, seq_size, win_size):
            mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(seq[j:i-1], True)
            window = f'{j}-{i}'
            if window not in windows:
                windows[window] = {}
            windows[window][k] = {'mh': mh}
            j = i
    print(f'{datetime.datetime.now()}: Split sequences into {len(windows.items())}')
    return seqs, windows

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
    return squareform(mat), ids

def createClusters(data, threshold, start, stage, outdir):
    mat, ids = computeMatrix(data)
    Z = linkage(mat, method='complete')
    clusters = cut_tree(Z, height=threshold).flatten()
    max_cluster = start
    for i, id in enumerate(ids):
        c = int(clusters[i]) + start + 1
        data[id]['cluster'] = c
        max_cluster = max(max_cluster, c)
    try:
        os.makedirs(f'{outdir}/windows/', exist_ok=True)
        plt.figure(figsize=(10, 5))
        dendrogram(Z, labels=ids, orientation='left')
        plt.axvline(x=threshold, color='r', linestyle='--')
        plt.title(f'{stage}')
        plt.xlabel('Distance')
        plt.ylabel('Samples')
        plt.savefig(f'{outdir}/windows/{stage}.png')
        plt.close()
    except Exception as e:
        print(f'Plot not made: {e}')
    return data, max_cluster

def assignClusters(data1, data2, data3, threshold):
    for k1, v1 in data1.items():
        for k2, v2 in data2.items():
            if v1['mh'].containment_ani(v2['mh']).dist < threshold:
                v1['cluster'] = v2['cluster']
                data3[k1] = v1
    return data3

def clusterSeqs(data, max_cluster, threshold, round, start, outdir):
    if start != 0:
        print(f'{datetime.datetime.now()}: Round {round} - Continuing from cluster {start}')
    if len(data) > 1:
        main, rem = mainAndRemainder(data, max_cluster, round)
        print(f'{datetime.datetime.now()}: Round {round} - Clustering {len(main)} sequences')
        clusters, last_cluster = createClusters(main, threshold, start, f'main_{round}', outdir)
        print(f'{datetime.datetime.now()}: Round {round} - Clustered into {last_cluster - start} clusters')
        if rem:
            reps = selectReps(clusters)
            print(f'{datetime.datetime.now()}: Round {round} - Assigning {len(rem)} sequences')
            clusters = assignClusters(rem, reps, clusters, threshold)
        loose = {k: rem[k] for k in rem if k not in clusters}
    else:
        loose, last_cluster = {}, start + 1
        for v in data.values(): v['cluster'] = last_cluster
        clusters = data
    return clusters, loose, last_cluster

def dict2Json(data, taxon, segment, filename):
    """
    Restructure flat list of data entries into nested JSON:
      taxon -> segment -> sequence_id -> metadata fields
    """
    print(f'{datetime.datetime.now()}: Writing clusters to {filename}')
    nested = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for entry in data:
        nested[taxon][segment][entry['seq']] = { 'cluster': entry['cluster'] }

    with open(filename, 'w') as f:
        json.dump(nested, f, indent=4)

def main():
    version = "1.1"

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, required=True, help="Path to FASTA file.")
    parser.add_argument("--taxon", default='null', help="Taxon name or 'version'")
    parser.add_argument("--segment", default='null', help="Segment name")
    parser.add_argument("--dist", default=0.05, type=float, help="Distance threshold (1 - ANI/100)")
    parser.add_argument("--max_cluster", default=1000, type=int, help="Max sequences in initial clustering round")
    parser.add_argument("--ksize", default=31, type=int, help="K-mer size for Sourmash")
    parser.add_argument("--scaled", default=1000, type=int, help="Scaled value for Sourmash")
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
    seqs, windows = loadSeqs(args.fasta, args.ksize, args.scaled, args.window_size)
    # Create clusters
    results = {}
    for window, data in windows.items():
        if len(data) == len(seqs):
            print(f'{datetime.datetime.now()}: Clustering window {window}')
            start, round, win_res = 0, 1, {}
            while data:
                clusters, data, start = clusterSeqs(data, args.max_cluster, args.dist, f'{window}_{round}', start, args.outdir)
                win_res.update(clusters)
                round += 1
            for sample, info in win_res.items():
                results[sample] = results.get(sample, '') + ':' + str(info['cluster'])

    unique, new_id = {}, 1
    for v in results.values():
        if v not in unique:
            unique[v] = new_id
            new_id += 1
    final = [{'seq': k, 'cluster': unique[v]} for k, v in results.items()]
    dict2Json(final, args.taxon, args.segment, f'{args.outdir}/clusters.json')

    end = f"""
    Total clusters {new_id - 1}
    Results saved to {args.outdir}
    """
    print(textwrap.dedent(end))

if __name__ == "__main__":
    main()