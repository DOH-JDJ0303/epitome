#!/usr/bin/env python3

# epitome-cluster.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

version = 1.0
start = f"""

epitome-cluster.py v{version}

Written by Jared Johnson
jared.johnson@doh.wa.gov
"""
print(start)

import sourmash
import screed
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
import random
import csv
import random
import datetime
import argparse
import concurrent.futures
import matplotlib.pyplot as plt
import itertools



#----- ARGS -----#
parser = argparse.ArgumentParser()
parser.add_argument("--fasta", dest="fasta", type=str, help="Path to FASTA file.")
parser.add_argument("--dist_threshold", dest="dist_threshold",  default=0.05, type=float, help="Distance thresholds used to create clusters (1-%ANI/100)(default: 0.05)")
parser.add_argument("--max_cluster", dest="max_cluster",  default=1000, type=int, help="Maximum number of sequences to compare in the initial clustering round")
parser.add_argument("--ksize", dest="ksize",  default=31, type=int, help="kmer size used by Sourmash")
parser.add_argument("--scaled", dest="scaled",  default=1000, type=int, help="scaled value used by Sourmash")
parser.add_argument("--outdir", dest="outdir",  default='./', type=str, help="Name of outputfile.")
args = parser.parse_args()

#----- FUNCTIONS -----#
# Function for splitting a dictionary into a random subset and its remainder
def mainAndRemainder(data, threshold, round):
    if len(data.items()) > threshold:
        # Get a random subset of keys
        keys = list(data.keys())
        random_keys = random.sample(keys, threshold)
        
        # Create the random subset dictionary
        subset = {key: data[key] for key in random_keys}
        for key, value in subset.items():
            value['source'] = f'main_{round}'
        
        # Create the remainder dictionary
        remainder = {key: data[key] for key in data if key not in random_keys}
        for key, value in remainder.items():
            value['source'] = f'assigned_{round}'
    else:
        for key, value in data.items():
            value['source'] = f'main_{round}'
        subset    = data
        remainder = {}
    
    return subset, remainder

# Function for selecting representatives of each cluster
def selectReps(data):
    result = {}
    clusters = []
    for key, value in data.items():
        cluster = value['cluster']
        if cluster not in clusters:
            result[key] = value
            clusters.append(cluster)
    return result

# Function to compute the ANI distance for a pair of items
def computeDistance(pair, data):
    id1, id2 = pair
    x = data[id1]['mh'].containment_ani(data[id2]['mh'])
    return { 'id1': id1, 'id2': id2, 'dist': x.dist }

# Function to parallelize the computation of the matrix
def computeMatrix(data):
    dist_long = []

    # Get pairwise comparisons
    ids = list(set([ key for key, value in data.items() ]))
    ids.sort()
    pairs = list(itertools.combinations(ids, 2))
    for pair in pairs:
        dist_long.append(computeDistance(pair, data))

    # Extract unique samples 
    ids = list(set(d['id1'] for d in dist_long).union(d['id2'] for d in dist_long)) 
    ids.sort()
    # Initialize the matrix with zeros 
    mat = [[0 for _ in ids] for _ in ids]

    # Fill the matrix with distances from the long format data 
    for entry in dist_long:
        i = ids.index(entry['id1'])
        j = ids.index(entry['id2'])
        mat[i][j] = entry['dist'] 
        mat[j][i] = entry['dist']
    
    mat = np.array(mat)
    np.fill_diagonal(mat, 0)
    mat = squareform(mat)
    
    return mat, ids

# Function for assigning clusters using a distance threshold
def createClusters(data, threshold, start, stage):
    result = data
    # Compute pairwise ANI
    mat, ids = computeMatrix(data)
   
    # Compute complete-linkage of pairwise ANI
    Z = linkage(mat, method='complete')

    # Cut dendrogram at distance threshold
    max_cluster = start
    clusters = cut_tree(Z, height=threshold).flatten()
    for index, id in enumerate(ids):
        cluster = int(clusters[index]) + int(start) + 1
        result[id]['cluster'] = cluster
        # Calculate the largest cluster number
        if cluster > int(max_cluster):
            max_cluster = cluster
    
    try:
        # Plot the dendrogram
        plt.figure(figsize=(10, 5))
        dendro = dendrogram(Z, labels=ids, orientation='left')

        # Draw a horizontal line at the threshold
        plt.axvline(x=threshold, color='r', linestyle='--')
        plt.title(f'{stage}')
        plt.xlabel('Samples')
        plt.ylabel('Distance')

        # Save the plot to a file
        plt.savefig(f'{args.outdir}/{stage}.png', format='png')
        plt.close()
    except:
        print('Plot not made!')
        
    return result, max_cluster

# Function for assigning clusters based on existing clusters
def assignClusters(data1, data2, data3, threshold):
    for key1, value1 in data1.items():
        for key2, value2 in data2.items():
            x = value1['mh'].containment_ani(value2['mh'])
            if x.dist < threshold:
                value1['cluster'] = value2['cluster']
                data3[key1] = value1
    return data3

# Function for peforming a clustering round
def clusterSeqs(data, max_cluster, threshold, round, start):
    if start != 0:
        print(f'{datetime.datetime.now()}: Round {round} - Continuing from cluster {start}')
    if len(data.items()) > 1:
        # Subset hashes based on the max allowed initial cluster dataset size
        main, remainder = mainAndRemainder(data, max_cluster, round)
        # Clustering Step 1 (Main)
        print(f'{datetime.datetime.now()}: Round {round} - Clustering {len(main.items())} sequences')
        clusters, last_cluster = createClusters(main, threshold, start, f'main_{round}')
        print(f'{datetime.datetime.now()}: Round {round} - Clustered {len(clusters.items())} sequences into {last_cluster - start} clusters')

        # Check if any further clustering is required
        if len(remainder.items()) > 0:
            # Assign clusters to the remaining sequences using representatives from the main clustering
            reps = selectReps(clusters)
            print(f'{datetime.datetime.now()}: Round {round} - Assigning {len(remainder.items())} sequences')
            pre_assigned = len(clusters.items())
            clusters = assignClusters(remainder, reps, clusters, threshold)
            post_assigned = len(clusters.items())
            print(f'{datetime.datetime.now()}: Round {round} - Assigned {post_assigned - pre_assigned} sequences')
        
        # Identify loose-ends
        looseends = {}
        looseends = {key: remainder[key] for key in remainder if key not in clusters.keys()}
        print(f'{datetime.datetime.now()}: Round {round} - Done; {len(looseends.items())} sequences remain')
    else:
        looseends = {}
        last_cluster = start + 1
        for key, value in data.items():
            value['cluster'] = last_cluster
        clusters = data

    return clusters, looseends, last_cluster

# Function for writing a dictionary to CSV
def dictToCsv(dicts, filename):
    # Get all unique keys
    keys = set()
    for d in dicts:
        keys.update(d.keys())
    keys = list(keys)

    # Write to CSV
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(dicts)

#----- MAIN -----#
# Create dictionary of hashes
minhashes = {}
for record in screed.open(args.fasta):
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)
    mh.add_sequence(record.sequence, True)
    minhashes[record.name] = { 'mh': mh }
print(f'{datetime.datetime.now()}: Loaded {len(minhashes.items())} sequences from {args.fasta}')

data = minhashes
start = 0
round = 1
results = {}
while len(data.items()) > 0:
    clusters, data, start = clusterSeqs(data, args.max_cluster, args.dist_threshold, round, start)
    results.update(clusters)
    round += 1

# Write results to CSV
print(f'{datetime.datetime.now()}: Writing file.')
final = []
for key, value in results.items():
    value['seq'] = key
    del value['mh']
    final.append(value)

dictToCsv(final, f'{args.outdir}/clusters.csv')

end = f"""

Total clusters {start}
Results saved to {args.outdir}

"""
print(end)