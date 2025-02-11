#!/usr/bin/env python3

# epitome-condense.py
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
parser.add_argument("--clusters", dest="clusters", type=str, help="Path to the cluster info.")
parser.add_argument("--dist_threshold", dest="dist_threshold",  default=0.05, type=float, help="Distance thresholds used to create clusters (1-%ANI/100)(default: 0.05)")
parser.add_argument("--ksize", dest="ksize",  default=31, type=int, help="kmer size used by Sourmash")
parser.add_argument("--outdir", dest="outdir",  default='./', type=str, help="Name of outputfile.")
args = parser.parse_args()

#----- FUNCTIONS -----#
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
    
    return mat, ids, dist_long

# Function for assigning clusters using a distance threshold
def createClusters(data, threshold, start, stage):
    result = data
    # Compute pairwise ANI
    mat, ids, dist_long = computeMatrix(data)
   
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
        
    return result, max_cluster, dist_long
# Function for condensing sequences
def condense(data):
    # group by cluster
    clusters = {}
    for key, value in data.items():
        if value['cluster'] in clusters:
            clusters[value['cluster']].append(key)
        else:
            clusters[value['cluster']] = [ key ]
    # select references for each cluster
    refs = []
    for key, value in clusters.items():
        if len(value) > 1:
            max_n      = 0
            max_len    = 0
            select_ref = []
            for id in value:
                if int(data[id]['n']) > max_n:
                    select_ref = [ id ]
                    max_n = int(data[id]['n'])
                elif int(data[id]['n']) == max_n:
                    select_ref.append(id)
            if len(select_ref) > 1:
                select_ref2 = []
                for id in select_ref:
                    if int(data[id]['length']) > max_len:
                        select_ref2 = [ id ]
                        max_len = int(data[id]['length'])
                    elif int(data[id]['length'])== max_len:
                        select_ref2.append(id)
                select_ref = select_ref2
            refs.append(select_ref[0])
        else:
            select_ref = [ value[0] ]
        refs.append(select_ref[0])
        for id in value:
            if id in refs:
                data[id]['condensed'] = ''
            else:
                data[id]['condensed'] = select_ref[0]
    return refs, data

# Function for finding the next closests sequence based on ANI
def findNeighbor(inclusions, exclusions, dists, data):
    sub_dist = []
    for pair in dists:
        if pair['id1'] not in exclusions and pair['id2'] not in exclusions:
                sub_dist.append(pair)
    for id in inclusions:
        neighbor      = []
        neighbor_dist = 1
        for pair in sub_dist:
            if pair['id1'] == id:
                other_id = pair['id2']
                dist = pair['dist']
            elif pair['id2'] == id:
                other_id = pair['id1']
                dist = pair['dist']
            else:
                continue
            if dist < neighbor_dist:
                neighbor = [ other_id ]
                neighbor_dist = dist
            elif dist == neighbor_dist:
                neighbor.append(other_id)
        neighbor = list(set(neighbor))
        neighbor_ani = round(100*(1-neighbor_dist))
        data[id]['neighbor_ani'] = neighbor_ani
        if neighbor_ani == 0:
            neighbor = []
            
        data[id]['neighbor'] = neighbor
    return data

# Function for writing a dictionary to CSV
def dictToCsv(dicts, filename):
    # Get all unique keys
    keys = ["taxon","segment","cluster","ref","condensed","n","length","neighbor_ani","neighbor"]

    # Write to CSV
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys, quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(dicts)

#----- MAIN -----#
# create a dictionary of sequence supporting each reference
seqs_per_cluster = {}
with open(args.clusters, newline='') as csvfile: 
    csvreader = csv.reader(csvfile)
    header = next(csvreader)
    cluster_index = header.index('cluster')
    for row in csvreader:
        cluster = row[cluster_index]
        if cluster in seqs_per_cluster:
            seqs_per_cluster[cluster] += 1
        else:
            seqs_per_cluster[cluster] = 1

# Create dictionary of hashes
minhashes = {}
ids = []
for record in screed.open(args.fasta):
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1000)
    mh.add_sequence(record.sequence, True)
    cluster = record.name.split('-')[2]
    minhashes[record.name] = { 'mh': mh, 'length': len(record.sequence), 'n': seqs_per_cluster[ cluster ] }
    ids.append(record.name)
ids = list(set(ids))
print(f'{datetime.datetime.now()}: Loaded {len(minhashes.items())} sequences from {args.fasta}')

# Cluster sequences
print(f'{datetime.datetime.now()}: Creating clusters')
data, max_cluster, dist_long = createClusters(minhashes, args.dist_threshold, 0, 'condense')
print(f'{datetime.datetime.now()}: {len(data.items())} clusters created.')
if len(data.items()) == max_cluster:
    print(f'{datetime.datetime.now()}: Nothing to condense!')
    refs = ids
    for key, value in data.items():
        value['condensed'] = ''
else:
    print(f'{datetime.datetime.now()}: Condensing references.')
    refs, data = condense(data)
    print

# Create sets of refs and condensed
refs = list(set(refs))
condensed = list(set([ id for id in ids if id not in refs ]))
print(f'{datetime.datetime.now()}: Finding next closest neighbors')
# Find next closest neighbor(s)
data = findNeighbor(refs, condensed, dist_long, data)
if condensed:
    data = findNeighbor(condensed, [], dist_long, data)
# Write results to CSV
print(f'{datetime.datetime.now()}: Writing file.')
results = []
for key, value in data.items():
    value['ref'] = key
    ref_name_split = key.split('-')
    value['taxon'] = ref_name_split[0]
    value['segment'] = ref_name_split[1]
    value['cluster'] = ref_name_split[2]
    del value['mh']
    results.append(value)

dictToCsv(results, f'{args.outdir}/condensed.csv')

# End
end = f"""
Results saved to {args.outdir}
"""
print(end)