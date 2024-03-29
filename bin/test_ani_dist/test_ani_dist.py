#!/usr/bin/env python

import os
import argparse
import subprocess
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from typing import List
from Bio import SeqIO
import math
import numpy as np


parser = ArgumentParser()
parser.add_argument("--input", dest="in_dir",  required=True, help="/path/to/input/seq/dir")
parser.add_argument("--min", dest="dist_min",  required=True, type=float, help="min ani difference to test")
parser.add_argument("--max", dest="dist_max",  required=True, type=float, help="max ani difference to test")
parser.add_argument("--step", dest="dist_step",  required=True, type=float, help="steps between min and max values of ANI distance to test")
parser.add_argument("--rep", dest="replicates",  required=True, type=int, help="number of replicates")
parser.add_argument("--outdir", dest="outdir",  help="number of replicates")


args = parser.parse_args() 
in_dir = args.in_dir
dist_min = args.dist_min
dist_max = args.dist_max
dist_step = args.dist_step
replicates = args.replicates
outdir = args.outdir


def num_point_mutations(filename:str, perc: float) -> int:
    if len(list(SeqIO.parse(filename, "fasta"))) != 1:
        raise Exception(f'{filename} is a multifasta, please split records \
                        in to individual fasta files')
    nucs = math.ceil(len(list(SeqIO.parse(filename, "fasta"))[0].seq)*perc)
    return nucs


def generate_outname(filename:str, replicate: int, ani: int, outdir:str) -> str:
    basename ='.'.join(filename.split(".")[:-1]).split("/")[-1]
    print(basename)
    ext = "." + filename.split(".")[-1]
    rename = basename + ".rep_" + str(replicate) + ".ani_" + str(ani) 
    outname = rename + ext
    if outdir:
        outname = os.path.join(outdir,outname)
    else:
        outname = os.path.join("./",outname)
    return outname


if __name__ == "__main__":
    # Generate test dataset
    if outdir:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    for filename in os.listdir(in_dir):
        filename = os.path.join(in_dir,filename)   
        for ani_dist in np.arange(dist_min, dist_max, dist_step):
            ani = str(1-ani_dist).replace(".", "-")
            point_muts = num_point_mutations(filename, ani_dist)
            for replicate in range(replicates):
                outname = generate_outname(filename, replicate, ani, outdir)
                subprocess.run(f'msbar {filename} \
                     -count {point_muts} -point 1 \
                     -block 0 -codon 0 \
                     -outseq {outname}', shell=True)
                print(outname)
                