#!/usr/bin/env python3
# epitome_consensus.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import gzip
import logging
import os
import random
import sys
from typing import Iterator, Tuple, TextIO

import numpy as np

from epitome_utils import sanitize_filename, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

SYMBOLS = ['-', 'A', 'T', 'C', 'G', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
IDX = {c: i for i, c in enumerate(SYMBOLS)}

# Which IUPAC codes contribute to A/T/C/G tallies
A_GROUP = ['A', 'R', 'W', 'M', 'D', 'H', 'V', 'N']
T_GROUP = ['T', 'Y', 'W', 'K', 'B', 'D', 'H', 'N']
C_GROUP = ['C', 'Y', 'S', 'M', 'B', 'H', 'V', 'N']
G_GROUP = ['G', 'R', 'S', 'K', 'B', 'D', 'V', 'N']

# -------------------------------
#  FUNCTIONS
# -------------------------------

def open_maybe_gz(path: str) -> TextIO:
    """Open a file path in text mode, supporting .gz transparently.

    Args:
        path: File path, optionally ending in '.gz'.

    Returns:
        A text-mode file handle (TextIO).
    """
    return gzip.open(path, 'rt', encoding='utf-8') if path.endswith('.gz') else open(path, 'r', encoding='utf-8')

def fasta_sequences(path: str) -> Iterator[str]:
    """Stream sequences from FASTA(.gz) as single uppercase strings (no headers).

    Args:
        path: Path to FASTA file (optionally gzipped).

    Yields:
        Uppercased sequence strings with whitespace removed.
    """
    with open_maybe_gz(path) as fh:
        seq_chunks = []
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if seq_chunks:
                    yield ''.join(seq_chunks).upper().replace(' ', '').replace('\t', '')
                    seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if seq_chunks:
            yield ''.join(seq_chunks).upper().replace(' ', '').replace('\t', '')

def build_lut() -> np.ndarray:
    """Build lookup table mapping ASCII code -> index in SYMBOLS (or -1 if illegal).

    Returns:
        NumPy array of shape (256,) dtype=int16 mapping ord(char) to SYMBOLS index.
    """
    lut = np.full(256, -1, dtype=np.int16)
    for c, i in IDX.items():
        lut[ord(c)] = i
    return lut

def tally_alignment(aln_path: str, rng: np.random.Generator) -> Tuple[np.ndarray, int]:
    """
    Stream over sequences in an alignment and tally per-site symbol counts.

    Args:
        aln_path: Path to alignment FASTA (optionally gzipped).
        rng: NumPy random Generator (unused here but kept for interface parity).

    Returns:
        counts: NumPy array of shape (L, 16) with per-site tallies for SYMBOLS.
        n_seq: Total number of sequences processed.
    """
    lut = build_lut()
    seq_iter = fasta_sequences(aln_path)
    try:
        first = next(seq_iter)
    except StopIteration:
        raise SystemExit("Error: empty alignment")

    L = len(first)
    if L == 0:
        raise SystemExit("Error: first sequence is empty")

    counts = np.zeros((L, len(SYMBOLS)), dtype=np.uint32)

    def process_seq(seq: str) -> None:
        nonlocal counts
        if len(seq) != L:
            raise SystemExit(f"Error: sequences have differing lengths (expected {L}, got {len(seq)})")
        arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        idx = lut[arr]
        if (idx < 0).any():
            bad = sorted(set(chr(c) for c in arr[idx < 0]))
            raise SystemExit(f"Error: the alignment contains illegal characters: {''.join(bad)}")
        counts[np.arange(L), idx] += 1

    process_seq(first)
    n_seq = 1
    for s in seq_iter:
        process_seq(s)
        n_seq += 1
        if n_seq % 1000 == 0:
            LOGGER.debug("Processed %d sequences...", n_seq)

    return counts, n_seq

def compute_site_calls(counts: np.ndarray, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Compute per-site max counts and consensus calls among -,A,T,C,G with random tie-breaking.

    Args:
        counts: NumPy array (L, 16) of per-symbol counts for each site.
        rng: NumPy random Generator used for tie-breaking.

    Returns:
        max_counts: NumPy array (L,) of maximum counts at each site among -,A,T,C,G.
        calls: NumPy uint8 array (L,) of ASCII codes for chosen calls.
    """
    # Grab per-symbol columns
    col = {c: counts[:, IDX[c]] for c in SYMBOLS}
    dash = col['-']
    A = sum(col[c] for c in A_GROUP)
    T = sum(col[c] for c in T_GROUP)
    C = sum(col[c] for c in C_GROUP)
    G = sum(col[c] for c in G_GROUP)

    # Per-site maximum among dash, A,T,C,G
    stacked = np.stack([dash, A, T, C, G], axis=1)  # shape (L, 5)
    max_counts = stacked.max(axis=1)

    # Tie-breaking: choose randomly among dash/A/T/C/G with the max
    mask = (stacked == max_counts[:, None])  # shape (L, 5)
    # Build choices array of same shape with codepoints: ['-','A','T','C','G'] -> bytes
    choices = np.array([ord('-'), ord('A'), ord('T'), ord('C'), ord('G')], dtype=np.int16)

    # For each row, pick one of the true positions uniformly
    L = counts.shape[0]
    calls = np.empty(L, dtype=np.uint8)
    for i in range(L):
        winners = np.flatnonzero(mask[i])
        pick = rng.integers(0, len(winners))
        calls[i] = choices[winners[pick]]
    return max_counts, calls

def write_site_list(site_csv_path: str, counts: np.ndarray, max_counts: np.ndarray, calls: np.ndarray) -> None:
    """Write per-site symbol tallies, max counts, and calls to CSV.

    Args:
        site_csv_path: Output CSV file path.
        counts: NumPy array (L, 16) of per-symbol counts.
        max_counts: NumPy array (L,) of maximum counts.
        calls: NumPy uint8 array (L,) of ASCII codes for calls.

    Returns:
        None
    """
    header = "site,-,A,T,C,G,R,Y,S,W,K,M,B,D,H,V,N,max,call\n"
    with open(site_csv_path, 'w', encoding='utf-8') as f:
        f.write(header)
        L = counts.shape[0]
        for i in range(L):
            row = [str(i+1)]
            row.extend(str(int(counts[i, IDX[c]])) for c in SYMBOLS)
            row.append(str(int(max_counts[i])))
            row.append(chr(int(calls[i])))
            f.write(','.join(row) + '\n')

def write_consensus_fasta(name: str, calls: np.ndarray, out_path: str) -> None:
    """Write a gzipped FASTA containing the consensus sequence (gaps removed).

    Args:
        name: FASTA header (sequence name).
        calls: NumPy uint8 array (L,) of ASCII codes for calls.
        out_path: Output .fa.gz path.

    Returns:
        None
    """
    seq = bytes([c for c in calls if c != ord('-')]).decode('ascii')
    with gzip.open(out_path, 'wt', encoding='utf-8') as f:
        f.write(f">{name}\n{seq}\n")

# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point: compute consensus from an alignment and write CSV/FASTA outputs.

    Parses arguments, tallies per-site counts, computes consensus calls with
    random tie-breaking among -,A,T,C,G, then writes a site-wise CSV and a
    gzipped FASTA consensus.

    Returns:
        None
    """
    version = "2.0"

    parser = argparse.ArgumentParser(description="Consensus caller for multiple sequence alignments.")
    parser.add_argument("--prefix", help="Prefix used to name the file and output fasta header.")
    parser.add_argument("--aln",  help="Alignment FASTA (.fa/.fasta[.gz]).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed.")
    parser.add_argument("--outdir", default='.', help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    rng = np.random.default_rng(args.seed)

    # Tally alignment
    counts, n_seq = tally_alignment(args.aln, rng)
    L = counts.shape[0]
    LOGGER.info("Alignment Length: %d", L)
    LOGGER.info("Number of Sequences: %d", n_seq)

    # Compute calls
    max_counts, calls = compute_site_calls(counts, rng)

    # Write site list CSV
    site_csv_path = os.path.join(args.outdir, f"{args.prefix}.pileup.csv")
    write_site_list(site_csv_path, counts, max_counts, calls)
    LOGGER.info("Wrote site-wise counts: %s", site_csv_path)

    # Write consensus FASTA
    consensus_path = os.path.join(args.outdir, f"{args.prefix}.fa.gz")
    write_consensus_fasta(args.prefix, calls, consensus_path)
    LOGGER.info("Wrote consensus FASTA: %s", consensus_path)

if __name__ == "__main__":
    main()
