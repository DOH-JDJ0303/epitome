#!/usr/bin/env python3

# epitome_inputs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import json
import logging
import gzip
import os
import sys
import re
from datetime import datetime
from collections import OrderedDict
from typing import List, Dict, Any, Tuple, Optional
import screed
import sourmash
import numpy as np
from sklearn.cluster import DBSCAN

from epitome_utils import sanitize_filename, detect_and_read, normalize_keys, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

STRING_KEYS = ['accession', 'segment', 'tax_id']
YEAR_KEYS   = ['collection_date']
MISSING_DATA = ["na", "n/a", "null", "none", "nan"]

# -------------------------------
#  FUNCTIONS
# -------------------------------

def _sanitize_str(s: str) -> str:
    """Keep letters, numbers, underscores, dashes, and dots; replace others with '_'."""
    s0 = s
    s = s.strip().replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9._-]', "_", s)
    s = re.sub(r'__+', "_", s)
    s = s.lower()
    if s != s0:
        LOGGER.debug("Sanitized string '%s' -> '%s'", s0, s)
    return s

def _extract_year(rec: Dict[str, Any]) -> None:
    """Extract a 4-digit year from YEAR_KEYS fields and replace the value with an int if found."""
    year_re = re.compile(r"^\d{4}")
    for k in YEAR_KEYS:
        v = rec.get(k, None)
        if v:
            m = year_re.search(str(v))
            if m:
                old = rec[k]
                rec[k] = int(m.group())
                if rec[k] != old:
                    LOGGER.debug("Extracted year for key '%s': %s -> %s", k, old, rec[k])

def _coerce_to_str(rec: Dict[str, Any]) -> None:
    """Convert certain keys in a record to strings."""
    for k in STRING_KEYS:
        if k in rec and rec[k] is not None and rec[k] != "" and not isinstance(rec[k], str):
            old = rec[k]
            rec[k] = str(rec[k])
            LOGGER.debug("Coerced key '%s' to str: %r -> %r", k, old, rec[k])

def _update_missing(rec: Dict[str, Any]) -> None:
    """Replace certain missing value markers with None."""
    for k, v in list(rec.items()):
        if isinstance(v, str) and v.lower() in MISSING_DATA:
            LOGGER.debug("Setting missing value for key '%s': %r -> None", k, v)
            rec[k] = None

def merge_records(metadata, fasta):
    LOGGER.info("Merging metadata with FASTA sequences")
    merged = []
    required_keys = ["accession"]
    accs = set()

    for src_path, recs in metadata:
        LOGGER.info("Source %s: %d metadata records", src_path, len(recs))
        for rec in recs:
            rec = normalize_keys(rec)
            _coerce_to_str(rec)
            _extract_year(rec)
            _update_missing(rec)

            if not all(k in rec for k in required_keys):
                raise ValueError(f'One or more required key ({required_keys}) is missing from record {rec}')

            acc = rec['accession']
            if acc not in fasta:
                LOGGER.warning("%s not found in supplied assembly files", acc)
                continue

            rec['sequence'] = fasta[acc]

            if acc in accs:
                raise ValueError(f'Multiple metadata entries for accession {acc}')
            accs.add(acc)

            merged.append(rec)

    LOGGER.info("Merged %d records with sequences (from %d sources)", len(merged), len(metadata))
    return merged

def drop_na(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove keys that are None for all records."""
    if not data:
        return data
    n_rec = len(data)
    na_counts: Dict[str, int] = {}
    rm_set = set()

    for rec in data:
        for k, v in rec.items():
            if v is None:
                na_counts[k] = na_counts.get(k, 0) + 1
                if na_counts[k] == n_rec:
                    rm_set.add(k)

    if rm_set:
        LOGGER.info("Dropping all-NA keys across %d records: %s", n_rec, ", ".join(sorted(rm_set)))
    for rec in data:
        for k in rm_set:
            if k in rec:
                rec.pop(k)
    return data

def group_segments(
    metadata: List[Dict[str, Any]],
    ksize: int,
    scaled: int,
    dist: float = 0.5,
    keep_singletons: bool = False,
) -> Dict[str, List[Dict[str, Any]]]:
    LOGGER.info("Grouping segments (ksize=%d, scaled=%d, dist=%.3f, keep_singletons=%s)",
                ksize, scaled, dist, keep_singletons)

    def _concat_seq(recs: List[Dict[str, Any]]) -> str:
        return "".join(rec.get("sequence", "") for rec in recs)

    def _mk_mh(seq: str, track_abundance: bool) -> sourmash.MinHash:
        mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=track_abundance)
        mh.add_sequence(seq, force=True)
        return mh

    def _cluster_with_noise(D: np.ndarray, ids: List[str], eps: float) -> Dict[int, List[str]]:
        db = DBSCAN(eps=float(eps), min_samples=1, metric="precomputed").fit(D)
        labels = db.labels_
        clusters: Dict[int, List[str]] = {}
        max_cluster = max(labels) if len(labels) else -1
        for i, lbl in enumerate(labels):
            if lbl == -1:
                max_cluster += 1
                lbl = max_cluster
            clusters.setdefault(int(lbl), []).append(ids[i])
        LOGGER.debug("DBSCAN produced %d clusters on %d ids", len(clusters), len(ids))
        return clusters

    def compute_matrix(
        mh1: Dict[str, sourmash.MinHash],
        mh2: Dict[str, sourmash.MinHash],
    ) -> Tuple[np.ndarray, List[str]]:
        ids = sorted(set(mh1.keys()) & set(mh2.keys()))
        n = len(ids)
        D = np.zeros((n, n), dtype=float)

        for i in range(n):
            ai = ids[i]
            for j in range(i + 1, n):
                aj = ids[j]
                d = 1.0 - mh1[ai].max_containment(mh2[aj])
                D[i, j] = D[j, i] = d

        np.fill_diagonal(D, 0.0)
        LOGGER.debug(f"Computed distance matrix of shape {D.shape} for {n} ids")
        return D, ids

    assigned: Dict[str, List[Dict[str, Any]]] = {}
    unassigned: Dict[str, Dict[str, Any]] = {}
    for rec in metadata:
        seg = rec.get("segment")
        if seg:
            seg = _sanitize_str(rec.get("segment"))
            assigned.setdefault(seg, []).append(rec)
        else:
            unassigned[rec["accession"]] = rec

    LOGGER.info(f"Assigned segments: {len(assigned)} groups; Unassigned records: {len(unassigned)} ({round(100*len(unassigned)/len(metadata), 1)}%)")

    group_full: Dict[str, sourmash.MinHash] = {}
    group_core: Dict[str, sourmash.MinHash] = {}

    for seg, records in assigned.items():
        seq = _concat_seq(records)
        if not seq:
            LOGGER.debug(f"Skipping empty sequence group '{seg}'")
            continue

        mh_full = _mk_mh(seq, track_abundance=True)
        group_full[seg] = mh_full

        mh_core = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=False)
        keep_all = (len(records) == 1)
        for h, count in mh_full.hashes.items():
            if keep_all or count > 1:
                mh_core.add_hash(h)
        group_core[seg] = mh_core

    D, ids = compute_matrix(group_full, group_core)
    clusters = _cluster_with_noise(D, ids, dist)

    final: Dict[str, List[Dict[str, Any]]] = {}
    for seg_list in clusters.values():
        if len(seg_list) == 1:
            seg = seg_list[0]
            final[seg] = assigned.get(seg, []).copy()
            continue

        common_name = max(
            seg_list,
            key=lambda s: (len(assigned.get(s, [])), -ord(s[0]) if s else 0),
        )
        merged: List[Dict[str, Any]] = []
        for s in seg_list:
            merged.extend(assigned.get(s, []))
        final[common_name] = merged

    final_mh: Dict[str, sourmash.MinHash] = {}
    for seg, recs in final.items():
        seq = _concat_seq(recs)
        if not seq:
            continue
        final_mh[seg] = _mk_mh(seq, track_abundance=False)

    leftovers: Dict[str, Tuple[Dict[str, Any], Optional[sourmash.MinHash]]] = {}
    for acc, rec in unassigned.items():
        seq = rec.get("sequence", "")
        if not seq or not final_mh:
            leftovers[acc] = (rec, None)
            continue

        mh = _mk_mh(seq, track_abundance=False)

        best_seg = None
        best_d = 1.0
        for seg, mh2 in final_mh.items():
            d = 1.0 - float(mh.contained_by(mh2))
            if d < best_d:
                best_d = d
                best_seg = seg

        if best_seg is not None and best_d < dist:
            final[best_seg].append(rec)
        else:
            leftovers[acc] = (rec, mh)

    LOGGER.info(f"Post-assignment groups: {len(final)}; leftovers: {len(leftovers)} ({round(100*len(leftovers)/len(metadata), 1)}%)")

    if not leftovers:
        return final
    if not keep_singletons:
        LOGGER.info("Discarding singletons")
        return final
    else:
        LOGGER.info("Naming singletons")
        singletons = {acc: mh for acc, (_, mh) in leftovers.items() if mh is not None}
        if len(singletons) == 1:
            (acc, _mh), = leftovers.items()
            final.setdefault("segment_group_0", []).append(leftovers[acc][0])
        elif len(singletons) > 1:
            D2, ids2 = compute_matrix(singletons, singletons)
            clusters2 = _cluster_with_noise(D2, ids2, dist)
            for lbl_ids in clusters2.values():
                grp = f"segment_group_{int(list(clusters2.keys())[0])}"
                for acc in lbl_ids:
                    final.setdefault(grp, []).append(leftovers[acc][0])
        else:
            for acc, (rec, _mh) in leftovers.items():
                final.setdefault("segment_group_0", []).append(rec)

        LOGGER.info(f"After singleton clustering, groups: {len(final)}")

    return final

def load_fastas(files):
    fasta = {}
    n_files = 0
    n_seqs = 0
    for f in files:
        n_files += 1
        LOGGER.info("Loading FASTA: %s", f)
        for rec in screed.open(f):
            name = rec.name.split()[0]
            if name in fasta:
                raise ValueError(f"Multiple sequences supplied for {name}")
            fasta[name] = rec.sequence
            n_seqs += 1
    LOGGER.info("Loaded %d sequences from %d FASTA files", n_seqs, n_files)
    return fasta

def write_group_files(group: str, records: List[Dict[str, Any]], taxon: str, outdir: str) -> Tuple[str, str]:
    """Write one metadata JSONL.GZ and one FASTA for a single group."""
    safe_taxon = sanitize_filename(taxon)
    safe_group = sanitize_filename(group)
    prefix = os.path.join(outdir, f"{safe_taxon}-{safe_group}")

    meta_path = prefix + ".metadata.jsonl.gz"
    fasta_path = prefix + ".fa.gz"

    # Write metadata JSONL.GZ
    n_rows = 0
    with gzip.open(meta_path, "wt", encoding="utf-8") as mf:
        for row in records:
            if isinstance(row, dict):
                row = {k: v for k, v in row.items() if k not in ['sequence']}
                mf.write(json.dumps(row, sort_keys=True) + "\n")
                n_rows += 1

    # Write FASTA (one record per accession with available sequence)
    n_seqs = 0
    with gzip.open(fasta_path, "wt", encoding="utf-8") as ff:
        for row in records:
            acc = row.get("accession")
            seq = row.get("sequence")
            if acc and seq:
                ff.write(f">{acc}\n{str(seq).strip()}\n")
                n_seqs += 1

    # LOGGER.info("Group '%s': wrote %d metadata rows, %d sequences", group, n_rows, n_seqs)
    # LOGGER.info("Paths: %s ; %s", meta_path, fasta_path)
    return meta_path, fasta_path

def save_outputs_by_group(groups: Dict[str, List[Dict[str, Any]]], taxon: str, outdir: str) -> Dict[str, Tuple[str, str]]:
    """Write per-group metadata and FASTA files. Returns dict[group] -> (meta_path, fasta_path)."""
    outputs: Dict[str, Tuple[str, str]] = {}
    LOGGER.info("Saving outputs for %d groups", len(groups))
    for grp, recs in groups.items():
        outputs[grp] = write_group_files(grp, recs, taxon, outdir)
    return outputs

def write_manifest(
    outputs: Dict[str, Tuple[str, str]],
    taxon: str,
    outdir: str,
    filename: str = "manifest.csv",
) -> str:
    path = os.path.join(outdir, filename)
    with open(path, "wt", encoding="utf-8") as f:
        f.write("taxon,segment,assembly,metadata\n")
        for seg, files in outputs.items():
            meta, asm = files
            f.write(f"{taxon},{seg},{os.path.basename(asm)},{os.path.basename(meta)}\n")

    LOGGER.info("Wrote manifest: %s", path)
    return path


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    """CLI entry point for merging JSON/CSV metadata into JSONL.GZ."""
    version = "1.0"

    parser = argparse.ArgumentParser(
        description="Merge JSON/CSV by 'accession' into flat JSONL.GZ metadata."
    )
    parser.add_argument("--taxon", default="null", help="Taxon name.")
    parser.add_argument("--assembly", nargs="+", help="Assembly FASTA(s).")
    parser.add_argument("--metadata", nargs="+", help="Metadata file(s).")
    parser.add_argument("--segmented", action="store_true", help="Cluster/merge by segment labels using MinHash/DBSCAN.")
    parser.add_argument("--singletons", action="store_true", help="Keep singletons")
    parser.add_argument("--scaled", type=int, default=10, help="MinHash scaled factor.")
    parser.add_argument("--ksize", type=int, default=16, help="MinHash k-mer size.")
    parser.add_argument("--containment", type=float, default=0.25, help="Minimum fraction of k-mers not contained for segment clustering")
    parser.add_argument("--outdir", default="./", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info("%s v%s", os.path.basename(__file__).replace('.py', ''), version)
    LOGGER.info("Author: Jared Johnson")
    LOGGER.info("Args: taxon=%s, segmented=%s, ksize=%d, scaled=%d, dist=%.3f, outdir=%s",
                args.taxon, args.segmented, args.ksize, args.scaled, args.containment, args.outdir)

    os.makedirs(args.outdir, exist_ok=True)

    # Load all inputs via detector
    metadata: List[Tuple[str, List[Dict[str, Any]]]] = []
    for p in args.metadata:
        try:
            src_path, recs = detect_and_read(p)
            LOGGER.info(f"Read {len(recs)} metadata records from {src_path}")
            metadata.append((src_path, recs))
        except Exception as e:
            LOGGER.exception(f"Error reading {p}: {e}")
            sys.exit(2)

    fastas = load_fastas(args.assembly)

    # Merge and write
    merged = merge_records(metadata, fastas)
    no_na  = drop_na(merged)

    if args.segmented:
        groups = group_segments(no_na, args.ksize, args.scaled, args.containment, args.singletons)
    else:
        groups = {'wg': no_na}
        LOGGER.info("Segmenting disabled: writing a single whole-genome group")

    for k, v in groups.items():
        c = {}
        for r in v:
            s = r.get('segment', 'none')
            if s not in c:
                c[s] = 0
            c[s] += 1

        smry = '; '.join([ f"{k2} ({v2})" for k2, v2 in c.items()])
        LOGGER.info(f"{k}: {smry}")
    outputs = save_outputs_by_group(groups, args.taxon, args.outdir)
    LOGGER.info(f"Completed. Wrote {len(outputs)} group outputs.")

    write_manifest(outputs, args.taxon, args.outdir)


if __name__ == "__main__":
    main()
