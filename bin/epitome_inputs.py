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
import sourmash as sm
import numpy as np
from sklearn.cluster import DBSCAN

from epitome_utils import sanitize_filename, DistanceCache, logging_config, build_full_minhash_map, _distance, compute_matrix, detect_and_read, normalize_keys

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
    dist: float = 0.25,
    keep_singletons: bool = False,
    unassigned_csv: str = 'unassigned.csv',
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Group sequences into segments based on MinHash containment and segment labels.

    Steps:
      1. Build per-accession MinHash map.
      2. Build initial segment assignments (from 'segment' field) and segment-level hash totals.
      3. Build segment-global MinHashes from those totals.
      4. Compute pairwise segment-global overlaps.
      5. Reassign members to the best-matching segment (if another segment explains them better).
      6. Rebuild segment-global hashes.
      7. Merge similar segments whose globals are highly overlapping.
      8. Rebuild segment-global hashes after merges.
      9. Try to assign previously unassigned accessions based on containment to merged segment globals.
     10. Rebuild segment-global hashes.
     11. Return final {segment: [metadata records]} mapping.

    If unassigned_csv is provided, a CSV of remaining unassigned accessions is written.
    """
    n_records = len(metadata)
    LOGGER.info(
        f"Grouping segments: {n_records} records (ksize={ksize}, scaled={scaled}, "
        f"dist={dist:.3f}, keep_singletons={keep_singletons})"
    )

    if not metadata:
        LOGGER.warning("No metadata records provided to group_segments; returning empty groups.")
        return {}

    # 1) Build MinHash map once
    LOGGER.info(f"Building MinHash map for {n_records} accessions")
    mh_map = build_full_minhash_map(
        {r["accession"]: r["sequence"] for r in metadata},
        ksize=ksize,
        scaled=scaled,
    )

    # 2) First pass: assign to segments (if labeled) and accumulate hashes
    assigned: Dict[str, Dict[str, Any]] = {}
    unassigned: Dict[str, Dict[str, Any]] = {}
    orig_seg_counts: Dict[str, int] = {}

    n_with_seg = 0
    n_without_seg = 0

    for rec in metadata:
        acc = rec["accession"]
        seg = rec.get("segment")

        if seg:
            n_with_seg += 1
            seg = _sanitize_str(seg)

            orig_seg_counts[seg] = orig_seg_counts.get(seg, 0) + 1

            if seg not in assigned:
                assigned[seg] = {"members": [], "hashes": {}}

            assigned[seg]["members"].append(acc)

            # Accumulate hashes for this segment
            _mh = mh_map[acc]
            seg_hashes = assigned[seg]["hashes"]
            for h, cnt in _mh.hashes.items():
                seg_hashes[h] = seg_hashes.get(h, 0) + cnt
        else:
            n_without_seg += 1
            unassigned[acc] = rec

    LOGGER.info(
        f"Initial segment assignment: {n_with_seg} with segment label, {n_without_seg} without"
    )
    LOGGER.info(f"Found {len(assigned)} initial labeled segments")

    def build_segment_globals(assigned_dict: Dict[str, Dict[str, Any]]) -> Dict[str, sm.MinHash]:
        """Helper to build segment-level MinHash objects from assigned members."""
        segment_globals: Dict[str, sm.MinHash] = {}
        for seg, data in assigned_dict.items():
            mh = sm.MinHash(
                n=0,
                scaled=scaled,
                ksize=ksize,
                seed=11,
                track_abundance=True,
            )
            mh.set_abundances(data["hashes"])
            segment_globals[seg] = mh
        return segment_globals

    # 3) Build initial segment-level MinHash objects
    LOGGER.info(f"Building segment-global MinHashes for {len(assigned)} segments")
    segment_globals = build_segment_globals(assigned)

    # 4) Precompute all pairwise overlaps (symmetric)
    segments = list(segment_globals.keys())
    global_overlaps: Dict[str, set] = {k: set() for k in segments}
    LOGGER.info(f"Computing pairwise segment global overlaps (threshold={dist:.3f})")

    for i, k1 in enumerate(segments):
        v1 = segment_globals[k1]
        for k2 in segments[i + 1 :]:  # upper triangle only
            v2 = segment_globals[k2]
            d12 = v1.contained_by_weighted(v2)
            d21 = v2.contained_by_weighted(v1)

            if d12 > dist:
                global_overlaps[k1].add(k2)
            if d21 > dist:
                global_overlaps[k2].add(k1)

    n_with_overlap = sum(1 for s, ovs in global_overlaps.items() if ovs)
    LOGGER.info(
        f"Segment overlap graph: {n_with_overlap} segments with overlaps, "
        f"{len(segments) - n_with_overlap} isolated segments"
    )

    # 5) Reassign members to best-matching segments
    reassignments: List[Tuple[str, str, str]] = []  # (acc, old_seg, new_seg)

    LOGGER.info("Reassigning members to best-matching segments based on containment")

    for seg, overlap in global_overlaps.items():
        if not overlap:
            continue

        candidates = list(overlap) + [seg]
        candidate_mh = {cand: segment_globals[cand] for cand in candidates}

        for acc in assigned[seg]["members"]:
            acc_mh = mh_map[acc]

            max_contain = 0.0
            orig_contain = 0.0
            best_assign = seg

            for cand_seg in candidates:
                d = candidate_mh[cand_seg].contained_by_weighted(acc_mh)
                if cand_seg == seg:
                    orig_contain = d
                if d > max_contain:
                    max_contain = d
                    best_assign = cand_seg

            if best_assign != seg and max_contain > orig_contain:
                reassignments.append((acc, seg, best_assign))

    LOGGER.info(f"Planned {len(reassignments)} member reassignments among segments")

    # Apply reassignments
    for acc, old_seg, new_seg in reassignments:
        LOGGER.info(f"{acc}: {old_seg} -> {new_seg}")
        assigned[old_seg]["members"].remove(acc)
        assigned[new_seg]["members"].append(acc)

    # 6) Rebuild segment-level hash totals after reassignments
    LOGGER.info("Rebuilding segment hash totals after reassignments")
    for seg in assigned:
        seg_hashes: Dict[int, int] = {}
        for acc in assigned[seg]["members"]:
            _mh = mh_map[acc]
            for h, cnt in _mh.hashes.items():
                seg_hashes[h] = seg_hashes.get(h, 0) + cnt
        assigned[seg]["hashes"] = seg_hashes

    segment_globals = build_segment_globals(assigned)

    # 7) Compare globals and merge segments that meet distance threshold
    segments = list(segment_globals.keys())
    merge_map: Dict[str, str] = {}  # source_seg -> target_seg
    LOGGER.info(f"Merging segments based on global containment (threshold={dist:.3f})")

    for i, k1 in enumerate(segments):
        if k1 in merge_map:
            continue  # already marked for merging

        v1 = segment_globals[k1]
        for k2 in segments[i + 1 :]:
            if k2 in merge_map:
                continue

            v2 = segment_globals[k2]
            d12 = v1.contained_by_weighted(v2)
            d21 = v2.contained_by_weighted(v1)

            # If either direction meets threshold, consider merging
            if d12 >= dist or d21 >= dist:
                # Original metadata-based counts
                k1_orig = orig_seg_counts.get(k1, 0)
                k2_orig = orig_seg_counts.get(k2, 0)

                if k1_orig > k2_orig:
                    survivor, merged = k1, k2
                elif k2_orig > k1_orig:
                    survivor, merged = k2, k1
                else:
                    # Tie-break: use current member counts
                    k1_curr = len(assigned[k1]["members"])
                    k2_curr = len(assigned[k2]["members"])
                    if k1_curr >= k2_curr:
                        survivor, merged = k1, k2
                    else:
                        survivor, merged = k2, k1

                merge_map[merged] = survivor
                LOGGER.info(
                    "Merging segment %s -> %s "
                    "(orig_counts: %s=%d, %s=%d; curr_sizes: %s=%d, %s=%d)",
                    merged, survivor,
                    k1, k1_orig, k2, k2_orig,
                    k1, len(assigned[k1]["members"]),
                    k2, len(assigned[k2]["members"]),
                )

                # If k1 got merged into k2, stop comparing k1 to later segments
                if merged == k1:
                    break

    # Apply merges
    LOGGER.info(f"Applying {len(merge_map)} segment merges")
    for source_seg, target_seg in merge_map.items():
        # Follow merge chain to find ultimate target
        final_target = target_seg
        while final_target in merge_map:
            final_target = merge_map[final_target]

        n_source = len(assigned[source_seg]["members"])
        LOGGER.info(
            f"Final merge path: {source_seg} -> {final_target} ({n_source} members moved)"
        )
        assigned[final_target]["members"].extend(assigned[source_seg]["members"])
        del assigned[source_seg]

    # 8) Rebuild segment-level hash totals again after merges
    LOGGER.info("Rebuilding segment hash totals after merges")
    for seg in assigned:
        seg_hashes: Dict[int, int] = {}
        for acc in assigned[seg]["members"]:
            _mh = mh_map[acc]
            for h, cnt in _mh.hashes.items():
                seg_hashes[h] = seg_hashes.get(h, 0) + cnt
        assigned[seg]["hashes"] = seg_hashes

    segment_globals = build_segment_globals(assigned)

    # 9) Assign unassigned sequences that meet distance threshold
    LOGGER.info(
        f"Attempting to assign {len(unassigned)} previously unassigned accessions to segments"
    )
    newly_assigned: List[Tuple[str, str, Dict[str, Any]]] = []

    for acc, rec in unassigned.items():
        acc_mh = mh_map[acc]
        best_seg = None
        max_contain = 0.0

        for seg, seg_mh in segment_globals.items():
            d = seg_mh.contained_by_weighted(acc_mh)
            if d >= dist and d > max_contain:
                max_contain = d
                best_seg = seg

        if best_seg:
            LOGGER.info(
                f"{acc}: unassigned -> {best_seg} (containment={max_contain:.3f})"
            )
            newly_assigned.append((acc, best_seg, rec))

    LOGGER.info(f"Newly assigned {len(newly_assigned)} previously unassigned accessions")

    # Apply new assignments and update unassigned set
    for acc, seg, rec in newly_assigned:
        if seg not in assigned:
            assigned[seg] = {"members": [], "hashes": {}}
        assigned[seg]["members"].append(acc)
        unassigned.pop(acc, None)

    # 10) Rebuild segment-level hash totals one last time
    LOGGER.info("Rebuilding segment hash totals after assigning unassigned sequences")
    for seg in assigned:
        seg_hashes: Dict[int, int] = {}
        for acc in assigned[seg]["members"]:
            _mh = mh_map[acc]
            for h, cnt in _mh.hashes.items():
                seg_hashes[h] = seg_hashes.get(h, 0) + cnt
        assigned[seg]["hashes"] = seg_hashes

    # 11) Optionally log and save remaining unassigned accessions
    n_unassigned_final = len(unassigned)
    if n_unassigned_final:
        LOGGER.info(
            f"After grouping, {n_unassigned_final} accessions remain unassigned to any segment"
        )
        if unassigned_csv:
            try:
                os.makedirs(os.path.dirname(unassigned_csv) or ".", exist_ok=True)
                with open(unassigned_csv, "w", encoding="utf-8") as f:
                    f.write("accession\n")
                    for acc in sorted(unassigned.keys()):
                        f.write(f"{acc}\n")
                LOGGER.info(f"Wrote unassigned accessions to {unassigned_csv}")
            except Exception:
                LOGGER.exception(f"Failed to write unassigned accessions CSV: {unassigned_csv}")
    else:
        LOGGER.info("All accessions assigned to segments")

    # Build final mapping of segment -> list of metadata records
    meta_map = {r["accession"]: r for r in metadata}
    final: Dict[str, List[Dict[str, Any]]] = {}

    for seg, data in assigned.items():
        for acc in data["members"]:
            final.setdefault(seg, []).append(meta_map[acc])

    LOGGER.info(f"Final grouping produced {len(final)} segments")
    for seg, recs in final.items():
        LOGGER.debug(f"Segment {seg}: {len(recs)} members")

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
