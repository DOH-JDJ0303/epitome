#!/usr/bin/env python3
# epitome_filter.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Tuple, Set, Optional, Any

from epitome_utils import (
    logging_config,
    read_jsonl,
    sanitize_filename,
)

LOGGER = logging_config()

REQUIRED_FIELDS: Tuple[str, ...] = ("taxon", "segment", "variant")
PRESETS: Tuple[str, ...] = ("asm5", "asm10", "asm20")


# ============================================================================
# RECORD FIELD ACCESS
# ============================================================================

def get_field(record: Dict[str, Any], key: str) -> Optional[str]:
    """Extract a field from a record (top-level or nested metadata)."""
    if value := record.get(key):
        return value
    if metadata := record.get("metadata"):
        if isinstance(metadata, dict):
            return metadata.get(key)
    return None


def header_id(record: Dict[str, Any], index: int) -> str:
    """Build a FASTA header from taxon|segment|variant plus a unique index.

    The trailing index keeps IDs unique so each alignment maps back to exactly
    one record; whitespace and '|' in fields are collapsed so the ID stays a
    single token.
    """
    parts = [sanitize_filename(str(get_field(record, k) or "NA")) for k in REQUIRED_FIELDS]
    return "|".join(parts) + f"|{index}"


# ============================================================================
# DATA LOADING
# ============================================================================

def load_records(
    input_paths: List[str],
) -> Tuple[Dict[int, Dict[str, Any]], Set[int]]:
    """
    Load JSONL input files, assigning each record a stable integer index.

    Returns:
        all_records: {index: record_dict} for every record
        with_sequence: set of indices that carry a non-empty sequence
    """
    all_records: Dict[int, Dict[str, Any]] = {}
    with_sequence: Set[int] = set()
    index = 0

    for path in input_paths:
        try:
            count = 0
            for rec in read_jsonl(path):
                all_records[index] = rec
                if rec.get("sequence"):
                    with_sequence.add(index)
                index += 1
                count += 1
            LOGGER.info("Loaded %d records from %s", count, path)
        except Exception as e:
            LOGGER.error("Failed to read %s: %s", path, e)
            sys.exit(2)

    if not all_records:
        raise ValueError("No records found in input - at least one required")

    return all_records, with_sequence


# ============================================================================
# FASTA GENERATION
# ============================================================================

def write_query_fasta(
    all_records: Dict[int, Dict[str, Any]],
    with_sequence: Set[int],
    fasta_path: str,
) -> int:
    """Write a FASTA of all records that carry a sequence. Returns count."""
    written = 0
    with open(fasta_path, "w", encoding="utf-8") as f:
        for index in sorted(with_sequence):
            record = all_records[index]
            f.write(f">{header_id(record, index)}\n{record['sequence']}\n")
            written += 1
    LOGGER.info("Wrote %d query sequences → %s", written, os.path.basename(fasta_path))
    return written


# ============================================================================
# MINIMAP2 ALIGNMENT
# ============================================================================

def run_minimap2(
    query_fasta: str,
    target: str,
    preset: str,
    threads: int,
    paf_path: str,
    extra_args: List[str],
) -> None:
    """Align the query FASTA against the target with minimap2 (PAF output)."""
    exe = shutil.which("minimap2")
    if not exe:
        LOGGER.error("minimap2 not found in PATH (install it, or pass --paf)")
        sys.exit(2)

    # minimap2 argument order is: [options] <target> <query>
    cmd = [exe, "-c", "-x", preset, "-t", str(threads), *extra_args, target, query_fasta]
    LOGGER.info("Running: %s", " ".join(cmd))

    with open(paf_path, "w", encoding="utf-8") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        LOGGER.error("minimap2 exited %d:\n%s", result.returncode, result.stderr)
        sys.exit(2)


# ============================================================================
# PAF PARSING AND COVERAGE
# ============================================================================

def merge_intervals(intervals: List[Tuple[int, int]]) -> int:
    """Total span of merged 0-based half-open [start, end) intervals."""
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    cur_start, cur_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def parse_paf(paf_path: str) -> Dict[int, Dict[str, Dict[str, Any]]]:
    """
    Parse a PAF into per-record, per-target alignment aggregates.

    Returns:
        {index: {target: {"qlen": int, "ivs": [(qs, qe)], "matches": int,
                          "block": int, "n_aln": int}}}
    """
    by_record: Dict[int, Dict[str, Dict[str, Any]]] = {}

    with open(paf_path, encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                LOGGER.warning("PAF line %d: fewer than 12 columns; skipping", line_no)
                continue
            try:
                index = int(fields[0].rsplit("|", 1)[1])
                qlen = int(fields[1])
                qstart, qend = int(fields[2]), int(fields[3])
                target = fields[5]
                matches, block = int(fields[9]), int(fields[10])
            except (IndexError, ValueError):
                LOGGER.warning("PAF line %d: unparseable; skipping", line_no)
                continue

            cell = by_record.setdefault(index, {}).setdefault(
                target, {"qlen": qlen, "ivs": [], "matches": 0, "block": 0, "n_aln": 0}
            )
            cell["ivs"].append((qstart, qend))
            cell["matches"] += matches
            cell["block"] += block
            cell["n_aln"] += 1

    return by_record


# ============================================================================
# FILTERING
# ============================================================================

def evaluate_records(
    by_record: Dict[int, Dict[str, Dict[str, Any]]],
    min_cov: float,
    min_id: float,
) -> Tuple[Set[int], Dict[int, Dict[str, Any]]]:
    """
    Decide which records align to a target above the coverage/identity cutoffs.

    A record is flagged for removal if, for ANY target, the merged query
    coverage >= min_cov AND the identity across that coverage >= min_id.

    Returns:
        to_remove: set of record indices to drop
        stats: {index: {"target", "coverage", "identity", "removed", "n_aln"}}
               describing the deciding (or best-covered) target per record
    """
    to_remove: Set[int] = set()
    stats: Dict[int, Dict[str, Any]] = {}

    for index, targets in by_record.items():
        best: Optional[Tuple[float, float, str]] = None      # highest coverage
        passing: Optional[Tuple[float, float, str]] = None    # passes both cutoffs
        n_aln = 0

        for target, cell in targets.items():
            n_aln += cell["n_aln"]
            qlen = cell["qlen"]
            covered = merge_intervals(cell["ivs"])
            coverage = 100.0 * covered / qlen if qlen else 0.0
            identity = 100.0 * cell["matches"] / cell["block"] if cell["block"] else 0.0

            if coverage >= min_cov and identity >= min_id:
                if passing is None or coverage > passing[0]:
                    passing = (coverage, identity, target)
            if best is None or coverage > best[0] or (
                coverage == best[0] and identity > best[1]
            ):
                best = (coverage, identity, target)

        chosen = passing if passing else best
        stats[index] = {
            "coverage": chosen[0],
            "identity": chosen[1],
            "target": chosen[2],
            "removed": passing is not None,
            "n_aln": n_aln,
        }
        if passing is not None:
            to_remove.add(index)

    return to_remove, stats


# ============================================================================
# OUTPUT GENERATION
# ============================================================================

def write_jsonl(path: str, records: List[Dict[str, Any]]) -> None:
    """Write records to a gzipped JSONL file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def write_summary_tsv(
    path: str,
    all_records: Dict[int, Dict[str, Any]],
    with_sequence: Set[int],
    stats: Dict[int, Dict[str, Any]],
    to_remove: Set[int],
) -> None:
    """Write a per-record decision table as TSV."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("index\ttaxon|segment|variant\tdecision\ttarget\t"
                "coverage_pct\tidentity_pct\tn_aln\n")
        for index in sorted(all_records):
            record = all_records[index]
            name = "|".join(
                sanitize_filename(str(get_field(record, k) or "NA"))
                for k in REQUIRED_FIELDS
            )
            decision = decision_label(index, with_sequence, stats, to_remove)
            st = stats.get(index)
            target = st["target"] if st else "-"
            cov = f"{st['coverage']:.2f}" if st else "-"
            idn = f"{st['identity']:.2f}" if st else "-"
            n_aln = st["n_aln"] if st else 0
            f.write(f"{index}\t{name}\t{decision}\t{target}\t{cov}\t{idn}\t{n_aln}\n")
    LOGGER.info("Wrote per-record summary → %s", os.path.basename(path))


def write_kept_records(
    output_path: str,
    all_records: Dict[int, Dict[str, Any]],
    to_remove: Set[int],
) -> int:
    """Write surviving (non-removed) records to a gzipped JSONL. Returns count."""
    kept = [all_records[i] for i in sorted(all_records) if i not in to_remove]
    write_jsonl(output_path, kept)
    LOGGER.info("Wrote %d kept records → %s", len(kept), os.path.basename(output_path))
    return len(kept)


# ============================================================================
# SUMMARY REPORTING
# ============================================================================

def decision_label(
    index: int,
    with_sequence: Set[int],
    stats: Dict[int, Dict[str, Any]],
    to_remove: Set[int],
) -> str:
    """Human-readable reason a record was kept or removed."""
    if index in to_remove:
        return "removed"
    if index in stats:
        return "kept (below threshold)"
    if index in with_sequence:
        return "kept (no alignment)"
    return "kept (no sequence)"


def log_summary(
    all_records: Dict[int, Dict[str, Any]],
    with_sequence: Set[int],
    stats: Dict[int, Dict[str, Any]],
    to_remove: Set[int],
    min_cov: float,
    min_id: float,
    kept: int,
) -> None:
    """Log the filter summary and the list of removed records."""
    aligned = set(stats)
    no_sequence = set(all_records) - with_sequence
    no_alignment = with_sequence - aligned
    kept_below = aligned - to_remove

    LOGGER.info("=" * 70)
    LOGGER.info("Filter summary (>=%.3g%% coverage & >=%.3g%% identity)", min_cov, min_id)
    LOGGER.info("-" * 70)
    LOGGER.info("  Input records        : %d", len(all_records))
    LOGGER.info("    with sequence      : %d", len(with_sequence))
    LOGGER.info("    without sequence   : %d", len(no_sequence))
    LOGGER.info("  Aligned to a target  : %d", len(aligned))
    LOGGER.info("  Removed              : %d", len(to_remove))
    LOGGER.info("  Kept                 : %d", kept)
    LOGGER.info("    below threshold    : %d", len(kept_below))
    LOGGER.info("    no alignment       : %d", len(no_alignment))
    LOGGER.info("    no sequence        : %d", len(no_sequence))
    LOGGER.info("-" * 70)

    if to_remove:
        LOGGER.info("Removed records:")
        for index in sorted(to_remove):
            record = all_records[index]
            name = "|".join(
                sanitize_filename(str(get_field(record, k) or "NA"))
                for k in REQUIRED_FIELDS
            )
            st = stats[index]
            LOGGER.info(
                "  %-35s → %-14s cov=%.2f%%  id=%.2f%%",
                name[:35], str(st["target"])[:14], st["coverage"], st["identity"],
            )
    else:
        LOGGER.info("Removed records: (none)")
    LOGGER.info("=" * 70)


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Filter JSONL sequence records that align to a reference "
                    "using minimap2 containment (coverage + identity)"
    )
    parser.add_argument("inputs", nargs="+", help="Input JSONL/JSONL.GZ files")
    parser.add_argument("-o", "--output", required=True, help="Output JSONL.GZ of kept records")
    parser.add_argument("--target", help="Reference FASTA or minimap2 .mmi index")
    parser.add_argument("--paf", help="Use an existing PAF instead of running minimap2")
    parser.add_argument("--preset", choices=PRESETS, default="asm5",
                        help="minimap2 divergence preset")
    parser.add_argument("--min-cov", type=float, default=10.0,
                        help="Coverage %% threshold to remove a record")
    parser.add_argument("--min-id", type=float, default=80.0,
                        help="Identity %% threshold to remove a record")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 threads")
    parser.add_argument("--summary", help="Write a per-record decision TSV to this path")
    parser.add_argument("--keep-fasta", help="Also save the generated query FASTA here")
    parser.add_argument("--keep-paf", help="Also save the PAF here")
    parser.add_argument("--mm-extra", default="",
                        help="Extra args passed to minimap2 as one quoted string")
    parser.add_argument("--version", action="version", version="1.0")
    args = parser.parse_args()

    # Log configuration
    LOGGER.info("=" * 70)
    LOGGER.info("epitome_filter v1.0")
    LOGGER.info("=" * 70)
    LOGGER.info("Configuration:")
    LOGGER.info("  Inputs: %d files", len(args.inputs))
    LOGGER.info("  Output: %s", args.output)
    LOGGER.info("  Target: %s", args.target or "(using --paf)")
    LOGGER.info("  Preset: %s", args.preset)
    LOGGER.info("  Min coverage: %.3f", args.min_cov)
    LOGGER.info("  Min identity: %.3f", args.min_id)
    LOGGER.info("  Threads: %d", args.threads)

    # Load data
    LOGGER.info("-" * 70)
    all_records, with_sequence = load_records(args.inputs)
    LOGGER.info("Loaded: %d records (%d with sequence)",
                len(all_records), len(with_sequence))

    workdir = tempfile.mkdtemp(prefix="epitome_filter_")
    try:
        # Obtain PAF: provided, or generated by minimap2
        LOGGER.info("-" * 70)
        if args.paf:
            paf_path = args.paf
            LOGGER.info("Using existing PAF: %s", paf_path)
        else:
            if not args.target:
                parser.error("need --target (FASTA or .mmi) or --paf")
            fasta_path = args.keep_fasta or os.path.join(workdir, "query.fasta")
            n_seq = write_query_fasta(all_records, with_sequence, fasta_path)
            paf_path = args.keep_paf or os.path.join(workdir, "alignments.paf")
            if n_seq:
                run_minimap2(fasta_path, args.target, args.preset, args.threads,
                             paf_path, shlex.split(args.mm_extra))
            else:
                LOGGER.warning("No sequences to align; producing empty PAF")
                open(paf_path, "w").close()

        # Evaluate coverage/identity and decide removals
        LOGGER.info("-" * 70)
        by_record = parse_paf(paf_path)
        to_remove, stats = evaluate_records(by_record, args.min_cov, args.min_id)
        LOGGER.info("Evaluated %d aligned records; %d flagged for removal",
                    len(stats), len(to_remove))

        # Write outputs
        LOGGER.info("-" * 70)
        kept = write_kept_records(args.output, all_records, to_remove)
        if args.summary:
            write_summary_tsv(args.summary, all_records, with_sequence, stats, to_remove)
    finally:
        if not (args.keep_fasta or args.keep_paf):
            shutil.rmtree(workdir, ignore_errors=True)

    # Summary
    log_summary(all_records, with_sequence, stats, to_remove,
                args.min_cov, args.min_id, kept)


if __name__ == "__main__":
    main()