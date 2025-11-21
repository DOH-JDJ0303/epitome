#!/usr/bin/env python3
import csv
import json
import logging
import sys
import argparse
import re
from typing import List, Dict, Any, Mapping, Set
import screed
import gzip

# ----------------------------
# Minimal logging
# ----------------------------
def setup_logging(level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("meta_clean")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    if not logger.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logger.addHandler(h)
    return logger

# ----------------------------
# Regex
# ----------------------------
RE_EPI_ISL = re.compile(r"\bEPI_ISL_\d+\b", re.IGNORECASE)

# ----------------------------
# IO helpers
# ----------------------------
def read_delim(path: str) -> List[Dict[str, Any]]:
    """Load a CSV/TSV file into a list of dicts; trims strings and converts empty to None.
       Auto-detects delimiter as ',' or '\t'.
    """
    recs: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        first_line = f.readline()
        delimiter = "\t" if first_line.count("\t") > first_line.count(",") else ","

    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header row")
        for row in reader:
            clean: Dict[str, Any] = {}
            for k, v in row.items():
                k = k.strip() if isinstance(k, str) else k
                if isinstance(v, str):
                    v = v.strip()
                    v = v if v != "" else None
                clean[k] = v
            recs.append(clean)
    return recs

def read_metadata(paths: List[str], logger: logging.Logger) -> List[Dict[str, Any]]:
    meta: List[Dict[str, Any]] = []
    for p in paths:
        chunk = read_delim(p)
        logger.info(f"Loaded {len(chunk)} rows from {p}")
        meta.extend(chunk)
    logger.info(f"Total metadata rows: {len(meta)}")
    return meta

def read_fastas(paths: List[str], logger: logging.Logger, extract_epi_isl: bool = False) -> Dict[str, str]:
    """Return dict of header->sequence. If extract_epi_isl, use first EPI_ISL_#### found in header as key."""
    fa: Dict[str, str] = {}
    total = normalized = no_match = collisions = 0
    for p in paths:
        added = 0
        with screed.open(p) as fh:
            for r in fh:
                total += 1
                hdr = r.name
                key = hdr
                if extract_epi_isl:
                    m = RE_EPI_ISL.search(hdr)
                    if m:
                        key = m.group(0)
                        normalized += 1
                    else:
                        no_match += 1  # kept as original header
                if key in fa and fa[key] != r.sequence:
                    collisions += 1  # later file overwrites earlier entry
                fa[key] = r.sequence
                added += 1
        logger.info(f"Loaded {added} sequences from {p}")
    if extract_epi_isl:
        logger.info(f"FASTA header normalization: total={total}, normalized={normalized}, no_epi_isl_found={no_match}, key_collisions={collisions}")
    logger.info(f"Total sequences indexed: {len(fa)}")
    return fa

# ----------------------------
# Domain transforms
# ----------------------------
def parse_flu(data: List[Dict[str, Any]], logger: logging.Logger) -> List[Dict[str, Any]]:
    """Flatten per-segment columns -> rows with accession and segment number (1-8)."""
    new_data: List[Dict[str, Any]] = []
    seg_map = {"PB2":"1","PB1":"2","PA":"3","HA":"4","NP":"5","NA":"6","MP":"7","NS":"8","HE":"4","P3":"3"}
    if not data:
        return new_data
    segment_cols = [k for k in data[0].keys() if isinstance(k, str) and k.endswith(" Segment_Id")]
    logger.info(f"Detected {len(segment_cols)} influenza segment columns: {segment_cols}")
    for rec in data:
        for col in segment_cols:
            gene_name = col.replace(" Segment_Id", "")
            seg_num = seg_map.get(gene_name)
            col_v = rec.get(col)
            if not col_v or not seg_num or not str(col_v).startswith('EPI'):
                continue
            acc = str(col_v).split("|")[0].strip()  # Accession is before first pipe
            rec_copy = dict(rec)  # avoid mutating original across multiple segments
            rec_copy["Accession ID"] = acc
            rec_copy["segment"] = seg_num
            new_data.append(rec_copy)
    logger.info(f"Influenza rows expanded to {len(new_data)} records")
    return new_data

def normalize_keys(d: Mapping[str, Any]) -> Dict[str, Any]:
    """Lowercase snake_case keys."""
    def to_snake(s: str) -> str:
        s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
        s = s.replace(" ", "_")
        return s.lower()
    return {to_snake(str(k)): v for k, v in d.items()}

# ----------------------------
# Main
# ----------------------------
def main():
    version = "2.6"

    parser = argparse.ArgumentParser(description="Filter metadata/FASTA to matched records and emit CSV+FASTA (strict accession==header)")
    parser.add_argument("--metadata", nargs="+", required=True, help="Path(s) to metadata CSV/TSV file(s)")
    parser.add_argument("--fasta", nargs="+", required=True, help="Path(s) to multi-FASTA file(s)")
    parser.add_argument("--out-csv", default='metadata.csv.gz', help="Output CSV with matched metadata")
    parser.add_argument("--out-fasta", default='sequences.fa.gz', help="Output FASTA with matched sequences")
    parser.add_argument("--species", help="Species name to add to each matched metadata record")
    parser.add_argument("--flu", action="store_true", help="Metadata is influenza (use segment mapping)")
    # accept both spellings
    parser.add_argument("--extract-epi-isl", dest="extract_epi_isl", action="store_true",
                        help="Extract EPI_ISL_#### from FASTA headers and use it as the header key for matching/output.")
    parser.add_argument("--extract_epi_isl", dest="extract_epi_isl", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    logger = setup_logging(args.log_level)
    logger.info(f"meta_clean v{version}")
    logger.info(f"Influenza mode: {args.flu}")
    if args.species:
        logger.info(f"Species tag: {args.species}")
    if args.extract_epi_isl:
        logger.info("FASTA headers will be normalized to EPI_ISL_####")
        if args.flu:
            logger.warning("Using --extract-epi-isl with --flu may prevent matches (flu metadata uses EPI##### segment accessions).")

    meta = read_metadata(args.metadata, logger)
    fa = read_fastas(args.fasta, logger, extract_epi_isl=args.extract_epi_isl)

    if not meta:
        logger.warning("No metadata rows provided.")
        return
    if not fa:
        logger.warning("No FASTA sequences provided.")
        return

    if args.flu:
        meta = parse_flu(meta, logger)
    else:
        sample_row = meta[0] if meta else {}
        if any(isinstance(k, str) and k.endswith(" Segment_Id") for k in sample_row.keys()):
            logger.warning("Metadata includes '* Segment_Id' columns; consider '--flu' mode.")

    matched_rows: List[Dict[str, Any]] = []
    matched_headers: Set[str] = set()
    missing_accession_idx: List[int] = []
    missing_in_fasta_vals: List[str] = []

    for idx, rec in enumerate(meta):
        r = normalize_keys(rec)

        # accession_id -> accession (normalized)
        acc = r.get("accession") or r.get("accession_id")
        if not acc:
            missing_accession_idx.append(idx)
            continue
        acc = str(acc).strip()
        if "accession" not in r:
            r["accession"] = acc
        if "location" in r and not r.get("geographic_region"):
            r["geographic_region"] = r["location"]

        # STRICT match: accession must exactly equal a FASTA key (possibly normalized EPI_ISL)
        if acc in fa:
            if args.species:
                r["species"] = args.species.strip()
            cleaned: Dict[str, Any] = {k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v)
                                       for k, v in r.items()}
            matched_rows.append(cleaned)
            matched_headers.add(acc)
        else:
            missing_in_fasta_vals.append(acc)

    if missing_in_fasta_vals:
        logger.warning(f"Preview of unmatched metadata accessions: {missing_in_fasta_vals[:5]}")

    # FASTA headers with no matching metadata
    unmatched_fasta_headers = sorted([h for h in fa.keys() if h not in matched_headers])
    if unmatched_fasta_headers:
        logger.warning(f"Preview of unmatched FASTA headers: {unmatched_fasta_headers[:5]}")

    # Summaries
    logger.info(f"Matched metadata rows: {len(matched_rows)}")
    logger.info(f"Matched FASTA headers: {len(matched_headers)}")
    logger.info(f"Metadata rows missing accession: {len(missing_accession_idx)}")
    logger.info(f"Metadata accessions not in FASTA: {len(missing_in_fasta_vals)}")
    logger.info(f"FASTA headers without metadata: {len(unmatched_fasta_headers)}")

    # ----------------------------
    # Emit CSV (normalized keys)
    # ----------------------------
    if matched_rows:
        # Union of keys; accession, segment, species first
        all_keys: Set[str] = set().union(*[set(r.keys()) for r in matched_rows])
        ordered = []
        for k in ("accession", "segment", "species"):
            if k in all_keys:
                ordered.append(k)
                all_keys.remove(k)
        ordered.extend(sorted(all_keys))
        with gzip.open(args.out_csv, "wt", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=ordered)
            w.writeheader()
            for r in matched_rows:
                w.writerow({k: r.get(k) for k in ordered})
        logger.info(f"Wrote matched metadata CSV: {args.out_csv} ({len(matched_rows)} rows)")
    else:
        logger.warning("No matched metadata rows; CSV not written.")

    # ----------------------------
    # Emit FASTA (only exact-header matches; preserve original read order via insertion order)
    # ----------------------------
    if matched_headers:
        with gzip.open(args.out_fasta, "wt", encoding="utf-8") as f:
            for hdr, seq in fa.items():
                if hdr in matched_headers:
                    f.write(f">{hdr}\n{seq}\n")
        logger.info(f"Wrote matched FASTA: {args.out_fasta} ({len(matched_headers)} sequences)")
    else:
        logger.warning("No matched FASTA headers; FASTA not written.")

if __name__ == '__main__':
    main()
