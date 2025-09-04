#!/usr/bin/env python3

# epitome_ncbi.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import logging
import sys
from pathlib import Path
import screed
import re
import os
from typing import Optional, Iterable, Mapping, Sequence

from epitome_utils import load_json_or_jsonl, logging_config

# -------------------------------
#  GLOBAL CONFIG
# -------------------------------

LOGGER = logging_config()

# -------------------------------
#  FUNCTIONS
# -------------------------------

def load_fasta_map(path: str) -> dict[str, str]:
    """Load FASTA(.gz ok via screed) and map both full and base accessions to sequences."""
    LOGGER.info(f"Loading FASTA sequences: %s", path)
    seq_by_acc: dict[str, str] = {}
    count = 0
    for rec in screed.open(path):
        header = rec.name
        acc_full = header.split()[0]
        seq = rec.sequence
        if acc_full not in seq_by_acc:
            seq_by_acc[acc_full] = seq
        acc_base = acc_full.split(".", 1)[0]
        if acc_base not in seq_by_acc:
            seq_by_acc[acc_base] = seq
        count += 1
        if count % 10000 == 0:
            LOGGER.debug("Processed %d FASTA records...", count)
    LOGGER.info("Loaded %d sequences into FASTA map (keys include full and base accessions)", len(seq_by_acc))
    return seq_by_acc


SUBTYPE_TARGETS = {
    "subtype", "serotype", "genotype", "serogroup", "variety",
    "biotype", "lineage", "clade", "subclade", "pathotype", "group", "subgroup"
}

def _format_segment_name(seg: Optional[str]) -> Optional[str]:
    """Normalize a raw segment name string."""
    if seg is None:
        return seg
    seg = seg.replace(r'\;', ';').strip()
    seg = re.sub(r'\s+', ' ', seg)
    seg = seg.lower()
    return seg


def parse_edirect_json_docsum(edirect_json_obj: Iterable[Mapping]) -> dict[str, dict[str, str]]:
    """Extract subtype-like attributes from EDirect esummary/efetch docsum JSON."""
    LOGGER.info("Parsing EDirect JSON docsum for subtype-like attributes")
    data: dict[str, dict[str, str]] = {}
    rows = 0
    for row in edirect_json_obj:
        rows += 1
        result = row.get("result", row)
        if not isinstance(result, dict):
            LOGGER.debug("Skipping non-dict result at row %d", rows)
            continue
        uids = result.get("uids")
        records = [result.get(u, {}) for u in uids] if isinstance(uids, list) else [result]
        for value in records:
            if not isinstance(value, dict):
                continue
            accession = value.get("accessionversion") or value.get("AccessionVersion")
            if not accession:
                continue
            subtype_keys = value.get("SubType") or value.get("subtype") or []
            subtype_vals = value.get("SubName") or value.get("subname") or []
            if isinstance(subtype_keys, str):
                subtype_keys = subtype_keys.split("|")
            if isinstance(subtype_vals, str):
                subtype_vals = subtype_vals.split("|")
            m: dict[str, str] = {}
            for i, k in enumerate(subtype_keys):
                k_norm = (k or "").strip().lower()
                if k_norm in SUBTYPE_TARGETS and i < len(subtype_vals):
                    v = str(subtype_vals[i]).strip()
                    if v:
                        m[k_norm] = v
            if m:
                data[accession] = m
    LOGGER.info("Extracted subtype data for %d accessions (from %d rows)", len(data), rows)
    return data


def extract_taxids(ds_taxa_json: list[dict]) -> dict[str, str]:
    """Extract a mapping from taxid to species name from NCBI datasets taxonomy JSON."""
    LOGGER.info("Extracting taxid -> species mapping")
    out: dict[str, str] = {}

    if not ds_taxa_json or not isinstance(ds_taxa_json[0], dict):
        LOGGER.warning("Unexpected taxonomy JSON structure; returning empty map")
        return out

    reports = ds_taxa_json[0].get("reports", [])
    for idx, row in enumerate(reports, 1):
        taxonomy = row.get("taxonomy", {})
        name = taxonomy.get("current_scientific_name", {}).get("name")
        primary_tid = taxonomy.get("tax_id")
        secondary_tids = taxonomy.get("secondary_tax_ids", [])

        if not name or not primary_tid:
            LOGGER.debug("No species info in report #%d", idx)
            continue

        tids = [primary_tid] + (secondary_tids if isinstance(secondary_tids, list) else [])
        for taxid in tids:
            if str(taxid) in out:
                LOGGER.debug("Overwriting existing mapping for taxid %s", taxid)
            out[str(taxid)] = name

    LOGGER.info("Extracted %d taxids", len(out))
    return out



def lineage_to_species_and_taxid(
    lineage_list: Optional[list],
    taxid_name_map: dict[object, str]
) -> tuple[Optional[str], Optional[object]]:
    """Map a lineage list to (species_name, taxid) using a taxid->name map."""
    if not isinstance(lineage_list, list):
        return None, None
    for item in lineage_list:
        taxid = item.get("taxId")
        if taxid is None:
            continue
        taxid = str(taxid)
        if taxid in taxid_name_map:
            return taxid_name_map[taxid], taxid
    return None, None


def write_segment_files(
    taxon: str,
    seg: str,
    seg_rows: Sequence[dict],
    fasta_map: dict[str, str],
    out_dir: Path
) -> list[str]:
    """Write per-segment JSON and FASTA files."""
    LOGGER.info("Writing files for segment '%s' (n=%d)", seg, len(seg_rows))
    file_prefix = f"{taxon}-{seg}"
    json_path = out_dir / f"{file_prefix}.json"
    fasta_path = out_dir / f"{file_prefix}.fa.gz"

    # JSON
    import json
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(seg_rows, f, ensure_ascii=False, indent=2)
    LOGGER.debug("Wrote JSON: %s", json_path)

    # FASTA
    import gzip
    n_fasta = 0
    seen: set[str] = set()
    with gzip.open(f"{fasta_path}", "wt", encoding="utf-8") as f:
        for r in seg_rows:
            acc = r.get("accession")
            if not acc:
                LOGGER.debug("Skipping row without accession in segment %s", seg)
                continue
            seq = fasta_map.get(acc) or fasta_map.get(acc.split(".", 1)[0])
            if not seq:
                LOGGER.warning("No sequence found for accession %s in FASTA map", acc)
                continue
            if acc in seen:
                LOGGER.debug("Duplicate accession in segment %s: %s (skipping)", seg, acc)
                continue
            seen.add(acc)
            f.write(f">{acc}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            n_fasta += 1
    LOGGER.info("Wrote FASTA: %s (records=%d)", fasta_path, n_fasta)

    return [taxon, seg, file_prefix]


def split_by_segment(
    taxon: str,
    rows: list[dict],
    fasta_map: dict[str, str],
    out_dir: Path
) -> None:
    """Split rows by segment, write per-segment JSON/FASTA, and a manifest CSV."""
    LOGGER.info("Splitting %d sequences into segments", len(rows))
    by_seg: dict[str, list[dict]] = {}
    for r in rows:
        seg = r.get("segment") or "wg"
        by_seg.setdefault(seg, []).append(r)

    man: list[list[str]] = [["taxon", "segment", "prefix"]]
    for seg, seg_rows in by_seg.items():
        man_row = write_segment_files(taxon, seg, seg_rows, fasta_map, out_dir)
        man.append(man_row)

    man_path = out_dir / f"{taxon}.manifest.csv"
    with man_path.open("w", newline="", encoding="utf-8") as f:
        out = '\n'.join([','.join(row) for row in man])
        f.write(out)

    LOGGER.info("Manifest written: %s (segments=%d)", man_path, len(by_seg))


# -------------------------------
#  MAIN
# -------------------------------

def main() -> None:
    version = "2.0"

    """CLI entry point: merge NCBI datasets & EDirect docsum; write per-segment outputs."""
    parser = argparse.ArgumentParser(
        description="Merge NCBI datasets & EDirect subtype data, then write segment-specific outputs."
    )
    parser.add_argument("--taxon", required=True, help="Taxon name.")
    parser.add_argument("--datasets_genome_fasta", required=True,
                        help="Path to genomic.fna from 'datasets download virus genome taxon'.")
    parser.add_argument("--datasets_genome_json", required=True,
                        help="Path to 'data_report.jsonl' from 'datasets download virus genome taxon'.")
    parser.add_argument("--datasets_taxonomy_json", required=True,
                        help="Path to JSON from 'datasets summary taxonomy taxon'.")
    parser.add_argument("--edirect_json", required=True,
                        help="Path to JSON or JSONL from EDirect esummary/efetch docsum.")
    parser.add_argument("--segmented", action="store_true",
                        help="Set if the virus is segmented (segment names used when present).")
    parser.add_argument("--segment_synonyms", default=None,
                        help="Optional semicolon list of segment synonyms 's|small;m|medium|middle;l|large'.")
    parser.add_argument("--out_dir", default=".", help="Output directory.")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    args = parser.parse_args()

    LOGGER.info(f"{os.path.basename(__file__).replace('.py', '')} v{version}")
    LOGGER.info(f"Author: Jared Johnson")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.debug("Output directory ensured: %s", out_dir)

    try:
        fasta_map = load_fasta_map(args.datasets_genome_fasta)
    except Exception as e:
        LOGGER.exception("Failed to load FASTA: %s", e)
        sys.exit(2)

    try:
        ds_genome_list = load_json_or_jsonl(args.datasets_genome_json)
        LOGGER.info("Loaded genome report (%s)", args.datasets_genome_json)
    except Exception as e:
        LOGGER.exception("Failed to load genome JSON/JSONL: %s", e)
        sys.exit(2)

    try:
        ds_taxa_json = load_json_or_jsonl(args.datasets_taxonomy_json)
        LOGGER.info("Loaded taxonomy summary (%s)", args.datasets_taxonomy_json)
    except Exception as e:
        LOGGER.exception("Failed to load taxonomy JSON/JSONL: %s", e)
        sys.exit(2)

    try:
        edirect_json_obj = load_json_or_jsonl(args.edirect_json)
        LOGGER.info("Loaded EDirect docsum (%s)", args.edirect_json)
    except Exception as e:
        LOGGER.exception("Failed to load EDirect JSON/JSONL: %s", e)
        sys.exit(2)

    taxid_name_map = extract_taxids(ds_taxa_json)
    LOGGER.debug("taxid->species entries: %d", len(taxid_name_map))

    subtypes_map = parse_edirect_json_docsum(edirect_json_obj)
    LOGGER.debug("subtype entries: %d", len(subtypes_map))

    seg_syn_map: dict[str, str] = {}
    if args.segment_synonyms:
        LOGGER.info("Processing segment synonyms: %s", args.segment_synonyms)
        groups = re.split(r'(?<!\\);', args.segment_synonyms)
        for grp in groups:
            token_list = [_format_segment_name(t) for t in grp.split("|") if t.strip()]
            token_set = set(token_list)
            canonical = token_list[0]
            for t in token_set:
                seg_syn_map[t] = canonical
        LOGGER.debug("Segment synonym map size: %d; sample: %s",
                    len(seg_syn_map),
                    dict(list(seg_syn_map.items())[:10]))

    rows: list[dict] = []
    missing_accession = 0
    missing_seg_syns = set()
    for i, rec in enumerate(ds_genome_list, 1):
        accession = rec.get("accession")
        if not accession:
            missing_accession += 1
            if missing_accession <= 5:
                LOGGER.warning(f"Record #{i} missing 'accession'; skipping")
            continue

        segment = rec.get("segment")
        if args.segmented:
            if not segment:
                LOGGER.warning(f"Record #{i} missing 'segment'; skipping")
                continue
            segment = _format_segment_name(segment)
            if segment not in seg_syn_map:
                missing_seg_syns.add(segment)
            segment = seg_syn_map.get(segment)
        else:
            segment = "wg"

        virus_lineage = (rec.get("virus") or {}).get("lineage") or []
        species_name, species_taxid = lineage_to_species_and_taxid(virus_lineage, taxid_name_map)

        host_lineage = rec.get("host", {}).get("lineage") or []
        host = host_lineage[-1].get("name") if host_lineage else None
        collection_date = rec.get("isolate", {}).get("collectionDate")
        geographic_region = rec.get("location", {}).get("geographicRegion")

        base = {
            "accession": accession,
            "taxon": args.taxon,
            "segment": segment,
            "species": species_name,
            "taxId": species_taxid,
            "host": host,
            "collectionDate": collection_date,
            "geographicRegion": geographic_region,
        }

        sub = (subtypes_map.get(accession) or {})
        rows.append({**base, **sub})

        if i % 1000 == 0:
            LOGGER.debug("Processed %d genome records...", i)
    
    if missing_seg_syns:
        LOGGER.warning(f"The following segment names do not match provided synonyms: {missing_seg_syns}")

    if missing_accession:
        LOGGER.warning("Skipped %d records without accession", missing_accession)

    unique_segments = sorted({r["segment"] for r in rows if "segment" in r and r["segment"]})
    LOGGER.info("Observed %d unique segment names: %s", len(unique_segments), ", ".join(unique_segments))

    split_by_segment(args.taxon, rows, fasta_map, out_dir)

    LOGGER.info("Processing complete for taxon '%s' (rows=%d, outdir=%s)", args.taxon, len(rows), out_dir)


if __name__ == "__main__":
    main()
