#!/usr/bin/env python3
"""
Parse GeneSeekr BLASTN TSV outputs and keep the best BLAST hit for each
separate genomic region in each sample.

Why this script exists:
    The older script kept only the highest percent_match hit across the whole
    genome file. That can drop real additional stx regions if a genome has more
    than one region. This script first groups BLAST hits into genomic regions
    based on overlapping query coordinates, then keeps the best hit within each
    region.

How "best BLAST hit" is chosen:
    1. Highest bit_score
    2. Lowest evalue
    3. Highest percent_match
    4. Longest alignment_length

Input:
    A directory containing files named: *_blastn_geneseekr.tsv

Output:
    parsed_outputs_best_hit_per_region.csv\
        
To run the script, provide the directory path containing the TSV files as an argument:
    python parse_geneseekr_best_hit_per_region.py /path/to/tsv_directory        
        
"""

import csv
import glob
import os
import sys
from typing import Dict, List, Tuple, Any

# Original GeneSeekr BLASTN TSV columns.
ORIGINAL_HEADER = [
    "query_id", "subject_id", "positives", "mismatches", "gaps", "evalue",
    "bit_score", "subject_length", "alignment_length", "query_start", "query_end",
    "subject_start", "subject_end", "percent_match", "query_sequence", "subject_sequence"
]

# Extra columns added by this script.
NEW_HEADER = [
    "sample_name",
    "region_number",
    "region_contig",
    "region_start",
    "region_end",
    "num_hits_in_region",
] + ORIGINAL_HEADER


def safe_float(value: str, default: float = 0.0) -> float:
    """Convert a value to float, returning default if conversion fails."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: str, default: int = 0) -> int:
    """Convert a value to int, accepting strings like '123.0' if needed."""
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default


def read_geneseekr_tsv(tsv_filename: str) -> List[Dict[str, str]]:
    """Read one GeneSeekr TSV file into a list of row dictionaries."""
    rows: List[Dict[str, str]] = []
    with open(tsv_filename, "r", newline="") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            # Skip incomplete/empty rows.
            if not row.get("query_id") or not row.get("subject_id"):
                continue
            rows.append(row)
    return rows


def row_query_interval(row: Dict[str, str]) -> Tuple[str, int, int]:
    """
    Return the query/genome interval for a BLAST hit.

    query_id is the contig name. query_start/query_end are coordinates on that
    contig. We normalize start/end because BLAST alignments can theoretically be
    in either orientation.
    """
    contig = row["query_id"]
    q_start = safe_int(row["query_start"])
    q_end = safe_int(row["query_end"])
    start = min(q_start, q_end)
    end = max(q_start, q_end)
    return contig, start, end


def cluster_hits_by_overlapping_region(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """
    Group BLAST hits into genomic regions based on overlapping query intervals.

    Hits are only grouped together if they are on the same contig and their
    query_start/query_end ranges overlap. This is better than grouping by exact
    coordinates because related BLAST hits for the same biological region can
    have slightly different starts/ends.
    """
    sorted_rows = sorted(rows, key=lambda r: row_query_interval(r))
    regions: List[Dict[str, Any]] = []

    for row in sorted_rows:
        contig, start, end = row_query_interval(row)

        if not regions:
            regions.append({"contig": contig, "start": start, "end": end, "rows": [row]})
            continue

        current = regions[-1]

        # Same contig and overlapping coordinates means same genomic region.
        if contig == current["contig"] and start <= current["end"]:
            current["end"] = max(current["end"], end)
            current["rows"].append(row)
        else:
            regions.append({"contig": contig, "start": start, "end": end, "rows": [row]})

    return regions


def best_blast_hit(rows: List[Dict[str, str]]) -> Dict[str, str]:
    """
    Pick the best BLAST hit from a list of hits from one genomic region.

    BLAST's strongest hit is normally based on bit_score/evalue, not simply
    percent_match. The tie-breakers help keep the choice deterministic.
    """
    return max(
        rows,
        key=lambda r: (
            safe_float(r.get("bit_score")),             # higher is better
            -safe_float(r.get("evalue"), default=1e300), # lower is better
            safe_float(r.get("percent_match")),         # higher is better
            safe_float(r.get("alignment_length")),      # higher is better
        ),
    )


def best_hits_per_region(tsv_filename: str) -> List[Dict[str, Any]]:
    """Return one best hit for each genomic region in one GeneSeekr TSV file."""
    rows = read_geneseekr_tsv(tsv_filename)
    regions = cluster_hits_by_overlapping_region(rows)

    output_rows: List[Dict[str, Any]] = []
    for i, region in enumerate(regions, start=1):
        best = best_blast_hit(region["rows"])
        output_rows.append({
            "region_number": i,
            "region_contig": region["contig"],
            "region_start": region["start"],
            "region_end": region["end"],
            "num_hits_in_region": len(region["rows"]),
            "best_hit": best,
        })

    return output_rows


def sample_name_from_filename(tsv_file: str) -> str:
    """Extract sample name from a filename like SAMPLE_blastn_geneseekr.tsv."""
    filename_base = os.path.basename(tsv_file)
    return filename_base.split("_blastn_geneseekr")[0]


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python parse_geneseekr_best_hit_per_region.py <TSV directory path>")
        sys.exit(1)

    tsv_dir = sys.argv[1]
    tsv_pattern = os.path.join(tsv_dir, "*_blastn_geneseekr.tsv")
    output_csv = os.path.join(tsv_dir, "parsed_outputs_best_hit_per_region.csv")

    tsv_files = sorted(glob.glob(tsv_pattern))
    if not tsv_files:
        print(f"No files found matching: {tsv_pattern}")
        sys.exit(1)

    with open(output_csv, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(NEW_HEADER)

        for tsv_file in tsv_files:
            sample_name = sample_name_from_filename(tsv_file)
            region_results = best_hits_per_region(tsv_file)

            for result in region_results:
                best = result["best_hit"]
                out_row = [
                    sample_name,
                    result["region_number"],
                    result["region_contig"],
                    result["region_start"],
                    result["region_end"],
                    result["num_hits_in_region"],
                ] + [best.get(col, "") for col in ORIGINAL_HEADER]
                writer.writerow(out_row)

    print(f"Done. Wrote: {output_csv}")


if __name__ == "__main__":
    main()
