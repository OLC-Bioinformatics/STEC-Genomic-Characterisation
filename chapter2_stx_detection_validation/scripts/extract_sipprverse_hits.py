#!/usr/bin/env python3
"""
Extract Sipprverse hits from the original wide Sipprverse CSV output.

This script is meant for the subtype-level validation workflow where:
- the original Sipprverse result file has one row per database target;
- useful hit information is stored in rows ending with "_match_details";
- genome/sample names are the column headers;
- cell values look like: 99.92% (107.28 +/- 24.69)
- sequence rows are ignored completely.

The output contains no sequences and no bracketed numbers. It keeps only:
Genome, Sipprverse_Hit, Subtype, Percent_Identity

Example:
    python extract_sipprverse_hits.py \
        --input sipprverse_output_0.98.csv \
        --output Sipprverse_0.98_extracted_hits_no_sequences.csv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


SUBTYPE_RE = re.compile(r"(Stx\d[a-z])", re.IGNORECASE)
PERCENT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*%")


def normalize_subtype(text: str) -> str:
    """Return subtype with consistent capitalization, e.g. stx2a -> Stx2a."""
    match = SUBTYPE_RE.search(str(text))
    if not match:
        return ""
    subtype = match.group(1).lower()
    return "Stx" + subtype.replace("stx", "")


def extract_percent(cell_value: str):
    """Extract only the percentage before the bracketed numbers."""
    match = PERCENT_RE.search(str(cell_value))
    if not match:
        return None
    return float(match.group(1))


def extract_hits(input_csv: Path) -> pd.DataFrame:
    """Parse Sipprverse CSV and return a long-format extracted hits table."""
    df = pd.read_csv(input_csv, dtype=str)

    if df.empty:
        raise ValueError("Input CSV is empty.")

    sample_col = df.columns[0]
    genome_columns = list(df.columns[1:])

    # Only the *_match_details rows contain the useful identity percentages.
    # The row immediately after each match_details row contains sequence text, which we ignore.
    details = df[df[sample_col].astype(str).str.endswith("_match_details", na=False)].copy()

    records = []

    for _, row in details.iterrows():
        details_name = str(row[sample_col])
        hit_name = details_name.replace("_match_details", "")
        subtype = normalize_subtype(hit_name)

        for genome in genome_columns:
            value = row.get(genome)

            # Empty cells and '-' mean Sipprverse did not report that hit for that genome.
            if pd.isna(value) or str(value).strip() in {"", "-"}:
                continue

            percent_identity = extract_percent(value)

            # If a non-empty cell does not contain a percent, skip it rather than pulling bracketed values.
            if percent_identity is None:
                continue

            records.append(
                {
                    "Genome": genome,
                    "Sipprverse_Hit": hit_name,
                    "Subtype": subtype,
                    "Percent_Identity": percent_identity,
                }
            )

    out = pd.DataFrame.from_records(records)

    if not out.empty:
        out = out.sort_values(["Genome", "Subtype", "Sipprverse_Hit", "Percent_Identity"]).reset_index(drop=True)

    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract Sipprverse hit names, subtypes, and percent identities from original wide CSV output."
    )
    parser.add_argument("--input", "-i", required=True, help="Original Sipprverse CSV file.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file for extracted hits.")
    args = parser.parse_args()

    input_csv = Path(args.input)
    output_csv = Path(args.output)

    if not input_csv.exists():
        print(f"ERROR: input file not found: {input_csv}", file=sys.stderr)
        return 1

    extracted = extract_hits(input_csv)
    extracted.to_csv(output_csv, index=False)

    print(f"Extracted {len(extracted)} Sipprverse hit rows.")
    print(f"Saved: {output_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
