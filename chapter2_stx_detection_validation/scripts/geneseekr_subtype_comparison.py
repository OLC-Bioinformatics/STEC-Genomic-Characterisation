#!/usr/bin/env python3
"""
Compare parsed GeneSeekr subtype calls with the subtype-level ground truth.

This version:
- uses every row in the parsed GeneSeekr CSV;
- applies no additional percent-identity cutoff;
- requires no separate validation workbook;
- maps GeneSeekr results using the Fasta column already present in the
  shared StxDB subtype-level ground-truth workbook;
- writes True/False as plain text rather than Excel booleans;
- formats the output to match GeneSeekr_Subtype_Comparison_Full.xlsx.

Special non-functional targets:
- 2025-SEQ-1796: Stx2a is interrupted by a transposon.
- SRR18191635: Stx2b contains a frameshift.

A non-functional subtype:
- cannot be a TP;
- is an FP when detected;
- does not become an FN when absent.

Example:
python geneseekr_subtype_comparison_no_cutoff.py \
  --groundtruth StxDB_groundtruth_subtype_level.xlsx \
  --geneseekr parsed_geneseekr_outputs_best_hit_per_region.csv \
  --output GeneSeekr_Subtype_Comparison.xlsx
"""

import argparse
import csv
import re
from artifact_tool import Blob, SpreadsheetFile, Workbook

TOTAL_SUBTYPES = 19

SUBTYPE_ORDER = [
    "Stx1a", "Stx1c", "Stx1d", "Stx1e",
    "Stx2a", "Stx2b", "Stx2c", "Stx2d", "Stx2e",
    "Stx2f", "Stx2g", "Stx2h", "Stx2i", "Stx2j",
    "Stx2k", "Stx2l", "Stx2m", "Stx2n", "Stx2o",
]

NONFUNCTIONAL_SUBTYPES = {
    "2025-SEQ-1796": {"Stx2a"},
    "SRR18191635": {"Stx2b"},
}

SUBTYPE_PATTERN = re.compile(r"Stx[12][a-z]", re.IGNORECASE)


def normalize_subtype(value):
    """Extract and standardize a subtype such as Stx2a."""
    if value is None:
        return None

    match = SUBTYPE_PATTERN.search(str(value))
    if not match:
        return None

    subtype = match.group(0)
    return f"Stx{subtype[3]}{subtype[4].lower()}"


def split_values(value):
    """Split semicolon-, comma-, or newline-separated values."""
    if value is None:
        return []

    return [
        item.strip()
        for item in re.split(r"[;,\n]+", str(value))
        if item.strip()
    ]


def sort_subtypes(values):
    """Sort subtypes according to the recognized 19-subtype universe."""
    return sorted(
        values,
        key=lambda subtype: (
            SUBTYPE_ORDER.index(subtype)
            if subtype in SUBTYPE_ORDER
            else 999,
            subtype,
        ),
    )


def join_subtypes(values):
    return "; ".join(sort_subtypes(values))


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else 0.0


def text_bool(value):
    """Return text so Excel displays True/False rather than icons or checkboxes."""
    return "True" if value else "False"


def load_ground_truth(path):
    """
    Read StxDB_groundtruth_subtype_level.xlsx using the Operon_vs_Subtype sheet with:
    Fastq, Fasta, Expected_Operons, Expected_Subtypes,
    Operon_Count, Unique_Subtype_Count.
    """
    workbook = SpreadsheetFile.import_xlsx(Blob.load(path))
    sheet = workbook.worksheets.get_item("Operon_vs_Subtype")
    rows = sheet.get_range("A1:F10000").values

    if not rows or not rows[0]:
        raise ValueError("The Operon_vs_Subtype sheet is empty.")

    headers = rows[0]
    required = {"Fastq", "Fasta", "Expected_Operons", "Expected_Subtypes"}
    missing = required - set(headers)

    if missing:
        raise ValueError(
            "Operon_vs_Subtype is missing required columns: "
            + ", ".join(sorted(missing))
        )

    records = []
    for row in rows[1:]:
        if not row or row[0] in (None, ""):
            continue

        record = dict(zip(headers, row))
        fastq = str(record["Fastq"]).strip()
        fasta = str(record["Fasta"]).strip()

        expected = {
            subtype
            for item in split_values(record.get("Expected_Subtypes"))
            if (subtype := normalize_subtype(item))
        }

        records.append({
            "Fastq": fastq,
            "Fasta": fasta,
            "Expected_Operons": record.get("Expected_Operons") or "",
            "Expected_Subtypes": expected,
            "Operon_Count": record.get("Operon_Count"),
            "Unique_Subtype_Count": record.get("Unique_Subtype_Count"),
        })

    return records


def load_geneseekr_hits(path):
    """
    Read every retained best-hit-per-region row from the parsed GeneSeekr CSV.

    No percent-identity cutoff is applied here. Percent identity is retained
    only for reporting.
    """
    hits_by_fasta = {}

    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        required = {
            "sample_name",
            "subject_id",
            "percent_match",
            "bit_score",
            "evalue",
            "alignment_length",
        }
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Parsed GeneSeekr CSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            fasta = str(row["sample_name"]).strip()
            subtype = normalize_subtype(row.get("subject_id"))

            if not fasta or subtype is None:
                continue

            try:
                percent_identity = float(row.get("percent_match"))
            except (TypeError, ValueError):
                percent_identity = None

            row["_subtype"] = subtype
            row["_percent_identity"] = percent_identity
            hits_by_fasta.setdefault(fasta, []).append(row)

    return hits_by_fasta


def best_hit_for_subtype(rows):
    """
    Choose the representative GeneSeekr hit for reporting.

    The parsed file already contains the best hit per genomic region. If the
    same subtype occurs in more than one region, retain the strongest regional
    hit using:
      1. highest bit score;
      2. lowest e-value;
      3. highest percent identity;
      4. longest alignment.
    """
    def numeric(value, default=0.0):
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    return sorted(
        rows,
        key=lambda row: (
            -numeric(row.get("bit_score")),
            numeric(row.get("evalue"), float("inf")),
            -numeric(row.get("percent_match")),
            -numeric(row.get("alignment_length")),
        ),
    )[0]


def write_workbook(summary_rows, detail_rows, ground_truth_rows, metric_values, output_path):
    """Create a workbook matching the supplied GeneSeekr comparison format."""
    workbook = Workbook.create()

    summary = workbook.worksheets.add("Summary")
    details = workbook.worksheets.add("Subtype_Details")
    ground_truth = workbook.worksheets.add("Ground_Truth")
    metrics = workbook.worksheets.add("Metrics")

    summary_headers = [
        "Fastq", "Fasta", "Expected_Subtypes", "Detected_Subtypes",
        "TP", "FP", "FN", "TN",
    ]

    detail_headers = [
        "Fastq", "Fasta", "Compared_Subtype", "Expected", "Detected",
        "Best_GeneSeekr_Hit", "Percent_Identity", "Classification",
    ]

    ground_truth_headers = [
        "Fastq", "Fasta", "Expected_Operons", "Expected_Subtypes",
        "Operon_Count", "Unique_Subtype_Count",
    ]

    metric_rows = [
        ["Metric", "Value"],
        ["TP", metric_values["TP"]],
        ["FP", metric_values["FP"]],
        ["FN", metric_values["FN"]],
        ["TN", metric_values["TN"]],
        ["Sensitivity", metric_values["Sensitivity"]],
        ["Specificity", metric_values["Specificity"]],
        ["Precision", metric_values["Precision"]],
        ["Accuracy", metric_values["Accuracy"]],
        ["F1 Score", metric_values["F1 Score"]],
    ]

    summary.get_range_by_indexes(
        0, 0, len(summary_rows) + 1, len(summary_headers)
    ).values = [summary_headers] + summary_rows

    details.get_range_by_indexes(
        0, 0, len(detail_rows) + 1, len(detail_headers)
    ).values = [detail_headers] + detail_rows

    ground_truth.get_range_by_indexes(
        0, 0, len(ground_truth_rows) + 1, len(ground_truth_headers)
    ).values = [ground_truth_headers] + ground_truth_rows

    metrics.get_range("A1:B10").values = metric_rows

    border = {
        "top": {"style": "thin", "color": "#D9D9D9"},
        "bottom": {"style": "thin", "color": "#D9D9D9"},
        "left": {"style": "thin", "color": "#D9D9D9"},
        "right": {"style": "thin", "color": "#D9D9D9"},
    }

    for sheet, rows, columns in [
        (summary, len(summary_rows) + 1, 8),
        (details, len(detail_rows) + 1, 8),
        (ground_truth, len(ground_truth_rows) + 1, 6),
        (metrics, 10, 2),
    ]:
        full_range = sheet.get_range_by_indexes(0, 0, rows, columns)
        full_range.format = {
            "font": {"name": "Calibri", "size": 11},
            "borders": border,
            "wrap_text": True,
            "vertical_alignment": "center",
        }

        header = sheet.get_range_by_indexes(0, 0, 1, columns)
        header.format = {
            "font": {"name": "Calibri", "size": 11, "bold": True},
            "borders": border,
            "wrap_text": True,
            "vertical_alignment": "center",
        }

        sheet.freeze_panes.freeze_rows(1)

    # Match the practical widths of the supplied workbook.
    for column, width in {
        "A": 18, "B": 18, "C": 27, "D": 27,
        "E": 9, "F": 9, "G": 9, "H": 9,
    }.items():
        summary.get_range(f"{column}:{column}").format.column_width = width

    for column, width in {
        "A": 18, "B": 18, "C": 18, "D": 12,
        "E": 12, "F": 27, "G": 16, "H": 15,
    }.items():
        details.get_range(f"{column}:{column}").format.column_width = width

    for column, width in {
        "A": 18, "B": 18, "C": 42,
        "D": 26, "E": 14, "F": 20,
    }.items():
        ground_truth.get_range(f"{column}:{column}").format.column_width = width

    metrics.get_range("A:A").format.column_width = 20
    metrics.get_range("B:B").format.column_width = 16

    # Percent identity is reported as a regular number, as in the template.
    if detail_rows:
        details.get_range(f"G2:G{len(detail_rows) + 1}").format.number_format = "0.00"

    # Performance metrics are percentages.
    metrics.get_range("B6:B10").format.number_format = "0.00%"

    SpreadsheetFile.export_xlsx(workbook).save(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="GeneSeekr subtype comparison with no additional cutoff."
    )
    parser.add_argument(
        "--groundtruth",
        required=True,
        help="StxDB_groundtruth_subtype_level.xlsx containing the Operon_vs_Subtype sheet.",
    )
    parser.add_argument(
        "--geneseekr",
        required=True,
        help="Parsed GeneSeekr best-hit-per-region CSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output comparison workbook.",
    )
    args = parser.parse_args()

    ground_truth = load_ground_truth(args.groundtruth)
    hits_by_fasta = load_geneseekr_hits(args.geneseekr)

    known_fastas = {record["Fasta"] for record in ground_truth}
    unmatched_fastas = sorted(set(hits_by_fasta) - known_fastas)

    if unmatched_fastas:
        preview = ", ".join(unmatched_fastas[:10])
        raise ValueError(
            f"{len(unmatched_fastas)} GeneSeekr sample name(s) were not found "
            f"in the ground truth Fasta column. First unmatched values: {preview}"
        )

    summary_rows = []
    detail_rows = []
    ground_truth_rows = []

    total_tp = total_fp = total_fn = total_tn = 0

    for record in ground_truth:
        fastq = record["Fastq"]
        fasta = record["Fasta"]
        all_curated = set(record["Expected_Subtypes"])
        nonfunctional = NONFUNCTIONAL_SUBTYPES.get(fastq, set())

        functional_expected = all_curated - nonfunctional

        genome_hits = hits_by_fasta.get(fasta, [])
        detected = {row["_subtype"] for row in genome_hits}

        true_positives = detected & functional_expected
        false_positives = detected - functional_expected
        false_negatives = functional_expected - detected

        true_negatives = (
            TOTAL_SUBTYPES
            - len(true_positives)
            - len(false_positives)
            - len(false_negatives)
        )

        if true_negatives < 0:
            raise ValueError(f"Negative TN count calculated for {fastq}.")

        total_tp += len(true_positives)
        total_fp += len(false_positives)
        total_fn += len(false_negatives)
        total_tn += true_negatives

        summary_rows.append([
            fastq,
            fasta,
            join_subtypes(functional_expected),
            join_subtypes(detected),
            len(true_positives),
            len(false_positives),
            len(false_negatives),
            true_negatives,
        ])

        ground_truth_rows.append([
            fastq,
            fasta,
            record["Expected_Operons"],
            join_subtypes(functional_expected),
            record["Operon_Count"],
            len(functional_expected),
        ])

        # TN rows are intentionally excluded, matching the supplied workbook.
        for subtype in sort_subtypes(
            true_positives | false_positives | false_negatives
        ):
            subtype_hits = [
                row for row in genome_hits
                if row["_subtype"] == subtype
            ]

            if subtype_hits:
                best = best_hit_for_subtype(subtype_hits)
                best_hit = best.get("subject_id", "")
                percent_identity = best.get("_percent_identity")
            else:
                best_hit = ""
                percent_identity = ""

            if subtype in true_positives:
                classification = "TP"
            elif subtype in false_positives:
                classification = "FP"
            else:
                classification = "FN"

            detail_rows.append([
                fastq,
                fasta,
                subtype,
                text_bool(subtype in functional_expected),
                text_bool(subtype in detected),
                best_hit,
                percent_identity,
                classification,
            ])

    metric_values = {
        "TP": total_tp,
        "FP": total_fp,
        "FN": total_fn,
        "TN": total_tn,
        "Sensitivity": safe_ratio(total_tp, total_tp + total_fn),
        "Specificity": safe_ratio(total_tn, total_tn + total_fp),
        "Precision": safe_ratio(total_tp, total_tp + total_fp),
        "Accuracy": safe_ratio(
            total_tp + total_tn,
            total_tp + total_fp + total_fn + total_tn,
        ),
        "F1 Score": safe_ratio(
            2 * total_tp,
            2 * total_tp + total_fp + total_fn,
        ),
    }

    write_workbook(
        summary_rows,
        detail_rows,
        ground_truth_rows,
        metric_values,
        args.output,
    )


if __name__ == "__main__":
    main()
