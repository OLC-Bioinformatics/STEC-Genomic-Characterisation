#!/usr/bin/env python3
"""
NCBI StxTyper subtype-level comparison using COMPLETE operons only.

Inputs
------
1. StxTyper combined-results TSV.
2. StxDB_groundtruth_subtype_level.xlsx, sheet: Operon_vs_Subtype.

Rules
-----
- Only rows whose operon column is exactly COMPLETE are used.
- PARTIAL, PARTIAL_CONTIG_END, COMPLETE_NOVEL, and all other statuses are ignored.
- COMPLETE calls are converted to recognized subtypes and duplicate calls of
  the same subtype within a genome are collapsed.
- If an expected functional subtype has no matching COMPLETE StxTyper call,
  it is counted as an FN.
- A COMPLETE subtype call absent from the functional ground truth is an FP.
- TN = 19 - TP - FP - FN per genome.
- 2025-SEQ-1796: Stx2a is interrupted and excluded from expected positives.
- SRR18191635: Stx2b has a frameshift and is excluded from expected positives.

Example
-------
python stxtyper_subtype_comparison.py \
  --groundtruth StxDB_groundtruth_subtype_level.xlsx \
  --stxtyper stxtyper_combined_results.tsv \
  --output StxTyper_Subtype_Comparison.xlsx
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

SUBTYPE_PATTERN = re.compile(r"^stx([12])([a-z])$", re.IGNORECASE)


def normalize_subtype(value):
    text = str(value or "").strip()
    match = SUBTYPE_PATTERN.fullmatch(text)
    if not match:
        return None
    return f"Stx{match.group(1)}{match.group(2).lower()}"


def split_values(value):
    return [
        item.strip()
        for item in re.split(r"[;,\n]+", str(value or ""))
        if item.strip()
    ]


def sort_subtypes(values):
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


def text_bool(value):
    return "True" if value else "False"


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_ground_truth(path):
    workbook = SpreadsheetFile.import_xlsx(Blob.load(path))
    sheet = workbook.worksheets.get_item("Operon_vs_Subtype")
    rows = sheet.get_range("A1:F10000").values

    if not rows:
        raise ValueError("The Operon_vs_Subtype sheet is empty.")

    headers = rows[0]
    required = {
        "Fasta", "Fastq", "Expected_Operons", "Expected_Subtypes",
        "Operon_Count", "Unique_Subtype_Count",
    }
    missing = required - set(headers)
    if missing:
        raise ValueError(
            "Ground truth is missing required columns: "
            + ", ".join(sorted(missing))
        )

    records = []
    for row in rows[1:]:
        if not row or row[0] in (None, ""):
            continue

        record = dict(zip(headers, row))
        expected = {
            subtype
            for item in split_values(record.get("Expected_Subtypes"))
            if (subtype := normalize_subtype(str(item).lower()))
        }

        records.append({
            "Fasta": str(record["Fasta"]).strip(),
            "Fastq": str(record["Fastq"]).strip(),
            "Expected_Operons": record.get("Expected_Operons") or "",
            "All_Expected_Subtypes": expected,
            "Operon_Count": record.get("Operon_Count"),
            "Unique_Subtype_Count": record.get("Unique_Subtype_Count"),
        })

    return records


def load_complete_stxtyper_calls(path):
    """
    Return only StxTyper rows whose operon status is exactly COMPLETE.
    """
    calls_by_genome = {}

    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "Genome", "stx_type", "operon", "identity", "target_contig",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "StxTyper TSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            if str(row.get("operon") or "").strip().upper() != "COMPLETE":
                continue

            genome = str(row.get("Genome") or "").strip()
            subtype = normalize_subtype(row.get("stx_type"))

            if not genome:
                continue
            if subtype is None:
                raise ValueError(
                    f"COMPLETE StxTyper row for {genome} does not contain "
                    f"a resolved subtype: {row.get('stx_type')!r}"
                )

            row["_subtype"] = subtype
            row["_identity"] = safe_float(row.get("identity"))
            calls_by_genome.setdefault(genome, []).append(row)

    return calls_by_genome


def best_call(rows):
    def identity_value(row):
        value = row.get("_identity")
        return value if value is not None else -1.0

    return sorted(
        rows,
        key=lambda row: (
            -identity_value(row),
            str(row.get("target_contig") or ""),
        ),
    )[0]


def build_workbook(ground_truth, calls_by_genome, output_path):
    known_fastas = {record["Fasta"] for record in ground_truth}
    unmatched = sorted(set(calls_by_genome) - known_fastas)

    if unmatched:
        preview = ", ".join(unmatched[:10])
        raise ValueError(
            f"{len(unmatched)} StxTyper Genome value(s) were not found in "
            f"the ground-truth Fasta column. First unmatched: {preview}"
        )

    summary_rows = []
    detail_rows = []
    ground_rows = []

    for record in ground_truth:
        fasta = record["Fasta"]
        fastq = record["Fastq"]

        all_expected = set(record["All_Expected_Subtypes"])
        nonfunctional = NONFUNCTIONAL_SUBTYPES.get(fastq, set())
        functional_expected = all_expected - nonfunctional

        complete_calls = calls_by_genome.get(fasta, [])
        detected = {row["_subtype"] for row in complete_calls}

        tp = detected & functional_expected
        fp = detected - functional_expected
        fn = functional_expected - detected
        tn = TOTAL_SUBTYPES - len(tp) - len(fp) - len(fn)

        if tn < 0:
            raise ValueError(f"Negative TN calculated for {fastq}.")

        summary_rows.append([
            fastq,
            fasta,
            join_subtypes(functional_expected),
            join_subtypes(detected),
            len(tp),
            len(fp),
            len(fn),
            tn,
        ])

        ground_rows.append([
            fasta,
            fastq,
            record["Expected_Operons"],
            join_subtypes(all_expected),
            join_subtypes(functional_expected),
            join_subtypes(nonfunctional),
            record["Operon_Count"],
            len(functional_expected),
        ])

        # Include TP, FP, and FN rows only. TN rows are omitted.
        for subtype in sort_subtypes(tp | fp | fn):
            matches = [
                row for row in complete_calls
                if row["_subtype"] == subtype
            ]

            if matches:
                retained = best_call(matches)
                stxtyper_call = retained.get("stx_type", "")
                identity = retained.get("_identity")
                operon_status = retained.get("operon", "")
                target_contig = retained.get("target_contig", "")
            else:
                stxtyper_call = ""
                identity = ""
                operon_status = ""
                target_contig = ""

            classification = (
                "TP" if subtype in tp
                else "FP" if subtype in fp
                else "FN"
            )

            detail_rows.append([
                fastq,
                fasta,
                subtype,
                text_bool(subtype in functional_expected),
                text_bool(subtype in detected),
                text_bool(subtype in nonfunctional),
                stxtyper_call,
                identity,
                operon_status,
                target_contig,
                classification,
            ])

    workbook = Workbook.create()
    summary = workbook.worksheets.add("Summary")
    details = workbook.worksheets.add("Subtype_Details")
    ground = workbook.worksheets.add("Ground_Truth")
    metrics = workbook.worksheets.add("Metrics")

    summary_headers = [
        "Fastq", "Fasta", "Expected_Subtypes", "Detected_Subtypes",
        "TP", "FP", "FN", "TN",
    ]
    detail_headers = [
        "Fastq", "Fasta", "Compared_Subtype", "Expected", "Detected",
        "Nonfunctional_Special_Case", "StxTyper_Call", "Identity",
        "Operon_Status", "Target_Contig", "Classification",
    ]
    ground_headers = [
        "Fasta", "Fastq", "Expected_Operons", "All_Curated_Subtypes",
        "Expected_Functional_Subtypes", "Nonfunctional_Subtypes",
        "Operon_Count", "Functional_Subtype_Count",
    ]

    summary.get_range_by_indexes(
        0, 0, len(summary_rows) + 1, len(summary_headers)
    ).values = [summary_headers] + summary_rows

    details.get_range_by_indexes(
        0, 0, len(detail_rows) + 1, len(detail_headers)
    ).values = [detail_headers] + detail_rows

    ground.get_range_by_indexes(
        0, 0, len(ground_rows) + 1, len(ground_headers)
    ).values = [ground_headers] + ground_rows

    metrics.get_range("A1:B10").values = [
        ["Metric", "Value"],
        ["TP", None],
        ["FP", None],
        ["FN", None],
        ["TN", None],
        ["Sensitivity", None],
        ["Specificity", None],
        ["Precision", None],
        ["Accuracy", None],
        ["F1 Score", None],
    ]

    last_summary_row = len(summary_rows) + 1
    metrics.get_range("B2:B10").formulas = [
        [f"=SUM(Summary!E2:E{last_summary_row})"],
        [f"=SUM(Summary!F2:F{last_summary_row})"],
        [f"=SUM(Summary!G2:G{last_summary_row})"],
        [f"=SUM(Summary!H2:H{last_summary_row})"],
        ["=IFERROR(B2/(B2+B4),0)"],
        ["=IFERROR(B5/(B5+B3),0)"],
        ["=IFERROR(B2/(B2+B3),0)"],
        ["=IFERROR((B2+B5)/(B2+B3+B4+B5),0)"],
        ["=IFERROR((2*B2)/(2*B2+B3+B4),0)"],
    ]

    border = {
        "top": {"style": "thin", "color": "#D9D9D9"},
        "bottom": {"style": "thin", "color": "#D9D9D9"},
        "left": {"style": "thin", "color": "#D9D9D9"},
        "right": {"style": "thin", "color": "#D9D9D9"},
    }

    for sheet, row_count, col_count in [
        (summary, len(summary_rows) + 1, len(summary_headers)),
        (details, len(detail_rows) + 1, len(detail_headers)),
        (ground, len(ground_rows) + 1, len(ground_headers)),
        (metrics, 10, 2),
    ]:
        sheet.get_range_by_indexes(0, 0, row_count, col_count).format = {
            "font": {"name": "Calibri", "size": 11},
            "borders": border,
            "wrap_text": True,
            "vertical_alignment": "center",
        }
        sheet.get_range_by_indexes(0, 0, 1, col_count).format = {
            "font": {"name": "Calibri", "size": 11, "bold": True},
            "borders": border,
            "wrap_text": True,
            "vertical_alignment": "center",
        }
        sheet.freeze_panes.freeze_rows(1)

    for column, width in {
        "A": 18, "B": 20, "C": 28, "D": 28,
        "E": 9, "F": 9, "G": 9, "H": 9,
    }.items():
        summary.get_range(f"{column}:{column}").format.column_width = width

    for column, width in {
        "A": 18, "B": 20, "C": 18, "D": 12, "E": 12,
        "F": 22, "G": 15, "H": 12, "I": 16, "J": 26, "K": 15,
    }.items():
        details.get_range(f"{column}:{column}").format.column_width = width

    for column, width in {
        "A": 20, "B": 18, "C": 42, "D": 26,
        "E": 28, "F": 24, "G": 14, "H": 22,
    }.items():
        ground.get_range(f"{column}:{column}").format.column_width = width

    metrics.get_range("A:A").format.column_width = 20
    metrics.get_range("B:B").format.column_width = 16
    metrics.get_range("B6:B10").format.number_format = "0.00%"

    if detail_rows:
        details.get_range(
            f"H2:H{len(detail_rows) + 1}"
        ).format.number_format = "0.00"

    SpreadsheetFile.export_xlsx(workbook).save(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Compare COMPLETE NCBI StxTyper subtype calls with the shared ground truth."
    )
    parser.add_argument("--groundtruth", required=True)
    parser.add_argument("--stxtyper", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ground_truth = load_ground_truth(args.groundtruth)
    calls_by_genome = load_complete_stxtyper_calls(args.stxtyper)
    build_workbook(ground_truth, calls_by_genome, args.output)


if __name__ == "__main__":
    main()
