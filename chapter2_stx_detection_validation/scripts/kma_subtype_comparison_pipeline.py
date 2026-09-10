#!/usr/bin/env python3
"""
KMA StxDB subtype-level validation pipeline 
================================================================

Project
-------
MSc Thesis - StxDB validation

Purpose
-------
This script compares KMA subtype-level calls against a collapsed ground truth
Excel table and generates Excel workbooks that summarize TP, FP, FN, TN, and
performance metrics for each KMA parameter set.

This annotated version is intended to be easier to read, reuse, and cite in a
reproducible analysis workflow.

Inputs
------
1. Ground truth Excel workbook
   Example: StxDB_groundtruth_subtype_level.xlsx

2. One or more KMA TSV result files
   Examples:
       kmaresults_70_conclave1.tsv
       kmaresults_70_conclave2.tsv
       kmaresults_75_conclave1.tsv
       kmaresults_75_conclave2.tsv
       kmaresults_80_conclave1.tsv
       kmaresults_80_conclave2.tsv
       kmaresults_85_conclave1.tsv
       kmaresults_85_conclave2.tsv
       kmaresults_90_conclave1.tsv
       kmaresults_90_conclave2.tsv
       kmaresults_95_conclave1.tsv
       kmaresults_95_conclave2.tsv

Outputs
-------
For each KMA TSV file, the script creates one Excel workbook containing:

1. Summary
   One row per genome with expected subtypes, detected subtypes, and counts for
   TP, FP, FN, and TN.

2. Subtype_Details
   One row per TP, FP, or FN subtype-level event.
   TN rows are intentionally not listed because TN is calculated as:

       TN = 19 - TP - FP - FN

3. Ground_Truth
   Copy of the collapsed ground truth table used for comparison.

4. Metrics
   Overall TP, FP, FN, TN, Sensitivity, Specificity, Precision, Accuracy,
   and F1 Score for that parameter set.

The script also creates one combined workbook:

    KMA_parameter_metrics_comparison.xlsx

This workbook has one row per KMA parameter set so the results can be compared
more easily.

Classification definitions
--------------------------
TP = subtype reported by KMA and present in the ground truth.
FP = subtype reported by KMA but absent from the ground truth.
FN = subtype present in the ground truth but not reported by KMA.
TN = 19 - TP - FP - FN, calculated per genome.

Important assumption
--------------------
Template identity is retained for reporting only. It is NOT used to decide
whether a result is TP, FP, FN, or TN.

Dependencies
------------
Install required Python packages with:

    pip install pandas openpyxl

Example command
---------------
python kma_subtype_comparison_pipeline.py \
  --groundtruth StxDB_groundtruth_subtype_level.xlsx \
  --kma kmaresults_70_conclave1.tsv kmaresults_70_conclave2.tsv \
        kmaresults_75_conclave1.tsv kmaresults_75_conclave2.tsv \
        kmaresults_80_conclave1.tsv kmaresults_80_conclave2.tsv \
        kmaresults_85_conclave1.tsv kmaresults_85_conclave2.tsv \
        kmaresults_90_conclave1.tsv kmaresults_90_conclave2.tsv \
        kmaresults_95_conclave1.tsv kmaresults_95_conclave2.tsv \
  --outdir kma_comparison_outputs
"""

from __future__ import annotations

###############################################################################
# Import Python libraries
###############################################################################

# argparse is used to let the script accept command-line arguments.
# This makes it easier to run the same script on different input files.
import argparse

# re is used for regular expressions, mainly to extract subtype names such as
# Stx1a or Stx2i from longer database/template names.
import re

# pathlib provides a safer and cleaner way to handle file paths.
from pathlib import Path

# These typing imports are not required for the script to run, but they make the
# expected inputs and outputs of functions clearer.
from typing import Dict, List, Optional, Sequence, Set, Tuple

# pandas is used to read, manipulate, and write tabular data.
import pandas as pd

# openpyxl is used for formatting the Excel workbooks after pandas writes them.
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

###############################################################################
# Global settings
###############################################################################

# There are 19 possible Stx subtypes in the comparison space.
# TN is calculated per genome as:
#
#     TN = 19 - TP - FP - FN
#
# This means TNs are not listed as individual rows in Subtype_Details, but they
# are still included in the Summary and Metrics sheets.
TOTAL_STX_SUBTYPES = 19

# This list is used only to make the Excel output easier to read.
# It controls the order in which subtypes appear in tables.
# It does NOT affect classification.
KNOWN_SUBTYPE_ORDER = [
    "Stx1a", "Stx1c", "Stx1d", "Stx1e",
    "Stx2a", "Stx2b", "Stx2c", "Stx2d", "Stx2e", "Stx2f", "Stx2g",
    "Stx2h", "Stx2i", "Stx2j", "Stx2k", "Stx2l", "Stx2m", "Stx2n", "Stx2o",
]

# Subtypes listed here are present as stx-like targets in the validation table,
# but they are not counted as functional expected toxin targets for validation.
# If a tool detects one of these subtypes in the corresponding genome, the call
# is classified as FP, not TP. If the tool misses it, it is not counted as FN.
NONFUNCTIONAL_SUBTYPES_BY_GENOME = {
    # GCA_018966965.1 / 2025-SEQ-1796 has functional Stx1a, but interrupted Stx2a.
    "2025-SEQ-1796": {"Stx2a"},

    # GCA_026134865.1 / SRR18191635 has functional Stx2d, but frameshifted Stx2b.
    "SRR18191635": {"Stx2b"},
}

# Regular expression used to find subtype strings.
# Examples it can detect:
#   Stx1a
#   stx1a
#   Stx2i
#   stx2i
# The script later normalizes the case to Stx1a, Stx2i, etc.
SUBTYPE_RE = re.compile(r"Stx[12][a-z]", re.IGNORECASE)

###############################################################################
# Helper functions for parsing and formatting values
###############################################################################

def split_multi_value_cell(value) -> List[str]:
    """
    Split cells that contain multiple values.

    Some ground truth cells may contain multiple operons or subtypes separated
    by semicolons, commas, or new lines. This function splits those cells into
    a clean Python list.

    Parameters
    ----------
    value : any
        Cell value from the Excel file.

    Returns
    -------
    list of str
        Clean list of values with blank entries removed.
    """
    if pd.isna(value):
        return []

    # Split on semicolon, newline, or comma.
    parts = re.split(r"[;\n,]+", str(value))

    # Remove blank entries and extra spaces.
    return [p.strip() for p in parts if p and p.strip()]


def normalize_subtype(value: str) -> Optional[str]:
    """
    Extract and normalize a subtype name.

    This function finds subtype text inside a larger string and normalizes its
    capitalization.

    Examples
    --------
    "StxOp1339_Stx2b_109" -> "Stx2b"
    "stx2i-O9-CB10366"    -> "Stx2i"

    Parameters
    ----------
    value : str
        Text that may contain a subtype.

    Returns
    -------
    str or None
        Normalized subtype, or None if no subtype was found.
    """
    if value is None or pd.isna(value):
        return None

    match = SUBTYPE_RE.search(str(value))
    if not match:
        return None

    raw = match.group(0)

    # Convert to canonical format:
    # Stx + toxin number + lowercase subtype letter.
    # Example: stx2I -> Stx2i
    return "Stx" + raw[3].upper() + raw[4].lower()


def subtype_sort_key(subtype: str) -> Tuple[int, int, str]:
    """
    Sort subtypes in a biologically readable order.

    Known Stx subtypes are ordered according to KNOWN_SUBTYPE_ORDER. Any unknown
    subtype is placed after the known subtypes.
    """
    if subtype in KNOWN_SUBTYPE_ORDER:
        return (0, KNOWN_SUBTYPE_ORDER.index(subtype), subtype)
    return (1, 999, subtype)


def parse_parameter_label(path: Path) -> str:
    """
    Create a readable parameter label from a KMA result filename.

    Example
    -------
    kmaresults_85_conclave1.tsv -> "85% ID, conclave 1"
    """
    stem = path.stem
    match = re.search(r"(\d+).*?conclave[_-]?(\d+)", stem, flags=re.IGNORECASE)
    if match:
        return f"{match.group(1)}% ID, conclave {match.group(2)}"
    return stem


def output_name_for_kma(path: Path) -> str:
    """
    Create the output workbook name for one KMA result file.

    Example
    -------
    kmaresults_85_conclave1.tsv -> KMA_85_conclave1_noTN.xlsx
    """
    stem = path.stem
    match = re.search(r"(\d+).*?conclave[_-]?(\d+)", stem, flags=re.IGNORECASE)
    if match:
        return f"KMA_{match.group(1)}_conclave{match.group(2)}_noTN.xlsx"

    # Fallback name if the file does not follow the expected naming pattern.
    safe = re.sub(r"[^A-Za-z0-9_\-]+", "_", stem)
    return f"{safe}_comparison_noTN.xlsx"

###############################################################################
# Ground truth parsing
###############################################################################

def collapse_groundtruth(gt_path: Path) -> pd.DataFrame:
    """
    Read and collapse the subtype-level ground truth workbook.

    The output table has one row per genome and contains:
    - Genome
    - Expected_Operons
    - Expected_Subtypes, meaning functional expected subtypes used for TP/FN calls
    - Nonfunctional_Subtypes, meaning interrupted/frameshifted subtypes that become FP if detected
    - All_Annotated_Subtypes, meaning all subtypes retained for transparency
    - Operon_Count
    - Unique_Subtype_Count

    The expected subtypes are converted into normalized names such as Stx1a,
    Stx2a, Stx2i, etc.

    Parameters
    ----------
    gt_path : pathlib.Path
        Path to the ground truth Excel workbook.

    Returns
    -------
    pandas.DataFrame
        Collapsed ground truth table with one row per genome.
    """
    gt = pd.read_excel(gt_path)

    # Flexible column detection:
    # The files used in this workflow usually have a Fastq column. If not,
    # the script uses the first column as the genome/sample identifier.
    fastq_col = "Fastq" if "Fastq" in gt.columns else gt.columns[0]

    # Expected operons are retained for reporting only.
    # They are not used for subtype-level classification.
    operon_col = "Expected_Operons" if "Expected_Operons" in gt.columns else None

    # Expected subtypes are required for this subtype-level comparison.
    subtype_col = "Expected_Subtypes" if "Expected_Subtypes" in gt.columns else None

    # If the exact Expected_Subtypes column name is not present, try to find a
    # column containing the word "subtype".
    if subtype_col is None:
        matches = [c for c in gt.columns if "subtype" in str(c).lower()]
        if not matches:
            raise ValueError(
                "Could not find an Expected_Subtypes column in the ground truth table."
            )
        subtype_col = matches[0]

    rows = []

    for _, row in gt.iterrows():
        genome = str(row[fastq_col]).strip()

        # Skip blank genome/sample names.
        if not genome or genome.lower() == "nan":
            continue

        # Operons are included in the output workbook for reference.
        expected_operons = split_multi_value_cell(row[operon_col]) if operon_col else []

        # Extract and normalize all annotated subtypes from the ground truth.
        # These include both functional and non-functional stx-like targets.
        all_annotated_subtypes = []
        for item in split_multi_value_cell(row[subtype_col]):
            subtype = normalize_subtype(item)
            if subtype and subtype not in all_annotated_subtypes:
                all_annotated_subtypes.append(subtype)

        # Apply the functional-toxin validation rule.
        # Non-functional subtypes remain visible in the output workbook, but are
        # removed from the expected positive set used for TP/FN calls.
        nonfunctional_set = NONFUNCTIONAL_SUBTYPES_BY_GENOME.get(genome, set())
        nonfunctional_subtypes = [
            subtype for subtype in all_annotated_subtypes if subtype in nonfunctional_set
        ]
        expected_subtypes = [
            subtype for subtype in all_annotated_subtypes if subtype not in nonfunctional_set
        ]

        rows.append(
            {
                "Genome": genome,
                "Expected_Operons": "; ".join(expected_operons),
                "Expected_Subtypes": "; ".join(expected_subtypes),
                "Nonfunctional_Subtypes": "; ".join(nonfunctional_subtypes),
                "All_Annotated_Subtypes": "; ".join(all_annotated_subtypes),
                "Operon_Count": len(expected_operons),
                "Unique_Subtype_Count": len(expected_subtypes),
            }
        )

    collapsed = pd.DataFrame(rows)

    # Keep one row per genome.
    collapsed = collapsed.drop_duplicates(subset=["Genome"], keep="first")

    # Sort output to make manual checking easier.
    collapsed = collapsed.sort_values("Genome").reset_index(drop=True)

    return collapsed

###############################################################################
# KMA result parsing
###############################################################################

def parse_kma(kma_path: Path) -> pd.DataFrame:
    """
    Read one KMA TSV result file.

    The script expects the KMA file to contain at least:
    - Sample
    - Template

    If Template_Identity is present, it is retained for reporting. If it is
    missing, a blank column is created.

    Classification does NOT depend on Template_Identity.

    Parameters
    ----------
    kma_path : pathlib.Path
        Path to one KMA TSV file.

    Returns
    -------
    pandas.DataFrame
        KMA results with an added Detected_Subtype column.
    """
    kma = pd.read_csv(kma_path, sep="\t")

    required_columns = {"Sample", "Template"}
    missing_columns = required_columns - set(kma.columns)
    if missing_columns:
        raise ValueError(
            f"{kma_path} is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}"
        )

    # Template identity is kept only for reporting.
    if "Template_Identity" not in kma.columns:
        kma["Template_Identity"] = pd.NA

    # Extract subtype from the KMA template name.
    kma["Detected_Subtype"] = kma["Template"].map(normalize_subtype)

    # Remove rows where no subtype could be extracted.
    kma = kma.dropna(subset=["Detected_Subtype"]).copy()

    # Clean sample names.
    kma["Sample"] = kma["Sample"].astype(str).str.strip()

    return kma


def top_hit_for_subtype(kma_sample: pd.DataFrame, subtype: str) -> Tuple[Optional[str], Optional[float]]:
    """
    Select the displayed top KMA hit for a subtype in one sample.

    This is only for reporting in the Subtype_Details sheet. If multiple KMA
    templates correspond to the same subtype, the template with the highest
    Template_Identity is displayed.

    Important:
    Template identity is NOT used to classify TP, FP, or FN.

    Parameters
    ----------
    kma_sample : pandas.DataFrame
        KMA rows for a single genome/sample.

    subtype : str
        Subtype being compared.

    Returns
    -------
    tuple
        (top template name, template identity)
    """
    hits = kma_sample[kma_sample["Detected_Subtype"] == subtype].copy()

    if hits.empty:
        return None, None

    # Convert identity to numeric so it can be sorted.
    hits["_identity_numeric"] = pd.to_numeric(hits["Template_Identity"], errors="coerce")

    # Highest identity first. Template name is used as a secondary sort key for
    # reproducible output if identities are tied.
    hits = hits.sort_values(["_identity_numeric", "Template"], ascending=[False, True])

    top = hits.iloc[0]
    identity = top.get("Template_Identity")

    if pd.isna(identity):
        identity = None

    return str(top["Template"]), identity

###############################################################################
# Core comparison logic
###############################################################################

def compare_one_kma(kma_path: Path, groundtruth: pd.DataFrame, outdir: Path) -> Dict[str, object]:
    """
    Compare one KMA result file to the ground truth table.

    For every genome:
    - expected_subtypes = subtypes listed in the ground truth
    - detected_subtypes = subtypes reported by KMA

    Classification:
    - TP = detected subtype is expected
    - FP = detected subtype is not expected
    - FN = expected subtype is not detected
    - TN = 19 - TP - FP - FN

    The Subtype_Details sheet only includes TP, FP, and FN rows.
    TN rows are intentionally excluded because they are implied by the formula.

    Parameters
    ----------
    kma_path : pathlib.Path
        Path to one KMA TSV file.

    groundtruth : pandas.DataFrame
        Collapsed ground truth table from collapse_groundtruth().

    outdir : pathlib.Path
        Directory where the output workbook should be saved.

    Returns
    -------
    dict
        Overall metric values for the parameter set. This dictionary is later
        used to build the combined metrics workbook.
    """
    kma = parse_kma(kma_path)

    # Group KMA rows by sample name for faster lookup during comparison.
    kma_by_sample = {sample: df for sample, df in kma.groupby("Sample")}

    summary_rows = []
    detail_rows = []

    for _, gt_row in groundtruth.iterrows():
        genome = gt_row["Genome"]

        # Functional expected subtype set for this genome.
        # Only these subtypes can become TP if detected or FN if missed.
        expected_subtypes: Set[str] = set(split_multi_value_cell(gt_row["Expected_Subtypes"]))

        # Non-functional/interrupted subtype set for this genome.
        # If detected, these are counted as FP; if missed, they are not FN.
        nonfunctional_subtypes: Set[str] = set(
            split_multi_value_cell(gt_row.get("Nonfunctional_Subtypes", ""))
        )

        # KMA hits for this genome. If the genome has no KMA hits, use an empty
        # DataFrame with the same columns.
        sample_hits = kma_by_sample.get(genome, pd.DataFrame(columns=kma.columns))

        # Detected subtype set for this genome.
        if not sample_hits.empty:
            detected_subtypes: Set[str] = set(sample_hits["Detected_Subtype"].dropna().unique())
        else:
            detected_subtypes = set()

        #######################################################################
        # Count TP, FP, FN, and TN for this genome
        #######################################################################

        # TP: expected and detected.
        tp = len(expected_subtypes & detected_subtypes)

        # FP: detected but not expected.
        fp = len(detected_subtypes - expected_subtypes)

        # FN: expected but not detected.
        fn = len(expected_subtypes - detected_subtypes)

        # TN: all remaining subtype categories in the 19-subtype comparison
        # space. TN is counted but not listed in Subtype_Details.
        tn = TOTAL_STX_SUBTYPES - tp - fp - fn

        #######################################################################
        # Add one row to the Summary sheet
        #######################################################################
        summary_rows.append(
            {
                "Genome": genome,
                "Expected operons": gt_row["Expected_Operons"],
                "Expected functional subtypes": "; ".join(sorted(expected_subtypes, key=subtype_sort_key)),
                "Nonfunctional subtypes": "; ".join(sorted(nonfunctional_subtypes, key=subtype_sort_key)),
                "Detected subtypes": "; ".join(sorted(detected_subtypes, key=subtype_sort_key)),
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "TN": tn,
            }
        )

        #######################################################################
        # Add rows to the Subtype_Details sheet
        #######################################################################
        # Only the union of expected and detected subtypes is reviewed.
        # This means the details sheet includes only TP, FP, and FN rows.
        # It does not include TN rows.
        for subtype in sorted(expected_subtypes | nonfunctional_subtypes | detected_subtypes, key=subtype_sort_key):
            expected = subtype in expected_subtypes
            nonfunctional = subtype in nonfunctional_subtypes
            detected = subtype in detected_subtypes

            if expected and detected:
                classification = "TP"
            elif nonfunctional and detected:
                classification = "FP"
            elif detected and not expected:
                classification = "FP"
            elif expected and not detected:
                classification = "FN"
            else:
                # This would be a TN, but TN rows are not included in the
                # Subtype_Details sheet.
                continue

            top_hit, identity = top_hit_for_subtype(sample_hits, subtype)

            detail_rows.append(
                {
                    "Genome": genome,
                    "Compared subtype": subtype,
                    "Expected functional?": "Yes" if expected else "No",
                    "Nonfunctional?": "Yes" if nonfunctional else "No",
                    "Detected?": "Yes" if detected else "No",
                    "Top KMA hit": top_hit if top_hit is not None else "",
                    "Template identity": identity if identity is not None else "",
                    "Classification": classification,
                }
            )

    summary = pd.DataFrame(summary_rows)
    details = pd.DataFrame(detail_rows)

    ###########################################################################
    # Calculate overall counts for the full parameter set
    ###########################################################################
    tp_total = int(summary["TP"].sum())
    fp_total = int(summary["FP"].sum())
    fn_total = int(summary["FN"].sum())
    tn_total = int(summary["TN"].sum())

    ###########################################################################
    # Calculate performance metrics
    ###########################################################################
    # Sensitivity / Recall:
    #     TP / (TP + FN)
    sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) else None

    # Specificity:
    #     TN / (TN + FP)
    specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) else None

    # Precision / Positive predictive value:
    #     TP / (TP + FP)
    precision = tp_total / (tp_total + fp_total) if (tp_total + fp_total) else None

    # Accuracy:
    #     (TP + TN) / (TP + FP + FN + TN)
    denominator = tp_total + fp_total + fn_total + tn_total
    accuracy = (tp_total + tn_total) / denominator if denominator else None

    # F1 Score:
    #     2 * Precision * Sensitivity / (Precision + Sensitivity)
    if precision is not None and sensitivity is not None and (precision + sensitivity):
        f1 = 2 * precision * sensitivity / (precision + sensitivity)
    else:
        f1 = None

    metrics = pd.DataFrame(
        {
            "Metric": [
                "TP",
                "FP",
                "FN",
                "TN",
                "Sensitivity",
                "Specificity",
                "Precision",
                "Accuracy",
                "F1 Score",
            ],
            "Value": [
                tp_total,
                fp_total,
                fn_total,
                tn_total,
                sensitivity,
                specificity,
                precision,
                accuracy,
                f1,
            ],
        }
    )

    ###########################################################################
    # Write the per-parameter Excel workbook
    ###########################################################################
    out_path = outdir / output_name_for_kma(kma_path)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        details.to_excel(writer, sheet_name="Subtype_Details", index=False)
        groundtruth.to_excel(writer, sheet_name="Ground_Truth", index=False)
        metrics.to_excel(writer, sheet_name="Metrics", index=False)

    # Apply formatting after writing the data.
    format_workbook(out_path)

    # Return the overall metrics so they can be added to the combined summary.
    return {
        "Parameter_Set": parse_parameter_label(kma_path),
        "Output_Workbook": out_path.name,
        "TP": tp_total,
        "FP": fp_total,
        "FN": fn_total,
        "TN": tn_total,
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Precision": precision,
        "Accuracy": accuracy,
        "F1 Score": f1,
    }

###############################################################################
# Excel formatting
###############################################################################

def format_workbook(path: Path) -> None:
    """
    Apply basic formatting to an Excel workbook.

    Formatting includes:
    - colored header row
    - bold white header font
    - borders
    - wrapped text
    - frozen header row
    - adjusted column widths
    - percentage formatting for metric values

    Parameters
    ----------
    path : pathlib.Path
        Path to the Excel workbook to format.
    """
    wb = load_workbook(path)

    header_fill = PatternFill("solid", fgColor="0F766E")
    header_font = Font(bold=True, color="FFFFFF")
    thin = Side(style="thin", color="D9E2E7")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for ws in wb.worksheets:
        # Freeze the first row so headers stay visible while scrolling.
        ws.freeze_panes = "A2"

        # Format the header row.
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

        # Format all cells with borders and wrapped text.
        for row in ws.iter_rows():
            for cell in row:
                cell.border = border
                cell.alignment = Alignment(vertical="top", wrap_text=True)

        # Adjust column widths based on the longest value in each column.
        for col_cells in ws.columns:
            letter = get_column_letter(col_cells[0].column)
            max_len = max(len(str(c.value)) if c.value is not None else 0 for c in col_cells)
            ws.column_dimensions[letter].width = min(max(max_len + 2, 10), 45)

        # Format metric values as percentages where appropriate.
        if ws.title == "Metrics":
            for row in range(2, ws.max_row + 1):
                metric = ws.cell(row=row, column=1).value
                if metric in {"Sensitivity", "Specificity", "Precision", "Accuracy", "F1 Score"}:
                    ws.cell(row=row, column=2).number_format = "0.00%"

        if ws.title == "Metrics_Comparison":
            headers = [cell.value for cell in ws[1]]
            for col_name in ["Sensitivity", "Specificity", "Precision", "Accuracy", "F1 Score"]:
                if col_name in headers:
                    col_idx = headers.index(col_name) + 1
                    for row in range(2, ws.max_row + 1):
                        ws.cell(row=row, column=col_idx).number_format = "0.00%"

    wb.save(path)

###############################################################################
# Combined metrics workbook
###############################################################################

def write_combined_metrics(metrics_rows: Sequence[Dict[str, object]], outdir: Path) -> Path:
    """
    Create one Excel workbook comparing all KMA parameter sets.

    The combined workbook contains:

    1. Metrics_Comparison
       One row per KMA parameter set.

    2. Definitions
       Definitions of TP, FP, FN, TN, and metric formulas.

    Parameters
    ----------
    metrics_rows : sequence of dict
        Metric dictionaries returned by compare_one_kma().

    outdir : pathlib.Path
        Directory where the combined workbook should be saved.

    Returns
    -------
    pathlib.Path
        Path to the combined metrics workbook.
    """
    combined = pd.DataFrame(metrics_rows)

    columns = [
        "Parameter_Set",
        "Output_Workbook",
        "TP",
        "FP",
        "FN",
        "TN",
        "Sensitivity",
        "Specificity",
        "Precision",
        "Accuracy",
        "F1 Score",
    ]
    combined = combined[columns]

    definitions = pd.DataFrame(
        {
            "Term": [
                "TP",
                "FP",
                "FN",
                "TN",
                "Sensitivity",
                "Specificity",
                "Precision",
                "Accuracy",
                "F1 Score",
                "Template identity",
                "Subtype_Details TN rows",
            ],
            "Definition / Formula": [
                "Subtype reported by KMA and present as a functional expected subtype in the ground truth.",
                "Subtype reported by KMA but absent from the functional expected subtype set, including interrupted or frameshifted non-functional targets.",
                "Functional expected subtype present in the ground truth but not reported by KMA. Non-functional targets are not counted as FN if missed.",
                "19 - TP - FP - FN per genome; summed across genomes for each parameter set.",
                "TP / (TP + FN)",
                "TN / (TN + FP)",
                "TP / (TP + FP)",
                "(TP + TN) / (TP + FP + FN + TN)",
                "2 × (Precision × Sensitivity) / (Precision + Sensitivity)",
                "Retained for reporting only; not used for TP/FP/FN classification.",
                "TN rows are not included in Subtype_Details because TN is calculated from the 19-subtype comparison space.",
            ],
        }
    )

    out_path = outdir / "KMA_parameter_metrics_comparison.xlsx"

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        combined.to_excel(writer, sheet_name="Metrics_Comparison", index=False)
        definitions.to_excel(writer, sheet_name="Definitions", index=False)

    format_workbook(out_path)

    return out_path

###############################################################################
# Command-line interface and workflow
###############################################################################

def main() -> None:
    """
    Run the full KMA subtype comparison workflow.

    Workflow
    --------
    1. Read command-line arguments.
    2. Create the output directory if it does not already exist.
    3. Load and collapse the ground truth table.
    4. For each KMA TSV file:
       - parse KMA results
       - compare detected subtypes to expected subtypes
       - write one per-parameter Excel workbook
       - collect overall metrics
    5. Write the combined metrics comparison workbook.
    """
    parser = argparse.ArgumentParser(
        description="Compare KMA subtype results to collapsed StxDB ground truth."
    )

    parser.add_argument(
        "--groundtruth",
        required=True,
        type=Path,
        help="Path to collapsed ground truth Excel workbook.",
    )

    parser.add_argument(
        "--kma",
        required=True,
        nargs="+",
        type=Path,
        help="One or more KMA result TSV files.",
    )

    parser.add_argument(
        "--outdir",
        default=Path("kma_comparison_outputs"),
        type=Path,
        help="Directory where output Excel workbooks will be saved.",
    )

    args = parser.parse_args()

    # Create output directory if needed.
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load and collapse ground truth once. The same ground truth is then used
    # for all KMA parameter sets.
    groundtruth = collapse_groundtruth(args.groundtruth)

    metrics_rows = []

    # Process each KMA parameter set.
    for kma_path in args.kma:
        print(f"Processing {kma_path}...")
        metrics_rows.append(compare_one_kma(kma_path, groundtruth, args.outdir))

    # Create one workbook that compares all parameter sets.
    combined_path = write_combined_metrics(metrics_rows, args.outdir)

    print("\nDone.")
    print(f"Combined metrics workbook: {combined_path}")
    print(f"Per-parameter workbooks saved in: {args.outdir}")


if __name__ == "__main__":
    main()
