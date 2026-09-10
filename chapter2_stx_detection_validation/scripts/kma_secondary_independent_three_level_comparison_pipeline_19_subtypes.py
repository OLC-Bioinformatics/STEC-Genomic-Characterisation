#!/usr/bin/env python3
"""
kma_secondary_independent_three_level_comparison.py

Secondary KMA validation analysis for StxDB.

This script compares KMA-reported Stx operon hits to the curated validation
truth table at three independent, non-hierarchical levels:

1) Operon level  = full StxOp ID, e.g. StxOp933_Stx2i_43
2) Variant level = subtype + allele, e.g. Stx2i_43
3) Subtype level = subtype only, e.g. Stx2i or Stx2o_like

This is NOT the older greedy/hierarchical one-to-one matching approach.
For each genome and each level, the script does a set comparison:

    TP = expected target detected by KMA
    FP = KMA target detected but not expected
    FN = expected target not detected by KMA
    TN = database universe size at that level - TP - FP - FN

The universe sizes are calculated directly from the Stx operon FASTA database:

    operon universe  = all unique full FASTA headers
    variant universe = all unique subtype+allele strings after removing StxOp###_
    subtype universe = the fixed set of 19 currently recognized Stx subtypes

Example:
    StxOp1352_Stx2o_like_1352
        operon  = StxOp1352_Stx2o_like_1352
        variant = Stx2o_like_1352
        subtype = Stx2o_like

Usage example:

python kma_secondary_independent_three_level_comparison_pipeline.py \
  --truth-xlsx 260609_StxDB_Validation.xlsx \
  --truth-sheet Final_Dataset \
  --db-fasta Stx_operon_20251114.fasta \
  --kma-files kmaresults_70_conclave1.tsv \
  --output-dir kma_secondary_outputs

Inputs expected:
- Truth Excel sheet must contain columns: Fastq, Fasta, StxOpDB
- KMA TSV files must contain columns: Sample and Template
  The Template column should contain IDs like StxOp933_Stx2i_43.

Dependencies:
    pandas, openpyxl
"""

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Set

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

LEVELS = ["operon", "variant", "subtype"]
STXOP_RE = re.compile(r"^(StxOp\d+)_(.+)$")

# Fixed subtype universe used for subtype-level TN calculations.
# These are the 19 currently recognized Stx subtypes represented in the
# benchmark analysis. Database labels describing chimeras, combined groups,
# or unusual aliases are not treated as additional recognized subtypes.
RECOGNIZED_STX_SUBTYPES = {
    "Stx1a", "Stx1c", "Stx1d", "Stx1e",
    "Stx2a", "Stx2b", "Stx2c", "Stx2d", "Stx2e", "Stx2f",
    "Stx2g", "Stx2h", "Stx2i", "Stx2j", "Stx2k", "Stx2l",
    "Stx2m", "Stx2n", "Stx2o",
}

# Special validation cases: these stx targets are present in the validation
# table, but are not considered functional toxin targets. Therefore, they are
# excluded from the expected-positive set. If KMA detects them, they are
# classified as FP rather than TP. If KMA does not detect them, they are not FN.
NONFUNCTIONAL_TARGETS_BY_SAMPLE = {
    "2025-SEQ-1796": {
        "operon": {"StxOp1089_Stx2a_1"},
        "variant": {"Stx2a_1"},
        "subtype": {"Stx2a"},
    },
    "SRR18191635": {
        "operon": {"StxOp3298_Stx2b_23"},
        "variant": {"Stx2b_23"},
        "subtype": {"Stx2b"},
    },
}


def parse_stx_id(value: object) -> Optional[Dict[str, str]]:
    """Parse a StxDB identifier into operon, variant, and subtype levels.

    This parser deliberately does not use cut -d'_' -f2 for subtype,
    because some subtype labels contain internal underscores, e.g. Stx2o_like.

    Instead:
    - remove only the StxOp###_ prefix to get the variant
    - remove only the final _allele field to get the subtype
    """
    if pd.isna(value):
        return None

    text = str(value).strip().lstrip(">")
    if not text:
        return None

    # If a FASTA header has extra text after the first whitespace, ignore it.
    text = text.split()[0]

    match = STXOP_RE.match(text)
    if not match:
        return None

    operon = text
    variant = match.group(2)
    subtype = re.sub(r"_[^_]+$", "", variant)

    return {
        "operon": operon,
        "variant": variant,
        "subtype": subtype,
    }


def load_database_universes(db_fasta: Path) -> Dict[str, Set[str]]:
    """Read the StxDB FASTA and define the universe for each comparison level.

    Operon and variant universes are derived directly from all unique database
    identifiers. The subtype universe is fixed at the 19 currently recognized
    Stx subtypes so that unusual database labels, combined subtype groups, and
    chimeric designations do not inflate subtype-level TN calculations.
    """
    universes = {level: set() for level in LEVELS}

    with db_fasta.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            parsed = parse_stx_id(line)
            if parsed is None:
                continue

            universes["operon"].add(parsed["operon"])
            universes["variant"].add(parsed["variant"])

    universes["subtype"] = set(RECOGNIZED_STX_SUBTYPES)

    return universes


def load_truth(truth_xlsx: Path, truth_sheet: str) -> Dict[str, dict]:
    """Load expected StxDB targets from the validation table.

    The sample key is the Fastq column because KMA was run on raw read IDs.
    If Fastq is blank for any row, the script falls back to Fasta.
    """
    df = pd.read_excel(truth_xlsx, sheet_name=truth_sheet, dtype=str)

    required = {"Fasta", "Fastq", "StxOpDB"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Truth table is missing columns: {sorted(missing)}")

    truth = defaultdict(lambda: {
        "Fasta": "",
        "Fastq": "",
        "operon": set(),
        "variant": set(),
        "subtype": set(),
    })

    for _, row in df.iterrows():
        parsed = parse_stx_id(row.get("StxOpDB"))
        if parsed is None:
            continue

        sample = str(row.get("Fastq") or row.get("Fasta") or "").strip()
        if not sample or sample.lower() == "nan":
            continue

        truth[sample]["Fastq"] = str(row.get("Fastq") or "").strip()
        truth[sample]["Fasta"] = str(row.get("Fasta") or "").strip()

        for level in LEVELS:
            truth[sample][level].add(parsed[level])

    return truth


def load_kma_hits(kma_tsv: Path) -> Dict[str, dict]:
    """Load detected KMA targets from one combined KMA TSV file."""
    df = pd.read_csv(kma_tsv, sep="\t", dtype=str)

    # Some older combined KMA files used File instead of Sample.
    sample_col = "Sample" if "Sample" in df.columns else "File"
    if sample_col not in df.columns or "Template" not in df.columns:
        raise ValueError(f"KMA file must contain Sample/File and Template columns: {kma_tsv}")

    hits = defaultdict(lambda: {
        "operon": set(),
        "variant": set(),
        "subtype": set(),
    })

    for _, row in df.iterrows():
        sample = str(row.get(sample_col) or "").strip().replace(".res", "")
        parsed = parse_stx_id(row.get("Template"))

        if not sample or parsed is None:
            continue

        for level in LEVELS:
            hits[sample][level].add(parsed[level])

    return hits


def safe_divide(numerator: int, denominator: int) -> Optional[float]:
    return numerator / denominator if denominator else None


def compare_one_kma_file(kma_tsv: Path, truth: Dict[str, dict], universes: Dict[str, Set[str]]):
    """Perform independent operon, variant, and subtype comparisons for one KMA file."""
    kma_hits = load_kma_hits(kma_tsv)
    samples = sorted(set(truth) | set(kma_hits))
    parameter_set = kma_tsv.stem

    summary_rows = []
    detail_rows = {level: [] for level in LEVELS}
    totals = {level: Counter() for level in LEVELS}

    for sample in samples:
        meta = truth.get(sample, {})
        base = {
            "Parameter_Set": parameter_set,
            "Genome": sample,
            "Fasta": meta.get("Fasta", ""),
            "Fastq": meta.get("Fastq", sample),
        }

        for level in LEVELS:
            # The truth table contains all curated stx targets, including two
            # targets that are interrupted/non-functional. For performance
            # calculations, only functional targets are allowed to contribute
            # to TP or FN. Non-functional targets are retained separately so
            # that detection of those targets is counted as FP.
            expected_all = set(truth.get(sample, {}).get(level, set()))
            nonfunctional = set(
                NONFUNCTIONAL_TARGETS_BY_SAMPLE.get(sample, {}).get(level, set())
            )
            expected = expected_all - nonfunctional
            detected = set(kma_hits.get(sample, {}).get(level, set()))

            tp = expected & detected
            fp = detected - expected
            fn = expected - detected
            detected_nonfunctional = detected & nonfunctional

            universe_size = len(universes[level])
            tn_count = universe_size - len(tp) - len(fp) - len(fn)

            row = {
                **base,
                "Level": level.title(),
                "Universe_Size": universe_size,
                "Expected_Count": len(expected),
                "Expected_All_Curated_Count": len(expected_all),
                "Nonfunctional_Count": len(nonfunctional),
                "Detected_Count": len(detected),
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "TN": tn_count,
                "Sensitivity": safe_divide(len(tp), len(tp) + len(fn)),
                "Specificity": safe_divide(tn_count, tn_count + len(fp)),
                "Precision": safe_divide(len(tp), len(tp) + len(fp)),
                "Accuracy": safe_divide(len(tp) + tn_count, universe_size),
                "F1": safe_divide(2 * len(tp), 2 * len(tp) + len(fp) + len(fn)),
                "Expected_Targets": "; ".join(sorted(expected)),
                "Expected_All_Curated_Targets": "; ".join(sorted(expected_all)),
                "Nonfunctional_Targets": "; ".join(sorted(nonfunctional)),
                "Detected_Targets": "; ".join(sorted(detected)),
                "Detected_Nonfunctional_Targets": "; ".join(sorted(detected_nonfunctional)),
                "TP_Targets": "; ".join(sorted(tp)),
                "FP_Targets": "; ".join(sorted(fp)),
                "FN_Targets": "; ".join(sorted(fn)),
            }
            summary_rows.append(row)

            for count_col in ["TP", "FP", "FN", "TN"]:
                totals[level][count_col] += row[count_col]

            # Do not write TN detail rows because that would create enormous sheets.
            # TN is still included in the per-genome Summary and aggregate Metrics tabs.
            for classification, targets in [("TP", tp), ("FP", fp), ("FN", fn)]:
                for target in sorted(targets):
                    detail_rows[level].append({
                        **base,
                        "Level": level.title(),
                        "Target": target,
                        "Expected": classification != "FP",
                        "Detected": classification != "FN",
                        "Nonfunctional_Target": target in nonfunctional,
                        "Classification": classification,
                        "Classification_Note": (
                            "Detected non-functional/interrupted target; counted as FP"
                            if classification == "FP" and target in nonfunctional
                            else ""
                        ),
                    })

    metrics_rows = []
    for level in LEVELS:
        universe_size = len(universes[level])
        total_comparisons = len(samples) * universe_size
        tp = totals[level]["TP"]
        fp = totals[level]["FP"]
        fn = totals[level]["FN"]
        tn = totals[level]["TN"]

        metrics_rows.append({
            "Parameter_Set": parameter_set,
            "Level": level.title(),
            "Genome_Count": len(samples),
            "Universe_Size_Per_Genome": universe_size,
            "Total_Comparisons": total_comparisons,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "TN": tn,
            "Sensitivity": safe_divide(tp, tp + fn),
            "Specificity": safe_divide(tn, tn + fp),
            "Precision": safe_divide(tp, tp + fp),
            "Accuracy": safe_divide(tp + tn, total_comparisons),
            "F1": safe_divide(2 * tp, 2 * tp + fp + fn),
        })

    return {
        "parameter_set": parameter_set,
        "summary": pd.DataFrame(summary_rows),
        "metrics": pd.DataFrame(metrics_rows),
        "details": {level: pd.DataFrame(rows) for level, rows in detail_rows.items()},
    }


def style_workbook(xlsx_path: Path):
    """Apply simple readable formatting: bold headers, no alternating row colours."""
    wb = load_workbook(xlsx_path)
    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(bold=True, color="FFFFFF")
    thin = Side(style="thin", color="D9E2F3")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            cell.border = border

        for row in ws.iter_rows(min_row=2):
            for cell in row:
                cell.border = border
                if isinstance(cell.value, str) and len(cell.value) > 80:
                    cell.alignment = Alignment(wrap_text=True, vertical="top")

        # Widths are capped so long target-list columns do not become enormous.
        for col in ws.columns:
            col_letter = col[0].column_letter
            header = str(col[0].value or "")
            if header.endswith("Targets"):
                ws.column_dimensions[col_letter].width = 45
            elif header in {"Target"}:
                ws.column_dimensions[col_letter].width = 30
            elif header in {"Parameter_Set", "Genome", "Fasta", "Fastq"}:
                ws.column_dimensions[col_letter].width = 24
            else:
                ws.column_dimensions[col_letter].width = 16

        for row in ws.iter_rows(min_row=2):
            for cell in row:
                if ws.cell(row=1, column=cell.column).value in {"Sensitivity", "Specificity", "Precision", "Accuracy", "F1"}:
                    cell.number_format = "0.00%"

    wb.save(xlsx_path)


def write_result_workbook(result: dict, universes: Dict[str, Set[str]], output_path: Path):
    """Write one output workbook for one KMA parameter file."""
    readme = pd.DataFrame([
        {"Item": "Analysis", "Value": "Independent KMA three-level secondary comparison"},
        {"Item": "Operon universe", "Value": len(universes["operon"])},
        {"Item": "Variant universe", "Value": len(universes["variant"])},
        {"Item": "Subtype universe", "Value": len(universes["subtype"])},
        {"Item": "TN formula", "Value": "Universe size - TP - FP - FN"},
        {"Item": "Functional-target rule", "Value": "Only functional expected targets can be TP or FN. Detection of specified interrupted/non-functional targets is FP."},
        {"Item": "Special case 1", "Value": "2025-SEQ-1796: StxOp1089_Stx2a_1 / Stx2a_1 / Stx2a is interrupted by a transposon and is counted as FP if detected."},
        {"Item": "Special case 2", "Value": "SRR18191635: StxOp3298_Stx2b_23 / Stx2b_23 / Stx2b has a frameshift and is counted as FP if detected."},
        {"Item": "Details tabs", "Value": "TP, FP, and FN rows only; TN rows are implicit."},
    ])

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        result["metrics"].to_excel(writer, sheet_name="Metrics", index=False)
        result["summary"].to_excel(writer, sheet_name="Summary_by_Genome", index=False)
        result["details"]["operon"].to_excel(writer, sheet_name="Operon_Details", index=False)
        result["details"]["variant"].to_excel(writer, sheet_name="Variant_Details", index=False)
        result["details"]["subtype"].to_excel(writer, sheet_name="Subtype_Details", index=False)

    style_workbook(output_path)


def write_combined_workbook(results: Iterable[dict], universes: Dict[str, Set[str]], output_path: Path):
    """Write one combined workbook across all KMA parameter files."""
    results = list(results)
    readme = pd.DataFrame([
        {"Item": "Analysis", "Value": "Combined independent KMA secondary comparison"},
        {"Item": "Operon universe", "Value": len(universes["operon"])},
        {"Item": "Variant universe", "Value": len(universes["variant"])},
        {"Item": "Subtype universe", "Value": len(universes["subtype"])},
        {"Item": "Functional-target rule", "Value": "Only functional expected targets can be TP or FN. Detection of specified interrupted/non-functional targets is FP."},
        {"Item": "Special case 1", "Value": "2025-SEQ-1796: StxOp1089_Stx2a_1 / Stx2a_1 / Stx2a is interrupted by a transposon and is counted as FP if detected."},
        {"Item": "Special case 2", "Value": "SRR18191635: StxOp3298_Stx2b_23 / Stx2b_23 / Stx2b has a frameshift and is counted as FP if detected."},
    ])

    combined_metrics = pd.concat([r["metrics"] for r in results], ignore_index=True)
    combined_summary = pd.concat([r["summary"] for r in results], ignore_index=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        combined_metrics.to_excel(writer, sheet_name="Combined_Metrics", index=False)
        combined_summary.to_excel(writer, sheet_name="Combined_Summary", index=False)

    style_workbook(output_path)


def main():
    parser = argparse.ArgumentParser(description="Independent secondary KMA comparison at operon, variant, and subtype levels.")
    parser.add_argument("--truth-xlsx", required=True, type=Path, help="Validation ground-truth Excel file.")
    parser.add_argument("--truth-sheet", default="Final_Dataset", help="Truth sheet name. Default: Final_Dataset")
    parser.add_argument("--db-fasta", required=True, type=Path, help="Stx operon database FASTA used to build the KMA index.")
    parser.add_argument("--kma-files", required=True, nargs="+", type=Path, help="One or more combined KMA TSV result files.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory where output Excel workbooks will be written.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    universes = load_database_universes(args.db_fasta)
    print("Database universes:")
    for level in LEVELS:
        print(f"  {level}: {len(universes[level])}")

    truth = load_truth(args.truth_xlsx, args.truth_sheet)
    print(f"Loaded truth for {len(truth)} genomes.")

    results = []
    for kma_file in args.kma_files:
        print(f"Comparing {kma_file.name}...")
        result = compare_one_kma_file(kma_file, truth, universes)
        results.append(result)

        output_path = args.output_dir / f"{result['parameter_set']}_secondary_three_level_comparison.xlsx"
        write_result_workbook(result, universes, output_path)
        print(f"  wrote {output_path}")

    combined_output = args.output_dir / "kma_secondary_three_level_combined_summary.xlsx"
    write_combined_workbook(results, universes, combined_output)
    print(f"Wrote combined summary: {combined_output}")


if __name__ == "__main__":
    main()
