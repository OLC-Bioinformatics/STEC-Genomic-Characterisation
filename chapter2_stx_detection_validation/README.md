# Chapter 2 — Validation of WGS approaches for stx detection and subtyping

Four workflows scored against a curated ground truth of 200 STEC genomes covering all 19
recognised Stx subtypes (293 operons, 136 unique operon sequences, 90 variants).

| Workflow | Input | Detection | Threshold used |
|---|---|---|---|
| KMA v1.6.8 | raw reads | k-mer alignment to STxOP_nt | 70% template identity, ConClave mode 1 |
| Sipprverse v0.2.46 | raw reads | read baiting + Bowtie2 | 98% identity |
| GeneSeekr v0.5.0 | SKESA assemblies | BLASTN against STxOP_nt | ≥ 90% percent match |
| StxTyper v1.0.45 | SKESA assemblies | translated BLAST, complete operons | defaults |

## The analysis chain

Each tool follows the same three steps: parse the raw output into a comparable form, score
it against the shared ground truth, then pool the results for plotting.

### 1. Parse

| Script | Reads | Writes |
|---|---|---|
| `parse_geneseekr_best_hit_per_region.py` | a directory of `*_blastn_geneseekr.tsv` | `data/geneseekr/parsed_geneseekr_outputs_best_hit_per_region.csv` |
| `extract_sipprverse_hits.py` | the wide Sipprverse CSV | `data/sipprverse/sipprverse_0.98_parsed_file.csv` |
| `combine_stxtyper_tsvs.py` | `data/stxtyper/stxtyper_results/` (200 per-genome TSVs) | `data/stxtyper/stxtyper_combined_results.tsv` |

KMA needs no parse step; its TSV output is read directly.

`parse_geneseekr_best_hit_per_region.py` groups BLAST hits into genomic regions by
overlapping query coordinates and keeps the best hit in each region, ranked by bit score,
then e-value, percent match and alignment length. Grouping by region rather than by genome
is what allows a genome carrying more than one stx locus to contribute more than one call.

### 2. Score against the ground truth

All four comparison scripts read `data/ground_truth/StxDB_groundtruth_subtype_level.xlsx`
(sheet `Operon_vs_Subtype`) and write a workbook with `Summary`, `Subtype_Details`,
`Ground_Truth` and `Metrics` sheets.

| Script | Output |
|---|---|
| `kma_subtype_comparison_pipeline.py` | `data/kma/comparison_outputs/` — one workbook per parameter set, plus `KMA_parameter_metrics_comparison.xlsx` (Table 2.4) |
| `sipprverse_subtype_comparison_pipeline.py` | `data/sipprverse/sipprverse_subtype_comparison_results_0.98.xlsx` |
| `geneseekr_subtype_comparison_cutoff.py` | `data/geneseekr/GeneSeekr_Subtype_Comparison_90.xlsx` — the ≥ 90% analysis reported in the thesis |
| `geneseekr_subtype_comparison.py` | the same comparison with no percent-match cutoff |
| `stxtyper_subtype_comparison.py` | `data/stxtyper/StxTyper_Subtype_Comparison.xlsx` — COMPLETE operons only |
| `kma_secondary_independent_three_level_comparison_pipeline_19_subtypes.py` | `data/kma/kmaresults_70_conclave1_secondary_three_level_comparison_19.xlsx` — operon, variant and subtype levels scored independently |

### 3. Pool and plot

`data/Subtype_Level_Tool_Metrics_Comparison.xlsx` and `data/all_tools_subtype_data_forR.xlsx`
collect the per-tool metrics. The three R scripts read the latter and produce the figures in
`figures/`: Figures 2.1 and 2.2 (`comparison_metrics_with_confidence_intervals_plot_allmetrics.R`),
Figure 2.3 (`radar_chart_F1scores.R`) and Figure 2.4
(`secondary_KMA70_conclave1_allmetrics_CI.R`). Confidence intervals are exact
Clopper–Pearson binomial intervals from `binom.test`.

## Scoring convention

A subtype is counted once per genome irrespective of copy number, and true negatives are
the remainder of the 19 recognised subtypes. A subtype encoded by a disrupted or
non-functional operon cannot be a true positive, is a false positive when detected, and is
not a false negative when absent. Two such targets are hard-coded in every comparison
script: `2025-SEQ-1796` (transposon insertion in stx2a subunit A) and `SRR18191635`
(frameshift in stx2b subunit A).

StxTyper is scored on rows whose operon status is exactly `COMPLETE`; `PARTIAL`,
`PARTIAL_CONTIG_END` and `COMPLETE_NOVEL` are ignored.

## Running the scripts

```bash
pip install pandas openpyxl

python scripts/kma_subtype_comparison_pipeline.py \
  --groundtruth data/ground_truth/StxDB_groundtruth_subtype_level.xlsx \
  --kma data/kma/raw_outputs/kmaresults_70_conclave1.tsv \
  --output KMA_70_conclave1_noTN.xlsx
```

Each script prints its own usage with `--help`.

> **Two scripts will not run as supplied.** `geneseekr_subtype_comparison.py` and
> `stxtyper_subtype_comparison.py` begin with `from artifact_tool import Blob,
> SpreadsheetFile, Workbook`. `artifact_tool` is not a public package and cannot be
> installed, so these two need their workbook-writing sections ported to `openpyxl` before
> they will execute. The other six Python scripts use only `pandas` and `openpyxl` and run
> as they are. `geneseekr_subtype_comparison_cutoff.py` — the version whose results the
> thesis reports — is among the working ones.

## Data

| Path | Contents |
|---|---|
| `data/ground_truth/` | subtype-level ground truth for the 200 benchmark genomes |
| `data/geneseekr/` | best-hit-per-region calls and the ≥ 90% comparison workbook |
| `data/kma/raw_outputs/` | 12 KMA TSVs — six identity thresholds × two ConClave modes |
| `data/kma/comparison_outputs/` | 12 comparison workbooks and the pooled parameter metrics |
| `data/sipprverse/` | parsed hits and the 98% comparison workbook |
| `data/stxtyper/` | 200 per-genome TSVs, the combined table (293 rows, 199 genomes) and the comparison workbook |
| `data/Subtype_Level_Tool_Metrics_Comparison.xlsx` | pooled metrics for all four workflows |
| `data/all_tools_subtype_data_forR.xlsx` | tidy input for the plotting scripts |

Two notes on the supplied files. In `data/stxtyper/stxtyper_results/`, 162 of the 200
filenames carry an `_FP` suffix; this is a local file-naming convention and has nothing to
do with false positives — genome identity comes from the `#name` field inside each file.
And `combine_stxtyper_tsvs.py` writes `stxtyper_combined.tsv`, while the file stored here
is named `stxtyper_combined_results.tsv`.

## Not included

The raw Sipprverse output (`sipprverse_output_0.98.csv`, 185 MB) exceeds GitHub's 100 MB
file limit. It is mostly sequence rows, which `extract_sipprverse_hits.py` discards; the
parsed 6.3 MB CSV it produces is stored here instead and is the direct input to the
comparison script. The raw file can be regenerated by rerunning Sipprverse at 98% identity.
