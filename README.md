# STEC Genomic Characterisation

Analysis code and derived data for the MSc thesis *Genomic Characterisation and Risk
Assessment of Shiga Toxin-Producing Escherichia coli* (N. Shubair, Carleton University,
in collaboration with the Canadian Food Inspection Agency).


## What this repository contains

Three analyses, one per thesis chapter, each with the scripts that produced the results
and the derived data files those scripts read. Raw sequencing reads and genome assemblies
are **not** stored here — they are in NCBI under BioProject
[PRJNA454819](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA454819), and per-isolate
accessions are listed in Supplementary Table S3.1.

```
chapter2_stx_detection_validation/    validation of four stx detection workflows
├── scripts/                          parsing, subtype comparison, metrics, figures
├── data/                             ground truth and per-tool comparison outputs
└── figures/                          Figures 2.1-2.4

chapter3_cfia_stec_characterisation/  497 CFIA STEC genomes
├── scripts/                          R scripts for Figures 3.1-3.3 and the summary tables
├── data/                             virulence presence/absence matrices and metadata
└── figures/                          Figures 3.1-3.3 and two exploratory heatmaps

chapter4_stx2i_lamb/                  66 stx2i-positive genomes
├── scripts/
├── data/                             curated metadata, virulence tables, chewBBACA output
└── figures/                          Figures 4.1-4.2

supplementary_tables/                 the electronic supplementary tables (xlsx)
```

## Chapter 2 — validation of stx detection and subtyping

Four workflows were compared against a curated ground truth built from 200 STEC genomes
spanning all 19 recognised Stx subtypes: **KMA** and **Sipprverse** on raw reads, and
**GeneSeekr** and **StxTyper** on SKESA assemblies. All four searched the STxOP_nt v1.0.0
reference database ([STxDB](https://github.com/OLC-Bioinformatics/STxDB)).

| Script | What it does |
|---|---|
| **Parse** | |
| `parse_geneseekr_best_hit_per_region.py` | Groups GeneSeekr BLASTN hits into genomic regions and keeps the best hit per region (bit score, then e-value, percent match, alignment length). |
| `extract_sipprverse_hits.py` | Reduces the wide Sipprverse CSV to Genome / Hit / Subtype / Identity, discarding sequence rows. |
| `combine_stxtyper_tsvs.py` | Concatenates the 200 per-genome StxTyper TSVs into one table. |
| **Score** | |
| `kma_subtype_comparison_pipeline.py` | Scores KMA across all 12 identity × ConClave parameter sets (Table 2.4). |
| `sipprverse_subtype_comparison_pipeline.py` | Scores Sipprverse at 98% identity. |
| `geneseekr_subtype_comparison_cutoff.py` | Scores GeneSeekr at ≥ 90% percent match — the analysis reported in the thesis. |
| `geneseekr_subtype_comparison.py` | The same comparison with no cutoff applied. |
| `stxtyper_subtype_comparison.py` | Scores StxTyper on COMPLETE operons only. |
| `kma_secondary_independent_three_level_comparison_pipeline_19_subtypes.py` | Secondary three-level validation: nucleotide operon, variant and subtype scored independently. |
| **Plot** | |
| `comparison_metrics_with_confidence_intervals_plot_allmetrics.R` | Figures 2.1 and 2.2; Clopper–Pearson 95% confidence intervals. |
| `radar_chart_F1scores.R` | Figure 2.3 radar plots (fmsb). |
| `secondary_KMA70_conclave1_allmetrics_CI.R` | Figure 2.4. |

Non-functional targets are handled explicitly: a subtype encoded by a disrupted operon
cannot be a true positive, counts as a false positive when detected, and does not become a
false negative when absent.

Two of these scripts import `artifact_tool`, which is not a public package and will not
install — see the chapter README.

## Chapter 3 — genomic characterisation of 497 CFIA STEC genomes

| Script | Output |
|---|---|
| `scatterpie_plot_subtypes_by_source.R` | Figure 3.1, Stx subtype distribution by isolation source |
| `bubble_plot_variants_by_source.R` | Figure 3.2, Stx variant distribution by isolation source |
| `chapter3_virulence_heatmap_497_corrected.R` | virulence gene presence/absence heatmap across all 497 genomes |
| `chapter3_virulence_heatmap_497_prevalence.R` | prevalence-scaled variant of the same heatmap |
| `chapter3_fig3.3_virulence_by_O157.R` | Figure 3.3, virulence gene prevalence in O157:H7 versus non-O157 |
| `chapter3_fig3.4_risk_levels.R` | FAO/WHO risk level assigned to each genome, and the counts underlying Table 3.6 |

`data/chapter3_497_virulence_presence_absence.tsv` is the binary genome × gene matrix that
the heatmap scripts read; `chapter3_497_virulence_gene_categories.tsv` assigns each gene to
a functional category, and `virulencefinder_497_identifier_mapping.tsv` reconciles the
SRR, GCA, CFIA and SEQ identifiers used across the source files.

## Chapter 4 — stx2i-positive STEC from lamb

Seven CFIA isolates from imported raw lamb, hybrid-assembled from Illumina and Oxford
Nanopore reads, analysed alongside 59 curated public genomes.

| Script | Output |
|---|---|
| `spanning_tree.R` | Figure 4.1 — cgMLST minimum-spanning tree annotated with serotype, isolation source, country of origin and MLST (ggtree) |
| `virulence_heatmap_with_source_annotation.R` | Figure 4.2 — virulence gene presence/absence across the 66 genomes, clustered by Ward D2 (pheatmap) |

`data/Stx2i_NCBI_CFIA_genomes_working_260730.xlsx` (sheet `NCBI_CFIA_dataset`) is the
metadata source for both scripts; `combined_virulence_results.tsv` is the raw
VirulenceFinder output, and the four derived tables alongside it are the heatmap script's
own outputs. `data/chewbbaca/` holds the chewBBACA allele calls, the 95% core schema and
the GrapeTree Newick file that `spanning_tree.R` reads.

## Software

| Tool | Version | Tool | Version |
|---|---|---|---|
| COWBAT | 0.5.0.23 | KMA | 1.6.8 |
| FastQC | 0.11.8 (COWBAT) / 0.11.6 (Ch. 4) | Sipprverse | 0.2.46 |
| BBTools (bbduk, Tadpole) | 38.22 (COWBAT) / 37.78 (Ch. 4) | GeneSeekr | 0.5.0 |
| ConFindr | 0.4.7 | StxTyper | 1.0.45 |
| SKESA | 2.3.0 | ECTyper | 2.0.0 |
| Pilon | 1.22 | VirulenceFinder | 3.2.1 (DB 2022-12-02) |
| Prodigal | 2.6.3 | ResFinder | 4.7.2 |
| Flye | 2.9.6-b1802 | MOB-recon | 1.4.0 |
| Unicycler | 0.4.8 | chewBBACA | 3.5.4 |
| Porechop / Filtlong | 0.2.4 / 0.3.1 | GrapeTree | 1.5.0 |
| Bowtie2 | 2.3.4-2.3.5 | SNVPhyl | 2.3.4 |
| STxOP_nt database | 1.0.0 | R | 4.5.3 |

R packages: tidyverse, readxl, readr, dplyr, tidyr, stringr, reshape2, ggplot2, pheatmap,
scatterpie, ggtree, RColorBrewer, scales, fmsb.

## Reproducing an analysis

```bash
git clone https://github.com/OLC-Bioinformatics/STEC-Genomic-Characterisation.git
cd STEC-Genomic-Characterisation/chapter3_cfia_stec_characterisation
Rscript scripts/chapter3_fig3.3_virulence_by_O157.R
```

Run each script from its chapter directory. Raw reads must be downloaded from
PRJNA454819 for any step that starts from sequence data; the scripts here begin from the
derived tables in `data/`.

## Citation

Shubair, N. (2026). *Genomic Characterisation and Risk Assessment of Shiga
Toxin-Producing Escherichia coli* [MSc thesis]. Carleton University.

## Related repositories

- [STxDB](https://github.com/OLC-Bioinformatics/STxDB) — the STxOP_nt reference database
- [GeneSeekr](https://github.com/OLC-Bioinformatics/GeneSeekr)
- [sipprverse](https://github.com/OLC-Bioinformatics/sipprverse)
- [COWBAT](https://github.com/OLC-Bioinformatics/COWBAT)
- [StxTyper](https://github.com/ncbi/stxtyper) (NCBI)
- [chewBBACA](https://github.com/B-UMMI/chewBBACA) (B-UMMI) — cgMLST/wgMLST schema creation and allele calling
- [GrapeTree](https://github.com/achtman-lab/GrapeTree) (Achtman lab) — MSTreeV2

## Contact

Noor Shubair — Carleton University / Canadian Food Inspection Agency
