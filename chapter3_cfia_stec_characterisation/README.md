# Chapter 3 — Genomic characterisation of STEC from Canadian food surveillance

497 STEC genomes recovered through CFIA food surveillance, characterised by serotype,
sequence type, Stx subtype and variant, virulence gene content, virulence genotype and
predicted FAO/WHO risk level.

## Scripts

| Script | Output |
|---|---|
| `scatterpie_plot_subtypes_by_source.R` | `figures/figure3.1_stx_subtype_by_source.jpg` — Figure 3.1, Stx subtype distribution by isolation source |
| `bubble_plot_variants_by_source.R` | `figures/figure3.2_stx_variant_by_source.jpg` — Figure 3.2, Stx variant distribution by isolation source |
| `chapter3_fig3.3_virulence_by_O157.R` | `figures/figure3.3_virulence_prevalence_O157_vs_nonO157.png` — Figure 3.3, virulence gene prevalence in O157:H7 versus non-O157 |
| `chapter3_virulence_heatmap_497_corrected.R` | `figures/virulence_gene_profiles_497_genomes.jpg` — presence/absence heatmap, all 497 genomes |
| `chapter3_virulence_heatmap_497_prevalence.R` | `figures/virulence_gene_prevalence_by_source_497.jpg` — gene prevalence within each isolation source |
| `chapter3_fig3.4_risk_levels.R` | risk level assigned to each genome, and the counts underlying Table 3.6 |

## Figures

`figures/` holds the rendered output of the scripts above.

Figures 3.1, 3.2 and 3.3 are the three figures that appear in Chapter 3 of the thesis. The
two 497-genome heatmaps are exploratory views of the same presence/absence matrix and are
**not** thesis figures: `virulence_gene_profiles_497_genomes.jpg` shows every genome as a
row, clustered by Ward D2, and `virulence_gene_prevalence_by_source_497.jpg` collapses the
genomes into their seven isolation sources and shades each cell by prevalence within that
source. Both are included because they informed how the gene set was narrowed for
Figure 3.3.

`chapter3_fig3.4_risk_levels.R` draws a stacked bar chart of the FAO/WHO risk levels that
was not used in the thesis. It is kept because the same script assigns a risk level to each
of the 497 genomes and prints the O157:H7 and non-O157 counts that populate **Table 3.6**.

## Data

| File | Contents |
|---|---|
| `CFIA_STEC.xlsx` | master metadata workbook; sheet `prjna454819` is the 497-genome list |
| `chapter3_497_virulence_presence_absence.tsv` | binary genome × virulence gene matrix |
| `chapter3_497_virulence_gene_categories.tsv` | functional category for each gene |
| `chapter3_497_virulence_genome_annotations.tsv` | per-genome annotations used for heatmap row labels |
| `chapter3_497_virulence_prevalence_by_source.tsv` | gene prevalence within each isolation source |
| `virulencefinder_497_identifier_mapping.tsv` | crosswalk between SRR, GCA, CFIA and SEQ identifiers |
| `virulencefinder_results_497_final.xlsx` | raw VirulenceFinder output |
| `Chapter3_analysis_tables.xlsx` | summary tables underlying Tables 3.1-3.6 |

## Risk classification

Risk levels follow FAO and WHO (2018) as operationalised by Zhang et al. (2025):
level 1, stx2a with eae or aggR; level 2, stx2d; level 3, stx2c with eae; level 4, stx1a
with eae; level 5, any other stx subtype or combination. Where a genome meets more than one
level, the highest is assigned. Genomes with no stx detected are recorded as not
classifiable.
