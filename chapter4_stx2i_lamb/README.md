# Chapter 4 — Sequence analysis of Stx2i-producing E. coli from lamb

Seven CFIA isolates recovered from imported raw lamb products, hybrid-assembled from
Illumina and Oxford Nanopore reads, analysed together with 59 curated publicly available
stx2i-positive genomes (66 genomes in total).

## Scripts

| Script | Output |
|---|---|
| `spanning_tree.R` | Figure 4.1 — cgMLST minimum-spanning tree. Reads the GrapeTree Newick file, attaches isolate metadata with `ggtree`, and draws tip points coloured by serotype and shaped by isolation source, with country of origin and MLST added as `gheatmap` panels. Also produces unrooted (daylight) and circular layouts. |
| `virulence_heatmap_with_source_annotation.R` | Figure 4.2 — virulence gene presence/absence heatmap. Builds the binary genome × gene matrix from the combined VirulenceFinder output, assigns each gene to a functional category, annotates rows by isolation source, clusters genomes by Ward D2, and writes the four derived tables listed below. |

## Figures

`figures/` holds the rendered output of both scripts, at the resolution used in the thesis.

| File | Figure |
|---|---|
| `figure4.1_core_genome_phylogeny_stx2i.jpg` | Figure 4.1 — core genome phylogeny of the 66 Stx2i-positive genomes, tips coloured by serotype and shaped by isolation source, with a country-of-origin panel |
| `figure4.2_virulence_gene_profiles_stx2i.jpg` | Figure 4.2 — presence/absence of 45 virulence genes across the 66 genomes |

## Data

| File | Contents |
|---|---|
| `Stx2i_NCBI_CFIA_genomes_working_260730.xlsx` | curation record for the assemblies retrieved from NCBI Pathogen Detection, with the MicroBIGG-E call, the KMA result and the retain/exclude outcome (Supplementary Table S4.2). Sheet `NCBI_CFIA_dataset` is the metadata source for both scripts. |
| `Stx2i_genomes_working_260730_forR_cgmlstfigures.xlsx` | metadata used to annotate the cgMLST tree (Figure 4.1) |
| `combined_virulence_results.tsv` | raw combined VirulenceFinder output for the 66 genomes (1,814 hits) |
| `virulence_presence_absence.tsv` | binary matrix, 66 genomes × 45 genes — written by the heatmap script |
| `virulence_gene_categories.tsv` | functional category assigned to each of the 45 genes |
| `virulence_gene_prevalence.tsv` | genome count and prevalence for each gene |
| `virulence_genome_source_annotations.tsv` | per-genome source, serotype, host and isolation type |

The four derived tables are outputs of `virulence_heatmap_with_source_annotation.R` and are
included so the figure can be inspected without re-running the analysis.

> **Not included:** `cgmlst95_MSTreeV2.nwk`, the GrapeTree Newick file that
> `spanning_tree.R` reads. It is produced by the chewBBACA → GrapeTree steps described
> below and must be added before that script will run.

## Analysis outline

1. Hybrid assembly — Flye v2.9.6 long-read assembly passed to Unicycler v0.4.8 with the
   processed Illumina reads, polished with Pilon v1.22.
2. Public genome retrieval — NCBI Pathogen Detection MicroBIGG-E and Isolates browsers,
   filtered on the element symbol `stx2i operon`; 99 assemblies retrieved, 59 retained.
3. Verification — blastn against STxOP_nt v1.0.0 and StxTyper v1.0.45 on the assemblies,
   with KMA v1.6.8 on the corresponding raw reads as the independent check.
4. cgMLST — chewBBACA v3.5.4, core schema at 95% locus presence (3,249 of 7,219 loci),
   minimum-spanning tree with GrapeTree v1.5.0 (MSTreeV2), annotated in R with ggtree.
5. Virulence gene content — VirulenceFinder v3.2.1 (90% identity, 60% coverage); 45 genes
   across seven functional categories.
6. Clonality — SNVPhyl v2.3.4; isolates within 20 hqSNVs of each other from the same
   sample were treated as clonal and represented by one strain.

## R packages

`tidyverse` · `readxl` · `pheatmap` · `RColorBrewer` · `ggplot2` · `ggtree` · `ggtreeExtra` ·
`treeio` · `ape` · `ggnewscale` · `janitor` · `paletteer`

`ggtree`, `ggtreeExtra` and `treeio` are Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))
```
