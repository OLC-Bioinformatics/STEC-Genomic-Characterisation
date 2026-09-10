# ==============================================================================
# Chapter 3 VirulenceFinder presence/absence heatmap
# 497 CFIA STEC genomes from the prjna454819 sheet
# ==============================================================================
#
# This script:
#   1. Reads the UPDATED VirulenceFinder output.
#   2. Reads the 497-genome list from the "prjna454819" sheet.
#   3. Uses the "alldata" sheet as an identifier crosswalk so VirulenceFinder
#      SEQIDs recorded as SRR, GCA/REFSEQ, CFIA, SEQ, or DAR-style IDs are all
#      mapped back to the REFSEQ identifier used in prjna454819.
#   4. Confirms that all 497 genomes are represented.
#   5. Creates a binary genome x virulence-gene matrix.
#   6. Adds broad isolation-source (Source1) annotation.
#   7. Orders genes by functional category and prevalence.
#   8. Saves the heatmap and supporting TSV files.
#
# Required packages:
# install.packages(c("tidyverse", "pheatmap", "readxl", "RColorBrewer"))
# ==============================================================================

library(tidyverse)
library(pheatmap)
library(readxl)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------------------------

# CHANGE THESE TWO PATHS to where the files are located on your computer.
virulence_file <- "/Users/noorshubair/Desktop/chapter3/virulencefinder_results_497_final.xlsx"
metadata_file  <- "/Users/noorshubair/Desktop/R_code_files/chapter3/CFIA_STEC.xlsx"

# Folder where output files will be written.
output_dir <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(virulence_file)) {
  stop("VirulenceFinder workbook not found: ", virulence_file)
}

if (!file.exists(metadata_file)) {
  stop("CFIA_STEC workbook not found: ", metadata_file)
}

# ------------------------------------------------------------------------------
# 2. READ VIRULENCEFINDER RESULTS
# ------------------------------------------------------------------------------

# Uses the first worksheet in the VirulenceFinder workbook.
virulence_results <- read_excel(virulence_file, sheet = 1)
colnames(virulence_results) <- str_trim(colnames(virulence_results))

required_vf_columns <- c("SEQID", "Virulence factor")
missing_vf_columns <- setdiff(required_vf_columns, colnames(virulence_results))

if (length(missing_vf_columns) > 0) {
  stop(
    "The following required VirulenceFinder column(s) are missing: ",
    paste(missing_vf_columns, collapse = ", ")
  )
}

virulence_results <- virulence_results %>%
  mutate(
    SEQID = str_trim(as.character(SEQID)),
    `Virulence factor` = str_trim(as.character(`Virulence factor`))
  )

message(
  "Unique SEQIDs in complete VirulenceFinder workbook: ",
  n_distinct(virulence_results$SEQID)
)

# ------------------------------------------------------------------------------
# 3. READ THE FINAL 497-GENOME DATASET
# ------------------------------------------------------------------------------

final_497 <- read_excel(metadata_file, sheet = "prjna454819")
colnames(final_497) <- str_trim(colnames(final_497))

required_final_columns <- c(
  "REFSEQ",
  "SRA accession",
  "Strain",
  "Source1",
  "Source2",
  "Serotype_HC list"
)

missing_final_columns <- setdiff(required_final_columns, colnames(final_497))

if (length(missing_final_columns) > 0) {
  stop(
    "The following required prjna454819 column(s) are missing: ",
    paste(missing_final_columns, collapse = ", ")
  )
}

final_497 <- final_497 %>%
  transmute(
    REFSEQ = str_trim(as.character(REFSEQ)),
    SRA_accession = str_trim(as.character(`SRA accession`)),
    Strain = str_trim(as.character(Strain)),
    Source1 = str_trim(as.character(Source1)),
    Source2 = str_trim(as.character(Source2)),
    Serotype = str_trim(as.character(`Serotype_HC list`)),
    STX = if ("STX" %in% colnames(final_497)) {
      str_trim(as.character(STX))
    } else {
      NA_character_
    }
  ) %>%
  filter(!is.na(REFSEQ), REFSEQ != "") %>%
  distinct(REFSEQ, .keep_all = TRUE)

if (nrow(final_497) != 497) {
  stop(
    "Expected exactly 497 unique REFSEQ genomes in prjna454819, but found ",
    nrow(final_497), "."
  )
}

message("Final Chapter 3 dataset contains 497 unique genomes.")

# ------------------------------------------------------------------------------
# 4. MAP VIRULENCEFINDER SEQIDS TO THE FINAL 497 REFSEQ IDS
# ------------------------------------------------------------------------------

# If the VirulenceFinder workbook already contains a Matched_REFSEQ column,
# use it directly. Otherwise, construct a crosswalk from prjna454819 + alldata.
if ("Matched_REFSEQ" %in% colnames(virulence_results)) {

  message("Matched_REFSEQ column detected in VirulenceFinder workbook; using it directly.")

  virulence_497 <- virulence_results %>%
    mutate(Matched_REFSEQ = str_trim(as.character(Matched_REFSEQ))) %>%
    filter(Matched_REFSEQ %in% final_497$REFSEQ)

  match_summary <- virulence_497 %>%
    distinct(Matched_REFSEQ, SEQID) %>%
    mutate(Identifier_type = "pre-matched") %>%
    left_join(final_497, by = c("Matched_REFSEQ" = "REFSEQ")) %>%
    arrange(Matched_REFSEQ)

} else {

  # alldata contains the alternate identifiers that were used as SEQID in the
  # original VirulenceFinder run. Restrict the crosswalk to identifier fields
  # that can genuinely identify a genome. This avoids treating shared metadata
  # such as BioProject or source values as genome IDs.
  alldata <- read_excel(metadata_file, sheet = "alldata")
  colnames(alldata) <- str_trim(colnames(alldata))

  if (!"REFSEQ" %in% colnames(alldata)) {
    stop("The alldata sheet does not contain a REFSEQ column.")
  }

  # Only these fields are used as possible genome identifiers.
  possible_id_columns <- c(
    "REFSEQ",
    "SEQID",
    "Strain",
    "SRA accession",
    "SRA accession.1",
    "Assembly",
    "Isolate identifiers"
  )

  id_columns <- intersect(possible_id_columns, colnames(alldata))

  # First restrict alldata to the 497 final REFSEQ records, then pivot only the
  # identifier columns. This preserves mappings such as:
  #   GCA_012149165.1 -> SRR7536940
  #   CFIAFB20180581  -> 2019-DAR-0014
  crosswalk_alldata <- alldata %>%
    mutate(REFSEQ = str_trim(as.character(REFSEQ))) %>%
    filter(REFSEQ %in% final_497$REFSEQ) %>%
    select(REFSEQ, all_of(setdiff(id_columns, "REFSEQ"))) %>%
    pivot_longer(
      cols = -REFSEQ,
      names_to = "Identifier_type",
      values_to = "Identifier"
    ) %>%
    mutate(Identifier = str_trim(as.character(Identifier))) %>%
    filter(
      !is.na(Identifier),
      Identifier != "",
      Identifier != "NA"
    )

  # Also include identifiers directly present in prjna454819.
  crosswalk_final <- bind_rows(
    final_497 %>%
      transmute(
        REFSEQ,
        Identifier_type = "REFSEQ",
        Identifier = REFSEQ
      ),
    final_497 %>%
      transmute(
        REFSEQ,
        Identifier_type = "SRA accession",
        Identifier = SRA_accession
      ),
    final_497 %>%
      transmute(
        REFSEQ,
        Identifier_type = "Strain",
        Identifier = Strain
      )
  ) %>%
    mutate(Identifier = str_trim(as.character(Identifier))) %>%
    filter(!is.na(Identifier), Identifier != "", Identifier != "NA")

  identifier_crosswalk <- bind_rows(crosswalk_alldata, crosswalk_final) %>%
    distinct(REFSEQ, Identifier, .keep_all = TRUE)

  # Some identifiers (for example shared sample-level identifiers) can refer to
  # more than one genome. They are excluded because only one-to-one identifiers
  # are safe for matching VirulenceFinder SEQIDs.
  ambiguous_identifiers <- identifier_crosswalk %>%
    distinct(REFSEQ, Identifier) %>%
    count(Identifier, name = "n_REFSEQ") %>%
    filter(n_REFSEQ > 1)

  if (nrow(ambiguous_identifiers) > 0) {
    message(
      nrow(ambiguous_identifiers),
      " shared identifier(s) were excluded from matching because they map to more than one REFSEQ genome."
    )
  }

  identifier_crosswalk_unique <- identifier_crosswalk %>%
    anti_join(ambiguous_identifiers, by = "Identifier") %>%
    distinct(Identifier, .keep_all = TRUE)

  virulence_497 <- virulence_results %>%
    inner_join(
      identifier_crosswalk_unique %>%
        select(Identifier, Matched_REFSEQ = REFSEQ, Identifier_type),
      by = c("SEQID" = "Identifier")
    )

  match_summary <- virulence_497 %>%
    distinct(Matched_REFSEQ, SEQID, Identifier_type) %>%
    left_join(final_497, by = c("Matched_REFSEQ" = "REFSEQ")) %>%
    arrange(Matched_REFSEQ)
}

# ------------------------------------------------------------------------------
# 5. VERIFY THAT ALL 497 GENOMES WERE MATCHED
# ------------------------------------------------------------------------------

matched_genomes <- sort(unique(virulence_497$Matched_REFSEQ))
missing_genomes <- setdiff(final_497$REFSEQ, matched_genomes)

message("Genomes matched to VirulenceFinder results: ", length(matched_genomes), " / 497")

if (length(missing_genomes) > 0) {
  write_tsv(
    tibble(REFSEQ = missing_genomes),
    file.path(output_dir, "UNMATCHED_GENOMES.tsv")
  )
  stop(
    "Only ", length(matched_genomes), " of 497 genomes matched. ",
    "See UNMATCHED_GENOMES.tsv in the output folder."
  )
}

message("All 497 genomes successfully matched.")

# Save the exact identifier mapping used.
write_tsv(
  match_summary,
  file.path(output_dir, "virulencefinder_497_identifier_mapping.tsv")
)

# ------------------------------------------------------------------------------
# 6. CREATE BINARY GENE PRESENCE/ABSENCE MATRIX
# ------------------------------------------------------------------------------

# Repeated hits to the same gene in the same genome are collapsed to one.
hits_unique <- virulence_497 %>%
  transmute(
    genome_id = Matched_REFSEQ,
    gene = str_trim(as.character(`Virulence factor`))
  ) %>%
  filter(
    !is.na(genome_id),
    genome_id != "",
    !is.na(gene),
    gene != ""
  ) %>%
  distinct(genome_id, gene)

presence_absence <- hits_unique %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = gene,
    values_from = present,
    values_fill = 0L
  )

# IMPORTANT: Include all 497 genomes even if a genome had zero retained
# VirulenceFinder hits (normally there should be at least one hit).
presence_absence <- final_497 %>%
  select(REFSEQ) %>%
  rename(genome_id = REFSEQ) %>%
  left_join(presence_absence, by = "genome_id") %>%
  mutate(across(-genome_id, ~ replace_na(.x, 0L))) %>%
  arrange(genome_id)

mat <- presence_absence %>%
  column_to_rownames("genome_id") %>%
  as.matrix()

storage.mode(mat) <- "numeric"

if (nrow(mat) != 497) {
  stop("Presence/absence matrix does not contain exactly 497 genomes.")
}

message("Genomes in heatmap matrix: ", nrow(mat))
message("Virulence genes represented: ", ncol(mat))

write_tsv(
  presence_absence,
  file.path(output_dir, "chapter3_497_virulence_presence_absence.tsv")
)

# ------------------------------------------------------------------------------
# 7. ASSIGN VIRULENCE GENES TO BROAD FUNCTIONAL CATEGORIES
# ------------------------------------------------------------------------------

# This preserves the same category logic used in the Chapter 4 stx2i script.
gene_cat <- tibble(gene = colnames(mat)) %>%
  mutate(
    Category = case_when(
      str_detect(gene, regex("^stx|^elt|^est|^ast|^sub|ehx", ignore_case = TRUE)) ~
        "Toxins",

      str_detect(
        gene,
        regex("fim|lpf|pap|afa|fdeC|iha|hra|tia|yeh|csg|eae|tir", ignore_case = TRUE)
      ) ~ "Adhesins",

      str_detect(
        gene,
        regex("fyu|irp|ire|sit|iuc|iut|chu|feo", ignore_case = TRUE)
      ) ~ "Iron acquisition",

      str_detect(
        gene,
        regex("gad|ter|hha|anr|aamR", ignore_case = TRUE)
      ) ~ "Stress response",

      str_detect(
        gene,
        regex("iss|traT|ompT|kps", ignore_case = TRUE)
      ) ~ "Immune evasion",

      str_detect(
        gene,
        regex("mch|mcm|cba|cia|cma|cea|colE", ignore_case = TRUE)
      ) ~ "Colicins / bacteriocins",

      TRUE ~ "Other"
    )
  )

write_tsv(
  gene_cat,
  file.path(output_dir, "chapter3_497_virulence_gene_categories.tsv")
)

ann_col <- gene_cat %>%
  column_to_rownames("gene")

ann_col <- ann_col[colnames(mat), , drop = FALSE]

# ------------------------------------------------------------------------------
# 8. CREATE GENOME (ROW) ANNOTATION
# ------------------------------------------------------------------------------

# Source1 contains the broad categories used in Chapter 3:
# Meat, Flour, Nuts/Seeds, Dairy, Environmental, Vegetables, Other.
genome_metadata <- final_497 %>%
  transmute(
    genome_id = REFSEQ,
    Source = case_when(
      is.na(Source1) | Source1 == "" ~ "Unknown",
      TRUE ~ Source1
    ),
    Source2 = Source2,
    Serotype = Serotype,
    STX = STX
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

missing_metadata <- setdiff(rownames(mat), genome_metadata$genome_id)

if (length(missing_metadata) > 0) {
  stop(
    "No metadata match was found for the following genome(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

ann_row <- genome_metadata %>%
  select(genome_id, Source) %>%
  column_to_rownames("genome_id")

ann_row <- ann_row[rownames(mat), , drop = FALSE]

write_tsv(
  genome_metadata,
  file.path(output_dir, "chapter3_497_virulence_genome_annotations.tsv")
)

message("Source groups represented:")
print(table(ann_row$Source, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 9. ANNOTATION COLOURS
# ------------------------------------------------------------------------------

# Change these colours if desired.
source_levels <- sort(unique(as.character(ann_row$Source)))

# Named palette for the expected Chapter 3 broad source categories.
source_colours <- c(
  "Meat" = "#D73027",
  "Flour" = "#E6AB02",
  "Nuts and Seeds" = "#A6761D",
  "Nuts/Seeds" = "#A6761D",
  "Dairy" = "#66A61E",
  "Environmental" = "#1B9E77",
  "Environment" = "#1B9E77",
  "Vegetables" = "#7570B3",
  "Other" = "grey50",
  "Unknown" = "grey80"
)

# Automatically assign colours if an unexpected Source1 category is present.
missing_source_colours <- setdiff(source_levels, names(source_colours))
if (length(missing_source_colours) > 0) {
  extra_cols <- setNames(
    grDevices::hcl.colors(length(missing_source_colours), palette = "Dark 3"),
    missing_source_colours
  )
  source_colours <- c(source_colours, extra_cols)
}

source_colours <- source_colours[source_levels]

annotation_colors <- list(
  Source = source_colours,
  Category = c(
    "Adhesins" = "#1B9E77",
    "Colicins / bacteriocins" = "#E7298A",
    "Immune evasion" = "#66A61E",
    "Iron acquisition" = "#E69F00",
    "Other" = "grey50",
    "Stress response" = "#7570B3",
    "Toxins" = "#D73027"
  )
)

# ------------------------------------------------------------------------------
# 10. ORDER GENES BY CATEGORY AND PREVALENCE
# ------------------------------------------------------------------------------

prevalence <- colSums(mat > 0)

category_order <- c(
  "Adhesins",
  "Colicins / bacteriocins",
  "Immune evasion",
  "Iron acquisition",
  "Other",
  "Stress response",
  "Toxins"
)

gene_order <- gene_cat %>%
  mutate(
    Category = factor(Category, levels = category_order),
    prevalence = prevalence[gene]
  ) %>%
  arrange(Category, desc(prevalence), gene) %>%
  pull(gene)

mat_ordered <- mat[, gene_order, drop = FALSE]
ann_col_ordered <- ann_col[colnames(mat_ordered), , drop = FALSE]

# ------------------------------------------------------------------------------
# 11. DRAW AND SAVE THE HEATMAP
# ------------------------------------------------------------------------------

heatmap_colours <- c("grey95", "#2166AC")
heatmap_breaks <- c(-0.5, 0.5, 1.5)

# For 497 genomes, individual genome labels are not readable at thesis-figure
# size. REFSEQ is still the row name used internally, but labels are hidden.
# Change show_rownames = FALSE to TRUE if you want a very large labelled figure.
draw_heatmap <- function() {
  pheatmap(
    mat_ordered,
    color = heatmap_colours,
    breaks = heatmap_breaks,
    annotation_col = ann_col_ordered,
    annotation_row = ann_row,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    clustering_method = "ward.D2",
    border_color = NA,
    show_rownames = FALSE,
    fontsize_col = 7,
    angle_col = 45,
    legend_breaks = c(0, 1),
    legend_labels = c("Absent", "Present"),
    main = expression(
      "Virulence gene profiles of 497 " * italic(E.~coli) * " genomes"
    )
  )
}

# JPEG output
jpeg(
  file.path(output_dir, "chapter3_497_virulence_heatmap.jpg"),
  width = 18,
  height = 14,
  units = "in",
  res = 300
)
draw_heatmap()
dev.off()

# PDF output (recommended for thesis-quality vector output)
pdf(
  file.path(output_dir, "chapter3_497_virulence_heatmap.pdf"),
  width = 18,
  height = 14
)
draw_heatmap()
dev.off()

# ------------------------------------------------------------------------------
# 12. SAVE GENE PREVALENCE SUMMARY
# ------------------------------------------------------------------------------

gene_prevalence <- gene_cat %>%
  mutate(
    genomes_positive = as.integer(prevalence[gene]),
    prevalence_percent = round(100 * genomes_positive / nrow(mat), 1)
  ) %>%
  arrange(desc(genomes_positive), gene)

write_tsv(
  gene_prevalence,
  file.path(output_dir, "chapter3_497_virulence_gene_prevalence.tsv")
)

# ------------------------------------------------------------------------------
# 13. FINAL CHECKS
# ------------------------------------------------------------------------------

message("------------------------------------------------------------")
message("Analysis complete.")
message("Genomes in final heatmap: ", nrow(mat))
message("Virulence genes in final heatmap: ", ncol(mat))
message("REFSEQ is used as the genome identifier/row name.")
message("Heatmap row labels are hidden because there are 497 genomes.")
message("Files written to: ", output_dir)
message("------------------------------------------------------------")
