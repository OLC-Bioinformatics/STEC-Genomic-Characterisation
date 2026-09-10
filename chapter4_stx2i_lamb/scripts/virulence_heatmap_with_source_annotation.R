# ==============================================================================
# VirulenceFinder heatmap with genome-source annotation
# 66 stx2i-positive Escherichia coli genomes
# ==============================================================================
# Inputs (stored in vir_dir and data_dir):
#   1. combined_virulence_results.tsv
#   2. Stx2i_NCBI_CFIA_genomes_working_260730.xlsx
#
# Outputs:
#   - virulence_presence_absence.tsv
#   - virulence_gene_categories.tsv
#   - virulence_gene_prevalence.tsv
#   - virulence_genome_source_annotations.tsv
#   - virulence_presence_absence_heatmap_with_source.jpg
#   - virulence_presence_absence_heatmap_with_source.pdf
#
# Required packages:
# install.packages(c("tidyverse", "pheatmap", "readxl"))
# ============================================================================== 

library(tidyverse)
library(pheatmap)
library(readxl)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------------------------

# Change this path if your files are stored elsewhere.
vir_dir <- "/Users/noorshubair/Desktop/chapter4/virulencefinder/virulencefinder_results"
data_dir <- "/Users/noorshubair/Desktop/chapter4/data/"

input_file <- file.path(vir_dir, "combined_virulence_results.tsv")
metadata_file <- file.path(data_dir, "Stx2i_NCBI_CFIA_genomes_working_260730.xlsx")

if (!file.exists(input_file)) {
  stop("VirulenceFinder input file not found: ", input_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata workbook not found: ", metadata_file)
}

# ------------------------------------------------------------------------------
# 2. Read the combined VirulenceFinder output
# ------------------------------------------------------------------------------

virulence_results <- read_tsv(
  input_file,
  show_col_types = FALSE,
  trim_ws = TRUE
)

required_columns <- c("Genome_ID", "Virulence factor")
missing_columns <- setdiff(required_columns, colnames(virulence_results))

if (length(missing_columns) > 0) {
  stop(
    "The following required column(s) are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

n_genomes <- n_distinct(virulence_results$Genome_ID)
message("Genomes represented in combined results: ", n_genomes)

if (n_genomes != 66) {
  warning("Expected 66 genomes, but found ", n_genomes, " genome IDs.")
}

# ------------------------------------------------------------------------------
# 3. Create a binary gene presence/absence matrix
# ------------------------------------------------------------------------------

# Repeated hits to the same gene in the same genome are collapsed to one.
hits_unique <- virulence_results %>%
  transmute(
    genome_id = str_trim(as.character(Genome_ID)),
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
  ) %>%
  arrange(genome_id)

mat <- presence_absence %>%
  column_to_rownames("genome_id") %>%
  as.matrix()

storage.mode(mat) <- "numeric"
message("Virulence genes represented: ", ncol(mat))

write_tsv(
  presence_absence,
  file.path(vir_dir, "virulence_presence_absence.tsv")
)

# ------------------------------------------------------------------------------
# 4. Assign genes to broad functional categories
# ------------------------------------------------------------------------------

gene_cat <- tibble(gene = colnames(mat)) %>%
  mutate(
    category = case_when(
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
  file.path(vir_dir, "virulence_gene_categories.tsv")
)

gene_cat <- gene_cat %>%
  rename(Category = category)

ann_col <- gene_cat %>%
  column_to_rownames("gene")

ann_col <- ann_col[colnames(mat), , drop = FALSE]

# ------------------------------------------------------------------------------
# 5. Read metadata and create row annotations
# ------------------------------------------------------------------------------

metadata_raw <- read_excel(
  metadata_file,
  sheet = "NCBI_CFIA_dataset"
)

# Remove leading/trailing spaces from workbook column names.
colnames(metadata_raw) <- str_trim(colnames(metadata_raw))

metadata_required <- c(
  "ASSEMBLY TO USE",
  "Host scientific name",
  "Isolation source",
  "Isolation type",
  "Source Type 2",
  "MLST",
  "EC-Typer"
)

metadata_missing <- setdiff(metadata_required, colnames(metadata_raw))

if (length(metadata_missing) > 0) {
  stop(
    "The following required metadata column(s) are missing: ",
    paste(metadata_missing, collapse = ", ")
  )
}

genome_metadata <- metadata_raw %>%
  transmute(
    genome_id = str_trim(as.character(`ASSEMBLY TO USE`)),
    host = str_trim(as.character(`Host scientific name`)),
    isolation_source = str_trim(as.character(`Isolation source`)),
    isolation_type = str_trim(as.character(`Isolation type`)),
    source_type = str_trim(as.character(`Source Type 2`)),
    ST = as.character(MLST),
    serotype = str_trim(as.character(`EC-Typer`))
  ) %>%
  filter(!is.na(genome_id), genome_id != "") %>%
  distinct(genome_id, .keep_all = TRUE) %>%
  mutate(
    source_text = str_c(
      replace_na(host, ""),
      replace_na(isolation_source, ""),
      replace_na(isolation_type, ""),
      replace_na(source_type, ""),
      sep = " | "
    ),
    Source = source_type,
  
    ST = if_else(
      is.na(ST) | ST == "" | ST == "NA",
      "Unknown",
      paste0("ST", ST)
    )
  )

# Retain only the 66 genomes shown in the heatmap.
genome_metadata <- genome_metadata %>%
  filter(genome_id %in% rownames(mat))

missing_metadata <- setdiff(rownames(mat), genome_metadata$genome_id)

if (length(missing_metadata) > 0) {
  stop(
    "No metadata match was found for the following heatmap genome(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (n_distinct(genome_metadata$genome_id) != nrow(mat)) {
  stop(
    "Metadata matching did not produce exactly one record per heatmap genome. ",
    "Please check duplicate or missing genome IDs."
  )
}

# Source is the main row annotation
ann_row <- genome_metadata %>%
  select(genome_id, Source) %>%
  column_to_rownames("genome_id")

ann_row <- ann_row[rownames(mat), , drop = FALSE]


#------------------------------------------------------------------------------

# Annotation colours

#------------------------------------------------------------------------------
sort(unique(genome_metadata$source_type))

annotation_colors <- list(
  Source = c(
    "Human Clinical" = "#D73027",
    "Ovis aries" = "#1A9850",
    "Oyster" = "lightblue",
    "Dairy" = "#FDDA24",
    "Environmental" = "#FF8200",
    "Raw Animal Feed" = "#7B3294",
    "Unknown Food" = "grey50" 
  ),
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

write_tsv(
  genome_metadata %>%
    select(
      genome_id,
      Source,
      serotype,
      host,
      isolation_source,
      isolation_type,
      source_type
    ),
  file.path(vir_dir, "virulence_genome_source_annotations.tsv")
)

message("Source groups represented:")
print(table(ann_row$Source, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 6. Order genes by category and prevalence
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
# 7. Draw and save the annotated heatmap
# ------------------------------------------------------------------------------

heatmap_colours <- c("grey95", "#2166AC")
heatmap_breaks <- c(-0.5, 0.5, 1.5)

# Use the same plot settings for both JPEG and PDF outputs.
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
    border_color = "grey",
    fontsize_row = 7,
    fontsize_col = 8,
    angle_col = 45,
    legend_breaks = c(0, 1),
    legend_labels = c("Absent", "Present"),
    main = expression(
      "Virulence gene profiles of " *
        italic(stx2i) * "-positive " *
        italic(E.~coli) * " genomes"
    )
  )
}

# JPEG output
dpi <- 300
jpeg("/Users/noorshubair/Desktop/R_code_files/chapter4/virulence_presence_absence_heatmap_with_source_rowcluster.jpg",
  width = 16 * dpi,
  height = 12 * dpi,
  res = dpi
)

draw_heatmap()
dev.off()

# PDF output
pdf(
  file.path(vir_dir, "virulence_presence_absence_heatmap_with_source.pdf"),
  width = 16,
  height = 12
)

draw_heatmap()
dev.off()

# ------------------------------------------------------------------------------
# 8. Save a gene-prevalence summary
# ------------------------------------------------------------------------------

gene_prevalence <- gene_cat %>%
  mutate(
    genomes_positive = as.integer(prevalence[gene]),
    prevalence_percent = round(100 * genomes_positive / nrow(mat), 1)
  ) %>%
  arrange(desc(genomes_positive), gene)

write_tsv(
  gene_prevalence,
  file.path(vir_dir, "virulence_gene_prevalence.tsv")
)

message("Analysis complete.")
message("Files written to: ", vir_dir)

