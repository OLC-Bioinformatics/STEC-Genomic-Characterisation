# ==============================================================================
# Figure 3.3  Virulence gene prevalence in O157:H7 versus non-O157 genomes
# ------------------------------------------------------------------------------
# Replaces the source-stratified prevalence heatmap. Two changes from the earlier
# version:
#   1. Rows are O157:H7 vs non-O157 rather than isolation source. The virulence
#      gene content of this collection is structured by lineage, not commodity.
#   2. Only genes present in 5-95% of genomes are shown. Six genes are in >95%
#      of genomes and 107 in <5%; neither group distinguishes anything. stx genes
#      are also dropped because they are covered by Figures 3.1 and 3.2.
# ==============================================================================

library(tidyverse)
library(pheatmap)
library(readxl)

metadata_file <- "/Users/noorshubair/Desktop/R_code_files/chapter3/CFIA_STEC.xlsx"
pa_file       <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap/chapter3_497_virulence_presence_absence.tsv"
cat_file      <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap/chapter3_497_virulence_gene_categories.tsv"
output_dir    <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap"

MIN_PREV <- 5     # drop near-absent genes
MAX_PREV <- 95    # drop near-universal genes

pa       <- read_tsv(pa_file, show_col_types = FALSE)
gene_cat <- read_tsv(cat_file, show_col_types = FALSE)
final497 <- read_excel(metadata_file, sheet = "prjna454819") %>%
  transmute(genome_id = str_trim(REFSEQ),
            Serotype  = str_trim(`Serotype_HC list`),
            Otype     = str_trim(Otype)) %>%
  mutate(Group = if_else(Otype == "O157", "O157:H7", "non-O157"))

mat <- pa %>% column_to_rownames("genome_id") %>% as.matrix()
storage.mode(mat) <- "numeric"

prev <- 100 * colSums(mat > 0) / nrow(mat)
keep <- names(prev)[prev >= MIN_PREV & prev <= MAX_PREV & !str_detect(names(prev), regex("^stx", ignore_case = TRUE))]
message("Genes retained: ", length(keep), " of ", ncol(mat))

grp <- final497$Group[match(rownames(mat), final497$genome_id)]
prev_by_group <- t(sapply(c("O157:H7", "non-O157"),
                          function(g) 100 * colMeans(mat[grp == g, keep, drop = FALSE] > 0)))
rownames(prev_by_group) <- c(sprintf("O157:H7 (n = %d)",  sum(grp == "O157:H7")),
                             sprintf("non-O157 (n = %d)", sum(grp == "non-O157")))

# order genes by functional category, then by how strongly they separate the groups
CAT_ORDER <- c("Adhesins","Toxins","Immune evasion","Iron acquisition",
               "Stress response","Colicins / bacteriocins","Other")
gene_order <- tibble(gene = keep) %>%
  left_join(gene_cat, by = "gene") %>%
  mutate(Category = factor(replace_na(Category, "Other"), levels = CAT_ORDER),
         diff = prev_by_group[1, gene] - prev_by_group[2, gene]) %>%
  arrange(Category, desc(diff), gene) %>%
  pull(gene)

prev_by_group <- prev_by_group[, gene_order, drop = FALSE]
ann_col <- gene_cat %>% filter(gene %in% gene_order) %>%
  mutate(Category = factor(replace_na(Category, "Other"), levels = CAT_ORDER)) %>%
  column_to_rownames("gene")
ann_col <- ann_col[gene_order, , drop = FALSE]

annotation_colors <- list(Category = c(
  "Adhesins" = "#1B9E77", "Toxins" = "#D73027", "Immune evasion" = "#66A61E",
  "Iron acquisition" = "#E69F00", "Stress response" = "#7570B3",
  "Colicins / bacteriocins" = "#E7298A", "Other" = "grey50"))

draw <- function() {
  pheatmap(prev_by_group,
           color = colorRampPalette(c("white","#DEEBF7","#9ECAE1","#3182BD","#08519C"))(100),
           breaks = seq(0, 100, length.out = 101),
           annotation_col = ann_col, annotation_colors = annotation_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, border_color = "white",
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 11, fontsize_col = 8, angle_col = 90,
           legend_breaks = c(0,25,50,75,100),
           legend_labels = c("0%","25%","50%","75%","100%"),
           main = "")
}

pdf(file.path(output_dir, "Figure_3.3_virulence_prevalence_O157.pdf"), width = 14, height = 3.6); draw(); dev.off()
jpeg(file.path(output_dir, "Figure_3.3_virulence_prevalence_O157.jpg"),
     width = 14, height = 3.6, units = "in", res = 300); draw(); dev.off()

write_tsv(as_tibble(prev_by_group, rownames = "Group"),
          file.path(output_dir, "Figure_3.3_underlying_values.tsv"))
message("Done.")
