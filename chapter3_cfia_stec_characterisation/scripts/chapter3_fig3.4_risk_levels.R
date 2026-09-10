# ==============================================================================
# Figure 3.4  FAO/WHO risk levels of the 497 CFIA STEC genomes
# ------------------------------------------------------------------------------
# Risk levels follow FAO and WHO (2018) as applied in Table 2 of Zhang et al.
# (2025):  1 = stx2a with eae or aggR;  2 = stx2d;  3 = stx2c with eae;
#          4 = stx1a with eae;  5 = any other stx subtype.
# Where a genome meets more than one criterion the highest (most severe) level
# is assigned, following the worked examples in Table 3 of that paper.
# aggR was not detected in any genome in this collection.
# ==============================================================================

library(tidyverse)
library(readxl)

metadata_file <- "/Users/noorshubair/Desktop/R_code_files/chapter3/CFIA_STEC.xlsx"
pa_file       <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap/chapter3_497_virulence_presence_absence.tsv"
output_dir    <- "/Users/noorshubair/Desktop/chapter3/virulence_heatmap"

final497 <- read_excel(metadata_file, sheet = "prjna454819") %>%
  transmute(genome_id = str_trim(REFSEQ),
            Otype     = str_trim(Otype),
            stx_raw   = str_trim(`Shiga Toxin Gene(s)`))

pa  <- read_tsv(pa_file, show_col_types = FALSE)
eae <- pa %>%
  transmute(genome_id,
            eae = as.integer(rowSums(across(starts_with("eae"))) > 0))

dat <- final497 %>%
  left_join(eae, by = "genome_id") %>%
  mutate(
    eae   = replace_na(eae, 0L),
    stx2a = str_detect(stx_raw, "Stx2a"),
    stx2c = str_detect(stx_raw, "Stx2c"),
    stx2d = str_detect(stx_raw, "Stx2d"),
    stx1a = str_detect(stx_raw, "Stx1a"),
    n_stx = str_count(stx_raw, "Stx[12][a-o]"),
    Group = if_else(Otype == "O157", "O157:H7", "non-O157"),
    Risk  = case_when(
      n_stx == 0            ~ "Not classifiable",
      stx2a & eae == 1      ~ "Level 1",
      stx2d                 ~ "Level 2",
      stx2c & eae == 1      ~ "Level 3",
      stx1a & eae == 1      ~ "Level 4",
      TRUE                  ~ "Level 5"))

LV  <- c("Level 1","Level 2","Level 3","Level 4","Level 5","Not classifiable")
COL <- c("Level 1"="#A50026","Level 2"="#F46D43","Level 3"="#FDAE61",
         "Level 4"="#FEE090","Level 5"="#4575B4","Not classifiable"="grey75")

plot_dat <- bind_rows(dat %>% mutate(Panel = Group),
                      dat %>% mutate(Panel = "All genomes")) %>%
  count(Panel, Risk) %>%
  group_by(Panel) %>% mutate(total = sum(n)) %>% ungroup() %>%
  mutate(Risk  = factor(Risk, levels = LV),
         Panel = factor(sprintf("%s\n(n = %d)", Panel, total),
                        levels = rev(unique(sprintf("%s\n(n = %d)",
                              c("O157:H7","non-O157","All genomes"),
                              c(sum(dat$Group=="O157:H7"), sum(dat$Group=="non-O157"), nrow(dat)))))))

p <- ggplot(plot_dat, aes(x = n, y = Panel, fill = Risk)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n / total > 0.045, n, "")),
            position = position_stack(vjust = 0.5), size = 3.4, fontface = "bold") +
  scale_fill_manual(values = COL, drop = FALSE) +
  labs(x = "Number of genomes", y = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(output_dir, "Figure_3.4_risk_levels.pdf"), p, width = 9.2, height = 3.3)
ggsave(file.path(output_dir, "Figure_3.4_risk_levels.jpg"), p, width = 9.2, height = 3.3, dpi = 300)

print(dat %>% count(Group, Risk) %>% pivot_wider(names_from = Group, values_from = n, values_fill = 0))
write_tsv(dat %>% select(genome_id, Otype, Group, stx_raw, eae, Risk),
          file.path(output_dir, "Figure_3.4_risk_assignments.tsv"))
message("Done.")
