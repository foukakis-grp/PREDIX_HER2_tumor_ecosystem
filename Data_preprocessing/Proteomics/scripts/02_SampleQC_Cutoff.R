

# Setup -------------------------------------------------------------------


# Libraries
library(vroom)
library(readxl)
library(janitor)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(factoextra)
library(reshape2)
library(ggpubr)


# Load functions
source("scripts/setup.R")

# Sample cutoff: 80th percentile of 187 samples with metadata: [quantile(seq(1, 187, 1), probs = 0.8)] = 149.8 -> 150 samples (minimal number of identified PGs: 5723)


# Load files --------------------------------------------------------------


# Load metdata
Metadata_V4 <- readRDS(file = "processed/Metadata_Predix_pretherapy_V4_MSruninfo.rds")

# Load quant data
Prot_Quant_wide_V1 <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V1.rds")


# Visualize identifications of precursors, peptides, proteins ---------------------------------------------


# MS identification data
MS_Ident_Vis <- Metadata_V4

# Ranked
MS_Ident_Vis$Sample_ID_MS <- forcats::fct_reorder(MS_Ident_Vis$Sample_ID_MS, MS_Ident_Vis$ProteinGroups_MS, .desc = TRUE)

Add_plot_title <- " (sorted by protein group number)"


# Bar plot for identified precursor number per sample
BarPrec_All <- ggplot(MS_Ident_Vis,
                      aes(x = Sample_ID_MS, y = Precursors_MS, fill = Sample_ID_MS)) +
  geom_bar(stat = "identity", fill = 'lightblue', color = 'black', linewidth = Borderwidth) +
  labs(title = paste0("Precursor identifications per sample", Add_plot_title),
       y = "No. of precursors",
       x = "Sample") +
  scale_y_continuous(breaks = seq(0, ceiling(max(MS_Ident_Vis$Precursors_MS)/10000)*10000, 10000), limits = c(0, ceiling(max(MS_Ident_Vis$Precursors_MS)/10000)*10000), expand = expansion(mult = c(0, 0.01))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = Axis_title),
        axis.text.y = element_text(size = Axis_text))

#BarPrec_All


# Bar plot for stripped peptide number per sample
BarPep_All <- ggplot(MS_Ident_Vis,
                     aes(x = Sample_ID_MS, y = Peptides_MS, fill = Sample_ID_MS)) +
  geom_bar(stat = "identity", fill = 'indianred', color = 'black', linewidth = Borderwidth) +
  labs(title = paste0("Stripped sequence identifications per sample", Add_plot_title),
       y = "No. of stripped sequences",
       x = "Sample") +
  scale_y_continuous(breaks = seq(0, ceiling(max(MS_Ident_Vis$Peptides_MS)/10000)*10000, 10000), limits = c(0, ceiling(max(MS_Ident_Vis$Peptides_MS)/10000)*10000), expand = expansion(mult = c(0, 0.01))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = Axis_title),
        axis.text.y = element_text(size = Axis_text))

#BarPep_All

# Bar plot for identified protein number per sample
BarProt_All <- ggplot(MS_Ident_Vis,
                      aes(x = Sample_ID_MS, y = ProteinGroups_MS, fill = Sample_ID_MS)) +
  geom_bar(stat = "identity", fill = 'palegreen3', color = 'black', linewidth = Borderwidth) +
  geom_hline(yintercept = 5723, linetype = "dashed", color = "black", linewidth = Linewidth) +
  labs(title = paste0("Protein group identifications per sample", Add_plot_title),
       y = "No. of protein groups",
       x = "Sample") +
  scale_y_continuous(breaks = seq(0, ceiling(max(MS_Ident_Vis$ProteinGroups_MS)/1000)*1000, 1000), limits = c(0, ceiling(max(MS_Ident_Vis$ProteinGroups_MS)/1000)*1000), expand = expansion(mult = c(0, 0.01))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = Axis_title),
        axis.text.y = element_text(size = Axis_text))

#BarProt_All

# Combined plot
Bar_Comb_All <- BarPrec_All + BarPep_All + BarProt_All +
  plot_layout(ncol = 1, axes = "collect")

Bar_Comb_All_mod <- Bar_Comb_All & labs(title = NULL, subtitle = NULL) & theme(plot.title = element_blank(), axis.title = element_text(size = Axis_title),
                                                                               axis.text.y = element_text(size = Axis_text), axis.text.x = element_blank())

#Bar_Comb_All_mod

#ggsave(plot = Bar_Comb_All_mod, paste0("figures/Sample_QC_Cutoff/MS_Entities_Ranked_PREDIX_pretherapy_Barplot", ".png"), device = "png", unit = "cm", width = 16, height = 13, dpi = 300)
#ggsave(plot = Bar_Comb_All_mod, paste0("figures/Sample_QC_Cutoff/MS_Entities_Ranked_PREDIX_pretherapy_Barplot", ".pdf"), device = "pdf", unit = "cm", width = 16, height = 13)


# Peptides for quantification ----------------------------------------------------------


# Load long-format Spectronaut data (measured pre-therapy samples)
Protein_Quant_Spec_long <- vroom(file = "data/Predix_directdia_GS_ENS106_BaselineSamples_withre-inj-Report_HJ_Protein_Quant_long (Normal).tsv",
                                 delim = "\t",
                                 col_select = c(R.FileName, R.ProteinGroupsIdentified, PG.NrOfStrippedSequencesUsedForQuantification))

# Filter samples according to clean metadata
Protein_Quant_long_temp <- Protein_Quant_Spec_long %>%
  filter(R.FileName %in% Metadata_V4$FileName_MS)

# Select columns of interest & rename columns
Peptide_used_quant <- Protein_Quant_long_temp %>%
  dplyr::select(R.FileName, PG.NrOfStrippedSequencesUsedForQuantification)

Peptide_used_quant <- Peptide_used_quant %>%
  dplyr::rename("Peptide_quant" = "PG.NrOfStrippedSequencesUsedForQuantification", "Sample" = "R.FileName")

# Count per sample how often 0/1/2/3 peptides were used for quantification
Peptide_count_table <- as.data.frame(table(Peptide_used_quant$Sample, Peptide_used_quant$Peptide_quant))
Peptide_count_table[is.na(Peptide_count_table)] <- 0

colnames(Peptide_count_table) <- c("Sample", "Peptide_Nr", "Frequency")

# Remove "0 peptide" entries
Peptide_count_table <- Peptide_count_table %>%
  filter(Peptide_Nr != 0)

# Arrange factor levels (samples sorted by protein group)
Peptide_count_table$Sample <- factor(Peptide_count_table$Sample, levels = levels(forcats::fct_reorder(Metadata_V4$FileName_MS, Metadata_V4$ProteinGroups_MS, .desc = TRUE)))
Peptide_count_table$Peptide_Nr <- factor(Peptide_count_table$Peptide_Nr, levels = rev(c(3, 2, 1)))

## Convert absolute numbers to proportions
Peptide_count_table <- Peptide_count_table %>%
  group_by(Sample) %>%
  mutate(Total_prot_Nr = sum(Frequency)) %>%
  ungroup() %>%
  mutate(Prop_of_total_prot = Frequency/Total_prot_Nr)

## Plot how often 1/2/3 peptides were used for quantification per sample (in proportion to total protein group number)
Plot_peptides_for_quant <- ggplot(Peptide_count_table,
                                  aes(x = Sample, y = Prop_of_total_prot*100, fill = Peptide_Nr)) +
  geom_bar(stat = "identity", width = 1.2) + #remove spaces between columns
  geom_vline(xintercept = 150.5, linetype = "dashed", linewidth = Linewidth, color = "black") + # insert sample cutoff line
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  labs(title = paste0("Proportion of peptides used for protein quantification per sample"),
       y = "Proportion of protein groups [%]",
       x = "Sample") +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = Axis_title),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = Axis_title),
        axis.text.y = element_text(size = Axis_text),
        legend.title = element_text(size = Axis_title),
        legend.text = element_text(size = Axis_title)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100.1)) +
  scale_x_discrete(expand = c(0, 0, 0, 0)) +
  coord_cartesian(ylim = c(0, 100), expand = FALSE) +
  guides(fill = guide_legend(title = "Peptides used \nfor quantification"))

#Plot_peptides_for_quant

#ggsave(plot = Plot_peptides_for_quant, paste0("figures/Sample_QC_Cutoff/Peptides_For_Quant_PREDIX_pretherapy_Barplot.png"), device = "png", unit = "cm", width = 16, height = 13, dpi = 300)
#ggsave(plot = Plot_peptides_for_quant, paste0("figures/Sample_QC_Cutoff/Peptides_For_Quant_PREDIX_pretherapy_Barplot.pdf"), device = "pdf", unit = "cm", width = 16, height = 13)


# Apply sample cutoff -----------------------------------------------------------


# Metdata

## Remove low-quality samples
Metadata_V5 <- Metadata_V4 %>%
  arrange(desc(ProteinGroups_MS)) %>% #arrange by identified protein number
  mutate(Rank = seq(1, nrow(.), 1)) %>% #rank by identified protein number
  filter(Rank <= 150) #filter samples according to the cutoff defined above

## Save
#saveRDS(Metadata_V5, file = "processed/Metadata_Predix_pretherapy_V5_HighQuality.rds")

# Quant data

## Remove low-quality samples
Prot_Quant_wide_V2 <- Prot_Quant_wide_V1 %>%
  dplyr::select(Metadata_V5$Sample_ID_MS)

## Save
#saveRDS(Prot_Quant_wide_V2, file = "processed/Protein_Quant_Predix_pretherapy_wide_V2_HighQuality.rds")


# Session End ------------------------------------------------------------

