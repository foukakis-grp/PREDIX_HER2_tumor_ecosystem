

# Setup -------------------------------------------------------------------


# Libraries
library(janitor)
library(tidyverse)
library(reshape2)
library(mixOmics)
library(RColorBrewer)
library(patchwork)


# Load functions
source("scripts/setup.R")

#Performed for quant data versions with & without batch effect correction


# Load files --------------------------------------------------------------


# Metadata
Metadata <- readRDS(file = "processed/Metadata_Predix_pretherapy_V5_HighQuality.rds")

# Quant files

## No batch effect correction
Prot_Quant_noCombat <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V4_Imp.rds")

## Batch effect correction
Prot_Quant_Combat <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V5_ComBat.rds")


# Remove biological replicates --------------------------------------------


# Remove biological replicate with lower number of identified proteins
Metadata_noRep <- Metadata %>%
  group_by(Number_MS) %>%
  filter(ProteinGroups_MS == max(ProteinGroups_MS)) %>%
  ungroup()

## Save
#saveRDS(Metadata_noRep, file = "processed/Metadata_Predix_pretherapy_V6_NoReps.rds")

# Adjust quant data

## Define input
Quant_list <- list(

  "noCombat" = Prot_Quant_noCombat,
  "Combat" = Prot_Quant_Combat

)

## Loop: remove replicates & center by median

Quant_output_list <- lapply(Quant_list, function(Quant_temp){

  # Remove replicates
  Quant_temp_noRep <- Quant_temp %>%
    dplyr::select(Metadata_noRep$Sample_ID_MS)

  Quant_temp_noRep <- as.matrix(Quant_temp_noRep)

  Quant_temp_noRep_shift <- Quant_temp_noRep - median(Quant_temp_noRep)

  Quant_output <- as.data.frame(Quant_temp_noRep_shift)

})

## Save output

saveRDS(Quant_output_list[["noCombat"]], file = "processed/Protein_Quant_Predix_pretherapy_wide_V6_1_NoComBat_NoReps_MedCent.rds")

saveRDS(Quant_output_list[["Combat"]], file = "processed/Protein_Quant_Predix_pretherapy_wide_V6_2_ComBat_NoReps_MedCent.rds")


# QC plot: batch effect correction -----------------------------------------------------------------

## Load files --------------------------------------------------------------


# Metadata
Metadata_noReps <- readRDS(file = "processed/Metadata_Predix_pretherapy_V6_NoReps.rds")

# Quant data
Prot_Quant_noCombat_noReps <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V6_1_NoComBat_NoReps_MedCent.rds")

Prot_Quant_Combat_noReps <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V6_2_ComBat_NoReps_MedCent.rds")


## Individual plots --------------------------------------------------------


# Input files
Quant_matrices <- c("Prot_Quant_noCombat_noReps", "Prot_Quant_Combat_noReps")

# Plot list to save ggplot2 objects per quant file
plot_list <- list()

# LOOP: adjustment procedure quality control

for(file_OI in Quant_matrices){ #file_OI <- "Prot_Quant_Combat_noReps"

  Quant_input_temp_log2 <- get(file_OI)

  # Plot name
  Status <- str_match(file_OI, "Prot_Quant_(\\w+)")[,2]

  cat("Processing quant file:", Status, "\n")

  # Ensure identical quant and metadata order
  Quant_input_temp_log2 <- Quant_input_temp_log2[, Metadata_noReps$Sample_ID_MS]

  # Total number of protein rows
  Total_proteins <- as.character(nrow(Quant_input_temp_log2))

  # Number of protein rows without missing values
  Proteins_wo_NAs <- as.character(nrow(Quant_input_temp_log2 %>% na.omit()))


  # Inter- vs intra-batch correlation

  ## Get correlation matrix of all samples (Spearman: log2 transformation (= monotonic) has no influence on correlation calculation)
  cor_samples <- cor(Quant_input_temp_log2, method = 'spearman', use = "pairwise.complete.obs")

  ## Get only one one half of the correlation matrix (to prevent duplicate values) and remove diagonal
  cor_samples_half <- cor_samples
  cor_samples_half[lower.tri(cor_samples_half, diag = TRUE)] <- NA

  ## Divide into batches & replicates

  ### Convert into long format
  cor_samples_long <- melt(cor_samples_half) %>%
    na.omit() #remove duplicate rows

  ### Add batch & replicate info from metadata
  cor_samples_long_info <- cor_samples_long %>%
    left_join(Metadata_noReps %>% dplyr::select(Sample_ID_MS, Plate_MS, Number_MS), by = c("Var1" = "Sample_ID_MS")) %>%
    left_join(Metadata_noReps %>% dplyr::select(Sample_ID_MS, Plate_MS, Number_MS), by = c("Var2" = "Sample_ID_MS"), suffix = c("_1", "_2")) %>%
    mutate(Replicates = case_when( #highlight replicates
      Number_MS_1 == Number_MS_2 ~ "Replicates",
      TRUE ~ "No Replicates"
    )) %>%
    mutate(Batches = case_when( #highlight batches
      Plate_MS_1 == Plate_MS_2 ~ "Same batch",
      TRUE ~ "Different batch"
    ))

  ### Wilcoxon rank-sum test
  Wilcoxon_test <- cor_samples_long_info %>%
    rstatix::wilcox_test(value ~ Batches) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_xy_position()

  ### Plot pairwise correlations grouped by category (same/different batch) in violin plot

  ### Adjusted plot for combined figure downstream

  Sample_correlations_violin_COMB <- ggplot(cor_samples_long_info,
                                            aes(x = Batches, y = value)) +
    geom_violin(width = 1, fill = "grey", linewidth = Borderwidth) +
    geom_boxplot(width = 0.1, outlier.size = Point_size, linewidth = Borderwidth) +
    stat_summary(fun = median, geom = "label", aes(label = sprintf("%.3f", after_stat(y))), vjust = -0.75, size = Axis_text/ggplot2::.pt, color = "black") +
    scale_x_discrete(labels = scales::label_wrap(6)) +
    scale_y_continuous(breaks = scales::breaks_width(0.1)) +
    labs(title = paste0("Pairwise spearman correlation of ", Status, " samples"),
         subtitle = paste0("Total protein rows: ", Total_proteins, " / Protein rows without NAs: ", Proteins_wo_NAs),
         y = "Correlation coefficient",
         x = "") +
    theme_bw() +
    theme(legend.position = "none") +
    ggpubr::stat_pvalue_manual(Wilcoxon_test,
                               label = "p-value = {signif(p, digits = 3)}",
                               label.size = Axis_text/ggplot2::.pt)

  #Sample_correlations_violin_COMB


  # Non-Linear Iterative Square (NIPALS) PCA
  result.pca.multi <- mixOmics::pca(t(Quant_input_temp_log2), scale = FALSE)

  Metadata_noReps$Plate_MS <- factor(Metadata_noReps$Plate_MS, levels = stringr::str_sort(unique(Metadata_noReps$Plate_MS), numeric = TRUE, locale = "en_UK"))


  ## Adjusted plot for combined figure downstream
  set.seed(1234) #seed for reproducibility of color palette

  PCA_with_NAs_COMB <- plotIndiv(result.pca.multi,
                                 style = "ggplot2",
                                 col.per.group = (sample(colorRampPalette(brewer.pal(8, "Set1"))(10))),
                                 comp = c(1,2),
                                 cex = 1.5,
                                 ind.names = TRUE,
                                 group = Metadata_noReps$Plate_MS,
                                 pch = Metadata_noReps$Plate_MS,
                                 legend.position = "bottom",
                                 title = "",
                                 legend = TRUE, legend.title = "Preparation batch")

  ### ggplot object can be saved & modified
  PCA_with_NAs_ggplot_COMB <- PCA_with_NAs_COMB$graph +
    theme_bw() +
    theme(legend.title = element_text(size = Axis_title),
          legend.text = element_text(size = Axis_title),
          strip.background = element_blank(),
          strip.text = element_blank())

  ### Change point size post-plotting
  for(n in seq(1,length(PCA_with_NAs_ggplot_COMB$layers),1)){

    PCA_with_NAs_ggplot_COMB$layers[[n]]$aes_params$size <- 1.5

    PCA_with_NAs_ggplot_COMB$layers[[n]]$aes_params$stroke <- 0.75

  }

  ### Save ggplot2 objects for combined plot
  plot_list[[paste0(Status, "_Sample_Cor")]] <- Sample_correlations_violin_COMB
  plot_list[[paste0(Status, "_NIPALS_PCA")]] <- PCA_with_NAs_ggplot_COMB

}


## Patchwork plot ----------------------------------------------------------


# Patchwork plot
QC_Imp <- plot_list[["noCombat_noReps_NIPALS_PCA"]] + plot_list[["noCombat_noReps_Sample_Cor"]] +
  plot_list[["Combat_noReps_NIPALS_PCA"]] + plot_list[["Combat_noReps_Sample_Cor"]] +
  plot_layout(ncol = 2, axes = "collect", guides = "collect")

## Modify plot aesthetics
QC_Imp_Mod <- QC_Imp & labs(title = NULL, subtitle = NULL) &
  theme(axis.title = element_text(size = Axis_title), axis.text = element_text(size = Axis_text), legend.position = "right")

## Correct aesthetics in plot "Sample correlations"
QC_Imp_Mod[[2]] <- QC_Imp_Mod[[2]] + coord_cartesian(ylim = c(-0.3, 0.75))
QC_Imp_Mod[[4]] <- QC_Imp_Mod[[4]] + coord_cartesian(ylim = c(-0.3, 0.75))

## Save
#ggsave(plot = QC_Imp_Mod, paste0("figures/BatchEffect_Correction/Batches_QC_Plots_PREDIX_pretherapy_Patchwork", ".png"), device = "png", units = "cm", width = 18, height = 13, dpi = 300)
#ggsave(plot = QC_Imp_Mod, paste0("figures/BatchEffect_Correction/Batches_QC_Plots_PREDIX_pretherapy_Patchwork", ".pdf"), device = "pdf", units = "cm", width = 18, height = 13)


# Session End ------------------------------------------------------------

