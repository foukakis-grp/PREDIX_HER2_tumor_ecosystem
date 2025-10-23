

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
library(limma)


# Load functions
source("scripts/setup.R")


# Differential abundance analysis ---------------------------------------------------------

## Load files --------------------------------------------------------------


# Metadata (137x high-quality samples)
Metadata <- readRDS(file = "processed/Metadata_Predix_pretherapy_V6_NoReps.rds")

## Rename "T-DM1" arm for Limma below
Metadata <- Metadata %>%
  mutate(Arm = case_when(
    Arm == "T-DM1" ~ "TDM1",
    TRUE ~ Arm
  ))

# Quant data (137x high-quality samples, no ComBat correction -> batches included as covariates in Limma model)
Prot_Quant <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V6_1_NoComBat_NoReps_MedCent.rds")


## Limma analysis ------------------------------------------------------------

### Define input ------------------------------------------------------------


# Sample sub-setting (as a list for loop later)

## Response stratified by treatment arm
DHP_meta <- Metadata %>% filter(Arm == "DHP")
DHP_quant <- Prot_Quant %>% dplyr::select(DHP_meta$Sample_ID_MS)

TDM1_meta <- Metadata %>% filter(Arm == "TDM1")
TDM1_quant <- Prot_Quant %>% dplyr::select(TDM1_meta$Sample_ID_MS)

# Contrasts of interest
Contrast_list <- list(

  # All samples
  list(name = "ALL_pCR_RD_Covar_Treat",
       meta_data = Metadata,
       quant_data = Prot_Quant,
       biol_group = "Response",
       biol_covariate_1 = "Arm",
       contrast = "pCR-RD"),
  list(name = "ALL_pCR_RD_Covar_Treat_ER",
       meta_data = Metadata,
       quant_data = Prot_Quant,
       biol_group = "Response",
       biol_covariate_1 = "Arm",
       biol_covariate_2 = "ER",
       contrast = "pCR-RD"),

  # Samples stratified by treatment arm
  list(name = "DHP_pCR_RD",
       meta_data = DHP_meta,
       quant_data = DHP_quant,
       biol_group = "Response",
       contrast = "pCR-RD"),
  list(name = "TDM1_pCR_RD",
       meta_data = TDM1_meta,
       quant_data = TDM1_quant,
       biol_group = "Response",
       contrast = "pCR-RD"),

  # Samples stratified by treatment arm (adjusted for ER status)
  list(name = "DHP_pCR_RD_Covar_ER",
       meta_data = DHP_meta,
       quant_data = DHP_quant,
       biol_group = "Response",
       biol_covariate_1 = "ER",
       contrast = "pCR-RD"),
  list(name = "TDM1_pCR_RD_Covar_ER",
       meta_data = TDM1_meta,
       quant_data = TDM1_quant,
       biol_group = "Response",
       biol_covariate_1 = "ER",
       contrast = "pCR-RD")

)


### Calculations ------------------------------------------------------------


# Loop over contrasts defined above
for(contrast in Contrast_list){ #contrast <- Contrast_list[[3]]

  Status <- contrast$name

  cat("In progress:", Status, "\n")

  if(sum(str_detect(names(contrast), "biol_")) == 1){ #one biological variable

    Contrast_OI <- as.factor(pull(contrast$meta_data, contrast$biol_group))
    Batch <- as.factor(pull(contrast$meta_data, "Plate_MS"))

    ## Define model matrix (means model = no intersect term)
    design = model.matrix(~0 + Contrast_OI + Batch)

    colnames(design) <- ifelse(str_detect(colnames(design), "Contrast_OI"), str_extract(colnames(design), "(?<=Contrast_OI).+"), colnames(design))

  } else if(sum(str_detect(names(contrast), "biol_")) == 2){ #two biological variables

    Contrast_OI <- as.factor(pull(contrast$meta_data, contrast$biol_group))
    Covariate_1 <- as.factor(pull(contrast$meta_data, contrast$biol_covariate_1))
    Batch <- as.factor(pull(contrast$meta_data, "Plate_MS"))

    ## Define model matrix (means model = no intersect term)
    design = model.matrix(~0 + Contrast_OI + Covariate_1 + Batch)

    colnames(design) <- ifelse(str_detect(colnames(design), "Contrast_OI"), str_extract(colnames(design), "(?<=Contrast_OI).+"), colnames(design))
    colnames(design) <- ifelse(str_detect(colnames(design), "Covariate_1"), str_extract(colnames(design), "(?<=Covariate_1).+"), colnames(design))

  } else if(sum(str_detect(names(contrast), "biol_")) == 3){ #three biological variables

    Contrast_OI <- as.factor(pull(contrast$meta_data, contrast$biol_group))
    Covariate_1 <- as.factor(pull(contrast$meta_data, contrast$biol_covariate_1))
    Covariate_2 <- as.factor(pull(contrast$meta_data, contrast$biol_covariate_2))
    Batch <- as.factor(pull(contrast$meta_data, "Plate_MS"))

    ## Define model matrix (means model = no intersect term)
    design = model.matrix(~0 + Contrast_OI + Covariate_1 + Covariate_2 + Batch)

    colnames(design) <- ifelse(str_detect(colnames(design), "Contrast_OI"), str_extract(colnames(design), "(?<=Contrast_OI).+"), colnames(design))
    colnames(design) <- ifelse(str_detect(colnames(design), "Covariate_1"), str_extract(colnames(design), "(?<=Covariate_1).+"), colnames(design))
    colnames(design) <- ifelse(str_detect(colnames(design), "Covariate_2"), str_extract(colnames(design), "(?<=Covariate_2).+"), colnames(design))

  }

  ## Define contrast matrix for all contrasts of interest
  ALL_contrasts <- makeContrasts(contrasts = contrast$contrast, levels = colnames(design))

  ## Limma procedure
  fit1 <- lmFit(contrast$quant_data, design = design)

  fit2 <- contrasts.fit(fit1, contrasts = ALL_contrasts)

  fit3 <- eBayes(fit2)

  Limma_res_temp <- topTable(fit3, coef = contrast$contrast, number = Inf, adjust.method = "BH", sort.by = "P")

  Limma_res_temp$gene <- rownames(Limma_res_temp)

  # Save Limma output
  write_tsv(Limma_res_temp, paste0("processed/Limma_output/", "Limma_output_", Status, ".tsv"))

}


# Agreement plot: RNA & Protein -------------------------------------------

## Load files --------------------------------------------------------------


# RNAseq-based transcriptomics data (from Kang Wang, KI)

RNAseq_DEG_path <- "data/DEG_All_RNAseq_KW_240415.xlsx"

# List of RNAseq data (used as input to load & combine with proteomics data)

RNA_contrast_list <- list(

  # RNA DEG results
  list(name = "ALL_pCR_RD_Covar_Treat",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_total_adjArm"),
  list(name = "DHP_pCR_RD",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_DHP"),
  list(name = "TDM1_pCR_RD",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_TDM1"),

  # ER-adjusted RNA DEG results
  list(name = "ALL_pCR_RD_Covar_Treat_ER",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_total_adjArm_ER"),
  list(name = "DHP_pCR_RD_Covar_ER",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_DHP_adjER"),
  list(name = "TDM1_pCR_RD_Covar_ER",
       Excel_file = RNAseq_DEG_path, Excel_sheet = "DEG_TDM1_adjER")

)


## Plots --------------------------------------------------------------------


#Two agreement plot versions: with and without pvalues highlighted

# LOOP over contrasts of interest (defined above)

for(RNA_contrast in RNA_contrast_list){ #RNA_contrast <- RNA_contrast_list[[2]]

  Status <- RNA_contrast$name

  cat("In progress:", Status, "\n")

  # RNAseq-based transcriptomics data
  RNAseq_DEG_temp <- read_xlsx(RNA_contrast$Excel_file, sheet = RNA_contrast$Excel_sheet, range = "A1:F15731")
  colnames(RNAseq_DEG_temp) <- paste0(colnames(RNAseq_DEG_temp), ".RNA")

  # MS-based proteomics data
  Prot_DEG_temp <- read_tsv(paste0("processed/Limma_output/", "Limma_output_", Status, ".tsv"))
  colnames(Prot_DEG_temp) <- paste0(colnames(Prot_DEG_temp), ".Prot")

  # Combine grouped DEG results in data frame (some measured proteins are not found in RNAseq data & vice versa)
  RNA_Prot_DEG_temp <- full_join(RNAseq_DEG_temp, Prot_DEG_temp, by = c("gene.RNA" = "gene.Prot")) %>%
    dplyr::rename("Gene_name" = "gene.RNA")

  # Define logFC cutoff for agreement plots (aligned to Kang Wang's cutoff)
  logFC_up <- 0.5
  logFC_down <- -0.5

  # p-value threshold
  pvalue <- 0.05 #to be highlighted
  pvalue_label <- paste0(pvalue*100, "p")

  # Unique plot aesthetics
  Label_size <- 2

  Point_size <- 1.5

  # Remove missing logFC values
  RNA_Prot_DEG_temp_woNA <- RNA_Prot_DEG_temp %>%
    filter(!is.na(log2FoldChange.RNA)) %>%
    filter(!is.na(logFC.Prot))

  # Create hierarchical categorical variables for RNA-Protein DEG comparison
  RNA_Prot_DEG_temp_woNA <- RNA_Prot_DEG_temp_woNA %>%
    mutate(

      # Version 1: Using both fold change thresholds AND p-value significance
      # (RNA: adjusted p-values, Protein: nominal p-values)
      category_sig = case_when(

        # Concordant: Both significant and same direction
        padj.RNA < pvalue & P.Value.Prot < pvalue &
          ((log2FoldChange.RNA > logFC_up & logFC.Prot > logFC_up) |
             (log2FoldChange.RNA < logFC_down & logFC.Prot < logFC_down)) ~ "Concordant",

        # Opposing: Both significant but opposite directions
        padj.RNA < pvalue & P.Value.Prot < pvalue &
          ((log2FoldChange.RNA > logFC_up & logFC.Prot < logFC_down) |
             (log2FoldChange.RNA < logFC_down & logFC.Prot > logFC_up)) ~ "Opposing",

        # RNA only: RNA significant, protein not (or not meeting FC threshold)
        padj.RNA < pvalue & abs(log2FoldChange.RNA) > logFC_up ~ "RNA_only",

        # Protein only: Protein significant, RNA not (or not meeting FC threshold)
        P.Value.Prot < pvalue & abs(logFC.Prot) > logFC_up ~ "Protein_only",

        # Not significant: Everything else
        TRUE ~ "Not_significant"
      ),

      # Version 2: Using only fold change thresholds (no p-values)
      category_fc_only = case_when(

        # Concordant: Both exceed FC threshold in same direction
          ((log2FoldChange.RNA > logFC_up & logFC.Prot > logFC_up) |
             (log2FoldChange.RNA < logFC_down & logFC.Prot < logFC_down)) ~ "Concordant",

        # Opposing: Both exceed FC threshold but opposite directions
          ((log2FoldChange.RNA > logFC_up & logFC.Prot < logFC_down) |
             (log2FoldChange.RNA < logFC_down & logFC.Prot > logFC_up)) ~ "Opposing",

        # RNA only: Only RNA exceeds FC threshold
        abs(log2FoldChange.RNA) > logFC_up ~ "RNA_only",

        # Protein only: Only protein exceeds FC threshold
        abs(logFC.Prot) > logFC_up ~ "Protein_only",

        # Not significant: Neither exceeds FC threshold
        TRUE ~ "Not_significant"
      )
    )


  ## Convert to factors with specific level ordering for plotting
  RNA_Prot_DEG_temp_woNA <- RNA_Prot_DEG_temp_woNA %>%
    mutate(
      category_sig = factor(category_sig,
                            levels = c("Not_significant", "RNA_only", "Protein_only",
                                       "Concordant", "Opposing")),
      category_fc_only = factor(category_fc_only,
                                levels = c("Not_significant", "RNA_only", "Protein_only",
                                           "Concordant", "Opposing"))
    )


  # Plot: considers both logFC & pvalue

  ## Arrange to enable hierarchy within the plot layers
  plot_data_sig <- RNA_Prot_DEG_temp_woNA %>%
    arrange(category_sig)

  ## Priority-based labeling
  label_data <- RNA_Prot_DEG_temp_woNA %>%
    filter(category_sig %in% c("Concordant", "Opposing"))

  ## Plot
  Comp_RNA_Prot_Scatter_sig <- ggplot(plot_data_sig, aes(x = log2FoldChange.RNA, y = logFC.Prot)) +

    geom_vline(xintercept = c(logFC_down, logFC_up), color = "black", linetype = "dashed", linewidth = Linewidth) +
    geom_hline(yintercept = c(logFC_down, logFC_up), color = "black", linetype = "dashed", linewidth = Linewidth) +

    geom_point(aes(fill = category_sig), shape = 21, size = Point_size,
               alpha = 1, stroke = Borderwidth, color = "white") +
    scale_fill_manual(values = c("Not_significant" = "grey",
                                 "RNA_only" = "gold3",
                                 "Protein_only" = "blue",
                                 "Concordant" = "green3",
                                 "Opposing" = "black")) +

    geom_label_repel(aes(label = `Gene_name`), data = label_data,
                     color = "black", size = Label_size,
                     max.overlaps = 100, #high value since labels were pre-filtered
                     seed = 1234, force = 1, max.iter = 10000) +

    # geom_smooth(method ='lm', formula = y~x, se = FALSE, color = "slategray4", linetype = "dashed", linewidth = Linewidth) +

    # stat_cor(aes(x = log2FoldChange.RNA, y = logFC.Prot), method = "spearman", label.x.npc = 0.8, label.y.npc = 0.1, inherit.aes = FALSE,
    #          color = "slategray4", size = Axis_title/ggplot2::.pt) +

    scale_y_continuous(breaks = scales::breaks_width(1)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +

    labs(title = paste0("DEG analysis: transcriptomics vs proteomics fold changes: ", Status),
         subtitle = paste0("Highlighted genes: logFC below/above ", round(logFC_up, 2), " & adj (RNA) or unadj (protein) pvalue < ", pvalue, " (gold: RNA / blue: protein / green: both / black: opposing) \n",
                           "All overlapping genes / ", "n = ", nrow(plot_data_sig), " / grey: Spearman correlation coefficient"),
         x = "log2FC RNA",
         y = "log2FC Proteins") +

    theme_classic() +
    theme(axis.text = element_text(size = Axis_text),
          axis.title = element_text(size = Axis_title),

          plot.title = element_blank(),
          plot.subtitle = element_blank(),

          legend.position = "none",
          legend.title = element_text(size = Axis_title),
          legend.text = element_text(size = Axis_title),

          panel.border = element_rect(fill = NA, colour = "black", linewidth = Linewidth),
          axis.line = element_blank())

  #Comp_RNA_Prot_Scatter_sig

  #Save (5x5 inches & smaller)
  # ggsave(plot = Comp_RNA_Prot_Scatter_sig, paste0("figures/Response_Groups/Agreement_plot/", "DEG_RNA_Prot_", Status, "_sig", pvalue_label, "_PREDIX_pretherapy_Scatterplot", ".png"),
  #        device = "png", unit = "in", width = 5, height = 5, dpi = 300)
  #
  # ggsave(plot = Comp_RNA_Prot_Scatter_sig, paste0("figures/Response_Groups/Agreement_plot/", "DEG_RNA_Prot_", Status, "_sig", pvalue_label, "_PREDIX_pretherapy_Scatterplot", ".pdf"),
  #        device = "pdf", unit = "in", width = 5, height = 5)

  # Plot: considers logFC only

  ## Arrange to enable hierarchy within the plot layers
  plot_data_fc <- RNA_Prot_DEG_temp_woNA %>%
    arrange(category_fc_only)

  ## Priority-based labeling
  label_data <- RNA_Prot_DEG_temp_woNA %>%
    filter(category_fc_only %in% c("Concordant", "Opposing"))

  ## Plot
  Comp_RNA_Prot_Scatter_logFC <- ggplot(plot_data_fc, aes(x = log2FoldChange.RNA, y = logFC.Prot)) +

    geom_vline(xintercept = c(logFC_down, logFC_up), color = "black", linetype = "dashed", linewidth = Linewidth) +
    geom_hline(yintercept = c(logFC_down, logFC_up), color = "black", linetype = "dashed", linewidth = Linewidth) +

    geom_point(aes(fill = category_fc_only), shape = 21, size = Point_size,
               alpha = 1, stroke = Borderwidth, color = "white") +
    scale_fill_manual(values = c("Not_significant" = "grey",
                                 "RNA_only" = "gold3",
                                 "Protein_only" = "blue",
                                 "Concordant" = "green3",
                                 "Opposing" = "black")) +

    geom_label_repel(aes(label = `Gene_name`), data = label_data,
                     color = "black", size = Label_size,
                     max.overlaps = 100, #high value since labels were pre-filtered
                     seed = 1234, force = 1, max.iter = 10000) +

    # geom_smooth(method ='lm', formula = y~x, se = FALSE, color = "slategray4", linetype = "dashed", linewidth = Linewidth) +

    # stat_cor(aes(x = log2FoldChange.RNA, y = logFC.Prot), method = "spearman", label.x.npc = 0.8, label.y.npc = 0.1, inherit.aes = FALSE,
    #          color = "slategray4", size = Axis_title/ggplot2::.pt) +

    scale_y_continuous(breaks = scales::breaks_width(1)) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +

    labs(title = paste0("DEG analysis: transcriptomics vs proteomics fold changes: ", Status),
         subtitle = paste0("Highlighted genes: logFC below/above ", round(logFC_up, 2), " (gold: RNA / blue: protein / green: both / black: opposing) \n",
                           "All overlapping genes / ", "n = ", nrow(plot_data_fc), " / grey: Spearman correlation coefficient"),
         x = "log2FC RNA",
         y = "log2FC Proteins") +

    theme_classic() +
    theme(axis.text = element_text(size = Axis_text),
          axis.title = element_text(size = Axis_title),

          plot.title = element_blank(),
          plot.subtitle = element_blank(),

          legend.position = "none",
          legend.title = element_text(size = Axis_title),
          legend.text = element_text(size = Axis_title),

          panel.border = element_rect(fill = NA, colour = "black", linewidth = Linewidth),
          axis.line = element_blank())

  #Comp_RNA_Prot_Scatter_logFC

  #Save (5x5 inches & smaller)
  # ggsave(plot = Comp_RNA_Prot_Scatter_logFC, paste0("figures/Response_Groups/Agreement_plot/", "DEG_RNA_Prot_", Status, "_logFC_PREDIX_pretherapy_Scatterplot", ".png"),
  #        device = "png", unit = "in", width = 5, height = 5, dpi = 300)
  #
  # ggsave(plot = Comp_RNA_Prot_Scatter_logFC, paste0("figures/Response_Groups/Agreement_plot/", "DEG_RNA_Prot_", Status, "_logFC_PREDIX_pretherapy_Scatterplot", ".pdf"),
  #        device = "pdf", unit = "in", width = 5, height = 5)

}


# Session End ------------------------------------------------------------

