

# Setup -------------------------------------------------------------------


# Libraries
library(janitor)
library(tidyverse)


# Load functions
source("scripts/setup.R")


# Load files --------------------------------------------------------------


# Quant data
Prot_Quant_V2 <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V2_HighQuality.rds")


# Normalization ---------------------------------------------------

## Input -------------------------------------------------------------------


Matrix_Quant <- as.matrix(Prot_Quant_V2)


## Mean-centering (ratios): features -------------------------------------------------


Matrix_Quant_mean_cent <- t(apply(Matrix_Quant, 1, function(x) x/mean(x, na.rm = TRUE)))


## Median normalization: samples ----------------------------------------------------


## Step 1: Calculate Median Protein Intensity for Each Sample
Sample_medians <- apply(Matrix_Quant_mean_cent, 2, function(x) median(x, na.rm = TRUE))

## Step 2: Calculate the Mean of All Median Intensities Across Samples
Mean_sample_medians <- mean(Sample_medians, na.rm = TRUE)

## Step 3: Divide Each Sample-Specific Median by the Mean Median Across Samples
Sample_ratios <- Sample_medians / Mean_sample_medians

## Step 4: Use the Ratio to Divide All Protein Intensities in the Sample
Matrix_Quant_mean_cent_median_norm <- t(t(Matrix_Quant_mean_cent) / Sample_ratios)


## Log2-transformation -----------------------------------------------------


Matrix_Quant_mean_cent_median_norm_log2 <- log2(Matrix_Quant_mean_cent_median_norm)


## Save mean-centered, mean-normalized & log2-transformed data -------------


Prot_Quant_Norm <- Matrix_Quant_mean_cent_median_norm_log2 %>%
  as.data.frame()

#saveRDS(Prot_Quant_Norm, file = "processed/Protein_Quant_Predix_pretherapy_wide_V3_Norm.rds")


# Session End ------------------------------------------------------------

