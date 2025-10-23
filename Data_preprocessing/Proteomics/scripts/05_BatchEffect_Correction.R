

# Setup -------------------------------------------------------------------


# Libraries
library(janitor)
library(tidyverse)
library(factoextra)
library(reshape2)
library(sva)


# Load functions
source("scripts/setup.R")


# Load files --------------------------------------------------------------


# Metadata
Metadata_Input <- readRDS(file = "processed/Metadata_Predix_pretherapy_V5_HighQuality.rds")

# Quant data
Prot_Quant_Input <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V4_Imp.rds")

# Check if all measured samples are included in metadata and in the same order
cat("Same sample order in quant & meta table:", all(colnames(Prot_Quant_Input) == Metadata_Input$Sample_ID_MS))


# ComBat Correction -------------------------------------------------------


# Define batches to be corrected
batch <- Metadata_Input$Plate_MS

# Perform ComBat (without biological covariates)
Quant_ComBat_noCoV <- ComBat(dat = Prot_Quant_Input, batch = batch, par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE)

Quant_ComBat_noCoV <- Quant_ComBat_noCoV %>%
  as.data.frame()

# Save
#saveRDS(Quant_ComBat_noCoV, file = "processed/Protein_Quant_Predix_pretherapy_wide_V5_ComBat.rds")


# Session End ------------------------------------------------------------

