

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


# PREDIX Metadata ----------------------------------------------------------------

## Cleaning ----------------------------------------------------------------


# PREDIX patients (provided by Kang Wang, KI)
Metadata_Predix <- read_tsv("data/Metadata_PREDIX_KW_20240318.txt")

# Proteomics samples (from Mahshid Zarrineh, KI)
Metadata_MS <- read_xlsx("data/Patient_IDs_samples_Predix_20231128.xlsx")

# Combine metadata

## Extract sample acquisition time, biological replicate info, biopsy number
Metadata_MS <- Metadata_MS %>%
  mutate(Time_manual = case_when( #Extract sample time from patient code if not already provided
    is.na(Time) ~ str_extract(Metadata_MS$Code, "(?<=-).*"),
    TRUE ~ as.character(Time)
  )) %>%
  mutate(Biol_Status = case_when( #P10 samples are biological replicates from the same time patient
    str_detect(Runs, "P10") ~ "Duplicate",
    TRUE ~ "Original"
  )) %>%
  mutate(Biopsy = case_when( #New column for "(2)" = second biopsy from a patient at the same time point
    str_detect(Code, "\\(2\\)") ~ "Second",
    TRUE ~ "First"
  ))

## Clean sample acquisition time column
Metadata_MS <- Metadata_MS %>%
  mutate(Sample_Time = case_when( #Remove "(2)" indicator from "Time_manual" column entries
    str_detect(Time_manual, "\\(2\\)") ~ str_extract(Metadata_MS$Time_manual, ".*(?=\\(2\\))"),
    TRUE ~ as.character(Time_manual)
  ))

## Double check sample time
Metadata_MS$Time == Metadata_MS$Time_manual

## Remove redundant information
Metadata_MS <- Metadata_MS %>%
  dplyr::select(-c(Time, Time_manual))

## Insert sample time description
Metadata_MS <- Metadata_MS %>%
  mutate(Sample_Time_Descr = case_when(
    Sample_Time == "0" ~ "Baseline",
    Sample_Time == "2" ~ "During therapy",
    Sample_Time == "OP" ~ "At surgery",
    TRUE ~ "No info"
  ))

## Tag MS metadata columns
colnames(Metadata_MS) <- paste0(colnames(Metadata_MS), "_MS")

## Extract patientID (for merging)
Metadata_MS <- Metadata_MS %>%
  mutate(patientID = str_extract(Code_MS, "^\\d{4}")) #first 4 digits

## Convert columns to character
Metadata_MS$Number_MS <- as.character(Metadata_MS$Number_MS)

Metadata_Predix$patientID <- as.character(Metadata_Predix$patientID)

## Change MS sample names to be compatible with downstream tools

### Insert leading zero to MS run IDs
Metadata_MS$Sample_ID_MS <- ifelse(str_detect(Metadata_MS$Runs_MS, "P"), Metadata_MS$Runs_MS, as.character(sprintf("%03d", as.numeric(Metadata_MS$Runs_MS))))

### Add prefix "S" for "Sample"
Metadata_MS$Sample_ID_MS <- ifelse(str_detect(Metadata_MS$Sample_ID_MS, "P"), Metadata_MS$Sample_ID_MS, paste0("S", Metadata_MS$Sample_ID_MS))

## Patient info but no MS run info
setdiff(Metadata_Predix$patientID, Metadata_MS$patientID)

## MS run info but no patient info
setdiff(Metadata_MS$patientID, Metadata_Predix$patientID)

## Combine tables (removes unique patientIDs)
Metadata_Predix_longit_clean <- inner_join(Metadata_Predix, Metadata_MS, by = "patientID")

## Retain pre-therapy samples
Metadata_Predix_pretherapy_clean <- Metadata_Predix_longit_clean %>%
  filter(Sample_Time_Descr_MS == "Baseline")

## Save clean longitudinal PREDIX metadata
#saveRDS(Metadata_Predix_pretherapy_clean, file = "processed/Metadata_Predix_pretherapy_V1.rds")


## Add MS preparation batches --------------------------------------------------


Metadata_V1 <- readRDS(file = "processed/Metadata_Predix_pretherapy_V1.rds")

# Read in batch data (from Mahshid Zarrineh, KI)
Batches_PREDIX <- read_xlsx("data/Sample_batches_Predix_20231204.xlsx", sheet = 1, range = "A1:F352")

## Clean & adapt to metadata format
colnames(Batches_PREDIX) <- paste0(colnames(Batches_PREDIX), "_MS")

Batches_PREDIX$Number_MS <- as.character(Batches_PREDIX$Number_MS)

Batches_PREDIX <- Batches_PREDIX %>%
  dplyr::select(-Time_MS) %>% #time for P10 duplicates is missing
  mutate(Plate_MS = case_when(
    Plate_MS == "Additional" ~ "10",
    TRUE ~ as.character(Plate_MS)
  ))

## Combine with metadata
Metadata_V2 <- left_join(Metadata_V1, Batches_PREDIX, by = c("Runs_MS", "Number_MS", "Code_MS"))

Metadata_V2$Plate_MS <- factor(Metadata_V2$Plate_MS, levels = as.character(seq(1,10,1)))

#saveRDS(Metadata_V2, file = "processed/Metadata_Predix_pretherapy_V2_BatchInfo.rds")


# PREDIX MS identifications ------------------------------------------

## Cleaning ----------------------------------------------------------------


# Load long-format Spectronaut data (measured pre-therapy samples)
Protein_Quant_Spec_long <- vroom(file = "data/Predix_directdia_GS_ENS106_BaselineSamples_withre-inj-Report_HJ_Protein_Quant_long (Normal).tsv",
                                 delim = "\t",
                                 col_select = c(R.FileName, R.PrecursorsIdentified, R.ModifiedSequencesIdentified, R.StrippedSequencesIdentified, R.ProteinGroupsIdentified))

# Extract sample IDs
MS_Ident <- Protein_Quant_Spec_long %>%
  distinct() %>%
  mutate(Sample_ID_MS = R.FileName)

MS_Ident$Sample_ID_MS <- str_extract(MS_Ident$Sample_ID_MS, "(?<=MZ_PREDIX_).+(?=_500ng)")
MS_Ident$Sample_ID_MS <- ifelse(str_detect(MS_Ident$Sample_ID_MS, "PRNr"), paste0("S", str_extract(MS_Ident$Sample_ID_MS, "\\d+")), MS_Ident$Sample_ID_MS)

## Keep original MS file names
MS_Ident <- dplyr::rename(MS_Ident, FileName_MS = R.FileName)


### Sample overlap between PREDIX metadata and MS measurements --------


## Load PREDIX metadata table
Metadata_V2 <- readRDS(file = "processed/Metadata_Predix_pretherapy_V2_BatchInfo.rds")

## How many patients with available metadata were measured?
table(Metadata_V2$Sample_ID_MS %in% MS_Ident$Sample_ID_MS)

#31x patients with metadata were not measured

### Update metadata (removal of non-measured patients)
Non_meas_samples <- as.vector(setdiff(Metadata_V2$Sample_ID_MS, MS_Ident$Sample_ID_MS))
Non_meas_samples <- Non_meas_samples[!is.na(Non_meas_samples)]

Metadata_V3 <- Metadata_V2 %>%
  filter(!Sample_ID_MS %in% Non_meas_samples)

#### Any re-injected samples?
table(duplicated(Metadata_V3$Sample_ID_MS))

#### Save file
#saveRDS(Metadata_V3, file = "processed/Metadata_Predix_pretherapy_V3_MeasuredOnly.rds")

### How many measured samples in metadata?
table(MS_Ident$Sample_ID_MS %in% Metadata_V3$Sample_ID_MS)

setdiff(MS_Ident$Sample_ID_MS, Metadata_V3$Sample_ID_MS)

#2x measured samples without metadata: removed above by merging patient and MS metadata tables without allowing unique patientIDs


#### Add MS run information to metadata ---------------------------------


MS_Ident_comb <- dplyr::rename(MS_Ident, c(Precursors_MS = R.PrecursorsIdentified, ModPeptides_MS = R.ModifiedSequencesIdentified,
                                           Peptides_MS = R.StrippedSequencesIdentified, ProteinGroups_MS = R.ProteinGroupsIdentified))

Metadata_V4 <- Metadata_V3 %>%
  left_join(MS_Ident_comb, by = "Sample_ID_MS")

#saveRDS(Metadata_V4, file = "processed/Metadata_Predix_pretherapy_V4_MSruninfo.rds")


# PREDIX Quant Data -------------------------------------------------------

## Cleaning ----------------------------------------------------------------


# Load long-format Spectronaut data (measured pre-therapy samples)
Protein_Quant_Spec_long <- vroom(file = "data/Predix_directdia_GS_ENS106_BaselineSamples_withre-inj-Report_HJ_Protein_Quant_long (Normal).tsv",
                                 delim = "\t")

# Load clean metadata
Metadata_V4 <- readRDS(file = "processed/Metadata_Predix_pretherapy_V4_MSruninfo.rds")

# Prepare quant matrix (rows: proteins / columns: samples)

## Pick first gene name for protein groups with multiple entries (semicolon-separated)
Protein_Quant_Spec_long$PG.Genes <- str_extract(Protein_Quant_Spec_long$PG.Genes, "[^;]+")

## Select columns of interest for wide format
Protein_Quant_long_sel <- Protein_Quant_Spec_long %>%
  dplyr::select(R.FileName, PG.Genes, PG.MS1Quantity)

## Replace "NaN" values with true NA values
Protein_Quant_long_sel[Protein_Quant_long_sel == "NaN"] <- NA

## Remove missing values
Protein_Quant_long_clean <- drop_na(Protein_Quant_long_sel, everything())

## Pivot into wide format
Protein_Quant_wide <- Protein_Quant_long_clean %>%
  pivot_wider(names_from = R.FileName, values_from = PG.MS1Quantity) %>%
  column_to_rownames(var = "PG.Genes")

## Correct column names
colnames(Protein_Quant_wide) <- str_extract(colnames(Protein_Quant_wide), "(?<=MZ_PREDIX_).+(?=_500ng)")
colnames(Protein_Quant_wide) <- ifelse(str_detect(colnames(Protein_Quant_wide), "PRNr"), paste0("S", str_extract(colnames(Protein_Quant_wide), "\\d+")), colnames(Protein_Quant_wide))

# Check quant file

## Check if there are semicolons left in "PG.Genes" column
sum(str_detect(rownames(Protein_Quant_wide), ";"), na.rm = TRUE)

## Filter samples according to metadata
Prot_Quant_wide_clean <- Protein_Quant_wide[, Metadata_V4$Sample_ID_MS]

all(colnames(Prot_Quant_wide_clean) == Metadata_V4$Sample_ID_MS)

## Remove protein rows without measurements (= only NA values)
Prot_Quant_wide_clean <- Prot_Quant_wide_clean %>%
  filter(rowSums(!is.na(.)) > 0) #remove proteins with only NA values

# Save clean wide protein quantification table
#saveRDS(Prot_Quant_wide_clean, file = "processed/Protein_Quant_Predix_pretherapy_wide_V1.rds")


# Session End ------------------------------------------------------------

