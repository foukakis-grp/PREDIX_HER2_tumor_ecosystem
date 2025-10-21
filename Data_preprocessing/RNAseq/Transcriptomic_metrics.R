# ============================================================
# Transcriptomic metrics pipeline
# - GGI / PIK3CA-GS (genefu)
# - ssGSEA (IOBR: KEGG sets + custom signature collection)
# - Simple mRNA panel features
# - PAM50 (SSP) labels (from Excel exports)
# - HER2DX component + composite scores (heuristic from paper)
# - Taxane response score (Mitotic − Ceramide)
# - TIDE results (precomputed)
# - Danaher cell-type scores
# - Merge -> save RDS + TSV
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(genefu)
  library(org.Hs.eg.db)
  library(IOBR)
  library(readxl)
})

# -----------------------------
# Paths (edit to your environment)
# -----------------------------
tpm_path      <- "E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds"
pam50_long_xl <- "E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_UF-3016_rnaseq.3.3_longitudinal_PAM50.xlsx"
pam50_pilot_xl<- "E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_UC-2913_rnaseq.3.3_pilot_PAM50.xlsx"
tide_path     <- "E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TIDE.txt"
danaher_gannot<- "E:/Projects/PREDIX_HER2/Multimodal/Resource/Danaher/gannot-cip-junecelltype.csv"

out_rds       <- "E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds"
out_tsv       <- "E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt"

# -----------------------------
# Load TPM (genes x samples)
# -----------------------------
tpm <- readRDS(tpm_path)
stopifnot(is.matrix(tpm) || is.data.frame(tpm))
tpm <- as.matrix(tpm)
# Ensure SYMBOL rownames
if (is.null(rownames(tpm))) stop("TPM matrix must have gene symbols as rownames.")
samples <- colnames(tpm)

# -----------------------------
# Helper: safe left_join by "sampleID"
# -----------------------------
join_by_sample <- function(x, y) left_join(x, y, by = "sampleID")

# -----------------------------
# 1) genefu: GGI and PIK3CA-GS
# -----------------------------
message("[genefu] Mapping symbols to Entrez and computing GGI / PIK3CA-GS...")
genes <- tibble(probe = rownames(tpm))
genes$EntrezGene.ID <- mapIds(org.Hs.eg.db,
                              keys = genes$probe,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first")
rownames(genes) <- genes$probe

data(sig.ggi); data(sig.pik3cags); data(sig.genius)  # (genefu ships these)
ggi_res    <- ggi(data = t(tpm), annot = genes, do.mapping = TRUE)
pik3ca_res <- pik3cags(data = t(tpm), annot = genes, do.mapping = TRUE)

genefu_df <- tibble(
  sampleID  = rownames(ggi_res$score),
  ggi_sig   = as.numeric(ggi_res$score),
  pik3ca_sig= as.numeric(pik3ca_res)
)

# -----------------------------
# 2) IOBR ssGSEA (KEGG subset + custom signatures)
#    Requires: objects `kegg` and `signature_collection` in scope
# -----------------------------
message("[IOBR] ssGSEA for KEGG + custom signature collection...")
if (!exists("kekk") && !exists("kegg")) {
  stop("Object `kegg` not found. Load KEGG signature list before running.")
}
if (!exists("signature_collection"))
  stop("Object `signature_collection` not found. Load your custom signature collection.")

input_expr <- round(tpm, 2)  # matrix genes x samples

ABC <- calculate_sig_score(
  eset = input_expr,
  signature = kegg,
  method = "ssgsea",
  parallel.size = 16
)

# Select KEGG pathways of interest and rename
ABC <- as_tibble(ABC) |>
  select(ID, KEGG_ABC_TRANSPORTERS, KEGG_APOPTOSIS, KEGG_LYSOSOME,
         KEGG_ENDOCYTOSIS, KEGG_OXIDATIVE_PHOSPHORYLATION,
         KEGG_PURINE_METABOLISM, KEGG_CITRATE_CYCLE_TCA_CYCLE,
         KEGG_GLUTATHIONE_METABOLISM, KEGG_FATTY_ACID_METABOLISM) |>
  rename(
    sampleID                 = ID,
    ABC_transporter          = KEGG_ABC_TRANSPORTERS,
    Apoptosis                = KEGG_APOPTOSIS,
    Lysosome                 = KEGG_LYSOSOME,
    Endocytosis              = KEGG_ENDOCYTOSIS,
    Oxidative_phosphorylation= KEGG_OXIDATIVE_PHOSPHORYLATION,
    Purine_metabolism        = KEGG_PURINE_METABOLISM,
    Citrate_cycles           = KEGG_CITRATE_CYCLE_TCA_CYCLE,
    Glutathione_metabolism   = KEGG_GLUTATHIONE_METABOLISM,
    Fatty_acid_metabolism    = KEGG_FATTY_ACID_METABOLISM
  )

sig <- calculate_sig_score(
  eset = input_expr,
  signature = signature_collection,
  method = "ssgsea",
  parallel.size = 16
) |>
  as_tibble() |>
  select(ID, Glycolysis, Nature_metabolism_Hypoxia, EMT2, CMLS_Review_Exosome) |>
  rename(
    sampleID = ID,
    Hypoxia  = Nature_metabolism_Hypoxia,
    EMT      = EMT2,
    Exosome  = CMLS_Review_Exosome
  )

Sig <- join_by_sample(ABC, sig)

# -----------------------------
# 3) Simple mRNA panel
# -----------------------------
panel_genes <- c("ESR1","ERBB2","CD8A","PTPRC","CD274","MKI67","PGR","MUC4")
missing_panel <- setdiff(panel_genes, rownames(tpm))
if (length(missing_panel)) warning("Panel genes missing from TPM: ", paste(missing_panel, collapse = ", "))

mRNA <- t(tpm[intersect(panel_genes, rownames(tpm)), , drop = FALSE]) |>
  as.data.frame() |>
  rownames_to_column("sampleID") |>
  as_tibble()

# prefix columns with mRNA-
mRNA <- mRNA |>
  rename_with(~ paste0("mRNA-", .x), -sampleID) |>
  mutate(`mRNA-ESR1_PGR_mean` = rowMeans(across(c(`mRNA-ESR1`, `mRNA-PGR`)), na.rm = TRUE))

# -----------------------------
# 4) PAM50 (SSP) from Excel outputs
# -----------------------------
message("[PAM50] Loading SSP labels from Excel exports...")
longitudinal <- read_excel(pam50_long_xl) |>
  select(sid, subtypecd_sspbc.subtype) |>
  rename(sampleID = sid, sspbc.subtype = subtypecd_sspbc.subtype)

pilot <- read_excel(pam50_pilot_xl) |>
  select(sid, subtypecd_sspbc.subtype.pilot) |>
  rename(sampleID = sid, sspbc.subtype = subtypecd_sspbc.subtype.pilot) |>
  mutate(sampleID = str_split_fixed(sampleID, "_", 3)[, 1])

sspbc <- bind_rows(longitudinal, pilot) |>
  mutate(
    sspbc.subtype = recode(sspbc.subtype, H2 = "Her2", LA = "LumA", LB = "LumB", BL = "Basal"),
    sspbc.subtype = factor(sspbc.subtype, levels = c("LumA","LumB","Her2","Basal"))
  ) |>
  filter(sampleID %in% samples)

# -----------------------------
# 5) HER2DX components and composite scores
#    (Heuristic replication from the paper’s component lists)
# -----------------------------
message("[HER2DX] Computing component scores...")
her2dx_genes <- list(
  IGG = c("CD27","CD79A","HLA-C","IGJ","IGKC","IGLC6","IGLV3-25","IL2RG",
          "CXCL8","LAX1","NTN3","PIM2","POU2AF1","TNFRSF17"),
  prolif = c("EXO1","ASPM","NEK2","KIF23"),
  luminal = c("BCL2","DNAJC12","AGR3","AFF3","ESR1"),
  her2amp = c("ERBB2","GRB7","STARD3","TCAP")
)

get_sum <- function(genes) {
  g <- intersect(genes, rownames(tpm))
  if (length(g) == 0) return(rep(NA_real_, ncol(tpm)))
  colSums(tpm[g, , drop = FALSE])
}

HER2DX <- tibble(
  sampleID             = samples,
  HER2DX_IGG           = get_sum(her2dx_genes$IGG),
  HER2DX_prolif        = get_sum(her2dx_genes$prolif),
  HER2DX_luminal       = get_sum(her2dx_genes$luminal),
  HER2DX_HER2_amplicon = get_sum(her2dx_genes$her2amp)
) |>
  mutate(
    HER2DX_pCR_likelihood_score =
      HER2DX_IGG/14 + HER2DX_prolif/4 + HER2DX_HER2_amplicon/5 - HER2DX_luminal/4,
    HER2DX_risk_score =
      HER2DX_IGG/14 + HER2DX_luminal/4 - HER2DX_prolif/4
  )

# -----------------------------
# 6) Taxane signature: (Mitotic − Ceramide)
# -----------------------------
message("[Taxane] Computing Mitotic − Ceramide signature...")
mitotic_genes  <- c("BUB1B","CDK1","AURKB","TTK")
ceramide_genes <- c("UGCG","COL4A3BP")

mit_m <- tpm[intersect(mitotic_genes,  rownames(tpm)), , drop = FALSE]
cer_m <- tpm[intersect(ceramide_genes, rownames(tpm)), , drop = FALSE]

Taxane <- tibble(
  sampleID = samples,
  Taxane_response = rowMeans(t(mit_m), na.rm = TRUE) - rowMeans(t(cer_m), na.rm = TRUE)
)

# -----------------------------
# 7) TIDE (precomputed) — make colnames syntactic and prefixed
# -----------------------------
message("[TIDE] Importing TIDE scores...")
TIDE <- fread(tide_path)
# Expect first column to be sample ID (V1)
if (!"V1" %in% names(TIDE)) stop("Unexpected TIDE columns; ensure first column is sample ID (V1).")
TIDE$sampleID <- TIDE$V1

tide_cols <- c("TIDE","IFNG","MSI Score","Dysfunction","Exclusion","MDSC","CAF","TAM M2","CTL")
stopifnot(all(tide_cols %in% names(TIDE)))

TIDE <- TIDE[, c("sampleID", tide_cols), with = FALSE]
names(TIDE) <- c("sampleID", paste0("TIDE_", make.names(tide_cols)))

# -----------------------------
# 8) Danaher cell-type scores
# -----------------------------
message("[Danaher] Computing cell-type scores...")
gannot <- read.csv(danaher_gannot, row.names = 1, check.names = FALSE)
celltypegenes <- rownames(gannot)[!is.na(gannot$Cell.Type.TCGA2)]
celltypes     <- setdiff(unique(gannot$Cell.Type.TCGA2), NA)

celltypematrix <- matrix(0, nrow = length(celltypegenes), ncol = length(celltypes),
                         dimnames = list(celltypegenes, celltypes))
for (ct in celltypes) {
  ct_genes <- rownames(gannot)[gannot$Cell.Type.TCGA2 == ct]
  ct_genes <- intersect(ct_genes, rownames(celltypematrix))
  celltypematrix[ct_genes, ct] <- 1
}

# Intersect with TPM rows
common_genes <- intersect(rownames(tpm), rownames(celltypematrix))
df <- t(tpm[common_genes, , drop = FALSE])

# Simple averaging scheme: multiply by normalized indicator matrix
denom <- colSums(celltypematrix[common_genes, , drop = FALSE])
denom[denom == 0] <- 1
Danaher <- as.data.frame(df[, colnames(celltypematrix), drop = FALSE] %*%
                           t(t(celltypematrix[common_genes, , drop = FALSE]) / denom)) |>
  as_tibble(rownames = "sampleID")

names(Danaher) <- paste0("Danaher-", gsub(" ", "-", names(Danaher)))
# Build TILs (mean of cell types highly correlated with CD45)
corr <- suppressWarnings(cor(Danaher |> select(-sampleID)))
cd45_col <- which(colnames(Danaher) == "Danaher-CD45")
if (length(cd45_col) == 1) {
  cells_to_use <- names(corr[, cd45_col])[corr[, cd45_col] > 0.7]
  cells_to_use <- setdiff(cells_to_use, "Danaher-CD45")
  Danaher <- Danaher |>
    mutate(`Danaher-TILs` = if (length(cells_to_use)) rowMeans(across(all_of(cells_to_use)), na.rm = TRUE) else NA_real_)
}

# -----------------------------
# 9) (Optional) hacksig
#    Requires: object `hacksig` already computed in your session
# -----------------------------
if (!exists("hacksig")) {
  warning("Object `hacksig` not found; skipping hacksig join.")
  hacksig <- tibble(sampleID = character())
}

# -----------------------------
# 10) Assemble final table
# -----------------------------
message("[Assemble] Merging all metrics...")
transcriptomic <- tibble(sampleID = samples) |>
  join_by_sample(mRNA)     |>
  join_by_sample(Sig)      |>
  join_by_sample(sspbc)    |>
  join_by_sample(Taxane)   |>
  join_by_sample(HER2DX)   |>
  join_by_sample(genefu_df)|>
  join_by_sample(hacksig)  |>
  join_by_sample(TIDE)     |>
  join_by_sample(Danaher)

# Derive patientID from sampleID (positions 9–12)
transcriptomic <- transcriptomic |>
  mutate(patientID = substr(sampleID, 9, 12)) |>
  relocate(patientID, .before = sampleID)

# Save
message("[Save] Writing: RDS + TSV")
saveRDS(transcriptomic, out_rds)
write.table(transcriptomic, file = out_tsv, quote = FALSE, row.names = FALSE, sep = "\t")

message("[Done] Rows: ", nrow(transcriptomic), "  Cols: ", ncol(transcriptomic))









