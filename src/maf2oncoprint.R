# Dependencies
library(tidyverse)
library(data.table)

# -----------------------------------------------------------------------------
# maf2oncoprint
# Build an oncoprint-ready matrix (genes x samples) of semicolon-tagged events
# from a MAF-like data.frame/data.table. Each cell is a concatenation of event
# labels (e.g., "Missense;Truncating;Amplification;"), or a single space " "
# if no events were found for that gene in that sample.
#
# Expected columns in `maf`:
#   - Tumor_Sample_Barcode
#   - Hugo_Symbol
#   - Variant_Classification
#
# Event collapsing:
#   * Missense_Mutation                     -> "Missense"
#   * Frame_Shift_Del/Ins, Splice_Site,
#     Nonsense_Mutation, Nonstop_Mutation,
#     Translation_Start_Site                -> "Truncating"
#   * In_Frame_Del/Ins                      -> "In.frame"
#   * Copy-number classes (Amplification, Deletion) kept as-is
#   * Silent/Non-coding categories collapsed into "Silent" and included,
#     but you can set `include_silent = FALSE` to drop them.
#
# Returns:
#   data.frame with rownames = genes and columns = samples
# -----------------------------------------------------------------------------
maf2oncoprint <- function(maf, include_silent = TRUE) {
  # Validate input
  required_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  stopifnot(all(required_cols %in% colnames(maf)))

  # Work with data.table for speed
  maf <- as.data.table(maf)

  # Define classes to drop into "Silent"
  silent_classes <- c(
    "3'UTR","5'UTR","3'Flank","Targeted_Region","Silent","Intron","RNA","IGR",
    "Splice_Region","5'Flank","lincRNA","De_novo_Start_InFrame","De_novo_Start_OutOfFrame",
    "Start_Codon_Ins","Start_Codon_SNP","Stop_Codon_Del"
  )

  # Canonicalize Variant_Classification into a small set of display labels
  maf[, Variant_Classification := fcase(
    Variant_Classification %in% "Missense_Mutation",                               "Missense",
    Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site",
                                  "Nonsense_Mutation","Nonstop_Mutation",
                                  "Translation_Start_Site"),                      "Truncating",
    Variant_Classification %in% c("In_Frame_Del","In_Frame_Ins"),                  "In.frame",
    Variant_Classification %in% silent_classes,                                     "Silent",
    Variant_Classification %in% c("Amplification","Deletion"),                      Variant_Classification,
    default = NA_character_
  )]

  # Optionally drop "Silent"
  if (!include_silent) {
    maf <- maf[Variant_Classification != "Silent"]
  }

  # Keep only rows with a recognized class
  maf <- maf[!is.na(Variant_Classification)]

  # Cache unique gene and sample orders for consistent shaping
  gene_names <- unique(maf$Hugo_Symbol)
  samples    <- unique(maf$Tumor_Sample_Barcode)

  # Helper: build a wide sample x gene matrix for one class label (e.g., "Missense")
  # Cells are either "Label;" if >=1 event, or " " (single space) if none.
  build_class_matrix <- function(dt, class_label) {
    sub <- dt[Variant_Classification == class_label]
    if (nrow(sub) == 0L) {
      # No events of this class: return a full " " matrix
      out <- matrix(" ", nrow = length(samples), ncol = length(gene_names),
                    dimnames = list(samples, gene_names))
      return(as.data.frame(out, stringsAsFactors = FALSE))
    }

    # Count events per sample-gene, then cast
    sub_counts <- sub[, .N, by = .(Tumor_Sample_Barcode, Hugo_Symbol)]
    wide <- dcast(
      sub_counts,
      Tumor_Sample_Barcode ~ Hugo_Symbol,
      value.var = "N",
      fill = 0
    )

    # Ensure all genes/samples present; fill missing with 0
    if (!all(gene_names %in% colnames(wide))) {
      missing_genes <- setdiff(gene_names, colnames(wide))
      wide[, (missing_genes) := 0]
    }
    if (!all(samples %in% wide$Tumor_Sample_Barcode)) {
      missing_samples <- setdiff(samples, wide$Tumor_Sample_Barcode)
      add <- data.table(Tumor_Sample_Barcode = missing_samples)
      for (g in gene_names) add[[g]] <- 0L
      wide <- rbindlist(list(wide, add), use.names = TRUE, fill = TRUE)
    }

    # Order rows/cols and convert counts to tagged strings
    setDF(wide)
    rownames(wide) <- wide$Tumor_Sample_Barcode
    wide <- wide[samples, gene_names, drop = FALSE]
    wide[] <- ifelse(wide[] > 0, paste0(class_label, ";"), " ")
    as.data.frame(wide, stringsAsFactors = FALSE)
  }

  # Classes to include in the oncoprint (order matters for concatenation)
  classes <- c("Missense", "Truncating", "In.frame", "Silent", "Amplification", "Deletion")

  # Build matrices per class and combine by string concatenation
  class_mats <- lapply(classes, function(cl) build_class_matrix(maf, cl))

  # Combine class layers: paste strings per cell in the defined order
  # Start with a matrix of blanks, then paste each layer
  combined <- matrix(" ", nrow = length(samples), ncol = length(gene_names),
                     dimnames = list(samples, gene_names))
  for (m in class_mats) {
    combined <- paste0(combined, as.matrix(m))
  }

  # Remove any accidental double spaces; keep a single space for empty
  combined[combined == paste0(rep(" ", 1), collapse = "")] <- " "

  # Return genes x samples (match your original output orientation)
  out <- t(combined)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  return(out)
}
