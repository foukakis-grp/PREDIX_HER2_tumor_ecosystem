# Dependencies
library(tableone)
library(tidyverse)
library(caret)
library(foreach)
library(stats)
library(lmtest)
library(data.table)

# -----------------------------
# Helpers
# -----------------------------

# Parse an "OR [LCI, UCI]" cell from ShowRegTable into numeric vector
.parse_ci <- function(ci_cell) {
  parts <- strsplit(ci_cell, " ")
  OR  <- as.numeric(sapply(parts, function(x) x[1]))
  LCI <- as.numeric(gsub("[[,]", "", sapply(parts, function(x) x[2])))
  UCI <- as.numeric(gsub("]",    "", sapply(parts, function(x) x[3])))
  c(OR = OR, LCI = LCI, UCI = UCI)
}

# Safely coerce a p-value to numeric
.as_num <- function(x) as.numeric(x)

# -----------------------------
# 1) Logistic regression (adjusted), with arm-specific fits and interaction test
# -----------------------------
Logistic_batch_adjER <- function(data, y, arm, variable, covariable) {
  results <- foreach(i = variable, .combine = rbind) %dopar% {
    cat(i, "...\n")

    # Whole cohort (biomarker + covariate)
    whole      <- glm(as.numeric(get(y)) ~ get(i) + get(covariable), family = "binomial", data = data)
    whole_ER   <- glm(as.numeric(get(y)) ~               get(covariable), family = "binomial", data = data)
    whole_ci   <- ShowRegTable(whole)[2]
    whole_lr_p <- lrtest(whole, whole_ER)$Pr[2]

    # T-DM1 arm
    experimental      <- glm(as.numeric(get(y)) ~ get(i) + get(covariable), family = "binomial",
                             data = data[data$Arm == "T-DM1", ])
    experimental_ER   <- glm(as.numeric(get(y)) ~               get(covariable), family = "binomial",
                             data = data[data$Arm == "T-DM1", ])
    experimental_ci   <- ShowRegTable(experimental)[2]
    experimental_lr_p <- lrtest(experimental, experimental_ER)$Pr[2]

    # DHP arm
    standard      <- glm(as.numeric(get(y)) ~ get(i) + get(covariable), family = "binomial",
                         data = data[data$Arm == "DHP", ])
    standard_ER   <- glm(as.numeric(get(y)) ~               get(covariable), family = "binomial",
                         data = data[data$Arm == "DHP", ])
    standard_ci   <- ShowRegTable(standard)[2]
    standard_lr_p <- lrtest(standard, standard_ER)$Pr[2]

    # Interaction test (biomarker × arm), adjusted for covariate
    interaction_1 <- glm(as.numeric(get(y)) ~ get(i) + get(arm) + get(covariable),
                         family = "binomial", data = data)
    interaction_2 <- glm(as.numeric(get(y)) ~ get(i) * get(arm) + get(covariable),
                         family = "binomial", data = data)
    # The interaction term is the new coefficient introduced in interaction_2
    interaction_coef_cell <- ShowRegTable(interaction_2)[5]
    interaction_lr_p      <- lrtest(interaction_1, interaction_2)$Pr[2]

    # Extract OR/CI for each model
    whole_vals <- .parse_ci(whole_ci)
    DHP_vals   <- .parse_ci(standard_ci)
    TDM1_vals  <- .parse_ci(experimental_ci)
    inter_vals <- .parse_ci(interaction_coef_cell)

    c(
      biomarker = i,
      whole_OR  = whole_vals["OR"],  Whole_LCI = whole_vals["LCI"],  Whole_UCI = whole_vals["UCI"],
      whole_lr_p = .as_num(whole_lr_p),

      DHP_OR    = DHP_vals["OR"],    DHP_LCI   = DHP_vals["LCI"],    DHP_UCI   = DHP_vals["UCI"],
      DHP_lr_p  = .as_num(standard_lr_p),

      TDM1_OR   = TDM1_vals["OR"],   TDM1_LCI  = TDM1_vals["LCI"],   TDM1_UCI  = TDM1_vals["UCI"],
      TDM1_lr_p = .as_num(experimental_lr_p),

      interaction_coefficient      = inter_vals["OR"],
      interaction_coefficient_LCI  = inter_vals["LCI"],
      interaction_coefficient_UCI  = inter_vals["UCI"],
      interaction_lr_p             = .as_num(interaction_lr_p)
    )
  }
  return(results)
}

# -----------------------------
# 2) Logistic regression (unadjusted), with arm-specific fits and interaction test
# -----------------------------
Logistic_batch_uni <- function(data, y, arm, variable) {
  results <- foreach(i = variable, .combine = rbind) %dopar% {
    cat(i, "...\n")

    # Whole cohort (biomarker only)
    whole      <- glm(as.numeric(get(y)) ~ get(i), family = "binomial", data = data)
    whole_null <- glm(as.numeric(get(y)) ~ 1,      family = "binomial", data = data)
    whole_ci   <- ShowRegTable(whole)[2]
    whole_lr_p <- lrtest(whole, whole_null)$Pr[2]

    # T-DM1 arm
    experimental      <- glm(as.numeric(get(y)) ~ get(i), family = "binomial",
                             data = data[data$Arm == "T-DM1", ])
    experimental_null <- glm(as.numeric(get(y)) ~ 1,      family = "binomial",
                             data = data[data$Arm == "T-DM1", ])
    experimental_ci   <- ShowRegTable(experimental)[2]
    experimental_lr_p <- lrtest(experimental, experimental_null)$Pr[2]

    # DHP arm
    standard      <- glm(as.numeric(get(y)) ~ get(i), family = "binomial",
                         data = data[data$Arm == "DHP", ])
    standard_null <- glm(as.numeric(get(y)) ~ 1,      family = "binomial",
                         data = data[data$Arm == "DHP", ])
    standard_ci   <- ShowRegTable(standard)[2]
    standard_lr_p <- lrtest(standard, standard_null)$Pr[2]

    # Interaction test (biomarker × arm), unadjusted
    interaction_1 <- glm(as.numeric(get(y)) ~ get(i) + get(arm), family = "binomial", data = data)
    interaction_2 <- glm(as.numeric(get(y)) ~ get(i) * get(arm), family = "binomial", data = data)
    interaction_coef_cell <- ShowRegTable(interaction_2)[3]  # interaction term cell
    interaction_lr_p      <- lrtest(interaction_1, interaction_2)$Pr[2]

    # Extract OR/CI
    whole_vals <- .parse_ci(whole_ci)
    DHP_vals   <- .parse_ci(standard_ci)
    TDM1_vals  <- .parse_ci(experimental_ci)
    inter_vals <- .parse_ci(interaction_coef_cell)

    c(
      biomarker = i,
      whole_OR  = whole_vals["OR"],  Whole_LCI = whole_vals["LCI"],  Whole_UCI = whole_vals["UCI"],
      whole_lr_p = .as_num(whole_lr_p),

      DHP_OR    = DHP_vals["OR"],    DHP_LCI   = DHP_vals["LCI"],    DHP_UCI   = DHP_vals["UCI"],
      DHP_lr_p  = .as_num(standard_lr_p),

      TDM1_OR   = TDM1_vals["OR"],   TDM1_LCI  = TDM1_vals["LCI"],   TDM1_UCI  = TDM1_vals["UCI"],
      TDM1_lr_p = .as_num(experimental_lr_p),

      interaction_coefficient      = inter_vals["OR"],
      interaction_coefficient_LCI  = inter_vals["LCI"],
      interaction_coefficient_UCI  = inter_vals["UCI"],
      interaction_lr_p             = .as_num(interaction_lr_p)
    )
  }
  return(results)
}

# -----------------------------
# 3) Logistic regression for subgroup variable (biomarker + Arm)
#    Returns OR/CI/p for biomarker effect in whole cohort adjusted for Arm.
# -----------------------------
Logistic_batch_for_subgroup <- function(data, y, variable) {
  results <- foreach(i = variable, .combine = rbind) %dopar% {
    cat(i, "...\n")
    fit       <- glm(as.numeric(get(y)) ~ get(i) + Arm, family = "binomial", data = data)
    fit_ci    <- ShowRegTable(fit)[2]
    fit_p     <- ShowRegTable(fit)[2, 2]
    fit_vals  <- .parse_ci(fit_ci)

    c(
      biomarker = i,
      whole_OR  = fit_vals["OR"],
      Whole_LCI = fit_vals["LCI"],
      Whole_UCI = fit_vals["UCI"],
      whole_lr_p = .as_num(fit_p)
    )
  }
  return(results)
}

# -----------------------------
# 4) Binary subgroup analysis
#    For each binary feature (e.g., mutation), fit pCR ~ Arm + ER within WT and Mut groups.
# -----------------------------
Logistic_batch_bin_subgroup <- function(data, subgroup) {
  results <- foreach(i = subgroup, .combine = rbind) %dopar% {
    cat(i, "...\n")

    # Mutant subgroup
    Mut    <- glm(as.numeric(pCR) ~ as.factor(Arm) + as.factor(ER),
                  family = "binomial", data = data %>% dplyr::filter(get(i) == "1"))
    Mut_ci <- ShowRegTable(Mut)[2]
    Mut_p  <- ShowRegTable(Mut)[2, 2]

    # Wild-type subgroup
    WT    <- glm(as.numeric(pCR) ~ as.factor(Arm) + as.factor(ER),
                 family = "binomial", data = data %>% dplyr::filter(get(i) == "0"))
    WT_ci <- ShowRegTable(WT)[2]
    WT_p  <- ShowRegTable(WT)[2, 2]

    WT_vals  <- .parse_ci(WT_ci)
    Mut_vals <- .parse_ci(Mut_ci)

    c(
      biomarker = i,
      WT_OR  = WT_vals["OR"],  WT_LCI = WT_vals["LCI"],  WT_UCI = WT_vals["UCI"],  WT_p  = .as_num(WT_p),
      Mut_OR = Mut_vals["OR"], Mut_LCI = Mut_vals["LCI"], Mut_UCI = Mut_vals["UCI"], Mut_p = .as_num(Mut_p)
    )
  }
  return(results)
}

# -----------------------------
# 5) Continuous subgroup analysis
#    For each continuous feature that has been discretized into "High"/"Low",
#    fit pCR ~ Arm + ER within each stratum and report OR/CI/p for Arm.
# -----------------------------
Logistic_batch_continuous_subgroup <- function(data, subgroup) {
  results <- foreach(i = subgroup, .combine = rbind) %dopar% {
    cat(i, "...\n")

    High    <- glm(as.numeric(pCR) ~ as.factor(Arm) + as.factor(ER),
                   family = "binomial", data = data %>% dplyr::filter(get(i) == "High"))
    Low     <- glm(as.numeric(pCR) ~ as.factor(Arm) + as.factor(ER),
                   family = "binomial", data = data %>% dplyr::filter(get(i) == "Low"))
    High_ci <- ShowRegTable(High)[2]; High_p <- ShowRegTable(High)[2, 2]
    Low_ci  <- ShowRegTable(Low)[2];  Low_p  <- ShowRegTable(Low)[2,  2]

    Low_vals  <- .parse_ci(Low_ci)
    High_vals <- .parse_ci(High_ci)

    c(
      biomarker = i,
      Low_OR  = Low_vals["OR"],  Low_LCI  = Low_vals["LCI"],  Low_UCI  = Low_vals["UCI"],  Low_p  = .as_num(Low_p),
      High_OR = High_vals["OR"], High_LCI = High_vals["LCI"], High_UCI = High_vals["UCI"], High_p = .as_num(High_p)
    )
  }
  return(results)
}

# -----------------------------
# 6) Percentile slicing utility
#    Create High/Low factors across a sequence of percentile cutoffs for given columns.
#    Example: slice_metric(df, c("geneA","geneB"), start_per = 10, end_per = 90, delta = 10)
#    will add columns like geneA_per_10, geneA_per_20, ..., each labeled "Low"/"High".
# -----------------------------
slice_metric <- function(df, sub_variable, start_per, end_per, delta) {
  stopifnot(start_per >= 0, end_per <= 100, delta > 0, start_per < end_per)

  for (col in sub_variable) {
    if (!is.numeric(df[[col]])) {
      warning(sprintf("Column '%s' is not numeric; skipping.", col))
      next
    }
    for (p in seq(from = start_per, to = end_per, by = delta)) {
      cut_off <- as.numeric(quantile(df[[col]], p / 100, na.rm = TRUE))
      new_col <- paste0(col, "_per_", p)
      df[[new_col]] <- cut(df[[col]],
                           breaks = c(-Inf, cut_off, Inf),
                           labels = c("Low", "High"),
                           include.lowest = TRUE, right = TRUE)
    }
  }
  df
}
