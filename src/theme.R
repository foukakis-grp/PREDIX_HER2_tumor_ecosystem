# =============================================================================
# Publication theme & palettes for ggplot2
# =============================================================================
# Dependencies (load once here; avoid library() inside functions)
library(ggplot2)
library(ggpubr)   # if you use stat_compare_means, etc.
library(ggthemes) # for theme_foundation()
library(scales)   # for discrete_scale()
library(grid)     # for unit()

# -----------------------------------------------------------------------------
# Theme: theme_manuscript()
# A clean, journal-friendly theme built on ggthemes::theme_foundation.
# -----------------------------------------------------------------------------
theme_manuscript <- function(base_size = 12, base_family = "Arial") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.background= element_rect(fill = NA, colour = NA),
      panel.border    = element_rect(colour = NA),
      axis.title.y    = element_text(angle = 90, vjust = 2, size = base_size),
      axis.title.x    = element_text(vjust = -0.2, size = base_size),
      axis.text       = element_text(size = base_size),
      axis.line       = element_line(colour = "black"),
      axis.ticks      = element_line(),
      panel.grid.major= element_line(colour = "#f0f0f0"),
      panel.grid.minor= element_blank(),
      legend.key      = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction= "horizontal",
      legend.key.size = unit(0.5, "cm"),
      legend.spacing  = unit(0, "cm"),
      strip.background= element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text      = element_text(size = base_size)
    )
}

# -----------------------------------------------------------------------------
# Palettes (named for reproducibility)
# -----------------------------------------------------------------------------
.pal_Arm         <- c(DHP  = "#8491B4FF", `T-DM1` = "#91D1C2FF")
.pal_pCR_RD      <- c(pCR  = "#20A39E",   RD     = "#fdb462")
# PAM50 order: LumA, LumB, HER2, Basal, Normal
.pal_PAM50       <- c(LumA = "#1f78b4", LumB = "#a6cee3", HER2 = "#fb9a99",
                      Basal = "#e31a1c", Normal = "#33a02c")
.pal_ER          <- c(ERneg = "#5ac6e9", ERpos = "#e97371")
.pal_modality    <- c(Proteomics = "#8c564b", RNA = "#1f77b4", DNA = "#ff7f0e",
                      Imaging = "#2ca02c", Immune = "#d62728", Other = "#9467bd")

# -----------------------------------------------------------------------------
# Discrete fill/color scales
# (You can rename “colour”→“color” if you prefer US spelling; aliases provided.)
# -----------------------------------------------------------------------------
scale_fill_Arm <- function(...) {
  discrete_scale("fill", "Arm", manual_pal(values = .pal_Arm), ...)
}
scale_color_Arm  <- function(...) { discrete_scale("colour", "Arm",  manual_pal(values = .pal_Arm), ...) }
scale_colour_Arm <- scale_color_Arm

scale_colour_pCR_RD <- function(...) {
  discrete_scale("colour", "pCR_RD", manual_pal(values = .pal_pCR_RD), ...)
}
scale_color_pCR_RD  <- scale_colour_pCR_RD

scale_colour_PAM50 <- function(...) {
  discrete_scale("colour", "PAM50", manual_pal(values = .pal_PAM50), ...)
}
scale_color_PAM50 <- scale_colour_PAM50

scale_fill_ER <- function(...) {
  discrete_scale("fill", "ER", manual_pal(values = .pal_ER), ...)
}
scale_color_ER  <- function(...) { discrete_scale("colour", "ER", manual_pal(values = .pal_ER), ...) }
scale_colour_ER <- scale_color_ER

scale_fill_modality <- function(...) {
  discrete_scale("fill", "modality", manual_pal(values = .pal_modality), ...)
}
scale_color_modality  <- function(...) { discrete_scale("colour", "modality", manual_pal(values = .pal_modality), ...) }
scale_colour_modality <- scale_color_modality