

# Setup -------------------------------------------------------------------


# Libraries
library(janitor)
library(tidyverse)


# Load functions
source("scripts/setup.R")


# Load files --------------------------------------------------------------


Prot_Quant_V3 <- readRDS(file = "processed/Protein_Quant_Predix_pretherapy_wide_V3_Norm.rds")


# Imputation --------------------------------------------------------------

## Protein rank visualization --------------------------------------------


Prot_Quant_Vis <- Prot_Quant_V3

# Remove rows containing solely NAs
Prot_Quant_Vis <- Prot_Quant_Vis %>%
  filter(!rowSums(is.na(.)) == ncol(.))

# Overlap = how often does a protein appear proportionally in all samples
Protein_overlap <- Prot_Quant_Vis %>%
  mutate(Overlap = (ncol(Prot_Quant_Vis)-rowSums(is.na(Prot_Quant_Vis)))/ncol(Prot_Quant_Vis))

# Plot: quantified protein overlap across samples

## Rank proteins based on appearance across samples
Protein_overlap_rank <- Protein_overlap %>%
  rownames_to_column(var = "Gene.name") %>%
  dplyr::select(Overlap) %>%
  arrange(desc(Overlap)) %>%
  mutate(Rank = seq(1, nrow(Protein_overlap), 1)) %>%
  mutate(Total_sample_Nr = Overlap * ncol(Prot_Quant_Vis))

## Highest protein rank across samples per overlap threshold
Max_rank_overlap <- Protein_overlap_rank %>%
  group_by(Overlap) %>%
  filter(Rank == max(Rank)) %>%
  ungroup()

## Define overlap cutoffs of interest and get number of proteins with equal or higher overlap
Cutoffs <- c(0.6)
Intervals_OI <- numeric()
Interval_labels <- as.character()

for(i in Cutoffs){ #i <- 0.75 #for test purposes
  Interval_temp <- Max_rank_overlap$Rank[max(which(Max_rank_overlap$Overlap >= i))]
  Interval_label_temp <- paste0(as.character(Interval_temp), " proteins, >= ", as.character(i*100), "% overlap")

  Intervals_OI <- append(Intervals_OI, Interval_temp)
  Interval_labels <- append(Interval_labels, Interval_label_temp)

}

## Plot ranked overlap results
Plot_prot_overlap_rank <- ggplot(Protein_overlap_rank,
                                 aes(x = Rank, y = Total_sample_Nr)) +
  geom_point(size = Point_size) +
  geom_vline(xintercept = Intervals_OI, linetype = "dashed", linewidth = Linewidth) +
  annotate("text", x = Intervals_OI, y = median(unique(Protein_overlap_rank$Total_sample_Nr))-20, label = Interval_labels, angle = 90, vjust = -1, size = Axis_title/ggplot2::.pt) +
  scale_y_continuous(breaks = c(1, seq(10, floor(max(Protein_overlap_rank$Total_sample_Nr)/10)*10, 10), max(Protein_overlap_rank$Total_sample_Nr)), limits = c(0, max(Protein_overlap_rank$Total_sample_Nr)), expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(breaks = c(1, seq(1000, ceiling(max(Protein_overlap_rank$Rank)/1000)*1000, 1000)), limits = c(min(Protein_overlap_rank$Rank), ceiling(max(Protein_overlap_rank$Rank)/1000)*1000), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = paste0("Overlap of quantified proteins across samples"),
       subtitle = paste0("Sample number = ", as.character(max(Protein_overlap_rank$Total_sample_Nr)), ", Total number of proteins = ", as.character(nrow(Prot_Quant_Vis))),
       y = "Samples containing quantified protein",
       x = "Protein rank based on quantification frequency across samples") +
  theme_bw() +
  coord_cartesian(xlim = c(1,10200)) +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.title.y = element_text(size = Axis_title, margin = margin(r = 7)),
        axis.title.x = element_text(size = Axis_title, margin = margin(t = 7)),
        axis.text = element_text(size = Axis_text))

Plot_prot_overlap_rank

#ggsave(plot = Plot_prot_overlap_rank, paste0("figures/Imputation/Rank_quantified_protein_overlap_PREDIX_pretherapy_Dotplot", ".png"), device = "png", unit = "cm", width = 14, height = 13, dpi = 300)
#ggsave(plot = Plot_prot_overlap_rank, paste0("figures/Imputation/Rank_quantified_protein_overlap_PREDIX_pretherapy_Dotplot", ".pdf"), device = "pdf", unit = "cm", width = 14, height = 13)


## Apply protein completeness cutoff ---------------------------------------


# Protein cutoff: minimally 60% coverage across samples

INPUT_quant <- Prot_Quant_V3 %>%
  filter(rowMeans(is.na(.)) <= 0.40)


## Perform Perseus Imputation ----------------------------------------------


# Define output
df <- df1 <- INPUT_quant

# Define scaled standard deviation ("width") and downshifted mean ("downshift") of random draw distribution
width <- 0.3
downshift <- 1.8

# Perform imputation per sample
for(i in 1:ncol(df1)){ #temp <- df1[[1]]
  temp <- df1[[i]]
  if(sum(is.na(temp))>0){
    temp.sd <- width * stats::sd(temp, na.rm = TRUE)
    temp.mean <- mean(temp, na.rm = TRUE) - downshift * stats::sd(temp, na.rm = TRUE)
    n.missing <- sum(is.na(temp))
    set.seed(1234) #ensures reproducibility while getting different values per column since temp.mean & temp.sd are different per iteration
    temp[is.na(temp)] <- stats::rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    df[[i]]<-temp
  }
}


Prot_Quant_Imp <- df

# Save output
#saveRDS(Prot_Quant_Imp, file = "processed/Protein_Quant_Predix_pretherapy_wide_V4_Imp.rds")


# Session End ------------------------------------------------------------

