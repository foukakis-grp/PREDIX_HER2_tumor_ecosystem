##########AUC##########
library(data.table);library(ggplot2);library(Polychrome)
# All
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/TREAT_both_rank_features_fused_proteomics_07_ROCarray.csv")
rename <- c("tpr","fpr")

names(Int) <- rename
Int[, data := "Integrated with proteomic - ROC AUC: 0.87 ± 0.07"]

All <- rbind(Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Integrated with proteomic - ROC AUC: 0.87 ± 0.07"))

ROC_all <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_all

ggsave(ROC_all, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/All_roc_curves.pdf", width = 6.2, height = 6.2)

# T-DM1
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/TREAT_1_rank_features_fused_proteomics_07_ROCarray.csv")
rename <- c("tpr","fpr")

names(Int) <- rename
Int[, data := "Integrated with proteomic - ROC AUC: 0.81 ± 0.07"]

All <- rbind(Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Integrated with proteomic - ROC AUC: 0.81 ± 0.07"))

ROC_all <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_all

ggsave(ROC_all, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/TDM1_roc_curves.pdf", width = 6.2, height = 6.2)

# DHP
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/TREAT_0_rank_features_fused_proteomics_07_ROCarray.csv")
rename <- c("tpr","fpr")

names(Int) <- rename
Int[, data := "Integrated with proteomic - ROC AUC: 0.81 ± 0.05"]

All <- rbind(Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Integrated with proteomic - ROC AUC: 0.81 ± 0.05"))

ROC_all <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_all

ggsave(ROC_all, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/DHP_roc_curves.pdf", width = 6.2, height = 6.2)

##########Frequency of selected features##########
library(data.table);library(ggplot2);library(forcats);library(tidyverse)
All=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/All_rank_features.csv")
All=All[-1,];names(All)=c("Feature","Frequency");All=All[All$Frequency>30,]
All$group=sub("_.*", "",All$Feature) 
DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/DHP_rank_features.csv")
DHP=DHP[-1,];names(DHP)=c("Feature","Frequency");DHP=DHP[DHP$Frequency>=30,]
DHP$group=sub("_.*", "",DHP$Feature) 
TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/TDM1_rank_features.csv")
TDM1=TDM1[-1,];names(TDM1)=c("Feature","Frequency");TDM1=TDM1[TDM1$Frequency>30,]
TDM1$group=sub("_.*", "",TDM1$Feature) 
group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0", 
  WSI = "#00743f",
  Prot= "#d62728"
)
#All
All$Feature=factor(All$Feature,levels =c("RNA_mRNA-ESR1","RNA_mRNA-PGR","Prot_SLC12A2","WSI_Cell_Interaction",    
           "Prot_RAB11B","Prot_RPL19","Prot_CDK12","RNA_HER2DX_HER2_amplicon","DNA_RAB11FIP1_CNA","Prot_EEA1",
           "RNA_ABC_transporter","RNA_sspbc.subtype","Prot_ERBB2_PG","Prot_ERBB2"))
All$Feature <- fct_rev(All$Feature)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
All_feature=ggplot(data=All,aes(Feature,Frequency,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Frequency (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 100)) # 设置X轴范围为0到100
All_feature
ggsave(All_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/All_feature_prot_rank.pdf", width = 8.3, height = 5)
library(IOBR);library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")
df=as.data.frame(df)
df=df[df$Prot_ERBB2!= "Unknown",]
RNA=c("RNA_mRNA-ESR1","RNA_mRNA-PGR","Prot_SLC12A2","WSI_Cell_Interaction",    
      "Prot_RAB11B","Prot_RPL19","Prot_CDK12","RNA_HER2DX_HER2_amplicon","DNA_RAB11FIP1_CNA","Prot_EEA1",
      "RNA_ABC_transporter","Prot_ERBB2","pCR")
data <- df[,RNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = RNA
)
res
#DHP
DHP$Feature=factor(DHP$Feature,levels =c("Prot_CDK12","Prot_MIEN1","Prot_RPL19","RNA_mRNA-ERBB2",             
                                         "Prot_ERBB2_PG","Prot_FCGR3A","DNA_MED1_CNA","DNA_COSMIC.Signature.10",    
                                         "WSI_Cell_Interaction","Prot_PPFIA1","DNA_TCRA.tcell.fraction.adj",
                                         "DNA_BRCA2_CNA","RNA_Purine_metabolism"))
DHP$Feature <- fct_rev(DHP$Feature)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
DHP_feature=ggplot(data=DHP,aes(Feature,Frequency,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Frequency (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 100)) 
DHP_feature
ggsave(DHP_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/DHP_feature_rank.pdf", width = 7.8, height = 5)


Protein=c("Prot_CDK12","Prot_MIEN1","Prot_RPL19",            
      "Prot_FCGR3A",    
      "Prot_PPFIA1","pCR")
RNA=c("RNA_mRNA-ERBB2","RNA_Purine_metabolism","pCR")
DNA=c("DNA_MED1_CNA","DNA_COSMIC.Signature.10","DNA_TCRA.tcell.fraction.adj",
      "DNA_BRCA2_CNA","pCR")
WSI=c("WSI_Cell_Interaction","pCR")
data <- df[df$Clin_Arm=="DHP",WSI] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = WSI
)
wilcox.test(WSI_Cell_Interaction~pCR,data)

#T-DM1
TDM1$Feature=factor(TDM1$Feature,levels =c("Clin_ER","RNA_mRNA-PGR","Prot_SLC12A2","RNA_mRNA-ESR1",                  
                                           "RNA_sspbc.subtype","DNA_RAB11FIP1_CNA","RNA_ABC_transporter",
                                           "RNA_Exosome","RNA_HER2DX_pCR_likelihood_score","DNA_CX2"))
TDM1$Feature <- fct_rev(TDM1$Feature)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
TDM1_feature=ggplot(data=TDM1,aes(Feature,Frequency,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values = group_colors) +
  labs(x="",y="Frequency (%)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "top")+
  coord_flip()+
  scale_y_continuous(limits = c(0, 100)) 
TDM1_feature
ggsave(TDM1_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/protein_subgroup/TDM1_feature_rank.pdf", width = 8, height = 5)

RNA=c("RNA_mRNA-PGR","RNA_mRNA-ESR1",                  
      "RNA_ABC_transporter",
      "RNA_Exosome","RNA_HER2DX_pCR_likelihood_score","DNA_CX2","pCR")
DNA=c("DNA_RAB11FIP1_CNA","pCR")
Protein=c("Prot_SLC12A2","pCR")
data <- df[df$Clin_Arm=="T-DM1",RNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = RNA
)
data <- df[df$Clin_Arm=="T-DM1",Protein] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = Protein
)
data <- df[df$Clin_Arm=="T-DM1",]
wilcox.test(Prot_SLC12A2~pCR,data)

chisq.test(data$pCR[data$Clin_ANYNODES!="Unknown"],
           data$Clin_ANYNODES[data$Clin_ANYNODES!="Unknown"]) 
#"RNA_chr17q12_fusion","RNA_ERBB2_fusion" "Clin_ER","Clin_ANYNODES"
