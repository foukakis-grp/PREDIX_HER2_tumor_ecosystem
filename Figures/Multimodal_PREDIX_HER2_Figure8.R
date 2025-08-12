##########predict switching treatment##########
library(data.table);library(ggplot2);library(ggpubr);library(tidyverse)
data1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/Switching_treat/soft_probs_trained_treatment__1__tested__0.csv")
data2=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/Switching_treat/soft_probs_trained_treatment__0__tested__1.csv")
data=rbind(data1,data2);colnames(data)[3]="switch_trt_response_prob"
data=data[,c("patientID","switch_trt_response_prob")]
data$switch_trt_response="No response"
data$switch_trt_response[data$switch_trt_response_prob>0.5]="Response"
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
data=left_join(data,clin)
DHP=data[data$Arm=="DHP",]
TDM1=data[data$Arm=="T-DM1",]

mytable=DHP %>%
  group_by(Response,switch_trt_response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
plot_DHP=ggbarplot(mytable, "Response", "freq",
          fill = "switch_trt_response", color = "switch_trt_response", 
          palette = c("#fdb462","#00A087FF"))+coord_flip()

mytable=TDM1 %>%
  group_by(Response,switch_trt_response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
plot_TDM1=ggbarplot(mytable, "Response", "freq",
          fill = "switch_trt_response", color = "switch_trt_response", 
          palette = c("#fdb462","#00A087FF"))+coord_flip()

ggarrange(plot_DHP, plot_TDM1, 
          ncol = 1, nrow = 2, 
          common.legend = TRUE, legend = "right")
##########All##########
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TREAT_BOTH_FUSED/TREAT_BOTH_rank_features_fusion_07_ROCarray.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TREAT_BOTH_clinical/TREAT_BOTH_rank_features_clinical_06_ROCarray.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TREAT_BOTH_DNA/TREAT_BOTH_rank_features_DNA_06_ROCarray.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TREAT_BOTH_RNA/TREAT_BOTH_rank_features_RNA_07_ROCarray.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TREAT_BOTH_image/TREAT_BOTH_rank_features_image_06_ROCarray.csv")

rename <- c("tpr","fpr")

names(Int) <- rename
names(Clin) <- rename
names(DNA) <- rename
names(RNA) <- rename
names(WSI) <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.64 ± 0.09"]
DNA[, data := "Genomic Only - ROC AUC: 0.65 ± 0.05"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.77 ± 0.09"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.65 ± 0.09"]
Int[, data := "Integrated - ROC AUC: 0.81 ± 0.06"]


All <- rbind(Clin,DNA,RNA,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.64 ± 0.09","Genomic Only - ROC AUC: 0.65 ± 0.05",
                                    "Transcriptomic Only - ROC AUC: 0.77 ± 0.09","Digital Pathology Only - ROC AUC: 0.65 ± 0.09",
                                    "Integrated - ROC AUC: 0.81 ± 0.06"))

ROC_all <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#00743f","#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_all

ggsave(ROC_all, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/All_roc_curves.pdf", width = 6.2, height = 6.2)


##########T-DM1##########
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_ROCarray.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_clinical/TREAT_1_rank_features_clinical_06_ROCarray.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_DNA/TREAT_1_rank_features_DNA_06_ROCarray.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_RNA/TREAT_1_rank_features_RNA_07_ROCarray.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_image/TREAT_1_rank_features_image_06_ROCarray.csv")

rename <- c("tpr","fpr")

names(Int) <- rename
names(Clin) <- rename
names(DNA) <- rename
names(RNA) <- rename
names(WSI) <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.64 ± 0.08"]
DNA[, data := "Genomic Only - ROC AUC: 0.65 ± 0.07"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.74 ± 0.07"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.65 ± 0.09"]
Int[, data := "Integrated - ROC AUC: 0.75 ± 0.06"]


All <- rbind(Clin,DNA,RNA,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.64 ± 0.08","Genomic Only - ROC AUC: 0.65 ± 0.07",
     "Transcriptomic Only - ROC AUC: 0.74 ± 0.07","Digital Pathology Only - ROC AUC: 0.65 ± 0.09",
     "Integrated - ROC AUC: 0.75 ± 0.06"))

ROC_TDM1 <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#00743f","#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_TDM1

ggsave(ROC_TDM1, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/TDM1_roc_curves.pdf", width = 6.2, height = 6.2)
##########DHP##########
library(data.table);library(ggplot2);library(Polychrome)
Int=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_ROCarray.csv")
Clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_clinical/TREAT_0_rank_features_clinical_06_ROCarray.csv")
DNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_DNA/TREAT_0_rank_features_DNA_065_ROCarray.csv")
RNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_RNA/TREAT_0_rank_features_RNA_075_ROCarray.csv")
WSI=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_image/TREAT_0_rank_features_image_06_ROCarray.csv")

rename <- c("tpr","fpr")

names(Int) <- rename
names(Clin) <- rename
names(DNA) <- rename
names(RNA) <- rename
names(WSI) <- rename

Clin[, data := "Clinic Only - ROC AUC: 0.64 ± 0.09"]
DNA[, data := "Genomic Only - ROC AUC: 0.68 ± 0.06"]
RNA[, data := "Transcriptomic Only - ROC AUC: 0.77 ± 0.06"]
WSI[, data := "Digital Pathology Only - ROC AUC: 0.66 ± 0.1"]
Int[, data := "Integrated - ROC AUC: 0.79 ± 0.07"]


All <- rbind(Clin,DNA,RNA,WSI,Int)
All=All[!is.na(All$tpr),]
All$data=factor(All$data,levels = c("Clinic Only - ROC AUC: 0.64 ± 0.09","Genomic Only - ROC AUC: 0.68 ± 0.06",
                                    "Transcriptomic Only - ROC AUC: 0.77 ± 0.06","Digital Pathology Only - ROC AUC: 0.66 ± 0.1",
                                    "Integrated - ROC AUC: 0.79 ± 0.07"))

ROC_DHP <- ggplot(All, aes(y = tpr, x = fpr, colour = data, group = data)) + 
  geom_line() +
  theme_bw() + 
  geom_abline(linetype = "dashed") +
  ylab("True Positive Rate") +
  xlab("False Positive Rate") +
  labs(col = "Feature set") +
  scale_colour_manual(values = c("#434279", "#f2a104","#72a2c0","#00743f","#d62728")) +# protein "#9467bd"
  theme(legend.position = c(0.725, 0.18))+
  theme(text = element_text(family = "ArialMT"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12),  # 修改 x 轴标签字符大小
        axis.text.y = element_text(size = 12))
ROC_DHP 
ggsave(ROC_DHP , file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/DHP_roc_curves.pdf", width = 6.2, height = 6.2)

##########Frequency of selected features##########
library(data.table);library(ggplot2);library(forcats);library(tidyverse)
All=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/All_rank_features.csv")
All=All[-1,];names(All)=c("Feature","Frequency");All=All[All$Frequency>30,]
All$group=sub("_.*", "",All$Feature) 
DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/DHP_rank_features.csv")
DHP=DHP[-1,];names(DHP)=c("Feature","Frequency");DHP=DHP[DHP$Frequency>=30,]
DHP$group=sub("_.*", "",DHP$Feature) 
TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/data/TDM1_rank_features.csv")
TDM1=TDM1[-1,];names(TDM1)=c("Feature","Frequency");TDM1=TDM1[TDM1$Frequency>30,]
TDM1$group=sub("_.*", "",TDM1$Feature) 
group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0", 
  WSI = "#00743f"
)
#All
All$Feature=factor(All$Feature,levels =c("RNA_mRNA-ESR1","RNA_mRNA-PGR","RNA_mRNA-ERBB2","RNA_Exosome",                    
                                        "RNA_HER2DX_pCR_likelihood_score","DNA_COSMIC.Signature.6","WSI_Cell_Interaction",            
                                        "DNA_CX1","RNA_FCGR3B","RNA_sspbc.subtype"))
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
ggsave(All_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/All_feature_rank.pdf", width = 8.3, height = 5)
library(IOBR);library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")
df=as.data.frame(df)
RNA=c("RNA_mRNA-ESR1","RNA_mRNA-PGR","RNA_mRNA-ERBB2","RNA_Exosome","RNA_HER2DX_pCR_likelihood_score","RNA_FCGR3B","RNA_G0scores","pCR")
DNA=c("DNA_COSMIC.Signature.6","DNA_CX1","pCR")
WSI=c("WSI_Cell_Interaction","pCR")
data <- df[,RNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = RNA
)
data <- df[,DNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = DNA
)
#DHP
DHP$Feature=factor(DHP$Feature,levels =c("RNA_mRNA-ESR1","RNA_mRNA-ERBB2","RNA_mRNA-PGR","RNA_HER2DX_HER2_amplicon",
                                         "RNA_Purine_metabolism","RNA_sspbc.subtype",         
                                         "RNA_Th2 cells","RNA_Glutathione_metabolism","RNA_FCGR3A",
                                         "RNA_pik3ca_sig","Clin_ANYNODES","RNA_ABC_transporter","DNA_subclone_per"))
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
ggsave(DHP_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/DHP_feature_rank.pdf", width = 7.8, height = 5)

    
RNA=c("RNA_mRNA-ESR1","RNA_mRNA-ERBB2","RNA_mRNA-PGR","RNA_HER2DX_HER2_amplicon",
      "RNA_Purine_metabolism",   
      "RNA_Th2 cells","RNA_Glutathione_metabolism","RNA_FCGR3A",
      "RNA_pik3ca_sig","RNA_ABC_transporter","DNA_subclone_per","pCR")
data <- df[df$Clin_Arm=="DHP",RNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = RNA
)
#T-DM1
TDM1$Feature=factor(TDM1$Feature,levels =c("RNA_mRNA-PGR","Clin_ER","RNA_mRNA-CD8A","RNA_chr17q12_fusion","DNA_PPP1R1B_CNA",
                                           "RNA_mRNA-MKI67","RNA_Mast-cells","RNA_mRNA-ERBB2","RNA_Glycolysis","RNA_Hypoxia",
                                           "RNA_ERBB2_fusion","Clin_ANYNODES","RNA_Exosome","RNA_FCGR3B","DNA_LOH_Del_burden", 
                                           "RNA_ABC_transporter"))
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
ggsave(TDM1_feature, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/TDM1_feature_rank.pdf", width = 7, height = 5)

RNA=c("RNA_mRNA-PGR","RNA_mRNA-CD8A",
      "RNA_mRNA-MKI67","RNA_Mast-cells","RNA_mRNA-ERBB2","RNA_Glycolysis","RNA_Hypoxia",
      "RNA_Exosome","RNA_FCGR3B","RNA_ABC_transporter","pCR")
DNA=c("DNA_PPP1R1B_CNA","DNA_LOH_Del_burden","pCR")
data <- df[df$Clin_Arm=="T-DM1",RNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = RNA
)
data <- df[df$Clin_Arm=="T-DM1",DNA] %>% mutate(across(everything(), as.numeric))
data=na.omit(data)
res=batch_wilcoxon(
  data,
  target = "pCR",
  feature = DNA
)
data <- df[df$Clin_Arm=="T-DM1",]
chisq.test(data$pCR[data$Clin_ANYNODES!="Unknown"],
           data$Clin_ANYNODES[data$Clin_ANYNODES!="Unknown"]) 
#"RNA_chr17q12_fusion","RNA_ERBB2_fusion" "Clin_ER","Clin_ANYNODES"








library(data.table);library(tidyverse)
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
#image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image$V1=NULL
image$TREAT=NULL
image$ERPRdic=NULL
image$pCR=NULL
image=image[,c("patientID","V1","CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")]
zero_counts <- sapply(image, function(x) sum(x == 0))
cols_with_few_zeros <- names(zero_counts[zero_counts < 20])
length(cols_with_few_zeros)
varibale=intersect(cols_with_few_zeros,colnames(image)[1:ncol(image)-1])
#cols_to_keep <- grep("^SPA(?!.*(Max|Min))", colnames(image), perl = TRUE, value = TRUE)
#varibale=intersect(cols_with_few_zeros,cols_to_keep)
varibale=colnames(image)[1:ncol(image)-1]
image=image[,c("patientID",varibale)]
data=left_join(image,clin,by="patientID")%>%as.data.frame()
data[,varibale]=data[,varibale] %>% mutate(across(where(is.numeric), scale))
na_counts <- sapply(data[,varibale], function(x) sum(is.na(x)))
table(na_counts)
varibale=varibale[na_counts==0]
results=Logistic_batch_adjER(data,"pCR","Arm",varibale,"ER")%>%as.data.frame()
library(openxlsx)
write.xlsx(results, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/Image_FULL_Logistic_regression.xlsx")
write.xlsx(results, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure7/Image_multiplex_like_Logistic_regression.xlsx")


p.adjust(results$DHP_lr_p,method ="BH" )%>%min()
p.adjust(results$TDM1_lr_p,method ="BH" )%>%min()




library(data.table);library(tidyverse)
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
#image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image$V1=NULL
image$TREAT=NULL
image$ERPRdic=NULL
image$pCR=NULL
image=image[,c("patientID","CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")]
data=left_join(image,clin,by="patientID")%>%as.data.frame()
varibale=c("CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")
data=data[,c("patientID","pCR","Arm","ER","CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")]
data[,varibale]=data[,varibale] %>% mutate(across(where(is.numeric), scale))
results=Logistic_batch_adjER(data,"pCR","Arm",varibale,"ER")%>%as.data.frame()

str(data)






