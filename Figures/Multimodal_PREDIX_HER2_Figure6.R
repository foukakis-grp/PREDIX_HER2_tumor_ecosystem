########Figure6B########
#T celll , CAF, meanHEAD/Neoantigen  3D PCA plot
library(data.table);library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr);library(sm)
mydata=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")%>%as.data.frame()
mydata=mydata[,c("patientID","TCR_clonality","TIDE_CAF","meanHED")]
mydata=na.omit(mydata)%>%as.data.frame()
mydata[,c("TCR_clonality","TIDE_CAF","meanHED")]=scale(mydata[,c("TCR_clonality","TIDE_CAF","meanHED")])
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin=clin[,c("patientID","Response","pCR","ER","Arm")]
data=left_join(mydata,clin,by="patientID")

model <- glm(as.numeric(pCR)~TIDE_CAF,family="binomial",data=data)
ShowRegTable(model)

outputDir="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/"
pdf(paste0(outputDir,"Fig6b_meanHED_density.pdf"),width=(5.38)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)

par(mfrow=c(2,1), oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("RD","pCR")){
  meanHED   <- data[data$Response %in% q,"meanHED"]
  TIDE_CAF <- data[data$Response %in% q,"TIDE_CAF"]
  y     <- cbind(`HLA evolutionary divergence` = meanHED,`CAF`=TIDE_CAF)
  sm.density(y,display="image",panel=F,ylim=c(-3,3),xlim=c(-3,3),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=median(data$meanHED),lty="dashed")
  abline(h=median(data$TIDE_CAF),lty="dashed")
}
mtext("CAF",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("HLA evolutionary divergence",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()
table(data$Response[data$meanHED>median(data$meanHE)&data$TIDE_CAF>median(data$TIDE_CAF)])
table(data$Response[data$meanHE>median(data$meanHE)&data$TIDE_CAF<median(data$TIDE_CAF)])
table(data$Response[data$meanHE<median(data$meanHE)&data$TIDE_CAF<median(data$TIDE_CAF)])
table(data$Response[data$meanHE<median(data$meanHE)&data$TIDE_CAF>median(data$TIDE_CAF)])


library(data.table);library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr);library(sm)
mydata=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")%>%as.data.frame()
mydata=mydata[,c("patientID","TIDE_CAF")]
mydata=na.omit(mydata)%>%as.data.frame()
image=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
image=image[,c("patientID","Immune_Cell_prop")]
data=inner_join(image,mydata,by="patientID")%>%as.data.frame()
data[,c("TIDE_CAF","Immune_Cell_prop")]=scale(data[,c("TIDE_CAF","Immune_Cell_prop")])
data$patientID=as.character(data$patientID)
colnames(data)=c("patientID","Immune_cell_proportion","CAF")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
clin=clin[,c("patientID","Response","pCR")]
data=left_join(data,clin,by="patientID")
outputDir="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/"
pdf(paste0(outputDir,"Fig6b_Immune_density.pdf"),width=(5.38)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)

par(mfrow=c(2,1), oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("RD","pCR")){
  Immune_cell_proportion   <- data[data$Response %in% q,"Immune_cell_proportion"]
  TIDE_CAF <- data[data$Response %in% q,"CAF"]
  y     <- cbind(Immune_cell_proportion= Immune_cell_proportion,CAF=TIDE_CAF)
  sm.density(y,display="image",panel=F,ylim=c(-3,3),xlim=c(-3,3),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=median(data$Immune_cell_proportion),lty="dashed")
  abline(h=median(data$CAF),lty="dashed")
}
mtext("CAF",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("Immune cell proportion",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()

model <- glm(as.numeric(pCR)~CAF+Immune_cell_proportion,family="binomial",data=data)
ShowRegTable(model)

table(data$Response[data$Immune_cell_proportion>median(data$Immune_cell_proportion)&data$CAF>median(data$CAF)])
table(data$Response[data$Immune_cell_proportion>median(data$Immune_cell_proportion)&data$CAF<median(data$CAF)])
table(data$Response[data$Immune_cell_proportion<median(data$Immune_cell_proportion)&data$CAF<median(data$CAF)])
table(data$Response[data$Immune_cell_proportion<median(data$Immune_cell_proportion)&data$CAF>median(data$CAF)])
########################
########Figure6D########
########################
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#RNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
colnames(data)=gsub("Danaher-","",colnames(data))
colnames(data)=gsub("TIDE_","",colnames(data))
data1=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune.rds")
data1=data1[,c("patientID","Th2 cells","MHC.I_19272155")]
data=left_join(data1,data,by="patientID")%>%as.data.frame()
variable=c("B-cells","DC","Macrophages","T-cells","CD8-T-cells",     
           "Neutrophils","Cytotoxic-cells","Treg",
           "Mast-cells","NK-cells","CD45","TILs","Dysfunction",        
           "Exclusion","CAF","TAM_M2","Th2 cells","MHC.I_19272155",
           "TCR_clonality","BCR_clonality","FCGR3A","FCGR3B")
data[,variable]=scale(data[,variable])
data=data[!is.na(data$CAF),]
rna_results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
#DNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
variable=c("A01","A02","A03","B44","B07","meanHED","lohhla","TCRA.tcell.fraction.adj","Neoantigen")
data$Neoantigen[is.na(data$Neoantigen)]=0
data[,c("meanHED","TCRA.tcell.fraction.adj")]=scale(data[,c("meanHED","TCRA.tcell.fraction.adj")])
data=data[!is.na(data$TCRA.tcell.fraction.adj),]
dna_results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()

#image
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
varibale=c("Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin=clin[,c("patientID","Response","pCR","ER","Arm")]
data=left_join(data,clin,by="patientID")
image_results=Logistic_batch_adjER(data,"pCR","Arm",varibale,"ER")%>%as.data.frame()
# merge
results=rbind(rna_results,dna_results)%>%
       rbind(image_results)
write.csv(results, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure5/Figure5d.csv")

TDM1=results[,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Signature","OR","Pvalue","group")
colnames(DHP)=c("Signature","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group=factor(df$group,levels = c("TDM1","DHP"))
unique(df$Signature)
df$Signature=factor(df$Signature,levels =rev(c("Neoantigen","meanHED","lohhla","MHC.I_19272155","A01","A03","B07","B44",
                                               "B-cells","DC","Macrophages","T-cells","CD8-T-cells","Cytotoxic-cells","NK-cells","CD45","TILs","Immune_Cell_prop",
                                               "TCR_clonality","BCR_clonality","TCRA.tcell.fraction.adj","FCGR3A","FCGR3B",
                                               "Dysfunction","Exclusion","CAF","Neutrophils","Mast-cells","Treg","TAM_M2","Th2 cells",
                                               "Distance_tumor_immune","Cell_Interaction")))
df=df[order(df$Signature,decreasing = F),]

baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Signature,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip()
# P 5.5X7
########################
########Figure6E########
########################
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
## continuous variable ##
# Distance_tumor_immune,Cell_Interaction,
# Th2 cells,FCGR3B Mast-cells Neutrophils  CAF
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
#RNA
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
colnames(data)=gsub("Danaher-","",colnames(data))
colnames(data)=gsub("TIDE_","",colnames(data))
data1=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune.rds")
data1=data1[,c("patientID","Th2 cells","MHC.I_19272155")]
data=left_join(data1,data,by="patientID")%>%as.data.frame()
variable=c("FCGR3B","CAF","Neutrophils","Mast-cells","Th2 cells")
data[,variable]=scale(data[,variable])
data=data[!is.na(data$CAF),]
immune=slice_metric(data,variable,20,80,5)%>%as.data.frame()
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(immune), value = TRUE)
immune=immune[,c("pCR","Arm","ER",selected_columns)]
Interact_result=Logistic_batch_adjER(immune,"pCR","Arm",selected_columns,"ER")%>%as.data.frame()
immune$Arm=factor(immune$Arm,levels = c("DHP","T-DM1"))
res=Logistic_batch_continuous_subgroup(immune,selected_columns)%>%as.data.frame()
res_continuous=res[res$biomarker%in%c("FCGR3B_per_20","CAF_per_80","Mast-cells_per_35","Neutrophils_per_30","Th2 cells_per_60"),]
Interact_continuous=Interact_result[Interact_result$biomarker%in%c("FCGR3B_per_20","CAF_per_80","Mast-cells_per_35","Neutrophils_per_30","Th2 cells_per_60"),]
continuous=cbind(res_continuous,Interact_continuous)
#Digital no good results
#HLA supertype  Manually add
# merge
library(openxlsx);library(tableone)
list_of_datasets <- list("immune" = continuous)
write.xlsx(list_of_datasets, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Predictive_biomarker.xlsx")
table(genomic$Hypoxia_per_80,genomic$Arm)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$A01=="No",])
ShowRegTable(whole)
whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$A01=="Yes",])
ShowRegTable(whole)
table(data$A01,data$Arm)
# forest plot #
library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure6/Subgroup_Forestplot.csv")
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 20), collapse = " ")
a=df$`OR (95% CI)`;b=df$`P for interaction`
df$`OR (95% CI)`=NULL;df$`P for interaction`=NULL
df$`OR (95% CI)`=a
df$`P for interaction`=b

tm <- forest_theme(base_size = 10,
                   refline_col = "black",
                   arrow_type = "closed",
                   footnote_col = "blue")

p <- forest(df[,c(1:3,7:9)],
            est = df$OR,
            lower = df$Low, 
            upper = df$High,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("DHP Better", "T-DM1 Better"),
            xlim = c(0, 4),
            ticks_at = c(0,0.5, 1, 2, 3),
            theme = tm)

# Print plot
plot(p)

#7X5







