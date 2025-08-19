# PREDIX_HER2_tumor_ecosystem

## Figures: notebook code for each figure.
## scr: source code.
## data: processed data.


# Multi-omics metrics


Malignant epithelial cells are the most heterogeneous cell type with almost every patient forming a separate cluster. Here, we present a method, AI-EPI (**A**tlas-level **I**ntegrated **E**pithelial **P**rogram **I**dentification), which identify patient-shared and patient-specific gene modules (GM) simultaneously and efficiently.The
method mainly contains two steps:

-   source code 
-   resouce files
-   code for data pre-processing
-   code for figure

<div align=center> 
<img src="./inst/workflow.png" alt="workflow.png">
</div> 

nf-core workflow
------------

The nf-core pipelines can be set up based on: https://nf-co.re/pipelines/

```
# Launch the RNAseq pipeline
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --output ./results/ \
    --genome GRCh38 \
    -profile singularity
```

Quick start
-----------

Here, we provide an example data of [GBC_epithelial](http://lifeome.net/supp/gbc/Single%20cell/Epithelial/adata_adeno_P.h5ad) from 10X Genomics. Users can download it and run following scripts to understand the workflow of AIEPI.

## Step1: gene module identification

For malignant epithelial cell number vary among patients, we sample the same cell number from every patient so that they are equally weighted.

```
import AI_EPI
sample_number = 500
weighted_sample = downsampling(epithelial_adata.obs,sample_number)
epithelial_downsample_adata = epithelial_adata[weighted_sample.index,]
```
AI-EPI identifies gene modules by consensus non-negative matrix factorization (cNMF). You can select a appropriate pragram number by the curve of stability and error at each choice of K. The down-sampled data of the exampled data can be obtained in [downsampled_data](http://lifeome.net/supp/gbc/Single%20cell/Epithelial/adata_adeno_p_sample.h5ad).The detailed gene module identification code based on down-sampled data can be obtained in [GM_identification.ipynb](https://github.com/JulieBaker1/AIEPI/blob/main/code/1.GM_identification.ipynb).

<div align=center> 
<img src="./inst/Epithelial.k_selection.png" width = "300" alt="Epithelial.k_selection.png">
</div> 


## Step2: gene module classification

In the second step, we distinguish the patient-shared GM from the patient-specific GM by a permutation test p-value.  

```
source("code/2.GM_classification.R")
patient_GM_score_23GM_100genes = read.csv("./data/patient_GM_score_23GM_100genes.csv",row.names = 1)
GM_classification_result = GM_classification(patient_GM_score_23GM_100genes)
```
<div align=center> 
<img src="./inst/IQR.png" width = "300"  alt="IQR.png">
</div> 

## Downstream analysis

We can define the state of each cell by the GM with highest GM score. The code can be obtained in [cell_state_assignment.ipynb](https://github.com/JulieBaker1/AIEPI/blob/main/code/3.cell_state_assignment.ipynb).


<div align=center> 
<img src="./inst/cell_state.png" width = "300" alt="cell_state.png">
</div> 





