# Integrating Machine Learning and Differential Expression Analysis to Identify Key Biomarkers in Lymphoid Leukemia: Pathways, Predictions and Future Directions

#### Authors (@slack): Chairunnisa Amanda (@AmandaC), Bassam Elhamsa (@Bassam_Elhamsa), Chioma Onyido (@Omabekee), Chukwuemeka Nwachuya (@Emzy247), Emmanuel Afolayemi (@EA(NOVACHRONO)), Hayford Osei Offei (@Hayford72), Oluwatobi Ogundepo (@Oluwatobi), Ibrahim Fangary (@IbrahimFangary)

#### GitHub Code: [DGE-ML-biomarker-analysis R Script](https://github.com/Omabekee/hackbio-cancer-internship/blob/main/Stage-3/code/DGE-ML-Analysis.R)

---

## Table of Contents
1. [Introduction](#1-introduction-to-lymphoid-leukemias)
2. [Dataset and Data Preprocessing](#2-description-of-dataset-and-data-preprocessing-steps)
3. [Methodology for Biomarker Discovery](#3-methodology-for-biomarker-discovery)
4. [Methodology for Machine Learning Analysis](#4-methodology-for-machine-learning-analysis)
5. [Result and Interpretation](#5-result-and-interpretation-of-model-performance)
6. [Conclusion and Future Directions for Research](#6-conclusion-and-future-directions-for-research)
7. [References](#references)

---
## 1. Introduction to Lymphoid Leukemias

Lymphoid leukemias (LL) are blood cancers resulting from the abnormal growth of lymphoid cells- B, T, or NK cells. They are categorised into acute lymphoblastic leukaemia (ALL), common in children, and chronic lymphocytic leukaemia (CLL), more common in adults, especially the elderly (Chennamadhavuni _et al_., 2023).


### 1.1 Project Aim:

To identify the key biomarkers associated with primary and recurrent LL samples using differential gene expression analysis and machine learning.

## 2. Description of Dataset and Data Preprocessing Steps

The data was downloaded from The Cancer Genome Atlas (TCGA) database via the GDC data portal. We selected 25 primary and 25 recurrent samples for analysis.

### 2.1 Handling Missing Values

The LL dataset was preprocessed by checking for missing or blank values using the `is.na` function, to reduce redundancies and understand what is unexpressed/undetected in our data.

### 2.2 Normalisation and Filtering

Normalisation and filtering were performed using the `TCGAnalyze_Normalization`, `TCGAnalyze_Filtering` and `betweenlaneNormalization` functions from the TCGAbiolinks and EDASeq R packages to adjust for gene length and sequencing depth (Mounir _et al._, 2019).

## 3. Methodology for Biomarker Discovery

### 3.1 Differential Gene Expression Analysis (DGE)

The analysis was conducted using the `TCGAanalyze_DEA` function from the TCGAbiolinks R package. Comparison was made between "primary" and "recurrent" samples, filtering results by an adjusted p-value < 0.05 and log2 fold change > 1.

<figure>  
  <img src="figures/fig2.png" alt="Figure 2:  Volcano plot showing the significant genes between primary and recurrent LL samples" width="800">  
  <figcaption>Figure 1:  Volcano plot showing the significant genes between primary and recurrent LL samples</figcaption>  
</figure>



### 3.2 Functional Enrichment Analysis

Functional enrichment was performed on 104 upregulated and 1,305 downregulated genes using the `enrichGO` function in R.

### 3.3 Pathway Visualisation

The steps involved filtering pathways based on p-value and q-value, calculating the gene ratio and rich factor for each pathway, and visualising the top 20 enriched pathways using `ggplot`. 


<figure>  
  <img src="figures/fig3.png" alt="Figure 3:  Key Pathways and Biological Processes Associated With Upregulated Genes" width="600">  
  <figcaption>Figure 2:  Key Pathways and Biological Processes Associated With Upregulated Genes</figcaption> 
</figure>


<figure>  
  <img src="figures/fig4.png" alt="Figure 4:  Key Pathways and Biological Processes Associated With Downregulated Genes" width="600">  
  <figcaption>Figure 3:  Key Pathways and Biological Processes Associated With Downregulated Genes</figcaption> 
</figure>


## 4. Methodology for Machine Learning Analysis

### 4.1 Feature Extraction

After performing differential gene expression analysis, genes were filtered by selecting those with a LogFC > 1 and adjusted p-values < 0.05 which reduced the total number of genes from **14,898** to **1,502**. Subsequently, near-zero variance filtering was applied with no genes removed. Finally, to avoid collinearity, genes with a correlation higher than 80% were removed, resulting in 364 genes.

### 4.2 Model Training and Testing.

A random forest classification model was built to classify sample type—either primary or recurrent— using the feature-selected training dataset consisting of 364 genes and 20 samples (10 primary and 10 recurrent). The model was configured with 500 trees in the forest (`ntree = 500`) and 27 genes considered at each split (`mtry = 27`). Model testing and validation were performed on an independent set of 10 samples (5 primary and 5 recurrent).

## 5. Result and Interpretation of Model Performance

After the analysis, the model achieved a prediction accuracy of 80%. Out of the 10 samples used for the testing, the model accurately predicted 4/5 Primaries and 4/5 Recurrent cancer samples

<figure>  
  <img src="figures/fig5.png" alt="Figure 5:  Confusion matrices (CM) summarising the performance of the model on the test data to give insight into the precision of the classification approach" width="600">  
  <figcaption>Figure 4:  Confusion matrices (CM) summarising the performance of the model on the test data to give insight into the precision of the classification approach</figcaption> 
</figure>

<figure>  
  <img src="figures/fig6.png" alt="Figure 6:  Top 20 genes by gini importance that helped the model performance" width="600">  
  <figcaption>Figure 5:  Top 20 genes by gini importance that helped the model performance</figcaption> 
</figure>


## 6. Conclusion and Future Directions for Research

This project combined machine learning and differential expression analysis to identify key biomarkers in LL. Through DEG and functional enrichment analysis, we were able to identify specific molecular changes between the two stages- primary and recurrent. The analysis revealed important genes like AKR1C3, ARHGEF11 and AHNAK which show potential as therapeutic and diagnostic biomarkers for early detection and more effective therapies. Although our random forest classifier achieved 80% prediction accuracy, some primary and recurrent samples were misclassified, indicating areas for improvement. Future directions include refining the model for higher accuracy, expanding sample sizes and integrating additional data types to provide a more comprehensive understanding of LL.

---

## References

1. Chennamadhavuni, A., Lyengar, V., Mukkamalla, S.K.R. & Shimanovsky, A. (2023). Leukemia. In: StatPearls [Internet]. Treasure Island (FL): StatPearls Publishing. Available from: https://www.ncbi.nlm.nih.gov/books/NBK560490/ 
2. Mounir, M., Lucchetta, M., Silva, T. C., Olsen, C., Bontempi, G., Chen, X., ... & Papaleo, E. (2019). New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx. *PLoS computational biology*, *15*(3), e1006701.  
   
