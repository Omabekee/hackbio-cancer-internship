## Integrating Machine Learning and Differential Expression Analysis to Identify Key Biomarkers in Lymphoid Leukemia: Pathwways, Predictions and Future Directions

This project involves the identification of key biomarkers associated with primary and recurrent lymphoid leukemia (LL) samples, using differential gene expression analysis (DGE) and machine learning (Random forest). The analysis was performed by the OncoHackers team using several R libraries for data preprocessing, visualization, DGE analysis and ML classification.

### **Libraries Used**

1. **tidyverse**: A collection of R packages for data manipulation, organization and visualization.  
2. **cowplot**: for combining multiple plots into a single figure, enhancing visual presentations.  
3. **randomForest**: for implementing random forest algorithms for classification and regression in ML tasks.  
4. **TCGAbiolinks**: for accessing and analyzing data from The Cancer Genome Atlas (TCGA) for bioinformatics research.  
5. **caret**: for streamlining the process of training ML models, including cross-validation and tuning.  
6. **DOSE**: for performing disease ontology-based enrichment analysis, useful for interpreting gene expression results.  
7. **forcats**: makes it easier to work with categorical variables (factors) in R.  
8. **enrichplot**: for visualizing functional enrichment analysis results, such as Gene Ontology (GO) and pathway data.  
9. **clusterProfiler**: for statistical analysis and visualization of functional profiles for genes and gene clusters.  
10. **org.Hs.eg.db**: gene annotation package which provides mappings between human gene IDs and biological information.  
11. **AnnotationDbi**: for linking different annotation databases with R objects for easier access and usage. 
12. **EnhancedVolcano**: for generating enhanced volcano plots to visualize differential gene expression analysis.  
13. **ggplot2**: for creating custom, high-quality visualizations in R using a layered approach.  
14. **FactoMineR**: for conducting multivariate data analysis, including PCA and clustering.  
15. **factoextra**: for providing tools for visualizing multivariate data analyses, such as PCA and clustering.  
16. **biomaRt**: for querying BioMart databases (such as Ensembl) for retrieving gene annotations.  
17. **dplyr**: for data manipulation tasks like filtering, selecting, and summarizing data.  
18. **readr**: for reading and writing tabular data more efficiently.  
19. **EDASeq**: for normalizing and adjusting RNA-Seq data for sequencing depth and gene length.  
20. **SummarizedExperiment**: for managing experimental datasets with metadata and matrix-based assay data.

**Code:** To install packages, use `install.packages(“packagename”)`  
To load libraries after installation, use `library(packagename)` 

### Analysis Workflow

<figure>  
  <img src="Stage-3/DGE-ML-Roadmap.png" alt="Figure 7: My Biomarker Discovery Roadmap" width="600">  
  <figcaption>Figure 7: Biomarker Discovery Roadmap</figcaption>  
</figure>


The analysis workflow, including all steps, detailed methodology,  data preprocessing and results, are comprehensively documented within the [project report](https://github.com/Omabekee/hackbio-cancer-internship/blob/main/Stage-3/report/DGE-ML-for-biomarker-discovery.md).

### Directory Structure

**/code**: This directory contains all the R codes used for analysis for each step in the workflow.    
**/data**: In this folder, you will find the .csv format of the gene expression data downloaded from TCGA. 
**/figures**: This directory holds the figures generated during the analysis, including the volcano plot, confusion matrix and lollipop plots of the top pathways.    
**/ouput**: The output directory contains results from the project such as .csv files of the filtered significant genes obtained after analysis.    
**/report**: This is where you’ll find the project report documenting the steps, methodology and interpretations of the analysis.
