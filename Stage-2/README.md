## Gene Expression Analysis and Visualisation of Glioblastoma Data

This project involves a walk-through of gene expression analysis of glioblastoma data conducted by the OncoHackers team, using various R libraries and visualisation techniques. Below are the steps and functions used throughout the analysis.

### Libraries Used

1. **gplots**: for generating heatmaps as well as other plots.  
2. **ggplot2**: for creating various types of visualisations.  
3. **RColorBrewer**: for diverging and sequential colour palettes in the heatmaps.  
4. **dplyr**: for data manipulation such as filtering, selecting and summarising data.  
5. **biomaRt**: for querying BioMart databases (such as Ensembl) for retrieving gene annotations.
6. **tidyverse**: is a set of R packages that includes ggplot2 for making plots, dplyr for handling data, tidyr for organizing data, and more.  
   

**Code:** To install packages, use *install.packages(“packagename”)*  
To load libraries after installation, use  *library(packagename)* 

### Analysis Workflow

1. **Reading in Glioblastoma Data:** The glioblastoma dataset contains the expression profiles of glioblastoma samples, where each value shows a gene's expression level in a specific sample. To begin analysis, load the gene expression data for glioblastoma samples.  
     
   **R Function**: *read.csv()* is a function used to read a CSV file into a dataframe.  
     
2. **Heatmap Generation**: Heatmaps were generated to visualise the gene expression patterns in the samples using either diverging colour palette or sequential colour palette. This approach aids in showing variations between upregulated and downregulated genes as well as gradual changes in expression levels.  

3. **Clustering Methods in the Heatmap**: Performed clustering by rows (genes), columns (samples) and both. This approach helps in understanding which samples are similar and which genes have similar expression patterns across these samples.

     i). **Clustering Genes (Rows)**: helps identify gene sets with similar behaviour, which 
     may indicate co-regulated genes or genes involved in the same biological processes.

     ii). **Clustering Samples (Columns)**: helps to find subtypes within the dataset based 
     on their gene expression patterns.

     iii). **Clustering Both Genes and Samples (Rows and Columns)**: provided a 
     comprehensive view of the data, showing how gene groups and sample profiles relate to 
     each other.
  
     
     **R Function**: *heatmap.2()* is a function used to generate heatmaps.  
     *colorRampPalette* used to create custom colour palettes to visually distinguish data 
     patterns.  
     
5. **Identification of Significant Genes:** Using the clusters generated from figure above, the samples were divided into two groups. Subsequently, log fold change was calculated by taking the log2 difference between the mean gene expression of Group B and Group A, with a small constant added for stability.  
   **LogFC formula**: *log2(groupB\_mean \+ 0.5) \- log2(groupA\_mean \+ 0.5)*  
     
   The p-value was calculated using the Wilcoxon test to compare the mean gene counts between Group A and Group B.  
   **p-value formula:** *wilcox.test(mean gene counts in groupA, mean gene counts in groupB*  
     
   Next, we determined significantly upregulated and downregulated genes based on log fold change (\> 1.5 or \< \-1.5) and p-values (\<0.05).  
     
   The Gene IDs of both up- and downregulated genes were mapped to their gene symbols from Ensembl.  
     
   **R Function**: *useMart()* is used to connect to a BioMart database for querying.  
   *getBM()* is used to retrieve data from the connected BioMart database.  
     
6. **Functional Enrichment Analysis**: Used ShinyGO (with default parameters and GO biological process database) to identify pathways enriched in the upregulated genes. The results obtained from ShinyGO were saved in a .csv file and used for visualisation (find file in *output* directory).  
     
7. **Visualisation of Pathways**: Visualised the top 5 upregulated pathways using the results obtained from the enrichment analysis to generate a lollipop plot (scaled by the \-log10 of the p-value to indicate pathway significance) for better clarity and interpretation of the enriched pathways.  
   **R Function**: *ggplot()* is used to create various plots in R, such as scatter plots, bar plots and lollipop plots in R, with customizations for layers, scales and  themes.

### Biological Significance of Top Upregulated Pathways
In summary, the top three pathways all point to the adaptive ability of tumour cells to ensure their survival and proliferation, leading to metastasis.
   
### Directory Structure

* **/code**: This directory contains all the R codes used for data analysis and visualisation, for each step in the workflow.  
* **/data**: In this folder, you will find the raw gene expression data used for analysis.  
* **/figures**: This directory holds the visualisations generated during the analysis, including the heatmaps and the lollipop plot of the top pathways.  
* **/ouput**: The output directory contains results from the project such as tables with the filtered significant genes and ShinyGO results obtained after enrichment analysis.  
* **/report**: This is where you’ll find the project report documenting the steps and interpretations of the analysis.

#### Visit live site
You can view the live version of this project [here](https://omabekee.github.io/hackbio-cancer-internship/)

  

