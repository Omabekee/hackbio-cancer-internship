## ========================================

#             Differential expression and machine learning in discovery of
#             Lymphoid leukemia biomarkers  

# ========================================


# ----------------------------------------

#             accessing all data and downloading Lymphoid leukemia data

# ----------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)

# Query for Lymphoid leukemia data from TARGET-ALL-P2 project
query_expression <- GDCquery(
  project = "TARGET-ALL-P2",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Blood Derived Cancer - Bone Marrow", "Recurrent Blood Derived Cancer - Bone Marrow")
)

# Extracting 25 Primary Blood Derived Cancer - Bone Marrow samples from the total query of the project

# creating a TRUE, FALSE array to index the 25 Primary Blood Derived Cancer - Bone Marrow samples
prim_index <- c()
# prim is a holder that will be used in a for loop to get just 25 TRUE 
prim <- 1
# For loop to extract the indicies of the first 25 25 Primary Blood Derived Cancer - Bone Marrow samples of our query
for (i in 1:nrow(query_expression$results[[1]])){
  
  if (query_expression$results[[1]]$sample_type[i] == "Primary Blood Derived Cancer - Bone Marrow" & prim <= 25) {
    
    prim_index[i] <- TRUE
    prim = prim+1
  } else {
    prim_index[i] <- FALSE
  }
  
}  

# primary is a dataframe has 25 Primary Blood Derived Cancer - Bone Marrow samples 
primary <- query_expression$results[[1]][prim_index,]


# Extracting 25 Recurrent Blood Derived Cancer - Bone Marrow samples from the total query of the project

# creating a TRUE, FALSE array to index the 25 Recurrent Blood Derived Cancer - Bone Marrow samples
rec_index <-  c()
# rec is a holder that will be used in a for loop to get just 25 TRUE 
rec <- 1
# For loop to extract the indicies of the first 25 Recurrent Blood Derived Cancer - Bone Marrow samples of our query
for (i in 1:nrow(query_expression$results[[1]])){
  
  if (query_expression$results[[1]]$sample_type[i] == "Recurrent Blood Derived Cancer - Bone Marrow" & rec <= 25) {
    
    rec_index[i] <- TRUE
    rec = rec+1
  } else {
    rec_index[i] <- FALSE
  }
  
}

# recurrent is a dataframe has 25 Recurrent Blood Derived Cancer - Bone Marrow samples 
recurrent <- query_expression$results[[1]][rec_index,]
# binding both 25 Recurrent Blood Derived Cancer - Bone Marrow samples and 25 Primary Blood Derived Cancer - Bone Marrow 
res <- rbind(primary,recurrent)

# update our query to get just 25 Recurrent Blood Derived Cancer - Bone Marrow samples and 25 Primary Blood Derived Cancer - Bone Marrow samples
query_expression$results[[1]] <- res


# Query for metadata
query_metadata <- GDCquery(
  project = "TARGET-ALL-P2",
  data.category = "Clinical",
  data.type = "Clinical Supplement"
)


GDCdownload(query_expression)

# Download the metadata
GDCdownload(query_metadata)

# Prepare the expression data for analysis
LL_expression_Data <- GDCprepare(query_expression)

# Prepare the metadata for analysis
LL_metadata <- GDCprepare(query_metadata)

# Gene expression dataframe 
gene_expression <- assay(LL_expression_Data)

write.csv(gene_expression, file = "gene_expression.csv", row.names = TRUE)

# ----------------------------------------

#              pre-processing and normalization

# ----------------------------------------


# This is an extension of "Accessing_ALL_data.R" code

# Load the libraries
library(readr)
library(TCGAbiolinks)
library(EDASeq)


# samplesheet_LL.tsv was downloaded from TCGA database directly 
samplesheet_LL <- read_tsv("samplesheet_LL.tsv")

gene_expression<- read_tsv("gene_expression.csv")

# Matching the IDs in 'gene_expression' data to the sample type (whether primary or recurrent)
# Extract the desired part (up to the third hyphen) from the column names
gene_expression_prefix <- sub("^(([^-]+-){3}[^-]+)-.*$", "\\1", colnames(gene_expression))

# Check the result
print(gene_expression_prefix)


# Match the gene expression prefixes with Sample ID in samplesheet_LL
matched_metadata <- samplesheet_LL[samplesheet_LL$`Sample ID` %in% gene_expression_prefix, ]

# Filter for Primary and Recurrent Blood Derived Cancer - Bone Marrow samples
primary_metadata <- matched_metadata[matched_metadata$`Sample Type` == "Primary Blood Derived Cancer - Bone Marrow", ]
recurrent_metadata <- matched_metadata[matched_metadata$`Sample Type` == "Recurrent Blood Derived Cancer - Bone Marrow", ]

write.csv(primary_metadata, "primary_metadata.csv", row.names=FALSE)
write.csv(recurrent_metadata, "recurrent_metadata.csv", row.names=FALSE)


# NEXT- Create a mapping from primary and recurrent metadata

# Combine primary and recurrent metadata into one table for easy mapping
combined_metadata <- rbind(primary_metadata, recurrent_metadata)

# Create a map from Sample ID to Sample Type ("Primary" or "Recurrent")
sample_type_map <- setNames(combined_metadata$`Sample Type`, combined_metadata$`Sample ID`)

# Replace gene expression column names with their corresponding Sample Type ("Primary" or "Recurrent")
# Match the extracted Sample IDs to the `Sample Type` in the combined metadata
new_colnames <- sample_type_map[gene_expression_ids]

# Assign the new column names back to the gene expression data
colnames(gene_expression) <- new_colnames


# Write the updated gene expression data to a CSV file
write.csv(gene_expression, file = "updated_gene_expression.csv", row.names = TRUE)


#NEXT- Data Preprocessing
LL_exp_data <- read.csv("updated_gene_expression.csv", row.names = 1)

#Check for missing values in gene expression data
any(is.na(gene_expression))    #this returned false indicating that there is no missing values
#check for zeros in the dataframe
sum(LL_exp_data == 0)      #this returned 1276958



# Normalize the data using gene length normalization
normData <- TCGAanalyze_Normalization(tabDF = LL_exp_data,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")


# Next, to correct for sequencing depth, apply another normalization method-
# upper quartile normalization.
# Convert normData to a matrix
normDataMatrix <- as.matrix(normData)


# Upper quartile normalization
normDataDepthAdjusted <- EDASeq::betweenLaneNormalization(normDataMatrix, which = "upper")

# Convert back to a data frame if needed
normDataDepthAdjusted <- as.data.frame(normDataDepthAdjusted)

# Now filter out low counts
filteredData <- TCGAanalyze_Filtering(tabDF = normDataDepthAdjusted, 
                                      method = "quantile", 
                                      qnt.cut = 0.25)

#still contains zero values, reduced to 24464 genes
# Check for any remaining zero rows
# Remove rows with at least one zero
filtered_out_zeros <- LL_exp_data[rowSums(LL_exp_data == 0) == 0, ]

# Check the result
print(dim(filtered_out_zeros)) #reduced to 14898 genes

#write to csv row.names = true
write.csv(filtered_out_zeros, file = "updated_preprocessed_data.csv", row.names = TRUE)




# ----------------------------------------

#             differential expression analysis

# ----------------------------------------

# Load the necessary library
library(TCGAbiolinks)
library(dplyr)

# Step 1: Load preprocessed data
LL_exp_data <- read.csv("updated_preprocessed_data.csv", row.names = 1)

# Check the column names
head(colnames(LL_exp_data))

# Step 2: Subset based on the column names
# Find column names that contain "Primary" or "Recurrent"
primary_samples <- grep("Primary", colnames(LL_exp_data), value = TRUE)
recurrent_samples <- grep("Recurrent", colnames(LL_exp_data), value = TRUE)

# Select the first 20 "Primary" and 20 "Recurrent" samples
primary_subset <- primary_samples[1:20]
recurrent_subset <- recurrent_samples[1:20]

# Combine selected samples
selected_samples <- c(primary_subset, recurrent_subset)

# Subset the expression data to only include the selected samples
subset_exp_data <- LL_exp_data[, selected_samples]

# Step 3: Prepare a condition/group vector
group_condition <- ifelse(grepl("Primary", selected_samples), "Primary", "Recurrent")

# Step 4: Run DEA using TCGAbiolinks
DEresults <- TCGAanalyze_DEA(
  mat1 = subset_exp_data[, group_condition == "Primary"], # Primary samples matrix
  mat2 = subset_exp_data[, group_condition == "Recurrent"], # Recurrent samples matrix
  Cond1type = "Primary", # Condition 1: Primary
  Cond2type = "Recurrent", # Condition 2: Recurrent
  method = "glmLRT"  # Method to use, "glmLRT" is recommended for RNA-Seq
)

# Step5: Filter significant results (adjusted p-value < 0.05 and log2 fold change > 1)
significant_genes <- DEresults %>%
  filter(FDR < 0.05 & abs(logFC) > 1)

write.csv(significant_genes, "DE_sig_genes.csv")

#----
# Get Gene Symbol from Ensembl
library(biomaRt)

# Step 1: Clean Ensembl IDs by removing the version numbers
cleaned_ensembl_ids <- sub("\\..*", "", rownames(DEresults))

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Step 2: Connect to Ensembl
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = cleaned_ensembl_ids,
  mart = ensembl)

# Step3: Check for gene IDs without symbols
missing_symbols <- cleaned_ensembl_ids[!cleaned_ensembl_ids %in% gene_annotations$ensembl_gene_id]

# The gene IDs without symbols are from archive
# Connect to the GRCh37 archive and query for missing symbols:
ensembl_grch37 <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "grch37.ensembl.org"
)

# Query biomaRt to get gene symbols for GRCh37
gene_annotations_grch37 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = missing_symbols,
  mart = ensembl_grch37
)

# Step 4: Combine the annotations
combined_annotations <- rbind(gene_annotations, gene_annotations_grch37)

DEresults$Ensembl_ID <- cleaned_ensembl_ids

#------------------------------------------------------------------
#           Merge the annotations with the DEA results by Ensembl ID
# ------------------------------------------------------------------
DEresults_with_symbols <- merge(DEresults, combined_annotations, 
                                by.x = "Ensembl_ID", 
                                by.y = "ensembl_gene_id", 
                                all.x = TRUE)

# Step 1: Remove rows with no hgnc_symbol (i.e., NA values in the hgnc_symbol column)
DEresults_with_symbols <- DEresults_with_symbols[!is.na(DEresults_with_symbols$hgnc_symbol), ]

# Step 2: Remove the columns "gene_name" and "gene_type" if they exist
# Check if the columns "gene_name" and "gene_type" are present, then remove them
if ("gene_name" %in% colnames(DEresults_with_symbols)) {
  DEresults_with_symbols <- DEresults_with_symbols[ , !(colnames(DEresults_with_symbols) %in% c("gene_name", "gene_type"))]
}

# Step 3: Make hgnc_symbol unique
DEresults_with_symbols$hgnc_symbol <- make.unique(as.character(DEresults_with_symbols$hgnc_symbol))

# Step 4: Set "hgnc_symbol" as row names and remove the "hgnc_symbol" column
rownames(DEresults_with_symbols) <- DEresults_with_symbols$hgnc_symbol
DEresults_with_symbols$hgnc_symbol <- NULL  # Remove the column after setting rownames

# Step 5: Save the cleaned result
write.csv(DEresults_with_symbols, "DE_results.csv")


#--------------------------------------------------

#                         Visualisation

#--------------------------------------------------
res <- DEresults_with_symbols
#or
res <- read.csv("DE_results.csv", row.names = 1)

# Volcano Plot

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

png("DE_volcano_plot.png", width = 1280, height = 720)
EnhancedVolcano(
  res, 
  lab = rownames(res), 
  x = 'logFC', 
  y = 'FDR', 
  pCutoff = 10e-3, 
  FCcutoff = 1.5, 
  title = 'Volcano plot of DE results',
  xlim = c(-3, 3) 
)
dev.off()


# Extract significant genes (you can adjust the filtering criteria as needed)
significant_genes <- read.csv("DE_sig_genes.csv", row.names = 1)

top50_genes <- significant_genes %>%
  arrange(desc(abs(logFC))) %>%
  head(50)

# Extract the expression matrix for the significant genes
expression_matrix <- subset_exp_data[rownames(subset_exp_data) %in% rownames(top50_genes), ]

# Create annotation for samples (group: Primary or Recurrent)
annotation <- data.frame(Group = group_condition)
rownames(annotation) <- colnames(expression_matrix)

# Generate the heatmap
pheatmap(expression_matrix, 
         annotation_col = annotation,  # Add the group annotation
         scale = "row",                # Scale rows to make expression levels comparable
         show_colnames = FALSE,
         show_rownames = FALSE, # Option to hide gene names if too many
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

# Perform PCA using normalized expression data
# Assuming expression_data is a matrix with rows as genes and columns as samples
pca_res <- prcomp(t(subset_exp_data), scale. = TRUE)

# Prepare a data frame for ggplot2
pca_df <- as.data.frame(pca_res$x)
pca_df$group <- factor(c(rep("Primary", X), rep("Recurrent", Y)))  # Assuming X primary and Y recurrent samples

# Create a PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +                                    # Set point size
  theme_minimal() +                                         # Set theme
  labs(title = "PCA of Gene Expression Data",               # Add title
       x = paste0("PC1: ", round(100 * summary(pca_res)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca_res)$importance[2, 2], 1), "% variance")) +
  theme(legend.position = "right")                          # Position the legend on the right


# Step 1: Load necessary libraries for PCA
library(FactoMineR)
library(factoextra)

# Step 2: Transpose the data matrix so that genes are columns, and samples are rows
# Assuming `subset_exp_data` is the expression data matrix
pca_data <- t(subset_exp_data)

# Step 3: Perform PCA
pca_result <- PCA(pca_data, graph = FALSE)

# Step 4: Visualize PCA (2D plot)
fviz_pca_ind(pca_result, 
             geom.ind = "point",  # Show points only (individuals)
             col.ind = group_condition,  # Color by groups
             palette = c("#00AFBB", "#FC4E07"),  # Color palette for primary and recurrent
             addEllipses = TRUE,  # Add confidence ellipses
             legend.title = "Groups",
             title = "PCA - Primary vs Recurrent Samples")



# ----------------------------------------

#             pathway enrichment analysis

# ----------------------------------------

# loading packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggplot2)

DE_results<- read_csv('DE_results.csv')

#looking at the data 
head(DE_results)

# extracting up-regulated genes 
UPgenes <-  DE_results %>%
  filter(logFC > 1.5 & FDR < 10e-3)
dim(UPgenes)
####there are 104 up regualted genes based on the filtering criteria

# extractingdown-regulated genes 
downgenes <-  DE_results %>%
  filter(logFC < -1.5 & FDR < 10e-3)
dim(downgenes)
#### we have 1305downregulated genes


# extracting gene names 
UPgenes_n <- UPgenes['Ensembl_ID']# the Ensembl_ID corresponds to the coulmn name of the genes ID
downgenes_n <-downgenes['Ensembl_ID'] 

# Check for missing or duplicated values
any(is.na(downgenes_n$Ensembl_ID))   # Should return FALSE
any(duplicated(downgenes_n$Ensembl_ID))  # Should return FALSE
# Check for missing or duplicated values
any(is.na(UPgenes_n$Ensembl_ID))   
any(duplicated(UPgenes_n$Ensembl_ID)) 

# removing dublucated names
downgenes_n <-downgenes_n[!duplicated(downgenes_n), ]# here we need to remove dublicated gene ENSEBML ID
# extracting the genes ID
downgenes_n <- as.vector(downgenes_n$Ensembl_ID)
UPgenes_n <- as.vector(UPgenes_n$Ensembl_ID)


GO_result_up<- enrichGO(gene= UPgenes_n,
                        OrgDb= 'org.Hs.eg.db', #Hs for human annotation 
                        keyType = 'ENSEMBL', # i used the ENSEBML gene ID
                        ont= 'Bp' ) # for biological processes

as.data.frame(GO_result_up) #looking at what we have
# dump plot
plot(barplot(GO_result_up, showCategory = 10))



GO_result_down<- enrichGO(gene=downgenes_n,
                          OrgDb= 'org.Hs.eg.db', #Hs for human annotation 
                          keyType = 'ENSEMBL', # i used the ENSEBML gene ID
                          ont= 'Bp' ) # for biological processes

as.data.frame(GO_result_down) #looking at what we have
# dump plot
plot(barplot(GO_result_down, showCategory = 5))

#visualization of the result

#visualziation of up-regulated genes
modified_result_up<- filter(GO_result_up, p.adjust < .05, qvalue < 0.2) # filtering pathways


library(DOSE)
library(forcats)
library(enrichplot)
modified_result_up= mutate(modified_result_up, geneRatio = parse_ratio(GeneRatio)) %>% #converting to numeric value ratio 
  arrange(desc(geneRatio)) # sorting the dataframe

# k/M
# Calculate richFactor by dividing Count by the numerator of BgRatio
y <- mutate(modified_result_up, richFactor = Count / as.numeric(sub("/\d+", "", BgRatio)))
ggplot(y, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Up regulated pathways")

#visualization of down-regulation genes genes
modified_result_down<- filter(GO_result_down, p.adjust < .05, qvalue < 0.2)# filtering pathways

modified_result_down= mutate(modified_result_down, geneRatio = parse_ratio(GeneRatio)) %>% #converting to numeric value ratio
  arrange(desc(geneRatio))# sorting the dataframe

# k/M
# Calculate richFactor by dividing Count by the numerator of BgRatio
y_ <- mutate(modified_result_down, richFactor = Count / as.numeric(sub("/\d+", "", BgRatio)))
ggplot(y_, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Down regulated pathways")


# ------------------------------------------

#             machine learning: Random Forest

# ------------------------------------------

# loading libraries 
library(tidyverse)
library(cowplot)
library(randomForest)
library(TCGAbiolinks)
library(caret)

# reading data and split it for training and testing data
training_data<- read.csv("updated_preprocessed_data.csv",row.names = 1)
test_data <- training_data[,c(21:25,46:50)]
training_data <- training_data[,c(1:20,26:45)]

# creating vectors for samples classes (Amanda code)
group_condition_training <- ifelse(grepl("Primary", colnames(training_data)), "Primary", "Recurrent")
group_condition_testing <- ifelse(grepl("Primary", colnames(test_data)), "Primary", "Recurrent")

# perform differential expression analysis (Amanda code)
DEresults <- TCGAanalyze_DEA(
  mat1 = training_data[, group_condition_training == "Primary"], # Primary samples matrix
  mat2 = training_data[, group_condition_training == "Recurrent"], # Recurrent samples matrix
  Cond1type = "Primary", # Condition 1: Primary
  Cond2type = "Recurrent", # Condition 2: Recurrent
  method = "glmLRT"  # Method to use, "glmLRT" is recommended for RNA-Seq
)

# filter significant upregu;ated and downregulated genes
significant_genes <- DEresults %>%
  filter(FDR < 0.05 & abs(logFC) > 1.5)

# Vector of significant genes names
significant_genes_names <- row.names(significant_genes)

# preparing training and test data
FS_training_data <- training_data[significant_genes_names,] # Extracting significant genes out of training data
FS_training_data <- t(FS_training_data) # Transpose training data
FS_training_data <- data.frame(FS_training_data) # Turn it back to dataframe 
FS_training_data <-  add_column(FS_training_data,State=group_condition_training) # Adding State column which include samples types (Primary,Recurrent)
FS_training_data$State <- as.factor(FS_training_data$State) # Convert the State type to factor  

#remove near zero variation i.e same thing accross board (Tobi code)
all.zero <- preProcess(FS_training_data, method = 'nzv', uniqueCut = 15)
FS_training_data <- predict(all.zero, FS_training_data)
#remove highly correlated - to avoid colinearity (Tobi code)
all.corr <- preProcess(FS_training_data, method = 'corr', cutoff = 0.8)  
FS_training_data <- predict(all.corr, FS_training_data)


FS_testing_data <- test_data[significant_genes_names,]
FS_testing_data <- t(FS_testing_data)
FS_testing_data <- data.frame(FS_testing_data)
FS_testing_data <-  add_column(FS_testing_data,State=group_condition_testing) 

# setting seed for reproducibility 
set.seed(42)
# Building random forest classification model  
model <- randomForest(State~.,data = FS_training_data,proximity=TRUE,ntree=500,mtry=27)

# This commented code is for Hyper tuning of parameters ntree and mtry
#obb_error_df <- data.frame(
#  Trees= rep(1:nrow(model$err.rate),times=3),
#  Type= rep(c("OOB","Primary","Recurrent"),each=nrow(model$err.rate)),
#  Error=c(model$err.rate[,"OOB"],
#          model$err.rate[,"Primary"],
#          model$err.rate[,"Recurrent"])
#)

#ggplot(obb_error_df,aes(Trees,Error))+
#  geom_line(aes(color=Type))

#oob_values <- vector(length = 100)
#for (i in 1:100){
# tempmodel <- randomForest(State~.,data = FS_training_data,mtry=i,ntree=500)
#oob_values[i] <- tempmodel$err.rate[nrow(tempmodel$err.rate),1]
#}

#which(oob_values==min(oob_values))

predict(model,FS_testing_data)


#Evaluating model accuracy1
type_pred = predict(model, FS_testing_data)
FS_testing_data$type_pred = type_pred

#building confusion matrix
CFM = table(FS_testing_data$State, FS_testing_data$type_pred)
CFM

Classificaatio_Accuracy = sum(diag(CFM)/sum(CFM))
Classificaatio_Accuracy #= 0.8


conf_matrix <- confusionMatrix(as.factor(FS_testing_data$type_pred), as.factor(FS_testing_data$State))
conf_df <- as.data.frame(conf_matrix$table)


# Calculate the total number of instances
total <- sum(conf_df$Freq)

# Add a percentage column to the data frame
conf_df$Percentage <- (conf_df$Freq / total) * 100


# Visualization


ggplot(data = conf_df, aes(x = Prediction, y = Reference)) +
  geom_tile(aes(fill = Percentage), color = "white") +  # Use percentage for fill
  scale_fill_gradient(low = "#f4845c", high = "#1a122b", name = "Percentage") +  # Set color gradient for percentage
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = 1, color = "white", size = 5) +  # Display percentages inside tiles
  labs(title = "Confusion Matrix Heatmap (Percentage)",
       x = "Predicted Label",
       y = "Actual Label") +
  theme_minimal() +  # Clean theme
  theme(plot.title = element_text(hjust = 0.5, size = 16),  # Center title
        axis.text = element_text(size = 12))  # Resize axis labels



# Extract variable importance from the caret model
gini_importance <- varImp(model, scale = FALSE)  # rf.ll is your train model

# Convert to a data frame for plotting
gini_importance_df <- data.frame(Feature = rownames(gini_importance),
                                 GiniImportance = gini_importance$Overall)

# Sort by Gini importance for better visualization
gini_importance_df <- gini_importance_df[order(gini_importance_df$GiniImportance, decreasing = TRUE), ]

# Keep only the top 20 features based on Gini importance
gini_importance_top20 <- gini_importance_df[1:20, ]


# Create a lollipop plot for the top 20 features using ggplot2
ggplot(gini_importance_top20, aes(x = reorder(Feature, GiniImportance), y = GiniImportance)) +
  geom_point(color = "#f4845c", size = 4) +  # Create the lollipop "head"
  geom_segment(aes(xend = Feature, yend = 0), color = "#1a122b", size = 1.5) +  # Create the lollipop "stick"
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = "Top 20 Genes by Gini Importance (Random Forest)",
       x = "Genes",
       y = "Gini Importance (Mean Decrease)") +
  theme_minimal() +  # Apply a clean theme
  theme(plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
        axis.text.x = element_text(size = 12),  # Resize x-axis text
        axis.text.y = element_text(size = 12))  # Resize y-axis text

