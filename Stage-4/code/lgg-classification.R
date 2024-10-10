## =============================================================

#             DATASET EXTRACTION FROM TCGA

# ==============================================================


# Install required packages if not already installed
BiocManager::install("TCGAbiolinks", force = TRUE)
BiocManager::install("SummarizedExperiment", force = TRUE)

# Load necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# Step 1: Query the Gene Expression Data for LGG (Gliomas)
query_exp <- GDCquery(
  project = "TCGA-LGG",                          # Lower Grade Gliomas project
  data.category = "Transcriptome Profiling",     # Transcriptome Profiling category
  data.type = "Gene Expression Quantification",  # Expression Quantification type
  workflow.type = "STAR - Counts"                 # Choose workflow type (FPKM, counts, etc.)
)

# Step 2: Download the data
GDCdownload(query_exp)

# Step 3: Prepare the expression data
exp_data <- GDCprepare(query_exp)

# Step 4: Check the first few rows of the expression data
head(assay(exp_data))  # View the expression matrix (gene expression levels)

# Step 5: Query for Clinical Data to Obtain IDH Status
query_clin <- GDCquery(
  project = "TCGA-LGG",                          # Same project to align with expression data
  data.category = "Clinical",                    # Clinical data category
  data.type = "Clinical Supplement",             # Clinical data type
  data.format = "BCR XML"                        # XML format for detailed clinical info
)

# Step 6: Download and Prepare the Clinical Data
GDCdownload(query_clin)
clin_data <- GDCprepare_clinic(query_clin, clinical.info = "patient")  # Prepare clinical data

# Step 7: Create a column for IDH status in the clinical data
clin_data$idh_status <- ifelse(clin_data$ldh1_mutation_found == "Yes", "Mutant", "Wildtype")

# Select only patient barcodes and IDH status
idh_status <- clin_data[, c("bcr_patient_barcode", "idh_status")]

# Step 8: Merge Expression Data with IDH Status
# Extract colData as a data frame from expression data
exp_data_col <- as.data.frame(colData(exp_data))

# Add the barcode to the expression column data
exp_data_col$barcode <- rownames(exp_data_col)

# Merge based on the barcode
merged_data <- merge(exp_data_col, idh_status, by.x = "barcode", by.y = "bcr_patient_barcode", all.x = TRUE)

# Check for first few rows of the merged dataset
head(merged_data)

# save the expression data to a CSV file:
write.csv(assay(exp_data), "expression_data.csv", row.names = TRUE)

# Check the structure of merged_data to identify list columns
str(merged_data)

# Convert list columns to character format
merged_data[] <- lapply(merged_data, function(x) {
  if (is.list(x)) {
    return(sapply(x, function(y) if (is.null(y)) NA else paste(y, collapse = ";")))  # Convert lists to strings
  } else {
    return(x)  # Keep other types unchanged
  }
})

# Write to CSV
write.csv(merged_data, "merged_expression_clinical_data.csv", row.names = TRUE)


## =============================================================

#             DATA PREPROCESSING

# ==============================================================


#Check for zeroes
any(is.na(selected_data)) #returns FALSE meaning there is no N/A in the sample
#Check for zeroes
sum(selected_data == 0) #returns 13760483

# Normalize the data using gene length normalization
normData <- TCGAanalyze_Normalization(tabDF = selected_data,
                                      geneInfo = geneInfoHT,  # Check if this annotation is right for your data
                                      method = "geneLength")

# Next, to correct for sequencing depth, apply another normalization method-
# upper quartile normalization.
# Convert normData to a matrix for between-lane normalization
normDataMatrix <- as.matrix(normData)

# Apply upper quartile normalization
normDataDepthAdjusted <- EDASeq::betweenLaneNormalization(normDataMatrix, which = "upper")

# Convert back to a data frame if needed
normDataDepthAdjusted <- as.data.frame(normDataDepthAdjusted)

# Filter out low counts using quantile filtering
filteredData <- TCGAanalyze_Filtering(tabDF = normDataDepthAdjusted, 
                                      method = "quantile", 
                                      qnt.cut = 0.25)
#reduced to 34528 genes
#check number of zeroes present
sum(filteredData == 0)      # returns 2229859

#Remove rows with at least one zero
filtered_out_zeros <- filteredData[rowSums(filteredData == 0) == 0, ]
#reduced to 17100 genes

#save filtered data containing zeroes
write.csv(filteredData, file = "filteredData_with_zeroes.csv", row.names = TRUE)

# Write the preprocessed data to a CSV file
write.csv(filtered_out_zeros, file = "preprocessed_exp_data.csv", row.names = TRUE)

#Get the barcodes for the wild type from the metadata
wt_barcodes <- filtered_IDM$barcode[filtered_IDM$IDH_status == "WT"]
#Get the barcodes for the mutant type
mutant_barcodes <- filtered_IDM$barcode[filtered_IDM$IDH_status == "Mutant"]



## =============================================================

#             DIFFERENTIAL GENE EXPRESSION ANALYSIS

# ==============================================================


#Run DEA using TCGAbiolinks
DEresults <- TCGAanalyze_DEA(
  mat1 = filtered_out_zeros[, wt_barcodes], # wild type samples matrix
  mat2 = filtered_out_zeros[, mutant_barcodes], # mutant samples matrix
  Cond1type = "WT", # Condition 1: wild type
  Cond2type = "Mutant", # Condition 2: Mutant
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT"  # Method to use, "glmLRT" is recommended for RNA-Seq
)

# Write DEresults to CSV
write.csv(DEresults, file = "DEresults.csv", row.names = TRUE)

#DE analysis with IDH status
DEresults.level <- TCGAanalyze_LevelTab(DEresults, "WT", "Mutant",
                                        filtered_out_zeros[, wt_barcodes],
                                        filtered_out_zeros[, mutant_barcodes])

#Run a simple volcano plot
plot(DEresults$logFC, -log10(DEresults$FDR))

#Get upregulated and downregulated genes
upreg.genes <- rownames(subset(DEresults.level, logFC > 1)) # 412 genes
dnreg.genes <- rownames(subset(DEresults.level, logFC < 1)) # 1417 genes


# Get gene IDS and hgnc symbol from Ensembl using biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = upreg.genes,
                     mart = mart)$hgnc_symbol

dnreg.genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = dnreg.genes,
                     mart = mart)$hgnc_symbol


#Perform enrichment analysis for both
up.EA <- TCGAanalyze_EAcomplete(TFname="Upregulated", upreg.genes) 
dn.EA <- TCGAanalyze_EAcomplete(TFname="Downregulated", dnreg.genes) 



#Visualize the results for upreg and downreg
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP),
                        GOBPTab = up.EA$ResBP,
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat,
                        nRGTab = upreg.genes,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP),
                        GOBPTab = dn.EA$ResBP,
                        GOCCTab = dn.EA$ResCC,
                        GOMFTab = dn.EA$ResMF,
                        PathTab = dn.EA$ResPat,
                        nRGTab = dnreg.genes,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Enchanced Volcano plot
#Visualisation
res <- DEresults.level

library(EnhancedVolcano)

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




## =============================================================

#             K-NEAREST NEIGHBOR PREDICTION OF
#             LOW-GRADE GLIOBLASTOMA BIOMARKERS BY IDH STATUS 

# ==============================================================


# Step 1: Load necessary packages
library(dplyr)

# Step 2: Read the CSV file
clinical_data <- read.csv("merged_expression_clinical_data.csv", stringsAsFactors = FALSE)

# Step 3: Create a new dataframe with 'sample' and 'paper_IDH.status' columns
metadata <- clinical_data %>%
  select(sample, paper_IDH.status)

# Step 4: Rename 'paper_IDH.status' to 'IDH-status'
metadata <- metadata %>%
  rename(IDH.status = paper_IDH.status)

# Step 5: View the new dataframe
head(metadata)


lgg.data <- read.csv("preprocessed_exp_data.csv", row.names=1)


##=======MACHINE LEARNING STEP==========
#predicting IDH status as either wild type or mutant

install.packages("caret")
install.packages("DALEX")
install.packages("pROC")

library(caret)
library(DALEX)
library(pROC)
set.seed(34567)


#Unsupervised ML=KNN
table(metadata$IDH.status)  #419 mutant and 94 WT


#preview filtered normalized data
boxplot(lgg.data[.1:50], las=2)
par(oma = c(10,0,0,0))  #to visualize better
boxplot(log10(lgg.data[,1:100]+1),
        ylim = c(0,10),
        las = 2)   #log transform  

# Rename the sample columns so that it matches with meta data
colnames(lgg.data) <- gsub("(TCGA)\\.(\\w{2})\\.(\\w{4})\\.(\\d{2}[A-Z]).*", "\\1-\\2-\\3-\\4", colnames(lgg.data))


# Check the updated column names
colnames(lgg.data)



##========DATA PREPROCESSING========
# select top variable genes

top.var.genes <- data.frame(t(lgg.data))  # transpose the data
SDs = apply(top.var.genes, 2, sd)
topPreds = order(SDs, decreasing = T)[1:1000]
top.var.genes = top.var.genes[, topPreds]
dim(top.var.genes)  # after selecitng top 1000 genes


# adding meta data
top.var.genes[1:5,1:5]

# confirm if the row names are same
head(metadata)  
head(top.var.genes)  # not the same, metadata has a different row name from top.var.genes


# Reset row names
rownames(metadata) <- NULL

# Set the 'sample' column as row names
library(tibble)
new.metadata <- metadata %>% 
  column_to_rownames(var = "sample")

# Check the new metadata
head(new.metadata)

# Merge
merged.data <- merge(top.var.genes, new.metadata, by = "row.names")

dim(merged.data)

rownames(merged.data) <- merged.data$Row.names #make the sample ids row names again
merged.data[1:4,1:4]

merged.data <- merged.data[,-1]
merged.data[1:4,1:4]



# Further preprocessing: remove highly correlated genes
# 1. Remove near zero variation
all.zero <- preProcess(merged.data, method = 'nzv', uniqueCut = 15)
merged.data <- predict(all.zero, merged.data)

#2. center -I skipped this step because I noticed that it
# negated some of the values,
# all.center <- preProcess(lgg.exp, method = 'center')
# lgg.exp <- predict(all.center, lgg.exp)

#3. remove highly correlated genes to avoid collinearity
all.corr <- preProcess(merged.data, method = 'corr', cutoff = 0.5)
merged.data <- predict(all.corr, merged.data)   # resulted in 124 columns from 1000


###============= preprocessing done =================



###============= TRAIN THE MODEL ====================
# splitting daatset using the standard 70:30 split
intrain <- createDataPartition(y = merged.data$IDH.status, p =0.7)[[1]]
length(intrain)     #375

# separate the test and training
train.lgg <- merged.data[intrain,]
test.lgg <- merged.data[-intrain,]

dim(train.lgg)   #375 124  
dim(test.lgg)    #159 124  

train.lgg$IDH.status
test.lgg$IDH.status

# let's train
# remove missing values

sum(is.na(train.lgg))   #15

train.lgg <- na.omit(train.lgg)


ctrl.lgg <- trainControl(method = 'cv', number = 5)

install.packages("randomForest")
library(randomForest)
knn.lgg <- train(IDH.status~.,
                 data = train.lgg,
                 method = "knn",
                 trControl = ctrl.lgg,
                 tuneGrid = data.frame(k=1:20))

# the best k is:
knn.lgg$bestTune    # k = 1

# predict
trainPred <- predict(knn.lgg, newdata = train.lgg)
testPred <- predict(knn.lgg, newdata = test.lgg)


###========== TESTING AND INTERPREATION ===========

# Confusion Matrix

confusionMatrix(trainPred, train.lgg$IDH.status)
confusionMatrix(testPred, test.lgg$IDH.status)

# this brough out an error: 

# Check the levels of the predictions and the actual values
levels(trainPred)              # "Mutant" "WT"    
levels(train.lgg$IDH.status)   #  NULL

levels(testPred)               # "Mutant" "WT"    
levels(test.lgg$IDH.status)    #  NULL

# Convert IDH.status to factors in both train and test datasets
train.lgg$IDH.status <- factor(train.lgg$IDH.status, levels = c("Mutant", "WT"))
test.lgg$IDH.status <- factor(test.lgg$IDH.status, levels = c("Mutant", "WT"))
levels(train.lgg$IDH.status)    # "Mutant" "WT"  


# Rerun Confusion Matrix

# Confusion Matrix for Training Data
train_cm <- confusionMatrix(trainPred, train.lgg$IDH.status)  #Accuracy: 100%, 294 acc predicted for mutant and 66 for WT, 0 incorrectly predicted
print(train_cm)

# Confusion Matrix for Test Data
test_cm <- confusionMatrix(testPred, test.lgg$IDH.status)     #Accuracy: 91%, 123 acc predicted mutant, 2 incorrectly predicted; 17 accurately predicted WT, 11 inaccurately predicted
print(test_cm)

# Extracting the table from confusion matrix for visualization for training data
train_cm_table <- as.data.frame(train_cm$table)


# Visualization
ggplot(train_cm_table, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "lightcoral", high = "darkslateblue") + 
  geom_text(aes(label = Freq), vjust = 1, color = "white", size = 5) +  # to display frequencies
  theme_minimal() +
  labs(title = "KNN Confusion Matrix (Training Data)", 
       x = "Actual IDH Status", 
       y = "Predicted IDH Status")

# Visualization for test data
test_cm_table <- as.data.frame(test_cm$table)

# Visualization 
ggplot(test_cm_table, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "lightcoral", high = "darkslateblue") + 
  geom_text(aes(label = Freq), vjust = 1, color = "white", size = 5) +  # to display frequencies
  theme_minimal() +
  labs(title = "KNN Confusion Matrix (Test Data)", 
       x = "Actual IDH Status", 
       y = "Predicted IDH Status")


###========== VARIABLE IMPORTANCE ===========

# Determine variable importance

explainer.lgg <- explain(knn.lgg, label = 'knn',
                         data = train.lgg,
                         y = as.numeric(train.lgg$IDH.status))

importance.lgg <- feature_importance(explainer.lgg,
                                     n_sample=50,
                                     type = 'difference')

View(importance.lgg)
plot(importance.lgg)

head(importance.lgg$variable) 

tail(importance.lgg$variable)   #check the least expressed variables (genes)

# highly expressed genes = "ENSG00000113140.11", "ENSG00000133048.13" "ENSG00000109846.9" ,
#"ENSG00000198695.2"  "ENSG00000154978.13"

# least expressed genes = "ENSG00000161642.18" "ENSG00000185920.16", "ENSG00000198734.12"
# "ENSG00000136002.20"



## =============================================================

#             RANDOM FOREST CLASSIFICATION OF
#             LOW-GRADE GLIOBLASTOMA BIOMARKERS  

# ==============================================================


llg.data <- read.csv("preprocessed_exp_data.csv")
llg.data.metadata <- read.csv("merged_expression_clinical_data.csv")

rownames(llg.data) <- llg.data$X #make the sample ids row names
llg.data <- llg.data[,-1] #remove the extra column


#change the - separated names in the metadata 
llg.data.metadata$barcode1 <- gsub("-", ".", llg.data.metadata$barcode)


#install.packages("caret")
#install.packages("DALEX")
#install.packages("pROC")
#install.packages("ranger")
#install.packages("reshape2")

library(reshape2)
library("caret")
library(pROC)
library(DALEX)
library("ranger")
library(tidyverse)



set.seed(34567)

#Supervised machine learning


#transpose the data so features are column
llg.data <- data.frame(t(llg.data))


#create meta list
meta <- llg.data.metadata %>% select(barcode1, paper_IDH.status, )

rownames(meta) <- meta$barcode1 #make the sample ids row names
meta <- meta %>% select(-barcode1) #remove the barcode1


#how many sampls do we have under each subset
table(meta$paper_IDH.status) #- Mutant; 419     WT; 94

dim(llg.data) #[1]   534 16949

#select top variable genes
SDs = apply(llg.data, 2, sd)
topPreds = order(SDs, decreasing = T)[1:1000]
llg.data1 = llg.data[, topPreds]
dim(llg.data1) #534 1000


#ourFinal Dataset for ML
llg.data2 <- merge(llg.data1, meta, by = "row.names")
dim(llg.data2) #534 1000
ll.data[1:5,1:5]
rownames(llg.data2) <- llg.data2$Row.names #make the sample ids row names again
llg.data2[1:5,1:5]
llg.data2 <- llg.data2[,-1]
llg.data2[1:5,1:5]


#we need to perform some preprocessing
#1 remove near zero variation i.e same thing accross board
all.zero <- preProcess(llg.data2, method = 'nzv', uniqueCut = 15)
llg.data2 <- predict(all.zero, llg.data2)

##2. center - to have a mean of 1 and sd of 0? --- didnt do this
#all.center <- preProcess(llg.data2, method = 'center')
#llg.data2 <- predict(all.center, llg.data2)

#3. remove highly correlated - to avoid colinearity
all.corr <- preProcess(llg.data2, method = 'corr', cutoff = 0.5)
llg.data2 <- predict(all.corr, llg.data2)

dim(llg.data2) #534 124

#count NA in paper_IDH.status column
na_count <- sum(is.na(llg.data2$paper_IDH.status))

# Print the count
print(na_count) # = 21

# Remove rows where 'gene_expression' has NA
llg.data3 <- llg.data2 %>% filter(!is.na(paper_IDH.status))

dim(llg.data3) # =  513 124



#splitting the dataset (70:30)
intrain <- createDataPartition(y = llg.data3$paper_IDH.status, p = 0.7)[[1]]

#separate the test and training
train.lgg <- llg.data3[intrain,]
test.lgg <- llg.data3[-intrain,]

dim(train.lgg) # 360 124
dim(test.lgg) #  153 124

#let's train
#random forest

rf.ctrl <- trainControl(method = 'cv') #bootstrap can be used also

rf.ll <- train(paper_IDH.status~.,
               data = train.lgg,
               method = 'ranger', #random forest
               trControl = rf.ctrl,
               importance = 'permutation',
               tuneGrid = data.frame(mtry=100,
                                     min.node.size = 1,
                                     splitrule="gini"
               ))
rf.ll$finalModel$prediction.error # = 0.03

plot(varImp(rf.ll), top = 10)


#Evaluating model accuracy1
type_pred = predict(rf.ll, test.lgg)
test.lgg$type_pred = type_pred

#building confusion matrix
CFM = table(test.lgg$paper_IDH.status, test.lgg$type_pred)
CFM

Classificaatio_Accuracy = sum(diag(CFM)/sum(CFM))
Classificaatio_Accuracy #= 0.90


conf_matrix <- confusionMatrix(as.factor(test.lgg$type_pred), as.factor(test.lgg$paper_IDH.status))
conf_df <- as.data.frame(conf_matrix$table)


# Calculate the total number of instances
total <- sum(conf_df$Freq)

# Add a percentage column to the data frame
conf_df$Percentage <- (conf_df$Freq / total) * 100

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
gini_importance <- varImp(rf.ll, scale = FALSE)  # rf.ll is your train model

# Convert to a data frame for plotting
gini_importance_df <- data.frame(Feature = rownames(gini_importance$importance),
                                 GiniImportance = gini_importance$importance$Overall)

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

print(conf_matrix)

library(gridExtra)
library(grid)



# Assuming you've already created your confusion matrix
# conf_matrix <- confusionMatrix(predictions, actual)

# Convert confusion matrix to a data frame
conf_matrix_df <- as.data.frame(as.table(conf_matrix$table))

# Create tables for metrics
overall_metrics_df <- as.data.frame(conf_matrix$overall)
class_metrics_df <- as.data.frame(conf_matrix$byClass)

# Create table plots
conf_matrix_plot <- tableGrob(conf_matrix_df)
overall_metrics_plot <- tableGrob(overall_metrics_df)
class_metrics_plot <- tableGrob(class_metrics_df)

# Add titles using textGrob
conf_matrix_title <- textGrob("Confusion Matrix", gp = gpar(fontsize = 14, fontface = "bold"))
overall_metrics_title <- textGrob("Overall Performance Metrics", gp = gpar(fontsize = 14, fontface = "bold"))
class_metrics_title <- textGrob("Class-level Performance Metrics", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange the tables and titles in a grid
grid_plot <- grid.arrange(
  conf_matrix_title, conf_matrix_plot,
  overall_metrics_title, overall_metrics_plot,
  class_metrics_title, class_metrics_plot,
  ncol = 1, heights = c(0.5, 2, 0.5, 2, 0.5, 2)  # Adjust heights to make room for titles
)

# Save the grid as a JPEG file
ggsave("confusion_matrix_and_metrics_with_titles.jpg", plot = grid_plot, width = 8, height = 10, dpi = 300)
