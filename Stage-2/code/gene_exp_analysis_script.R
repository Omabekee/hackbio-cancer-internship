# Calling libraries 
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(biomaRt)

# STEP 1: Read in data and assign to variable

gbm_data <- read.csv("glioblastoma.csv", row.names=1)
colnames(gbm_data) <- sapply(strsplit(colnames(gbm_data),"\\."), function(x) paste(x[c(1,3,6)],collapse = "."))


# STEP 2: Generate heatmaps and visualize clusters

# a- Define a color palette for heatmap
# Diverging color palettes

colors_diverging <- colorRampPalette(c("Blue", "White", "Brown"))(100)
heatmap.2(as.matrix(gbm_data), trace = 'none',
          scale='row', dendogram = 'col',
          Colv = TRUE, Rowv = FALSE,
          cexCol = 0.6,
          col = colors_diverging,
          main = "Diverging color scale")

# Sequential color palettes
heatmap.2(as.matrix(gbm_data), trace = 'none',
          scale='row', dendogram = 'col',
          Colv = TRUE, Rowv = FALSE,
          cexCol = 0.6,
          col = colorRampPalette(brewer.pal(11, "Blues"))(100),
          main = "Sequential color scale")

# b- Cluster by genes, samples and both
# Cluster by both genes and samples
heatmap.2(as.matrix(gbm_data), trace = 'none',
          scale='row', dendrogram = 'both',
          Colv = TRUE, Rowv = TRUE,
          cexCol = 0.6,
          col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
          main = "Clustering by genes and samples")

# Cluster by rows (genes)
heatmap.2(as.matrix(gbm_data), trace = 'none',
          scale='row', dendrogram = 'row',
          Colv = FALSE, Rowv = TRUE,
          cexCol = 0.6,
          col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
          main = "Clustering by genes")

# Cluster by columns (samples)
heatmap.2(as.matrix(gbm_data), trace = 'none',
          scale='row', dendrogram = 'col',
          Colv = TRUE, Rowv = FALSE,
          cexCol = 0.6,
          col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
          main = "Clustering by samples")


# STEP 3: Subset genes that are significantly upregulated and downregulated

# a- First, calculate fold change and pvalue

# Divide the groups based on the clusters from the heatmap
# Select groups by index positions, group A and group B
groupA <- c(1,2,3,4,5)
groupB <- c(6,7,8,9,10)

groupA_data <- gbm_data[, groupA]
groupB_data <- gbm_data[, groupB]

# Get means of both groups to calculate fold change
groupA_mean <- rowMeans(groupA_data)
groupB_mean <- rowMeans(groupB_data)

# Calculate log fold change
logFC <- log2(groupB_mean+0.5) - log2(groupA_mean+0.5)

# Calculate p value
pvalues <- apply(gbm_data, 1, function(row) {
  wilcox.test(row[1:5], row[6:10])$p.value
})
# Create a dataframe that has the logFC and pvalue for each gene
DEG <- data.frame(logFC, pvalues)
rownames(DEG) <- rownames(gbm_data)

# Drawing volcano plot of Glioblastoma DEGs 
with(DEG, plot(logFC, -log10(pvalues), pch=19, main="Volcano plot of Glioblastoma DEGs"))
with(subset(DEG, pvalues<=0.05 & (logFC)>=1.5), points(logFC, -log10(pvalues), pch=19, col="red"))
with(subset(DEG, pvalues<=0.05 & (logFC)<= -1.5), points(logFC, -log10(pvalues), pch=19, col="blue"))
legend(x=-9,y=2.3,c("upregulated","downgulated"), cex=.8, bty="n", col=c("red","blue"),pch=19)


# b- Subset upregulated and downregulated genes
upreg_genes <- DEG %>% filter(logFC >= 1.5 & pvalues<= 0.05)

downreg_genes <- DEG %>% filter(logFC <= -1.5 & pvalues<= 0.05)

# writing the up and downregulated genes 
write.csv(upreg_genes,"upregulated_genes.csv",row.names=TRUE)  
write.csv(downreg_genes,"downregulated_genes.csv",row.names=TRUE)

# c- Get gene symbol from Ensembl
# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(DEG)

gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",                     
  values = ensembl_ids,                           
  mart = ensembl                                    
)

# Create a list of gene annotations
keys <- gene_annotations$ensembl_gene_id
values <- gene_annotations$hgnc_symbol

l <- list()
for (i in 1:length(keys)){
  l[keys[i]] <- values[i]
}

# add new column "symbol" to upreg_genes_with_symbol
upreg_genes_with_symbol <- upreg_genes
upreg_genes_with_symbol$symbol <- unlist(l[rownames(upreg_genes_with_symbol)], use.names = FALSE)

# add new column "symbol" to downreg_genes_with_symbol
downreg_genes_with_symbol <- downreg_genes
downreg_genes_with_symbol$symbol <- unlist(l[rownames(downreg_genes_with_symbol)], use.names = FALSE)

write.csv(upreg_genes_with_symbol, "upregu_genes_with_symb.csv", row.names=TRUE)
write.csv(downreg_genes_with_symbol, "downreg_genes_with_symb.csv", row.names=TRUE)


# STEP 4: Pathway visualization 
upreg_pathways <- read.csv("UpReg_Top_Pathways.csv") %>%  # Import upregulated enriched pathways
  add_column(Pathways_names=trimws(sub("GO:........","",upreg_pathways$Pathway)))    # getting the name of the pathways only in new column 


# Lollipop plot of the top 5 upregulated enriched pathways color scaling the points according to - log 10 of the p_val
ggplot(upreg_pathways[c(1:5),], aes(x = Pathways_names, y = nGenes)) + # defining the data and x,y axes
  geom_segment(aes(x = Pathways_names, xend = Pathways_names,
                   y = 0, yend = nGenes),
               color="black") +  # drawing segments of the lollipop in black color
  geom_point(aes(color=-log10(Enrichment.FDR)),
             size=8) + # Drawing the points of the lollipop and color scale it by - log 10 of the p_val
  scale_color_gradient(low = "blue",high = "red")+ # Specify the color scaling from blue to red
  coord_flip()+ # Flipping the axes for better visualization 
  theme_grey()+ # Using grey theme 
  labs(title = "Top 5 upregulated pathways: Gene Number and Significance Levels",
       x = "pathway",
       y = "gene_num",
       color="-log10 p_val")+    # specifying the labels for each element in the plot
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) # Adjusting text sizes in the plot

