# Loading packages
library(airway)
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)


# Loading scale_rows function
scale_rows <- function (x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}


# Loading airway data
data(airway)


# Creating a variable with annotation data
anno <- as.data.frame(colData(airway))
anno <- data.frame(cond = anno$dex, 
                   row.names = rownames(anno))


# Creating a variable with count data
counts <- assay(airway)


# Creating a DESeq object
ddsHTSeq <- DESeqDataSetFromMatrix(countData = counts, 
                                   colData = anno, 
                                   design = ~ cond)


# Sorting the factors
colData(ddsHTSeq)$cond <- factor(colData(ddsHTSeq)$cond, levels=c("untrt", "trt"))


# Performing DE analysis
dds <- DESeq(ddsHTSeq)


# Normalizing the data
vsd <- varianceStabilizingTransformation(dds)


# Getting DE results
res <- results(dds, contrast=c("cond","untrt","trt"), alpha = 0.05)
summary(res)


# Plotting PCA
plotPCA(vsd, intgroup="cond")


# Selecting DE genes
selected_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange)>1.5))


# Heatmap
tab_heatmap <- assay(vsd)[selected_genes, ]
tab_heatmap <- as.data.frame(tab_heatmap)
tab_heatmap <- scale_rows(tab_heatmap)

tab_heatmap <- tab_heatmap[, rownames(anno)]

ha1 = HeatmapAnnotation(df = anno, col = list(cond = c("untrt" = "black", "trt" = "green")))

breaks <- seq(-2,2, by= 0.1)
ht <- Heatmap(tab_heatmap, 
              top_annotation = ha1, 
              name = "zscore", column_title = "", width = 1, 
              show_row_names = F, show_column_names = F,
              cluster_rows = T, cluster_columns = F,
              clustering_distance_columns = "pearson", clustering_method_columns = "ward.D2",
              clustering_distance_rows = "pearson", clustering_method_rows = "ward.D2", 
              show_column_dend = F, show_row_dend = T, 
              row_names_gp = gpar(fontsize = 2))
print(ht)

######################################################
#################### CHALLENGE!!! ####################
######################################################
# Sort anno dataframe based on conditions and plot the 
# heatmap with the grouped columns

######################################################
######################################################
######################################################



# Functional Enrichment
group <- groupGO(gene = selected_genes, 
                 OrgDb = org.Hs.eg.db, 
                 keyType = "ENSEMBL", 
                 ont = "BP", 
                 readable = T)
group_df <- as.data.frame(group)


enrich <- enrichGO(gene = selected_genes, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENSEMBL", 
                   ont = "BP", 
                   readable = T)
enrich_df <- as.data.frame(enrich)

barplot(enrich)

dotplot(enrich)

emapplot(enrich)

heatplot(enrich)




######################################################
#################### CHALLENGE!!! ####################
######################################################
# Do the same analysis as we did above using the counting 
# table located in /ref/counts.txt and taking into account only comparing 
# CellType in the table with the annotations located in /ref/SampleInfo.txt.
# Use org.Mm.eg.db package to perform the enrichment

######################################################
######################################################
######################################################