# Load required packages
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)


# load in the data
counts <- read.delim("C:/Users/user/Documents/Bioinformatics Projects/Sickle Cell Project/GSE254951_geo_rawcounts.txt/GSE254951_geo_rawcounts.txt", row.names = 1)
metadata <- read.csv("C:/Users/user/Documents/Bioinformatics Projects/Sickle Cell Project/Metadata.csv", stringsAsFactors = TRUE)

metadata$Sample.ID <- NULL
rownames(metadata) <- metadata$Sample
metadata <- metadata[, -1]

# Check alignment

all(colnames(counts) %in% rownames(metadata))
counts <- counts[, rownames(metadata)]


# Create dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),  # DESeq2 expects integers
                              colData = metadata,
                              design = ~ Condition)

# Pre-filtering low counts
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "HU-MTD", "pre-HU"))
resultsNames(dds)

BiocManager::install("apeglm")
res <- lfcShrink(dds, coef = "Condition_pre.HU_vs_HU.MTD", type = "apeglm")


# Save as CSV

if (!dir.exists("results")) {
  dir.create("results")
}

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, "results/deseq2_results.csv", row.names = FALSE)

#Visualizations

#Volcano Plot
png("results/volcano_plot.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(res_df,
                lab = res_df$gene_id,
                x = "log2FoldChange",
                y = "pvalue",
                title = "Hydroxyurea Response: HU-MTD vs Pre-HU",
                pCutoff = 0.05,
                FCcutoff = 1.5)
dev.off()


# PCA PLOT

vsd <- vst(dds, blind = FALSE)

png("results/pca_plot.png", width = 1200, height = 1000, res = 150)
plotPCA(vsd, intgroup = "Condition")
dev.off()


#Functional Enrichment Analysis

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

up_genes <- readLines("C:/Users/user/Documents/results/upregulated_genes.txt")
down_genes <- readLines("C:/Users/user/Documents/results/downregulated_genes.txt")


# Remove version numbers
up_genes <- sub("\\..*", "", up_genes)

# Remove NAs and duplicates
up_genes <- unique(na.omit(up_genes))

# Convert Ensembl IDs to gene symbols
gene_map <- mapIds(org.Hs.eg.db,
                   keys = up_genes,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")

# Clean final list
up_gene_symbols <- na.omit(gene_map)


ego_up <- enrichGO(gene = up_gene_symbols,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

#kegg_up <- enrichKEGG(gene = up_gene_symbols,
                     # organism = "hsa",
                      #keyType = "kegg",
                      #pvalueCutoff = 0.05)

head(ego_up)
barplot(ego_up, showCategory = 10, title = "GO Biological Processes - Upregulated")


# Remove version numbers
down_genes <- sub("\\..*", "", down_genes)

# Remove NAs and duplicates
down_genes <- unique(na.omit(down_genes))

# Convert Ensembl IDs to gene symbols
gene_map <- mapIds(org.Hs.eg.db,
                   keys = down_genes,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")

# Clean final list
down_gene_symbols <- na.omit(gene_map)


ego_down <- enrichGO(gene = down_gene_symbols,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

#kegg_up <- enrichKEGG(gene = up_gene_symbols,
# organism = "hsa",
#keyType = "kegg",
#pvalueCutoff = 0.05)

head(ego_down)
barplot(ego_down, showCategory = 10, title = "GO Biological Processes - Downregulated")
