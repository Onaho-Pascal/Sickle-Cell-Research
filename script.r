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
