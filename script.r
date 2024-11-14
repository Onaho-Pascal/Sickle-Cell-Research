################################################
## Visualization of RNA-seq data of SCD Patients
################################################


install.packages("readxl")

library(readxl)
data <- read_excel("C:/Users/user/Documents/30-day Challenge/Geo rawcounts.xlsx", sheet = "GSE254951_geo_rawcounts")

View(data)

Scd_data <- as.data.frame(data)
View(Scd_data)

rownames(Scd_data) <- Scd_data$ensembl_id
Scd_data$ensembl_id <- NULL

# Filter out genes with low counts (e.g., counts < 10 in at least half of samples)
keep_genes <- rowSums(Scd_data >= 10) >= (ncol(Scd_data) / 2)
filtered_counts <- Scd_data[keep_genes, ]

View(filtered_counts)


library(gplots)
library(RColorBrewer)
div_color_palette <- rev(brewer.pal(11, "RdBu"))
seq_color_palette <- brewer.pal(9, "Blues")

boxplot(filtered_counts, xlab = "samples", ylab = "counts", las = 2, col = "lightblue") # las = 2 is to rotate the x-axis labels 

hist(filtered_counts[, "HU9"], 
     main = "distribution of Raw Counts for Sample 9 (Normal)", 
     xlab = "counts", 
     col = "lightgreen", 
     breaks = 50)

log_scd_data <- log2(filtered_counts + 1)
boxplot(log_scd_data, xlab = "samples", ylab = "counts", las = 2, col = "lightblue") # las = 2 is to rotate the x-axis labels 
hist(log_scd_data[, "HU9"], 
     main = "distribution of Raw Counts for Sample 9 (Log)", 
     xlab = "counts",
     col = "lightgreen", 
     breaks = 50)

z_scd_data <- t(scale(t(log_scd_data)))
boxplot(z_scd_data, xlab = "samples", ylab = "counts", las = 2, col = "blue")
hist(z_scd_data[, "HU9"], 
     main = "distribution of Raw Counts for Sample 9 (Zeta)", 
     xlab = "counts", 
     col = "purple", 
     breaks = 50)


heatmap.2(as.matrix(filtered_counts),
          col = seq_color_palette,
          Rowv = F, Colv = F, dendrogram = "none",
          sepcolor = "black",
          trace = "none",
          key = T,
          key.title = "Expression",
          density.info = "none",
          main = "Heatmap of top differentially expressed genes in Sickle Cell Disease",
          cexRow = 0.9, cexCol = 0.7, margins = c(11,10))
