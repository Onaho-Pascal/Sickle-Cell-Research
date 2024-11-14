# Evaluating Hydroxyurea’s Non-HbF-Related Impacts on Gene Expression and Pathways in Pediatric Sickle Cell Disease: An RNA-Seq Analysis
N.B: This study has been carried out before and the data was retrieved from the GEO database
This study specifically examined HC’s non-HbF (fetal hemoglobin) related impacts, providing insight into how HC affects genes and pathways critical to SCD pathophysiology in pediatric patients.  
Here’s a breakdown of my process so far:  
🧬 Data Exploration and Preprocessing  
* Dataset Loading and Preparation:
Imported the raw RNA-Seq counts data into R.
Assigned row names as gene identifiers for easier manipulation and analysis.
Filtered out genes with low expression (less than 10 counts in at least half the samples) to focus on genes with significant read counts.
* Data Normalization and Transformation:
Applied log2 transformation to stabilize variance, an essential step for RNA-Seq data visualization and downstream analysis.
Standardized the data further with Z-scores for sample comparability, making it ready for in-depth comparisons.
📊 Visualizations:
1) Boxplots:
Created boxplots to inspect count distribution across samples before and after normalization, helping identify any potential outliers or sample inconsistencies.

2) Histograms:
Generated histograms of raw counts and log-transformed counts for individual samples, giving a clear view of distribution shifts post-normalization.
3) Heatmap:
Visualized the expression of top genes through a heatmap, capturing gene expression patterns across all samples and highlighting potentially significant expression clusters affected by HC.
![Rplot26](https://github.com/user-attachments/assets/ff8c9977-7ffb-4f7b-ad79-b4dc7f6f59b2)
![Rplot27](https://github.com/user-attachments/assets/d9151abe-e89f-427e-b845-ac6910b196c6)
![Rplot28](https://github.com/user-attachments/assets/04778117-0d24-4858-bf6d-85bca08a75c5)
![Rplot25](https://github.com/user-attachments/assets/742ab741-72e5-4044-944a-2d5dcdd495f1)

