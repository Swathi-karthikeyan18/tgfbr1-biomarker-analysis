# 1. Load Libraries
library(affy)
library(limma)
library(sva)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(pheatmap)

# 2. Setup Annotation & Batch Info
# Based on your study description of 7 datasets
# Ensure your sample_info has a 'Batch' column corresponding to the 7 GSE IDs
sample_info <- data.frame(
  SampleType = factor(c(rep("Non_TNBC", 129), rep("TNBC", 544))),
  Batch = factor(c(rep("GSE95700", 57), rep("GSE58812", 107), 
                   rep("GSE83937", 131), rep("GSE157284", 82),
                   rep("GSE167213", 124), rep("GSE65194", 130), 
                   rep("GSE43358", 63))) # Adjusted to match your total 673
)
rownames(sample_info) <- colnames(expr_matrix)

# 3. Batch Correction for VISUALIZATION (PCA/Heatmap only)
# Study used removeBatchEffect for visualization purposes
expr_viz <- limma::removeBatchEffect(as.matrix(expr_matrix), 
                                     batch = sample_info$Batch,
                                     design = model.matrix(~SampleType, data=sample_info))

# 4. Differential Expression Analysis (The Proper Way)
# The study accounts for batch as a covariate within the linear model
design <- model.matrix(~ 0 + SampleType + Batch, data = sample_info)
colnames(design) <- make.names(colnames(design)) # Clean names for limma

# Rename columns for clarity in contrasts
colnames(design)[1:2] <- c("Non_TNBC", "TNBC")

fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(TNBC - Non_TNBC, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, adjust = "fdr", number = Inf)

# 5. TGF-beta Target Identification
tgfb_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>%
  pull(gene_symbol)

# 1. Stricter Thresholds (As per your study text)
logFC_threshold <- 1.5
fdr_threshold <- 0.05

# 2. Identify DEGs with the new |log2FC| > 1.5 filter
deg_results$diffexpressed <- "NO"
deg_results$diffexpressed[deg_results$logFC > logFC_threshold & deg_results$adj.P.Val < fdr_threshold] <- "UP"
deg_results$diffexpressed[deg_results$logFC < -logFC_threshold & deg_results$adj.P.Val < fdr_threshold] <- "DOWN"

deg_sig <- deg_results[deg_results$diffexpressed != "NO", ]

# 3. Extract TGF-beta Pathway Genes for Visualisation
# Using the Hallmark set as specified
tgfb_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>%
  pull(gene_symbol)

# Get the overlap between your DEGs and the TGF-beta pathway
overlap_genes <- intersect(rownames(deg_sig), tgfb_hallmark)

# 4. Prepare for DAVID (Export Gene Symbols)
# DAVID usually requires a simple list of gene symbols or Entrez IDs
write.table(rownames(deg_sig), "DEGs_for_DAVID.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 5. Prepare for GSEA (v4.3.2)
# GSEA desktop requires a .rnk file or a full expression matrix (.gct) and phenotype (.cls)
# Creating a ranked list based on Signal-to-Noise (or t-statistic as a proxy)
gsea_rank <- data.frame(Gene = rownames(deg_results), stat = deg_results$t)
gsea_rank <- gsea_rank[order(-gsea_rank$stat),]
write.table(gsea_rank, "TNBC_vs_NonTNBC.rnk", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

tgfb_results <- deg_results[rownames(deg_results) %in% tgfb_genes, ]