# --- 1. SETTINGS & LIBRARIES ---
library(pheatmap)
library(RColorBrewer)

save_pheat <- function(x, filename, width=10, height=12) {
  # Opens a png device with high resolution (300 DPI)
  png(filename, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# --- 2. HALLMARK TGF-BETA HEATMAP ---
tgfb_genes_overlap <- intersect(tgfb_genes, rownames(expr_matrix))
heatmap_tgfb_data <- expr_matrix[tgfb_genes_overlap, ]
heatmap_tgfb_scaled <- t(scale(t(heatmap_tgfb_data)))

# Generate Plot
p_tgfb <- pheatmap(heatmap_tgfb_scaled, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   annotation_col = sample_info,
                   show_colnames = FALSE,
                   main = "Hallmark TGF-Beta Signaling: TNBC vs Non-TNBC",
                   color = colorRampPalette(c("green", "black", "firebrick3"))(100),
                   fontsize_row = 8,
                   border_color = NA,
                   silent = TRUE) # Don't show in RStudio yet, just capture


save_pheat(p_tgfb, "Heatmap_Hallmark_TGFB.png")

# --- 3. HALLMARK EMT HEATMAP ---
emt_plot_genes <- intersect(emt_genes, rownames(deg_sig))
heatmap_emt_data <- expr_matrix[emt_plot_genes, ]
heatmap_emt_scaled <- t(scale(t(heatmap_emt_data)))

p_emt <- pheatmap(heatmap_emt_scaled, 
                  cluster_rows = TRUE, 
                  cluster_cols = TRUE, 
                  annotation_col = sample_info,
                  show_colnames = FALSE,
                  main = "Hallmark EMT Pathway: TNBC vs Non-TNBC",
                  color = colorRampPalette(c("green", "black", "red"))(100),
                  fontsize_row = 5,
                  border_color = NA,
                  silent = TRUE)

save_pheat(p_emt, "Heatmap_Hallmark_EMT.png")