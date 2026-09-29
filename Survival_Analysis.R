# =============================================================================
# TCGA-BRCA overall-survival analysis  (TGF-betaR1 and other TGF-beta genes)
# =============================================================================

getwd()
setwd("D:\\TNBC\\Code_correction")

# 1. LIBRARIES -----------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)

GENES        <- c("TGFBR1", "TGFB1", "SMAD3", "JUNB")  # Figure 2 A, B, C, D
PRIMARY_GENE <- "TGFBR1"
ASSAY_NAME   <- Sys.getenv("ASSAY_NAME", unset = "fpkm_unstrand")  # see deviation note

# 2. QUERY & DOWNLOAD (cached) --------------------------------------------------
CACHE <- "data_se.rds"
if (file.exists(CACHE)) {
  message("Loading cached SummarizedExperiment from ", CACHE)
  data_se <- readRDS(CACHE)
} else {
  query_TCGA <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query_TCGA)          # ~GBs, one-off
  data_se <- GDCprepare(query_TCGA)
  saveRDS(data_se, CACHE)
}

# 3. EXPRESSION MATRIX (collapse duplicate gene symbols) ------------------------
assay_name <- if (ASSAY_NAME %in% assayNames(data_se)) ASSAY_NAME else assayNames(data_se)[1]
message("Using assay: ", assay_name)
exp_matrix <- assay(data_se, assay_name)

symbols <- rowData(data_se)$gene_name
keep <- !is.na(symbols)
exp_matrix <- exp_matrix[keep, , drop = FALSE]; symbols <- symbols[keep]
row_order  <- order(rowMeans(exp_matrix), decreasing = TRUE)          # highest first
exp_matrix <- exp_matrix[row_order, , drop = FALSE]; symbols <- symbols[row_order]
dup <- duplicated(symbols)
exp_matrix <- exp_matrix[!dup, , drop = FALSE]
rownames(exp_matrix) <- symbols[!dup]

# 4. CLINICAL TABLE (survival endpoint + covariates) ---------------------------
cd <- as.data.frame(colData(data_se))
# pick covariate columns defensively (names vary slightly across GDC releases)
pick <- function(df, cands) { hit <- cands[cands %in% names(df)]; if (length(hit)) hit[1] else NA }
age_col   <- pick(cd, c("age_at_index", "age_at_diagnosis", "age_at_initial_pathologic_diagnosis"))
stage_col <- pick(cd, c("ajcc_pathologic_stage", "tumor_stage", "ajcc_pathologic_tumor_stage"))
node_col  <- pick(cd, c("ajcc_pathologic_n", "pathologic_n"))

clinical <- data.frame(
  barcode      = cd$barcode,
  vital_status = cd$vital_status,
  days_to_death          = cd$days_to_death,
  days_to_last_follow_up = cd$days_to_last_follow_up,
  age   = if (!is.na(age_col))   suppressWarnings(as.numeric(as.character(cd[[age_col]]))) else NA,
  stage = if (!is.na(stage_col)) as.character(cd[[stage_col]]) else NA,
  node  = if (!is.na(node_col))  as.character(cd[[node_col]])  else NA,
  stringsAsFactors = FALSE
)
# tidy covariates: collapse stage to I/II/III/IV, node to N0 vs N+
clinical$stage <- gsub("[ABC]$", "", sub("^Stage ", "", clinical$stage))
clinical$stage <- factor(clinical$stage, levels = c("I", "II", "III", "IV"))
clinical$node  <- ifelse(grepl("^N0", clinical$node), "N0",
                         ifelse(is.na(clinical$node) | clinical$node %in% c("NX", ""), NA, "N+"))
clinical$node  <- factor(clinical$node, levels = c("N0", "N+"))

# 5. PER-GENE KM ANALYSIS (Figure 2 A-D) ---------------------------------------
# 5. PER-GENE KM ANALYSIS (Figure 2 A-D) ---------------------------------------
run_km <- function(gene) {
  if (!gene %in% rownames(exp_matrix)) { warning("Gene ", gene, " absent; skipping."); return(NULL) }
  
  df <- data.frame(barcode = colnames(exp_matrix), expr = exp_matrix[gene, ])
  df <- merge(df, clinical, by = "barcode") |>
    dplyr::mutate(time   = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
                  status = ifelse(vital_status == "Dead", 1L, 0L)) |>
    dplyr::filter(!is.na(time) & time > 0)
  
  df$group <- factor(ifelse(df$expr > median(df$expr), "High", "Low"), levels = c("Low", "High"))
  
  fit <- survfit(Surv(time, status) ~ group, data = df)
  cox <- coxph(Surv(time, status) ~ group, data = df)
  
  hr <- summary(cox)$conf.int[1]
  ci <- summary(cox)$conf.int[3:4]
  pv <- survminer::surv_pvalue(fit, data = df)$pval
  
  # Manuscript & Publication-ready statistical formatting
  pv_txt <- ifelse(pv < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", pv)))
  hr_txt <- paste0("HR = ", sprintf("%.2f", hr), " (95% CI: ", sprintf("%.2f", ci[1]), "–", sprintf("%.2f", ci[2]), ")")
  
  # Generate Kaplan-Meier plot object
  g <- ggsurvplot(
    fit, data = df,
    palette = c("#377EB8", "#E41A1C"), 
    size = 1,
    censor.shape = 124, 
    censor.size = 2, 
    legend = "top",
    legend.title = gene, 
    legend.labs = c("Low (Ref)", "High (Risk)"),
    pval = FALSE,                        # Disabled built-in pval to avoid curve collision
    risk.table = TRUE, 
    risk.table.height = 0.22,
    risk.table.y.text = FALSE,           # Clean, non-redundant risk table labels
    ggtheme = theme_classic(base_size = 13)
  )
  
  # Place statistics cleanly in the top-right open area (x = 5200)
  g$plot <- g$plot +
    annotate("text", x = 5200, y = 0.90, label = pv_txt, size = 4.2, hjust = 0) +
    annotate("text", x = 5200, y = 0.82, label = hr_txt, fontface = "bold", size = 4.2, hjust = 0) +
    labs(x = "Time (Days)", y = "Overall Survival Probability")
  
  # 1. Render in RStudio Plots tab
  print(g)
  
  # 2. Save grid object safely via ggsave without graphics device locking
  res_plot <- survminer::arrange_ggsurvplots(list(g), print = FALSE, ncol = 1, nrow = 1)
  ggplot2::ggsave(
    filename = paste0("Survival_", gene, ".png"),
    plot = res_plot,
    width = 7,
    height = 7,
    dpi = 300
  )
  
  data.frame(
    gene = gene, 
    n = nrow(df), 
    HR = round(hr, 2), 
    CI_low = round(ci[1], 2), 
    CI_high = round(ci[2], 2), 
    logrank_p = pv
  )
}

results <- do.call(rbind, lapply(GENES, run_km))
print(results)
write.csv(results, "Survival_summary.csv", row.names = FALSE)

# 6. COVARIATE-ADJUSTED COX + SCHOENFELD (methods requirement) -----------------
adj_df <- data.frame(barcode = colnames(exp_matrix), expr = exp_matrix[PRIMARY_GENE, ]) |>
  merge(clinical, by = "barcode") |>
  dplyr::mutate(time   = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
                status = ifelse(vital_status == "Dead", 1L, 0L),
                group  = factor(ifelse(expr > median(expr), "High", "Low"), levels = c("Low","High"))) |>
  dplyr::filter(!is.na(time) & time > 0)

# only include covariates that are actually present/non-empty
covars <- c("group",
            if (any(!is.na(adj_df$age)))   "age"   else NULL,
            if (any(!is.na(adj_df$stage))) "stage" else NULL,
            if (any(!is.na(adj_df$node)))  "node"  else NULL)
fml <- as.formula(paste("Surv(time, status) ~", paste(covars, collapse = " + ")))
cox_adj <- coxph(fml, data = adj_df)
cat("\n=== Covariate-adjusted Cox model for ", PRIMARY_GENE, " ===\n", sep = "")
print(summary(cox_adj))

zph <- survival::cox.zph(cox_adj)   # Schoenfeld residual test of PH assumption
cat("\n=== Schoenfeld residual (proportional-hazards) test ===\n")
print(zph)
message("PH assumption ", if (zph$table["GLOBAL", "p"] > 0.05) "holds" else "VIOLATED",
        " (global Schoenfeld p = ", signif(zph$table["GLOBAL", "p"], 3), ")")

message("Done. Wrote Survival_{", paste(GENES, collapse=","), "}.png + Survival_summary.csv")

# Visualise TGFBR1 (or any gene) directly in the RStudio Plots tab
PRIMARY_GENE_KM <- run_km(PRIMARY_GENE)

