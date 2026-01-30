# 1. LOAD LIBRARIES
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)

# 2. QUERY & DOWNLOAD DATA
# Specifically targeting RNA-Seq counts for TCGA-BRCA
query_TCGA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_TCGA)
data_se <- GDCprepare(query_TCGA)

# 3. DATA PROCESSING
# Extract FPKM-UQ or FPKM (study preferred unstranded FPKM)
exp_matrix <- assay(data_se, "fpkm_unstrand") 
rownames(exp_matrix) <- rowData(data_se)$gene_name

# Clinical data extraction
clinical <- as.data.frame(colData(data_se)) %>%
  select(barcode, vital_status, days_to_death, days_to_last_follow_up)

# Define Gene of Interest (Target)
target_gene <- "JUNB" 

# Prepare final dataframe
survival_data <- as.data.frame(t(exp_matrix[target_gene, , drop = FALSE]))
survival_data$barcode <- rownames(survival_data)

final_df <- merge(survival_data, clinical, by = "barcode") %>%
  mutate(
    time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time) & time > 0)

# Stratify into High/Low groups based on median expression
final_df$group <- factor(ifelse(final_df[[target_gene]] > median(final_df[[target_gene]]), "High", "Low"), 
                         levels = c("Low", "High"))

# 4. STATISTICAL CALCULATIONS
# Fit survival curve
fit <- survfit(Surv(time, status) ~ group, data = final_df)

# Calculate Cox Hazard Ratio (HR)
cox_fit <- coxph(Surv(time, status) ~ group, data = final_df)
hr_value <- round(summary(cox_fit)$conf.int[1], 2)
p_val_raw <- surv_pvalue(fit)$pval
p_val_formatted <- paste0("p = ", formatC(p_val_raw, format = "e", digits = 1))

# 5. FINAL VISUALIZATION (Publication Quality)
p <- ggsurvplot(
  fit, 
  data = final_df,
  palette = c("#377EB8", "#E41A1C"), # Blue (Low) and Red (High)
  size = 1,                 
  censor.shape = 124,       # Vertical bar ticks
  censor.size = 2,
  legend = "top",
  legend.title = target_gene, 
  legend.labs = c("Low (Ref)", "High (Risk)"),
  pval = p_val_formatted,   
  pval.coord = c(800, 0.25), 
  risk.table = TRUE,
  risk.table.height = 0.20,
  risk.table.title = "Number at risk",
  ggtheme = theme_classic(base_size = 14) 
)

# Add Hazard Ratio annotation
p$plot <- p$plot + 
  annotate("text", x = 800, y = 0.15, 
           label = paste0("HR = ", hr_value), 
           fontface = "bold", size = 5, hjust = 0) +
  labs(x = "Time (Days)", y = "Overall Survival Probability")

print(p)