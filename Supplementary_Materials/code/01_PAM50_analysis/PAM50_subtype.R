# analysis_PAM50_subtype.R
# PAM50 subtype stratification + eRNA risk score independent prognostic analysis

# 0. Environment setup
cat(">>> [0/6] Initializing environment...\n")

setwd(".")

out_dir <- "Results/PAM50_Subtype"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pkgs <- c("data.table", "survival", "survminer", "ggplot2",
          "dplyr", "patchwork", "cowplot", "dunn.test", "ggrastr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Global color palette
col_high <- "#D62728"
col_low  <- "#1F77B4"

pam50_colors <- c(
  "Luminal A"     = "#E64B35",
  "Luminal B"     = "#4DBBD5",
  "HER2-enriched" = "#00A087",
  "Basal-like"    = "#3C5488",
  "Normal-like"   = "#F39B7F"
)

# 1. Load data
cat("\n>>> [1/6] Loading data files...\n")

# eRNA expression matrix
data <- fread("Data_Source/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv")

# Standardize column names to TCGA format
original_col <- colnames(data)
new_col <- ifelse(
  grepl("_tumor", original_col), gsub("_tumor", "-01A", original_col),
  ifelse(grepl("_normal", original_col), gsub("_normal", "-11A", original_col),
         original_col)
)
colnames(data) <- new_col

# Define tumor sample columns
tumor_col <- colnames(data)[which(substr(colnames(data), 14, 14) == 0)]
data <- as.data.frame(data)
rownames(data) <- data[, 1]
data <- data[, -1]
cat(sprintf("  eRNA matrix: %d eRNAs x %d samples\n", nrow(data), ncol(data)))

# Survival data (UTF-16 encoding)
survival <- read.delim("Data_Source/TCGA-BRCA.survival.tsv", fileEncoding = "UTF-16")
survival <- unique(survival[, c("X_PATIENT", "OS.time", "OS")])
survival <- as.data.frame(survival)
rownames(survival) <- survival[, 1]
survival <- survival[, -1]

# Convert survival time to months
survival$OS.time <- round(as.numeric(survival$OS.time) / 30, 2)
colnames(survival)[1] <- "OS.time/month"
survival$OS <- as.numeric(survival$OS)
cat(sprintf("  Survival data: %d patients\n", nrow(survival)))

# Clinical data
data_clinical <- fread("Data_Source/TCGA-BRCA.clinical.tsv")
cat(sprintf("  Clinical data: %d rows\n", nrow(data_clinical)))

# PAM50 subtype data
pam50_raw <- fread("Data_Source/TCGA-BRCA.clinicalMatrix.tsv",
                   select = c("sampleID", "PAM50Call_RNAseq"))
pam50_raw <- pam50_raw[!is.na(PAM50Call_RNAseq) &
                         PAM50Call_RNAseq != "null" &
                         PAM50Call_RNAseq != ""]
cat(sprintf("  PAM50 data: %d valid samples\n", nrow(pam50_raw)))

# 2. eRNA matrix preprocessing
cat("\n>>> [2/6] Preprocessing eRNA matrix...\n")

# Filter low-expression features
keep <- rowSums(data >= 0.3) >= 500
data.filtered <- data[keep, ]
cat(sprintf("  After filtering: %d eRNAs\n", nrow(data.filtered)))

# log2 transformation
data.filtered <- log2(data.filtered + 1)

# Sanitize row names
rownames(data.filtered) <- gsub(":", "_", rownames(data.filtered))

# Subset tumor samples and truncate to 12-character patient IDs
data.filtered.tumor <- data.filtered[, tumor_col]
colnames(data.filtered.tumor) <- substr(colnames(data.filtered.tumor), 1, 12)

# 3. Build analysis dataset
cat("\n>>> [3/6] Building analysis dataset...\n")

# 10-eRNA prognostic signature
selected_genes <- c("chr1_155158995", "chr3_11236700", "chr8_22624675",
                    "chr3_138070534", "chr9_71398719", "chr9_114689796",
                    "chr10_5531356",  "chr9_71398939", "chr12_13371038",
                    "chr10_5528926")

# Validate eRNA presence
missing_genes <- selected_genes[!selected_genes %in% rownames(data.filtered.tumor)]
if (length(missing_genes) > 0) {
  stop(sprintf("Missing eRNAs: %s", paste(missing_genes, collapse = ", ")))
} else {
  cat("  All 10 eRNAs found\n")
}

# Intersect with survival data
overlap_combined <- intersect(rownames(survival), colnames(data.filtered.tumor))
cat(sprintf("  Overlapping samples: %d\n", length(overlap_combined)))

combined_data <- cbind(
  survival[overlap_combined, ],
  t(data.filtered.tumor)[overlap_combined, selected_genes]
)
combined_data <- as.data.frame(combined_data)

# Clean column names
colnames(combined_data) <- gsub(":", "_", colnames(combined_data))
colnames(combined_data) <- gsub("/month", "_month", colnames(combined_data))

# 4. Compute risk scores
cat("\n>>> [4/6] Computing risk scores...\n")

combined_data$risk_score <-
  (-0.12338077) * combined_data[, "chr1_155158995"] +
  ( 0.21888531) * combined_data[, "chr3_11236700"]  +
  (-0.24173782) * combined_data[, "chr8_22624675"]  +
  (-0.85296049) * combined_data[, "chr3_138070534"] +
  (-0.20550119) * combined_data[, "chr9_71398719"]  +
  (-0.20987822) * combined_data[, "chr9_114689796"] +
  ( 0.21679753) * combined_data[, "chr10_5531356"]  +
  (-0.03772779) * combined_data[, "chr9_71398939"]  +
  ( 0.35946524) * combined_data[, "chr12_13371038"] +
  ( 0.06413942) * combined_data[, "chr10_5528926"]

median_risk <- median(combined_data$risk_score, na.rm = TRUE)
combined_data$risk_group <- ifelse(combined_data$risk_score > median_risk, "High", "Low")
combined_data$risk_group <- factor(combined_data$risk_group, levels = c("Low", "High"))
cat(sprintf("  Median risk score: %.4f\n", median_risk))
cat(sprintf("  High: %d, Low: %d\n",
            sum(combined_data$risk_group == "High"),
            sum(combined_data$risk_group == "Low")))

# 5. Integrate PAM50 data
cat("\n>>> [5/6] Integrating PAM50 data...\n")

pam50_df <- as.data.frame(pam50_raw)
pam50_df$patient_id <- substr(pam50_df$sampleID, 1, 12)

# Recode to standard subtype labels
pam50_df$PAM50_subtype <- dplyr::recode(pam50_df$PAM50Call_RNAseq,
                                        "LumA"   = "Luminal A",
                                        "LumB"   = "Luminal B",
                                        "Her2"   = "HER2-enriched",
                                        "Basal"  = "Basal-like",
                                        "Normal" = "Normal-like"
)
pam50_df$PAM50_subtype <- factor(pam50_df$PAM50_subtype,
                                 levels = names(pam50_colors))
pam50_clean <- pam50_df[!duplicated(pam50_df$patient_id),
                        c("patient_id", "PAM50_subtype")]

# Merge PAM50 into combined_data
combined_data$patient_id <- rownames(combined_data)
df <- merge(combined_data, pam50_clean, by = "patient_id", all.x = FALSE)
df <- df[!is.na(df$PAM50_subtype), ]
df <- df[!is.na(df$OS.time_month) & df$OS.time_month > 0, ]

# Apply 150-month administrative censoring
df$OS <- as.numeric(df$OS)
df$OS <- ifelse(df$OS.time_month >= 150, 0, df$OS)
df$OS.time_month <- ifelse(df$OS.time_month >= 150, 150, df$OS.time_month)

cat(sprintf("  Final analysis samples: %d\n", nrow(df)))
cat("  PAM50 subtype distribution:\n")
print(table(df$PAM50_subtype))
cat("  Risk group x PAM50:\n")
print(table(df$risk_group, df$PAM50_subtype))

# 6. Generate figures
cat("\n>>> [6/6] Generating figures...\n")

# Figure 1: KM curves stratified within each PAM50 subtype
# Risk groups defined by subtype-specific median risk score
cat("  [6-1] KM curves (subtype-specific median cutoff)\n")

subtypes <- levels(df$PAM50_subtype)
km_plots <- list()

for (st in subtypes) {

  df_sub <- df[df$PAM50_subtype == st, ]
  n_total <- nrow(df_sub)

  # Subtype-specific median cutoff for balanced High/Low split
  median_st <- median(df_sub$risk_score, na.rm = TRUE)
  df_sub$risk_group_st <- factor(
    ifelse(df_sub$risk_score > median_st, "High", "Low"),
    levels = c("Low", "High")
  )

  n_high <- sum(df_sub$risk_group_st == "High")
  n_low  <- sum(df_sub$risk_group_st == "Low")

  cat(sprintf("    %s: n=%d, median_score=%.4f, High=%d, Low=%d\n",
              st, n_total, median_st, n_high, n_low))

  if (n_total < 20) {
    km_plots[[st]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = sprintf("%s\nn=%d\n(insufficient samples)", st, n_total),
               size = 4, color = "grey50") +
      theme_void()
    next
  }

  fit_sub <- survfit(Surv(OS.time_month, OS) ~ risk_group_st, data = df_sub)

  p_km <- ggsurvplot(
    fit_sub,
    data              = df_sub,
    pval              = TRUE,
    pval.method       = TRUE,
    pval.size         = 4,
    pval.method.size  = 4,
    pval.coord        = c(10, 0.18),
    pval.method.coord = c(10, 0.25),
    conf.int          = FALSE,
    palette           = c(col_low, col_high),
    risk.table        = FALSE,
    xlim              = c(0, 150),
    break.time.by     = 30,
    legend.title      = "Risk Group",
    legend.labs       = c(sprintf("Low (n=%d)",  n_low),
                          sprintf("High (n=%d)", n_high)),
    legend            = c(0.80, 0.90),
    title             = st,
    xlab              = "Time (Months)",
    ylab              = "Survival Probability",
    font.title        = c(12, "bold"),
    font.x            = c(11, "bold"),
    font.y            = c(11, "bold"),
    font.tickslab     = c(10, "plain"),
    font.legend       = c(10, "plain"),
    ggtheme = theme_classic2() +
      theme(
        panel.border      = element_rect(colour = "black", fill = NA),
        axis.line         = element_line(linewidth = 0.8),
        legend.title      = element_text(size = 10, face = "bold"),
        legend.text       = element_text(color = "black"),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        plot.title        = element_text(
          hjust = 0.5,
          face  = "bold",
          color = pam50_colors[st],
          size  = 12
        )
      )
  )
  km_plots[[st]] <- p_km$plot
}

# Arrange KM plots in 2x3 grid
km_combined <-
  (km_plots[["Luminal A"]]  | km_plots[["Luminal B"]]  | km_plots[["HER2-enriched"]]) /
  (km_plots[["Basal-like"]] | km_plots[["Normal-like"]] | plot_spacer()) +
  plot_annotation(
    title    = "Kaplan-Meier Survival Curves Stratified by PAM50 Subtype",
    subtitle = "High vs. Low eRNA risk group within each PAM50 subtype (cutoff: subtype-specific median risk score)",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

out_km <- file.path(out_dir, "Fig_KM_by_PAM50_subtype.svg")
ggsave(out_km, km_combined, width = 16, height = 12, device = "svg")
cat(sprintf("  Saved: %s\n", out_km))

# Figure 2: Risk score boxplot across PAM50 subtypes
cat("  [6-2] Risk score boxplot\n")

# Kruskal-Wallis global test
kw_test  <- kruskal.test(risk_score ~ PAM50_subtype, data = df)
kw_label <- sprintf("Kruskal-Wallis: p = %.2e", kw_test$p.value)

# Dunn post-hoc test with Bonferroni correction
dunn_res <- dunn.test::dunn.test(
  df$risk_score,
  g      = df$PAM50_subtype,
  method = "bonferroni",
  altp   = TRUE
)

# Sample size labels per subtype
n_labels <- df %>%
  group_by(PAM50_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", n),
         y_pos = min(df$risk_score) - 0.12)

p_box <- ggplot(df, aes(x = PAM50_subtype, y = risk_score,
                        fill = PAM50_subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55,
               color = "grey30", linewidth = 0.4) +
  ggrastr::geom_point_rast(
    aes(color = PAM50_subtype),
    position = position_jitter(width = 0.18, seed = 42),
    size = 1.0, alpha = 0.45
  ) +
  geom_text(data = n_labels,
            aes(x = PAM50_subtype, y = y_pos, label = label),
            size = 3.2, color = "grey40", inherit.aes = FALSE) +
  scale_fill_manual(values  = pam50_colors) +
  scale_color_manual(values = pam50_colors) +
  labs(
    title    = "eRNA Risk Score Distribution across PAM50 Subtypes",
    subtitle = kw_label,
    x        = "PAM50 Subtype",
    y        = "eRNA Risk Score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle      = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x        = element_text(angle = 20, hjust = 1, size = 11),
    axis.title         = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

out_box <- file.path(out_dir, "Fig_RiskScore_by_PAM50.svg")
ggsave(out_box, p_box, width = 10, height = 7, device = "svg")
cat(sprintf("  Saved: %s\n", out_box))

# Figure 3: Multivariate Cox regression forest plot
cat("  [6-3] Multivariate Cox regression\n")

# Extract and clean staging variables from clinical data
clinical_ref <- data_clinical %>%
  dplyr::select(
    patient_id = submitter_id,
    age        = `age_at_earliest_diagnosis_in_years.diagnoses.xena_derived`,
    T_stage    = `ajcc_pathologic_t.diagnoses`,
    N_stage    = `ajcc_pathologic_n.diagnoses`,
    M_stage    = `ajcc_pathologic_m.diagnoses`
  ) %>%
  mutate(patient_id = substr(patient_id, 1, 12)) %>%
  distinct(patient_id, .keep_all = TRUE)

# Clean stage labels and encode as numeric
df_cox <- df %>%
  left_join(clinical_ref, by = "patient_id") %>%
  mutate(
    T_stage = gsub("[abc]", "", T_stage),
    N_stage = gsub("[abc]|\\(mol\\+\\)|mi|\\(i[-+]\\)", "", N_stage),
    M_stage = gsub("[abc]", "", M_stage),
    age     = as.numeric(age)
  ) %>%
  filter(
    T_stage %in% c("T1", "T2", "T3", "T4"),
    N_stage %in% c("N0", "N1", "N2", "N3"),
    M_stage %in% c("M0", "M1"),
    !is.na(age)
  ) %>%
  mutate(
    T_stage_num = as.numeric(as.factor(T_stage)),
    N_stage_num = as.numeric(as.factor(N_stage)),
    M_stage_num = as.numeric(as.factor(M_stage))
  )

# Set reference levels: Low risk, Luminal A subtype
df_cox$risk_group    <- factor(df_cox$risk_group,    levels = c("Low", "High"))
df_cox$PAM50_subtype <- relevel(df_cox$PAM50_subtype, ref = "Luminal A")

cat(sprintf("  Cox model samples: %d\n", nrow(df_cox)))

# Fit full model, fall back to simplified model on error
cox_fit <- tryCatch(
  coxph(Surv(OS.time_month, OS) ~
          risk_group + PAM50_subtype +
          age + T_stage_num + N_stage_num + M_stage_num,
        data = df_cox),
  error = function(e) {
    cat(sprintf("  Full model failed: %s\n  Trying simplified model...\n", e$message))
    coxph(Surv(OS.time_month, OS) ~ risk_group + PAM50_subtype + age,
          data = df_cox)
  }
)

cat("  Cox model summary:\n")
print(summary(cox_fit))

# Extract HR and confidence intervals
cox_sum  <- summary(cox_fit)
cox_coef <- as.data.frame(cox_sum$conf.int)
cox_coef$pval <- cox_sum$coefficients[, "Pr(>|z|)"]
cox_coef$term <- rownames(cox_coef)
colnames(cox_coef)[1:4] <- c("HR", "HR_inv", "CI_low", "CI_high")

# Map variable names to readable labels
label_map <- c(
  "risk_groupHigh"             = "Risk Group: High vs Low",
  "PAM50_subtypeLuminal B"     = "PAM50: Luminal B vs A",
  "PAM50_subtypeHER2-enriched" = "PAM50: HER2-enriched vs A",
  "PAM50_subtypeBasal-like"    = "PAM50: Basal-like vs A",
  "PAM50_subtypeNormal-like"   = "PAM50: Normal-like vs A",
  "age"                        = "Age",
  "T_stage_num"                = "T Stage",
  "N_stage_num"                = "N Stage",
  "M_stage_num"                = "M Stage"
)
cox_coef$label <- dplyr::recode(cox_coef$term, !!!label_map)
cox_coef$label <- factor(cox_coef$label, levels = rev(cox_coef$label))

# Color significant terms in red
cox_coef$dot_color <- ifelse(cox_coef$pval < 0.05, col_high, "grey50")
cox_coef$sig <- dplyr::case_when(
  cox_coef$pval < 0.001 ~ "***",
  cox_coef$pval < 0.01  ~ "**",
  cox_coef$pval < 0.05  ~ "*",
  TRUE                  ~ "ns"
)

# Forest plot on log-scale x-axis
p_forest <- ggplot(cox_coef, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.25, color = "black", linewidth = 0.5) +
  geom_point(aes(color = dot_color), size = 4, shape = 16) +
  geom_text(aes(label = sprintf("HR=%.3f [%.3f-%.3f]\np=%s",
                                HR, CI_low, CI_high,
                                ifelse(pval < 0.001, "< 0.001",
                                       sprintf("%.3f", pval)))),
            hjust = -0.08, size = 2.8, color = "grey30",
            lineheight = 0.9) +
  scale_color_identity() +
  scale_x_continuous(
    trans  = "log",
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25", "0.5", "1", "2", "4", "8"),
    expand = expansion(mult = c(0.05, 0.50))
  ) +
  labs(
    title    = "Multivariate Cox Regression",
    subtitle = "Reference: Low risk group; Luminal A subtype",
    x        = "Hazard Ratio (HR) with 95% CI",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.title.x       = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 10, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank()
  )

out_forest <- file.path(out_dir, "Fig_Multivariate_Cox_PAM50.svg")
ggsave(out_forest, p_forest, width = 10, height = 8, device = "svg")
cat(sprintf("  Saved: %s\n", out_forest))

# Figure 4: Stacked bar chart of PAM50 subtype distribution by risk group
cat("  [6-4] Subtype distribution by risk group\n")

dist_df <- df %>%
  group_by(risk_group, PAM50_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(risk_group) %>%
  mutate(pct = n / sum(n) * 100)

p_dist <- ggplot(dist_df,
                 aes(x = risk_group, y = pct, fill = PAM50_subtype)) +
  geom_bar(stat = "identity", width = 0.6,
           color = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(pct >= 5,
                               sprintf("%s\n%.1f%%", PAM50_subtype, pct), "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = pam50_colors, name = "PAM50 Subtype") +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "PAM50 Subtype Distribution by Risk Group",
    x     = "Risk Group",
    y     = "Proportion (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.position    = "right",
    axis.text.x        = element_text(size = 12, face = "bold"),
    axis.title         = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

out_dist <- file.path(out_dir, "Fig_HighRisk_Subtype_Distribution.svg")
ggsave(out_dist, p_dist, width = 7, height = 6, device = "svg")
cat(sprintf("  Saved: %s\n", out_dist))

# Export Cox results to CSV
cat("  [6-5] Exporting Cox results\n")

cox_table_out <- cox_coef %>%
  dplyr::select(label, HR, CI_low, CI_high, pval, sig) %>%
  dplyr::rename(
    Variable     = label,
    CI_lower_95  = CI_low,
    CI_upper_95  = CI_high,
    P_value      = pval,
    Significance = sig
  ) %>%
  dplyr::arrange(P_value)

out_csv <- file.path(out_dir, "Table_Multivariate_Cox_PAM50.csv")
write.csv(cox_table_out, out_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", out_csv))

# Summary
cat("\nAll analyses complete.\n")
cat(sprintf("  Final sample size: %d\n", nrow(df)))
output_files <- c(
  "Fig_KM_by_PAM50_subtype.svg",
  "Fig_RiskScore_by_PAM50.svg",
  "Fig_Multivariate_Cox_PAM50.svg",
  "Fig_HighRisk_Subtype_Distribution.svg",
  "Table_Multivariate_Cox_PAM50.csv"
)
for (f in output_files) {
  fp <- file.path(out_dir, f)
  cat(sprintf("  [%s] %s\n", ifelse(file.exists(fp), "OK", "MISSING"), f))
}
