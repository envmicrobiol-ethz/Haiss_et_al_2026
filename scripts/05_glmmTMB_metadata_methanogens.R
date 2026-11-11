# ==============================================================================
# FINAL MULTIVARIABLE GLMM ANALYSIS (RAW TR COUNTS + glmmTMB)
# ==============================================================================

# Required Libraries
library(dplyr)
library(stringr)
library(glmmTMB)     # Primary engine for Negative Binomial GLMMs
library(car)         # For Type II Wald Chisq Tests
library(broom.mixed) # For tidy model summaries and random effects
library(DHARMa)      # For simulated residual diagnostics
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# 0. SETUP & DATA CLEANING
# ------------------------------------------------------------------------------

# Read in data
df <- read.csv2("2026-07-21_RefSeq_methanogens_metadata.csv", header = TRUE, sep = ",")

# Ensure numeric conversions
df$TR_counts <- as.numeric(df$TR_counts)
df$genome_size <- as.numeric(df$genome_size)
df$ph_optimum_predicted <- as.numeric(df$ph_optimum_predicted)
df$salinity_optimum <- as.numeric(df$salinity_optimum)
df$temperature_optimum <- as.numeric(df$temperature_optimum)

# Clean dataset globally
df_clean <- df %>%
  filter(!is.na(TR_counts), !is.na(genome_size), genome_size > 0) %>%
  mutate(
    ecosystem_group = ifelse(is.na(ecosystem_group) | ecosystem_group == "", "Unknown", ecosystem_group),
    geographic_location = ifelse(is.na(geographic_location) | geographic_location == "", "Unknown", geographic_location),
    ecosystem_general = ifelse(is.na(ecosystem_general) | ecosystem_general == "", "Unknown", ecosystem_general),    
    Ecosystem = ifelse(is.na(Ecosystem) | Ecosystem == "", "Unknown", Ecosystem),
    oxygen = ifelse(is.na(oxygen) | oxygen == "", "Unknown", oxygen)
  )

# ------------------------------------------------------------------------------
# STEP 1: DEFINE CANDIDATE METADATA COLUMNS
# ------------------------------------------------------------------------------

candidate_cols <- c(
  "ph_optimum_predicted",
  "salinity_optimum",
  "temperature_optimum",
  "Methylotrophic.methanogenesis..MeOH.",
  "Methylotrophic.methanogenesis..TMA.",
  "Acetoclastic.methanogenesis..2.3.",
  "Hydrogenotrophic.methanogenesis"
)

# ------------------------------------------------------------------------------
# STEP 2: PREPROCESSING, TAXONOMY PARSING & DENSITY CALCULATION
# ------------------------------------------------------------------------------

analysis_data <- df_clean %>%
  mutate(
    # Extract Genus from Organism_name
    genus = Organism_name %>%
      str_remove("^uncultured\\s+") %>%
      str_remove("^Candidatus\\s+") %>%
      str_extract("^\\w+"),
    
    # Pre-calculate genome size metrics & explicit log term
    genome_size_mbp = genome_size / 1e6,
    log_genome_size_mbp = log(genome_size_mbp),
    genome_size_100kb = genome_size / 1e5,
    tr_density = TR_counts / genome_size_100kb
  ) %>%
  # Reclassify rare ecosystems (N < 6) to "other"
  group_by(ecosystem_group) %>%
  mutate(ecosystem_group = ifelse(n() < 6, "other", ecosystem_group)) %>%
  ungroup() %>%
  # Reclassify rare ecosystems (N < 3) to "other"
  group_by(ecosystem_general) %>%
  mutate(ecosystem_general = ifelse(n() < 6, "other", ecosystem_general)) %>%
  ungroup() %>%
  # Convert logical/boolean columns to factors
  mutate(across(where(is.logical), as.factor))

# ------------------------------------------------------------------------------
# STEP 3: SCREEN AND FILTER OUT CONSTANT COLUMNS
# ------------------------------------------------------------------------------

valid_covariates <- c()

for (col in candidate_cols) {
  if (col %in% colnames(analysis_data)) {
    non_na_vals <- na.omit(analysis_data[[col]])
    if (length(unique(non_na_vals)) > 1) {
      valid_covariates <- c(valid_covariates, col)
    } else {
      cat("Dropping constant/single-value column:", col, "\n")
    }
  }
}

# Retain complete cases across valid covariates
complete_model_data <- analysis_data %>%
  filter(if_all(all_of(c("TR_counts", "log_genome_size_mbp", valid_covariates)), ~ !is.na(.)))

cat("\nGenomes in final model:", nrow(complete_model_data), "\n")
cat("Retained non-constant covariates:", paste(valid_covariates, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# STEP 4: FIT GLMM (glmmTMB with Negative Binomial family)
# ------------------------------------------------------------------------------

fixed_terms <- c("log_genome_size_mbp", "as.factor(ecosystem_general)") # ecosystem_general / ecosystem_group
if (length(valid_covariates) > 0) {
  formatted_covs <- paste0("`", valid_covariates, "`")
  fixed_terms <- c(fixed_terms, formatted_covs)
}

glmm_formula <- as.formula(
  paste0("TR_counts ~ ", paste(fixed_terms, collapse = " + "), " + (1 | genus)")
)

cat("\nFitting Negative Binomial GLMM using glmmTMB...\n")
final_model <- glmmTMB(
  formula = glmm_formula, 
  data = complete_model_data, 
  family = nbinom2
)

# ------------------------------------------------------------------------------
# STEP 5: OUTPUT MODEL RESULTS & SCALING HYPOTHESIS TEST
# ------------------------------------------------------------------------------

# TYPE II ANCOVA TABLE (Wald Chisq Test)
cat("\n--- TYPE II ANCOVA TABLE (Wald Chisq) ---\n")
print(car::Anova(final_model, type = "II"))

# FIXED EFFECT COEFFICIENTS (Log-scale)
cat("\n--- FIXED EFFECT COEFFICIENTS (Log-scale Estimates) ---\n")
fixed_summary <- tidy(final_model, effects = "fixed", exponentiate = FALSE, conf.int = TRUE) %>%
  mutate(p.value_formatted = ifelse(p.value < 0.001, "< 0.001", as.character(round(p.value, 4)))) %>%
  arrange(p.value)
print(fixed_summary, n = 50)

# FIXED EFFECT COEFFICIENTS (Incidence Rate Ratios - IRR)
cat("\n--- FIXED EFFECT COEFFICIENTS (Incidence Rate Ratios - IRR) ---\n")
tidy(final_model, effects = "fixed", exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value_formatted = ifelse(p.value < 0.001, "< 0.001", as.character(round(p.value, 4)))) %>%
  arrange(p.value) %>%
  print(n = 50)

# DIRECT HYPOTHESIS TEST FOR TR DENSITY HYPER-SCALING (H0: beta = 1)
cat("\n--- HYPOTHESIS TEST: IS TR DENSITY CONSTANT ACROSS GENOME SIZES? (H0: beta = 1) ---\n")
beta_row <- fixed_summary %>% filter(term == "log_genome_size_mbp")
beta_hat <- beta_row$estimate
se_hat   <- beta_row$std.error

# Calculate Wald z-score and 2-tailed p-value against H0: beta = 1
z_stat <- (beta_hat - 1) / se_hat
p_val_density <- 2 * (1 - pnorm(abs(z_stat)))

cat(sprintf("Beta estimate: %.4f (SE: %.4f)\n", beta_hat, se_hat))
cat(sprintf("Wald z-statistic (against H0: beta = 1): %.4f\n", z_stat))
cat(sprintf("p-value: %s\n", ifelse(p_val_density < 0.001, "< 0.001", sprintf("%.4f", p_val_density))))

# TAXONOMIC LINEAGE VARIANCE
cat("\n--- TAXONOMIC LINEAGE VARIANCE ---\n")
ran_pars <- tidy(final_model, effects = "ran_pars")
var_genus <- (ran_pars %>% filter(group == "genus", term == "sd__(Intercept)"))$estimate^2
cat(sprintf("Genus Random Effect Variance (log scale): %.4f\n", var_genus))

# ==============================================================================
# STEP 6: TEST TAXONOMY & LINEAGE ENRICHMENT
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 EXTRACT GENUS RANDOM EFFECTS (DEVIATIONS FROM SCALING LAW)
# ------------------------------------------------------------------------------

genus_effects <- tidy(final_model, effects = "ran_vals") %>%
  filter(group == "genus") %>%
  rename(genus = level) %>%
  mutate(
    rate_factor = exp(estimate),
    ci_lower = estimate - (1.96 * std.error),
    ci_upper = estimate + (1.96 * std.error),
    significance = case_when(
      ci_lower > 0 ~ "Significantly Higher",
      ci_upper < 0 ~ "Significantly Lower",
      TRUE         ~ "Not Significantly Different"
    ),
    significance = factor(significance, levels = c("Significantly Higher", "Significantly Lower", "Not Significantly Different"))
  ) %>%
  arrange(desc(estimate))

cat("\n--- TOP 10 HIGHEST TR RATE GENERA (Relative to Genome Size Expectation) ---\n")
print(head(genus_effects %>% select(genus, estimate, rate_factor, std.error), 10))

cat("\n--- TOP 10 LOWEST TR RATE GENERA (Relative to Genome Size Expectation) ---\n")
print(tail(genus_effects %>% select(genus, estimate, rate_factor, std.error), 10))

cat("\n--- GENERA WITH SIGNIFICANTLY HIGHER TR RATES ---\n")
print(genus_effects %>% filter(significance == "Significantly Higher") %>% select(genus, estimate, ci_lower, ci_upper))

cat("\n--- GENERA WITH SIGNIFICANTLY LOWER TR RATES ---\n")
print(genus_effects %>% filter(significance == "Significantly Lower") %>% select(genus, estimate, ci_lower, ci_upper))

# ------------------------------------------------------------------------------
# 6.2 LIKELIHOOD RATIO TEST FOR GENUS RANDOM EFFECT
# ------------------------------------------------------------------------------

cat("\n--- LIKELIHOOD RATIO TEST FOR GENUS RANDOM EFFECT ---\n")

# Extract fixed formula without random effect
fixed_formula <- update(glmm_formula, . ~ . - (1 | genus))

# Fit fixed-effects Negative Binomial model via glmmTMB
null_model <- glmmTMB(fixed_formula, data = complete_model_data, family = nbinom2)

# LRT comparison
lrt_genus <- anova(final_model, null_model)
print(lrt_genus)

# ------------------------------------------------------------------------------
# 6.3 WILCOXON RANK-SUM TEST: Methanosarcina vs. Other Genera
# ------------------------------------------------------------------------------

complete_model_data <- complete_model_data %>%
  mutate(is_methanosarcina = ifelse(genus == "Methanosarcina", "Methanosarcina", "Other Genera"))

summary_stats <- complete_model_data %>%
  group_by(is_methanosarcina) %>%
  summarise(
    N = n(),
    Median_TR_Density = median(tr_density, na.rm = TRUE),
    IQR_TR_Density = IQR(tr_density, na.rm = TRUE)
  )

cat("\n--- SUMMARY STATISTICS (TR Density: TRs / 100 kb) ---\n")
print(summary_stats)

wilcox_res <- wilcox.test(tr_density ~ is_methanosarcina, data = complete_model_data)

cat("\n--- WILCOXON TEST RESULTS ---\n")
print(wilcox_res)

# ==============================================================================
# STEP 7: MODEL SELECTION & RESIDUAL DIAGNOSTICS
# ==============================================================================

# Fit minimal reduced model (Genome size + Genus random effect)
reduced_model <- glmmTMB(
  TR_counts ~ log_genome_size_mbp + (1 | genus), 
  data = complete_model_data, 
  family = nbinom2
)

cat("\n--- MODEL SELECTION (AIC COMPARISON) ---\n")
print(AIC(final_model, reduced_model))

# Simulated residual diagnostics using DHARMa
cat("\nRunning DHARMa residual diagnostics...\n")
simulationOutput <- simulateResiduals(fittedModel = final_model, plot = TRUE)










# ==============================================================================
# SUBSAMPLING METADATA PARAMETER SIGNIFICANCE ANALYSIS
# ==============================================================================

library(glmmTMB)
library(dplyr)
library(broom.mixed)
library(car)
library(ggplot2)

set.seed(42)

# Subsample sizes and iterations
target_sizes <- c(10, 20, 40, 80, 120, 140)
n_total <- nrow(complete_model_data)
n_reps <- 25  # Number of random draws per subsample size

all_terms_list <- list()

# ------------------------------------------------------------------------------
# HELPER FUNCTION: FIT GLMM, FIXED ANCOVA & RANDOM EFFECT LRT
# ------------------------------------------------------------------------------
fit_and_extract_anova <- function(sub_df, size, rep_id) {
  tryCatch({
    # 1. Fit Full Model (with random effect)
    fit_full <- glmmTMB(formula = glmm_formula, data = sub_df, family = nbinom2)
    
    # Extract Fixed Effects via Type II Wald Test
    anova_tab <- as.data.frame(car::Anova(fit_full, type = "II"))
    anova_tab$term <- rownames(anova_tab)
    p_col <- grep("Pr\\(>", colnames(anova_tab), value = TRUE)[1]
    
    fixed_res <- data.frame(
      sample_size = size,
      iteration = rep_id,
      term = anova_tab$term,
      p_value = anova_tab[[p_col]],
      is_significant = anova_tab[[p_col]] < 0.05,
      fit_success = TRUE,
      stringsAsFactors = FALSE
    )
    
    # 2. Evaluate Random Effect (1 | genus) via Likelihood Ratio Test (LRT)
    no_genus_formula <- as.formula(paste0("TR_counts ~ ", paste(fixed_terms, collapse = " + ")))
    fit_reduced <- glmmTMB(formula = no_genus_formula, data = sub_df, family = nbinom2)
    
    lrt <- anova(fit_reduced, fit_full)
    p_val_genus <- lrt$`Pr(>Chisq)`[2]
    
    random_res <- data.frame(
      sample_size = size,
      iteration = rep_id,
      term = "(1 | genus) [Random Effect]",
      p_value = p_val_genus,
      is_significant = !is.na(p_val_genus) & (p_val_genus < 0.05),
      fit_success = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Combine fixed terms + random term
    bind_rows(fixed_res, random_res)
    
  }, error = function(e) {
    # Fail-safe wrapper for tiny subsamples where GLMM fails to fit
    data.frame(
      sample_size = size,
      iteration = rep_id,
      term = "FAILED_FIT",
      p_value = NA,
      is_significant = FALSE,
      fit_success = FALSE,
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# LOOP OVER SUBSAMPLE SIZES
# ------------------------------------------------------------------------------
cat("Running parameter and random effect subsampling...\n")

for (s in target_sizes) {
  if (s >= n_total) next
  cat(sprintf("Sampling N = %d (%d iterations)...\n", s, n_reps))
  
  for (r in 1:n_reps) {
    sub_df <- complete_model_data %>% sample_n(size = s)
    res <- fit_and_extract_anova(sub_df, s, r)
    all_terms_list[[length(all_terms_list) + 1]] <- res
  }
}

# Add Full Dataset (N = total, 1 iteration)
cat(sprintf("Sampling Full Dataset (N = %d)...\n", n_total))
full_res <- fit_and_extract_anova(complete_model_data, n_total, 1)
all_terms_list[[length(all_terms_list) + 1]] <- full_res

# Combine all results
term_results <- bind_rows(all_terms_list) %>%
  filter(fit_success == TRUE)

# ------------------------------------------------------------------------------
# DISPLAY SUMMARY TABLE: % DETECTED AS SIGNIFICANT (p < 0.05)
# ------------------------------------------------------------------------------
cat("\n--- PARAMETER & RANDOM EFFECT DETECTION RATE ACROSS SAMPLE SIZES (% p < 0.05) ---\n")

param_summary <- term_results %>%
  group_by(sample_size, term) %>%
  summarise(
    Total_Runs = n(),
    Pct_Significant = round(mean(is_significant, na.rm = TRUE) * 100, 1),
    Median_P_Value = round(median(p_value, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(term, sample_size)

print(param_summary, n = 100)

# ------------------------------------------------------------------------------
# PLOT: PARAMETER & RANDOM EFFECT SIGNIFICANCE HEATMAP
# ------------------------------------------------------------------------------
param_summary_clean <- param_summary %>%
  mutate(
    term_clean = gsub("as.factor\\(|\\)|`", "", term),
    # Bold highlight for random effect in axis labels
    is_random = ifelse(grepl("Random Effect", term), "Random Effect", "Fixed Effect")
  )

plot_heatmap <- ggplot(param_summary_clean, aes(x = factor(sample_size), y = reorder(term_clean, Pct_Significant))) +
  geom_tile(aes(fill = Pct_Significant), color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(Pct_Significant, "%")), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#e74c3c", 
    mid = "#f1c40f", 
    high = "#2ecc71", 
    midpoint = 50,
    limits = c(0, 100),
    name = "% Significant\n(p < 0.05)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Covariate & Taxonomic Random Effect Stability Across Subsamples",
    subtitle = "Fixed effects tested via Type II Wald tests; Genus random effect tested via Likelihood Ratio Test (LRT)",
    x = "Subsample Size (N)",
    y = "Model Parameter"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(plot_heatmap)

# Save heatmap figure
ggsave("Figure_Subsampling_Parameter_and_RandomEffect_Heatmap.pdf", plot_heatmap, width = 9, height = 6.5)







# ==============================================================================
# STRATIFIED & BALANCED BOOTSTRAP ANALYSIS FOR GLMM
# ==============================================================================

library(glmmTMB)
library(dplyr)
library(car)
library(ggplot2)
library(patchwork)

set.seed(42)

# Number of bootstrap iterations
B <- 100  
min_group_size <- 8  # Size of smallest ecosystem group

# Store bootstrap results
boot_results_list <- list()

# ------------------------------------------------------------------------------
# HELPER FUNCTION: FIT GLMM, FIXED ANCOVA & RANDOM EFFECT LRT
# ------------------------------------------------------------------------------
fit_boot_model <- function(sub_df, boot_type, rep_id) {
  tryCatch({
    # 1. Fit full model
    fit_full <- glmmTMB(formula = glmm_formula, data = sub_df, family = nbinom2)
    
    # CHECK CONVERGENCE & HESSIAN MATRIX
    # If pdHess is FALSE or NA, the model did not properly converge
    is_pdHess <- !is.null(fit_full$sdr) && isTRUE(fit_full$sdr$pdHess)
    if (!is_pdHess) {
      stop("Non-positive-definite Hessian matrix encountered.")
    }
    
    # 2. Type II Wald test for fixed effects
    anova_tab <- as.data.frame(car::Anova(fit_full, type = "II"))
    anova_tab$term <- rownames(anova_tab)
    p_col <- grep("Pr\\(>", colnames(anova_tab), value = TRUE)[1]
    
    fixed_res <- data.frame(
      boot_type = boot_type,
      iteration = rep_id,
      term = anova_tab$term,
      p_value = anova_tab[[p_col]],
      is_significant = anova_tab[[p_col]] < 0.05,
      fit_success = TRUE,
      stringsAsFactors = FALSE
    )
    
    # 3. Likelihood Ratio Test for (1 | genus)
    no_genus_formula <- as.formula(paste0("TR_counts ~ ", paste(fixed_terms, collapse = " + ")))
    fit_reduced <- glmmTMB(formula = no_genus_formula, data = sub_df, family = nbinom2)
    
    lrt <- anova(fit_reduced, fit_full)
    p_val_genus <- lrt$`Pr(>Chisq)`[2]
    
    random_res <- data.frame(
      boot_type = boot_type,
      iteration = rep_id,
      term = "(1 | genus) [Random Effect]",
      p_value = p_val_genus,
      is_significant = !is.na(p_val_genus) & (p_val_genus < 0.05),
      fit_success = TRUE,
      stringsAsFactors = FALSE
    )
    
    bind_rows(fixed_res, random_res)
    
  }, error = function(e) {
    # Silently catch failed or non-converged fits
    data.frame(
      boot_type = boot_type,
      iteration = rep_id,
      term = "FAILED_FIT",
      p_value = NA,
      is_significant = FALSE,
      fit_success = FALSE,
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# RUN BOOTSTRAP LOOPS
# ------------------------------------------------------------------------------
cat(sprintf("Starting %d Bootstrap Iterations...\n", B))

for (b in 1:B) {
  if (b %% 20 == 0) cat(sprintf("  Completed %d/%d iterations...\n", b, B))
  
  # METHOD 1: Stratified Bootstrap (Preserving original group sizes)
  df_stratified <- complete_model_data %>%
    group_by(ecosystem_general) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    ungroup()
  
  res_strat <- fit_boot_model(df_stratified, "1. Stratified (Original Ratios)", b)
  
  # METHOD 2: Balanced Bootstrap (Equal N = 8 per ecosystem group)
  df_balanced <- complete_model_data %>%
    group_by(ecosystem_general) %>%
    slice_sample(n = min_group_size, replace = TRUE) %>%
    ungroup()
  
  res_bal <- fit_boot_model(df_balanced, "2. Balanced (Equal N = 8/Group)", b)
  
  boot_results_list[[length(boot_results_list) + 1]] <- res_strat
  boot_results_list[[length(boot_results_list) + 1]] <- res_bal
}

# Combine results
boot_results <- bind_rows(boot_results_list) %>%
  filter(fit_success == TRUE)

# ------------------------------------------------------------------------------
# SUMMARY TABLE: % DETECTION RATE (p < 0.05) ACROSS BOOTSTRAPS
# ------------------------------------------------------------------------------
cat("\n--- BOOTSTRAP PARAMETER SIGNIFICANCE SUMMARY (% p < 0.05) ---\n")

boot_summary <- boot_results %>%
  group_by(boot_type, term) %>%
  summarise(
    Successful_Fits = n(),
    Pct_Significant = round(mean(is_significant, na.rm = TRUE) * 100, 1),
    Median_P_Value = round(median(p_value, na.rm = TRUE), 4),
    IQR_P_Value = round(IQR(p_value, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(term, boot_type)

print(boot_summary, n = 100)

# ------------------------------------------------------------------------------
# COMPARISON PLOT: STRATIFIED VS. BALANCED BOOTSTRAP DETECTION RATES
# ------------------------------------------------------------------------------
boot_summary_clean <- boot_summary %>%
  mutate(term_clean = gsub("as.factor\\(|\\)|`", "", term))

plot_boot_comparison <- ggplot(
  boot_summary_clean, 
  aes(x = reorder(term_clean, Pct_Significant), y = Pct_Significant, fill = boot_type)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  geom_text(
    aes(label = paste0(Pct_Significant, "%")), 
    position = position_dodge(width = 0.8), 
    vjust = -0.4, size = 3.2, fontface = "bold"
  ) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("1. Stratified (Original Ratios)" = "#2b5c8f", "2. Balanced (Equal N = 8/Group)" = "#e67e22")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Parameter Significance under Stratified Bootstrapping",
    subtitle = "Comparing original ecosystem proportions vs. equalized group sizes (N = 8/group)",
    x = "Model Term",
    y = "% of Bootstrap Replicates with p < 0.05",
    fill = "Bootstrap Strategy"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(plot_boot_comparison)

# Save figure
ggsave("Figure_Bootstrap_Parameter_Stability.pdf", plot_boot_comparison, width = 10, height = 6)
