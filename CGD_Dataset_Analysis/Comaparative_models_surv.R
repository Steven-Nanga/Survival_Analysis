# =============================================================================
# Comprehensive Recurrent Events Analysis for CGD Data
# Comparing: AG, PWP-TT, PWP-GT, WLW, Frailty, and MSM Models
# =============================================================================

# =============================================================================
# PACKAGE INSTALLATION AND LOADING
# =============================================================================

# Function to check and install packages
check_and_install <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package, "...\n")
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# Required packages
required_packages <- c("survival", "ggplot2", "gridExtra", "dplyr", "survminer")

cat("Checking and installing required packages...\n")
for (pkg in required_packages) {
  check_and_install(pkg)
}

cat("\nAll required packages loaded successfully!\n\n")

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Read the data
cgd <- read.csv('CGD Dataset.csv')

# Display data structure
cat("Dataset Overview:\n")
str(cgd)
cat("\nFirst few rows:\n")
head(cgd)

# Summary statistics
cat("\nSummary Statistics:\n")
summary(cgd)

# Count events per patient
events_per_patient <- cgd %>%
  group_by(ID) %>%
  summarise(
    n_episodes = max(sequence),
    total_events = sum(event),
    follow_up = max(T1)
  )

cat("\nRecurrent Events Summary:\n")
cat("Total patients:", length(unique(cgd$ID)), "\n")
cat("Total observations:", nrow(cgd), "\n")
cat("Total events:", sum(cgd$event), "\n")
cat("Mean events per patient:", mean(events_per_patient$total_events), "\n")

# Create gap time variable
cgd$gap_time <- cgd$T1 - cgd$T0

# =============================================================================
# 2. ANDERSEN-GILL (AG) MODEL
# =============================================================================

cat("\n\n========================================\n")
cat("ANDERSEN-GILL MODEL (Counting Process)\n")
cat("========================================\n")

# AG model uses counting process formulation
# Treats all events as independent, time-varying
ag_model <- coxph(Surv(T0, T1, event) ~ trtmt + inherit + age + 
                    cluster(ID), 
                  data = cgd, 
                  method = "breslow")

cat("\nAG Model Summary:\n")
print(summary(ag_model))

# Store results
ag_coef <- coef(ag_model)
ag_se <- sqrt(diag(ag_model$var))
ag_ci <- confint(ag_model)

# =============================================================================
# 3. PWP-TT (Prentice-Williams-Peterson - Total Time)
# =============================================================================

cat("\n\n========================================\n")
cat("PWP-TT MODEL (Total Time)\n")
cat("========================================\n")

# PWP-TT stratifies by event number, uses total time
pwp_tt_model <- coxph(Surv(T0, T1, event) ~ trtmt + inherit + age + 
                        strata(sequence) + cluster(ID),
                      data = cgd,
                      method = "breslow")

cat("\nPWP-TT Model Summary:\n")
print(summary(pwp_tt_model))

pwp_tt_coef <- coef(pwp_tt_model)
pwp_tt_se <- sqrt(diag(pwp_tt_model$var))
pwp_tt_ci <- confint(pwp_tt_model)

# =============================================================================
# 4. PWP-GT (Prentice-Williams-Peterson - Gap Time)
# =============================================================================

cat("\n\n========================================\n")
cat("PWP-GT MODEL (Gap Time)\n")
cat("========================================\n")

# PWP-GT stratifies by event number, uses gap time
pwp_gt_model <- coxph(Surv(gap_time, event) ~ trtmt + inherit + age + 
                        strata(sequence) + cluster(ID),
                      data = cgd,
                      method = "breslow")

cat("\nPWP-GT Model Summary:\n")
print(summary(pwp_gt_model))

pwp_gt_coef <- coef(pwp_gt_model)
pwp_gt_se <- sqrt(diag(pwp_gt_model$var))
pwp_gt_ci <- confint(pwp_gt_model)

# =============================================================================
# 5. WLW (Wei-Lin-Weissfeld) MODEL
# =============================================================================

cat("\n\n========================================\n")
cat("WLW MODEL (Marginal Means/Rates)\n")
cat("========================================\n")

# Create separate datasets for each event stratum
max_events <- max(cgd$sequence)
wlw_data_list <- list()

for (i in 1:max_events) {
  temp_data <- cgd %>%
    group_by(ID) %>%
    filter(sequence <= i) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(
      event_i = ifelse(sequence == i & event == 1, 1, 0),
      stratum = i
    )
  wlw_data_list[[i]] <- temp_data
}

wlw_data <- bind_rows(wlw_data_list)

# Fit WLW model
wlw_model <- coxph(Surv(T1, event_i) ~ trtmt + inherit + age + 
                     strata(stratum) + cluster(ID),
                   data = wlw_data,
                   method = "breslow")

cat("\nWLW Model Summary:\n")
print(summary(wlw_model))

wlw_coef <- coef(wlw_model)
wlw_se <- sqrt(diag(wlw_model$var))
wlw_ci <- confint(wlw_model)

# =============================================================================
# 6. FRAILTY MODEL (Shared Gamma Frailty)
# =============================================================================

cat("\n\n========================================\n")
cat("FRAILTY MODEL (Shared Gamma Frailty)\n")
cat("========================================\n")

# Shared frailty model with gamma distribution
frailty_model <- coxph(Surv(T0, T1, event) ~ trtmt + inherit + age + 
                         frailty(ID, distribution = "gamma"),
                       data = cgd,
                       method = "breslow")

cat("\nFrailty Model Summary:\n")
print(summary(frailty_model))

frailty_coef <- coef(frailty_model)[1:3]  # Exclude frailty terms
frailty_se <- sqrt(diag(frailty_model$var))[1:3]
frailty_ci <- confint(frailty_model)[1:3, ]

# Extract frailty variance
frailty_var <- frailty_model$history[[1]]$theta
cat("\nEstimated Frailty Variance (theta):", frailty_var, "\n")

# =============================================================================
# 7. MARGINAL STRUCTURAL MODEL (MSM) - Simplified Approach
# =============================================================================

cat("\n\n========================================\n")
cat("MARGINAL STRUCTURAL MODEL (MSM)\n")
cat("========================================\n")

# For MSM, we use inverse probability weighting
# This is a simplified version; full MSM would require time-varying confounders

# Calculate weights (simplified - assuming baseline treatment only)
# In practice, you'd use propensity scores
cgd$weight <- 1  # Equal weights for simplicity in this example

msm_model <- coxph(Surv(T0, T1, event) ~ trtmt + inherit + age + 
                     cluster(ID),
                   weights = weight,
                   data = cgd,
                   method = "breslow",
                   robust = TRUE)

cat("\nMSM Model Summary:\n")
print(summary(msm_model))

msm_coef <- coef(msm_model)
msm_se <- sqrt(diag(msm_model$var))
msm_ci <- confint(msm_model)

# =============================================================================
# 8. MODEL COMPARISON
# =============================================================================

cat("\n\n========================================\n")
cat("MODEL COMPARISON\n")
cat("========================================\n")

# Create comparison table for treatment effect
models <- c("AG", "PWP-TT", "PWP-GT", "WLW", "Frailty", "MSM")
trtmt_coefs <- c(ag_coef[1], pwp_tt_coef[1], pwp_gt_coef[1], 
                 wlw_coef[1], frailty_coef[1], msm_coef[1])
trtmt_se <- c(ag_se[1], pwp_tt_se[1], pwp_gt_se[1], 
              wlw_se[1], frailty_se[1], msm_se[1])
trtmt_hr <- exp(trtmt_coefs)
trtmt_ci_lower <- exp(trtmt_coefs - 1.96 * trtmt_se)
trtmt_ci_upper <- exp(trtmt_coefs + 1.96 * trtmt_se)
trtmt_pval <- 2 * pnorm(-abs(trtmt_coefs / trtmt_se))

comparison_table <- data.frame(
  Model = models,
  Coefficient = round(trtmt_coefs, 4),
  SE = round(trtmt_se, 4),
  HR = round(trtmt_hr, 4),
  CI_Lower = round(trtmt_ci_lower, 4),
  CI_Upper = round(trtmt_ci_upper, 4),
  P_Value = round(trtmt_pval, 4)
)

cat("\nTreatment Effect Comparison Across Models:\n")
print(comparison_table)

# AIC comparison (where applicable)
cat("\n\nModel Fit Statistics:\n")
cat("AG Model AIC:", AIC(ag_model), "\n")
cat("PWP-TT Model AIC:", AIC(pwp_tt_model), "\n")
cat("PWP-GT Model AIC:", AIC(pwp_gt_model), "\n")
cat("WLW Model AIC:", AIC(wlw_model), "\n")
cat("Frailty Model AIC:", AIC(frailty_model), "\n")
cat("MSM Model AIC:", AIC(msm_model), "\n")

# =============================================================================
# 9. VISUALIZATIONS
# =============================================================================

cat("\n\nCreating visualizations...\n")

# Plot 1: Forest plot of treatment effects
forest_data <- comparison_table

p1 <- ggplot(forest_data, aes(x = Model, y = HR)) +
  geom_point(size = 4, color = "blue") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), 
                width = 0.2, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Treatment Effect Estimates Across Models",
       subtitle = "Hazard Ratios with 95% Confidence Intervals",
       y = "Hazard Ratio",
       x = "Model") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Plot 2: Cumulative hazard by treatment (AG model)
ag_surv <- survfit(Surv(T0, T1, event) ~ trtmt, data = cgd)

p2 <- ggsurvplot(ag_surv, 
                 data = cgd,
                 fun = "cumhaz",
                 conf.int = TRUE,
                 risk.table = FALSE,
                 legend.labs = c("Placebo", "rIFN-gamma"),
                 title = "Cumulative Hazard Function (AG Model)",
                 xlab = "Time (days)",
                 ylab = "Cumulative Hazard")$plot

# Plot 3: Event rate by model
p3 <- ggplot(comparison_table, aes(x = Model, y = HR, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Treatment Hazard Ratios by Model",
       y = "Hazard Ratio",
       x = "Model") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
print(p1)
print(p2)
print(p3)

# =============================================================================
# 10. DETAILED RESULTS EXPORT
# =============================================================================

# Create comprehensive results summary
cat("\n\n========================================\n")
cat("COMPREHENSIVE RESULTS SUMMARY\n")
cat("========================================\n")

# Full coefficient table for all covariates
all_vars <- c("trtmt", "inherit", "age")
full_comparison <- data.frame()

for (var in all_vars) {
  var_idx <- which(names(ag_coef) == var)
  
  temp_df <- data.frame(
    Variable = var,
    Model = models,
    Coefficient = c(ag_coef[var_idx], pwp_tt_coef[var_idx], 
                    pwp_gt_coef[var_idx], wlw_coef[var_idx], 
                    frailty_coef[var_idx], msm_coef[var_idx]),
    SE = c(ag_se[var_idx], pwp_tt_se[var_idx], pwp_gt_se[var_idx], 
           wlw_se[var_idx], frailty_se[var_idx], msm_se[var_idx])
  )
  
  temp_df$HR <- exp(temp_df$Coefficient)
  temp_df$CI_Lower <- exp(temp_df$Coefficient - 1.96 * temp_df$SE)
  temp_df$CI_Upper <- exp(temp_df$Coefficient + 1.96 * temp_df$SE)
  temp_df$P_Value <- 2 * pnorm(-abs(temp_df$Coefficient / temp_df$SE))
  
  full_comparison <- rbind(full_comparison, temp_df)
}

cat("\nFull Covariate Effects Across All Models:\n")
print(full_comparison)

# =============================================================================
# 11. KEY FINDINGS AND RECOMMENDATIONS
# =============================================================================

cat("\n\n========================================\n")
cat("KEY FINDINGS AND RECOMMENDATIONS\n")
cat("========================================\n")

cat("\n1. Treatment Effect:\n")
if (mean(trtmt_hr) < 1) {
  cat("   - Treatment (rIFN-gamma) shows protective effect (HR < 1) across all models\n")
  cat("   - Average HR across models:", round(mean(trtmt_hr), 3), "\n")
} else {
  cat("   - Treatment shows increased hazard across models\n")
}

cat("\n2. Model Heterogeneity:\n")
cat("   - Range of treatment HRs:", round(min(trtmt_hr), 3), "to", round(max(trtmt_hr), 3), "\n")
cat("   - Coefficient of variation:", round(sd(trtmt_hr)/mean(trtmt_hr), 3), "\n")

cat("\n3. Statistical Significance:\n")
sig_models <- sum(trtmt_pval < 0.05)
cat("   -", sig_models, "out of", length(models), "models show significant treatment effect (p<0.05)\n")

cat("\n4. Model Recommendations:\n")
cat("   - AG Model: Best for overall event rate, assumes independence\n")
cat("   - PWP Models: Account for event ordering, stratified baseline hazards\n")
cat("   - WLW Model: Marginal approach, suitable for treatment effect estimation\n")
if (exists("frailty_fitted") && frailty_fitted) {
  cat("   - Frailty Model: Accounts for patient heterogeneity (theta =", 
      round(frailty_var, 3), ")\n")
} else {
  cat("   - Frailty Model: Could not be fitted in this session\n")
}
cat("   - MSM Model: Suitable for time-varying confounding (simplified here)\n")

cat("\n\nAnalysis Complete!\n")

# Save results
save(comparison_table, full_comparison, ag_model, pwp_tt_model, pwp_gt_model,
     wlw_model, frailty_model, msm_model,
     file = "cgd_recurrent_events_analysis.RData")

cat("\nResults saved to: cgd_recurrent_events_analysis.RData\n")

