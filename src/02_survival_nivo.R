# =============================================================================
# 0. LOAD REQUIRED LIBRARIES
# =============================================================================
library(survival)    # For Kaplan–Meier estimation
library(flexsurv)    # For parametric/spline survival models
library(IPDfromKM)   # For reconstructing IPD from Kaplan–Meier curves
library(ggplot2)     # For visualizations
library(dplyr)       # For data manipulation
library(gridExtra)   # For arranging multiple plots
library(bbplot)      # For custom themes (e.g., theme_fivethirtyeight)
library(cowplot)     # For arranging multiple plots in a grid

# =============================================================================
# 1. READ AND PREPARE DIGITIZED SURVIVAL DATA
# =============================================================================

# Custom function to enforce non-increasing survival (monotonicity)
enforce_monotonicity <- function(df) {
  for (i in 2:nrow(df)) {
    min_previous <- min(df$percentage[1:(i - 1)])
    if (df$percentage[i] > min_previous) {
      df$percentage[i] <- min_previous
    }
  }
  return(df)
}

# Read survival event data
PFS <- read.csv("data/nivo 10yr_PFsurvival.csv", 
                header = FALSE, col.names = c("time", "percentage"))  %>% 
  arrange(time) %>%
  mutate(percentage = percentage/100)

OS  <- read.csv("data/nivo 10yr_survival.csv", 
                header = FALSE, col.names = c("time", "percentage"))  %>% 
  arrange(time) %>%
  mutate(percentage = percentage/100) 

# Enforce monotonic survival
OS  <- enforce_monotonicity(OS)
PFS <- enforce_monotonicity(PFS)

# Create data frames for IPD reconstruction
os_dat  <- OS  %>% 
  select(time, sur = percentage) %>% 
  arrange(time)

pfs_dat <- PFS %>% 
  select(time, sur = percentage) %>% 
  arrange(time)

# =============================================================================
# 2. DEFINE NUMBER-AT-RISK DATA
# =============================================================================

patients_at_risk <- list(
  os = list(
    timepoints = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 
                   60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120),
    nivo   = c(316, 265, 231, 201, 181, 171, 158, 145, 141, 137, 134, 130, 126, 
                   123, 118, 107, 102, 98, 96, 92, 77)
  ),
  pfs = list(
    timepoints = c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 
                   60, 66, 72, 78, 84, 90, 96, 102, 108, 114, 120),
    nivo   = c(316, 151, 120, 106, 97, 84, 78, 73, 69, 66, 62, 58,
                55, 52, 49, 45, 42, 40, 38, 24, 12)
  )
)




# =============================================================================
# 3. RECONSTRUCT INDIVIDUAL PATIENT DATA (IPD)
# =============================================================================

# Extract number-at-risk information
trisk_os  <- patients_at_risk$os$timepoints
nrisk_os  <- patients_at_risk$os$nivo
trisk_pfs <- patients_at_risk$pfs$timepoints
nrisk_pfs <- patients_at_risk$pfs$nivo

# Preprocess data (for IPDfromKM package)
pre_os  <- preprocess(dat = os_dat[, c("time", "sur")],
                      trisk = trisk_os,
                      nrisk = nrisk_os,
                      maxy  = 1)

pre_pfs <- preprocess(dat = pfs_dat[, c("time", "sur")],
                      trisk = trisk_pfs,
                      nrisk = nrisk_pfs,
                      maxy  = 1)

# Reconstruct IPD
est_os  <- getIPD(prep = pre_os)
est_pfs <- getIPD(prep = pre_pfs)

ipd_os  <- est_os$IPD  # Contains: time, status (1 = event, 0 = censored)
ipd_pfs <- est_pfs$IPD

# =============================================================================
# 4. PLOT RECONSTRUCTED KAPLAN–MEIER CURVES
# =============================================================================

# Helper function to prepare survival data for ggplot
prepare_km_data_with_ci <- function(km_fit, dataset_name) {
  data.frame(
    time    = km_fit$time,
    surv    = km_fit$surv,
    lower   = km_fit$lower,
    upper   = km_fit$upper,
    status  = ifelse(km_fit$n.censor > 0, 0, 1),  # 0 = Censored, 1 = Event
    dataset = dataset_name
  )
}

# Create Kaplan–Meier fits
km_os  <- survfit(Surv(time, status) ~ 1, data = ipd_os)
km_pfs <- survfit(Surv(time, status) ~ 1, data = ipd_pfs)

# Prepare data with confidence intervals
os_km_data  <- prepare_km_data_with_ci(km_os,  "Overall Survival (OS)")
pfs_km_data <- prepare_km_data_with_ci(km_pfs, "Progression-Free Survival (PFS)")
km_combined <- bind_rows(os_km_data, pfs_km_data)

# Plot the Kaplan–Meier curves with confidence intervals
ggplot(km_combined, aes(x = time, y = surv, color = dataset)) +
  geom_step(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dataset), alpha = 0.03) +
  scale_color_manual(values = c("blue", "orange")) +
  scale_fill_manual(values = c("blue", "orange")) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  geom_vline(xintercept = c(36, 60, 120), linetype = "dashed", color = "gray") +
  labs(title  = "Kaplan–Meier Survival Curves",
       x      = "Time (Months)",
       y      = "Survival Probability",
       color  = "Survival Type",
       fill   = "Survival Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# =============================================================================
# 5. FIT MODELS FOR OVERALL SURVIVAL (OS)
# =============================================================================

# 5A. Fit Parametric Models
param_models_os <- list(
  "Exponential"  = flexsurvreg(Surv(time, status) ~ 1, data = ipd_os, dist = "exp"),
  "Weibull"      = flexsurvreg(Surv(time, status) ~ 1, data = ipd_os, dist = "weibull"),
  "Log-Logistic" = flexsurvreg(Surv(time, status) ~ 1, data = ipd_os, dist = "llogis"),
  "Log-Normal"   = flexsurvreg(Surv(time, status) ~ 1, data = ipd_os, dist = "lnorm"),
  "Gompertz"     = flexsurvreg(Surv(time, status) ~ 1, data = ipd_os, dist = "gompertz")
)

# 5B. Fit Royston–Parmar Spline Models (k = 1 to 10)
spline_models_os <- list()
for (k_val in 1:6) {
  model_name <- paste0("Spline_k", k_val)
  spline_models_os[[model_name]] <- flexsurvspline(
    Surv(time, status) ~ 1,
    data  = ipd_os,
    k     = k_val,
    scale = "hazard"
  )
}

# 5C. Compare AIC for OS Models (Parametric + Spline)
all_os_models <- c(param_models_os, spline_models_os)
aic_os <- sapply(all_os_models, function(m) m$AIC)
aic_os <- sort(aic_os)
cat("AIC (OS) for all models (sorted):\n")
print(aic_os)

best_os_name  <- names(aic_os)[1]
best_os_model_n <- all_os_models[[best_os_name]]
cat("Best OS Model by AIC:", best_os_name, "\n")

# =============================================================================
# 6. PLOT KM CURVE & BEST OS MODEL
# =============================================================================

# Recreate the KM fit for OS (for plotting)
km_os <- survfit(Surv(time, status) ~ 1, data = ipd_os)
km_data_os <- data.frame(time = km_os$time, surv = km_os$surv)

# Define a time grid for prediction
max_time_plot <- max(ipd_os$time)  # or use max(c(ipd_os$time, ipd_pfs$time))
time_grid <- seq(0, max_time_plot, length.out = 200)

# Generate predicted survival from the best OS model
best_os_summary <- summary(best_os_model_n,
                           newdata = data.frame(time = time_grid),
                           type = "survival")

best_os_df <- as.data.frame(best_os_summary)

colnames(best_os_df) <- c("time", "surv", "lower", "upper")

ggplot() +
  geom_step(data = km_data_os, aes(x = time, y = surv),
            color = "black", size = 1) +
  geom_line(data = best_os_df, aes(x = time, y = surv),
            color = "red", linetype = "dashed", size = 1) +
  labs(title = paste("OS: KM vs. Best Model (", best_os_name, ")", sep = ""),
       x     = "Time (Months)",
       y     = "Survival Probability") +
  theme_minimal()

# =============================================================================
# 7. FIT MODELS FOR PROGRESSION-FREE SURVIVAL (PFS)
# =============================================================================

# 7A. Fit Parametric Models
param_models_pfs <- list(
  "Exponential"  = flexsurvreg(Surv(time, status) ~ 1, data = ipd_pfs, dist = "exp"),
  "Weibull"      = flexsurvreg(Surv(time, status) ~ 1, data = ipd_pfs, dist = "weibull"),
  "Log-Logistic" = flexsurvreg(Surv(time, status) ~ 1, data = ipd_pfs, dist = "llogis"),
  "Log-Normal"   = flexsurvreg(Surv(time, status) ~ 1, data = ipd_pfs, dist = "lnorm"),
  "Gompertz"     = flexsurvreg(Surv(time, status) ~ 1, data = ipd_pfs, dist = "gompertz")
)

# 7B. Fit a Spline Model for PFS
spline_models_pfs <- list()
for (k_val in 1:2) {
  model_name <- paste0("Spline_k", k_val)
  spline_models_pfs[[model_name]] <- flexsurvspline(
    Surv(time, status) ~ 1,
    data  = ipd_pfs,
    k     = k_val,
    scale = "hazard"
  )
}

# 7C. Compare AIC for PFS Models
all_pfs_models <- c(param_models_pfs, spline_models_pfs)
aic_pfs <- sapply(all_pfs_models, function(m) m$AIC)
aic_pfs <- sort(aic_pfs)
cat("\nAIC (PFS) for all models (sorted):\n")
print(aic_pfs)

best_pfs_name  <- names(aic_pfs)[1]
best_pfs_model_n <- all_pfs_models[[best_pfs_name]]
cat("Best PFS Model by AIC:", best_pfs_name, "\n")

# =============================================================================
# 8. PLOT KM CURVE & BEST PFS MODEL
# =============================================================================

km_pfs <- survfit(Surv(time, status) ~ 1, data = ipd_pfs)
km_data_pfs <- data.frame(time = km_pfs$time, surv = km_pfs$surv)

best_pfs_summary <- summary(best_pfs_model_n,
                            newdata = data.frame(time = time_grid),
                            type = "survival")

best_pfs_df <- as.data.frame(best_pfs_summary)

colnames(best_pfs_df) <- c("time", "surv", "lower", "upper")

ggplot() +
  geom_step(data = km_data_pfs, aes(x = time, y = surv),
            color = "black", size = 1) +
  geom_line(data = best_pfs_df, aes(x = time, y = surv),
            color = "red", linetype = "dashed", size = 1) +
  labs(title = paste("PFS: KM vs. Best Model (", best_pfs_name, ")", sep = ""),
       x     = "Time (Months)",
       y     = "Survival Probability") +
  theme_minimal()



# =============================================================================
# 9. DEFINE CUSTOM THEME FOR PUBLICATION & COMBINE PLOTS
# =============================================================================
library(ggthemes)
custom_theme <- theme_fivethirtyeight() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 12),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.key       = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")
  )

# Combined Plot A: OS and PFS models overlay
C <- ggplot() +
  geom_step(data = km_data_os, aes(x = time, y = surv),
            color = "grey40", size = 1) +
  geom_line(data = best_os_df, aes(x = time, y = surv),
            color = "#008Fd5", linetype = "dashed", size = 1) +
  geom_step(data = km_data_pfs, aes(x = time, y = surv),
            color = "grey40", size = 1) +
  geom_line(data = best_pfs_df, aes(x = time, y = surv),
            color = "#008Fd5", linetype = "dashed", size = 1) +
  scale_color_fivethirtyeight() +  # Use bbplot's color scale
  labs(x = "Time (months)",
       y = "Survival Probability",
       color = "Model") +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 1)) +
  custom_theme

# =============================================================================
# 10. VALIDATION AGAINST PUBLISHED DATA
#     (Timepoints: 36, 60, 120; plus median)
# =============================================================================

timepoints <- c(36, 60, 120)

# Generate model predictions at specified timepoints
os_preds_at  <- data.frame(summary(best_os_model_n, t = timepoints, type = "survival"))
pfs_preds_at <- data.frame(summary(best_pfs_model_n, t = timepoints, type = "survival"))

# Published data for comparison
timepoint_estimates_df <- data.frame(
  time     = rep(c(36, 60, 120), 2),
  endpoint = rep(c("OS Model Estimate", "PFS Model Estimate"), each = 3),
  estimate = c(0.51, 0.44, 0.37, 0.32, 0.29, 0.23),
  ci_lower = c(0.45, 0.39, 0.32, 0.26, 0.24, 0.18),
  ci_upper = c(0.56, 0.50, 0.43, 0.37, 0.35, 0.29)
)

# Convert time to factor for plotting purposes
os_preds_at$time <- as.factor(os_preds_at$time)
pfs_preds_at$time <- as.factor(pfs_preds_at$time)
timepoint_estimates_df$time <- as.factor(timepoint_estimates_df$time)

# Add endpoint column to model estimates
os_preds_at$endpoint  <- "OS Model Estimate"
pfs_preds_at$endpoint <- "PFS Model Estimate"

# Rename columns to match
colnames(os_preds_at)  <- c("time", "estimate", "ci_lower", "ci_upper", "endpoint")
colnames(pfs_preds_at) <- c("time", "estimate", "ci_lower", "ci_upper", "endpoint")

# Combine OS and PFS model estimates
model_preds <- bind_rows(os_preds_at, pfs_preds_at)

D <- ggplot() +
  # Published estimates (shifted left)
  geom_point(data = timepoint_estimates_df, 
             aes(x = time, y = estimate, color = "Published Data"), 
             size = 3, shape = 17, position = position_dodge(width = 1)) +
  geom_errorbar(data = timepoint_estimates_df, 
                aes(x = time, ymin = ci_lower, ymax = ci_upper, color = "Published Data"), 
                width = 0.2, position = position_dodge(width = 1)) +
  geom_point(data = model_preds, 
             aes(x = time, y = estimate, color = endpoint), 
             size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(data = model_preds, 
                aes(x = time, ymin = ci_lower, ymax = ci_upper, color = endpoint), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  facet_wrap(~endpoint, ncol = 1) +
  scale_color_manual(values = c("OS Model Estimate" = "#008Fd5", 
                                "PFS Model Estimate" = "#008Fd5", 
                                "Published Data"     = "black")) +
  labs(x = "Time (months)",
       y = "",
       color = "") +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 1)) +
  custom_theme +
  theme(legend.position = "none")

library(cowplot)

# Create individual titles
title_top <- ggdraw() + draw_label("Nivolumab + Ipilimumab", fontface = 'bold', hjust = 0.5)
title_bottom <- ggdraw() + draw_label("Nivolumab", fontface = 'bold', hjust = 0.5)

# Arrange the plots with titles
combined_top <- plot_grid(title_top, plot_grid(A, B, labels = c("A", "B"), ncol = 2, rel_widths = c(1.25, 1)), ncol = 1, rel_heights = c(0.1, 1))
combined_bottom <- plot_grid(title_bottom, plot_grid(C, D, labels = c("C", "D"), ncol = 2, rel_widths = c(1.25, 1)), ncol = 1, rel_heights = c(0.1, 1))

# Combine everything
combined_plot <- plot_grid(combined_top, combined_bottom, ncol = 1)

# Display the final plot
print(combined_plot)

ggsave("survival_curves.png", combined_plot, width = 12, height = 10, dpi = 300, bg = "white")
