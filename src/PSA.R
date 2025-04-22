########################################################################
# PARTITIONED SURVIVAL MODEL - PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
# Using existing survival curves
########################################################################

library(flexsurv)
library(ggplot2)

# 1) Time Horizon & Discounting
########################################################################
max_months <- 120    # 10 years
cycle_length <- 1
time_vec <- seq(0, max_months, by = cycle_length)
disc_rate <- 0.03
discount_factor <- function(m) {
  # discrete monthly discount factor for cycle m
  (1 + disc_rate)^(-(m / 12))
}

# 2) Use your existing survival curves
########################################################################
# The user already has:
#   best_os_model_ni, best_pfs_model_ni  => NIVO+IPI
#   best_os_model_n,  best_pfs_model_n   => NIVO alone

get_surv <- function(model, times) {
  # Summaries survival from a flexsurv model
  s_list <- summary(model, t = times, type = "survival")
  sapply(s_list, function(x) x$est)
}

# OS/PFS for combination - use the existing curves
os_ni  <- get_surv(best_os_model_ni, time_vec)
pfs_ni <- get_surv(best_pfs_model_ni, time_vec)

# OS/PFS for NIVO alone - use the existing curves
os_n  <- get_surv(best_os_model_n, time_vec)
pfs_n <- get_surv(best_pfs_model_n, time_vec)

# 3) PSA Setup
########################################################################
n_simulations <- 1000  # Number of PSA iterations

# 4) Define Parameter Distributions
########################################################################
# Drug costs (gamma distributions)
cost_nivolumab_mean <- 14975
cost_nivolumab_se <- (21417 - 9703) / (2 * 1.96)  # Estimated from 95% CI
cost_nivolumab_shape <- (cost_nivolumab_mean^2) / (cost_nivolumab_se^2)
cost_nivolumab_scale <- (cost_nivolumab_se^2) / cost_nivolumab_mean

cost_combo_mean <- 26425
cost_combo_se <- (37662 - 17089) / (2 * 1.96)  # Estimated from 95% CI
cost_combo_shape <- (cost_combo_mean^2) / (cost_combo_se^2)
cost_combo_scale <- (cost_combo_se^2) / cost_combo_mean

# Disease costs
cost_stable_mean <- 2166
cost_stable_se <- (3098 - 1397) / (2 * 1.96)  # Estimated from 95% CI
cost_stable_shape <- (cost_stable_mean^2) / (cost_stable_se^2)
cost_stable_scale <- (cost_stable_se^2) / cost_stable_mean

# Second-line treatment costs - assuming 20% variation
cost_prog_combo_mean <- 8908
cost_prog_combo_cv <- 0.2  # Coefficient of variation
cost_prog_combo_se <- cost_prog_combo_mean * cost_prog_combo_cv
cost_prog_combo_shape <- (cost_prog_combo_mean^2) / (cost_prog_combo_se^2)
cost_prog_combo_scale <- (cost_prog_combo_se^2) / cost_prog_combo_mean

cost_prog_nivo_mean <- 12093
cost_prog_nivo_cv <- 0.2
cost_prog_nivo_se <- cost_prog_nivo_mean * cost_prog_nivo_cv
cost_prog_nivo_shape <- (cost_prog_nivo_mean^2) / (cost_prog_nivo_se^2)
cost_prog_nivo_scale <- (cost_prog_nivo_se^2) / cost_prog_nivo_mean

# Palliative care cost - assuming 20% variation
cost_palliative_death_mean <- 15957
cost_palliative_death_cv <- 0.2
cost_palliative_death_se <- cost_palliative_death_mean * cost_palliative_death_cv
cost_palliative_death_shape <- (cost_palliative_death_mean^2) / (cost_palliative_death_se^2)
cost_palliative_death_scale <- (cost_palliative_death_se^2) / cost_palliative_death_mean

# Utilities (beta distributions)
util_stable_base_mean <- 0.754
util_stable_base_se <- 0.03  # Assuming SE of 0.03
util_stable_base_alpha <- util_stable_base_mean * ((util_stable_base_mean * (1 - util_stable_base_mean) / (util_stable_base_se^2)) - 1)
util_stable_base_beta <- (1 - util_stable_base_mean) * ((util_stable_base_mean * (1 - util_stable_base_mean) / (util_stable_base_se^2)) - 1)

util_prog_decrement_mean <- 0.180
util_prog_decrement_se <- 0.02  # Assuming SE of 0.02
util_prog_decrement_alpha <- util_prog_decrement_mean * ((util_prog_decrement_mean * (1 - util_prog_decrement_mean) / (util_prog_decrement_se^2)) - 1)
util_prog_decrement_beta <- (1 - util_prog_decrement_mean) * ((util_prog_decrement_mean * (1 - util_prog_decrement_mean) / (util_prog_decrement_se^2)) - 1)

# AE penalties (beta distributions)
ae_penalty_combo_mean <- 0.012
ae_penalty_combo_se <- 0.003  # Assuming SE of 0.003
ae_penalty_combo_alpha <- ae_penalty_combo_mean * ((ae_penalty_combo_mean * (1 - ae_penalty_combo_mean) / (ae_penalty_combo_se^2)) - 1)
ae_penalty_combo_beta <- (1 - ae_penalty_combo_mean) * ((ae_penalty_combo_mean * (1 - ae_penalty_combo_mean) / (ae_penalty_combo_se^2)) - 1)

ae_penalty_nivo_mean <- 0.0036
ae_penalty_nivo_se <- 0.0009  # Assuming SE of 0.0009
ae_penalty_nivo_alpha <- ae_penalty_nivo_mean * ((ae_penalty_nivo_mean * (1 - ae_penalty_nivo_mean) / (ae_penalty_nivo_se^2)) - 1)
ae_penalty_nivo_beta <- (1 - ae_penalty_nivo_mean) * ((ae_penalty_nivo_mean * (1 - ae_penalty_nivo_mean) / (ae_penalty_nivo_se^2)) - 1)

# 5) Pre-calculate health state partitioning
########################################################################
# Using the fixed survival curves
stable_ni <- pfs_ni
progress_ni <- pmax(os_ni - pfs_ni, 0)  # Avoid negative values
dead_ni <- 1 - os_ni

stable_n <- pfs_n
progress_n <- pmax(os_n - pfs_n, 0)  # Avoid negative values
dead_n <- 1 - os_n

# Fixed therapy parameters
max_therapy_months <- 24

# 6) Run PSA simulations
########################################################################
# Prepare to store results
psa_results <- data.frame(
  iter = 1:n_simulations,
  cost_combo = numeric(n_simulations),
  qaly_combo = numeric(n_simulations),
  cost_nivo = numeric(n_simulations),
  qaly_nivo = numeric(n_simulations),
  inc_cost = numeric(n_simulations),
  inc_qaly = numeric(n_simulations),
  icer = numeric(n_simulations)
)

# Run simulations
for (i in 1:n_simulations) {
  # 1) Sample parameter values for this iteration
  
  # Cost parameters (from gamma distributions)
  cost_nivolumab <- rgamma(1, shape = cost_nivolumab_shape, scale = cost_nivolumab_scale)
  cost_combo <- rgamma(1, shape = cost_combo_shape, scale = cost_combo_scale)
  cost_stable <- rgamma(1, shape = cost_stable_shape, scale = cost_stable_scale)
  cost_prog_combo <- rgamma(1, shape = cost_prog_combo_shape, scale = cost_prog_combo_scale)
  cost_prog_nivo <- rgamma(1, shape = cost_prog_nivo_shape, scale = cost_prog_nivo_scale)
  cost_palliative_death <- rgamma(1, shape = cost_palliative_death_shape, scale = cost_palliative_death_scale)
  
  # Utility parameters (from beta distributions)
  util_stable_base <- rbeta(1, shape1 = util_stable_base_alpha, shape2 = util_stable_base_beta)
  util_prog_decrement <- rbeta(1, shape1 = util_prog_decrement_alpha, shape2 = util_prog_decrement_beta)
  util_prog <- util_stable_base - util_prog_decrement
  if (util_prog < 0) util_prog <- 0  # Ensure utility isn't negative
  
  ae_penalty_combo <- rbeta(1, shape1 = ae_penalty_combo_alpha, shape2 = ae_penalty_combo_beta)
  ae_penalty_nivo <- rbeta(1, shape1 = ae_penalty_nivo_alpha, shape2 = ae_penalty_nivo_beta)
  
  # 2) Calculate costs and QALYs for both arms using sampled parameters
  
  # Combo arm calculation
  monthly_cost_combo <- numeric(length(time_vec))
  monthly_qaly_combo <- numeric(length(time_vec))
  
  for (m in seq_along(time_vec)) {
    # Distribution at month m
    fs <- stable_ni[m]
    fp <- progress_ni[m]
    
    # Drug cost
    if (m <= max_therapy_months) {
      drug_c <- fs * cost_combo
    } else {
      drug_c <- 0
    }
    
    # Progressed cost
    prog_c <- fp * cost_prog_combo
    
    # Stable disease cost
    stable_c <- fs * cost_stable
    
    cost_m <- drug_c + prog_c + stable_c
    
    # Utility
    if (m <= max_therapy_months) {
      u_stable <- util_stable_base - ae_penalty_combo
    } else {
      u_stable <- util_stable_base
    }
    
    u_prog <- util_prog
    util_m <- fs * u_stable + fp * u_prog
    
    # Discount
    dfact <- discount_factor(m-1)
    monthly_cost_combo[m] <- cost_m * dfact
    monthly_qaly_combo[m] <- util_m * (1/12) * dfact
  }
  
  # Sum cost & qaly
  total_cost_combo <- sum(monthly_cost_combo)
  total_qaly_combo <- sum(monthly_qaly_combo)
  
  # Palliative care costs
  nd_sum <- 0
  for(m in 2:length(time_vec)) {
    newly_dead <- (dead_ni[m] - dead_ni[m-1])
    if(newly_dead < 0) newly_dead <- 0
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead * cost_palliative_death * dfact
  }
  total_cost_combo <- total_cost_combo + nd_sum
  
  # Nivo alone arm calculation
  monthly_cost_nivo <- numeric(length(time_vec))
  monthly_qaly_nivo <- numeric(length(time_vec))
  
  for (m in seq_along(time_vec)) {
    # Distribution at month m
    fs <- stable_n[m]
    fp <- progress_n[m]
    
    # Drug cost
    if (m <= max_therapy_months) {
      drug_c <- fs * cost_nivolumab
    } else {
      drug_c <- 0
    }
    
    # Progressed cost
    prog_c <- fp * cost_prog_nivo
    
    # Stable disease cost
    stable_c <- fs * cost_stable
    
    cost_m <- drug_c + prog_c + stable_c
    
    # Utility
    if (m <= max_therapy_months) {
      u_stable <- util_stable_base - ae_penalty_nivo
    } else {
      u_stable <- util_stable_base
    }
    
    u_prog <- util_prog
    util_m <- fs * u_stable + fp * u_prog
    
    # Discount
    dfact <- discount_factor(m-1)
    monthly_cost_nivo[m] <- cost_m * dfact
    monthly_qaly_nivo[m] <- util_m * (1/12) * dfact
  }
  
  # Sum cost & qaly
  total_cost_nivo <- sum(monthly_cost_nivo)
  total_qaly_nivo <- sum(monthly_qaly_nivo)
  
  # Palliative care costs
  nd_sum <- 0
  for(m in 2:length(time_vec)) {
    newly_dead <- (dead_n[m] - dead_n[m-1])
    if(newly_dead < 0) newly_dead <- 0
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead * cost_palliative_death * dfact
  }
  total_cost_nivo <- total_cost_nivo + nd_sum
  
  # Calculate incremental values
  inc_cost <- total_cost_combo - total_cost_nivo
  inc_qaly <- total_qaly_combo - total_qaly_nivo
  
  # Calculate ICER
  if (inc_qaly > 0) {
    icer <- inc_cost / inc_qaly
  } else if (inc_cost < 0) {
    icer <- NA  # Dominant (more effective, less costly)
  } else {
    icer <- NA  # Dominated (less effective, more costly)
  }
  
  # Store results
  psa_results$cost_combo[i] <- total_cost_combo
  psa_results$qaly_combo[i] <- total_qaly_combo
  psa_results$cost_nivo[i] <- total_cost_nivo
  psa_results$qaly_nivo[i] <- total_qaly_nivo
  psa_results$inc_cost[i] <- inc_cost
  psa_results$inc_qaly[i] <- inc_qaly
  psa_results$icer[i] <- icer
  
  # Progress indicator
  if (i %% 100 == 0) {
    cat(sprintf("Completed %d of %d simulations\n", i, n_simulations))
  }
}

# 7) Analyze and visualize PSA results
########################################################################

# Summary statistics
psa_summary <- data.frame(
  Metric = c("Cost Combo", "QALY Combo", "Cost Nivo", "QALY Nivo", "Inc Cost", "Inc QALY", "ICER"),
  Mean = c(
    mean(psa_results$cost_combo), 
    mean(psa_results$qaly_combo),
    mean(psa_results$cost_nivo),
    mean(psa_results$qaly_nivo),
    mean(psa_results$inc_cost),
    mean(psa_results$inc_qaly),
    mean(psa_results$icer, na.rm = TRUE)
  ),
  Median = c(
    median(psa_results$cost_combo), 
    median(psa_results$qaly_combo),
    median(psa_results$cost_nivo),
    median(psa_results$qaly_nivo),
    median(psa_results$inc_cost),
    median(psa_results$inc_qaly),
    median(psa_results$icer, na.rm = TRUE)
  ),
  LCI = c(
    quantile(psa_results$cost_combo, 0.025), 
    quantile(psa_results$qaly_combo, 0.025),
    quantile(psa_results$cost_nivo, 0.025),
    quantile(psa_results$qaly_nivo, 0.025),
    quantile(psa_results$inc_cost, 0.025),
    quantile(psa_results$inc_qaly, 0.025),
    quantile(psa_results$icer, 0.025, na.rm = TRUE)
  ),
  UCI = c(
    quantile(psa_results$cost_combo, 0.975), 
    quantile(psa_results$qaly_combo, 0.975),
    quantile(psa_results$cost_nivo, 0.975),
    quantile(psa_results$qaly_nivo, 0.975),
    quantile(psa_results$inc_cost, 0.975),
    quantile(psa_results$inc_qaly, 0.975),
    quantile(psa_results$icer, 0.975, na.rm = TRUE)
  )
)

print(psa_summary)

# Calculate CE probabilities at different thresholds
calculate_ce_probability <- function(psa_results, threshold) {
  n_dominant <- sum(psa_results$inc_qaly > 0 & psa_results$inc_cost < 0, na.rm = TRUE)
  n_cost_effective <- sum(psa_results$inc_qaly > 0 & psa_results$inc_cost > 0 & 
                            psa_results$inc_cost / psa_results$inc_qaly < threshold, na.rm = TRUE)
  prob <- (n_dominant + n_cost_effective) / nrow(psa_results)
  return(prob)
}

thresholds <- c(50000, 100000, 150000, 200000)
ce_probs <- sapply(thresholds, function(t) calculate_ce_probability(psa_results, t))
ce_table <- data.frame(
  Threshold = thresholds,
  Probability = ce_probs
)
print(ce_table)

# Create Cost-Effectiveness Plane
ce_plane <- ggplot(psa_results, aes(x = inc_qaly, y = inc_cost)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 100000, intercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Cost-Effectiveness Plane",
    x = "Incremental QALYs (Combo vs. Nivo)",
    y = "Incremental Cost (Combo vs. Nivo)",
    caption = "Red dotted line: $100,000/QALY threshold"
  ) +
  theme_minimal()

print(ce_plane)

# Create Cost-Effectiveness Acceptability Curve (CEAC)
wtp_values <- seq(0, 300000, by = 10000)
ceac_probs <- sapply(wtp_values, function(t) calculate_ce_probability(psa_results, t))

ceac_df <- data.frame(
  WTP = wtp_values,
  Probability = ceac_probs
)

ceac_plot <- ggplot(ceac_df, aes(x = WTP, y = Probability)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(labels = scales::dollar_format()) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Cost-Effectiveness Acceptability Curve",
    x = "Willingness to Pay Threshold ($/QALY)",
    y = "Probability Cost-Effective"
  ) +
  theme_minimal()

print(ceac_plot)

# Save results
save(psa_results, file = "psa_results.RData")
ggsave("ce_plane.png", ce_plane, width = 8, height = 6)
ggsave("ceac_plot.png", ceac_plot, width = 8, height = 6)