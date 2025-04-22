########################################################################
# PARTITIONED SURVIVAL MODEL FOR ADVANCED MELANOMA
# Based on CheckMate 067 trial data with 10-year follow-up
########################################################################

library(flexsurv)

# 1) Time Horizon & Discounting
########################################################################
max_months <- 120    # 10 years (matches trial follow-up)
cycle_length <- 1    # Monthly cycles
time_vec <- seq(0, max_months, by = cycle_length)
disc_rate <- 0.03    # Standard 3% annual discount rate
discount_factor <- function(m) {
  # Discrete monthly discount factor for cycle m
  (1 + disc_rate)^(-(m / 12))
}

# 2) Survival Curves (from your best-fit parametric models)
########################################################################
# The user has:
#   best_os_model_ni, best_pfs_model_ni  => NIVO+IPI (Combination therapy)
#   best_os_model_n,  best_pfs_model_n   => NIVO alone

get_surv <- function(model, times) {
  # Summarize survival from a flexsurv model
  s_list <- summary(model, t = times, type = "survival")
  sapply(s_list, function(x) x$est)
}

# OS/PFS for combination
os_ni  <- get_surv(best_os_model_ni, time_vec)
pfs_ni <- get_surv(best_pfs_model_ni, time_vec)

# OS/PFS for NIVO alone
os_n  <- get_surv(best_os_model_n, time_vec)
pfs_n <- get_surv(best_pfs_model_n, time_vec)

# 3) Partition: Stable, Progressed, Dead
########################################################################
# Ensure no negative values in state proportions
stable_ni   <- pfs_ni
progress_ni <- pmax(os_ni - pfs_ni, 0)  # Prevent negative values
dead_ni     <- 1 - os_ni

stable_n    <- pfs_n
progress_n  <- pmax(os_n - pfs_n, 0)    # Prevent negative values
dead_n      <- 1 - os_n

# 4) Parameter Inputs (based on trial data and literature)
########################################################################
# DRUG COSTS (per cycle)
cost_nivolumab       <- 14975   # from table, (95% CI: 9703–21417)
cost_combo           <- 26425   # from table, (17089–37662)

# TREATMENT DURATION
# Based on the trial protocol that continued until progression/toxicity
# and common clinical practice of limiting treatment to 2 years maximum
max_therapy_months   <- 24

# DISCONTINUATION DUE TO ADVERSE EVENTS
# From trial: 44.4% of combo and 15.7% of nivolumab patients discontinued due to toxicity
disc_ae_combo        <- 0.444   # Proportion discontinuing due to AEs (combo)
disc_ae_nivo         <- 0.157   # Proportion discontinuing due to AEs (nivolumab) 
disc_ae_period       <- 3       # Assume most AE discontinuations occur in first 3 months

# SUBSEQUENT TREATMENT (from trial data)
# 36% of combo patients received subsequent systemic therapy
# 50% of nivolumab patients received subsequent systemic therapy
prop_2ndline_combo   <- 0.36    # From trial
prop_2ndline_nivo    <- 0.50    # From trial

# Base second-line costs (to be adjusted by proportion receiving)
cost_prog_combo_base <- 8908
cost_prog_nivo_base  <- 12093

# Adjusted second-line costs (applied to all progressed patients)
cost_prog_combo      <- cost_prog_combo_base * prop_2ndline_combo
cost_prog_nivo       <- cost_prog_nivo_base * prop_2ndline_nivo

# SUBSEQUENT LOCAL THERAPY (from trial data)
# 46% of combo and 56% of nivolumab patients received local therapy
prop_local_combo     <- 0.46    # From trial
prop_local_nivo      <- 0.56    # From trial
cost_local_therapy   <- 5000    # Estimated cost for local therapy (e.g., radiation, surgery)

# Apply local therapy as one-time cost upon progression (see implementation in main function)

# ADVERSE EVENT DISUTILITIES
# From trial: Weighted average disutilities based on AE incidence and duration
ae_penalty_combo     <- 0.012
ae_penalty_nivo      <- 0.0036

# AE TIMING (based on trial data showing front-loaded AE occurrence)
ae_early_months      <- 3       # First 3 months: highest impact
ae_medium_months     <- 6       # Months 4-6: medium impact

# DISEASE COSTS
cost_stable          <- 2166    # From table, (95% CI: 1397–3098)
cost_palliative_death<- 15957   # One-time lumpsum cost at death

# UTILITY
util_stable_base     <- 0.754   # Base utility for stable disease
util_prog_decrement  <- 0.180   # Utility decrement for progression
util_prog            <- util_stable_base - util_prog_decrement  # = 0.574
util_death           <- 0.0

# 5) Partitioned Survival Function
########################################################################
calc_partition_arm <- function(stable_vec, prog_vec, dead_vec, arm = "Combo") {
  n_cycles <- length(stable_vec)
  monthly_cost <- numeric(n_cycles)
  monthly_qaly <- numeric(n_cycles)
  
  # Track first progression events for local therapy costs
  prog_events <- numeric(n_cycles)
  if (n_cycles > 1) {
    for (m in 2:n_cycles) {
      # New progressions = decrease in stable population - deaths from stable
      # This is an approximation assuming stable patients don't die directly
      prog_events[m] <- max(0, (stable_vec[m-1] - stable_vec[m]) - 
                              (dead_vec[m] - dead_vec[m-1]) * (stable_vec[m-1]/(stable_vec[m-1] + prog_vec[m-1])))
    }
  }
  
  for (m in seq_len(n_cycles)) {
    # Distribution at month m
    fs <- stable_vec[m]
    fp <- prog_vec[m]
    
    # 1) DRUG COST - with discontinuation due to AEs
    if (arm == "Combo") {
      if (m <= max_therapy_months) {
        if (m <= disc_ae_period) {
          # In first 3 months, apply monthly discontinuation rate
          monthly_disc_rate <- disc_ae_combo / disc_ae_period
          # Proportion still on therapy after accounting for AE discontinuations
          prop_on_therapy <- 1.0 - (monthly_disc_rate * (m-1))
          # Apply cost only to those still on therapy
          drug_c <- fs * cost_combo * prop_on_therapy
        } else {
          # After disc_ae_period, all remaining patients continue to max_therapy_months
          prop_remaining <- 1.0 - disc_ae_combo
          drug_c <- fs * cost_combo * prop_remaining
        }
      } else {
        drug_c <- 0
      }
      
      # Progressed => cost_prog_combo (already adjusted for proportion receiving)
      prog_c <- fp * cost_prog_combo
      
    } else { # "NivoAlone"
      if (m <= max_therapy_months) {
        if (m <= disc_ae_period) {
          monthly_disc_rate <- disc_ae_nivo / disc_ae_period
          prop_on_therapy <- 1.0 - (monthly_disc_rate * (m-1))
          drug_c <- fs * cost_nivolumab * prop_on_therapy
        } else {
          prop_remaining <- 1.0 - disc_ae_nivo
          drug_c <- fs * cost_nivolumab * prop_remaining
        }
      } else {
        drug_c <- 0
      }
      
      # Progressed => cost_prog_nivo (already adjusted for proportion receiving)
      prog_c <- fp * cost_prog_nivo
    }
    
    # 2) DISEASE COST
    stable_c <- fs * cost_stable
    
    # 3) LOCAL THERAPY COSTS (one-time at progression)
    local_c <- 0
    if (m > 1) {
      if (arm == "Combo") {
        local_c <- prog_events[m] * cost_local_therapy * prop_local_combo
      } else {
        local_c <- prog_events[m] * cost_local_therapy * prop_local_nivo
      }
    }
    
    # Total costs for this cycle
    cost_m <- drug_c + prog_c + stable_c + local_c
    
    # 4) UTILITY WITH TIME-DEPENDENT AE IMPACT
    if (arm == "Combo") {
      if (m <= max_therapy_months) {
        if (m <= ae_early_months) {
          # First 3 months: full impact (most AEs occur and are most severe)
          ae_factor <- 1.0
        } else if (m <= ae_medium_months) {
          # Months 4-6: reduced impact (some AEs resolved, others improving)
          ae_factor <- 0.6
        } else {
          # Months 7-24: minimal impact (mainly endocrine AEs persist)
          ae_factor <- 0.2
        }
        
        u_stable <- util_stable_base - (ae_penalty_combo * ae_factor)
      } else {
        u_stable <- util_stable_base
      }
    } else { # "NivoAlone"
      if (m <= max_therapy_months) {
        if (m <= ae_early_months) {
          ae_factor <- 1.0
        } else if (m <= ae_medium_months) {
          ae_factor <- 0.6
        } else {
          ae_factor <- 0.2
        }
        
        u_stable <- util_stable_base - (ae_penalty_nivo * ae_factor)
      } else {
        u_stable <- util_stable_base
      }
    }
    
    # Progressed utility (same for both arms)
    u_prog <- util_prog
    
    # Total utility for this cycle
    util_m <- fs * u_stable + fp * u_prog
    
    # Apply discounting
    dfact <- discount_factor(m-1)
    monthly_cost[m] <- cost_m * dfact
    monthly_qaly[m] <- util_m * (1/12) * dfact  # Convert to QALYs (1/12 of a year)
  }
  
  # Sum costs & QALYs
  total_cost <- sum(monthly_cost)
  total_qaly <- sum(monthly_qaly)
  
  # 5) PALLIATIVE CARE AT DEATH
  # Calculate newly dead at each cycle
  nd_sum <- 0
  for(m in 2:n_cycles) {
    newly_dead <- (dead_vec[m] - dead_vec[m-1])
    if(newly_dead < 0) newly_dead <- 0  # Safety check
    # Discount cost at cycle midpoint
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead * cost_palliative_death * dfact
  }
  
  total_cost <- total_cost + nd_sum
  
  # Return results
  return(list(
    arm = arm,
    total_cost = total_cost,
    total_qaly = total_qaly,
    monthly_cost = monthly_cost,
    monthly_qaly = monthly_qaly
  ))
}

# 6) Model Validation Function
########################################################################
validate_model <- function(os_model, pfs_model, observed_os, observed_pfs, timepoints) {
  # Get predicted survival at key timepoints
  predicted_os <- summary(os_model, t = timepoints, type = "survival")
  predicted_os <- sapply(predicted_os, function(x) x$est)
  
  predicted_pfs <- summary(pfs_model, t = timepoints, type = "survival")
  predicted_pfs <- sapply(predicted_pfs, function(x) x$est)
  
  # Create validation table
  validation <- data.frame(
    Timepoint = timepoints,
    OS_Predicted = predicted_os,
    OS_Observed = observed_os,
    OS_Diff = predicted_os - observed_os,
    PFS_Predicted = predicted_pfs,
    PFS_Observed = observed_pfs,
    PFS_Diff = predicted_pfs - observed_pfs
  )
  
  return(validation)
}

# 7) RUN CALCULATIONS
########################################################################
res_combo <- calc_partition_arm(stable_ni, progress_ni, dead_ni, arm = "Combo")
res_nivo <- calc_partition_arm(stable_n, progress_n, dead_n, arm = "NivoAlone")

# 8) RESULTS REPORTING
########################################################################
print_results <- function(res_combo, res_nivo) {
  cat("=== PARTITIONED SURVIVAL RESULTS (BASE CASE) ===\n\n")
  
  cat("Nivolumab + Ipilimumab:\n")
  cat(sprintf("  Total Cost:  $%.2f\n", res_combo$total_cost))
  cat(sprintf("  Total QALYs: %.3f\n", res_combo$total_qaly))
  
  cat("\nNivolumab Alone:\n")
  cat(sprintf("  Total Cost:  $%.2f\n", res_nivo$total_cost))
  cat(sprintf("  Total QALYs: %.3f\n", res_nivo$total_qaly))
  
  icost <- res_combo$total_cost - res_nivo$total_cost
  iqaly <- res_combo$total_qaly - res_nivo$total_qaly
  
  cat("\nIncremental (Combo vs. Nivo):\n")
  cat(sprintf("  Inc Cost: $%.2f\n", icost))
  cat(sprintf("  Inc QALY: %.3f\n", iqaly))
  
  if(iqaly > 0) {
    cat(sprintf("  ICER:     $%.2f per QALY\n", icost/iqaly))
  } else {
    cat("  Combo is dominated (IQALY <= 0)\n")
  }
  
  # Additional model statistics
  cat("\nKey Model Statistics:\n")
  cat(sprintf("  Time Horizon: %d months (%.1f years)\n", max_months, max_months/12))
  cat(sprintf("  Discount Rate: %.1f%%\n", disc_rate*100))
  cat(sprintf("  Undiscounted Life Years (Combo): %.2f\n", sum(os_ni)/12))
  cat(sprintf("  Undiscounted Life Years (Nivo): %.2f\n", sum(os_n)/12))
  cat(sprintf("  10-year OS (Combo): %.1f%%\n", os_ni[121]*100))
  cat(sprintf("  10-year OS (Nivo): %.1f%%\n", os_n[121]*100))
}

# Print results
print_results(res_combo, res_nivo)

# 9) MODEL CALIBRATION
########################################################################
# Observed survival data from CheckMate 067 trial (replace with actual values)
# Key timepoints in months: 12, 24, 36, 60, 120 (10 years)
timepoints <- c(12, 24, 36, 60, 120)

# From the trial paper: Overall survival at 10 years was 43% (combo), 37% (nivo)
os_observed_combo <- c(0.73, 0.64, 0.58, 0.52, 0.43)  # Example - replace with actual data
os_observed_nivo <- c(0.71, 0.59, 0.52, 0.44, 0.37)   # Example - replace with actual data

pfs_observed_combo <- c(0.50, 0.43, 0.39, 0.36, 0.33) # Example - replace with actual data
pfs_observed_nivo <- c(0.37, 0.32, 0.29, 0.28, 0.27)  # Example - replace with actual data

# Validate models
validation_combo <- validate_model(best_os_model_ni, best_pfs_model_ni, 
                                   os_observed_combo, pfs_observed_combo, timepoints)
validation_nivo <- validate_model(best_os_model_n, best_pfs_model_n, 
                                  os_observed_nivo, pfs_observed_nivo, timepoints)

# Optionally, print validation results
cat("\nModel Calibration (Combo):\n")
print(validation_combo)

cat("\nModel Calibration (Nivo):\n")
print(validation_nivo)

# 10) SAVE RESULTS
########################################################################
# Optionally save detailed results for further analysis
results <- list(
  combo = res_combo,
  nivo = res_nivo,
  icost = res_combo$total_cost - res_nivo$total_cost,
  iqaly = res_combo$total_qaly - res_nivo$total_qaly,
  icer = (res_combo$total_cost - res_nivo$total_cost) / 
    (res_combo$total_qaly - res_nivo$total_qaly),
  parameters = list(
    max_months = max_months,
    disc_rate = disc_rate,
    max_therapy_months = max_therapy_months,
    prop_2ndline_combo = prop_2ndline_combo,
    prop_2ndline_nivo = prop_2ndline_nivo,
    disc_ae_combo = disc_ae_combo,
    disc_ae_nivo = disc_ae_nivo
  ),
  validation_combo = validation_combo,
  validation_nivo = validation_nivo
)

# Uncomment to save results
# save(results, file = "psa_base_case_results.RData")