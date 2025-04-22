############################################################################
# PARTITIONED SURVIVAL MODEL
# Using updated parameters & constant monthly AE penalty
############################################################################

library(flexsurv)

# 1) Time Horizon & Discounting
############################################################################
max_months <- 120    # 10 years
cycle_length <- 1
time_vec <- seq(0, max_months, by = cycle_length)

disc_rate <- 0.03

discount_factor <- function(m) {
  # discrete monthly discount factor for cycle m
  (1 + disc_rate)^(-(m / 12))
}

# 2) Survival Curves (from your best-fit parametric/spline models)
############################################################################
# The user has:
#   best_os_model_ni, best_pfs_model_ni  => NIVO+IPI
#   best_os_model_n,  best_pfs_model_n   => NIVO alone

get_surv <- function(model, times) {
  # Summaries survival from a flexsurv model
  s_list <- summary(model, t = times, type = "survival")
  sapply(s_list, function(x) x$est)
}

# OS/PFS for combination
os_ni  <- get_surv(best_os_model_ni,  time_vec)
pfs_ni <- get_surv(best_pfs_model_ni, time_vec)

# OS/PFS for NIVO alone
os_n  <- get_surv(best_os_model_n,  time_vec)
pfs_n <- get_surv(best_pfs_model_n, time_vec)

# 3) Partition: Stable, Progressed, Dead
############################################################################
stable_ni   <- pfs_ni
progress_ni <- os_ni - pfs_ni
dead_ni     <- 1 - os_ni

stable_n    <- pfs_n
progress_n  <- os_n - pfs_n
dead_n      <- 1 - os_n

# 4) Parameter Inputs (from your updated table)
############################################################################
# DRUG COSTS (per cycle)
cost_nivolumab       <- 14975   # from table,  (95% CI: 9703–21417)
cost_combo           <- 26425   # from table, (17089–37662)

# For an up-to-24-month time-limited therapy
max_therapy_months   <- 24

# SECOND-LINE TREATMENT (progressed) monthly or lumpsum?
# The table shows $8,908 for second-line after combo, $12,093 after NIVO alone.
# We'll treat them as monthly costs while progressed, or lumpsum?
# Often, it is simpler to treat them as the *progressed state cost* if we assume
# all progressed eventually get second-line. We'll do that here:

cost_prog_combo      <- 8908
cost_prog_nivo       <- 12093

# DRUG-TOXICITY Disutilities
# Weighted monthly penalty: 
# e.g. 0.012 for combo, 0.0036 for NIVO
ae_penalty_combo     <- 0.012
ae_penalty_nivo      <- 0.0036

# DISEASE COSTS
cost_stable          <- 2166   # e.g. from table, 95% CI: 1397–3098
cost_palliative_death<- 15957  # one-time lumpsum if you want to charge at death
# We'll incorporate that below.

# UTILITY
# stable disease utility: 0.754, progressed decrement: 0.180 => 0.754-0.180=0.574?
util_stable_base     <- 0.754
util_prog_decrement  <- 0.180
util_prog            <- 0.754 - 0.180   # =0.574
util_death           <- 0.0

# 5) Partitioned Survival Function
############################################################################
# Each cycle i, we compute costs/QALYs for stable & progressed.
# - If stable & i<=24 => apply drug cost (combo or nivo)
# - Add disease cost (2,166 stable, progressed?)
# - For progressed, we might add second-line cost
# - Utility for stable or progressed
# - Add AE penalty each cycle for stable *only if i<=24*

calc_partition_arm <- function(stable_vec, prog_vec, dead_vec,
                               arm = "Combo") {
  
  n_cycles <- length(stable_vec)
  monthly_cost <- numeric(n_cycles)
  monthly_qaly <- numeric(n_cycles)
  
  for (m in seq_len(n_cycles)) {
    # distribution at month m
    fs <- stable_vec[m]
    fp <- prog_vec[m]
    
    # 1) DRUG COST
    if (arm == "Combo") {
      # stable & within 24 months => cost_combo
      if (m <= max_therapy_months) {
        drug_c <- fs * cost_combo
      } else {
        drug_c <- 0
      }
      # progressed => cost_prog_combo
      prog_c <- fp * cost_prog_combo
      
    } else { # "NivoAlone"
      if (m <= max_therapy_months) {
        drug_c <- fs * cost_nivolumab
      } else {
        drug_c <- 0
      }
      # progressed => cost_prog_nivo
      prog_c <- fp * cost_prog_nivo
    }
    
    # 2) DISEASE COST
    # stable = cost_stable, progressed => NA? Or do we keep a separate cost_stable?
    # The table suggests stable=2,166. The second-line might be "progressed cost"
    # We'll incorporate that 2,166 for stable, and 4,000 is not shown in your snippet now. 
    # But let's assume the second-line cost is the entire cost for progressed? We'll do that.
    stable_c <- fs * cost_stable
    
    cost_m <- drug_c + prog_c + stable_c
    
    # 3) UTILITY
    # stable => 0.754, progressed => 0.574
    # But we also have AE penalty each cycle for stable if on therapy
    if (arm == "Combo") {
      # stable utility w/ AE penalty if m <=24
      if (m <= max_therapy_months) {
        u_stable <- util_stable_base - ae_penalty_combo
      } else {
        u_stable <- util_stable_base
      }
    } else { # "NivoAlone"
      if (m <= max_therapy_months) {
        u_stable <- util_stable_base - ae_penalty_nivo
      } else {
        u_stable <- util_stable_base
      }
    }
    # progressed utility
    u_prog <- util_stable_base - util_prog_decrement
    # total utility
    util_m <- fs * u_stable + fp * u_prog
    
    # discount
    dfact <- discount_factor(m-1)
    monthly_cost[m] <- cost_m * dfact
    monthly_qaly[m] <- util_m * (1/12) * dfact
  }
  
  # sum cost & qaly
  total_cost <- sum(monthly_cost)
  total_qaly <- sum(monthly_qaly)
  
  # 4) OPTIONAL: PALLIATIVE CARE AT DEATH
  # If we want a lumpsum for each death event, we can approximate:
  # newly dead at cycle m is (dead_vec[m] - dead_vec[m-1]) if m>1
  # let's do a quick approach
  nd_sum <- 0
  for(m in 2:n_cycles) {
    newly_dead <- (dead_vec[m] - dead_vec[m-1])
    if(newly_dead < 0) newly_dead <- 0
    # discount cost at midpoint? We'll discount by cycle factor
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead*cost_palliative_death*dfact
  }
  
  total_cost <- total_cost + nd_sum
  
  return(list(arm=arm,
              total_cost=total_cost,
              total_qaly=total_qaly))
}


# 6) RUN CALCULATIONS
############################################################################
res_combo <- calc_partition_arm(stable_ni, progress_ni, dead_ni, arm="Combo")
res_nivo  <- calc_partition_arm(stable_n,  progress_n,  dead_n,  arm="NivoAlone")

cat("=== PARTITIONED SURVIVAL RESULTS (BASE CASE) ===\n\n")

cat("Nivolumab + Ipilimumab:\n")
cat(sprintf("  Total Cost:  $%.2f\n", res_combo$total_cost))
cat(sprintf("  Total QALYs: %.3f\n",  res_combo$total_qaly))

cat("\nNivolumab Alone:\n")
cat(sprintf("  Total Cost:  $%.2f\n", res_nivo$total_cost))
cat(sprintf("  Total QALYs: %.3f\n",  res_nivo$total_qaly))

icost <- res_combo$total_cost - res_nivo$total_cost
iqaly <- res_combo$total_qaly - res_nivo$total_qaly
cat("\nIncremental (Combo vs. Nivo):\n")
cat(sprintf("  Inc Cost: $%.2f\n",  icost))
cat(sprintf("  Inc QALY: %.3f\n",  iqaly))
if(iqaly>0) {
  cat(sprintf("  ICER:     $%.2f per QALY\n", icost/iqaly))
} else {
  cat("  Combo is dominated (IQALY <= 0)\n")
}