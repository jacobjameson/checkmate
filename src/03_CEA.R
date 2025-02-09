

##################################################################
# 0) Setup and Key Parameters
##################################################################
library(survival)
library(flexsurv)
library(dplyr)

# We assume you have two survival models:
#   best_os_model_ni_ni  (flexsurvreg or flexsurvspline object for OS)
#   best_pfs_model_ni (flexsurvreg or flexsurvspline object for PFS)
#
# Example: best_os_model_ni <- ...
#          best_pfs_model_ni <- ...

# Markov model states: Stable, Progressed, Dead
state_names <- c("Stable","Progressed","Dead")

# 10-year horizon, monthly cycle => 120 cycles
max_months <- 120
cycle_length <- 1   # 1 month per cycle
n_cycles <- max_months / cycle_length

# Annual discount rate = 3%, so monthly discount factor = (1+0.03)^( - month/12 ).
disc_rate_annual <- 0.03

# We'll construct a vector of months from 0..120
time_vec <- seq(0, max_months, by = cycle_length)

##################################################################
# 1) Transition Probability Functions from OS/PFS
##################################################################
# We'll replicate the typical "ratio" approach. Each month 'm':
#   stable->dead:  1 - (S_os[m+1]/S_os[m])
#   stable->stable:  (S_pfs[m+1]/S_pfs[m])
#   stable->progressed: remainder
#   progressed->dead:   1 - ( (S_os[m+1]-S_pfs[m+1]) / (S_os[m]-S_pfs[m]) )
#   progressed->progressed: remainder
#
# We'll define a helper to get OS and PFS at these discrete months.

get_surv <- function(model, times) {
  # summary(..., type="survival") => list of length(times)
  s_list <- summary(model, t= times, type = "survival")
  sapply(s_list, function(x) x$est)
}


os_surv <- get_surv(best_os_model_ni,  time_vec)  # length=121
pfs_surv <- get_surv(best_pfs_model_ni, time_vec)

##################################################################
# 2) COSTS & UTILITIES (Base Case)
##################################################################
# As per the paper, they used monthly costs (2017 USD) from Table 1,
# plus monthly utilities. We pick some placeholders. 
# In your real code, fill in actual numbers from your table(s).

### Example: Standard therapy arm
cost_std_drug       <- 3000   # monthly cost of standard drug
cost_std_infusion   <- 136    # from the paper, e.g. "Chemo/immunotherapy infusion"
cost_std_toxicity   <- 4947   # from Table 1 for standard therapy toxicity
prob_std_toxicity   <- 0.15   # example probability of G3-G4 toxicity each month
util_std_tox_penalty<- 0.0069 # monthly utility penalty if toxicity occurs
# Utility stable / progressed per month:
util_stable_m <- 0.517  # from Table 1 (0.517 per month if stable)
util_prog_m   <- 0.280  # from Table 1
cost_stable_m <- 400    # from the example or your table
cost_prog_m    <- 100
cost_dead_m    <- 0
# palliative cost if you die that month (in the last month):
cost_palliative <- 9836  # from Table 1

### Example: Nivolumab arm
cost_nivo_drug       <- 13432  # monthly cost from Table 1
cost_nivo_infusion   <- 136
cost_nivo_toxicity   <- 1677   # from Table 1 for nivolumab toxicity
prob_nivo_toxicity   <- 0.12   # placeholder
util_nivo_tox_penalty<- 0.0070

##################################################################
# 3) Function to run Markov for a single "arm"
##################################################################
run_markov_one_arm <- function(os_vec, pfs_vec,
                               cost_drug, cost_infusion, cost_toxicity,
                               prob_tox_month, util_tox_penalty,
                               cost_stable, cost_prog, cost_dead,
                               cost_palliative,
                               util_stable, util_prog, util_dead=0,
                               disc_rate=0.03) {
  # os_vec, pfs_vec are length(n_cycles+1) survival probabilities at each month
  # We'll store Markov trace as a matrix (n_cycles+1) x 3 states
  trace <- matrix(0, nrow=n_cycles+1, ncol=3, dimnames=list(0:n_cycles, state_names))
  # Everyone starts in stable at t=0
  trace[1,] <- c(1,0,0)
  
  # Store cost & QALYs each cycle
  cycle_cost <- numeric(n_cycles)
  cycle_qaly <- numeric(n_cycles)
  
  for (m in 1:n_cycles) {
    # OS/PFS at start month m
    S_os_t   <- os_vec[m]   
    S_os_t1  <- os_vec[m+1] 
    S_pfs_t  <- pfs_vec[m]
    S_pfs_t1 <- pfs_vec[m+1]
    
    # stable->dead
    p_SD_D <- ifelse(S_os_t>0, 1 - S_os_t1/S_os_t, 0)
    # stable->stable
    p_SD_SD <- ifelse(S_pfs_t>0, S_pfs_t1/S_pfs_t, 0)
    # stable->progressed
    p_SD_PD <- 1 - p_SD_SD - p_SD_D
    if (p_SD_PD<0) p_SD_PD <- 0
    
    # progressed->dead
    denom_prog <- (S_os_t - S_pfs_t)
    numer_prog <- (S_os_t1 - S_pfs_t1)
    ratio_prog_alive <- ifelse(denom_prog>0, numer_prog/denom_prog, 0)
    p_PD_D <- 1 - ratio_prog_alive
    p_PD_PD <- ratio_prog_alive
    
    # We'll get the distribution at cycle m
    dist_m <- trace[m,]
    # Next cycle distribution
    dist_mplus1 <- numeric(3)
    # stable
    dist_mplus1[1] <- dist_m[1]*p_SD_SD
    dist_mplus1[2] <- dist_m[1]*p_SD_PD + dist_m[2]*p_PD_PD
    dist_mplus1[3] <- dist_m[1]*p_SD_D   + dist_m[2]*p_PD_D + dist_m[3]
    
    trace[m+1,] <- dist_mplus1
    
    # Now compute cost and QALYs for cycle m
    # We consider the "average" fraction in each state during this cycle.
    # A simpler approach is to assume the entire cycle is spent in the distribution at the start of the cycle.
    frac_stable    <- dist_m[1]
    frac_progress  <- dist_m[2]
    frac_dead      <- dist_m[3]
    
    # Because the therapy is given if you're alive (stable or progressed):
    # cost of drug + infusion if alive
    cycle_cost_alive <- (frac_stable + frac_progress)*(cost_drug + cost_infusion)
    # cost of disease management
    cycle_cost_disease <- frac_stable*cost_stable + frac_progress*cost_prog + frac_dead*cost_dead
    # cost of palliative if they die *this cycle*? 
    # Strictly we'd do dist_mplus1[3] - dist_m[3], i.e. new deaths. Let's do a quick approximation:
    newly_dead <- dist_mplus1[3] - dist_m[3]
    if(newly_dead<0) newly_dead <- 0
    cycle_cost_palliative <- newly_dead * cost_palliative
    
    # Toxicity cost: if you are alive, there's a prob(prob_tox_month) of G3-G4 toxicity => cost_toxicity
    # The expected cost of tox = cost_toxicity * (frac_alive * prob_tox_month)
    frac_alive <- frac_stable + frac_progress
    cycle_cost_tox <- cost_toxicity * (frac_alive * prob_tox_month)
    
    # sum cost
    cost_m <- cycle_cost_alive + cycle_cost_disease + cycle_cost_palliative + cycle_cost_tox
    
    # QALYs: 
    #   stable => 0.517 monthly? or 0.517 is a *yearly* utility? 
    #   Actually the paper states 0.517 as "per year." 
    #   For a 1-month cycle, we do 0.517/12. We'll do that conversion in the code:
    util_stable_cycle <- util_stable/12
    util_prog_cycle   <- util_prog/12
    # Death = 0
    # Then subtract toxicity penalty if it occurs: 
    #   If you have an AE, you lose "util_tox_penalty" for that cycle. 
    #   expected AE = frac_alive * prob_tox_month => apply that penalty in expectation
    #   We do it *in addition* or we subtract it from the stable/prog utility?
    #   The paper subtracts an absolute decrement from stable or progressed. We'll do:
    
    # baseline utility
    qaly_base <- frac_stable*util_stable_cycle + frac_progress*util_prog_cycle
    # toxicity penalty
    tox_penalty <- (frac_alive * prob_tox_month)*util_tox_penalty
    
    qaly_m <- qaly_base - tox_penalty
    
    # discount for month m => discount factor = 1/(1+disc_rate)^(m/12)
    # approximate using continuous or discrete. We'll do discrete:
    month_in_years <- m/12
    discount_factor <- (1+disc_rate)^( - month_in_years )
    
    cycle_cost[m] <- cost_m * discount_factor
    cycle_qaly[m] <- qaly_m * discount_factor
  }
  
  total_cost <- sum(cycle_cost)
  total_qalys <- sum(cycle_qaly)
  
  list(
    trace = trace,
    total_cost = total_cost,
    total_qalys = total_qalys
  )
}

##################################################################
# 4) Run the Model for Two Arms: (1) Standard Therapy, (2) Nivolumab
##################################################################
res_std <- run_markov_one_arm(
  os_vec = os_surv, 
  pfs_vec= pfs_surv,
  cost_drug = cost_std_drug,
  cost_infusion = cost_std_infusion,
  cost_toxicity = cost_std_toxicity,
  prob_tox_month= prob_std_toxicity,
  util_tox_penalty= util_std_tox_penalty,
  cost_stable = cost_stable_m,
  cost_prog   = cost_prog_m,
  cost_dead   = cost_dead_m,
  cost_palliative = cost_palliative,
  util_stable = util_stable_m,
  util_prog   = util_prog_m,
  util_dead   = 0,
  disc_rate   = disc_rate_annual
)

res_nivo <- run_markov_one_arm(
  os_vec = os_surv,
  pfs_vec= pfs_surv,
  cost_drug = cost_nivo_drug,
  cost_infusion = cost_nivo_infusion,
  cost_toxicity = cost_nivo_toxicity,
  prob_tox_month= prob_nivo_toxicity,
  util_tox_penalty= util_nivo_tox_penalty,
  cost_stable = cost_stable_m,
  cost_prog   = cost_prog_m,
  cost_dead   = cost_dead_m,
  cost_palliative = cost_palliative,
  util_stable = util_stable_m,
  util_prog   = util_prog_m,
  disc_rate   = disc_rate_annual
)

# Summaries
cat("\n--- Standard Therapy ---\n")
cat(sprintf("Total Cost:  $%0.2f\n", res_std$total_cost))
cat(sprintf("Total QALYs: %0.3f\n", res_std$total_qalys))

cat("\n--- Nivolumab ---\n")
cat(sprintf("Total Cost:  $%0.2f\n", res_nivo$total_cost))
cat(sprintf("Total QALYs: %0.3f\n", res_nivo$total_qalys))

# Incremental
inc_cost  <- res_nivo$total_cost  - res_std$total_cost
inc_qalys <- res_nivo$total_qalys - res_std$total_qalys
icer <- inc_cost / inc_qalys

cat("\n--- Incremental (Nivolumab vs. Standard) ---\n")
cat(sprintf("Incremental Cost:  $%0.2f\n", inc_cost))
cat(sprintf("Incremental QALYs: %0.3f\n", inc_qalys))
cat(sprintf("ICER ($/QALY):     $%0.2f\n", icer))
