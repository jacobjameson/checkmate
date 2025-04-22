#########################################################
# 1) Prepare the four OS/PFS curves for subgroups
#########################################################

# "Truly Responder on single-agent" => same survival as combination, 
#   but single-agent cost. We'll store them in vectors:
os_singleResp  <- os_comb  # or just define a copy
pfs_singleResp <- pfs_comb

# "Truly Non-responder on single-agent" => your derived "resistant" curves:
os_resistant  
pfs_resistant 
  
  # We'll re-use os_comb,pfs_comb for "Truly R + combo" and "Truly NR + combo" 
  # We'll define cost/utility parameters for single-agent vs. combo.
  
  #########################################################
# 2) Partition Surv function
#########################################################
# We'll define a simpler function:
partition_ce <- function(os_vec, pfs_vec, 
                         cost_stable, cost_progressed,
                         drug_cost, # cost of therapy while stable & <=24 mo
                         ae_penalty, # monthly AE penalty while on therapy
                         base_stable_util, base_prog_util,
                         horizon=120) {
  # returns a list(cost=..., qaly=...)
  
  # derive stable/prog/dead from OS/PFS
  stable_vec   <- pfs_vec
  progress_vec <- os_vec - pfs_vec
  dead_vec     <- 1 - os_vec
  
  # monthly accumulators
  monthly_cost <- numeric(horizon)
  monthly_qaly <- numeric(horizon)
  
  for (m in seq_len(horizon)) {
    fs <- stable_vec[m]
    fp <- progress_vec[m]
    
    # if still within 24 months => drug cost
    if (m <= 24) {
      # cost while stable
      c_drug <- fs * drug_cost
      # AE penalty
      u_stable <- base_stable_util - ae_penalty
    } else {
      c_drug <- 0
      u_stable <- base_stable_util
    }
    
    c_stable <- fs * cost_stable
    c_prog   <- fp * cost_progressed
    
    cost_m <- c_drug + c_stable + c_prog
    
    # utility
    util_m <- fs*u_stable + fp*base_prog_util
    
    # discount
    dfact <- (1+0.03)^(-((m-1)/12))
    monthly_cost[m] <- cost_m * dfact
    monthly_qaly[m] <- util_m*(1/12)*dfact
  }
  
  list(cost=sum(monthly_cost), qaly=sum(monthly_qaly))
}

#########################################################
# 3) Subgroup-level calculations
#########################################################

# We'll define a function that calculates cost & QALYs for each subgroup:
calc_subgroup <- function(trulyR=TRUE, predictedR=TRUE) {
  # decide which OS/PFS to use and which cost structure
  if (predictedR) {
    # single agent
    cost_stable <- 2166
    cost_prog   <- 4000 
    drug_cost   <- 15000 # single-agent nivolumab, for example
    if (trulyR) {
      # OS/PFS = single-agent "responder" => same as combo
      os_vec <- os_singleResp
      pfs_vec<- pfs_singleResp
      ae_pen <- 0.0036  # from your table for single agent
    } else {
      # truly non-responder
      os_vec <- os_resistant
      pfs_vec<- pfs_resistant
      ae_pen <- 0.0036
    }
  } else {
    # predicted non-responder => combo
    cost_stable <- 2166
    cost_prog   <- 4000
    drug_cost   <- 26425 # combination cost
    # OS/PFS always combo
    os_vec <- os_comb
    pfs_vec<- pfs_comb
    ae_pen <- 0.012  # from table for combo
  }
  
  # run partition
  partition_ce(os_vec = os_vec, 
               pfs_vec= pfs_vec,
               cost_stable=cost_stable,
               cost_progressed=cost_prog,
               drug_cost=drug_cost,
               ae_penalty=ae_pen,
               base_stable_util=0.754,
               base_prog_util=0.754-0.180,
               horizon=120)
}

#########################################################
# 4) Loop over sensitivity & specificity
#########################################################
sens_values <- seq(0,1,by=0.1)
spec_values <- seq(0,1,by=0.1)

# store results
results <- data.frame(Sens=double(), Spec=double(),
                      Cost=double(), QALY=double())

for (Se in sens_values) {
  for (Sp in spec_values) {
    
    # subpop fractions
    # we define piR=0.6, piNR=0.4
    piR  <- 0.6
    piNR <- 0.4
    
    # 4 subgroups:
    # A) truly R + test says R => fraction= piR*Sp
    # B) truly R + test says NR => fraction= piR*(1-Sp)
    # C) truly NR + test says R => fraction= piNR*(1-Se)
    # D) truly NR + test says NR=> fraction= piNR*Se
    
    # compute cost & QALYs for each subgroup
    # A) trulyR=TRUE, predictedR=TRUE
    resA <- calc_subgroup(trulyR=TRUE, predictedR=TRUE)
    # B) trulyR=TRUE, predictedR=FALSE
    resB <- calc_subgroup(trulyR=TRUE, predictedR=FALSE)
    # C) trulyR=FALSE, predictedR=TRUE
    resC <- calc_subgroup(trulyR=FALSE, predictedR=TRUE)
    # D) trulyR=FALSE, predictedR=FALSE
    resD <- calc_subgroup(trulyR=FALSE, predictedR=FALSE)
    
    # Weighted average
    total_cost <- resA$cost*(piR*Sp) + resB$cost*(piR*(1-Sp)) +
      resC$cost*(piNR*(1-Se)) + resD$cost*(piNR*Se)
    total_qaly <- resA$qaly*(piR*Sp) + resB$qaly*(piR*(1-Sp)) +
      resC$qaly*(piNR*(1-Se)) + resD$qaly*(piNR*Se)
    
    results <- rbind(results,
                     data.frame(Sens=Se, Spec=Sp,
                                Cost=total_cost, QALY=total_qaly))
  }
}

# 'results' now has a grid of (Se,Sp) and the resulting cost & QALY.
# You can then plot or examine the pairs.


# create a heatmap where the axes are sensitivity and specificity and the color represents the cost
library(ggplot2)
library(viridis)
ggplot(results, aes(x=Sens, y=Spec, fill=Cost)) +
  geom_tile() +
  # make the color virdis
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="Cost heatmap", x="Sensitivity", y="Specificity")

# create a heatmap where the axes are sensitivity and specificity and the color represents the QALY
ggplot(results, aes(x=Sens, y=Spec, fill=QALY)) +
  geom_tile() +
  # make the color red
  scale_fill_continuous(low="white", high="red", breaks = 1000) +
  theme_minimal() +
  labs(title="QALY heatmap", x="Sensitivity", y="Specificity") 


base_cost  <- 493596.84    # from your snippet
base_qalys <- 3.121


results$base_cost <- base_cost
results$base_qalys <- base_qalys

# calculate the incremental cost and QALYs
results$incr_cost <- results$Cost - results$base_cost
results$incr_qalys <- results$QALY - results$base_qalys

results$ICER <- results$incr_cost / results$incr_qalys

# create a heatmap where the axes are sensitivity and specificity and the color represents the ICER
ggplot(results, aes(x=Sens, y=Spec, fill=ICER)) +
  geom_tile() +
  scale_fill_continuous(low="white", high="red") +
  theme_minimal() +
  labs(title="ICER heatmap", x="Sensitivity", y="Specificity")



# your "resistant" OS/PFS
os_resistant  <- pmax( pmin( (os_n - 0.6*os_ni)/0.4 , 1 ), 0 )
pfs_resistant <- pmax( pmin( (pfs_n - 0.6*pfs_ni)/0.4 , 1 ), 0 )


library(ggplot2)
library(dplyr)

# Suppose we already know the baseline cost/QALY of Nivolumab alone:
base_cost  <- 493597    # from your snippet
base_qalys <- 3.121


############################################################################
# FULL SCRIPT: PARTITIONED SURVIVAL + SUBGROUP TEST LOGIC
# with lumpsum palliative cost in partition_ce()
############################################################################

library(flexsurv)
library(dplyr)
library(ggplot2)
library(viridis)

# 1) TIME HORIZON & DISCOUNT
max_months <- 120
time_vec   <- seq(0, max_months, by=1)  # monthly cycles
disc_rate  <- 0.03

discount_factor <- function(m) {
  # discrete monthly discount factor for cycle m
  (1 + disc_rate)^(-(m/12))
}

# 2) SURVIVAL CURVES: COMBO & NIVO
#    - best_os_model_ni, best_pfs_model_ni => combo
#    - best_os_model_n,  best_pfs_model_n  => nivo alone
# You must have these flexsurv objects loaded
get_surv <- function(model, times) {
  s_list <- summary(model, t=times, type="survival")
  sapply(s_list, function(x) x$est)
}

# OS/PFS for combination
os_ni  <- get_surv(best_os_model_ni, time_vec)
pfs_ni <- get_surv(best_pfs_model_ni, time_vec)

# OS/PFS for single-agent
os_n  <- get_surv(best_os_model_n,  time_vec)
pfs_n <- get_surv(best_pfs_model_n, time_vec)

# partition for each
stable_ni   <- pfs_ni
progress_ni <- os_ni - pfs_ni
dead_ni     <- 1 - os_ni

stable_n    <- pfs_n
progress_n  <- os_n - pfs_n
dead_n      <- 1 - os_n

############################################################################
# 3) "RESISTANT" SINGLE-AGENT CURVES
#    define your derived curves. Example:
# os_resistant  <- pmax( pmin( (os_n - 0.6*os_ni)/0.4 , 1 ), 0 )
# pfs_resistant <- pmax( pmin( (pfs_n - 0.6*pfs_ni)/0.4 , 1 ), 0 )
############################################################################

#   os_resistant  <- ...
#   pfs_resistant <- ...

############################################################################
# 4) COST AND UTIL PARAMS
############################################################################

# DRUG COSTS (per cycle)
cost_nivolumab <- 14975
cost_combo     <- 26425
max_therapy_months <- 24  # 2 years

# PROG COSTS (includes second line?), per cycle:
cost_prog_combo <- 8908
cost_prog_nivo  <- 12093

# Weighted AE disutilities
ae_penalty_combo <- 0.012
ae_penalty_nivo  <- 0.0036

# Disease cost stable
cost_stable <- 2166

# lumpsum palliative cost at death
cost_palliative_death <- 15957

# utilities
util_stable_base <- 0.754
util_prog_decrement <- 0.180
util_prog <- util_stable_base - util_prog_decrement  # 0.574
util_dead  <- 0.0

############################################################################
# 5) PARTITION FUNCTION (NOW WITH LUMPSUM PALLIATIVE)
############################################################################

partition_ce <- function(os_vec, pfs_vec,
                         cost_stable, cost_progressed,
                         drug_cost, # cost of therapy while stable & <=24 mo
                         ae_penalty, # monthly AE penalty while on therapy
                         base_stable_util, base_prog_util,
                         horizon=120) {
  
  stable_vec   <- pfs_vec
  progress_vec <- os_vec - pfs_vec
  dead_vec     <- 1 - os_vec
  
  monthly_cost <- numeric(horizon)
  monthly_qaly <- numeric(horizon)
  
  for (m in seq_len(horizon)) {
    fs <- stable_vec[m]
    fp <- progress_vec[m]
    
    # if still within 24 months => drug cost
    if (m <= max_therapy_months) {
      c_drug <- fs * drug_cost
      # AE penalty
      u_stable <- base_stable_util - ae_penalty
    } else {
      c_drug <- 0
      u_stable <- base_stable_util
    }
    
    c_stable <- fs * cost_stable
    c_prog   <- fp * cost_progressed
    
    cost_m <- c_drug + c_stable + c_prog
    util_m <- fs*u_stable + fp*base_prog_util
    
    dfact <- discount_factor(m-1)
    monthly_cost[m] <- cost_m * dfact
    monthly_qaly[m] <- util_m*(1/12)*dfact
  }
  
  total_cost <- sum(monthly_cost)
  total_qaly <- sum(monthly_qaly)
  
  # lumpsum palliative cost at death
  n_cycles <- length(stable_vec)
  nd_sum <- 0
  for (m in 2:n_cycles) {
    newly_dead <- (dead_vec[m] - dead_vec[m-1])
    if(newly_dead < 0) newly_dead <- 0
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead*cost_palliative_death*dfact
  }
  total_cost <- total_cost + nd_sum
  
  list(cost=total_cost, qaly=total_qaly)
}

############################################################################
# 6) Reproduce Base 2-arm Analysis
############################################################################

calc_partition_arm <- function(stable_vec, prog_vec, dead_vec, arm="Combo") {
  if (arm=="Combo") {
    partition_ce(stable_vec, prog_vec,
                 cost_stable=2166,
                 cost_progressed=8908,
                 drug_cost=26425,
                 ae_penalty=0.012,
                 base_stable_util=0.754,
                 base_prog_util=0.574,
                 horizon=120)
  } else {
    # single agent
    partition_ce(stable_vec, prog_vec,
                 cost_stable=2166,
                 cost_progressed=12093,
                 drug_cost=14975,
                 ae_penalty=0.0036,
                 base_stable_util=0.754,
                 base_prog_util=0.574,
                 horizon=120)
  }
}

res_combo <- calc_partition_arm(stable_ni, progress_ni, dead_ni, arm="Combo")
res_nivo  <- calc_partition_arm(stable_n, progress_n, dead_n, arm="NivoAlone")

cat("\n=== PARTITIONED SURVIVAL RESULTS (BASE CASE) ===\n")
cat("Nivolumab + Ipilimumab:\n")
cat(sprintf("  Total Cost:  $%.2f\n", res_combo$cost))
cat(sprintf("  Total QALYs: %.3f\n",  res_combo$qaly))

cat("\nNivolumab Alone:\n")
cat(sprintf("  Total Cost:  $%.2f\n", res_nivo$cost))
cat(sprintf("  Total QALYs: %.3f\n",  res_nivo$qaly))

icost <- res_combo$cost - res_nivo$cost
iqaly <- res_combo$qaly - res_nivo$qaly
cat("\nIncremental (Combo vs. Nivo):\n")
cat(sprintf("  Inc Cost: $%.2f\n",  icost))
cat(sprintf("  Inc QALY: %.3f\n",  iqaly))
if(iqaly>0) {
  cat(sprintf("  ICER:     $%.2f per QALY\n", icost/iqaly))
} else {
  cat("  Dominated or extendedly dominated\n")
}

############################################################################
# 7) Subgroup logic for test-based approach
############################################################################

# single-agent "responder" OS/PFS = same as combo
os_singleResp  <- os_ni
pfs_singleResp <- pfs_ni

# single-agent "resistant" OS/PFS => define or load:
# e.g.:
# os_resistant  <- pmax(pmin( (os_n  - 0.6*os_ni)/0.4 ,1),0)
# pfs_resistant <- pmax(pmin( (pfs_n - 0.6*pfs_ni)/0.4,1),0)

calc_subgroup <- function(trulyR=TRUE, predictedR=TRUE) {
  if (predictedR) {
    # single agent
    cost_stable <- 2166
    cost_prog   <- 12093
    drug_cost   <- 14975
    if (trulyR) {
      os_vec <- os_singleResp
      pfs_vec<- pfs_singleResp
      ae_pen <- 0.0036
    } else {
      os_vec <- os_resistant
      pfs_vec<- pfs_resistant
      ae_pen <- 0.0036
    }
  } else {
    # combo
    cost_stable <- 2166
    cost_prog   <- 8908
    drug_cost   <- 26425
    os_vec <- os_ni
    pfs_vec<- pfs_ni
    ae_pen <- 0.012
  }
  
  partition_ce(os_vec, pfs_vec,
               cost_stable   = cost_stable,
               cost_progressed = cost_prog,
               drug_cost     = drug_cost,
               ae_penalty    = ae_pen,
               base_stable_util=0.754,
               base_prog_util=0.574,
               horizon=120)
}

############################################################################
# 8) Loop over Se, Sp
############################################################################

sens_values <- seq(0,1,by=0.1)
spec_values <- seq(0,1,by=0.1)

results <- data.frame(Sens=double(), Spec=double(),
                      Cost=double(), QALY=double())

for (Se in sens_values) {
  for (Sp in spec_values) {
    
    piR  <- 0.60
    piNR <- 0.40
    
    # four sub-subgroups:
    # A) trulyR + test=R => fraction= piR*Sp
    # B) trulyR + test=NR=> fraction= piR*(1-Sp)
    # C) trulyNR+test=R => fraction= piNR*(1-Se)
    # D) trulyNR+test=NR=> fraction= piNR*Se
    
    A <- calc_subgroup(trulyR=TRUE,  predictedR=TRUE)
    B <- calc_subgroup(trulyR=TRUE,  predictedR=FALSE)
    C <- calc_subgroup(trulyR=FALSE, predictedR=TRUE)
    D <- calc_subgroup(trulyR=FALSE, predictedR=FALSE)
    
    total_cost <- A$cost*(piR*Sp) + B$cost*(piR*(1-Sp)) +
      C$cost*(piNR*(1-Se)) + D$cost*(piNR*Se)
    total_qaly <- A$qaly*(piR*Sp) + B$qaly*(piR*(1-Sp)) +
      C$qaly*(piNR*(1-Se)) + D$qaly*(piNR*Se)
    
    results <- rbind(results,
                     data.frame(Sens=Se, Spec=Sp,
                                Cost=total_cost, QALY=total_qaly))
  }
}

# Check corner extremes:
cat("\nCheck extremes:\n")
cat("Se=1,Sp=0 => everyone predicted non-responder => entire pop on combo\n")
cornerA <- subset(results, Sens==1 & Spec==0)
print(cornerA)

cat("\nSe=0,Sp=1 => everyone predicted responder => entire pop on single agent\n")
cornerB <- subset(results, Sens==0 & Spec==1)
print(cornerB)

############################################################################
# 9) Plot Heatmaps
############################################################################
ggplot(results, aes(x=Sens, y=Spec, fill=Cost)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="Cost Heatmap: Biomarker Strategy",
       fill="Total Cost",
       x="Test Sensitivity", y="Test Specificity")

ggplot(results, aes(x=Sens, y=Spec, fill=QALY)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="QALY Heatmap: Biomarker Strategy",
       fill="Total QALYs",
       x="Test Sensitivity", y="Test Specificity")

############################################################################
# DONE
############################################################################








############################################################################
# FINAL INTEGRATED SCRIPT
#  - Reproduce base 2-arm results (Combo vs. Nivo)
#  - Extend to subgroups for predictive test with (Se,Sp) loop
############################################################################

library(flexsurv)
library(dplyr)
library(ggplot2)
library(viridis)

# 1) TIME & DISCOUNT
max_months <- 120
time_vec   <- seq(0, max_months, by=1)
disc_rate  <- 0.03

discount_factor <- function(m) {
  (1 + disc_rate)^(-(m/12))
}

# 2) LOAD SURVIVAL CURVES FROM flexsurv
#    best_os_model_ni, best_pfs_model_ni => combo
#    best_os_model_n,  best_pfs_model_n  => nivo alone
#    define them or load them before running

get_surv <- function(model, times) {
  s_list <- summary(model, t=times, type="survival")
  sapply(s_list, function(x) x$est)
}

# OS/PFS for combo
os_ni  <- get_surv(best_os_model_ni,  time_vec)
pfs_ni <- get_surv(best_pfs_model_ni, time_vec)

# OS/PFS for single
os_n   <- get_surv(best_os_model_n,   time_vec)
pfs_n  <- get_surv(best_pfs_model_n,  time_vec)

stable_ni   <- pfs_ni
progress_ni <- os_ni - pfs_ni
dead_ni     <- 1 - os_ni

stable_n    <- pfs_n
progress_n  <- os_n - pfs_n
dead_n      <- 1 - os_n

# 3) "RESISTANT" single-agent
# define or load them from your formula:
# os_resistant  <- ...
# pfs_resistant <- ...

############################################################################
# 4) PARAMETERS (matching final table)
############################################################################
cost_nivolumab <- 14975     # single-agent monthly cost
cost_combo     <- 26425
max_therapy_months <- 24

cost_prog_combo <- 8908
cost_prog_nivo  <- 12093

cost_stable <- 2166
cost_palliative_death <- 15957

ae_penalty_combo <- 0.012
ae_penalty_nivo  <- 0.0036

util_stable_base   <- 0.754
util_prog_decrement<- 0.180
util_prog          <- util_stable_base - util_prog_decrement # => 0.574
util_dead          <- 0.0

############################################################################
# 5) PARTITION FUNCTION (unified lumpsum approach)
############################################################################
partition_ce <- function(stable_vec, prog_vec, dead_vec,
                         drug_cost, cost_stable, cost_prog,
                         ae_penalty,
                         horizon=120) {
  
  # stable=fs, progressed=fp, dead=?
  monthly_cost <- numeric(horizon)
  monthly_qaly <- numeric(horizon)
  
  for (m in seq_len(horizon)) {
    fs <- stable_vec[m]
    fp <- prog_vec[m]
    
    # if stable & within 24 mo => pay drug cost
    if (m <= max_therapy_months) {
      c_drug   <- fs*drug_cost
      u_stable <- util_stable_base - ae_penalty
    } else {
      c_drug   <- 0
      u_stable <- util_stable_base
    }
    
    c_stab <- fs*cost_stable
    c_prog <- fp*cost_prog
    
    cost_m <- c_drug + c_stab + c_prog
    util_m <- fs*u_stable + fp*util_prog
    
    # discount
    dfact <- discount_factor(m-1)
    monthly_cost[m] <- cost_m*dfact
    monthly_qaly[m] <- util_m*(1/12)*dfact
  }
  
  total_cost <- sum(monthly_cost)
  total_qaly <- sum(monthly_qaly)
  
  # lumpsum palliative cost for newly dead
  nd_sum <- 0
  n_cycles <- length(stable_vec)
  for (m in 2:n_cycles) {
    newly_dead <- (dead_vec[m] - dead_vec[m-1])
    if (newly_dead<0) newly_dead <- 0
    dfact <- discount_factor(m-1)
    nd_sum <- nd_sum + newly_dead*cost_palliative_death*dfact
  }
  
  total_cost <- total_cost + nd_sum
  
  list(cost=total_cost, qaly=total_qaly)
}

############################################################################
# 6) REPLICATE 2-ARM BASE-CASE
############################################################################

calc_partition_arm <- function(stable_vec, prog_vec, dead_vec, arm="Combo") {
  # Just picks the right param for single vs. combo
  if (arm=="Combo") {
    partition_ce(stable_vec, prog_vec, dead_vec,
                 drug_cost = cost_combo,
                 cost_stable = cost_stable,
                 cost_prog   = cost_prog_combo,
                 ae_penalty  = ae_penalty_combo,
                 horizon=120)
  } else {
    # Nivo alone
    partition_ce(stable_vec, prog_vec, dead_vec,
                 drug_cost   = cost_nivolumab,
                 cost_stable = cost_stable,
                 cost_prog   = cost_prog_nivo,
                 ae_penalty  = ae_penalty_nivo,
                 horizon=120)
  }
}

res_combo <- calc_partition_arm(stable_ni, progress_ni, dead_ni, arm="Combo")
res_nivo  <- calc_partition_arm(stable_n,  progress_n,  dead_n,  arm="NivoAlone")

cat("\n=== BASE 2-ARM RESULTS ===\n")
cat("Nivo+IPI:\n")
cat(sprintf("  Cost:  $%.2f\n", res_combo$cost))
cat(sprintf("  QALY:  %.3f\n",  res_combo$qaly))

cat("\nNivo alone:\n")
cat(sprintf("  Cost:  $%.2f\n", res_nivo$cost))
cat(sprintf("  QALY:  %.3f\n",  res_nivo$qaly))

icost <- res_combo$cost - res_nivo$cost
iqaly <- res_combo$qaly - res_nivo$qaly
cat("\nIncremental (Combo vs. Nivo):\n")
cat(sprintf("  IncCost: $%.2f\n", icost))
cat(sprintf("  IncQaly: %.3f\n", iqaly))
if (iqaly>0) {
  cat(sprintf("  ICER: $%.2f per QALY\n", icost/iqaly))
} else {
  cat("  Dominated or extended domin.\n")
}

############################################################################
# 7) SUBGROUP APPROACH
############################################################################
# single-agent "responder" => same OS/PFS as combo
os_singleResp  <- os_ni
pfs_singleResp <- pfs_ni

# single-agent "resistant"
#  define them or load them
# os_resistant  <- ...
# pfs_resistant <- ...

calc_subgroup <- function(trulyR=TRUE, predictedR=TRUE) {
  if (predictedR) {
    # single-agent
    drug_c <- cost_nivolumab
    c_stab <- cost_stable
    c_prog <- cost_prog_nivo
    aep    <- ae_penalty_nivo
    
    if (trulyR) {
      # "responder" => same OS as combo
      os_vec <- os_singleResp
      pfs_vec<- pfs_singleResp
    } else {
      # "resistant"
      os_vec <- os_resistant
      pfs_vec<- pfs_resistant
    }
  } else {
    # predicted non-responder => combo
    drug_c <- cost_combo
    c_stab <- cost_stable
    c_prog <- cost_prog_combo
    aep    <- ae_penalty_combo
    
    os_vec <- os_ni
    pfs_vec<- pfs_ni
  }
  
  partition_ce(stable_vec=pfs_vec,
               prog_vec=(os_vec - pfs_vec),
               dead_vec=(1 - os_vec),
               drug_cost   = drug_c,
               cost_stable = c_stab,
               cost_prog   = c_prog,
               ae_penalty  = aep,
               horizon=120)
}

############################################################################
# 8) Loop over (Se,Sp)
############################################################################
sens_values <- seq(0,1,by=0.1)
spec_values <- seq(0,1,by=0.1)

results <- data.frame(Sens=numeric(), Spec=numeric(),
                      Cost=numeric(), QALY=numeric())

for (Se in sens_values) {
  for (Sp in spec_values) {
    # subpop fractions
    piR  <- 0.60
    piNR <- 0.40
    
    # four sub-subgroups
    # A) truly R + test=R => fraction= piR*Sp
    # B) truly R + test=NR=> fraction= piR*(1-Sp)
    # C) truly NR+test=R => fraction= piNR*(1-Se)
    # D) truly NR+test=NR=> fraction= piNR*Se
    A <- calc_subgroup(trulyR=TRUE,  predictedR=TRUE)
    B <- calc_subgroup(trulyR=TRUE,  predictedR=FALSE)
    C <- calc_subgroup(trulyR=FALSE, predictedR=TRUE)
    D <- calc_subgroup(trulyR=FALSE, predictedR=FALSE)
    
    total_cost <- A$cost*(piR*Sp) + B$cost*(piR*(1-Sp)) +
      C$cost*(piNR*(1-Se)) + D$cost*(piNR*Se)
    total_qaly <- A$qaly*(piR*Sp) + B$qaly*(piR*(1-Sp)) +
      C$qaly*(piNR*(1-Se)) + D$qaly*(piNR*Se)
    
    results <- rbind(results,
                     data.frame(Sens=Se, Spec=Sp,
                                Cost=total_cost,
                                QALY=total_qaly))
  }
}

cat("\nCheck extremes:\n")
cat("Se=1,Sp=0 => entire pop predicted 'NR' => all combo\n")
cornerA <- subset(results, Sens==1 & Spec==0)
print(cornerA)

cat("\nSe=0,Sp=1 => entire pop predicted 'R' => all single agent\n")
cornerB <- subset(results, Sens==0 & Spec==1)
print(cornerB)

############################################################################
# 9) Plot the heatmaps
############################################################################
ggplot(results, aes(x=Sens, y=Spec, fill=Cost)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="Cost Heatmap: Biomarker-based strategy",
       fill="Total Cost ($)",
       x="Test Sensitivity", y="Test Specificity")

ggplot(results, aes(x=Sens, y=Spec, fill=QALY)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="QALY Heatmap: Biomarker-based strategy",
       fill="QALYs",
       x="Test Sensitivity", y="Test Specificity")

############################################################################
# DONE
############################################################################



base_cost <-  492021.16
base_qaly <-  3.105

results$base_cost <- base_cost
results$base_qaly <- base_qaly

results$icost <- results$Cost - base_cost
results$iqaly <- results$QALY - base_qaly

results$icer <- results$icost / results$iqaly

# heat map of postive icers

ggplot(subset(results, icer > 0), aes(x=Sens, y=Spec, fill=icer)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title="ICER Heatmap: Biomarker-based strategy",
       fill="ICER ($/QALY)",
       x="Test Sensitivity", y="Test Specificity") +
  # make the color scale values more readable
  scale_fill_viridis(labels=scales::dollar_format(scale=1e-1))

# show 
