# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define a common time vector (in months)
time_vec <- seq(0, 120, by = 1)  # From 0 to 120 months (10 years)

0.3759776
0.2432601
# -----------------------------------------------
# Step 1: Get Predictions from the Survival Models
# -----------------------------------------------
# For combination therapy (nivo+ipi)
os_comb_summary <- summary(best_os_model_ni, t = time_vec, type = "survival")
pfs_comb_summary <- summary(best_pfs_model_ni, t = time_vec, type = "survival")

# For monotherapy (nivo only)
os_mon_summary <- summary(best_os_model_n, t = time_vec, type = "survival")
pfs_mon_summary <- summary(best_pfs_model_n, t = time_vec, type = "survival")

# Extract survival estimates as vectors using the same time grid.
os_comb <- sapply(os_comb_summary, function(x) x$est)
pfs_comb <- sapply(pfs_comb_summary, function(x) x$est)
os_mon  <- sapply(os_mon_summary, function(x) x$est)
pfs_mon <- sapply(pfs_mon_summary, function(x) x$est)

# -----------------------------------------------
# Step 2: Compute Counterfactual Survival for Intrinsically Resistant Patients
# -----------------------------------------------
# We assume that the monotherapy survival is a weighted average:
# OS_mon = 0.8 * OS_comb + 0.2 * OS_resistant
# Solving for OS_resistant:
os_resistant <- (os_mon - (0.6 * os_comb)) / 0.4

# Similarly, for PFS:
pfs_resistant <- (pfs_mon - (0.6 * pfs_comb)) / 0.4
pfs_resistant <- pmax(pmin(pfs_resistant, 1), 0)

# -----------------------------------------------
# Step 3: Create Data Frames for Plotting
# -----------------------------------------------
df_OS <- data.frame(
  time = time_vec,
  Nivolumab = os_mon,
  Nivolumab.Ipilimumab = os_comb,
  Resistant   = os_resistant
)

df_PFS <- data.frame(
  time = time_vec,
  Nivolumab = pfs_mon,
  Nivolumab.Ipilimumab = pfs_comb,
  Resistant   = pfs_resistant
)

# Reshape data to long format for ggplot
df_OS_long <- pivot_longer(df_OS, cols = c("Nivolumab", "Nivolumab.Ipilimumab", "Resistant"),
                           names_to = "Group", values_to = "surv")

df_PFS_long <- pivot_longer(df_PFS, cols = c("Nivolumab", "Nivolumab.Ipilimumab", "Resistant"),
                            names_to = "Group", values_to = "surv")

df_OS_long$type <- "OS"
df_PFS_long$type <- "PFS"
# -----------------------------------------------
# Step 4: Plot the Survival Curves
# -----------------------------------------------

df <- rbind(df_OS_long, df_PFS_long)

# replace values in the Group column
df %>%
  mutate(Group = recode(Group, 
                        "Nivolumab" = "Anti-PD-1 (Nivolumab) Model",
                        "Nivolumab.Ipilimumab" = "Anti-PD-1 Responders",
                        "Resistant" = "Anti-PD-1 Resistant")) -> df


# Overall Survival Plot
ggplot(df, aes(x = time, y = surv, color = Group)) +
  geom_line(size = 1.15) +
  labs(title = "Counterfactual OS and PFS for Anti-PD-1 Responders/Non-Responders",
       x = "Time (months)",
       y = "Survival Probability",
       color = '') +
  facet_wrap(~type) +
  scale_linetype_manual(values = c("OS" = "solid", "PFS" = "dashed")) +
  scale_color_manual(values = c("Anti-PD-1 (Nivolumab) Model" = "#008Fd5", 
                                "Anti-PD-1 Responders" = "#77AB43", 
                                "Anti-PD-1 Resistant" = "#FF2700")) +
  custom_theme +
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     limits = c(0,1)) +
  theme(legend.position = "top") 

ggsave("counterfactual_survival.png", width = 10, height = 6, dpi = 300) 




