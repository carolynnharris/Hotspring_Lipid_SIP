cat("\014")
rm(list = ls())
graphics.off()

################################################################################
# 06_Generate_IPL_BP_summary_table.R
#
# This script:
# 1. Loads the cleaned merged LH-SIP dataset
# 2. Filters to IPL-derived biphytane data only
# 3. Summarizes biological replicates by Site x Type x TimePoint
# 4. Calculates mean ± SEM values for:
#      - BP relative abundance (%)
#      - BP Ring Index
#      - compound-specific BP d2H values
#      - abundance-weighted mean BP d2H
#      - mean ∆d2H/ring
# 5. Formats the output to generate Supplementary Table S1.
#
# Output:
# - summary_tables/TableS1_IPL_BP_summary_raw.csv
# - summary_tables/TableS1_IPL_BP_summary_formatted.csv
################################################################################


#### packages ####
library(tidyverse)

#### paths ####
clean_dir <- "data_inputs"
summary_dir <- "summary_tables"


#### load data ####
dat <- read.csv(file.path(clean_dir, "LHSIP_All_Merge.csv")) %>%
  select(-any_of(c("X", "...1", "Unnamed: 0"))) %>%
  filter(Fraction == "IPL") %>%
  mutate(
    TimePoint = as.numeric(TimePoint)
  )

#### helper functions ####
sem <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

# Use replicate SEM when there is more than one replicate
# If there is only one rep, use the already-propagated analytical SEM column
group_sem <- function(value_col, analytical_sem_col) {
  n_vals <- sum(!is.na(value_col))
  
  if (n_vals > 1) {
    sem(value_col)
  } else if (n_vals == 1) {
    analytical_sem_col[which(!is.na(value_col))[1]]
  } else {
    NA_real_
  }
}

format_mean_sem <- function(mean, sem, digits_mean = 0, digits_sem = 1) {
  case_when(
    is.na(mean) ~ "-",
    is.na(sem) ~ as.character(round(mean, digits_mean)),
    TRUE ~ paste0(
      round(mean, digits_mean),
      " ± ",
      round(sem, digits_sem)
    )
  )
}

#### summarize IPL BP data ####
BP_table_raw <- dat %>%
  group_by(Site, Type, Incubation_days) %>%
  summarise(
    N = n(),
    
    BP0_RelAbund_mean = mean(BP0_RelAbund, na.rm = TRUE),
    BP0_RelAbund_sem  = sem(BP0_RelAbund),
    BP1_RelAbund_mean = mean(BP1_RelAbund, na.rm = TRUE),
    BP1_RelAbund_sem  = sem(BP1_RelAbund),
    BP2_RelAbund_mean = mean(BP2_RelAbund, na.rm = TRUE),
    BP2_RelAbund_sem  = sem(BP2_RelAbund),
    BP3_RelAbund_mean = mean(BP3_RelAbund, na.rm = TRUE),
    BP3_RelAbund_sem  = sem(BP3_RelAbund),
    BP4_RelAbund_mean = mean(BP4_RelAbund, na.rm = TRUE),
    BP4_RelAbund_sem  = sem(BP4_RelAbund),
    
    RI_BP_mean = mean(RI_BP, na.rm = TRUE),
    RI_BP_sem  = sem(RI_BP),
    
    BP0_d2H_mean_mean = mean(BP0_d2H_mean, na.rm = TRUE),
    BP0_d2H_mean_sem  = group_sem(BP0_d2H_mean, BP0_d2H_sem),
    
    BP1_d2H_mean_mean = mean(BP1_d2H_mean, na.rm = TRUE),
    BP1_d2H_mean_sem  = group_sem(BP1_d2H_mean, BP1_d2H_sem),
    
    BP2_d2H_mean_mean = mean(BP2_d2H_mean, na.rm = TRUE),
    BP2_d2H_mean_sem  = group_sem(BP2_d2H_mean, BP2_d2H_sem),
    
    d2H_wt_mean_mean = mean(d2H_wt_mean, na.rm = TRUE),
    d2H_wt_mean_sem  = group_sem(d2H_wt_mean, d2H_wt_mean_sem),
    
    RingDiff_all_mean = mean(RingDiff_all, na.rm = TRUE),
    RingDiff_all_sem  = group_sem(RingDiff_all, RingDiff_all_sem),
    
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

#### order rows ####
BP_table_raw <- BP_table_raw %>%
  mutate(
    Site = factor(Site, levels = c("Beryl", "ETAT-3")),
    Type = factor(Type, levels = c("Control-Kill", "Control-Live", "Experimental"))
  ) %>%
  arrange(Site, Type, Incubation_days)

#### format table  ####
BP_table_formatted <- BP_table_raw %>%
  transmute(
    Site = as.character(Site),
    Type = as.character(Type),
    `Time Point (days)` = Incubation_days,
    N,
    
    `BP-0 Rel. Abund. (%)` = format_mean_sem(BP0_RelAbund_mean, BP0_RelAbund_sem, 0, 1),
    `BP-1 Rel. Abund. (%)` = format_mean_sem(BP1_RelAbund_mean, BP1_RelAbund_sem, 0, 1),
    `BP-2 Rel. Abund. (%)` = format_mean_sem(BP2_RelAbund_mean, BP2_RelAbund_sem, 0, 1),
    `BP-3 Rel. Abund. (%)` = format_mean_sem(BP3_RelAbund_mean, BP3_RelAbund_sem, 0, 1),
    `BP-4 Rel. Abund. (%)` = format_mean_sem(BP4_RelAbund_mean, BP4_RelAbund_sem, 0, 1),
    
    `Ring Index` = format_mean_sem(RI_BP_mean, RI_BP_sem, 2, 2),
    
    `BP-0 δ2H (‰)` = format_mean_sem(BP0_d2H_mean_mean, BP0_d2H_mean_sem, 0, 1),
    `BP-1 δ2H (‰)` = format_mean_sem(BP1_d2H_mean_mean, BP1_d2H_mean_sem, 0, 1),
    `BP-2 δ2H (‰)` = format_mean_sem(BP2_d2H_mean_mean, BP2_d2H_mean_sem, 0, 1),
    `Wt. Mean δ2H (‰)` = format_mean_sem(d2H_wt_mean_mean, d2H_wt_mean_sem, 0, 1),
    
    `∆δ2H/ring (‰)` = format_mean_sem(RingDiff_all_mean, RingDiff_all_sem, 0, 1)
  )

#### save outputs ####
write.csv(
  BP_table_raw,
  file.path(summary_dir, "TableS1_IPL_BP_summary_raw.csv"),
  row.names = FALSE
)

write.csv(
  BP_table_formatted,
  file.path(summary_dir, "TableS1_IPL_BP_summary_formatted.csv"),
  row.names = FALSE
)

BP_table_formatted
