cat("\014")
rm(list=ls()) 
graphics.off()

################################################################################
# 02_generate_summary_tables.R
#
# This script:
# 1. Loads the cleaned merged LH-SIP dataset
# 2. Subsets to IPL-derived biphytane data only
# 3. Summarizes biological replicates by Site x Type x incubation_days
# 4. Calculates Δ2H enrichment in ppm relative to the day-0 value.
# 5. Flags values above the 2 ppm assay detection threshold
# 6. Calculates growth rate and apparent generation time using:
#       alpha = 0.66 ± 0.10
# 7. Propagates uncertainty in F0, Ft, tracer water, and alpha
# 8. Calculates integrated 3+14 day growth estimates for experimental samples
# 9. Saves summary tables for downstream tables and figures
# 10. Generates Table 2
# Outputs:
# - summary_tables/Table2_growth_summary_by_timepoint.csv
################################################################################

#### load packages ####
library("tidyverse")

#### load data ####
dat <- read.csv("data_inputs/LHSIP_All_Merge.csv") 

#### define assimilation efficiency & detection limit ####
alpha_mean <- 0.66
alpha_sd   <- 0.10

detection_limit_ppm <- 2.0011 # calculated elsewhere

#### subset for IPL lipids only ####
dat <- dat %>%
  filter(Fraction == "IPL") %>%
  mutate(
    Incubation_days = round(Incubation_days)
  )


#### helper functions ####

sem_replicates <- function(x) {
  n_vals <- sum(!is.na(x))
  
  if (n_vals <= 1) {
    return(NA_real_)
  } else {
    return(sd(x, na.rm = TRUE) / sqrt(n_vals))
  }
}

combine_mean_sem <- function(values, sems) {
  n_vals <- sum(!is.na(values))
  
  if (n_vals == 0) {
    return(NA_real_)
  }
  
  if (n_vals == 1) {
    return(sems[!is.na(values)][1])
  }
  
  sqrt(sum(sems[!is.na(values)]^2, na.rm = TRUE) / n_vals^2)
}

calc_growth <- function(Ft, F0, FL, t_days, Ft_sem, F0_sem, FL_sem,
                        alpha = alpha_mean, alpha_err = alpha_sd) {
  
  A <- alpha * FL - Ft
  B <- alpha * FL - F0
  
  if (is.na(A) | is.na(B) | is.na(t_days) | t_days <= 0) {
    return(tibble(
      mu_day = NA_real_, mu_day_sem = NA_real_,
      mu_year = NA_real_, mu_year_sem = NA_real_,
      T_G_years = NA_real_, T_G_years_sem = NA_real_
    ))
  }
  
  if (A <= 0 | B <= 0) {
    return(tibble(
      mu_day = NA_real_, mu_day_sem = NA_real_,
      mu_year = NA_real_, mu_year_sem = NA_real_,
      T_G_years = NA_real_, T_G_years_sem = NA_real_
    ))
  }
  
  mu_day <- -(1 / t_days) * log(A / B)
  
  if (mu_day <= 0) {
    return(tibble(
      mu_day = NA_real_, mu_day_sem = NA_real_,
      mu_year = NA_real_, mu_year_sem = NA_real_,
      T_G_years = NA_real_, T_G_years_sem = NA_real_
    ))
  }
  
  # Partial derivatives for propagated uncertainty
  dmu_dFt    <-  (1 / t_days) * (1 / A)
  dmu_dF0    <- -(1 / t_days) * (1 / B)
  dmu_dFL    <- -(alpha / t_days) * ((1 / A) - (1 / B))
  dmu_dalpha <- -(FL / t_days) * ((1 / A) - (1 / B))
  
  mu_day_sem <- sqrt(
    (dmu_dFt    * Ft_sem)^2 +
      (dmu_dF0    * F0_sem)^2 +
      (dmu_dFL    * FL_sem)^2 +
      (dmu_dalpha * alpha_err)^2
  )
  
  mu_year <- mu_day * 365
  mu_year_sem <- mu_day_sem * 365
  
  T_G_years <- log(2) / mu_year
  T_G_years_sem <- abs(log(2) / mu_year^2) * mu_year_sem
  
  tibble(
    mu_day = mu_day,
    mu_day_sem = mu_day_sem,
    mu_year = mu_year,
    mu_year_sem = mu_year_sem,
    T_G_years = T_G_years,
    T_G_years_sem = T_G_years_sem
  )
}




#### summary table ####

vars_with_sems <- c(
  BP0_d2H_mean = "BP0_d2H_sem",
  BP1_d2H_mean = "BP1_d2H_sem",
  BP2_d2H_mean = "BP2_d2H_sem",
  BP3_d2H_mean = "BP3_d2H_sem",
  BP4_d2H_mean = "BP4_d2H_sem",
  d2H_wt_mean = "d2H_wt_mean_sem",
  RingDiff_all = "RingDiff_all_sem",
  BP0_D_ppm_mean = "BP0_D_ppm_mean_sem",
  BP1_D_ppm_mean = "BP1_D_ppm_mean_sem",
  BP2_D_ppm_mean = "BP2_D_ppm_mean_sem",
  BP3_D_ppm_mean = "BP3_D_ppm_mean_sem",
  BP4_D_ppm_mean = "BP4_D_ppm_mean_sem",
  D_ppm_wt_mean = "D_ppm_wt_mean_sem",
  Tracer_Water_D_ppm = "Tracer_Water_D_ppm_sem"
)

vars_without_sems <- c(
  "BP0_RelAbund",
  "BP1_RelAbund",
  "BP2_RelAbund",
  "BP3_RelAbund",
  "BP4_RelAbund",
  "RI_BP",
  "BP_ug_sediment_g"
)

dat_summary <- dat %>%
  group_by(Site, Type, Incubation_days) %>%
  summarise(
    n = n(),
    
    across(
      all_of(names(vars_with_sems)),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sem = ~ combine_mean_sem(
          values = .x,
          sems = pick(all_of(vars_with_sems[[cur_column()]]))[[1]]
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    
    across(
      all_of(vars_without_sems),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sem = ~ sem_replicates(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

#### calculate Δ2H relative to day 0 ####

vars_to_delta <- c(
  "D_ppm_wt_mean",
  "BP0_D_ppm_mean",
  "BP1_D_ppm_mean",
  "BP2_D_ppm_mean"
)

for (var_base in vars_to_delta) {
  
  mean_col <- paste0(var_base, "_mean")
  sem_col  <- paste0(var_base, "_sem")
  
  delta_col <- paste0("Delta_", var_base)
  delta_sem_col <- paste0(delta_col, "_sem")
  
  day0_vals <- dat_summary %>%
    filter(Incubation_days == 0) %>%
    select(
      Site,
      Type,
      day0_mean = all_of(mean_col),
      day0_sem = all_of(sem_col)
    )
  
  dat_summary <- dat_summary %>%
    left_join(day0_vals, by = c("Site", "Type")) %>%
    mutate(
      !!delta_col := case_when(
        is.na(.data[[mean_col]]) | is.na(day0_mean) ~ NA_real_,
        Incubation_days == 0 ~ 0,
        TRUE ~ .data[[mean_col]] - day0_mean
      ),
      
      !!delta_sem_col := case_when(
        is.na(.data[[sem_col]]) | is.na(day0_sem) ~ NA_real_,
        Incubation_days == 0 ~ 0,
        TRUE ~ sqrt(.data[[sem_col]]^2 + day0_sem^2)
      )
    ) %>%
    select(-day0_mean, -day0_sem)
}


#### flag 2H uptake above detection limit ####

dat_summary <- dat_summary %>%
  mutate(
    above_detection_limit_wt_mean =
      Delta_D_ppm_wt_mean >= detection_limit_ppm
  )

#### calculate growth rates and generation times ####

# Add day-0 baseline values
growth_input <- dat_summary %>%
  filter(Type == "Experimental") %>%
  select(
    Site,
    Type,
    Incubation_days,
    D_ppm_wt_mean_mean,
    D_ppm_wt_mean_sem,
    Tracer_Water_D_ppm_mean,
    Tracer_Water_D_ppm_sem,
    starts_with("BP0_D_ppm_mean"),
    starts_with("BP1_D_ppm_mean"),
    starts_with("BP2_D_ppm_mean"),
    starts_with("Delta_")
  )

baseline <- growth_input %>%
  filter(Incubation_days == 0) %>%
  select(
    Site,
    F0_wt = D_ppm_wt_mean_mean,
    F0_wt_sem = D_ppm_wt_mean_sem,
    F0_BP0 = BP0_D_ppm_mean_mean,
    F0_BP0_sem = BP0_D_ppm_mean_sem,
    F0_BP1 = BP1_D_ppm_mean_mean,
    F0_BP1_sem = BP1_D_ppm_mean_sem,
    F0_BP2 = BP2_D_ppm_mean_mean,
    F0_BP2_sem = BP2_D_ppm_mean_sem
  )

# summarize growth results
growth_results <- growth_input %>%
  filter(Incubation_days > 0) %>%
  left_join(baseline, by = "Site") %>%
  rowwise() %>%
  mutate(
    growth_wt = list(calc_growth(
      Ft = D_ppm_wt_mean_mean,
      F0 = F0_wt,
      FL = Tracer_Water_D_ppm_mean,
      t_days = Incubation_days,
      Ft_sem = D_ppm_wt_mean_sem,
      F0_sem = F0_wt_sem,
      FL_sem = Tracer_Water_D_ppm_sem
    )),
    
    growth_BP0 = list(calc_growth(
      Ft = BP0_D_ppm_mean_mean,
      F0 = F0_BP0,
      FL = Tracer_Water_D_ppm_mean,
      t_days = Incubation_days,
      Ft_sem = BP0_D_ppm_mean_sem,
      F0_sem = F0_BP0_sem,
      FL_sem = Tracer_Water_D_ppm_sem
    )),
    
    growth_BP1 = list(calc_growth(
      Ft = BP1_D_ppm_mean_mean,
      F0 = F0_BP1,
      FL = Tracer_Water_D_ppm_mean,
      t_days = Incubation_days,
      Ft_sem = BP1_D_ppm_mean_sem,
      F0_sem = F0_BP1_sem,
      FL_sem = Tracer_Water_D_ppm_sem
    )),
    
    growth_BP2 = list(calc_growth(
      Ft = BP2_D_ppm_mean_mean,
      F0 = F0_BP2,
      FL = Tracer_Water_D_ppm_mean,
      t_days = Incubation_days,
      Ft_sem = BP2_D_ppm_mean_sem,
      F0_sem = F0_BP2_sem,
      FL_sem = Tracer_Water_D_ppm_sem
    ))
  ) %>%
  ungroup() %>%
  unnest_wider(growth_wt, names_sep = "_") %>%
  unnest_wider(growth_BP0, names_sep = "_") %>%
  unnest_wider(growth_BP1, names_sep = "_") %>%
  unnest_wider(growth_BP2, names_sep = "_")


#### integrate growth estimates over 3 + 14 day time points ####
# Growth rates are weighted by incubation duration
# ETAT-3 is excluded here because no 2H uptake was detected
# ETAT-3 minimum generation time is calculated separately from the assay detection limit

integrated_growth <- growth_results %>%
  filter(Site == "Beryl") %>% # Exclude ETAT-3 because it showed no uptake
  group_by(Site) %>%
  summarise(
    
    weights = list(Incubation_days),
    
    Type = "Experimental",
    Incubation_days = NA_real_,
    Time_Point = "Integrated_3_14_days",
    
    Delta_D_ppm_wt_mean_mean =
      mean(Delta_D_ppm_wt_mean, na.rm = TRUE),
    
    Delta_D_ppm_wt_mean_sem =
      sqrt(sum(Delta_D_ppm_wt_mean_sem^2, na.rm = TRUE)) / n(),
    
    
    #### weighted-average growth rates ####
    
    mu_year_wt =
      weighted.mean(
        growth_wt_mu_year,
        w = unlist(weights),
        na.rm = TRUE
      ),
    
    mu_year_wt_sem =
      sqrt(sum((unlist(weights) * growth_wt_mu_year_sem)^2, na.rm = TRUE)) /
      sum(unlist(weights), na.rm = TRUE),
    
    T_G_years_wt = log(2) / mu_year_wt,
    
    T_G_years_wt_sem =
      abs(log(2) / mu_year_wt^2) * mu_year_wt_sem,
    
    
    #### BP0 ####
    
    mu_year_BP0 =
      weighted.mean(
        growth_BP0_mu_year,
        w = unlist(weights),
        na.rm = TRUE
      ),
    
    mu_year_BP0_sem =
      sqrt(sum((unlist(weights) * growth_BP0_mu_year_sem)^2, na.rm = TRUE)) /
      sum(unlist(weights), na.rm = TRUE),
    
    T_G_years_BP0 = log(2) / mu_year_BP0,
    
    T_G_years_BP0_sem =
      abs(log(2) / mu_year_BP0^2) * mu_year_BP0_sem,
    
    
    #### BP1 ####
    
    mu_year_BP1 =
      weighted.mean(
        growth_BP1_mu_year,
        w = unlist(weights),
        na.rm = TRUE
      ),
    
    mu_year_BP1_sem =
      sqrt(sum((unlist(weights) * growth_BP1_mu_year_sem)^2, na.rm = TRUE)) /
      sum(unlist(weights), na.rm = TRUE),
    
    T_G_years_BP1 = log(2) / mu_year_BP1,
    
    T_G_years_BP1_sem =
      abs(log(2) / mu_year_BP1^2) * mu_year_BP1_sem,
    
    
    #### BP2 ####
    
    mu_year_BP2 =
      weighted.mean(
        growth_BP2_mu_year,
        w = unlist(weights),
        na.rm = TRUE
      ),
    
    mu_year_BP2_sem =
      sqrt(sum((unlist(weights) * growth_BP2_mu_year_sem)^2, na.rm = TRUE)) /
      sum(unlist(weights), na.rm = TRUE),
    
    T_G_years_BP2 = log(2) / mu_year_BP2,
    
    T_G_years_BP2_sem =
      abs(log(2) / mu_year_BP2^2) * mu_year_BP2_sem,
    
    .groups = "drop"
    
  ) %>%
  select(-weights)


#### join growth results back to main summary table ####

dat_summary_growth <- dat_summary %>%
  left_join(
    growth_results %>%
      select(
        Site,
        Type,
        Incubation_days,
        starts_with("growth_")
      ),
    by = c("Site", "Type", "Incubation_days")
  )


#### save tables ####

dir.create("summary_tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  dat_summary_growth,
  file = "summary_tables/LHSIP_Summary_Table.csv",
  row.names = FALSE
)

write.csv(
  growth_results,
  file = "summary_tables/LHSIP_Growth_Rates_by_Timepoint.csv",
  row.names = FALSE
)

write.csv(
  integrated_growth,
  file = "summary_tables/LHSIP_Growth_Rates_Integrated_3_14_days.csv",
  row.names = FALSE
)



#### Table 2 - summary of 2H-uptake, growth rate, TG - Beryl data ####

table2_beryl_timepoints <- dat_summary_growth %>%
  filter(
    Site == "Beryl",
    Type == "Experimental",
    Incubation_days %in% c(3, 14)
  ) %>%
  transmute(
    Site,
    `Time Point` = case_when(
      Incubation_days == 3  ~ "3 days*",
      Incubation_days == 14 ~ "14 days"
    ),
    `Lipid Δ²H (ppm)` = paste0(
      round(Delta_D_ppm_wt_mean, 2),
      " ± ",
      round(Delta_D_ppm_wt_mean_sem, 2)
    ),
    `μ (year⁻¹)` = paste0(
      round(growth_wt_mu_year, 3),
      " ± ",
      round(growth_wt_mu_year_sem, 2)
    ),
    `TG (years)` = paste0(
      round(growth_wt_T_G_years),
      " ± ",
      round(growth_wt_T_G_years_sem)
    )
  )


table2_beryl_average <- integrated_growth %>%
  transmute(
    Site,
    `Time Point` = "Average",
    `Lipid Δ²H (ppm)` = paste0(
      round(Delta_D_ppm_wt_mean_mean, 2),
      " ± ",
      round(Delta_D_ppm_wt_mean_sem, 2)
    ),
    `μ (year⁻¹)` = paste0(
      round(mu_year_wt, 3),
      " ± ",
      round(mu_year_wt_sem, 2)
    ),
    `TG (years)` = paste0(
      round(T_G_years_wt),
      " ± ",
      round(T_G_years_wt_sem)
    )
  )

table2_etat3 <- dat_summary_growth %>% # calculated elsewhere based on detection limits
  filter(
    Site == "ETAT-3",
    Type == "Experimental"
  ) %>%
  summarise(
    Site = "ETAT-3",
    `Time Point` = "Average",
    `Lipid Δ²H (ppm)` = "0",
    `μ (year⁻¹)` = paste0(
      "< ",
      round(max(growth_wt_mu_year, na.rm = TRUE), 4),
      " ± ",
      round(max(growth_wt_mu_year_sem, na.rm = TRUE), 3)
    ),
    `TG (years)` = paste0(
      "> ",
      round(35),
      " ± ",
      round(5)
    ),
    .groups = "drop"
  )

table2_summary <- bind_rows(
  table2_beryl_timepoints,
  table2_beryl_average,
  table2_etat3
)

table2_summary

write.csv(
  table2_summary,
  "summary_tables/Table2_growth_summary_by_timepoint.csv",
  row.names = FALSE
)
