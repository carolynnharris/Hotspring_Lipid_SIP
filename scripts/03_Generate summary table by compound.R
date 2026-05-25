cat("\014")
rm(list=ls()) 
graphics.off()

################################################################################
# 03_Generate_summary_table_by_compound.R
#
# This script:
# 1. Pulls compound-specific IPL BP standing stocks for Beryl
# 2. Pulls compound-specific growth rates TG for BP-0, BP-1, and BP-2 (time-integrated)
# 3. Calculates annual biomass production rates for each BP
# 4. Calculates the amount of new biomass produced during incubation
# 5. Calculates the percent increase in biomass over the incubation.
# 6. Calculates abundance-weighted mean growth and production metrics across BP pools
# 7. Formats and saves the Table 3 summary table
# 9. Plots BP ring number vs. TG (Fig S4) 
#
# Outputs:
# - summary_tables/Table3_beryl_production_raw.csv
# - summary_tables/Table3_beryl_production_formatted.csv
# - figures/FigureS4_BP_ring_number_generation_time.png
################################################################################

library(tidyverse)

#### settings ####
t_days <- 14
t_years <- t_days / 365
summary_dir <- "summary_tables"

#### load data ####
integrated_growth <- read.csv(
  file.path(summary_dir, "LHSIP_Growth_Rates_Integrated_3_14_days.csv")
)

lhsip_merge <- read.csv(
  file.path("data_inputs", "LHSIP_All_Merge.csv")
)

#### define Beryl initial IPL-BP standing stock ####

initial_bp_conc <- lhsip_merge %>%
  filter(ID_Fraction == "B01_IPL") %>%
  summarise(
    BP0_total = BP_ug_sediment_g * BP0_RelAbund/100,
    BP1_total = BP_ug_sediment_g * BP1_RelAbund/100,
    BP2_total = BP_ug_sediment_g * BP2_RelAbund/100
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "BP",
    values_to = "initial_conc_ug_g"
  ) %>%
  mutate(
    Compound = case_when(
      BP == "BP0_total" ~ "BP-0",
      BP == "BP1_total" ~ "BP-1",
      BP == "BP2_total" ~ "BP-2"
    ),
  ) %>%
  select(Compound, initial_conc_ug_g)



#### reshape integrated growth results ####
beryl_growth_long <- integrated_growth %>%
  filter(Site == "Beryl") %>%
  transmute(
    Site,
    BP0_mu = mu_year_BP0,
    BP0_mu_sem = mu_year_BP0_sem,
    BP0_TG = T_G_years_BP0,
    BP0_TG_sem = T_G_years_BP0_sem,
    
    BP1_mu = mu_year_BP1,
    BP1_mu_sem = mu_year_BP1_sem,
    BP1_TG = T_G_years_BP1,
    BP1_TG_sem = T_G_years_BP1_sem,
    
    BP2_mu = mu_year_BP2,
    BP2_mu_sem = mu_year_BP2_sem,
    BP2_TG = T_G_years_BP2,
    BP2_TG_sem = T_G_years_BP2_sem
  ) %>%
  pivot_longer(
    cols = -Site,
    names_to = c("BP", ".value"),
    names_pattern = "(BP[0-2])_(.*)"
  ) %>%
  mutate(
    Compound = case_when(
      BP == "BP0" ~ "BP-0",
      BP == "BP1" ~ "BP-1",
      BP == "BP2" ~ "BP-2",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Site, Compound, mu, mu_sem, TG, TG_sem)

#### calculate compound-specific production ####
#### calculate compound-specific production ####
BP_new_production <- beryl_growth_long %>%
  left_join(initial_bp_conc, by = "Compound") %>%
  mutate(
    production_rate_ug_g_yr = initial_conc_ug_g * mu,
    production_rate_ug_g_yr_sem = initial_conc_ug_g * mu_sem,
    
    production_rate_ng_g_yr = production_rate_ug_g_yr * 1000,
    production_rate_ng_g_yr_sem = production_rate_ug_g_yr_sem * 1000,
    
    new_biomass_ug_g =
      initial_conc_ug_g * (exp(mu * t_years) - 1),
    
    new_biomass_ug_g_sem =
      initial_conc_ug_g * t_years * exp(mu * t_years) * mu_sem,
    
    new_biomass_ng_g = new_biomass_ug_g * 1000,
    new_biomass_ng_g_sem = new_biomass_ug_g_sem * 1000,
    
    percent_increase =
      100 * new_biomass_ug_g / initial_conc_ug_g,
    
    percent_increase_sem =
      100 * new_biomass_ug_g_sem / initial_conc_ug_g
  )

#### calculate abundance-weighted mean row ####
new_production_wt_mean <- BP_new_production %>%
  summarise(
    Site = "Beryl",
    Compound = "Wt. Average",
    
    mu = weighted.mean(mu, w = initial_conc_ug_g, na.rm = TRUE),
    mu_sem = sqrt(sum((initial_conc_ug_g * mu_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    TG = log(2) / mu,
    TG_sem = abs(log(2) / mu^2) * mu_sem,
    
    production_rate_ug_g_yr =
      weighted.mean(production_rate_ug_g_yr, w = initial_conc_ug_g, na.rm = TRUE),
    production_rate_ug_g_yr_sem =
      sqrt(sum((initial_conc_ug_g * production_rate_ug_g_yr_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    production_rate_ng_g_yr =
      weighted.mean(production_rate_ng_g_yr, w = initial_conc_ug_g, na.rm = TRUE),
    production_rate_ng_g_yr_sem =
      sqrt(sum((initial_conc_ug_g * production_rate_ng_g_yr_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    new_biomass_ug_g =
      weighted.mean(new_biomass_ug_g, w = initial_conc_ug_g, na.rm = TRUE),
    new_biomass_ug_g_sem =
      sqrt(sum((initial_conc_ug_g * new_biomass_ug_g_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    new_biomass_ng_g =
      weighted.mean(new_biomass_ng_g, w = initial_conc_ug_g, na.rm = TRUE),
    new_biomass_ng_g_sem =
      sqrt(sum((initial_conc_ug_g * new_biomass_ng_g_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    percent_increase =
      weighted.mean(percent_increase, w = initial_conc_ug_g, na.rm = TRUE),
    percent_increase_sem =
      sqrt(sum((initial_conc_ug_g * percent_increase_sem)^2, na.rm = TRUE)) /
      sum(initial_conc_ug_g, na.rm = TRUE),
    
    .groups = "drop"
  )

#### combine and format ####
table3_raw <- bind_rows(
  BP_new_production,
  new_production_wt_mean
)

table3_formatted <- table3_raw %>%
  transmute(
    Compound,
    `μ (year⁻¹)` = paste0(
      sprintf("%.3f", mu),
      " ± ",
      sprintf("%.2f", mu_sem)
    ),
    `TG (years)` = paste0(
      round(TG),
      " ± ",
      round(TG_sem)
    ),
    `Production rate (ng g⁻¹ sed year⁻¹)` = paste0(
      sprintf("%.1f", production_rate_ng_g_yr),
      " ± ",
      sprintf("%.1f", production_rate_ng_g_yr_sem)
    ),
    `New biomass (ng g⁻¹ sed)` = paste0(
      sprintf("%.2f", new_biomass_ng_g),
      " ± ",
      sprintf("%.2f", new_biomass_ng_g_sem)
    ),
    `% Increase` = paste0(
      sprintf("%.2f", percent_increase),
      " ± ",
      sprintf("%.2f", percent_increase_sem)
    )
  )

table3_formatted

#### save Table 3 #### 
write.csv(
  table3_raw,
  file.path(summary_dir, "Table3_beryl_production_raw.csv"),
  row.names = FALSE
)

write.csv(
  table3_formatted,
  file.path(summary_dir, "Table3_beryl_production_formatted.csv"),
  row.names = FALSE
)


#### Figure S4: BP ring number vs. generation time ####

BP_ring_df <- BP_new_production %>%
  filter(Compound %in% c("BP-0", "BP-1", "BP-2")) %>%
  mutate(
    BP_ring_number = case_when(
      Compound == "BP-0" ~ 0,
      Compound == "BP-1" ~ 1,
      Compound == "BP-2" ~ 2
    )
  )

BP_ring_lm <- lm(TG ~ BP_ring_number, data = BP_ring_df)

BP_ring_stats <- broom::tidy(BP_ring_lm)
BP_ring_glance <- broom::glance(BP_ring_lm)

BP_ring_slope <- BP_ring_stats %>%
  filter(term == "BP_ring_number") %>%
  pull(estimate)

BP_ring_slope_se <- BP_ring_stats %>%
  filter(term == "BP_ring_number") %>%
  pull(std.error)

BP_ring_r2 <- BP_ring_glance$r.squared

BP_ring_label <- paste0(
  "atop(",
  "'Slope = ", round(BP_ring_slope, 1), " ± ", round(BP_ring_slope_se, 1), " years ring'^-1,",
  "R^2 == ", round(BP_ring_r2, 2),
  ")"
)

p_BP_ring_lm <- ggplot(BP_ring_df, aes(x = BP_ring_number, y = TG)) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = TG - TG_sem, ymax = TG + TG_sem),
    width = 0,
    color = "green3",
    linewidth = 0.6
  ) +
  geom_point(
    color = "green3",
    size = 3
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 37,
    label = BP_ring_label,
    hjust = 0,
    vjust = 1,
    size = 3.6,
    parse = TRUE
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    limits = c(-0.35, 2.35)
  ) +
  scale_y_continuous(
    limits = c(0, 40),
    breaks = c(0, 20, 40)
  ) +
  theme_bw() +
  labs(
    x = "BP Ring Number",
    y = "Generation Time (years)"
  )

p_BP_ring_lm

ggsave(
  filename = "figures/FigureS4_BP_ring_number_generation_time.png",
  plot = p_BP_ring_lm,
  width = 5,
  height = 4,
  dpi = 300
)
