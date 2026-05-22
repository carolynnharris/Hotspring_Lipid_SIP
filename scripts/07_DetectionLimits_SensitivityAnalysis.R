cat("\014")
rm(list = ls())
graphics.off()

################################################################################
# 07_DetectionLimits_SensitivityAnalysis.R
#
# This script:
# 1. Calculates LH-SIP assay detection-limits
# 2. Uses propagated uncertainty to identify the minimum lipid 2H enrichment
#    where growth rate is distinguishable from zero at the 2σ level
# 3. Creates Figure S6: lipid 2H enrichment vs. μ / σ_μ
# 4. Creates the Figure S5 sensitivity analysis showing incubation time
#    needed to detect growth across generation times, label strengths,
#    assimilation efficiencies, and lipid 2H enrichments
#
# Assumptions:
# - Mean assimilation efficiency:
#       alpha = 0.66
# - Assimilation efficiency uncertainty:
#       sigma_alpha = 0.10
# - VSMOW D/H reference ratio: RVSMOW = 0.00015576
#
# Outputs:
# - summary_tables/FigureS6_detection_limit_data.csv
# - summary_tables/FigureS6_detection_limit_value.csv
# - summary_tables/FigureS5_sensitivity_analysis_data.csv
# - figures/FigureS6_detection_limits.png
# - figures/FigureS5_sensitivity_analysis.png
#
# Notes:
# - Figure S6 identifies the minimum detectable lipid enrichment under the
#   14-day, 0.37 at% 2H incubation conditions used in this study
# - Figure S5 evaluates how stronger labeling solutions and longer incubation
#   durations improve sensitivity to slow microbial growth
################################################################################

library(tidyverse)

#### paths ####
figure_dir <- "figures"
summary_dir <- "summary_tables"



#### constants ####
RVSMOW <- 0.00015576
alpha_mean <- 0.66
alpha_sigma <- 0.10

#### approximate uncertainties in F values ####
# Approximate analytical uncertainty in lipid F values ~ +/-3 permil
sigma_F_lipid <- 0.47

# Approximate tracer-water uncertainty ~ +/-500 permil
sigma_FL <- 100

#### helper functions ####
delta_to_F_ppm <- function(delta_permil) {
  R <- RVSMOW * (delta_permil / 1000 + 1)
  F <- R / (1 + R)
  F * 1e6
}

F_ppm_to_delta <- function(F_ppm) {
  F <- F_ppm / 1e6
  R <- F / (1 - F)
  ((R / RVSMOW) - 1) * 1000
}

# growth rate
growth_rate_day <- function(Ft_ppm, F0_ppm, FL_ppm, t_days, alpha) {
  A <- alpha * FL_ppm - Ft_ppm
  B <- alpha * FL_ppm - F0_ppm
  
  ifelse(
    A > 0 & B > 0 & Ft_ppm > F0_ppm,
    -(1 / t_days) * log(A / B),
    NA_real_
  )
}

# growth rate error
growth_rate_error_day <- function(Ft_ppm, F0_ppm, FL_ppm, t_days,
                                  alpha, sigma_Ft, sigma_F0,
                                  sigma_FL, alpha_sigmalpha) {
  
  A <- alpha * FL_ppm - Ft_ppm
  B <- alpha * FL_ppm - F0_ppm
  
  dmu_dFt <-  (1 / t_days) * (1 / A)
  dmu_dF0 <- -(1 / t_days) * (1 / B)
  dmu_dFL <- -(alpha / t_days) * ((1 / A) - (1 / B))
  dmu_da  <- -(FL_ppm / t_days) * ((1 / A) - (1 / B))
  
  sqrt(
    (dmu_dFt * sigma_Ft)^2 +
      (dmu_dF0 * sigma_F0)^2 +
      (dmu_dFL * sigma_FL)^2 +
      (dmu_da  * alpha_sigmalpha)^2
  )
}


#### Detection Limits ####

# LHSIP assay conditions in this study
F0_ppm <- 95 # Baseline natural-abundance lipid value
FL_detection_ppm <- 3700 # 0.37 at % 2H
t_detection_days <- 14 # 14-day incubations

detection_df <- tibble(
  Ft_minus_F0_ppm = seq(0.1, 100, length.out = 1000)
) %>%
  mutate(
    F0_ppm = F0_ppm,
    Ft_ppm = F0_ppm + Ft_minus_F0_ppm,
    FL_ppm = FL_detection_ppm,
    t_days = t_detection_days,
    alpha = alpha_mean,
    
    mu_day = growth_rate_day(
      Ft_ppm = Ft_ppm,
      F0_ppm = F0_ppm,
      FL_ppm = FL_ppm,
      t_days = t_days,
      alpha = alpha
    ),
    
    mu_day_sem = growth_rate_error_day(
      Ft_ppm = Ft_ppm,
      F0_ppm = F0_ppm,
      FL_ppm = FL_ppm,
      t_days = t_days,
      alpha = alpha,
      sigma_Ft = sigma_F_lipid,
      sigma_F0 = sigma_F_lipid,
      sigma_FL = sigma_FL,
      alpha_sigmalpha = alpha_sigma
    ),
    
    mu_sigma_ratio = mu_day / mu_day_sem,
    
    mu_year = mu_day * 365,
    mu_year_sem = mu_day_sem * 365,
    
    T_G_years = log(2) / mu_year,
    T_G_years_sem = abs(log(2) / mu_year^2) * mu_year_sem
  )

detection_limit <- detection_df %>%
  filter(mu_sigma_ratio >= 2) %>%
  slice(1)

detection_limit_ppm <- detection_limit$Ft_minus_F0_ppm

# convert ppm atom-fraction enrichment to delta units
F0_delta <- F_ppm_to_delta(F0_ppm)

Ft_detection_delta <- F_ppm_to_delta(
  F0_ppm + detection_limit_ppm
)

detection_limit_delta_permil <-
  Ft_detection_delta - F0_delta

cat("Minimum detectable lipid 2H enrichment:",
    round(detection_limit_ppm, 2), "ppm\n")

cat("Minimum detectable lipid δ2H  enrichment:",
    round(detection_limit_delta_permil, 2), "permil\n")

cat("Growth rate at detection limit:",
    round(detection_limit$mu_year, 4), "±",
    round(detection_limit$mu_year_sem, 4), "year^-1\n")

cat("Generation time at detection limit:",
    round(detection_limit$T_G_years, 1), "±",
    round(detection_limit$T_G_years_sem, 1), "years\n")

# save files
write.csv(
  detection_df,
  file.path(summary_dir, "FigureS6_detection_limit_data.csv"),
  row.names = FALSE
)

write.csv(
  detection_limit,
  file.path(summary_dir, "FigureS6_detection_limit_value.csv"),
  row.names = FALSE
)


# plot detection limit 
p_detection_limit <- ggplot(
  detection_df,
  aes(x = Ft_minus_F0_ppm, y = mu_sigma_ratio)
) +
  geom_line(
    color = "red",
    linewidth = 0.8
  ) +
  
  # 95% CI lines
  geom_segment(
    aes(
      x = 0.1,
      xend = detection_limit_ppm,
      y = 2,
      yend = 2
    ),
    color = "black",
    linewidth = 0.7,
    linetype = 2
  ) +
  geom_segment(
    aes(
      x = detection_limit_ppm,
      xend = detection_limit_ppm,
      y = 0.1,
      yend = 2
    ),
    color = "black",
    linewidth = 0.7,
    linetype = 2
  ) +
  
  # detection limit point
  geom_point(
    aes(
      x = detection_limit_ppm,
      y = 0.1
    ),
    shape = 18,
    size = 5,
    color = "red"
  ) +
  coord_cartesian(clip = "off") +
  
  annotate(
    "text",
    x = detection_limit_ppm * 1.18,
    y = 0.1,
    label = "Detection limit",
    color = "red",
    hjust = 0,
    size = 5
  ) +
  
  annotate(
    "text",
    x = 0.13,
    y = 2.2,
    label = "95% CI",
    hjust = 0,
    size = 5
  ) +
  
  scale_x_log10(
    limits = c(0.1, 100),
    breaks = c(0.1, 1, detection_limit_ppm, 10, 100)
  ) +
  
  scale_y_log10(
    limits = c(0.1, 10),
    breaks = c(0.1, 1, 2, 10)
  ) +
  
  theme_bw() +
  
  # major gridlines
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  ) +
  
  labs(
    title = "Detection Limits",
    x = expression(
      Delta^2*H~"(ppm) ="~F[t]-F[0]
    ),
    y = expression(mu / sigma[mu])
  )

p_detection_limit

ggsave(
  file.path(figure_dir, "FigureS6_detection_limits.png"),
  plot = p_detection_limit,
  width = 5,
  height = 4,
  dpi = 300
)


#### Sensitivity Analysis ####
# calculate detectable generation times using stronger labeling solutions (FL)
# show 3-14 day incubations on the plots

label_strengths <- tibble(
  label_panel = factor(
    c("0.37 at% ²H", "1 at% ²H", "10 at% ²H"),
    levels = c("0.37 at% ²H", "1 at% ²H", "10 at% ²H")
  ),
  FL_ppm = c(3700, 10000, 100000)
)

alpha_scenarios <- tibble(
  alpha_label = factor(
    c("0.56", "0.76"),
    levels = c("0.56", "0.76")
  ),
  alpha = c(0.56, 0.76)
)

deltaH_scenarios <- tibble(
  deltaH_permil = c(10, 100, 1000, 10000),
  deltaH_ppm = delta_to_F_ppm(deltaH_permil) - delta_to_F_ppm(0)
)

generation_times_years <- 10^seq(-3, 3, length.out = 600)

sensitivity_df <- expand_grid(
  label_strengths,
  alpha_scenarios,
  deltaH_scenarios,
  generation_time_years = generation_times_years
) %>%
  mutate(
    mu_day = log(2) / (generation_time_years * 365),
    
    # Solve LH-SIP equation for incubation time:
    # ΔF = (1 - exp(-μt)) * (alpha * FL - F0)
    # t = -ln(1 - ΔF / (alpha * FL - F0)) / μ
    detectable_fraction = deltaH_ppm / (alpha * FL_ppm - F0_ppm),
    
    incubation_time_days = ifelse(
      detectable_fraction > 0 & detectable_fraction < 1,
      -log(1 - detectable_fraction) / mu_day,
      NA_real_
    ),
    
    TG_error_years = NA_real_
  )

write.csv(
  sensitivity_df,
  file.path(summary_dir, "FigureS5_sensitivity_analysis_data.csv"),
  row.names = FALSE
)

p_sensitivity <- ggplot(
  sensitivity_df,
  aes(
    x = generation_time_years,
    y = incubation_time_days,
    color = factor(deltaH_permil),
    linetype = alpha_label
  )
) +
  geom_hline(
    yintercept = 3,
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = 14,
    color = "black",
    linewidth = 0.6
  ) +
  
  annotate(
    "text",
    x = 0.01,
    y = 5.5,
    label = "3 to 14 days",
    color = "black",
    size = 4
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ label_panel, ncol = 1) +
  scale_x_log10(
    limits = c(1e-3, 1e3),
    breaks = c(1e-3, 1e-1, 1e1, 1e3),
    labels = scales::label_math(.x)
  ) +
  scale_y_log10(
    limits = c(1e-2, 1e2),
    breaks = c(1e-2, 1e-1, 1, 10, 100),
    labels = scales::label_math(.x)
  ) +
  scale_color_manual(
    values = c(
      "10" = "blue",
      "100" = "purple",
      "1000" = "magenta",
      "10000" = "red"
    ),
    name = expression(Delta^2 * H ~ "(\u2030)")
  ) +
  scale_linetype_manual(
    values = c("0.56" = "solid", "0.76" = "dashed"),
    name = expression(alpha)
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Generation Time (years)",
    y = "Incubation Time (days)"
  )

p_sensitivity

ggsave(
  file.path(figure_dir, "FigureS5_sensitivity_analysis.png"),
  plot = p_sensitivity,
  width = 6,
  height = 7,
  dpi = 300
)