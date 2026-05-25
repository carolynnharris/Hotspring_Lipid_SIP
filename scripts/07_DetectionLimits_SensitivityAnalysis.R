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
# 4. Creates the Figure S5: sensitivity analysis showing incubation time
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
# - summary_tables/FigureS5_detection_limit_data.csv
# - summary_tables/FigureS5_detection_limit_value.csv
# - summary_tables/FigureS6_sensitivity_analysis_data.csv
# - figures/FigureS5_detection_limits.png
# - figures/FigureS6_sensitivity_analysis.png
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
sigma_F_lipid <- 0.5

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
  Ft_minus_F0_ppm = c(
    seq(0.1, 4.9, by = 0.005),
    seq(5, 100, by = 1)
  )
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

#### print detection limits ####
cat("Minimum detectable lipid 2H enrichment:",
    round(detection_limit_ppm, 4), "ppm\n")

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
  file.path(summary_dir, "FigureS5_detection_limit_data.csv"),
  row.names = FALSE
)

write.csv(
  detection_limit,
  file.path(summary_dir, "FigureS5_detection_limit_value.csv"),
  row.names = FALSE
)

###### plot detection limits  ###### 
 
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
  
  annotate(
    "text",
    x = detection_limit_ppm * 1.18,
    y = 0.1,
    label = "Detection limit",
    color = "red",
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = 0.13,
    y = 2.4,
    label = "95% CI",
    hjust = 0,
    size = 4
  ) +
  
  scale_x_log10(
    limits = c(0.1, 100),
    breaks = c(0.1, 1, round(detection_limit_ppm,1), 10, 100)
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
  file.path(figure_dir, "FigureS5_detection_limits.png"),
  plot = p_detection_limit,
  width = 5,
  height = 4,
  dpi = 300
)


#### Sensitivity Analysis ####
# calculate detectable generation times using stronger labeling solutions (FL)
# show 3-14 day incubations on the plots

###### define error assumptions ######
sigma_FL <- 100
sigma_F0 <- 0.5
sigma_Ft <- 0.5
sigma_alpha <- 0.1


###### define scenarios, calculate growth rate, T_G, and errors ######
sensitivity_df <- expand_grid(
  incubation_time_days = c(
    seq(1, 9, by = 1),
    seq(10, 1000, by = 10)
  ),
  alpha = c(0.56, 0.76),
  FL_ppm = c(3700, 10000, 100000),
  deltaH_ppm = c(1.56, 15.6, 156, 1555)
) %>%
  mutate(
    F0_ppm = 95,
    Ft_ppm = F0_ppm + deltaH_ppm,
    
    growth_rate_day = ifelse(
      Ft_ppm < alpha * FL_ppm,
      -(1 / incubation_time_days) *
        log((Ft_ppm - alpha * FL_ppm) / (F0_ppm - alpha * FL_ppm)),
      NA_real_
    ),
    
    growth_rate_year = growth_rate_day * 365,
    
    TG_days = log(2) / growth_rate_day,
    TG_years = TG_days / 365,
    
    growth_rate_error_day = sqrt(
      (
        sigma_F0^2 * (alpha * FL_ppm - Ft_ppm)^2 +
          sigma_Ft^2 * (F0_ppm - alpha * FL_ppm)^2 +
          (F0_ppm - Ft_ppm)^2 *
          (FL_ppm^2 * sigma_alpha^2 + alpha^2 * sigma_FL^2)
      ) /
        (
          incubation_time_days^2 *
            (F0_ppm - alpha * FL_ppm)^2 *
            (alpha * FL_ppm - Ft_ppm)^2
        )
    ),
    
    growth_rate_error_year = growth_rate_error_day * 365,
    
    TG_error_days = (log(2) / growth_rate_day^2) *
      growth_rate_error_day,
    
    TG_error_years = TG_error_days / 365,
    
    # Signal-to-error ratio
    growth_rate_SNR = growth_rate_year /
      growth_rate_error_year,
    
    # Detection limit classification
    Detection_Limit = case_when(
      growth_rate_SNR >= 2 ~ "Above",
      growth_rate_SNR < 2 ~ "BDL",
      TRUE ~ NA_character_
    )
  )

# save results
write.csv(
  sensitivity_df,
  file.path(summary_dir, "FigureS6_sensitivity_analysis_data.csv"),
  row.names = FALSE
)


##### plot sensitivity analysis ####

# make subplot labels
panel_labels <- tibble(
  FL_ppm = c(3700, 10000, 100000),
  label_panel = factor(
    c("0.37 at% ²H", "1 at% ²H", "10 at% ²H"),
    levels = c("0.37 at% ²H", "1 at% ²H", "10 at% ²H")
  ),
  panel_letter = c("A", "B", "C"),
  x = 840,
  y = 0.002
)

# add plotting labels
sensitivity_df_plot <- sensitivity_df %>%
  mutate(
    label_panel = case_when(
      FL_ppm == 3700 ~ "0.37 at% ²H",
      FL_ppm == 10000 ~ "1 at% ²H",
      FL_ppm == 100000 ~ "10 at% ²H"
    ),
    label_panel = factor(
      label_panel,
      levels = c("0.37 at% ²H", "1 at% ²H", "10 at% ²H")
    ),
    deltaH_label = paste0(round(deltaH_ppm, 1), " ppm"),
    alpha_label = as.character(alpha)
  )

# make plot
p_sensitivity <- ggplot(
  sensitivity_df_plot,
  aes(
    x = incubation_time_days,
    y = TG_years,
    color = factor(deltaH_ppm),
    linetype = alpha_label
  )
) +
  geom_rect(
    aes(
      xmin = 3,
      xmax = 14,
      ymin = 0.001,
      ymax = 1000
    ),
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.3
  ) +
  annotate(
    "text",
    x = 7,
    y = 0.002,
    label = "3 to 14 days",
    color = "black",
    size = 3.6
  ) +
  geom_text(
    data = panel_labels,
    aes(
      x = x,
      y = y,
      label = panel_letter
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 6
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ label_panel, nrow = 1) +
  scale_x_log10(
    limits = c(1, 1000),
    breaks = c(1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  scale_y_log10(
    limits = c(1e-3, 1e3),
    breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    labels = scales::label_math(.x)
  ) +
  scale_color_manual(
    values = c(
      "1.56" = "blue",
      "15.6" = "purple",
      "156" = "magenta",
      "1555" = "red"
    ),
    name = expression(Delta^2 * H),
    labels = c(
      "1.56" = "1.56 ppm (10 ‰ )",
      "15.6" = "15.6 ppm (100 ‰ )",
      "156" = "156 ppm (1,000 ‰ )",
      "1555" = "1555 ppm (10,000 ‰ )"
    )
  ) +
  scale_linetype_manual(
    values = c("0.56" = "solid", "0.76" = "dashed"),
    name = expression("Assimilation Efficiency ("*alpha*")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = expression("Sensitivity analysis across "^{2}*H[2]*"O labeling solutions"),
    x = "Incubation Time (days)",
    y = "Generation Time (years)"
  )

p_sensitivity

ggsave(
  file.path(figure_dir, "FigureS6_sensitivity_analysis.png"),
  plot = p_sensitivity,
  width = 11,
  height = 4.25,
  dpi = 300
)


###### print stats for fig caption ######

# Panel A = 0.37 at% 2H (FL = 3700)
# Panel C = 10 at% 2H (FL = 100000)
# Restrict to 3-14 day incubation window and detectable results

caption_df <- sensitivity_df %>%
  filter(
    incubation_time_days >= 3,
    incubation_time_days <= 14,
    Detection_Limit == "Above"
  )

##### Panel A #####

# heterotrophic (alpha = 0.56)
A_hetero <- caption_df %>%
  filter(
    FL_ppm == 3700,
    alpha == 0.56
  )

A_hetero_min <- min(A_hetero$TG_years, na.rm = TRUE)
A_hetero_max <- max(A_hetero$TG_years, na.rm = TRUE)

# autotrophic (alpha = 0.76)
A_auto <- caption_df %>%
  filter(
    FL_ppm == 3700,
    alpha == 0.76
  )

A_auto_min <- min(A_auto$TG_years, na.rm = TRUE)
A_auto_max <- max(A_auto$TG_years, na.rm = TRUE)


##### Panel B #####

# heterotrophic (alpha = 0.56)
B_hetero <- caption_df %>%
  filter(
    FL_ppm == 10000,
    alpha == 0.56
  )

B_hetero_min <- min(B_hetero$TG_years, na.rm = TRUE)
B_hetero_max <- max(B_hetero$TG_years, na.rm = TRUE)

# autotrophic (alpha = 0.76)
B_auto <- caption_df %>%
  filter(
    FL_ppm == 10000,
    alpha == 0.76
  )

B_auto_min <- min(B_auto$TG_years, na.rm = TRUE)
B_auto_max <- max(B_auto$TG_years, na.rm = TRUE)


##### Panel C #####

# heterotrophic (alpha = 0.56)
C_hetero <- caption_df %>%
  filter(
    FL_ppm == 100000,
    alpha == 0.56
  )

C_hetero_min <- min(C_hetero$TG_years, na.rm = TRUE)
C_hetero_max <- max(C_hetero$TG_years, na.rm = TRUE)

# autotrophic (alpha = 0.76)
C_auto <- caption_df %>%
  filter(
    FL_ppm == 100000,
    alpha == 0.76
  )

C_auto_min <- min(C_auto$TG_years, na.rm = TRUE)
C_auto_max <- max(C_auto$TG_years, na.rm = TRUE)


##### convert short generation times to days #####

A_hetero_min_days <- A_hetero_min * 365
A_auto_min_days <- A_auto_min * 365

B_hetero_min_days <- B_hetero_min * 365
B_auto_min_days <- B_auto_min * 365

C_hetero_min_days <- C_hetero_min * 365
C_auto_min_days <- C_auto_min * 365


##### print results #####

cat(
  "\n--------------------------------------------------\n",
  "Figure caption values:\n",
  "--------------------------------------------------\n\n",
  
  "A: 0.37 at% ²H labeling solution\n",
  "  Heterotrophic (α = 0.56): detectable generation times = ~",
  round(A_hetero_min_days,0),
  " days to ~",
  round(A_hetero_max,0),
  " years\n",
  
  "  Autotrophic (α = 0.76): detectable generation times = ~",
  round(A_auto_min_days,0),
  " days to ~",
  round(A_auto_max,0),
  " years\n\n",
  
  "B: 1 at% ²H labeling solution\n",
  "  Heterotrophic (α = 0.56): detectable generation times = ~",
  round(B_hetero_min_days,0),
  " days to ~",
  round(B_hetero_max,0),
  " years\n",
  
  "  Autotrophic (α = 0.76): detectable generation times = ~",
  round(B_auto_min_days,0),
  " days to ~",
  round(B_auto_max,0),
  " years\n\n",
  
  "C: 10 at% ²H labeling solution\n",
  "  Heterotrophic (α = 0.56): detectable generation times = ~",
  round(C_hetero_min_days,0),
  " days to ~",
  round(C_hetero_max,0),
  " years\n",
  
  "  Autotrophic (α = 0.76): detectable generation times = ~",
  round(C_auto_min_days,0),
  " days to ~",
  round(C_auto_max,0),
  " years\n\n",
  
  "--------------------------------------------------\n"
)


