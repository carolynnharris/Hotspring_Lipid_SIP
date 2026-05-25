cat("\014")
rm(list = ls())
graphics.off()

################################################################################
# 04_Lipid_2H_Uptake_Figures.R
#
# This script:
# 1. Loads the summary tables created in scripts 02 and 03.
# 2. Plots combo results figure (Fig 4):
#    A-B: Δ²H weighted mean for experimental + controls
#    C-D: Δ²H for BP-0, BP-1, BP-2, and weighted mean
#    E-F: integrated generation times
# 3. Plots uptake results in delta space (Fig S3):
#    mean d²HBP through time for experimental + controls.
################################################################################

#### load packages ####
library(tidyverse)

#### paths ####
summary_dir <- "summary_tables"
figure_dir <- "figures"

#### load data ####
dat <- read.csv(file.path(summary_dir, "LHSIP_Summary_Table.csv"))

integrated_growth <- read.csv(
  file.path(summary_dir, "LHSIP_Growth_Rates_Integrated_3_14_days.csv")
)

#### symbology ####
dat <- dat %>%
  mutate(
    col_site = case_when(
      Type != "Experimental" ~ "grey50",
      Site == "Beryl" ~ "green3",
      Site == "ETAT-3" ~ "darkblue",
      TRUE ~ "black"
    ),
    pch_type = case_when(
      Type == "Experimental" ~ 16,
      Type == "Control-Live" ~ 24,
      Type == "Control-Kill" ~ 25,
      TRUE ~ 16
    )
  )

line_types <- c(
  "Experimental" = 1,
  "Control-Live" = 2,
  "Control-Kill" = 3
)

line_widths <- c(
  "Experimental" = 2.5,
  "Control-Live" = 1.5,
  "Control-Kill" = 1.5
)

compound_pch <- c(
  "Wt. Mean" = 16,
  "BP-0" = 21,
  "BP-1" = 24,
  "BP-2" = 22
)

#### helper: plot weighted mean Δ2H through time ####
plot_delta_wt <- function(site, panel_letter) {
  
  df_site <- dat %>%
    filter(Site == site) %>%
    arrange(desc(Type), Incubation_days)
  
  plot(
    df_site$Incubation_days,
    df_site$Delta_D_ppm_wt_mean,
    las = 1,
    xlab = "Incubation Time (days)",
    ylab = expression(Delta^2 * H ~ "(ppm)"),
    ylim = c(-2, 6),
    xlim = c(0, 15),
    cex = 1.5,
    lwd = 1,
    col = df_site$col_site,
    bg = "white",
    pch = df_site$pch_type,
    cex.axis = 0.8,
    cex.lab = 1.1,
    xaxt = "n",
    yaxt = "n"
  )
  # lower detection limit
  polygon(
    x = c(-1, 16, 16, -1),
    y = c(-1.5, -1.5, 1.5, 1.5),
    col = rgb(0, 0, 0, alpha = 0.1),
    border = NA
  )
  
  abline(h = 0)
  
  for (type_i in unique(df_site$Type)) {
    subset_data <- df_site %>%
      filter(Type == type_i) %>%
      arrange(Incubation_days)
    
    lines(
      subset_data$Incubation_days,
      subset_data$Delta_D_ppm_wt_mean,
      col = ifelse(type_i == "Experimental", df_site$col_site[1], "grey50"),
      lty = line_types[type_i],
      lwd = line_widths[type_i]
    )
  }
  
  arrows(
    x0 = df_site$Incubation_days,
    y0 = df_site$Delta_D_ppm_wt_mean - df_site$Delta_D_ppm_wt_mean_sem,
    y1 = df_site$Delta_D_ppm_wt_mean + df_site$Delta_D_ppm_wt_mean_sem,
    angle = 90,
    code = 0,
    col = df_site$col_site,
    lwd = 1.5
  )
  
  points(
    df_site$Incubation_days,
    df_site$Delta_D_ppm_wt_mean,
    cex = 1.5,
    lwd = 1,
    col = df_site$col_site,
    bg = "white",
    pch = df_site$pch_type
  )
  
  mtext(panel_letter, side = 3, line = -12.5, cex = 1, font = 2, adj = 0.05)
  mtext(site, side = 3, line = 0.4, cex = 1, font = 2)
  
  axis(1, at = seq(0, 15, 5), cex.axis = 0.8, tck = -0.035)
  axis(2, at = seq(-2, 6, 2), las = 1, cex.axis = 0.8, tck = -0.035)
  
  legend(
    "topleft",
    inset = c(0.01, 0.01),
    legend = c("Experimental", "Control (Live)", "Control (Kill)"),
    pch = c(16, 24, 25),
    lty = c(1, 2, 3),
    lwd = c(2, 1, 1),
    col = c("black", "grey50", "grey50"),
    pt.bg = c("white", "white", "white"),
    bty = "n",
    cex = 0.85
  )
}

#### helper: plot BP-specific Δ2H through time ####
plot_delta_bp <- function(site, panel_letter) {
  
  df_site <- dat %>%
    filter(Site == site, Type == "Experimental") %>%
    arrange(Incubation_days)
  
  site_col <- unique(df_site$col_site)
  
  plot(
    df_site$Incubation_days,
    df_site$Delta_D_ppm_wt_mean,
    las = 1,
    xlab = "Incubation Time (days)",
    ylab = expression(Delta^2 * H ~ "(ppm)"),
    ylim = c(-2, 6),
    xlim = c(0, 15),
    cex = 1.5,
    lwd = 1.5,
    col = site_col,
    bg = "white",
    pch = 16,
    cex.axis = 0.8,
    cex.lab = 1.1,
    xaxt = "n",
    yaxt = "n"
  )
  # lower detection limit
  polygon(
    x = c(-1, 16, 16, -1),
    y = c(-1.5, -1.5, 1.5, 1.5),
    col = rgb(0, 0, 0, alpha = 0.1),
    border = NA
  )
  
  abline(h = 0)
  
  bp_plot_list <- tribble(
    ~label,     ~mean_col,                 ~sem_col,                     ~pch, ~lwd,
    "Wt. Mean", "Delta_D_ppm_wt_mean",      "Delta_D_ppm_wt_mean_sem",     16,   2.5,
    "BP-0",     "Delta_BP0_D_ppm_mean",     "Delta_BP0_D_ppm_mean_sem",    21,   1.5,
    "BP-1",     "Delta_BP1_D_ppm_mean",     "Delta_BP1_D_ppm_mean_sem",    24,   1.5,
    "BP-2",     "Delta_BP2_D_ppm_mean",     "Delta_BP2_D_ppm_mean_sem",    22,   1.5
  )
  
  for (i in seq_len(nrow(bp_plot_list))) {
    mean_col <- bp_plot_list$mean_col[i]
    sem_col <- bp_plot_list$sem_col[i]
    
    lines(
      df_site$Incubation_days,
      df_site[[mean_col]],
      col = site_col,
      lty = 1,
      lwd = bp_plot_list$lwd[i]
    )
    
    arrows(
      x0 = df_site$Incubation_days,
      y0 = df_site[[mean_col]] - df_site[[sem_col]],
      y1 = df_site[[mean_col]] + df_site[[sem_col]],
      angle = 90,
      code = 0,
      col = site_col,
      lwd = 1.5
    )
    
    points(
      df_site$Incubation_days,
      df_site[[mean_col]],
      cex = 1.5,
      lwd = 1.5,
      col = site_col,
      bg = "white",
      pch = bp_plot_list$pch[i]
    )
  }
  
  mtext(panel_letter, side = 3, line = -12.5, cex = 1, font = 2, adj = 0.05)
  
  axis(1, at = seq(0, 15, 5), cex.axis = 0.8, tck = -0.035)
  axis(2, at = seq(-2, 6, 2), las = 1, cex.axis = 0.8, tck = -0.035)
  
  legend(
    "topleft",
    inset = c(0.01, 0.01),
    legend = c("Wt. Mean", "BP0", "BP1", "BP2"),
    pch = c(16, 21, 24, 22),
    lty = 1,
    col = "black",
    pt.bg = "white",
    lwd = c(2.5, 1, 1, 1),
    bty = "n",
    cex = 0.85
  )
}

#### prepare generation-time panel data ####
beryl_tg_vals <- integrated_growth %>%
  filter(Site == "Beryl") %>%
  slice(1)

beryl_tg <- tibble(
  Site = "Beryl",
  Compound = c("BP-0", "BP-1", "BP-2", "Wt. Mean"),
  TG = c(
    beryl_tg_vals$T_G_years_BP0,
    beryl_tg_vals$T_G_years_BP1,
    beryl_tg_vals$T_G_years_BP2,
    beryl_tg_vals$T_G_years_wt
  ),
  TG_sem = c(
    beryl_tg_vals$T_G_years_BP0_sem,
    beryl_tg_vals$T_G_years_BP1_sem,
    beryl_tg_vals$T_G_years_BP2_sem,
    beryl_tg_vals$T_G_years_wt_sem
  ),
  Col = "green3",
  pch = c(21, 24, 22, 16)
)

etat_tg <- tibble(
  Site = "ETAT-3",
  Compound = c("BP-0", "BP-1", "BP-2", "Wt. Mean"),
  TG = 42,
  TG_sem = 21,
  Col = "darkblue",
  pch = c(21, 24, 22, 16)
)

tg_dat <- bind_rows(beryl_tg, etat_tg) %>%
  mutate(
    Compound = factor(Compound, levels = c("BP-0", "BP-1", "BP-2", "Wt. Mean")),
    y_values = as.numeric(Compound)
  )

#### helper: plot generation time panels ####
plot_tg_panel <- function(site, panel_letter) {
  
  subset <- tg_dat %>%
    filter(Site == site)
  
  plot(
    subset$TG,
    subset$y_values,
    las = 1,
    log = "x",
    xlab = "Generation Time (years)",
    ylab = "",
    ylim = c(0.5, 4.5),
    xlim = c(1, 1000),
    cex = 1.5,
    lwd = 1.5,
    col = subset$Col,
    bg = "white",
    pch = subset$pch,
    cex.axis = 0.8,
    cex.lab = 1.1,
    xaxt = "n",
    yaxt = "n"
  )
  # detection limits 
  polygon(
    x = c(42, 2000, 2000, 42),
    y = c(0, 0, 10, 10),
    col = rgb(0, 0, 0, alpha = 0.08),
    border = NA
  )
  
  mtext("Compound", side = 2, line = 2.5, cex = 0.8)
  mtext("Detectable Growth", col = "grey30", side = 1, font = 3, line = -13, adj = 0.05, cex = 0.7)
  mtext("Undetectable Growth", col = "grey30", side = 1, font = 3, line = -13, adj = 0.95, cex = 0.7)
  
  if (site == "Beryl") {
    arrows(
      y0 = subset$y_values,
      x0 = subset$TG - subset$TG_sem,
      x1 = subset$TG + subset$TG_sem,
      angle = 90,
      code = 0,
      col = subset$Col,
      lwd = 1.5
    )
  } else {
    arrows(
      y0 = subset$y_values,
      x0 = subset$TG - subset$TG_sem,
      x1 = 1000,
      angle = 90,
      code = 0,
      col = subset$Col,
      lwd = 1.5
    )
  }
  
  points(
    subset$TG,
    subset$y_values,
    cex = 1.5,
    lwd = 1.5,
    col = subset$Col,
    bg = "white",
    pch = subset$pch
  )
  
  mtext(panel_letter, side = 3, line = -12.5, cex = 1, font = 2, adj = 0.05)
  
  axis(1, at = c(1, 10, 100, 1000), labels = TRUE, cex.axis = 0.8, tck = -0.035)
  axis(
    2,
    at = 1:4,
    labels = c("BP-0", "BP-1", "BP-2", "Wt. Mean"),
    las = 1,
    cex.axis = 0.8,
    tck = -0.035
  )
}

#### Figure 4 composite ####
png(
  filename = file.path(figure_dir, "Figure4_Lipid_2H_Uptake_Composite.png"),
  width = 7,
  height = 7,
  units = "in",
  res = 300
)

par(
  mfrow = c(3, 2),
  mar = c(2.5, 3.5, 1, 1),
  oma = c(0.5, 1, 1, 1),
  mgp = c(1.4, 0.5, 0)
)

plot_delta_wt("Beryl", "A")
plot_delta_wt("ETAT-3", "B")
plot_delta_bp("Beryl", "C")
plot_delta_bp("ETAT-3", "D")
plot_tg_panel("Beryl", "E")
plot_tg_panel("ETAT-3", "F")

dev.off()

#### Figure S3: mean δ2HBP through time ####
png(
  filename = file.path(figure_dir, "FigureS3_Lipid_2H_Uptake_deltaSpace.png"),
  width = 7 * 1.3,
  height = 3 * 1.2,
  units = "in",
  res = 300
)

par(
  mfrow = c(1, 2),
  mar = c(3, 3.5, 2, 1),
  oma = c(1, 1, 1, 1),
  mgp = c(1.9, 0.5, 0)
)

for (site in c("Beryl", "ETAT-3")) {
  
  df_site <- dat %>%
    filter(Site == site) %>%
    arrange(desc(Type), Incubation_days)
  
  ymin <- round(min(df_site$d2H_wt_mean_mean, na.rm = TRUE), -1) - 10
  ymax <- round(max(df_site$d2H_wt_mean_mean, na.rm = TRUE), -1) + 10
  
  plot(
    df_site$Incubation_days,
    df_site$d2H_wt_mean_mean,
    las = 1,
    xlab = "Incubation Time (days)",
    ylab = expression("Mean " * delta^2 * H["BP"] ~ "(" * "\u2030" * ")"),
    ylim = c(ymin, ymax),
    xlim = c(0, 15),
    cex = 1.5,
    lwd = 1,
    col = df_site$col_site,
    bg = "white",
    pch = df_site$pch_type,
    cex.axis = 0.8,
    cex.lab = 1.1,
    xaxt = "n"
  )
  
  for (type_i in unique(df_site$Type)) {
    subset_data <- df_site %>%
      filter(Type == type_i) %>%
      arrange(Incubation_days)
    
    lines(
      subset_data$Incubation_days,
      subset_data$d2H_wt_mean_mean,
      col = ifelse(type_i == "Experimental", df_site$col_site[1], "grey50"),
      lty = line_types[type_i],
      lwd = line_widths[type_i]
    )
  }
  
  arrows(
    x0 = df_site$Incubation_days,
    y0 = df_site$d2H_wt_mean_mean - df_site$d2H_wt_mean_sem,
    y1 = df_site$d2H_wt_mean_mean + df_site$d2H_wt_mean_sem,
    angle = 90,
    code = 0,
    col = df_site$col_site,
    lwd = 1.5
  )
  
  points(
    df_site$Incubation_days,
    df_site$d2H_wt_mean_mean,
    cex = 1.5,
    lwd = 1.5,
    col = df_site$col_site,
    bg = "white",
    pch = df_site$pch_type
  )
  
  mtext(ifelse(site == "Beryl", "A", "B"), side = 3, line = -10.5, cex = 1.1, font = 2, adj = 0.05)
  mtext(site, side = 3, line = 0.4, cex = 1, font = 2)
  
  axis(1, at = seq(0, 15, 5), cex.axis = 0.8, tck = -0.035)
  
  legend(
    "topleft",
    inset = c(0.01, 0.01),
    legend = c("Experimental", "Control (Live)", "Control (Kill)"),
    pch = c(16, 24, 25),
    lty = c(1, 2, 3),
    lwd = c(2, 1, 1),
    col = c("black", "grey50", "grey50"),
    pt.bg = c("white", "white", "white"),
    bty = "n",
    cex = 0.85
  )
}

dev.off()

