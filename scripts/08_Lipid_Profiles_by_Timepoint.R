cat("\014")
rm(list = ls())
graphics.off()

################################################################################
# 08_Lipid_Profiles_by_Timepoint.R
#
# This script:
# 1. Loads the merged LH-SIP dataset
# 2. Summarizes IPL-iGDGT and IPL-derived biphytane lipid distributions by
#    Site x Type x incubation time point
# 3. Calculates mean RI-GDGT and RI-BP values across biological replicates
# 4. Generates a multi-panel lipid profile figure showing:
#      - Ring Index through time
#      - iGDGT relative abundance distributions
#      - BP relative abundance distributions
#
# Output:
# - figures/FigureS2_lipid_profiles_by_timepoint.png
################################################################################

library(tidyverse)
library(patchwork)

#### paths ####
data_dir <- "data_inputs"
figure_dir <- "figures"

#### symbology ####
GDGT_bar_colors <- c(
  "#00EEEE", "#1DD2E1", "#3BB7D4", "#599BC8", "#7780BB",
  "grey50", "#9464AE", "#B249A2", "#D02D95", "#EE1289"
)

BP_bar_colors <- c(
  "#00EEEE", "#599BC8", "#9464AE", "#B249A2", "#EE1289"
)

#### load data ####
dat <- read.csv(file.path(data_dir, "LHSIP_All_Merge.csv")) %>%
  filter(Fraction == "IPL") %>%
  mutate(
    Incubation_days = round(Incubation_days),
    Type_label = case_when(
      Type == "Experimental" ~ "Exp.",
      Type == "Control-Live" ~ "Cnt-Live",
      Type == "Control-Kill" ~ "Cnt-Kill",
      TRUE ~ Type
    ),
    Type_label = factor(
      Type_label,
      levels = c("Exp.", "Cnt-Live", "Cnt-Kill")
    )
  )

#### summarize by site, type, and time point ####
lipid_summary <- dat %>%
  group_by(Site, Type_label, Incubation_days) %>%
  summarise(
    across(
      c(
        GDGT0_RelAbund, GDGT1_RelAbund, GDGT2_RelAbund,
        GDGT3_RelAbund, GDGT4_RelAbund, Cren_RelAbund,
        GDGT5_RelAbund, GDGT6_RelAbund, GDGT7_RelAbund,
        GDGT8_RelAbund,
        BP0_RelAbund, BP1_RelAbund, BP2_RelAbund,
        BP3_RelAbund, BP4_RelAbund,
        RI_GDGT, RI_BP
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .))) %>%
  mutate(
    sample_label = paste0(
      case_when(
        Type_label == "Exp." ~ "Exp",
        Type_label == "Cnt-Live" ~ "Cnt_L",
        Type_label == "Cnt-Kill" ~ "Cnt_K"
      ),
      "_",
      Incubation_days,
      "d"
    ),
    sample_label = factor(
      sample_label,
      levels = c(
        "Exp_0d", "Exp_3d", "Exp_14d",
        "Cnt_L_0d", "Cnt_L_14d",
        "Cnt_K_0d", "Cnt_K_14d"
      )
    )
  )

#### RI data ####
ri_long <- lipid_summary %>%
  select(Site, Type_label, Incubation_days, RI_GDGT, RI_BP) %>%
  pivot_longer(
    cols = c(RI_GDGT, RI_BP),
    names_to = "RI_type",
    values_to = "Ring_Index"
  ) %>%
  mutate(
    RI_type = recode(
      RI_type,
      RI_GDGT = "RI-GDGT",
      RI_BP = "RI-BP"
    ),
    RI_type = factor(RI_type, levels = c("RI-GDGT", "RI-BP"))
  )

#### GDGT data: faceted by Type, x-axis = Incubation_days ####
gdgt_long <- lipid_summary %>%
  select(
    Site, Type_label, Incubation_days,
    GDGT0_RelAbund, GDGT1_RelAbund, GDGT2_RelAbund,
    GDGT3_RelAbund, GDGT4_RelAbund, Cren_RelAbund,
    GDGT5_RelAbund, GDGT6_RelAbund, GDGT7_RelAbund,
    GDGT8_RelAbund
  ) %>%
  pivot_longer(
    cols = -c(Site, Type_label, Incubation_days),
    names_to = "Compound",
    values_to = "RelAbund"
  ) %>%
  mutate(
    RelAbund = replace_na(RelAbund, 0),
    Compound = recode(
      Compound,
      GDGT0_RelAbund = "GDGT-0",
      GDGT1_RelAbund = "GDGT-1",
      GDGT2_RelAbund = "GDGT-2",
      GDGT3_RelAbund = "GDGT-3",
      GDGT4_RelAbund = "GDGT-4",
      Cren_RelAbund = "Cren-4",
      GDGT5_RelAbund = "GDGT-5",
      GDGT6_RelAbund = "GDGT-6",
      GDGT7_RelAbund = "GDGT-7",
      GDGT8_RelAbund = "GDGT-8"
    ),
    Compound = factor(
      Compound,
      levels = c(
        "GDGT-0", "GDGT-1", "GDGT-2", "GDGT-3", "GDGT-4",
        "Cren-4", "GDGT-5", "GDGT-6", "GDGT-7", "GDGT-8"
      )
    )
  ) %>%
  group_by(Site, Type_label, Incubation_days) %>%
  mutate(
    RelAbund_plot = RelAbund / sum(RelAbund, na.rm = TRUE)
  ) %>%
  ungroup()

#### BP data: faceted by Type, x-axis = Incubation_days ####
bp_long <- lipid_summary %>%
  select(
    Site, Type_label, Incubation_days,
    BP0_RelAbund, BP1_RelAbund, BP2_RelAbund,
    BP3_RelAbund, BP4_RelAbund
  ) %>%
  pivot_longer(
    cols = -c(Site, Type_label, Incubation_days),
    names_to = "Compound",
    values_to = "RelAbund"
  ) %>%
  mutate(
    RelAbund = replace_na(RelAbund, 0),
    Compound = recode(
      Compound,
      BP0_RelAbund = "BP-0",
      BP1_RelAbund = "BP-1",
      BP2_RelAbund = "BP-2",
      BP3_RelAbund = "BP-3",
      BP4_RelAbund = "BP-4"
    ),
    Compound = factor(
      Compound,
      levels = c("BP-0", "BP-1", "BP-2", "BP-3", "BP-4")
    )
  ) %>%
  group_by(Site, Type_label, Incubation_days) %>%
  mutate(
    RelAbund_plot = 100 * RelAbund / sum(RelAbund, na.rm = TRUE)
  ) %>%
  ungroup()

#### helper theme ####
panel_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.tag = element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = c(0.06, 0.95),
    legend.position = "none"
  )

#### plotting functions ####
# make_ri_plot <- function(site_name, panel_letter) {
#   ggplot(
#     ri_long %>% filter(Site == site_name),
#     aes(
#       x = Incubation_days,
#       y = Ring_Index,
#       color = RI_type,
#       shape = Type_label,
#       linetype = Type_label,
#       group = interaction(RI_type, Type_label)
#     )
#   ) +
#     geom_line(linewidth = 0.7, na.rm = TRUE) +
#     geom_point(size = 2.8, fill = "white", stroke = 1.1, na.rm = TRUE) +
#     scale_color_manual(
#       values = c("RI-GDGT" = "#EF3B2C", "RI-BP" = "#40BFFF"),
#       name = NULL
#     ) +
#     scale_shape_manual(
#       values = c("Exp." = 21, "Cnt-Live" = 24, "Cnt-Kill" = 25),
#       name = NULL
#     ) +
#     scale_linetype_manual(
#       values = c("Exp." = "solid", "Cnt-Live" = "dashed", "Cnt-Kill" = "dotted"),
#       name = NULL
#     ) +
#     scale_x_continuous(
#       breaks = c(0, 3, 14),
#       limits = c(-0.5, 14.5)
#     ) +
#     scale_y_continuous(
#       limits = c(0, 3.35),
#       breaks = c(0, 1, 2, 3)
#     ) +
#     labs(
#       title = site_name,
#       x = "Incubation Time (days)",
#       y = "Ring Index",
#       tag = panel_letter
#     ) +
#     panel_theme
# }

make_ri_plot <- function(site_name, panel_letter) {
  ggplot(
    ri_long %>%
      filter(Site == site_name) %>%
      mutate(
        RI_fill = ifelse(RI_type == "RI-GDGT", Site, "RI-BP")
      ),
    aes(
      x = Incubation_days,
      y = Ring_Index,
      color = Site,
      fill = RI_fill,
      shape = Type_label,
      linetype = Type_label,
      group = interaction(RI_type, Type_label)
    )
  ) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 2.8, stroke = 1.1, na.rm = TRUE) +
    scale_color_manual(
      values = c("Beryl" = "green3", "ETAT-3" = "navyblue"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("Beryl" = "green3", "ETAT-3" = "navyblue", "RI-BP" = "white"),
      name = NULL
    ) +
    scale_shape_manual(
      values = c("Exp." = 21, "Cnt-Live" = 24, "Cnt-Kill" = 25),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("Exp." = "solid", "Cnt-Live" = "dashed", "Cnt-Kill" = "dotted"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = c(0, 3, 14),
      limits = c(-0.5, 14.5)
    ) +
    scale_y_continuous(
      limits = c(0, 3.35),
      breaks = c(0, 1, 2, 3)
    ) +
    labs(
      title = site_name,
      x = "Incubation Time (days)",
      y = "Ring Index",
      tag = panel_letter
    ) +
    panel_theme
}


make_gdgt_plot <- function(site_name, panel_letter) {
  ggplot(
    gdgt_long %>% filter(Site == site_name),
    aes(
      x = factor(Incubation_days),
      y = RelAbund_plot,
      fill = Compound
    )
  ) +
    geom_col(
      color = "white",
      linewidth = 0.25,
      position = position_stack(reverse = TRUE),
      width = 0.8
    ) +
    facet_grid(~ Type_label, scales = "free_x", space = "free_x", drop = TRUE) +
    scale_fill_manual(
      values = setNames(GDGT_bar_colors, levels(gdgt_long$Compound)),
      breaks = rev(levels(gdgt_long$Compound)),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1)
    ) +
    labs(
      x = "Incubation Time (days)",
      y = "Relative Abundance",
      tag = panel_letter
    ) +
    panel_theme +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid = element_blank()
    )
}



make_bp_plot <- function(site_name, panel_letter) {
  ggplot(
    bp_long %>% filter(Site == site_name),
    aes(
      x = factor(Incubation_days),
      y = RelAbund_plot,
      fill = Compound
    )
  ) +
    geom_col(
      color = "white",
      linewidth = 0.25,
      position = position_stack(reverse = TRUE),
      width = 0.8
    ) +
    facet_grid(~ Type_label, scales = "free_x", space = "free_x", drop = TRUE) +
    scale_fill_manual(
      values = setNames(BP_bar_colors, levels(bp_long$Compound)),
      breaks = rev(levels(bp_long$Compound)),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 50, 100),
      labels = c("0", "0.5", "1.0")
    ) +
    labs(
      x = "Incubation Time (days)",
      y = "Relative Abundance",
      tag = panel_letter
    ) +
    panel_theme +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid = element_blank()
    )
}

#### legend plots ####

ri_legend_df <- tibble(
  legend_group = c(
    "Ring Index", "Ring Index",
    "Treatment", "Treatment", "Treatment"
  ),
  label = c(
    "RI-GDGT", "RI-BP",
    "Experimental", "Control-Live", "Control-Killed"
  ),
  x = 1,
  # y = c(5, 4, 2.5, 1.5, 0.5),
  y = c(5.2, 4.2, 2.2, 1.2, 0.2),
  shape = c(21, 21, 21, 24, 25),
  fill = c("black", "white", "white", "white", "white")
)

ri_leg_plot <- ggplot(ri_legend_df, aes(x = x, y = y)) +
  geom_point(
    aes(shape = shape, fill = fill),
    color = "black",
    size = 3,
    stroke = 1.1,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = label),
    x = 1.25,
    hjust = 0,
    size = 4
  ) +
  annotate(
    "text",
    x = 1,
    y = 6,
    label = "Ring Index",
    hjust = 0,
    size = 4.2
  ) +
  annotate(
    "text",
    x = 1,
    y = 3,
    label = "Treatment",
    hjust = 0,
    size = 4.2
  ) +
  scale_shape_identity() +
  scale_fill_identity() +
  coord_cartesian(
    xlim = c(0.8, 3.2),
    ylim = c(0, 6.2),
    clip = "off"
  ) +
  theme_void()


gdgt_legend_plot <- ggplot(
  gdgt_long,
  aes(
    x = factor(Incubation_days),
    y = RelAbund_plot,
    fill = Compound
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = setNames(GDGT_bar_colors, levels(gdgt_long$Compound)),
    breaks = rev(levels(gdgt_long$Compound)),
    name = NULL
  ) +
  theme_void() +
  theme(legend.position = "right")


bp_legend_plot <- ggplot(
  bp_long,
  aes(
    x = factor(Incubation_days),
    y = RelAbund_plot,
    fill = Compound
  )
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = setNames(BP_bar_colors, levels(bp_long$Compound)),
    breaks = rev(levels(bp_long$Compound)),
    name = NULL
  ) +
  theme_void() +
  theme(legend.position = "right")

#### extract legends ####
gdgt_leg <- cowplot::get_legend(gdgt_legend_plot)
bp_leg <- cowplot::get_legend(bp_legend_plot)

gdgt_leg_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(gdgt_leg, x = -0.3)
bp_leg_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(bp_leg, x = -0.33)

#### make individual panels ####
p_A <- make_ri_plot("Beryl", "A")
p_B <- make_ri_plot("ETAT-3", "B")

p_C <- make_gdgt_plot("Beryl", "C")
p_D <- make_gdgt_plot("ETAT-3", "D")

p_E <- make_bp_plot("Beryl", "E")
p_F <- make_bp_plot("ETAT-3", "F")

#### combine final figure ####
p_lipid_profiles <-
  (p_A | p_B | ri_leg_plot) /
  (p_C | p_D | gdgt_leg_plot) /
  (p_E | p_F | bp_leg_plot) +
  plot_layout(
    widths = c(1, 1, 0.2),
    heights = c(1, 1, 1)
  )

p_lipid_profiles

ggsave(
  filename = file.path(figure_dir, "FigureS2_lipid_profiles_by_timepoint.png"),
  plot = p_lipid_profiles,
  width = 11,
  height = 7,
  dpi = 300
)