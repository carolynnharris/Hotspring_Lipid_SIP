cat("\014") #clears console
rm(list=ls()) 
graphics.off()

################################################################################
# 05_Lipid_Profiles_Summary_Figure.R
#
# This script:
# 1. Loads cleaned and formatted GDGT and BP abundance data
# 2. Calculates Ring Index values for GDGT and BP distributions
# 3. Creates lipid inventory & profile summary figure (Fig 3 in main text)
# Input:
# - data_inputs/LHSIP_Lipid_Summary_Plot.csv
#
# Output:
# - figures/Figure3_LipidInventorySummary.png
################################################################################



#### load relevant packages ####
library("tidyverse")

#### paths ####
clean_dir <- "data_inputs"
figure_dir <- "figures"


#### load data ####
dat <- read.csv("data_inputs/Lipid_Inventory_Plot.csv")

#### load formatted plotting data ####
lipid_plot <- dat %>%
  mutate(
    row_id = row_number(),
    Site = ifelse(is.na(Site) | Site == "", "Spacer", Site),
    Fraction = ifelse(is.na(Fraction) | Fraction == "", "Total", Fraction),
    row_type = case_when(
      rowSums(across(GDGT0:GDGT8)) > 0 ~ "iGDGT",
      rowSums(across(BP0:BP4)) > 0 ~ "BP",
      IPL > 0 | CL > 0 ~ "Total",
      TRUE ~ "Spacer"
    ),
    row_label = case_when(
      row_type == "iGDGT" ~ paste0("iGDGT-", Fraction),
      row_type == "BP" ~ paste0("BP-", Fraction),
      row_type == "Total" ~ "Total",
      TRUE ~ ""
    )
  )

#### symbology ####
GDGT_bar_colors <- c(
  "#00EEEE", "#1DD2E1", "#3BB7D4", "#599BC8", "#7780BB",
  "grey50", "#9464AE", "#B249A2", "#D02D95", "#EE1289"
)

BP_bar_colors <- c(
  "#00EEEE", "#599BC8", "#9464AE", "#B249A2", "#EE1289"
)

fraction_bar_colors <- c(
  IPL = "grey20",
  CL = "grey80"
)

site_cols <- c(
  Beryl = "green3",
  `ETAT-3` = "darkblue"
)


#### calculate Ring Index for GDGT and BP values ####
lipid_plot <- lipid_plot %>%
  mutate(
    RI_GDGT = (
      0 * GDGT0 +
        1 * GDGT1 +
        2 * GDGT2 +
        3 * GDGT3 +
        4 * (GDGT4 + Cren) +
        5 * GDGT5 +
        6 * GDGT6 +
        7 * GDGT7 +
        8 * GDGT8
    ) / (
      GDGT0 + GDGT1 + GDGT2 + GDGT3 +
        GDGT4 + Cren + GDGT5 + GDGT6 + GDGT7 + GDGT8
    ),
    
    RI_BP = (
      0 * BP0 +
        1 * BP1 +
        2 * BP2 +
        3 * BP3 +
        4 * BP4
    ) / (
      BP0 + BP1 + BP2 + BP3 + BP4
    )
  )

#### pull RI values for figure annotation ####
Beryl_IPL_RI_GDGT <- lipid_plot %>%
  filter(Site == "Beryl", Fraction == "IPL", row_type == "iGDGT") %>%
  pull(RI_GDGT) %>%
  round(2)

Beryl_CL_RI_GDGT <- lipid_plot %>%
  filter(Site == "Beryl", Fraction == "CL", row_type == "iGDGT") %>%
  pull(RI_GDGT) %>%
  round(2)

Beryl_IPL_RI_BP <- lipid_plot %>%
  filter(Site == "Beryl", Fraction == "IPL", row_type == "BP") %>%
  pull(RI_BP) %>%
  round(2)

Beryl_CL_RI_BP <- lipid_plot %>%
  filter(Site == "Beryl", Fraction == "CL", row_type == "BP") %>%
  pull(RI_BP) %>%
  round(2)

ETAT_IPL_RI_GDGT <- lipid_plot %>%
  filter(Site == "ETAT-3", Fraction == "IPL", row_type == "iGDGT") %>%
  pull(RI_GDGT) %>%
  round(2)

ETAT_CL_RI_GDGT <- lipid_plot %>%
  filter(Site == "ETAT-3", Fraction == "CL", row_type == "iGDGT") %>%
  pull(RI_GDGT) %>%
  round(2)

ETAT_IPL_RI_BP <- lipid_plot %>%
  filter(Site == "ETAT-3", Fraction == "IPL", row_type == "BP") %>%
  pull(RI_BP) %>%
  round(2)

ETAT_CL_RI_BP <- lipid_plot %>%
  filter(Site == "ETAT-3", Fraction == "CL", row_type == "BP") %>%
  pull(RI_BP) %>%
  round(2)


#### make matrix for base R stacked horizontal barplot ####
plot_cols <- c(
  "GDGT0", "GDGT1", "GDGT2", "GDGT3", "GDGT4", "Cren",
  "GDGT5", "GDGT6", "GDGT7", "GDGT8",
  "BP0", "BP1", "BP2", "BP3", "BP4",
  "IPL", "CL"
)

plot_matrix <- lipid_plot %>%
  select(all_of(plot_cols)) %>%
  as.matrix()

rownames(plot_matrix) <- lipid_plot$row_label

# Reverse row order for horizontal barplot display.
plot_matrix <- plot_matrix[nrow(plot_matrix):1, ]

bar_cols <- c(GDGT_bar_colors, BP_bar_colors, fraction_bar_colors)

#### save figure ####
png(
  filename = file.path(figure_dir, "Figure3_LipidInventorySummary.png"),
  width = 7.5,
  height = 4,
  units = "in",
  res = 300
)

par(
  mfrow = c(1, 1),
  mar = c(3, 9, 2, 9),
  oma = c(0, 0, 0, 0),
  mgp = c(1.6, 0.4, 0)
)

bar_positions <- barplot(
  t(plot_matrix),
  beside = FALSE,
  horiz = TRUE,
  col = bar_cols,
  border = "white",
  cex.axis = 0.8,
  cex.lab = 1,
  cex.names = 0.8,
  ylab = "",
  xlab = "Relative Abundance",
  las = 1,
  names.arg = rownames(plot_matrix)
)

#### site boxes ####
rect(-0.005, 7.3, 1.005, 13.3, border = site_cols["Beryl"], lwd = 3, xpd = TRUE)
rect(-0.005, 0.1, 1.005, 6.1, border = site_cols["ETAT-3"], lwd = 3, xpd = TRUE)

#### site labels and manuscript inventory circles ####
# Total iGDGT inventories are manuscript Table 1 values.
mtext("Beryl", col = site_cols["Beryl"], side = 3, line = -1.3,
      cex = 1, font = 2, adj = -0.43)

points(-0.35, 11, pch = 16, cex = 7, col = site_cols["Beryl"], xpd = TRUE)
text(-0.35, 11, labels = "80", col = "white", cex = 1, font = 2, xpd = TRUE)
text(-0.35, 9.2, labels = "µg iGDGT/", col = site_cols["Beryl"], cex = 0.7, font = 2, xpd = TRUE)
text(-0.35, 8.7, labels = "g dry sed", col = site_cols["Beryl"], cex = 0.7, font = 2, xpd = TRUE)

mtext("ETAT-3", col = site_cols["ETAT-3"], side = 3, line = -9.2,
      cex = 1, font = 2, adj = -0.46)

points(-0.35, 4.2, pch = 16, cex = 3, col = site_cols["ETAT-3"], xpd = TRUE)
text(-0.35, 4.2, labels = "12", col = "white", cex = 0.9, font = 2, xpd = TRUE)
text(-0.35, 3.2, labels = "µg iGDGT/", col = site_cols["ETAT-3"], cex = 0.7, font = 2, xpd = TRUE)
text(-0.35, 2.7, labels = "g dry sed", col = site_cols["ETAT-3"], cex = 0.7, font = 2, xpd = TRUE)

#### RI labels ####

RI_col <- "red"
RI_font <- 3

text(1.05, 13.7, labels = "RI", col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)

text(1.05, 12.7, labels = Beryl_IPL_RI_GDGT, col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 11.5, labels = Beryl_CL_RI_GDGT, col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 10.3, labels = Beryl_IPL_RI_BP,    col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 9.1,  labels = Beryl_CL_RI_BP,    col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)

text(1.05, 5.5, labels = ETAT_IPL_RI_GDGT, col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 4.3, labels = ETAT_CL_RI_GDGT, col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 3.1, labels = ETAT_IPL_RI_BP,    col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)
text(1.05, 2.0, labels = ETAT_CL_RI_BP,    col = RI_col, cex = 0.8, font = RI_font, xpd = TRUE)

#### legends ####
legend(
  "topright",
  inset = c(-0.35, 0),
  title = "iGDGT",
  legend = c("0", "1", "2", "3", "4", "Cren", "5", "6", "7", "8"),
  fill = GDGT_bar_colors,
  bty = "n",
  cex = 0.9,
  xpd = TRUE
)

legend(
  "topright",
  inset = c(-0.45, 0),
  title = "BP",
  legend = c("0", "1", "2", "3", "4"),
  fill = BP_bar_colors,
  bty = "n",
  cex = 0.9,
  xpd = TRUE
)

legend(
  "topright",
  inset = c(-0.32, 0.75),
  title = "Total",
  legend = c("IPL", "CL"),
  fill = fraction_bar_colors,
  bty = "n",
  cex = 0.9,
  xpd = TRUE
)

dev.off()


