cat("\014")
rm(list=ls()) 
graphics.off()

#### load packages ####
library("readxl")
library("dplyr")
library("tidyverse")
library("purrr")

#### load cleaned datasets ####
GDGT <- read.csv("01 Cleaned Data/Cleaned_GDGT_Abundance_data.csv")
BP <- read.csv("01 Cleaned Data/Cleaned_BP_Abundance_data.csv")
d2H <- read.csv("01 Cleaned Data/Cleaned_BP_d2H_data.csv")
metadata <- read.csv("01 Cleaned Data/Cleaned_LHSIP_Experimental_Metadata.csv")

#### merge datasets ####
# remove index columns before merging datasets
BP <- BP %>%
  select(-Site, -X)
GDGT <- GDGT %>%
  select(-X)
d2H <- d2H %>%
  select(-X)
metadata <- metadata %>%
  select(-X)

# Merge GDGT, BP, and d2H by "ID_Fraction"
merged_data <- reduce(list(GDGT, BP, d2H), 
                      full_join, 
                      by = c("ID_Fraction", "ID", "Fraction"))

merged_data <- merged_data %>%
  select(-Site)
# Merge with metadata by "ID"
all_data <- merged_data %>%
  left_join(metadata, by = c("ID"))  


#### calculate abundance weighted mean d2H  ####
all_data <- all_data %>%
  mutate(
    # Total weight (sum of relative abundances)
    Total_BP_Abund_in_d2H = BP0_RelAbund + BP1_RelAbund + BP2_RelAbund,
    
    # Abundance-weighted mean d2H
    d2H_wt_mean = (BP0_d2H_mean * BP0_RelAbund +
                     BP1_d2H_mean * BP1_RelAbund +
                     BP2_d2H_mean * BP2_RelAbund) / Total_BP_Abund_in_d2H,
    
    # Propagated SEM for the weighted mean
    d2H_wt_mean_sem = sqrt(
      (BP0_RelAbund^2 * BP0_d2H_sem^2 +
         BP1_RelAbund^2 * BP1_d2H_sem^2 +
         BP2_RelAbund^2 * BP2_d2H_sem^2) / Total_BP_Abund_in_d2H^2
    )
  ) 
# NOTE: even though we could not get a d2H measurement for the small amount of 
# BP-3 detected in some samples, the weighted mean values represent >=95% of the BP pool
# summary(all_data$Total_BP_Abund_in_d2H)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   94.78   98.28   98.94   98.61   99.17  100.00       2 


#### convert delta 2H values to fractional 2H abundance and ppm 2H ####
# Define the VSMOW standard ratio
R_std <- 0.00015576  # D/H ratio of VSMOW

d2H_map <- tribble(
  ~mean_col,            ~sem_col,
  "BP0_d2H_mean",       "BP0_d2H_sem",
  "BP1_d2H_mean",       "BP1_d2H_sem",
  "BP2_d2H_mean",       "BP2_d2H_sem",
  "BP3_d2H_mean",       "BP3_d2H_sem",
  "BP4_d2H_mean",       "BP4_d2H_sem",
  "d2H_wt_mean",        "d2H_wt_mean_sem",
  "Tracer_Water_d2H",   "Tracer_Water_d2H_sem"
)

for (i in seq_len(nrow(d2H_map))) {
  col <- d2H_map$mean_col[i]
  sem_col <- d2H_map$sem_col[i]
  
  F2H_col <- str_replace(col, "d2H", "F2H")
  F2H_sem_col <- paste0(F2H_col, "_sem")
  D_ppm_col <- str_replace(F2H_col, "F2H", "D_ppm")
  D_ppm_sem_col <- paste0(D_ppm_col, "_sem")
  
  all_data <- all_data %>%
    mutate(
      !!F2H_col := {
        delta <- .data[[col]]
        R <- R_std * (delta / 1000 + 1)
        R / (1 + R)
      },
      !!F2H_sem_col := {
        delta <- .data[[col]]
        sem <- .data[[sem_col]]
        R <- R_std * (delta / 1000 + 1)
        dF_dDelta <- (R_std / 1000) / (1 + R)^2
        abs(dF_dDelta) * sem
      },
      !!D_ppm_col := .data[[F2H_col]] * 1e6,
      !!D_ppm_sem_col := .data[[F2H_sem_col]] * 1e6
    )
}


#### Save worked up data ####
write.csv(all_data, file = "01 Cleaned Data/LHSIP_All_Merge.csv")
