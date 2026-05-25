# Hotspring LH-SIP Analysis in R

R scripts, cleaned input data, summary tables, and figures associated with the hot spring LH-SIP analysis for the manuscript "Lipid Hydrogen Stable Isotope Probing Reveals Decadal-Scale Generation Times for Archaea in Hot Spring Sediments" (Harris et al., submitted) (bioRxiv: prepint doi: https://doi.org/10.64898/2026.05.15.725266)

## Repository structure

- `data_inputs/` — cleaned input datasets
- `scripts/` — R scripts used for analysis and figure generation
- `summary_tables/` — generated summary tables
- `figures/` — generated manuscript and supplement figures
- `HotSpring_LHSIP.Rproj` — RStudio project file


## Script descriptions

### `01_Merge cleaned datasets.R`
Loads & merges cleaned GDGT abundance, biphytane abundance, compound-specific δ²H, and experimental metadata tables, merges them into(`LHSIP_All_Merge.csv`)
Calculates:
- abundance-weighted mean lipid δ²H values
- fractional ²H abundance (F²H)
- ²H concentration in ppm
- propagated analytical uncertainties

### `02_Generate summary table.R`
Calculates apparent archaeal growth rates (μ), generation times (T_G), and lipid ²H enrichments for individual incubations and abundance-weighted averages using LH-SIP isotope incorporation equations
Generates Table 2 summary

### `03_Generate summary table by compound.R`
Calculates compound-specific and abundance-weighted:
- growth rates
- generation times
- biomass production rates
- new biomass production during incubations
- percent biomass increase
Generates Table 3 summary

### `04_Lipid_2H_Uptake_Figs.R`
Generates multipane; figures showing:
- lipid ²H enrichment over time for experimental replicates and controls
- lipid ²H enrichment over time for individual biphytane compounds
- generation times
Generates Figure S4 regression of biphytane ring number vs. generation time.

### `05_Lipid_Profiles_Summary_Figure.R`
Generates lipid distribution and ring index summary figures for IPL and CL fractions from Beryl and ETAT-3 sediments.

### `06_Generate_IPL_BP_summary_table.R`
Generates Supplementary Table S1 summarizing:
- biphytane relative abundances
- BP ring indices
- compound-specific lipid δ²H values
- abundance-weighted mean lipid δ²H
- Δδ²H/ring values
Includes propagated analytical uncertainty and replicate variability.

### `07_DetectionLimits_SensitivityAnalysis.R`
Calculates LH-SIP assay detection limits and performs sensitivity analyses exploring:
- minimum detectable lipid ²H enrichment
- detectable growth rates and generation times
- effects of tracer enrichment
- assimilation efficiency
- incubation duration

Generates Figure S5 and Figure S6.

### `08_Lipid_Profiles_by_Timepoint.R`
Generates Supplementary Figure S2 showing no changes in IPL lipid distributions over the incubation. 
Calculates and visualizes:
- iGDGT relative abundance distributions
- BP relative abundance distributions
- GDGT and BP ring indices (RI-GDGT and RI-BP)
- average profiles across biological replicates for each Site × Type × Incubation time point
