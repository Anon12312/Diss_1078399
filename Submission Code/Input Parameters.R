#### SECTION 1 - SETUP AND LOAD DATA ####

# Clear the R environment by removing all objects
rm(list=ls())

# Load the dplyr library, which is used for data manipulation in this script
library(dplyr)

# Set the working directory to the specified path
setwd("/Users/harrystreet/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")

# Load the preprocessed data from a previous analysis
load("poverty_master.rda")

# Load the country-specific data from various CSV files

# Load the Gross National Income (GNI) per capita for 2021 and rename the relevant column
country_GNI <- read.csv("Country GNI Per Capita 2021 USD.csv") %>%
  select(Country.Code, X2021) %>%
  rename(GNI_per_capita = X2021)

# Load the out-of-pocket (OOP) health expenditure data for 2021, convert to a proportion, and rename the relevant column
country_OOP <- read.csv("Country OOP 2021.csv") %>%
  mutate(X2021 = X2021 / 100) %>%
  select(Country, X2021) %>%
  rename(oop = X2021)

# Load the life expectancy data and rename the relevant columns for further use
country_LE <- read.csv("Country Life Expectancy.csv") %>%
  select(SpatialDimValueCode, FactValueNumeric) %>%
  rename(ISO3 = SpatialDimValueCode, LE = FactValueNumeric)

# Load the hospital costs data, originally from the WHO 2008 CHOICE dataset, adjusted to 2021 USD
hosp_costs <- read.csv("hospital_costs.csv")

# Create a mapping of countries to their corresponding ISO3 codes
country_iso3_mapping <- data.frame(
  Country = c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi",
              "Cameroon", "Central African Republic", "Chad", "Comoros", "Congo",
              "Côte d'Ivoire", "Democratic Republic of the Congo", "Equatorial Guinea", "Eritrea", "Eswatini",
              "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea",
              "Guinea-Bissau", "Kenya", "Liberia", "Madagascar", "Malawi",
              "Mali", "Mauritania", "Mozambique", "Namibia", "Niger",
              "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal", "Sierra Leone",
              "South Africa", "South Sudan", "Togo", "Uganda", "United Republic of Tanzania",
              "Zambia", "Zimbabwe"),
  ISO3 = c("AGO", "BEN", "BWA", "BFA", "BDI", 
           "CMR", "CAF", "TCD", "COM", "COG", 
           "CIV", "COD", "GNQ", "ERI", "SWZ", 
           "ETH", "GAB", "GMB", "GHA", "GIN", 
           "GNB", "KEN", "LBR", "MDG", "MWI", 
           "MLI", "MRT", "MOZ", "NAM", "NER", 
           "NGA", "RWA", "STP", "SEN", "SLE", 
           "ZAF", "SSD", "TGO", "UGA", "TZA", 
           "ZMB", "ZWE")
)

# Update hospital costs data to ensure consistency in country names, filter to include only relevant countries, and merge with ISO3 codes
hosp_costs <- hosp_costs %>%
  mutate(Country = ifelse(Country == "Swaziland","Eswatini",Country)) %>%
  filter(Country %in% country_iso3_mapping$Country) %>%
  left_join(country_iso3_mapping, by = "Country")

# Calculate median values for columns of interest in the hospital costs data
median_values <- hosp_costs %>%
  summarise(across(c(Mean.value.from.sample, High.estimate..95..uncertainty.interval., Low.estimate..95..uncertainty.interval., SD), median, na.rm = TRUE))

# Assign median values to Zimbabwe where data is missing
hosp_costs <- hosp_costs %>%
  mutate(
    Mean.value.from.sample = ifelse(Country == "Zimbabwe" & is.na(Mean.value.from.sample), median_values$Mean.value.from.sample, Mean.value.from.sample),
    High.estimate..95..uncertainty.interval. = ifelse(Country == "Zimbabwe" & is.na(High.estimate..95..uncertainty.interval.), median_values$High.estimate..95..uncertainty.interval., High.estimate..95..uncertainty.interval.),
    Low.estimate..95..uncertainty.interval. = ifelse(Country == "Zimbabwe" & is.na(Low.estimate..95..uncertainty.interval.), median_values$Low.estimate..95..uncertainty.interval., Low.estimate..95..uncertainty.interval.),
    SD = ifelse(Country == "Zimbabwe" & is.na(SD), median_values$SD, SD)
  )

#### SECTION 2 - FUNCTION TO GENERATE MONTE CARLO PARAMETERS ####

# Set the number of simulations to be run in the Monte Carlo model
n_sims_montecarlo <- 2000

# Define a function to generate country-specific parameters for Monte Carlo simulations
generate_params_for_country <- function(country_name, country_iso3) {
  # Extract the GNI, OOP expenditure, and life expectancy for the specified country
  gni <- country_GNI %>% filter(Country.Code == country_iso3) %>% pull(GNI_per_capita)
  oop <- country_OOP %>% filter(Country == country_name) %>% pull(oop)
  le <- country_LE %>% filter(ISO3 == country_iso3) %>% pull(LE)
  
  # Calculate daily GNI based on the annual GNI
  cost_GNI_daily <- gni / 365
  
  # Extract the mean and standard deviation of hospital costs for the specified country
  hospital_cost_mean <- hosp_costs %>% filter(Country == country_name) %>% pull(Mean.value.from.sample)
  hospital_cost_sd <- hosp_costs %>% filter(Country == country_name) %>% pull(SD)
  
  # Generate Monte Carlo simulations for OOP expenditure
  oop_mean <- oop
  oop_sd <- (0.10 * oop_mean) / 1.96  # Calculate SD to achieve ±10% around the mean
  parvec_oop <- rnorm(n_sims_montecarlo, oop_mean, oop_sd)
  
  # Generate Monte Carlo simulations for clinical costs
  # Clinical cost estimates are based on a 2014 USD study and adjusted for 2021 inflation
  # https://gh.bmj.com/content/2/2/e000176
  diagnostic_test <- 1.22
  diagnostic_test_sd <- 0.47
  drug_cost <- 0.95
  drug_cost_sd <- 0.33
  hospital_cost <- hospital_cost_mean
  hospital_cost_sd <- hospital_cost_sd
  
  parvec_diagnostic_test <- rnorm(n_sims_montecarlo, diagnostic_test, diagnostic_test_sd)
  parvec_drug_cost <- rnorm(n_sims_montecarlo, drug_cost, drug_cost_sd)
  parvec_hospital_cost <- rnorm(n_sims_montecarlo, hospital_cost, hospital_cost_sd)
  
  # Generate Monte Carlo simulations for resistant cases
  # Clinical cost estimates are based on a 2014 USD study and adjusted for 2021 inflation
  # https://gh.bmj.com/content/2/2/e000176
  drug_cost_severe <- 10.26
  drug_cost_severe_sd <- 2.94
  hospital_cost_severe <- hospital_cost_mean
  hospital_cost_severe_sd <- hospital_cost_sd
  
  parvec_drug_cost_severe <- rnorm(n_sims_montecarlo, drug_cost_severe, drug_cost_severe_sd)
  parvec_hospital_cost_severe <- rnorm(n_sims_montecarlo, hospital_cost_severe, hospital_cost_severe_sd)
  
  # Generate Monte Carlo simulations for hospital duration parameters
  dur_ill_prehosp_mean <- 2
  dur_ill_prehosp_sd <- 0.5
  dur_hosp_clinical_survived_mean <- 3
  dur_hosp_clinical_survived_sd <- 1
  dur_hosp_clinical_survived_severe_mean <- 6
  dur_hosp_clinical_survived_severe_sd <- 1
  dur_hosp_resistant_survived_mean <- 6
  dur_hosp_resistant_survived_sd <- 1
  dur_hosp_died_mean <- 6
  dur_hosp_died_sd <- 1
  grief_mean <- 41.09
  grief_sd <- 42
  
  parvec_dur_ill_prehosp <- rnorm(n_sims_montecarlo, dur_ill_prehosp_mean, dur_ill_prehosp_sd)
  parvec_dur_hosp_clinical_survived <- rnorm(n_sims_montecarlo, dur_hosp_clinical_survived_mean, dur_hosp_clinical_survived_sd)
  parvec_dur_hosp_clinical_survived_severe <- rnorm(n_sims_montecarlo, dur_hosp_clinical_survived_severe_mean, dur_hosp_clinical_survived_severe_sd)
  parvec_dur_hosp_resistant_survived <- rnorm(n_sims_montecarlo, dur_hosp_resistant_survived_mean, dur_hosp_resistant_survived_sd)
  parvec_dur_hosp_died <- rnorm(n_sims_montecarlo, dur_hosp_died_mean, dur_hosp_died_sd)
  parvec_grief <- rnorm(n_sims_montecarlo, grief_mean, grief_sd)
  
  # Generate Monte Carlo simulations for disability weighting
  # Disability weights are based on data from the Global Burden of Disease (GBD) study 2019
  disutility_mod_mean <- 0.051
  disutility_mod_sd <- 0.0107
  disutility_severe_mean <- 0.133
  disutility_severe_sd <- 0.026
  
  parvec_disutility_mod <- rnorm(n_sims_montecarlo, disutility_mod_mean, disutility_mod_sd)
  parvec_disutility_severe <- rnorm(n_sims_montecarlo, disutility_severe_mean, disutility_severe_sd)
  
  # Generate Monte Carlo simulations for life expectancy parameters
  avg_years_left_to_live_mean <- le - 5
  avg_years_left_to_live_sd <- (0.10 * avg_years_left_to_live_mean) / 1.96  # Calculate SD to achieve ±10% around the mean
  parvec_avg_years <- rnorm(n_sims_montecarlo, avg_years_left_to_live_mean, avg_years_left_to_live_sd)
  
  
  parvec_avg_years <- rnorm(n_sims_montecarlo, avg_years_left_to_live_mean, avg_years_left_to_live_sd)
  
  # Combine all generated parameters into a data frame
  df_params_montecarlo <- data.frame(
    oop = parvec_oop,
    diagnostic_test = parvec_diagnostic_test,
    drug_cost = parvec_drug_cost,
    hospital_cost = parvec_hospital_cost,
    drug_cost_severe = parvec_drug_cost_severe,
    hospital_cost_severe = parvec_hospital_cost_severe,
    dur_ill_prehosp = parvec_dur_ill_prehosp,
    dur_hosp_clinical_survived = parvec_dur_hosp_clinical_survived,
    dur_hosp_clinical_survived_severe = parvec_dur_hosp_clinical_survived_severe,
    dur_hosp_resistant_survived = parvec_dur_hosp_resistant_survived,
    dur_hosp_died = parvec_dur_hosp_died,
    grief = parvec_grief,
    disutility_mod = parvec_disutility_mod,
    disutility_severe = parvec_disutility_severe,
    avg_years_left_to_live = parvec_avg_years
  )
  
  # Save the parameter data frame to a file, naming it based on the country ISO3 code
  save(df_params_montecarlo, file = paste0("params_montecarlo_", country_iso3, ".Rdata"))
}

#### SECTION 3 - LOOP THROUGH COUNTRIES AND SAVE PARAMETERS ####

# Loop through each country and its corresponding ISO3 code in the country_iso3_mapping data frame
for (i in 1:nrow(country_iso3_mapping)) {
  country_name <- country_iso3_mapping$Country[i]
  country_iso3 <- country_iso3_mapping$ISO3[i]
  
  # Generate and save Monte Carlo parameters for the current country
  generate_params_for_country(country_name, country_iso3)
}

