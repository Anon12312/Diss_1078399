# Clear the current workspace to ensure no residual objects affect the analysis
rm(list=ls())

#### SECTION 1 - SETUP ####

#### SECTION 1.1 - IMPORTING THE DATA ####

# Load the necessary libraries for data manipulation, visualisation, parallel processing, and geospatial analysis
library(dplyr)
library(ggplot2)
library(parallel)
library(scales)
library(pbmcapply)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggpubr)
library(tidyr)

# Set the working directory to the location where the data files are stored
setwd("/Users/X/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")

# Load pre-processed datasets that will be used for the analysis
load("all_malaria.rda") # Malaria case data for base and vaccine efficacy (VE) scenarios
load("daly_master.rda") # Monetised Disability-Adjusted Life Years (DALYs) data
load("poverty_master.rda") # Data on catastrophic health expenditure and impoverishment
load("vsl_master.rda") # Value of Statistical Life (VSL) estimates

# Load population data from a CSV file
population_data <- read.csv("Pop1_byCountry.csv") # Population data by country

# Define the number of Monte Carlo simulations to be conducted
n_sims_montecarlo <- 2000

# Create a mapping between country names and their corresponding ISO3 codes, which will be utilised in subsequent data operations
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

# Load data on Gross National Income (GNI) per capita for each country
country_GNI <- read.csv("Country GNI Per Capita 2021 USD.csv") %>%
  select(Country.Code, X2021) %>%
  rename(GNI_per_capita = X2021)

# Load data on the economic loss associated with the death of a child, adjusted for non-health GDP loss
gdp_loss <- read.csv("Child NonHealth GDP Loss 2021 USD.csv") %>%
  select(Country, Loss.Per.Child.Adjusted) %>%
  rename(GDP_Loss_Per_Child = Loss.Per.Child.Adjusted) %>%
  left_join(., country_iso3_mapping, by = "Country") %>%
  select(Country,ISO3,GDP_Loss_Per_Child)

# Define the discount rate to be applied in economic evaluations, and create a function for applying this discount
discRate <- 0.03
discount <- function(value, rate, year_k, base_year = 2021) {
  return(value * (1 / (1 + rate) ^ (year_k - base_year)))
}

#### SECTION 1.2 - FILL DATA GAPS ####

# Identify and fill in missing values in the Value of Statistical Life (VSL) data using appropriate estimates
na_vsl <- vsl_master %>% filter(is.na(vsl_estimate_2021_IntDoll))
na_vsl <- na_vsl %>% filter(GID_0 %in% country_iso3_mapping$ISO3)

# Impute missing VSL values based on available data, with specific estimates for certain countries
vsl_master <- vsl_master %>%
  mutate(vsl_estimate_2021_IntDoll = ifelse(GID_0 == "COD", 22729.16417, vsl_estimate_2021_IntDoll),
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "COG", 302547.4287, vsl_estimate_2021_IntDoll),
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "GMB", 64287.78445, vsl_estimate_2021_IntDoll),
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "SWZ", 966510.1256, vsl_estimate_2021_IntDoll),
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "ERI", 56096.62603, vsl_estimate_2021_IntDoll), # Using 2021 international USD low-income median
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "SSD", 56096.62603, vsl_estimate_2021_IntDoll), # Using 2021 international USD low-income median
         vsl_estimate_2021_IntDoll = ifelse(GID_0 == "STP", 56096.62603, vsl_estimate_2021_IntDoll)) # Using 2021 international USD low-income median

# Similarly, identify and impute missing values in the GNI data
na_gni <- country_GNI %>%
  filter(is.na(GNI_per_capita))
na_gni <- na_gni %>%
  filter(Country.Code %in% country_iso3_mapping$ISO3)

# Impute missing GNI values with appropriate estimates for certain countries
country_GNI <- country_GNI %>%
  mutate(GNI_per_capita = ifelse(Country.Code == "ERI", 640, GNI_per_capita),
         GNI_per_capita = ifelse(Country.Code == "SSD", 410, GNI_per_capita))


#### SECTION 2 - CREATING THE FUNCTIONS ####

#### SECTION 2.1 - BASE SCENARIO ####

calculate_base_scenario <- function(data, results, params, VSL, cost_daly, cost_GNI_daily, gdp_loss_per_child, prob_catastrophic, prob_impoverishment) {
  
  # Calculate the proportion of severe cases
  prop_severe <- 0.138
  
  # Calculate costs for mild and severe clinical and resistant cases
  clinical_cost_mild <- params$diagnostic_test + params$drug_cost + params$hospital_cost * params$dur_hosp_clinical_survived
  clinical_cost_severe <- params$diagnostic_test + params$drug_cost_severe + params$hospital_cost_severe * params$dur_hosp_clinical_survived_severe
  
  resistant_cost_mild <- params$diagnostic_test + params$drug_cost + params$hospital_cost * params$dur_hosp_resistant_survived
  resistant_cost_severe <- params$diagnostic_test + params$drug_cost_severe + params$hospital_cost_severe * params$dur_hosp_resistant_survived
  
  # Iterate over each year in the dataset to calculate various health and economic outcomes
  for (i in seq_along(data$year)) {
    year_k <- data$year[i]
    
    # Split clinical and resistant cases into mild and severe
    clinical_cases_mild <- data$cases_VE0[i] * (1 - prop_severe)
    clinical_cases_severe <- data$cases_VE0[i] * prop_severe
    
    resistant_cases_mild <- data$resistant_cases_VE0[i] * (1 - prop_severe)
    resistant_cases_severe <- data$resistant_cases_VE0[i] * prop_severe
    
    # Debugging checks
    if (length(resistant_cases_mild) == 0 || length(resistant_cases_severe) == 0) {
      stop(paste("Error: resistant_cases_mild or resistant_cases_severe has zero length for year", year_k))
    }
    
    if (length(resistant_cost_mild) == 0 || length(resistant_cost_severe) == 0) {
      stop(paste("Error: resistant_cost_mild or resistant_cost_severe has zero length for year", year_k))
    }
    
    # Calculate Years of Life Lost (YLL) and apply discounting
    YLL <- data$deaths_VE0[i] * params$avg_years_left_to_live
    YLL_disc <- discount(YLL, discRate, year_k)
    
    # Calculate Disability-Adjusted Life Years (DALYs)
    DALY_clinical_mild <- clinical_cases_mild * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived) * params$disutility_mod / 365
    DALY_clinical_severe <- clinical_cases_severe * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived_severe) * params$disutility_severe / 365
    DALY_clinical <- DALY_clinical_mild + DALY_clinical_severe
    
    DALY_resistant_mild <- resistant_cases_mild * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_mod / 365
    DALY_resistant_severe <- resistant_cases_severe * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_severe / 365
    DALY_resistant <- DALY_resistant_mild + DALY_resistant_severe
    
    DALY_death <- YLL
    DALY_death_disc <- YLL_disc
    DALY_total <- DALY_clinical + DALY_resistant + DALY_death
    DALY_total_disc <- DALY_clinical + DALY_resistant + DALY_death_disc
    
    # Calculate the costs associated with the loss of statistical life, discounted DALYs, clinical care, and productivity loss
    Cost_VSL <- data$deaths_VE0[i] * VSL
    Cost_VSLY <- DALY_death_disc * VSL / params$avg_years_left_to_live
    
    Cost_DALY_total_disc <- DALY_total_disc * cost_daly
    Cost_DALY_total_disc <- discount(Cost_DALY_total_disc, discRate, year_k)
    
    Cost_clin_gvt <- (clinical_cases_mild * clinical_cost_mild + clinical_cases_severe * clinical_cost_severe) * (1 - params$oop)
    Cost_clin_oop <- (clinical_cases_mild * clinical_cost_mild + clinical_cases_severe * clinical_cost_severe) * params$oop
    
    Cost_hosp_resistant_gvt <- (resistant_cases_mild * resistant_cost_mild + resistant_cases_severe * resistant_cost_severe) * (1 - params$oop)
    Cost_hosp_resistant_oop <- (resistant_cases_mild * resistant_cost_mild + resistant_cases_severe * resistant_cost_severe) * params$oop
    
    # Debugging checks for costs
    if (length(Cost_hosp_resistant_gvt) == 0 || length(Cost_hosp_resistant_oop) == 0) {
      stop(paste("Error: Cost_hosp_resistant_gvt or Cost_hosp_resistant_oop has zero length for year", year_k))
    }
    
    Cost_prod_hospital <- clinical_cases_mild * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived) +
      clinical_cases_severe * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived_severe)
    
    Cost_prod_resistant_hospital <- resistant_cases_mild * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) +
      resistant_cases_severe * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived)
    
    Cost_prod_death <- data$deaths_VE0[i] * (params$grief * cost_GNI_daily + gdp_loss_per_child)
    Cost_prod_total <- Cost_prod_hospital + Cost_prod_resistant_hospital + Cost_prod_death
    
    Cost_care_gvt <- Cost_clin_gvt + Cost_hosp_resistant_gvt
    Cost_care_oop <- Cost_clin_oop + Cost_hosp_resistant_oop
    
    # Apply discounting to the calculated costs
    Cost_clin_gvt_disc <- discount(Cost_clin_gvt, discRate, year_k)
    Cost_clin_oop_disc <- discount(Cost_clin_oop, discRate, year_k)
    Cost_hosp_resistant_gvt_disc <- discount(Cost_hosp_resistant_gvt, discRate, year_k)
    Cost_hosp_resistant_oop_disc <- discount(Cost_hosp_resistant_oop, discRate, year_k)
    Cost_prod_total_disc <- discount(Cost_prod_total, discRate, year_k)
    Cost_care_gvt_disc <- discount(Cost_care_gvt, discRate, year_k)
    Cost_care_oop_disc <- discount(Cost_care_oop, discRate, year_k)
    
    Cost_societal <- Cost_care_gvt + Cost_care_oop + Cost_prod_total
    Cost_societal_disc <- Cost_care_gvt_disc + Cost_care_oop_disc + Cost_prod_total_disc
    
    # Estimate the number of catastrophic and impoverishing cases due to out-of-pocket expenses
    Catastrophic_cases <- (data$cases_VE0[i] + data$resistant_cases_VE0[i]) * prob_catastrophic / 100
    Impoverished_cases <- (data$cases_VE0[i] + data$resistant_cases_VE0[i]) * prob_impoverishment / 100
    
    # Store the calculated results in the results dataframe
    results$YLL[i] <- YLL
    results$YLL_disc[i] <- YLL_disc
    results$DALY_clinical[i] <- DALY_clinical
    results$DALY_resistant[i] <- DALY_resistant
    results$DALY_death[i] <- DALY_death
    results$DALY_death_disc[i] <- DALY_death_disc
    results$DALY_total[i] <- DALY_total
    results$DALY_total_disc[i] <- DALY_total_disc
    results$Cost_VSL[i] <- Cost_VSL
    results$Cost_VSLY[i] <- Cost_VSLY
    results$Cost_DALY_total_disc[i] <- Cost_DALY_total_disc
    results$Cost_clin_gvt[i] <- Cost_clin_gvt
    results$Cost_clin_gvt_disc[i] <- Cost_clin_gvt_disc
    results$Cost_clin_oop[i] <- Cost_clin_oop
    results$Cost_clin_oop_disc[i] <- Cost_clin_oop_disc
    results$Cost_hosp_resistant_gvt[i] <- Cost_hosp_resistant_gvt
    results$Cost_hosp_resistant_gvt_disc[i] <- Cost_hosp_resistant_gvt_disc
    results$Cost_hosp_resistant_oop[i] <- Cost_hosp_resistant_oop
    results$Cost_hosp_resistant_oop_disc[i] <- Cost_hosp_resistant_oop_disc
    results$Cost_prod_hospital[i] <- Cost_prod_hospital
    results$Cost_prod_resistant_hospital[i] <- Cost_prod_resistant_hospital
    results$Cost_prod_death[i] <- Cost_prod_death
    results$Cost_prod_total[i] <- Cost_prod_total
    results$Cost_prod_total_disc[i] <- Cost_prod_total_disc
    results$Cost_care_gvt[i] <- Cost_care_gvt
    results$Cost_care_oop[i] <- Cost_care_oop
    results$Cost_societal[i] <- Cost_societal
    results$Cost_societal_disc[i] <- Cost_societal_disc
    results$Catastrophic_cases[i] <- Catastrophic_cases
    results$Impoverished_cases[i] <- Impoverished_cases
  }
  
  return(results)
}


#### SECTION 2.2 - VACCINATION SCENARIO FUNCTIONS ####

calculate_vaccination_scenario <- function(data, scenario, results, params, VSL, cost_daly, cost_GNI_daily, gdp_loss_per_child, prob_catastrophic, prob_impoverishment) {
  
  # Calculate the proportion of severe cases
  prop_severe <- 0.138
  
  # Calculate costs for mild and severe clinical and resistant cases
  clinical_cost_mild <- params$diagnostic_test + params$drug_cost + params$hospital_cost * params$dur_hosp_clinical_survived
  clinical_cost_severe <- params$diagnostic_test + params$drug_cost_severe + params$hospital_cost_severe * params$dur_hosp_clinical_survived_severe
  
  resistant_cost_mild <- params$diagnostic_test + params$drug_cost + params$hospital_cost * params$dur_hosp_resistant_survived
  resistant_cost_severe <- params$diagnostic_test + params$drug_cost_severe + params$hospital_cost_severe * params$dur_hosp_resistant_survived
  
  # Iterate over each year in the dataset to calculate various health and economic outcomes under the vaccination scenario
  for (i in seq_along(data$year)) {
    year_k <- data$year[i]
    
    # Split cases averted into mild and severe categories
    clinical_cases_averted_mild <- data[[paste0("cases_averted_", scenario)]][i] * (1 - prop_severe)
    clinical_cases_averted_severe <- data[[paste0("cases_averted_", scenario)]][i] * prop_severe
    
    resistant_cases_averted_mild <- data[[paste0("resistant_cases_averted_", scenario)]][i] * (1 - prop_severe)
    resistant_cases_averted_severe <- data[[paste0("resistant_cases_averted_", scenario)]][i] * prop_severe
    
    # Debugging checks
    if (length(resistant_cases_averted_mild) == 0 || length(resistant_cases_averted_severe) == 0) {
      stop(paste("Error: resistant_cases_averted_mild or resistant_cases_averted_severe has zero length for year", year_k, "scenario", scenario))
    }
    
    if (length(resistant_cost_mild) == 0 || length(resistant_cost_severe) == 0) {
      stop(paste("Error: resistant_cost_mild or resistant_cost_severe has zero length for year", year_k, "scenario", scenario))
    }
    
    # Calculate the Years of Life Lost (YLL) averted by the vaccination, applying discounting
    YLL <- data[[paste0("deaths_averted_", scenario)]][i] * params$avg_years_left_to_live
    YLL_disc <- discount(YLL, discRate, year_k)
    
    # Calculate the Disability-Adjusted Life Years (DALYs) averted by the vaccination
    DALY_clinical_mild <- clinical_cases_averted_mild * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived) * params$disutility_mod / 365
    DALY_clinical_severe <- clinical_cases_averted_severe * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived_severe) * params$disutility_severe / 365
    DALY_clinical <- DALY_clinical_mild + DALY_clinical_severe
    
    DALY_resistant_mild <- resistant_cases_averted_mild * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_mod / 365
    DALY_resistant_severe <- resistant_cases_averted_severe * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_severe / 365
    DALY_resistant <- DALY_resistant_mild + DALY_resistant_severe
    
    DALY_death <- YLL
    DALY_death_disc <- YLL_disc
    DALY_total <- DALY_clinical + DALY_resistant + DALY_death
    DALY_total_disc <- DALY_clinical + DALY_resistant + DALY_death_disc
    
    # Calculate the costs averted by the vaccination, including VSL, discounted DALYs, clinical care, and productivity loss
    Cost_VSL <- data[[paste0("deaths_averted_", scenario)]][i] * VSL
    Cost_VSLY <- DALY_death_disc * VSL / params$avg_years_left_to_live
    
    Cost_DALY_total_disc <- DALY_total_disc * cost_daly
    Cost_DALY_total_disc <- discount(Cost_DALY_total_disc, discRate, year_k)
    
    Cost_clin_gvt <- (clinical_cases_averted_mild * clinical_cost_mild + clinical_cases_averted_severe * clinical_cost_severe) * (1 - params$oop)
    Cost_clin_oop <- (clinical_cases_averted_mild * clinical_cost_mild + clinical_cases_averted_severe * clinical_cost_severe) * params$oop
    
    Cost_hosp_resistant_gvt <- (resistant_cases_averted_mild * resistant_cost_mild + resistant_cases_averted_severe * resistant_cost_severe) * (1 - params$oop)
    Cost_hosp_resistant_oop <- (resistant_cases_averted_mild * resistant_cost_mild + resistant_cases_averted_severe * resistant_cost_severe) * params$oop
    
    # Debugging checks for costs
    if (length(Cost_hosp_resistant_gvt) == 0 || length(Cost_hosp_resistant_oop) == 0) {
      stop(paste("Error: Cost_hosp_resistant_gvt or Cost_hosp_resistant_oop has zero length for year", year_k, "scenario", scenario))
    }
    
    Cost_prod_hospital <- clinical_cases_averted_mild * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived) +
      clinical_cases_averted_severe * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_clinical_survived_severe)
    
    Cost_prod_resistant_hospital <- resistant_cases_averted_mild * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) +
      resistant_cases_averted_severe * cost_GNI_daily * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived)
    
    Cost_prod_death <- data[[paste0("deaths_averted_", scenario)]][i] * (params$grief * cost_GNI_daily + gdp_loss_per_child)
    Cost_prod_total <- Cost_prod_hospital + Cost_prod_resistant_hospital + Cost_prod_death
    
    Cost_care_gvt <- Cost_clin_gvt + Cost_hosp_resistant_gvt
    Cost_care_oop <- Cost_clin_oop + Cost_hosp_resistant_oop
    
    # Apply discounting to the calculated costs under the vaccination scenario
    Cost_clin_gvt_disc <- discount(Cost_clin_gvt, discRate, year_k)
    Cost_clin_oop_disc <- discount(Cost_clin_oop, discRate, year_k)
    Cost_hosp_resistant_gvt_disc <- discount(Cost_hosp_resistant_gvt, discRate, year_k)
    Cost_hosp_resistant_oop_disc <- discount(Cost_hosp_resistant_oop, discRate, year_k)
    Cost_prod_total_disc <- discount(Cost_prod_total, discRate, year_k)
    Cost_care_gvt_disc <- discount(Cost_care_gvt, discRate, year_k)
    Cost_care_oop_disc <- discount(Cost_care_oop, discRate, year_k)
    
    Cost_societal <- Cost_care_gvt + Cost_care_oop + Cost_prod_total
    Cost_societal_disc <- Cost_care_gvt_disc + Cost_care_oop_disc + Cost_prod_total_disc
    
    # Estimate the number of catastrophic and impoverishing cases averted by the vaccination due to out-of-pocket expenses
    Catastrophic_cases <- (clinical_cases_averted_mild + clinical_cases_averted_severe + resistant_cases_averted_mild + resistant_cases_averted_severe) * prob_catastrophic / 100
    Impoverished_cases <- (clinical_cases_averted_mild + clinical_cases_averted_severe + resistant_cases_averted_mild + resistant_cases_averted_severe) * prob_impoverishment / 100
    
    # Store the calculated results in the results dataframe for the vaccination scenario
    results$YLL[i] <- YLL
    results$YLL_disc[i] <- YLL_disc
    results$DALY_clinical[i] <- DALY_clinical
    results$DALY_resistant[i] <- DALY_resistant
    results$DALY_death[i] <- DALY_death
    results$DALY_death_disc[i] <- DALY_death_disc
    results$DALY_total[i] <- DALY_total
    results$DALY_total_disc[i] <- DALY_total_disc
    results$Cost_VSL[i] <- Cost_VSL
    results$Cost_VSLY[i] <- Cost_VSLY
    results$Cost_DALY_total_disc[i] <- Cost_DALY_total_disc
    results$Cost_clin_gvt[i] <- Cost_clin_gvt
    results$Cost_clin_gvt_disc[i] <- Cost_clin_gvt_disc
    results$Cost_clin_oop[i] <- Cost_clin_oop
    results$Cost_clin_oop_disc[i] <- Cost_clin_oop_disc
    results$Cost_hosp_resistant_gvt[i] <- Cost_hosp_resistant_gvt
    results$Cost_hosp_resistant_gvt_disc[i] <- Cost_hosp_resistant_gvt_disc
    results$Cost_hosp_resistant_oop[i] <- Cost_hosp_resistant_oop
    results$Cost_hosp_resistant_oop_disc[i] <- Cost_hosp_resistant_oop_disc
    results$Cost_prod_hospital[i] <- Cost_prod_hospital
    results$Cost_prod_resistant_hospital[i] <- Cost_prod_resistant_hospital
    results$Cost_prod_death[i] <- Cost_prod_death
    results$Cost_prod_total[i] <- Cost_prod_total
    results$Cost_prod_total_disc[i] <- Cost_prod_total_disc
    results$Cost_societal[i] <- Cost_societal
    results$Cost_societal_disc[i] <- Cost_societal_disc
    results$Catastrophic_cases[i] <- Catastrophic_cases
    results$Impoverished_cases[i] <- Impoverished_cases
  }
  return(results)
}


#### SECTION 2.3 - COUNTRY FUNCTION ####

# This function processes data for a specific country, running simulations and calculating various outcomes
process_country <- function(country) {
  # Filter the extracted data to include only the selected country
  country_data <- extracted_data %>% filter(country == !!country)
  
  # Extract the Value of Statistical Life (VSL) for the country
  VSL <- vsl_master %>% filter(GID_0 == country_data$ISO3[1]) %>% select(vsl_estimate_2021_IntDoll) %>% pull()
  
  # Extract the GNI data for the country and calculate daily GNI
  gni <- country_GNI %>% filter(Country.Code == country_data$ISO3[1]) %>% pull(GNI_per_capita)
  cost_GNI_daily <- gni / 365
  
  # Extract the monetised DALY data for the country
  cost_daly <- daly_master %>%
    filter(GID_0 == country_data$ISO3[1]) %>%
    select(daly_value) %>%
    pull()
  
  # Extract the GDP loss per child data for the country
  gdp_loss_per_child <- gdp_loss %>%
    filter(ISO3 == country_data$ISO3[1]) %>%
    pull(GDP_Loss_Per_Child)
  
  # Extract the probability of catastrophic and impoverishing cases
  prob_catastrophic <- poverty_master %>%
    filter(GID_0 == country_data$ISO3[1]) %>%
    pull(prop_at_risk_catastrp_HCexp)
  
  prob_impoverishment <- poverty_master %>%
    filter(GID_0 == country_data$ISO3[1]) %>%
    pull(prop_at_risk_impoverishing)
  
  # Load the specific parameter file for the country to be used in Monte Carlo simulations
  load(paste0("params_montecarlo_", country_data$ISO3[1], ".Rdata"))
  
  # Initialize empty dataframes to store results
  results_base_all <- data.frame()
  all_results <- data.frame()
  
  # Run simulations for the base scenario
  for (sim in 1:n_sims_montecarlo) {
    # Extract parameters for the current simulation
    params <- df_params_montecarlo[sim, ]
    
    # Initialize a results dataframe for the current country and simulation
    results_base_sim <- data.frame(
      year = country_data$year,
      YLL = numeric(length(country_data$year)),
      YLL_disc = numeric(length(country_data$year)),
      DALY_clinical = numeric(length(country_data$year)),
      DALY_resistant = numeric(length(country_data$year)),
      DALY_death = numeric(length(country_data$year)),
      DALY_death_disc = numeric(length(country_data$year)),
      DALY_total = numeric(length(country_data$year)),
      DALY_total_disc = numeric(length(country_data$year)),
      Cost_VSL = numeric(length(country_data$year)),
      Cost_VSLY = numeric(length(country_data$year)),
      Cost_DALY_total_disc = numeric(length(country_data$year)),
      Cost_clin_gvt = numeric(length(country_data$year)),
      Cost_clin_gvt_disc = numeric(length(country_data$year)),
      Cost_clin_oop = numeric(length(country_data$year)),
      Cost_clin_oop_disc = numeric(length(country_data$year)),
      Cost_hosp_resistant_gvt = numeric(length(country_data$year)),
      Cost_hosp_resistant_gvt_disc = numeric(length(country_data$year)),
      Cost_hosp_resistant_oop = numeric(length(country_data$year)),
      Cost_hosp_resistant_oop_disc = numeric(length(country_data$year)),
      Cost_prod_hospital = numeric(length(country_data$year)),
      Cost_prod_resistant_hospital = numeric(length(country_data$year)),
      Cost_prod_death = numeric(length(country_data$year)),
      Cost_prod_total = numeric(length(country_data$year)),
      Cost_prod_total_disc = numeric(length(country_data$year)),
      Cost_societal = numeric(length(country_data$year)),
      Cost_societal_disc = numeric(length(country_data$year)),
      Catastrophic_cases = numeric(length(country_data$year)),
      Impoverished_cases = numeric(length(country_data$year)),
      Cases_averted = numeric(length(country_data$year)),  # New column to store cases averted
      Resistant_cases_averted = numeric(length(country_data$year)),  # New column to store resistant cases averted
      Deaths_averted = numeric(length(country_data$year))  # New column to store deaths averted
    )
    
    # Calculate the base scenario results using the function defined earlier
    scenario_results <- calculate_base_scenario(country_data, results_base_sim, params, VSL, cost_daly, cost_GNI_daily, gdp_loss_per_child, prob_catastrophic, prob_impoverishment)
    
    # Store the results with country information
    results_base_all <- rbind(results_base_all, cbind(country = country, ISO3 = country_data$ISO3[1], scenario_results))
  }
  
  # Run simulations for each vaccination efficacy scenario (VE1, VE2, VE3)
  for (scenario in c("VE1", "VE2", "VE3")) {
    for (sim in 1:n_sims_montecarlo) {
      # Extract parameters for the current simulation
      params <- df_params_montecarlo[sim, ]
      
      # Initialize a results dataframe for the current country and simulation
      results_sim <- data.frame(
        year = country_data$year,
        YLL = numeric(length(country_data$year)),
        YLL_disc = numeric(length(country_data$year)),
        DALY_clinical = numeric(length(country_data$year)),
        DALY_resistant = numeric(length(country_data$year)),
        DALY_death = numeric(length(country_data$year)),
        DALY_death_disc = numeric(length(country_data$year)),
        DALY_total = numeric(length(country_data$year)),
        DALY_total_disc = numeric(length(country_data$year)),
        Cost_VSL = numeric(length(country_data$year)),
        Cost_VSLY = numeric(length(country_data$year)),
        Cost_DALY_total_disc = numeric(length(country_data$year)),
        Cost_clin_gvt = numeric(length(country_data$year)),
        Cost_clin_gvt_disc = numeric(length(country_data$year)),
        Cost_clin_oop = numeric(length(country_data$year)),
        Cost_clin_oop_disc = numeric(length(country_data$year)),
        Cost_hosp_resistant_gvt = numeric(length(country_data$year)),
        Cost_hosp_resistant_gvt_disc = numeric(length(country_data$year)),
        Cost_hosp_resistant_oop = numeric(length(country_data$year)),
        Cost_hosp_resistant_oop_disc = numeric(length(country_data$year)),
        Cost_prod_hospital = numeric(length(country_data$year)),
        Cost_prod_resistant_hospital = numeric(length(country_data$year)),
        Cost_prod_death = numeric(length(country_data$year)),
        Cost_prod_total = numeric(length(country_data$year)),
        Cost_prod_total_disc = numeric(length(country_data$year)),
        Cost_societal = numeric(length(country_data$year)),
        Cost_societal_disc = numeric(length(country_data$year)),
        Catastrophic_cases = numeric(length(country_data$year)),
        Impoverished_cases = numeric(length(country_data$year)),
        Cases_averted = numeric(length(country_data$year)),  # New column for cases averted
        Resistant_cases_averted = numeric(length(country_data$year)),  # New column for resistant cases averted
        Deaths_averted = numeric(length(country_data$year))  # New column for deaths averted
      )
      
      # Calculate results for the current simulation and vaccination scenario
      scenario_results <- calculate_vaccination_scenario(country_data, scenario, results_sim, params, VSL, cost_daly, cost_GNI_daily, gdp_loss_per_child, prob_catastrophic, prob_impoverishment)
      
      # Store the results with country and scenario information
      all_results <- rbind(all_results, cbind(country = country, ISO3 = country_data$ISO3[1], scenario = scenario, scenario_results))
    }
  }
  
  return(list(results_base_all = results_base_all, all_results = all_results))
}

#### SECTION 3 - RESULTS ####

#### SECTION 3.1 - SIMULATE RESULTS ####

# Get a list of all countries from the extracted data
countries <- unique(extracted_data$country)

# Detect the number of available CPU cores for parallel processing
num_cores <- detectCores()

# Run the country-specific processing function across all countries in parallel
results <- pbmclapply(countries, process_country, mc.cores = num_cores)

# Combine the results from all countries into a single dataframe
results_base_all <- do.call(rbind, lapply(results, function(res) res$results_base_all))
all_results <- do.call(rbind, lapply(results, function(res) res$all_results))

# Define the vaccination scenarios and time horizons for summary statistics
scenarios <- c("VE1", "VE2", "VE3")
time_horizons <- c(3, 5, 10)

# Function to calculate cumulative statistics (e.g., mean DALYs, costs) for a given time horizon
calculate_cumulative_stats <- function(results, years) {
  results %>%
    filter(year <= 2021 + years) %>%
    group_by(country, ISO3, scenario) %>%
    summarise(
      Mean_YLL = mean(YLL, na.rm = TRUE),
      Mean_YLL_disc = mean(YLL_disc, na.rm = TRUE),
      Mean_DALY_clinical = mean(DALY_clinical, na.rm = TRUE),
      Mean_DALY_resistant = mean(DALY_resistant, na.rm = TRUE),
      Mean_DALY_death = mean(DALY_death, na.rm = TRUE),
      Mean_DALY_death_disc = mean(DALY_death_disc, na.rm = TRUE),
      Mean_DALY_total = mean(DALY_total, na.rm = TRUE),
      Mean_DALY_total_disc = mean(DALY_total_disc, na.rm = TRUE),
      Mean_Cost_VSL = mean(Cost_VSL, na.rm = TRUE),
      Mean_Cost_VSLY = mean(Cost_VSLY, na.rm = TRUE),
      Mean_Cost_DALY_total_disc = mean(Cost_DALY_total_disc, na.rm = TRUE),
      Mean_Cost_clin_gvt = mean(Cost_clin_gvt, na.rm = TRUE),
      Mean_Cost_clin_gvt_disc = mean(Cost_clin_gvt_disc, na.rm = TRUE),
      Mean_Cost_clin_oop = mean(Cost_clin_oop, na.rm = TRUE),
      Mean_Cost_clin_oop_disc = mean(Cost_clin_oop_disc, na.rm = TRUE),
      Mean_Cost_hosp_resistant_gvt = mean(Cost_hosp_resistant_gvt, na.rm = TRUE),
      Mean_Cost_hosp_resistant_gvt_disc = mean(Cost_hosp_resistant_gvt_disc, na.rm = TRUE),
      Mean_Cost_hosp_resistant_oop = mean(Cost_hosp_resistant_oop, na.rm = TRUE),
      Mean_Cost_hosp_resistant_oop_disc = mean(Cost_hosp_resistant_oop_disc, na.rm = TRUE),
      Mean_Cost_prod_hospital = mean(Cost_prod_hospital, na.rm = TRUE),
      Mean_Cost_prod_resistant_hospital = mean(Cost_prod_resistant_hospital, na.rm = TRUE),
      Mean_Cost_prod_death = mean(Cost_prod_death, na.rm = TRUE),
      Mean_Cost_prod_total = mean(Cost_prod_total, na.rm = TRUE),
      Mean_Cost_prod_total_disc = mean(Cost_prod_total_disc, na.rm = TRUE),
      Mean_Cost_societal = mean(Cost_societal, na.rm = TRUE),
      Mean_Cost_societal_disc = mean(Cost_societal_disc, na.rm = TRUE),
      Mean_Catastrophic_cases = mean(Catastrophic_cases, na.rm = TRUE),
      Mean_Impoverished_cases = mean(Impoverished_cases, na.rm = TRUE),
      Mean_Cases_averted = mean(Cases_averted, na.rm = TRUE),
      Mean_Resistant_cases_averted = mean(Resistant_cases_averted, na.rm = TRUE),
      Mean_Deaths_averted = mean(Deaths_averted, na.rm = TRUE)
    )
}

# Calculate cumulative statistics for different time horizons from the combined results
summary_stats_3_years <- calculate_cumulative_stats(all_results, 3)
summary_stats_5_years <- calculate_cumulative_stats(all_results, 5)
summary_stats_10_years <- calculate_cumulative_stats(all_results, 10)

# Save the summarised data to CSV files for each time horizon
output_dir_csv <- "/Users/X/Desktop/University of Oxford/3. Trinity Term (Placement)/10. Output Files"
write.csv(summary_stats_3_years, file.path(output_dir_csv, "summary_stats_3_years.csv"), row.names = FALSE)
write.csv(summary_stats_5_years, file.path(output_dir_csv, "summary_stats_5_years.csv"), row.names = FALSE)
write.csv(summary_stats_10_years, file.path(output_dir_csv, "summary_stats_10_years.csv"), row.names = FALSE)

#### SECTION 3.2 - RESULTS ANALYSIS ####

#### SECTION 3.2.1 - AFRICAN COUNTRY MAPPING ####
output_dir_graphs <- "/Users/X/Desktop/University of Oxford/3. Trinity Term (Placement)/11. Output Graphs"

# Load world map data and filter to include only African countries
world <- ne_countries(scale = "medium", returnclass = "sf")
africa <- world %>% filter(continent == "Africa")

# Define the function to calculate summary statistics for various outcomes (aggregate values)
calculate_summary_stats_aggregate <- function(all_results, time_horizon, scenario) {
  summary_stats <- all_results %>%
    filter(year <= 2021 + time_horizon, scenario == !!scenario) %>%
    group_by(country, ISO3, scenario) %>%
    summarise(
      Mean_DALYs_Averted = mean(DALY_total_disc, na.rm = TRUE),
      Lower_CI_DALYs_Averted = quantile(DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_DALYs_Averted = quantile(DALY_total_disc, 0.975, na.rm = TRUE),
      Mean_Cost_DALYs_Averted = mean(Cost_DALY_total_disc, na.rm = TRUE),
      Lower_CI_Cost_DALYs_Averted = quantile(Cost_DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_Cost_DALYs_Averted = quantile(Cost_DALY_total_disc, 0.975, na.rm = TRUE),
      Mean_Societal_Cost_Averted = mean(Cost_societal_disc, na.rm = TRUE),
      Lower_CI_Societal_Cost_Averted = quantile(Cost_societal_disc, 0.025, na.rm = TRUE),
      Upper_CI_Societal_Cost_Averted = quantile(Cost_societal_disc, 0.975, na.rm = TRUE),
      Mean_Catastrophic_Cases = mean(Catastrophic_cases, na.rm = TRUE),
      Lower_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.025, na.rm = TRUE),
      Upper_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.975, na.rm = TRUE),
      Mean_Impoverished_Cases = mean(Impoverished_cases, na.rm = TRUE),
      Lower_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.025, na.rm = TRUE),
      Upper_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.975, na.rm = TRUE)
    )
  return(summary_stats)
}

# Initialize global min and max values for the different outcome metrics
global_min_aggregate <- Inf
global_max_aggregate <- -Inf
global_min_cost_aggregate <- Inf
global_max_cost_aggregate <- -Inf
global_min_catastrophic_cases_aggregate <- Inf
global_max_catastrophic_cases_aggregate <- -Inf
global_min_impoverished_cases_aggregate <- Inf
global_max_impoverished_cases_aggregate <- -Inf
global_min_societal_cost_aggregate <- Inf
global_max_societal_cost_aggregate <- -Inf

# List to store summary statistics for different scenarios and time horizons
summary_stats_list_aggregate <- list()

# Calculate global min and max values across all time horizons and scenarios (for aggregate data)
for (time_horizon in time_horizons) {
  for (scenario in scenarios) {
    summary_stats_aggregate <- calculate_summary_stats_aggregate(all_results, time_horizon, scenario)
    summary_stats_list_aggregate[[paste0(time_horizon, "_", scenario)]] <- summary_stats_aggregate
    
    global_min_aggregate <- min(global_min_aggregate, min(summary_stats_aggregate$Mean_DALYs_Averted, na.rm = TRUE))
    global_max_aggregate <- max(global_max_aggregate, max(summary_stats_aggregate$Mean_DALYs_Averted, na.rm = TRUE))
    global_min_cost_aggregate <- min(global_min_cost_aggregate, min(summary_stats_aggregate$Mean_Cost_DALYs_Averted, na.rm = TRUE))
    global_max_cost_aggregate <- max(global_max_cost_aggregate, max(summary_stats_aggregate$Mean_Cost_DALYs_Averted, na.rm = TRUE))
    global_min_catastrophic_cases_aggregate <- min(global_min_catastrophic_cases_aggregate, min(summary_stats_aggregate$Mean_Catastrophic_Cases, na.rm = TRUE))
    global_max_catastrophic_cases_aggregate <- max(global_max_catastrophic_cases_aggregate, max(summary_stats_aggregate$Mean_Catastrophic_Cases, na.rm = TRUE))
    global_min_impoverished_cases_aggregate <- min(global_min_impoverished_cases_aggregate, min(summary_stats_aggregate$Mean_Impoverished_Cases, na.rm = TRUE))
    global_max_impoverished_cases_aggregate <- max(global_max_impoverished_cases_aggregate, max(summary_stats_aggregate$Mean_Impoverished_Cases, na.rm = TRUE))
    global_min_societal_cost_aggregate <- min(global_min_societal_cost_aggregate, min(summary_stats_aggregate$Mean_Societal_Cost_Averted, na.rm = TRUE))
    global_max_societal_cost_aggregate <- max(global_max_societal_cost_aggregate, max(summary_stats_aggregate$Mean_Societal_Cost_Averted, na.rm = TRUE))
  }
}

# Define a function to generate and return a ggplot object for aggregate data with larger legends
create_plot_aggregate <- function(data, fill_var, fill_label, title, min_val, max_val) {
  ggplot(data) +
    geom_sf(aes_string(fill = fill_var)) +
    scale_fill_gradientn(colors = c("navy", "turquoise", "green", "yellow", "orange", "red"), na.value = "grey", name = fill_label,
                         labels = scales::label_number(scale = 1, suffix = ""),
                         limits = c(min_val, max_val)) +
    labs(title = title,
         fill = fill_label) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), # Remove latitude and longitude lines
          legend.position = "right", # Place legend on the right
          legend.text = element_text(size = 14), # Increase the size of the legend text
          legend.title = element_text(size = 16)) + # Increase the size of the legend title
    annotation_scale(location = "bl", width_hint = 0.5) + # Add a scale bar in the bottom left
    annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) # Add a detailed north arrow in the top right
}

# Loop through each time horizon and create and save combined plots for all VE scenarios (for aggregate data)
for (time_horizon in time_horizons) {
  # Initialize an empty list to store plots for each outcome metric
  plot_lists_aggregate <- list("DALYs_Averted" = list(),
                               "Societal_Cost" = list(),
                               "Catastrophic_Cases" = list(),
                               "Impoverished_Cases" = list())
  
  for (scenario in scenarios) {
    # Get the summary statistics for the current scenario and time horizon
    summary_stats_aggregate <- summary_stats_list_aggregate[[paste0(time_horizon, "_", scenario)]]
    
    # Merge the map data with the summary data
    africa_data_aggregate <- merge(africa, summary_stats_aggregate, by.x = "iso_a3", by.y = "ISO3", all.x = TRUE)
    
    # Generate the DALYs Averted plot for the current VE scenario
    p_DALYs_aggregate <- create_plot_aggregate(africa_data_aggregate, "Mean_DALYs_Averted", "Total DALYs Averted",
                                               paste0("Scenario ", scenario, " (", time_horizon, " Years)"), global_min_aggregate, global_max_aggregate)
    plot_lists_aggregate$DALYs_Averted[[scenario]] <- p_DALYs_aggregate
    
    # Generate the Societal Cost Averted plot for the current VE scenario
    p_Societal_Cost_aggregate <- create_plot_aggregate(africa_data_aggregate, "Mean_Societal_Cost_Averted", "Societal Cost Averted (USD)",
                                                       paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                                       global_min_societal_cost_aggregate, global_max_societal_cost_aggregate)
    plot_lists_aggregate$Societal_Cost[[scenario]] <- p_Societal_Cost_aggregate
    
    # Generate the Catastrophic Cases Averted plot for the current VE scenario
    p_Catastrophic_Cases_aggregate <- create_plot_aggregate(africa_data_aggregate, "Mean_Catastrophic_Cases", "Total Catastrophic Cases Averted",
                                                            paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                                            global_min_catastrophic_cases_aggregate, global_max_catastrophic_cases_aggregate)
    plot_lists_aggregate$Catastrophic_Cases[[scenario]] <- p_Catastrophic_Cases_aggregate
    
    # Generate the Impoverished Cases Averted plot for the current VE scenario
    p_Impoverished_Cases_aggregate <- create_plot_aggregate(africa_data_aggregate, "Mean_Impoverished_Cases", "Total Impoverished Cases Averted",
                                                            paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                                            global_min_impoverished_cases_aggregate, global_max_impoverished_cases_aggregate)
    plot_lists_aggregate$Impoverished_Cases[[scenario]] <- p_Impoverished_Cases_aggregate
  }
  
  # Combine and save the DALYs Averted plots for Scenario 1, 2, and 3
  combined_DALYs_aggregate <- ggarrange(plot_lists_aggregate$DALYs_Averted[[1]], plot_lists_aggregate$DALYs_Averted[[2]], plot_lists_aggregate$DALYs_Averted[[3]], 
                                        ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_DALYs_Averted_", time_horizon, "_years.png")), plot = combined_DALYs_aggregate, width = 20, height = 7)
  
  # Combine and save the Societal Cost Averted plots for Scenario 1, 2, and 3
  combined_Societal_Cost_aggregate <- ggarrange(plot_lists_aggregate$Societal_Cost[[1]], plot_lists_aggregate$Societal_Cost[[2]], plot_lists_aggregate$Societal_Cost[[3]], 
                                                ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Societal_Cost_Averted_", time_horizon, "_years.png")), plot = combined_Societal_Cost_aggregate, width = 20, height = 7)
  
  # Combine and save the Catastrophic Cases Averted plots for Scenario 1, 2, and 3
  combined_Catastrophic_Cases_aggregate <- ggarrange(plot_lists_aggregate$Catastrophic_Cases[[1]], plot_lists_aggregate$Catastrophic_Cases[[2]], plot_lists_aggregate$Catastrophic_Cases[[3]], 
                                                     ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Catastrophic_Cases_Averted_", time_horizon, "_years.png")), plot = combined_Catastrophic_Cases_aggregate, width = 20, height = 7)
  
  # Combine and save the Impoverished Cases Averted plots for Scenario 1, 2, and 3
  combined_Impoverished_Cases_aggregate <- ggarrange(plot_lists_aggregate$Impoverished_Cases[[1]], plot_lists_aggregate$Impoverished_Cases[[2]], plot_lists_aggregate$Impoverished_Cases[[3]], 
                                                     ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Impoverished_Cases_Averted_", time_horizon, "_years.png")), plot = combined_Impoverished_Cases_aggregate, width = 20, height = 7)
}


#### SECTION 3.2.1.2 - PER 1000 CHILDREN GRAPHS ####

# Function to calculate summary statistics per 1000 children for various outcomes
calculate_summary_stats_per_1000 <- function(all_results, time_horizon, population_data, scenario) {
  summary_stats <- all_results %>%
    filter(year <= 2021 + time_horizon, scenario == !!scenario) %>%
    group_by(country, ISO3, scenario) %>%
    summarise(
      Mean_DALYs_Averted = mean(DALY_total_disc, na.rm = TRUE),
      Lower_CI_DALYs_Averted = quantile(DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_DALYs_Averted = quantile(DALY_total_disc, 0.975, na.rm = TRUE),
      Mean_Cost_DALYs_Averted = mean(Cost_DALY_total_disc, na.rm = TRUE),
      Lower_CI_Cost_DALYs_Averted = quantile(Cost_DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_Cost_DALYs_Averted = quantile(Cost_DALY_total_disc, 0.975, na.rm = TRUE),
      Mean_Societal_Cost_Averted = mean(Cost_societal_disc, na.rm = TRUE),
      Lower_CI_Societal_Cost_Averted = quantile(Cost_societal_disc, 0.025, na.rm = TRUE),
      Upper_CI_Societal_Cost_Averted = quantile(Cost_societal_disc, 0.975, na.rm = TRUE),
      Mean_Catastrophic_Cases = mean(Catastrophic_cases, na.rm = TRUE),
      Lower_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.025, na.rm = TRUE),
      Upper_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.975, na.rm = TRUE),
      Mean_Impoverished_Cases = mean(Impoverished_cases, na.rm = TRUE),
      Lower_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.025, na.rm = TRUE),
      Upper_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.975, na.rm = TRUE),
      Mean_Cases_Averted = mean(Cases_averted, na.rm = TRUE),
      Lower_CI_Cases_Averted = quantile(Cases_averted, 0.025, na.rm = TRUE),
      Upper_CI_Cases_Averted = quantile(Cases_averted, 0.975, na.rm = TRUE),
      Mean_Resistant_Cases_Averted = mean(Resistant_cases_averted, na.rm = TRUE),
      Lower_CI_Resistant_Cases_Averted = quantile(Resistant_cases_averted, 0.025, na.rm = TRUE),
      Upper_CI_Resistant_Cases_Averted = quantile(Resistant_cases_averted, 0.975, na.rm = TRUE),
    ) %>%
    left_join(population_data, by = c("ISO3" = "Country.Code")) %>%
    mutate(
      DALYs_Averted_per_1000 = (Mean_DALYs_Averted / Pop1) * 1000,
      Lower_CI_DALYs_Averted_per_1000 = (Lower_CI_DALYs_Averted / Pop1) * 1000,
      Upper_CI_DALYs_Averted_per_1000 = (Upper_CI_DALYs_Averted / Pop1) * 1000,
      Cost_DALYs_Averted_per_1000 = (Mean_Cost_DALYs_Averted / Pop1) * 1000,
      Lower_CI_Cost_DALYs_Averted_per_1000 = (Lower_CI_Cost_DALYs_Averted / Pop1) * 1000,
      Upper_CI_Cost_DALYs_Averted_per_1000 = (Upper_CI_Cost_DALYs_Averted / Pop1) * 1000,
      Societal_Cost_per_1000 = (Mean_Societal_Cost_Averted / Pop1) * 1000,
      Lower_CI_Societal_Cost_per_1000 = (Lower_CI_Societal_Cost_Averted / Pop1) * 1000,
      Upper_CI_Societal_Cost_per_1000 = (Upper_CI_Societal_Cost_Averted / Pop1) * 1000,
      Catastrophic_Cases_per_1000 = (Mean_Catastrophic_Cases / Pop1) * 1000,
      Lower_CI_Catastrophic_Cases_per_1000 = (Lower_CI_Catastrophic_Cases / Pop1) * 1000,
      Upper_CI_Catastrophic_Cases_per_1000 = (Upper_CI_Catastrophic_Cases / Pop1) * 1000,
      Impoverished_Cases_per_1000 = (Mean_Impoverished_Cases / Pop1) * 1000,
      Lower_CI_Impoverished_Cases_per_1000 = (Lower_CI_Impoverished_Cases / Pop1) * 1000,
      Upper_CI_Impoverished_Cases_per_1000 = (Upper_CI_Impoverished_Cases / Pop1) * 1000,
      Cases_Averted_per_1000 = (Mean_Cases_Averted / Pop1) * 1000,
      Lower_CI_Cases_Averted_per_1000 = (Lower_CI_Cases_Averted / Pop1) * 1000,
      Upper_CI_Cases_Averted_per_1000 = (Upper_CI_Cases_Averted / Pop1) * 1000,
      Resistant_Cases_Averted_per_1000 = (Mean_Resistant_Cases_Averted / Pop1) * 1000,
      Lower_CI_Resistant_Cases_Averted_per_1000 = (Lower_CI_Resistant_Cases_Averted / Pop1) * 1000,
      Upper_CI_Resistant_Cases_Averted_per_1000 = (Upper_CI_Resistant_Cases_Averted / Pop1) * 1000,
    )
  return(summary_stats)
}

# Initialize global min and max values for the different outcome metrics
global_min_per_1000 <- Inf
global_max_per_1000 <- -Inf
global_min_cost_per_1000 <- Inf
global_max_cost_per_1000 <- -Inf
global_min_catastrophic_cases_per_1000 <- Inf
global_max_catastrophic_cases_per_1000 <- -Inf
global_min_impoverished_cases_per_1000 <- Inf
global_max_impoverished_cases_per_1000 <- -Inf
global_min_societal_cost_per_1000 <- Inf
global_max_societal_cost_per_1000 <- -Inf

# List to store summary statistics for different scenarios and time horizons
summary_stats_list_per_1000 <- list()

# Calculate global min and max values across all time horizons and scenarios
for (time_horizon in time_horizons) {
  for (scenario in scenarios) {
    summary_stats_per_1000 <- calculate_summary_stats_per_1000(all_results, time_horizon, population_data, scenario)
    summary_stats_list_per_1000[[paste0(time_horizon, "_", scenario)]] <- summary_stats_per_1000
    global_min_per_1000 <- min(global_min_per_1000, min(summary_stats_per_1000$DALYs_Averted_per_1000, na.rm = TRUE))
    global_max_per_1000 <- max(global_max_per_1000, max(summary_stats_per_1000$DALYs_Averted_per_1000, na.rm = TRUE))
    global_min_cost_per_1000 <- min(global_min_cost_per_1000, min(summary_stats_per_1000$Cost_DALYs_Averted_per_1000, na.rm = TRUE))
    global_max_cost_per_1000 <- max(global_max_cost_per_1000, max(summary_stats_per_1000$Cost_DALYs_Averted_per_1000, na.rm = TRUE))
    global_min_catastrophic_cases_per_1000 <- min(global_min_catastrophic_cases_per_1000, min(summary_stats_per_1000$Catastrophic_Cases_per_1000, na.rm = TRUE))
    global_max_catastrophic_cases_per_1000 <- max(global_max_catastrophic_cases_per_1000, max(summary_stats_per_1000$Catastrophic_Cases_per_1000, na.rm = TRUE))
    global_min_impoverished_cases_per_1000 <- min(global_min_impoverished_cases_per_1000, min(summary_stats_per_1000$Impoverished_Cases_per_1000, na.rm = TRUE))
    global_max_impoverished_cases_per_1000 <- max(global_max_impoverished_cases_per_1000, max(summary_stats_per_1000$Impoverished_Cases_per_1000, na.rm = TRUE))
    global_min_societal_cost_per_1000 <- min(global_min_societal_cost_per_1000, min(summary_stats_per_1000$Societal_Cost_per_1000, na.rm = TRUE))
    global_max_societal_cost_per_1000 <- max(global_max_societal_cost_per_1000, max(summary_stats_per_1000$Societal_Cost_per_1000, na.rm = TRUE))
  }
}

# Define a function to generate and return a ggplot object with larger legends
create_plot <- function(data, fill_var, fill_label, title, min_val, max_val) {
  ggplot(data) +
    geom_sf(aes_string(fill = fill_var)) +
    scale_fill_gradientn(colors = c("navy", "turquoise", "green", "yellow", "orange", "red"), na.value = "grey", name = fill_label,
                         labels = scales::label_number(scale = 1, suffix = ""),
                         limits = c(min_val, max_val)) +
    labs(title = title,
         fill = fill_label) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), # Remove latitude and longitude lines
          legend.position = "right", # Place legend on the right
          legend.text = element_text(size = 14), # Increase the size of the legend text
          legend.title = element_text(size = 16)) + # Increase the size of the legend title
    annotation_scale(location = "bl", width_hint = 0.5) + # Add a scale bar in the bottom left
    annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) # Add a detailed north arrow in the top right
}

# Now, when you generate the plots, the legends will be larger.


# Loop through each time horizon and create and save combined plots for all VE scenarios
for (time_horizon in time_horizons) {
  # Initialize an empty list to store plots for each outcome metric
  plot_lists <- list("DALYs_Averted" = list(),
                     "Societal_Cost" = list(),
                     "Catastrophic_Cases" = list(),
                     "Impoverished_Cases" = list())
  
  for (scenario in scenarios) {
    # Get the summary statistics per 1000 children for the current scenario and time horizon
    summary_stats_per_1000 <- summary_stats_list_per_1000[[paste0(time_horizon, "_", scenario)]]
    
    # Merge the map data with the summary data
    africa_data_per_1000 <- merge(africa, summary_stats_per_1000, by.x = "iso_a3", by.y = "ISO3", all.x = TRUE)
    
    # Generate the DALYs Averted plot for the current VE scenario
    p_DALYs <- create_plot(africa_data_per_1000, "DALYs_Averted_per_1000", "DALYs Averted per 1000 Children",
                           paste0("Scenario ", scenario, " (", time_horizon, " Years)"), global_min_per_1000, global_max_per_1000)
    plot_lists$DALYs_Averted[[scenario]] <- p_DALYs
    
    # Generate the Societal Cost Averted plot for the current VE scenario
    p_Societal_Cost <- create_plot(africa_data_per_1000, "Societal_Cost_per_1000", "Societal Cost Averted per 1000 Children (USD)",
                                   paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                   global_min_societal_cost_per_1000, global_max_societal_cost_per_1000)
    plot_lists$Societal_Cost[[scenario]] <- p_Societal_Cost
    
    # Generate the Catastrophic Cases Averted plot for the current VE scenario
    p_Catastrophic_Cases <- create_plot(africa_data_per_1000, "Catastrophic_Cases_per_1000", "Catastrophic Cases Averted per 1000 Children",
                                        paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                        global_min_catastrophic_cases_per_1000, global_max_catastrophic_cases_per_1000)
    plot_lists$Catastrophic_Cases[[scenario]] <- p_Catastrophic_Cases
    
    # Generate the Impoverished Cases Averted plot for the current VE scenario
    p_Impoverished_Cases <- create_plot(africa_data_per_1000, "Impoverished_Cases_per_1000", "Impoverished Cases Averted per 1000 Children",
                                        paste0("Scenario ", scenario, " (", time_horizon, " Years)"),
                                        global_min_impoverished_cases_per_1000, global_max_impoverished_cases_per_1000)
    plot_lists$Impoverished_Cases[[scenario]] <- p_Impoverished_Cases
  }
  
  # Combine and save the DALYs Averted plots for Scenario 1, 2, and 3
  combined_DALYs <- ggarrange(plot_lists$DALYs_Averted[[1]], plot_lists$DALYs_Averted[[2]], plot_lists$DALYs_Averted[[3]], 
                              ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_DALYs_Averted_per_1000_", time_horizon, "_years.png")), plot = combined_DALYs, width = 20, height = 7)
  
  # Combine and save the Societal Cost Averted plots for Scenario 1, 2, and 3
  combined_Societal_Cost <- ggarrange(plot_lists$Societal_Cost[[1]], plot_lists$Societal_Cost[[2]], plot_lists$Societal_Cost[[3]], 
                                      ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Societal_Cost_Averted_per_1000_", time_horizon, "_years.png")), plot = combined_Societal_Cost, width = 20, height = 7)
  
  # Combine and save the Catastrophic Cases Averted plots for Scenario 1, 2, and 3
  combined_Catastrophic_Cases <- ggarrange(plot_lists$Catastrophic_Cases[[1]], plot_lists$Catastrophic_Cases[[2]], plot_lists$Catastrophic_Cases[[3]], 
                                           ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Catastrophic_Cases_Averted_per_1000_", time_horizon, "_years.png")), plot = combined_Catastrophic_Cases, width = 20, height = 7)
  
  # Combine and save the Impoverished Cases Averted plots for Scenario 1, 2, and 3
  combined_Impoverished_Cases <- ggarrange(plot_lists$Impoverished_Cases[[1]], plot_lists$Impoverished_Cases[[2]], plot_lists$Impoverished_Cases[[3]], 
                                           ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
  ggsave(filename = file.path(output_dir_graphs, paste0("Combined_Impoverished_Cases_Averted_per_1000_", time_horizon, "_years.png")), plot = combined_Impoverished_Cases, width = 20, height = 7)
}


#### SECTION 3.2.2 - AGGREGATE STATISTICS ####

#### SECTION 3.2.2.1 - SUMMARY STATISTICS FOR BASE SCENARIO ####

# Define the function to calculate summary statistics for the baseline scenario
calculate_summary_stats_base <- function(results) {
  results %>%
    group_by(country, ISO3) %>%
    summarise(
      Mean_Societal_Cost = mean(Cost_societal_disc, na.rm = TRUE),
      Lower_CI_Societal_Cost = quantile(Cost_societal_disc, 0.025, na.rm = TRUE),
      Upper_CI_Societal_Cost = quantile(Cost_societal_disc, 0.975, na.rm = TRUE),
      
      Mean_Resistant_OOP = mean(Cost_hosp_resistant_oop_disc, na.rm = TRUE),
      Lower_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.025, na.rm = TRUE),
      Upper_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.975, na.rm = TRUE),
      
      Mean_DALYs = mean(DALY_total_disc, na.rm = TRUE),
      Lower_CI_DALYs = quantile(DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_DALYs = quantile(DALY_total_disc, 0.975, na.rm = TRUE),
      
      Mean_Cost_VSL = mean(Cost_VSLY, na.rm = TRUE),
      Lower_CI_Cost_VSL = quantile(Cost_VSLY, 0.025, na.rm = TRUE),
      Upper_CI_Cost_VSL = quantile(Cost_VSLY, 0.975, na.rm = TRUE),
      
      Mean_YLL_Lost = mean(YLL_disc, na.rm = TRUE),
      Lower_CI_YLL_Lost = quantile(YLL_disc, 0.025, na.rm = TRUE),
      Upper_CI_YLL_Lost = quantile(YLL_disc, 0.975, na.rm = TRUE),
      
      Mean_Catastrophic_Cases = mean(Catastrophic_cases, na.rm = TRUE),
      Lower_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.025, na.rm = TRUE),
      Upper_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.975, na.rm = TRUE),
      
      Mean_Impoverished_Cases = mean(Impoverished_cases, na.rm = TRUE),
      Lower_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.025, na.rm = TRUE),
      Upper_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.975, na.rm = TRUE),
      
      Mean_Cost_DALYs_Disc = mean(Cost_DALY_total_disc, na.rm = TRUE),
      Lower_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.975, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Calculate summary statistics for the baseline scenario
summary_stats_base <- calculate_summary_stats_base(results_base_all)

# Aggregate statistics across all countries by summing
aggregate_stats_base <- summary_stats_base %>%
  summarise(
    Sum_Societal_Cost = sum(Mean_Societal_Cost, na.rm = TRUE),
    Lower_CI_Societal_Cost = sum(Lower_CI_Societal_Cost, na.rm = TRUE),
    Upper_CI_Societal_Cost = sum(Upper_CI_Societal_Cost, na.rm = TRUE),
    
    Sum_Resistant_OOP = sum(Mean_Resistant_OOP, na.rm = TRUE),
    Lower_CI_Resistant_OOP = sum(Lower_CI_Resistant_OOP, na.rm = TRUE),
    Upper_CI_Resistant_OOP = sum(Upper_CI_Resistant_OOP, na.rm = TRUE),
    
    Sum_DALYs = sum(Mean_DALYs, na.rm = TRUE),
    Lower_CI_DALYs = sum(Lower_CI_DALYs, na.rm = TRUE),
    Upper_CI_DALYs = sum(Upper_CI_DALYs, na.rm = TRUE),
    
    Sum_Cost_VSL = sum(Mean_Cost_VSL, na.rm = TRUE),
    Lower_CI_Cost_VSL = sum(Lower_CI_Cost_VSL, na.rm = TRUE),
    Upper_CI_Cost_VSL = sum(Upper_CI_Cost_VSL, na.rm = TRUE),
    
    Sum_YLL_Lost = sum(Mean_YLL_Lost, na.rm = TRUE),
    Lower_CI_YLL_Lost = sum(Lower_CI_YLL_Lost, na.rm = TRUE),
    Upper_CI_YLL_Lost = sum(Upper_CI_YLL_Lost, na.rm = TRUE),
    
    Sum_Catastrophic_Cases = sum(Mean_Catastrophic_Cases, na.rm = TRUE),
    Lower_CI_Catastrophic_Cases = sum(Lower_CI_Catastrophic_Cases, na.rm = TRUE),
    Upper_CI_Catastrophic_Cases = sum(Upper_CI_Catastrophic_Cases, na.rm = TRUE),
    
    Sum_Impoverished_Cases = sum(Mean_Impoverished_Cases, na.rm = TRUE),
    Lower_CI_Impoverished_Cases = sum(Lower_CI_Impoverished_Cases, na.rm = TRUE),
    Upper_CI_Impoverished_Cases = sum(Upper_CI_Impoverished_Cases, na.rm = TRUE),
    
    Sum_Cost_DALYs_Disc = sum(Mean_Cost_DALYs_Disc, na.rm = TRUE),
    Lower_CI_Cost_DALYs_Disc = sum(Lower_CI_Cost_DALYs_Disc, na.rm = TRUE),
    Upper_CI_Cost_DALYs_Disc = sum(Upper_CI_Cost_DALYs_Disc, na.rm = TRUE)
  )

write.csv(aggregate_stats_base,"aggregate_stats_base.csv")

# Create tables for each statistic and save as CSV
table_list_base <- list(
  "Societal_Cost_Base" = aggregate_stats_base %>%
    select(Sum_Societal_Cost, Lower_CI_Societal_Cost, Upper_CI_Societal_Cost),
  "Resistant_OOP_Base" = aggregate_stats_base %>%
    select(Sum_Resistant_OOP, Lower_CI_Resistant_OOP, Upper_CI_Resistant_OOP),
  "DALYs_Base" = aggregate_stats_base %>%
    select(Sum_DALYs, Lower_CI_DALYs, Upper_CI_DALYs),
  "Cost_VSL_Base" = aggregate_stats_base %>%
    select(Sum_Cost_VSL, Lower_CI_Cost_VSL, Upper_CI_Cost_VSL),
  "YLL_Lost_Base" = aggregate_stats_base %>%
    select(Sum_YLL_Lost, Lower_CI_YLL_Lost, Upper_CI_YLL_Lost),
  "Catastrophic_Cases_Base" = aggregate_stats_base %>%
    select(Sum_Catastrophic_Cases, Lower_CI_Catastrophic_Cases, Upper_CI_Catastrophic_Cases),
  "Impoverished_Cases_Base" = aggregate_stats_base %>%
    select(Sum_Impoverished_Cases, Lower_CI_Impoverished_Cases, Upper_CI_Impoverished_Cases),
  "Cost_DALYs_Disc_Base" = aggregate_stats_base %>%
    select(Sum_Cost_DALYs_Disc, Lower_CI_Cost_DALYs_Disc, Upper_CI_Cost_DALYs_Disc)
)

for (table_name in names(table_list_base)) {
  table <- table_list_base[[table_name]]
  write.csv(table, file.path(output_dir_csv, paste0(table_name, ".csv")), row.names = FALSE)
  print(knitr::kable(table, format = "html", table.attr = "class='table table-striped'", caption = paste("Table:", table_name)))
}

#### SECTION 3.2.2.2- REGION SEGREGATION ####

# Define the regions and their corresponding countries
regions <- list(
  "Northern Africa" = c("Algeria", "Egypt", "Libya", "Morocco", "Sudan", "Tunisia", "Western Sahara"),
  "Eastern Africa" = c("Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mauritius", "Mozambique", "Rwanda", "Seychelles", "Somalia", "South Sudan", "United Republic of Tanzania", "Uganda", "Zambia", "Zimbabwe"),
  "Middle Africa" = c("Angola", "Cameroon", "Central African Republic", "Chad", "Congo", "Democratic Republic of the Congo", "Equatorial Guinea", "Gabon", "Sao Tome and Principe"),
  "Southern Africa" = c("Botswana", "Eswatini", "Lesotho", "Namibia", "South Africa"),
  "Western Africa" = c("Benin", "Burkina Faso", "Cape Verde", "Côte d'Ivoire", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")
)

# Add region information to the data
results_base_all <- results_base_all %>%
  mutate(region = case_when(
    country %in% regions$`Northern Africa` ~ "Northern Africa",
    country %in% regions$`Eastern Africa` ~ "Eastern Africa",
    country %in% regions$`Middle Africa` ~ "Middle Africa",
    country %in% regions$`Southern Africa` ~ "Southern Africa",
    country %in% regions$`Western Africa` ~ "Western Africa",
    TRUE ~ "Other"
  ))

# Define the function to calculate summary statistics for the baseline scenario by region
calculate_summary_stats_base_region <- function(results) {
  results %>%
    group_by(region) %>%
    summarise(
      Mean_Societal_Cost = mean(Cost_societal_disc, na.rm = TRUE),
      Lower_CI_Societal_Cost = quantile(Cost_societal_disc, 0.025, na.rm = TRUE),
      Upper_CI_Societal_Cost = quantile(Cost_societal_disc, 0.975, na.rm = TRUE),
      
      Mean_Resistant_OOP = mean(Cost_hosp_resistant_oop_disc, na.rm = TRUE),
      Lower_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.025, na.rm = TRUE),
      Upper_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.975, na.rm = TRUE),
      
      Mean_DALYs = mean(DALY_total_disc, na.rm = TRUE),
      Lower_CI_DALYs = quantile(DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_DALYs = quantile(DALY_total_disc, 0.975, na.rm = TRUE),
      
      Mean_Cost_VSL = mean(Cost_VSLY, na.rm = TRUE),
      Lower_CI_Cost_VSL = quantile(Cost_VSLY, 0.025, na.rm = TRUE),
      Upper_CI_Cost_VSL = quantile(Cost_VSLY, 0.975, na.rm = TRUE),
      
      Mean_YLL_Lost = mean(YLL_disc, na.rm = TRUE),
      Lower_CI_YLL_Lost = quantile(YLL_disc, 0.025, na.rm = TRUE),
      Upper_CI_YLL_Lost = quantile(YLL_disc, 0.975, na.rm = TRUE),
      
      Mean_Catastrophic_Cases = mean(Catastrophic_cases, na.rm = TRUE),
      Lower_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.025, na.rm = TRUE),
      Upper_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.975, na.rm = TRUE),
      
      Mean_Impoverished_Cases = mean(Impoverished_cases, na.rm = TRUE),
      Lower_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.025, na.rm = TRUE),
      Upper_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.975, na.rm = TRUE),
      
      Mean_Cost_DALYs_Disc = mean(Cost_DALY_total_disc, na.rm = TRUE),
      Lower_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.025, na.rm = TRUE),
      Upper_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.975, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Calculate summary statistics for the baseline scenario by region
summary_stats_base_region <- calculate_summary_stats_base_region(results_base_all)

table_list_base_region <- list(
  "Societal_Cost_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Societal_Cost, Lower_CI_Societal_Cost, Upper_CI_Societal_Cost),
  "Resistant_OOP_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Resistant_OOP, Lower_CI_Resistant_OOP, Upper_CI_Resistant_OOP),
  "DALYs_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_DALYs, Lower_CI_DALYs, Upper_CI_DALYs),
  "Cost_VSL_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Cost_VSL, Lower_CI_Cost_VSL, Upper_CI_Cost_VSL),
  "YLL_Lost_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_YLL_Lost, Lower_CI_YLL_Lost, Upper_CI_YLL_Lost),
  "Catastrophic_Cases_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Catastrophic_Cases, Lower_CI_Catastrophic_Cases, Upper_CI_Catastrophic_Cases),
  "Impoverished_Cases_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Impoverished_Cases, Lower_CI_Impoverished_Cases, Upper_CI_Impoverished_Cases),
  "Cost_DALYs_Disc_Base_Region" = summary_stats_base_region %>%
    select(region, Mean_Cost_DALYs_Disc, Lower_CI_Cost_DALYs_Disc, Upper_CI_Cost_DALYs_Disc)
)

# Save tables as CSV and display
for (table_name in names(table_list_base_region)) {
  table <- table_list_base_region[[table_name]]
  write.csv(table, file.path(output_dir_csv, paste0(table_name, ".csv")), row.names = FALSE)
  print(knitr::kable(table, format = "html", table.attr = "class='table table-striped'", caption = paste("Table:", table_name)))
}


#### SECTION 3.2.2.3 - VACCINES ####

# Calculate summary statistics for each country and scenario
country_scenario_stats <- all_results %>%
  group_by(country, ISO3, scenario) %>%
  summarise(
    Mean_Societal_Cost = mean(Cost_societal_disc, na.rm = TRUE),
    Lower_CI_Societal_Cost = quantile(Cost_societal_disc, 0.025, na.rm = TRUE),
    Upper_CI_Societal_Cost = quantile(Cost_societal_disc, 0.975, na.rm = TRUE),
    
    Mean_Resistant_OOP = mean(Cost_hosp_resistant_oop_disc, na.rm = TRUE),
    Lower_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.025, na.rm = TRUE),
    Upper_CI_Resistant_OOP = quantile(Cost_hosp_resistant_oop_disc, 0.975, na.rm = TRUE),
    
    Mean_DALYs_Averted = mean(DALY_total_disc, na.rm = TRUE),
    Lower_CI_DALYs_Averted = quantile(DALY_total_disc, 0.025, na.rm = TRUE),
    Upper_CI_DALYs_Averted = quantile(DALY_total_disc, 0.975, na.rm = TRUE),
    
    Mean_Cost_VSL = mean(Cost_VSLY, na.rm = TRUE),
    Lower_CI_Cost_VSL = quantile(Cost_VSLY, 0.025, na.rm = TRUE),
    Upper_CI_Cost_VSL = quantile(Cost_VSLY, 0.975, na.rm = TRUE),
    
    Mean_YLL_Lost = mean(YLL_disc, na.rm = TRUE),
    Lower_CI_YLL_Lost = quantile(YLL_disc, 0.025, na.rm = TRUE),
    Upper_CI_YLL_Lost = quantile(YLL_disc, 0.975, na.rm = TRUE),
    
    Mean_Catastrophic_Cases = mean(Catastrophic_cases, na.rm = TRUE),
    Lower_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.025, na.rm = TRUE),
    Upper_CI_Catastrophic_Cases = quantile(Catastrophic_cases, 0.975, na.rm = TRUE),
    
    Mean_Impoverished_Cases = mean(Impoverished_cases, na.rm = TRUE),
    Lower_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.025, na.rm = TRUE),
    Upper_CI_Impoverished_Cases = quantile(Impoverished_cases, 0.975, na.rm = TRUE),
    
    Mean_Cost_DALYs_Disc = mean(Cost_DALY_total_disc, na.rm = TRUE),
    Lower_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.025, na.rm = TRUE),
    Upper_CI_Cost_DALYs_Disc = quantile(Cost_DALY_total_disc, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Aggregate statistics across all countries by summing
aggregate_stats <- country_scenario_stats %>%
  group_by(scenario) %>%
  summarise(
    Sum_Societal_Cost = sum(Mean_Societal_Cost, na.rm = TRUE),
    Lower_CI_Societal_Cost = sum(Lower_CI_Societal_Cost, na.rm = TRUE),
    Upper_CI_Societal_Cost = sum(Upper_CI_Societal_Cost, na.rm = TRUE),
    
    Sum_Resistant_OOP = sum(Mean_Resistant_OOP, na.rm = TRUE),
    Lower_CI_Resistant_OOP = sum(Lower_CI_Resistant_OOP, na.rm = TRUE),
    Upper_CI_Resistant_OOP = sum(Upper_CI_Resistant_OOP, na.rm = TRUE),
    
    Sum_DALYs_Averted = sum(Mean_DALYs_Averted, na.rm = TRUE),
    Lower_CI_DALYs_Averted = sum(Lower_CI_DALYs_Averted, na.rm = TRUE),
    Upper_CI_DALYs_Averted = sum(Upper_CI_DALYs_Averted, na.rm = TRUE),
    
    Sum_Cost_VSL = sum(Mean_Cost_VSL, na.rm = TRUE),
    Lower_CI_Cost_VSL = sum(Lower_CI_Cost_VSL, na.rm = TRUE),
    Upper_CI_Cost_VSL = sum(Upper_CI_Cost_VSL, na.rm = TRUE),
    
    Sum_YLL_Lost = sum(Mean_YLL_Lost, na.rm = TRUE),
    Lower_CI_YLL_Lost = sum(Lower_CI_YLL_Lost, na.rm = TRUE),
    Upper_CI_YLL_Lost = sum(Upper_CI_YLL_Lost, na.rm = TRUE),
    
    Sum_Catastrophic_Cases = sum(Mean_Catastrophic_Cases, na.rm = TRUE),
    Lower_CI_Catastrophic_Cases = sum(Lower_CI_Catastrophic_Cases, na.rm = TRUE),
    Upper_CI_Catastrophic_Cases = sum(Upper_CI_Catastrophic_Cases, na.rm = TRUE),
    
    Sum_Impoverished_Cases = sum(Mean_Impoverished_Cases, na.rm = TRUE),
    Lower_CI_Impoverished_Cases = sum(Lower_CI_Impoverished_Cases, na.rm = TRUE),
    Upper_CI_Impoverished_Cases = sum(Upper_CI_Impoverished_Cases, na.rm = TRUE),
    
    Sum_Cost_DALYs_Disc = sum(Mean_Cost_DALYs_Disc, na.rm = TRUE),
    Lower_CI_Cost_DALYs_Disc = sum(Lower_CI_Cost_DALYs_Disc, na.rm = TRUE),
    Upper_CI_Cost_DALYs_Disc = sum(Upper_CI_Cost_DALYs_Disc, na.rm = TRUE)
  )

# Create tables for each statistic and save as CSV
table_list <- list(
  "Societal_Cost" = aggregate_stats %>%
    select(scenario, Sum_Societal_Cost, Lower_CI_Societal_Cost, Upper_CI_Societal_Cost),
  "Resistant_OOP" = aggregate_stats %>%
    select(scenario, Sum_Resistant_OOP, Lower_CI_Resistant_OOP, Upper_CI_Resistant_OOP),
  "DALYs_Averted" = aggregate_stats %>%
    select(scenario, Sum_DALYs_Averted, Lower_CI_DALYs_Averted, Upper_CI_DALYs_Averted),
  "Cost_VSL" = aggregate_stats %>%
    select(scenario, Sum_Cost_VSL, Lower_CI_Cost_VSL, Upper_CI_Cost_VSL),
  "YLL_Lost" = aggregate_stats %>%
    select(scenario, Sum_YLL_Lost, Lower_CI_YLL_Lost, Upper_CI_YLL_Lost),
  "Catastrophic_Cases" = aggregate_stats %>%
    select(scenario, Sum_Catastrophic_Cases, Lower_CI_Catastrophic_Cases, Upper_CI_Catastrophic_Cases),
  "Impoverished_Cases" = aggregate_stats %>%
    select(scenario, Sum_Impoverished_Cases, Lower_CI_Impoverished_Cases, Upper_CI_Impoverished_Cases),
  "Cost_DALYs_Disc" = aggregate_stats %>%
    select(scenario, Sum_Cost_DALYs_Disc, Lower_CI_Cost_DALYs_Disc, Upper_CI_Cost_DALYs_Disc)
)

# Save tables as CSV and display
for (table_name in names(table_list)) {
  table <- table_list[[table_name]]
  write.csv(table, file.path(output_dir_csv, paste0(table_name, ".csv")), row.names = FALSE)
  print(knitr::kable(table, format = "html", table.attr = "class='table table-striped'", caption = paste("Table:", table_name)))
}

#### SECTION 4 - UNIVARIATE SENSITIVITY ANALYSIS ####

#### SECTION 4.1 - UNIVARIATE SENSITIVITY ANALYSIS USING 95% CI ####

# Function to process each country
process_country_sens <- function(country) {
  country_data <- extracted_data %>% filter(ISO3 == !!country)
  
  # Ensure ISO3 is not NA
  if (is.na(country_data$ISO3[1])) {
    message("Skipping country with NA ISO3 code")
    return(NULL)
  }
  
  # Extract VSL for the country
  VSL <- vsl_master %>% filter(GID_0 == country_data$ISO3[1]) %>% pull(vsl_estimate_2021_IntDoll)
  
  # Extract GNI data for the country
  gni <- country_GNI %>% filter(Country.Code == country_data$ISO3[1]) %>% pull(GNI_per_capita)
  cost_GNI_daily <- gni / 365
  
  # DALY Calculations and Parameters
  cost_daly <- daly_master %>% filter(GID_0 == country_data$ISO3[1]) %>% pull(daly_value)
  
  # GDP Loss per child for the country
  gdp_loss_per_child <- gdp_loss %>% filter(ISO3 == country_data$ISO3[1]) %>% pull(GDP_Loss_Per_Child)
  
  prob_catastrophic <- poverty_master %>% filter(GID_0 == country_data$ISO3[1]) %>% pull(prop_at_risk_catastrp_HCexp)
  
  prob_impoverishment <- poverty_master %>% filter(GID_0 == country_data$ISO3[1]) %>% pull(prop_at_risk_impoverishing)
  
  # Construct the file name
  param_file <- paste0("params_montecarlo_", country_data$ISO3[1], ".Rdata")
  
  # Check if the file exists before loading it
  if (!file.exists(param_file)) {
    message(paste("File not found for country:", country, "ISO3:", country_data$ISO3[1]))
    return(NULL)
  }
  
  # Load specific parameter file for the country
  load(param_file)
  
  results_all <- data.frame()
  
  # Define scenarios
  scenarios <- c("VE1", "VE2", "VE3")
  
  # Extract mean, 2.5% CI, and 97.5% CI for each parameter
  mean_params <- as.list(colMeans(df_params_montecarlo))
  lower_CI <- as.list(apply(df_params_montecarlo, 2, quantile, probs = 0.025))
  upper_CI <- as.list(apply(df_params_montecarlo, 2, quantile, probs = 0.975))
  
  # List of parameters to vary
  params_to_vary <- names(mean_params)
  
  # Loop over each parameter to vary it while keeping others constant
  for (param in params_to_vary) {
    for (CI_label in c("mean", "2.5CI", "97.5CI")) {
      
      # Initialize parameters to mean values
      params <- mean_params
      
      # Vary only the selected parameter
      if (CI_label == "mean") {
        params[[param]] <- mean_params[[param]]
      } else if (CI_label == "2.5CI") {
        params[[param]] <- lower_CI[[param]]
      } else if (CI_label == "97.5CI") {
        params[[param]] <- upper_CI[[param]]
      }
      
      for (scenario in scenarios) {
        # Initialise results dataframe for this country and scenario
        results_sim <- data.frame(
          year = country_data$year,
          DALYs = numeric(length(country_data$year)),
          Societal_cost_disc = numeric(length(country_data$year))
        )
        
        # Calculate results for this scenario
        scenario_results <- calculate_vaccination_scenario(
          country_data,
          scenario,
          results_sim,
          params,
          VSL,
          cost_daly,
          cost_GNI_daily,
          gdp_loss_per_child,
          prob_catastrophic,
          prob_impoverishment
        )
        
        # Store the results with country, scenario, parameter, and CI level information
        results_all <- rbind(
          results_all,
          data.frame(
            country = country,
            ISO3 = country_data$ISO3[1],
            scenario = scenario,
            parameter = param,
            CI_level = CI_label,
            DALYs = mean(scenario_results$DALY_total, na.rm = TRUE),
            Societal_cost_disc = mean(scenario_results$Cost_societal_disc, na.rm = TRUE)
          )
        )
      }
    }
  }
  
  return(results_all)
}

# Process all countries
countries <- unique(extracted_data$ISO3)

# Run the process_country_sens function for each country sequentially
results_sens <- lapply(countries, process_country_sens)

# Combine the results from all countries into a single dataframe
results_sens <- do.call(rbind, results_sens)

# Now we want to sum our results for the whol region:
summed_results <- results_sens %>%
  group_by(scenario, parameter, CI_level) %>%
  summarise(
    Total_DALYs = sum(DALYs, na.rm = TRUE),
    Total_Societal_cost_disc = sum(Societal_cost_disc, na.rm = TRUE)
  ) %>%
  ungroup()

# Then, we calculate Changes Relative to the Mean. We convert to millions to avoid false precision and to make the plots easier to read
tornado_data <- summed_results %>%
  pivot_wider(names_from = CI_level, values_from = c(Total_DALYs, Total_Societal_cost_disc)) %>%
  mutate(
    Change_DALYs_2.5 = (Total_DALYs_2.5CI - Total_DALYs_mean) / 1e6,  # Convert to millions
    Change_DALYs_97.5 = (Total_DALYs_97.5CI - Total_DALYs_mean) / 1e6,  # Convert to millions
    Change_Cost_2.5 = (Total_Societal_cost_disc_2.5CI - Total_Societal_cost_disc_mean) / 1e6,  # Convert to millions
    Change_Cost_97.5 = (Total_Societal_cost_disc_97.5CI - Total_Societal_cost_disc_mean) / 1e6  # Convert to millions
  )

# Apply a more readable parameter mapping with British English spelling
parameter_mapping <- c(
  "oop" = "OOP % CHE",
  "diagnostic_test" = "Diagnostic Test Costs",
  "drug_cost" = "Treatment Costs",
  "hospital_cost" = "Hospitalisation Cost",
  "diagnostic_test_res" = "Diagnostic Test Cost (Resistant Case)",
  "drug_cost_severe" = "Treatment Costs (Severe Case)",
  "hospital_cost_res" = "Hospitalisation Cost (Resistant Case)",
  "dur_ill_prehosp" = "Duration of Illness Before Hospital",
  "dur_hosp_clinical_survived" = "Duration of Hospitalisation (Clinical, Survived)",
  "dur_hosp_resistant_survived" = "Duration of Hospitalisation (Resistant, Survived)",
  "dur_hosp_clinical_survived_severe" = "Duration of Hospitalisation (Severe, Survived)",
  "dur_hosp_died" = "Duration of Hospitalisation (Died)",
  "grief" = "Grief Period",
  "disutility_mod" = "Disability Weighting - Moderate",
  "disutility_severe" = "Disability Weighting - Severe",
  "avg_years_left_to_live" = "Average Years Left to Live",
  "hospital_cost_severe" = "Hospitalisation Cost (Severe)"
)

# Update the tornado data with the new parameter names
tornado_data <- tornado_data %>%
  mutate(parameter = recode(parameter, !!!parameter_mapping))

# Step 3: Create Tornado Plots for DALYs
for (scenario_name in unique(tornado_data$scenario)) {
  scenario_data <- tornado_data %>%
    filter(scenario == scenario_name) %>%
    select(parameter, Change_DALYs_2.5, Change_DALYs_97.5) %>%
    arrange(desc(abs(Change_DALYs_97.5)))  # Arrange by the magnitude of change
  
  ggplot(scenario_data) +
    geom_bar(aes(x = reorder(parameter, abs(Change_DALYs_97.5)), y = Change_DALYs_2.5), stat = "identity", fill = "lightblue", width = 0.7) +
    geom_bar(aes(x = reorder(parameter, abs(Change_DALYs_97.5)), y = Change_DALYs_97.5), stat = "identity", fill = "blue", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = scales::comma) +  # Format y-axis labels with commas
    labs(title = paste("Tornado Plot for DALYs (in millions) -", scenario_name),
         x = "Parameter",
         y = "Change in DALYs (Millions)") +
    theme_minimal()
  
  ggsave(filename = file.path(output_dir_graphs, paste0("Tornado_Plot_DALYs_", scenario_name, ".png")))
}

# Step 4: Create Tornado Plots for Societal Cost
for (scenario_name in unique(tornado_data$scenario)) {
  scenario_data <- tornado_data %>%
    filter(scenario == scenario_name) %>%
    select(parameter, Change_Cost_2.5, Change_Cost_97.5) %>%
    arrange(desc(abs(Change_Cost_97.5)))  # Arrange by the magnitude of change
  
  ggplot(scenario_data) +
    geom_bar(aes(x = reorder(parameter, abs(Change_Cost_97.5)), y = Change_Cost_2.5), stat = "identity", fill = "lightblue", width = 0.7) +
    geom_bar(aes(x = reorder(parameter, abs(Change_Cost_97.5)), y = Change_Cost_97.5), stat = "identity", fill = "blue", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = scales::comma) +  # Format y-axis labels with commas
    labs(title = paste("Tornado Plot for Societal Cost (in millions) -", scenario_name),
         x = "Parameter",
         y = "Change in Societal Cost (Millions)") +
    theme_minimal()
  
  ggsave(filename = file.path(output_dir_graphs, paste0("Tornado_Plot_Cost_", scenario_name, ".png")))
}
