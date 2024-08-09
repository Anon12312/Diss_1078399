#### SECTION 1 - SET UP AND LOAD DATA ####

# Clear the workspace to ensure a clean environment
rm(list = ls())

# Load the necessary libraries for data manipulation, plotting, and parallel processing
library(dplyr)
library(ggplot2)
library(parallel)
library(scales)
library(pbmcapply)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# Set the working directory to the appropriate location and load the required datasets
setwd("/Users/harrystreet/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")
load("all_malaria.rda")
load("daly_master.rda")
load("poverty_master.rda")
load("vsl_master.rda")

# Set the number of simulations to be used in the Monte Carlo model
n_sims_montecarlo <- 2000

# Load country-specific data on Gross National Income (GNI) per capita and children's population numbers
country_GNI <- read.csv("Country GNI Per Capita 2021 USD.csv") %>%
  select(Country.Code, X2021) %>%
  rename(GNI_per_capita = X2021)

population_data <- read.csv("Pop1_byCountry.csv") # Children's population numbers from Hamilton et al. (2023) for consistency

# Define the output directory to save results and create it if it does not exist
output_dir <- "/Users/harrystreet/Desktop/University of Oxford/3. Trinity Term (Placement)/11. Output Graphs/Worst Case Scenario"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#### SECTION 2 - SETTING FUNCTIONS ####

# Define the discount rate to be applied and create a function for discounting values over time
discRate <- 0.03
discount <- function(value, rate, year_k, base_year = 2021) {
  return(value * (1 / (1 + rate) ^ (year_k - base_year)))
}

# Create a function to calculate costs for a specific scenario; this function focuses on the worst-case scenario by only considering costs
calculate_results <- function(data, scenario, results, params, cost_GNI_daily) {
  
  # Calculate the proportion of severe cases
  prop_severe <- 0.138
  
  # Calculate costs for mild and severe resistant cases
  resistant_cost_mild <- params$diagnostic_test + params$drug_cost + params$hospital_cost * params$dur_hosp_resistant_survived
  resistant_cost_severe <- params$diagnostic_test + params$drug_cost_severe + params$hospital_cost_severe * params$dur_hosp_resistant_survived
  
  for (i in seq_along(data$year)) {
    year_k <- data$year[i]
    
    # Split resistant cases averted into mild and severe categories
    resistant_cases_averted_mild <- data[[paste0("resistant_cases_averted_Ires2_", scenario)]][i] * (1 - prop_severe)
    resistant_cases_averted_severe <- data[[paste0("resistant_cases_averted_Ires2_", scenario)]][i] * prop_severe
    
    # Calculate Disability-Adjusted Life Years (DALYs) for mild and severe resistant cases
    DALY_resistant_mild <- resistant_cases_averted_mild * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_mod / 365
    DALY_resistant_severe <- resistant_cases_averted_severe * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived) * params$disutility_severe / 365
    DALY_resistant <- DALY_resistant_mild + DALY_resistant_severe
    
    DALY_total_disc <- DALY_resistant # Only considering resistant cases for DALY calculation
    
    # Calculate various costs, including governmental and out-of-pocket (OOP) costs for mild and severe resistant cases
    Cost_hosp_resistant_gvt_mild <- resistant_cases_averted_mild * resistant_cost_mild * (1 - params$oop)
    Cost_hosp_resistant_gvt_severe <- resistant_cases_averted_severe * resistant_cost_severe * (1 - params$oop)
    Cost_hosp_resistant_gvt <- Cost_hosp_resistant_gvt_mild + Cost_hosp_resistant_gvt_severe
    
    Cost_hosp_resistant_oop_mild <- resistant_cases_averted_mild * resistant_cost_mild * params$oop
    Cost_hosp_resistant_oop_severe <- resistant_cases_averted_severe * resistant_cost_severe * params$oop
    Cost_hosp_resistant_oop <- Cost_hosp_resistant_oop_mild + Cost_hosp_resistant_oop_severe
    
    Cost_prod_resistant_hospital_mild <- resistant_cases_averted_mild * cost_GNI_daily[i] * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived)
    Cost_prod_resistant_hospital_severe <- resistant_cases_averted_severe * cost_GNI_daily[i] * (params$dur_ill_prehosp + params$dur_hosp_resistant_survived)
    Cost_prod_resistant_hospital <- Cost_prod_resistant_hospital_mild + Cost_prod_resistant_hospital_severe
    
    # Apply discounting to the costs
    Cost_prod_total_disc <- discount(Cost_prod_resistant_hospital, discRate, year_k)
    
    Cost_care_gvt_disc <- discount(Cost_hosp_resistant_gvt, discRate, year_k)
    Cost_care_oop_disc <- discount(Cost_hosp_resistant_oop, discRate, year_k)
    
    Cost_societal_disc <- Cost_care_gvt_disc + Cost_care_oop_disc + Cost_prod_total_disc
    
    # Store the results in the respective columns
    results$DALY_resistant[i] <- DALY_resistant
    results$DALY_total_disc[i] <- DALY_total_disc
    results$Cost_hosp_resistant_gvt[i] <- Cost_hosp_resistant_gvt
    results$Cost_hosp_resistant_oop[i] <- Cost_hosp_resistant_oop
    results$Cost_prod_resistant_hospital[i] <- Cost_prod_resistant_hospital
    results$Cost_prod_total_disc[i] <- Cost_prod_total_disc
    results$Cost_care_gvt_disc[i] <- Cost_care_gvt_disc
    results$Cost_care_oop_disc[i] <- Cost_care_oop_disc
    results$Cost_societal_disc[i] <- Cost_societal_disc
  }
  
  return(results)
}


#### SECTION 3 - PROCESS RESULTS ####

# Create a function to process the simulation for each country, applying country-specific parameters such as GNI per capita
process_country <- function(country) {
  country_data <- extracted_data %>% filter(country == !!country)
  
  # Extract GNI data for the specific country
  gni <- country_GNI %>% filter(Country.Code == country_data$ISO3[1]) %>% pull(GNI_per_capita)
  cost_GNI_daily <- gni / 365
  
  # Load the Monte Carlo parameters specific to the country
  load(paste0("params_montecarlo_", country_data$ISO3[1], ".Rdata"))
  
  all_results <- data.frame()
  
  for (scenario in c("VE1")) { # Adjust for the worst-case scenario
    for (sim in 1:n_sims_montecarlo) {
      # Extract the parameters for the current simulation
      params <- df_params_montecarlo[sim, ]
      
      # Initialise a dataframe to store the results of the simulation
      results_sim <- data.frame(
        year = country_data$year,
        DALY_resistant = numeric(length(country_data$year)),
        DALY_total_disc = numeric(length(country_data$year)),
        Cost_hosp_resistant_gvt = numeric(length(country_data$year)),
        Cost_hosp_resistant_oop = numeric(length(country_data$year)),
        Cost_prod_resistant_hospital = numeric(length(country_data$year)),
        Cost_prod_total_disc = numeric(length(country_data$year)),
        Cost_care_gvt_disc = numeric(length(country_data$year)),
        Cost_care_oop_disc = numeric(length(country_data$year)),
        Cost_societal_disc = numeric(length(country_data$year))
      )
      
      # Calculate results for the current simulation and scenario
      scenario_results <- calculate_results(country_data, scenario, results_sim, params, cost_GNI_daily)
      
      # Combine results with country and scenario information
      all_results <- rbind(all_results, cbind(country = country, ISO3 = country_data$ISO3[1], scenario = scenario, scenario_results))
    }
  }
  
  return(all_results)
}

# Generate a list of countries to process
countries <- unique(extracted_data$country)

# Detect the number of available processor cores and use them all for parallel processing
num_cores <- detectCores()

# Process the simulation results for all countries in parallel
results <- pbmclapply(countries, process_country, mc.cores = num_cores)

#### SECTION 4 - RESULTS FORMATTING AND ANALYSIS ####

# Combine the results from all countries into a single dataframe
all_results <- do.call(rbind, results)

# Summarise and analyse the results, particularly focusing on resistant cases averted
summary_resistant_ve1 <- all_results %>%
  group_by(country, ISO3) %>%
  summarise(
    Mean_Resistant_Cases_Averted = mean(Cost_prod_resistant_hospital, na.rm = TRUE),
    Lower_CI_Resistant_Cases_Averted = quantile(Cost_prod_resistant_hospital, 0.025, na.rm = TRUE),
    Upper_CI_Resistant_Cases_Averted = quantile(Cost_prod_resistant_hospital, 0.975, na.rm = TRUE),
    Mean_OOP_Averted = mean(Cost_hosp_resistant_oop, na.rm = TRUE),
    Lower_CI_OOP_Averted = quantile(Cost_hosp_resistant_oop, 0.025, na.rm = TRUE),
    Upper_CI_OOP_Averted = quantile(Cost_hosp_resistant_oop, 0.975, na.rm = TRUE)
  )

# Save the summarised results as a CSV file
write.csv(summary_resistant_ve1, file.path(output_dir, "summary_resistant_scen4.csv"), row.names = FALSE)

# Calculate out-of-pocket (OOP) averted spending on resistant cases per 1000 children
summary_resistant_ve1 <- summary_resistant_ve1 %>%
  left_join(population_data %>% select(Country.Code, Pop1), by = c("ISO3" = "Country.Code")) %>%
  mutate(
    OOP_Averted_per_1000 = (Mean_OOP_Averted / Pop1) * 1000
  )

# Load map data for Africa
world <- ne_countries(scale = "medium", returnclass = "sf")
africa <- world %>% filter(continent == "Africa")

# Merge the map data with the summary results
africa_data_ve1_resistant <- merge(africa, summary_resistant_ve1, by.x = "iso_a3", by.y = "ISO3", all.x = TRUE)

# Plot the results on the map, focusing on resistant cases averted under the VE1 scenario
p_resistant_ve1 <- ggplot(africa_data_ve1_resistant) +
  geom_sf(aes(fill = Mean_Resistant_Cases_Averted)) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey", name = "Resistant Cases Averted",
                       labels = scales::label_number(scale = 1e-3, suffix = "k")) +
  labs(title = "Resistant Cases Averted for VE1 Scenario by Country",
       fill = "Resistant Cases Averted") +
  theme_minimal() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering)

# Save the plot to the output directory
ggsave(filename = file.path(output_dir, "Resistant_Cases_Averted_VE1.png"), plot = p_resistant_ve1)

# Plot the OOP averted spending per 1000 children for the VE1 scenario
p_oop_averted_per_1000 <- ggplot(africa_data_ve1_resistant) +
  geom_sf(aes(fill = OOP_Averted_per_1000)) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey", name = "OOP Averted per 1000",
                       limits = c(), labels = scales::label_number()) +
  labs(title = "Resistant OOP Averted Spending per 1000 Children for Scenario Four",
       fill = "OOP Averted per 1000") +
  theme_minimal() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering)

# Save the plot to the output directory
ggsave(filename = file.path(output_dir, "OOP_Averted_per_1000_VE1.png"), plot = p_oop_averted_per_1000)

# Plot the total OOP averted spending under the VE1 scenario
p_total_oop_averted <- ggplot(africa_data_ve1_resistant) +
  geom_sf(aes(fill = Mean_OOP_Averted)) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey", name = "Total OOP Averted",
                       labels = scales::label_dollar(scale = 1e-6, suffix = "M"), limits = c()) +
  labs(title = "Total Resistant OOP Averted Spending for Scenario Four",
       fill = "Total OOP Averted (USD)") +
  theme_minimal() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering)

# Save the plot to the output directory
ggsave(filename = file.path(output_dir, "Total_OOP_Averted_VE1.png"), plot = p_total_oop_averted)

# Calculate the aggregate OOP averted spending along with 95% uncertainty intervals
aggregate_oop_averted <- summary_resistant_ve1 %>%
  summarise(
    Total_OOP_Averted = sum(Mean_OOP_Averted, na.rm = TRUE),
    Lower_CI_Total_OOP_Averted = sum(Lower_CI_OOP_Averted, na.rm = TRUE),
    Upper_CI_Total_OOP_Averted = sum(Upper_CI_OOP_Averted, na.rm = TRUE)
  )

# Display the total and confidence interval values for OOP averted spending
sum(aggregate_oop_averted$Total_OOP_Averted)
sum(aggregate_oop_averted$Lower_CI_Total_OOP_Averted)
sum(aggregate_oop_averted$Upper_CI_Total_OOP_Averted)

