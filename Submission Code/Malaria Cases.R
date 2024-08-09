#### SECTION 1 - SETUP AND LOAD DATA ####

# Clear the R environment by removing all existing objects
rm(list=ls())

# Load necessary libraries for data manipulation and reading files
library(readxl)
library(dplyr)
library(tidyr)

# Set the working directory to the specified path
setwd("/Users/harrystreet/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")

# Load the ISO3 country codes and malaria model data
iso3 <- read_excel("Africa_ISO3.xlsx")
data <- read.csv("Malaria_MC.csv")

# Merge the ISO3 country codes with the malaria data and set the country column as the index
data <- merge(data, iso3, by = "country", all.x = TRUE)

# Convert all columns (excluding country, ISO3, and year) to represent values per thousand (e.g., cases per thousand)
data <- data %>%
  mutate(across(-c(country, ISO3, year), ~ . * 1000))

# Define the Vaccine Efficacy (VE) scenarios
scenarios_avd <- c("VE1", "VE2", "VE3")


#### SECTION 2 - CONVERT TO AVERTED CASES FOR VACCINE SCENARIOS ####

# Calculate the burden of disease averted for each Vaccine Efficacy (VE) scenario
for (ve in scenarios_avd) {
  data <- data %>%
    mutate(
      !!paste0("cases_averted_", ve) := I_VE0 - get(paste0("I_", ve)),
      !!paste0("cases_averted_", ve, "_min") := I_VE0_min - get(paste0("I_", ve, "_min")),
      !!paste0("cases_averted_", ve, "_max") := I_VE0_max - get(paste0("I_", ve, "_max")),
      !!paste0("resistant_cases_averted_", ve) := Ires_VE0 - get(paste0("Ires_", ve)),
      !!paste0("resistant_cases_averted_", ve, "_min") := Ires_VE0_min - get(paste0("Ires_", ve, "_min")),
      !!paste0("resistant_cases_averted_", ve, "_max") := Ires_VE0_max - get(paste0("Ires_", ve, "_max")),
      !!paste0("deaths_averted_", ve) := D_VE0 - get(paste0("D_", ve)),
      !!paste0("deaths_averted_", ve, "_min") := D_VE0_min - get(paste0("D_", ve, "_min")),
      !!paste0("deaths_averted_", ve, "_max") := D_VE0_max - get(paste0("D_", ve, "_max")),
      !!paste0("resistant_cases_averted_Ires2_", ve) := Ires2_VE0 - Ires2_VE1,
      !!paste0("resistant_cases_averted_Ires2_", ve, "_min") := Ires2_VE0_min - get(paste0("Ires2_", ve, "_min")),
      !!paste0("resistant_cases_averted_Ires2_", ve, "_max") := Ires2_VE0_max - get(paste0("Ires2_", ve, "_max"))
    )
}

# Define the columns to be selected for further analysis
selected_columns <- c("country", "ISO3", "year",
                      "I_VE0", "I_VE0_min", "I_VE0_max",
                      "Ires_VE0", "Ires_VE0_min", "Ires_VE0_max",
                      "D_VE0", "D_VE0_min", "D_VE0_max",
                      "Ires2_VE0", "Ires2_VE0_min", "Ires2_VE0_max")

# Add the newly calculated averted cases columns for each VE scenario to the selected columns
for (ve in scenarios_avd) {
  selected_columns <- c(selected_columns,
                        paste0("cases_averted_", ve), paste0("cases_averted_", ve, "_min"), paste0("cases_averted_", ve, "_max"),
                        paste0("resistant_cases_averted_", ve), paste0("resistant_cases_averted_", ve, "_min"), paste0("resistant_cases_averted_", ve, "_max"),
                        paste0("deaths_averted_", ve), paste0("deaths_averted_", ve, "_min"), paste0("deaths_averted_", ve, "_max"),
                        paste0("resistant_cases_averted_Ires2_", ve), paste0("resistant_cases_averted_Ires2_", ve, "_min"), paste0("resistant_cases_averted_Ires2_", ve, "_max"))
}

# Extract the selected columns from the dataset for further use
extracted_data <- data %>%
  select(all_of(selected_columns))


#### SECTION 3 - RENAME COLUMNS FOR BETTER UNDERSTANDING ####

# Rename columns for clarity, changing VE0 scenarios to 'baseline' and other descriptive labels
colnames(extracted_data) <- gsub("I_VE0", "cases_VE0", colnames(extracted_data))
colnames(extracted_data) <- gsub("I_VE0_min", "cases_VE0_min", colnames(extracted_data))
colnames(extracted_data) <- gsub("I_VE0_max", "cases_VE0_max", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires_VE0", "resistant_cases_VE0", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires_VE0_min", "resistant_cases_VE0_min", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires_VE0_max", "resistant_cases_VE0_max", colnames(extracted_data))
colnames(extracted_data) <- gsub("D_VE0", "deaths_VE0", colnames(extracted_data))
colnames(extracted_data) <- gsub("D_VE0_min", "deaths_VE0_min", colnames(extracted_data))
colnames(extracted_data) <- gsub("D_VE0_max", "deaths_VE0_max", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires2_VE0", "resistant_cases_averted_Ires2_VE0", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires2_VE0_min", "resistant_cases_averted_Ires2_VE0_min", colnames(extracted_data))
colnames(extracted_data) <- gsub("Ires2_VE0_max", "resistant_cases_averted_Ires2_VE0_max", colnames(extracted_data))

# Display the first few rows of the extracted data to verify the changes
print(head(extracted_data))

# Save the extracted data for use in subsequent economic modelling
save(extracted_data, file = "all_malaria.rda")
