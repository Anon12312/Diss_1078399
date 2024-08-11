#### SECTION 1 - SETUP AND LOAD DATA ####

# Clear the R environment by removing all objects
rm(list=ls())  

# Set the working directory to the specified path
setwd("/Users/X/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")  

# Configure options for displaying numbers: disable scientific notation and set the number of digits to display
options(scipen=999)
options(digits = 3)

# Load necessary libraries for data manipulation and reading Excel files
library(dplyr)
library(tidyr)
library(readxl)
library(purrr)

# Load population data for African countries from an Excel file
WorldPop_Africa <- read_excel("Population_Africa_WorldPopADM1.xlsx")

# Create a data frame that maps country names to their corresponding ISO3 codes
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
  GID_0 = c("AGO", "BEN", "BWA", "BFA", "BDI", 
            "CMR", "CAF", "TCD", "COM", "COG", 
            "CIV", "COD", "GNQ", "ERI", "SWZ", 
            "ETH", "GAB", "GMB", "GHA", "GIN", 
            "GNB", "KEN", "LBR", "MDG", "MWI", 
            "MLI", "MRT", "MOZ", "NAM", "NER", 
            "NGA", "RWA", "STP", "SEN", "SLE", 
            "ZAF", "SSD", "TGO", "UGA", "TZA", 
            "ZMB", "ZWE"))

# Load poverty data based on the threshold of $2.15 USD per day (2017 PPP) from a CSV file
below_2_15USD_pD_2017PPP <- read.csv("poverty line.csv")

# Filter the poverty data to include only the countries of interest
filtered_poverty_data <- below_2_15USD_pD_2017PPP %>%
  inner_join(country_iso3_mapping, by = c("Country.Code" = "GID_0"))

# Reshape the filtered poverty data to a long format, making the year variable explicit
long_format <- filtered_poverty_data %>%
  gather(key = "Year", value = "Poverty_Rate", X2015:X2023) %>%
  mutate(Year = as.integer(sub("X", "", Year)))  # Convert the year variable to integer format

# Filter the data to retain the most recent non-missing value for each country
poverty_master <- long_format %>%
  filter(!is.na(Poverty_Rate)) %>%
  group_by(Country.Name, Country.Code) %>%
  arrange(desc(Year)) %>%
  slice(1) %>%
  ungroup()

# Identify countries that are missing from the poverty data
missing_countries <- country_iso3_mapping %>%
  anti_join(poverty_master, by = c("GID_0" = "Country.Code"))

# Display the list of missing countries
print(missing_countries)

# Rename the Poverty_Rate variable to be more descriptive
poverty_master <- poverty_master %>%
  rename(percentage_below_2_15USD = Poverty_Rate)

# Calculate the median poverty rate across the available data
median_poverty_rate <- median(poverty_master$percentage_below_2_15USD, na.rm = TRUE)

# Manually add missing countries with specific values for the poverty rate
additional_data <- data.frame(
  Country.Name = c("Eritrea", "Comoros", "Congo", "Equatorial Guinea", "Madagascar", "South Africa"),
  Country.Code = c("ERI", "COM", "COG", "GNQ", "MDG", "ZAF"),
  Year = rep(2023, 6),  # Assuming the most recent year for consistency
  percentage_below_2_15USD = c(median_poverty_rate, 18.6, 35.4, 76.8, 80.7, 20.5)
)

# Combine the additional data with the existing poverty data and select relevant columns
poverty_master <- bind_rows(poverty_master, additional_data) %>%
  select(Country.Name, Country.Code, percentage_below_2_15USD) 

#### SECTION 2 - ESTIMATING CATASTROPHIC EXPENDITURE ####

# Background
# Study: PATRONAGE AND COST OF MALARIA TREATMENT IN PRIVATE HOSPITALS IN IBADAN NORTH L.G.A SOUTH WESTERN, NIGERIA  
# (Salawu et al., 2016)
# The cost estimates provided in the study are in local currency for the year 2016.
# The average total malaria treatment cost for children was N 10,371.49.

# Method for Cost Scaling to derive country-specific cost estimates
# This method is based on the tradable resources approach (see Turner et al., 2019; Value in Health, 22(9):1026–1032)
# Steps involved:
# 1. Take the original value expressed in local currency
# --> Average child malaria treatment cost: N 10,371.49

# 2. Convert to international dollars (I$) using the exchange rate at the time of the cost estimation (2016)
# Note: World Bank data was used (PPP_conv_factor.xlsx)
PPP_conv_factor <- read_excel("PPP_conv_factor.xlsx")
names(PPP_conv_factor)
PPP_conv_factor %>% filter(GID_0 == "NGA")
# --> PPP I$ conversion factor for Nigeria in 2016: 105.373917
# --> Calculation: 10,371.49 / 105.373917 = 98.4255 I$

# 3. Inflate using the US dollar inflation rate (divide GDP deflator)
# Note: World Bank data was used (GDP_deflator.xlsx)
GDP_deflator <- read_excel("GDP_deflator.xlsx")
names(GDP_deflator)
GDP_deflator %>% filter(GID_0 == "USA")
# --> GDP_2021 / GDP_2016 = 113.568894993064 / 101.002235480218 = 1.12442
# --> Calculation: 98.4255 I$ * 1.12442 = $I 110.671702
treat_cost_2021IntD <- 110.671702

# 4. Information on the share of out-of-pocket (OOP) expenditure
OOP_proportion <- read.csv("Country OOP 2021.csv") %>%
  select(Country, X2021) %>%
  left_join(., country_iso3_mapping, by = "Country") %>%
  select(Country,GID_0,X2021)
names(OOP_proportion)
OOP_proportion %>% filter(GID_0 == "NGA")
# Proportion of per capita healthcare expenditure paid OOP I$ PPP in 2021: 76.24337 (latest available estimate)

# 5. OOP expenditure adjustment factor:
# Calculation: 110.671602 / 0.7624337 = 145.1558372
treat_cost_2021IntD <- 145.1558372

# Example: Translating Nigeria study costs to Benin setting:
OOP_proportion %>% filter(GID_0 == "BEN")
# = 48.63362
# Calculation: Total Nigeria study costs in I$ * Benin OOP proportion
# 145.1558372 * (0.4863362) = 70.5945 I$

# The following section calculates catastrophic healthcare expenditure thresholds for various countries.
# These thresholds represent the proportion of households that would face financial catastrophe due to healthcare costs.
catastrophic_threshold_data <- data.frame(
  Country.Name = c("Angola", "Burundi", "Benin", "Burkina Faso", "Botswana",
                   "Central African Republic", "Cote d'Ivoire", "Cameroon", "Congo, Dem. Rep.",
                   "Congo", "Comoros", "Eritrea", "Ethiopia", "Gabon", "Ghana", "Guinea",
                   "Equatorial Guinea", "Gambia, The", "Guinea-Bissau", "Kenya", "Liberia", "Madagascar",
                   "Mali", "Mozambique", "Mauritania", "Malawi", "Namibia", "Niger", "Nigeria",
                   "Rwanda", "Senegal", "Sierra Leone", "South Sudan", "Sao Tome and Principe",
                   "Eswatini", "Chad", "Togo", "Tanzania", "Uganda", "South Africa", "Zambia", "Zimbabwe"),
  percentage_below_catastrophic_threshold = c(8.14, 15.63, 9.64, 6.14, 0, 61.22, 1.06, 31.98, 64.6,
                                              9.17, 19.55, 38.9, 11.14, 0.2, 7.99, 13.8, 85, 0.57,
                                              32.16, 2.46, 2.65, 58.95, 4.06, 4.65, 1.4, 3.25, 0,
                                              31.3, 52.08, 0.34, 5.86, 23.45, 32.91, 0.12, 0.33,
                                              36.47, 40.46, 5.22, 13.7, 0.02, 2.45, 0.09)
)

# Combine the catastrophic threshold data with the existing poverty data
poverty_master <- poverty_master %>%
  left_join(catastrophic_threshold_data, by = "Country.Name")

# Calculate the proportion of households at risk of catastrophic healthcare expenditure and impoverishment
poverty_master <- poverty_master %>% mutate(prop_at_risk_catastrp_HCexp = percentage_below_catastrophic_threshold)

poverty_master <- poverty_master %>%
  mutate(prop_at_risk_impoverishing = pmax(0, percentage_below_catastrophic_threshold - percentage_below_2_15USD)) %>%
  mutate(GID_0 = Country.Code) %>%
  mutate(NAME_0 = Country.Name)

# Reorder and select relevant columns for final output
names(poverty_master)
poverty_master <- poverty_master[, c("GID_0", "NAME_0", "prop_at_risk_catastrp_HCexp", "prop_at_risk_impoverishing")]

# Save the final dataset to a file
save(poverty_master, file='poverty_master.rda')  
