#### SECTION 1 - SETUP AND LOAD DATA ####

# Clear the R environment by removing all existing objects
rm(list=ls())

# Load the necessary libraries for data manipulation and reading Excel files
library(dplyr)
library(tidyr)
library(readxl)

# Set the working directory to the specified path
setwd("/Users/harrystreet/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")

# Load the population data for Africa from the WorldPop dataset
WorldPop_Africa <- read_excel(paste0("Population_Africa_WorldPopADM1.xlsx"))

# Display the number of rows in the dataset
nrow(WorldPop_Africa)

# Display the column names in the dataset
names(WorldPop_Africa)

# Display the number of unique countries (GID_0) in the dataset
length(unique(WorldPop_Africa$GID_0))

# Remove duplicate rows based on the GID_0 (country code) column
WorldPop_Africa <- WorldPop_Africa[!duplicated(WorldPop_Africa$GID_0), ]

# Load the GDP per capita data for 2021
gdp_pc_2021 <- read_excel(paste0("GDP_pC_for_R.xlsx"))

#### SECTION 2 - GENERATE MONETISED DALY DATA ####

# Merge the population data with the GDP per capita data using the GID_0 column as the key
daly_master <- left_join(WorldPop_Africa, gdp_pc_2021, by="GID_0")

# Identify and filter out any countries with missing GDP per capita data
daly_master %>%
  filter(is.na(GDP_pC_2021_currentUS)) %>%
  filter(!duplicated(GID_0))

# Manually impute missing GDP per capita values for Eritrea and South Sudan based on external sources
daly_master <- daly_master %>%
  mutate(GDP_pC_2021_currentUS = ifelse(GID_0 == "ERI", 614.26, GDP_pC_2021_currentUS)) %>%
  mutate(GDP_pC_2021_currentUS = ifelse(GID_0 == "SSD", 749.7, GDP_pC_2021_currentUS))

# Display the column names of the dataset
names(daly_master)

# Remove unnecessary columns from the dataset to simplify further analysis
daly_master$NAME_1 <- NULL
daly_master$Population <- NULL
daly_master$NAME_0.y <- NULL
daly_master$GID_1 <- NULL

# Rename the region variable (NAME_0.x) to a more appropriate name, 'country'
daly_master <- daly_master %>% 
  rename(NAME_0 = NAME_0.x)

# Load the Disability-Adjusted Life Year (DALY) estimates for the countries
daly_estimates <- read_excel(paste0("DALY_for_R.xlsx"))

# Merge the DALY estimates with the existing dataset using the country name as the key
daly_master <- left_join(daly_master, daly_estimates, by="NAME_0")

# Calculate the median percentage of GDP per capita allocated for DALYs across all countries
median_perc <- median(daly_master$perc_of_GDPpc, na.rm = TRUE)

# Impute specific percentages of GDP per capita for DALYs in certain countries based on external studies or the median value
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Angola", median_perc, perc_of_GDPpc))
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Central African Republic", median_perc, perc_of_GDPpc))
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Côte d'Ivoire", 0.19, perc_of_GDPpc)) # Value from Ochalke et al.
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Republic of the Congo", 0.15, perc_of_GDPpc)) 
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Eritrea", 0.27, perc_of_GDPpc)) # Value from Ochalke et al.
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Gambia", 0.69, perc_of_GDPpc))
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Equatorial Guinea", median_perc, perc_of_GDPpc))
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "Liberia", median_perc, perc_of_GDPpc)) # Value from Ochalke et al.
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "South Sudan", 0.16, perc_of_GDPpc)) # Value from Ochalke et al.
daly_master <- daly_master %>% mutate(perc_of_GDPpc = ifelse(NAME_0 == "São Tomé and Príncipe", median_perc, perc_of_GDPpc)) # Value from Ochalke et al.

# Remove the US_2015 column, which is no longer necessary
daly_master$US_2015 <- NULL  

# Calculate the monetised value of DALYs by multiplying GDP per capita by the percentage of GDP per capita allocated for DALYs
daly_master <- daly_master %>%
  mutate(daly_value = GDP_pC_2021_currentUS * perc_of_GDPpc)

# Save the processed data to a file for future use
save(daly_master, file="daly_master.rda")

# Display the final dataset to verify that the data has been processed as expected
daly_master 
