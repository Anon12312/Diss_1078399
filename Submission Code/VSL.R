#### SECTION 1 - SETUP AND LOAD DATA ####

# Clear the R environment by removing all existing objects
rm(list = ls())

# Set options to prevent scientific notation and to limit the number of digits displayed
options(scipen=999)
options(digits = 3)

# Load necessary libraries for data manipulation and reading Excel files
library(dplyr)
library(tidyr)
library(readxl)

# Set the working directory to the specified path
setwd("/Users/X/Desktop/University of Oxford/3. Trinity Term (Placement)/5. Data Sources")

# Load the WorldPop Africa dataset containing population data for African regions
WorldPop_Africa <- read_excel("Population_Africa_WorldPopADM1.xlsx")

# Display the number of rows and column names in the dataset for verification
nrow(WorldPop_Africa)
names(WorldPop_Africa)

# Check for and return the number of unique regions (GID_0 codes) in the dataset
length(unique(WorldPop_Africa$GID_0))

# Remove any duplicate entries based on the GID_0 (region code) to ensure each region is represented only once
WorldPop_Africa <- WorldPop_Africa[!duplicated(WorldPop_Africa$GID_0), ]

# Load the Value of Statistical Life (VSL) estimates for African countries
vsl <- read_excel("VSL_for_R.xlsx")

# Merge the WorldPop Africa dataset with the VSL dataset based on the country name (NAME_0)
vsl_master <- WorldPop_Africa %>%
  left_join(., vsl, by="NAME_0")

# Identify and display any regions (GID_0 codes) that have missing VSL estimates in the merged dataset
vsl_master %>%
  filter(is.na(vsl_estimate_2021_IntDoll)) %>%
  filter(!duplicated(GID_0))

# Save the merged dataset (vsl_master) as an R data file for future use
save(vsl_master, file="vsl_master.rda")
