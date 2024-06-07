###################Extracting outcome files###########################################
#This R script performed cancer outcome files extraction task from the raw UK Biobank data.
#df1 was Diagnoses - ICD10(field ID 41270) and df2 was Date of first in-patient diagnosis-ICD10(field ID 41280).
#The data tables were converted into character matrices (df1_matrix and df2_matrix). 
#This conversion was done to facilitate faster processing and compatibility with subsequent C++ functions.
#The use of C++ parallel processing, efficient data handling, and integration with C++ ensured that 
#the script could handle large datasets and perform complex computations quickly.
#Finally R loop was used through each code.For each ICD code in the list, the script called the 
#findDatesWithCode function to find dates associated with the current code and converted the 
#result to a data table, added sample IDs, and renamed columns for clarity.

#by: Luo wenjin
#date: 2024-01-12
#######################################################################################


library(data.table)
library(dplyr)
library(Rcpp)
library(RcppParallel)

# Remove all objects from the current environment to start fresh
rm(list=ls())

# Load necessary libraries
library(data.table) 
library(Rcpp) 

# Read two CSV files into data tables
df1 <- fread('~/data/UKB45050/icd_diagnosis/icd_diagnose_2.csv')
df2 <- fread('~/data/UKB45050/icd_diagnosis/icd_diagns_date_2.csv')

# Convert data tables to character matrices
df1_matrix <- as.matrix(df1, mode = "character")
df2_matrix <- as.matrix(df2, mode = "character")

# Read the first column from another CSV file into a vector
sample_ids <- fread("~/Documents/pso_ca/pso_time.csv", select = 1)[[1]]

# Remove the original data tables to free up memory
rm(df1, df2)

# Set options for parallel processing
setThreadOptions(numThreads = 64, stackSize = "auto")

# Source a C++ file containing a function definition for 'findDatesWithCode'
Rcpp::sourceCpp("~/mambaforge/envs/R403/lib/R/library/RcppParallel/include/rcpp_outcome.cpp")

# Set the working directory
setwd("~/data/PSO/cox/ca/")

# Read a file containing a list of codes
code <- fread('~/data/UKB45050/icd_diagnosis/ca_code.txt', header = F)

# Loop through each code in the list
for (i in code$V1) {
  # Call the C++ function to find dates with the current code
  result_matrix <- findDatesWithCode(df1_matrix, df2_matrix, i) %>%
    data.table %>%
    mutate('V1' = sample_ids) %>%
    rename(ID = V1, time = V2)
  
  # Write the result to a CSV file named after the current code
  fwrite(result_matrix, paste0(i, ".csv"), row.names = F)
}
