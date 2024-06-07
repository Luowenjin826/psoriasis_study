##########################Psoriasis_and_PRS_phewas##########################
#This R script was designed to perform comprehensive Cox and logistic regression analyses 
#using parallel processing to handle large data efficiently. Initially, it read covariate 
#data and preprocesses various input files, calculating time variables and merging with 
#covariates to prepare the datasets for analysis. The script defined functions for data 
#preprocessing, Cox regression for psoriasis, and logistic regression for PRS, ensuring 
#results were formatted with hazard ratios and confidence intervals. Specific lists of 
#files for different cancer were created to categorize and exclude certain data based on 
#case counts. The script defined multiple regression formulas tailored for various models 
#and subsets of data, including sex-specific and dietary factors. Finally, the script 
#orchestrated the entire analysis process, performing Cox regression and logistic 
#regression for each file, summarizing the results, and writing them to CSV files. 
#This approach leveraged parallel computing to efficiently analyze large datasets, 
#providing detailed summaries of survival and logistic regression results.

#by: Luo wenjin
#date: 2024-01-15
############################################################################


library(data.table)
library(plyr)
library(dplyr)
library(survival)
library(future)
library(future.apply)


##########################Psoriasis_Cox_regression##########################
# Remove all objects from the current environment to start fresh
rm(list=ls())
gc()
# Load necessary libraries
library(data.table) 
library(dplyr) 
library(survival) 
library(furrr) 
library(plyr) 

# Read the covariates data
covar <- fread('~/data/PSO/cox/covar.csv')

# Function to read and preprocess input data
input_dat <- function(file_name) {
  # Read the input file and rename the 'time' column to 'ca_date'
  dt1 <- fread(file_name) %>%
    dplyr::rename(ca_date = time)
  
  # Read the pso data, merge with covariates, and preprocess
  pso <- fread("~/data/PSO/cox/pso.csv") %>%
    left_join(covar, by = 'ID') %>%
    filter(!is.na(attend.0)) %>%
    left_join(dt1, by = 'ID') %>%
    dplyr::rename(pso_date = time) %>%
    mutate(
      'pso_time' = as.numeric(difftime(pso_date, attend.0, units = 'days')) / 365,
      'ca_time' = as.numeric(difftime(ca_date, attend.0, units = 'days')) / 365
    ) %>%
    filter(
      !(!is.na(pso_time) & 
          !is.na(ca_time) & 
          (pso_time > ca_time))
    ) %>%
    mutate('time' = ca_time)
  
  # Adjust 'time' for specific conditions
  ind <- which(
    !is.na(pso$pso_time) & 
      !is.na(pso$ca_time) & 
      (pso$ca_time > pso$pso_time)
  )
  pso$time[ind] <- pso$ca_time[ind] - pso$pso_time[ind]
  
  # Final data processing and outcome definition
  dt2 <- pso %>%
    filter(time > 0 | is.na(time)) %>%
    mutate(
      'pso_outcome' = ifelse(is.na(pso_time), 0, 1),
      'ca_outcome'  = ifelse(is.na(ca_time), 0, 1)
    ) %>%
    data.frame
  
  # Replace NA times with a large number (9999)
  dt2$time[which(is.na(dt2$time))] <- 9999
  
  # Convert columns 11 to 21 to factors
  for (i in 11:21) {
    dt2[, i] <- factor(dt2[, i])
  }
  
  return(dt2)
}

# Function to perform Cox regression and format the results
cox_reg <- function(fml, dat) {
  sm <- summary(coxph(fml, dat)) # Perform Cox regression
  dt1 <- data.frame(
    'hr' = sm$coefficients[1, 2],
    'lower' = sm$conf.int[1, 3],
    'upper' = sm$conf.int[1, 4],
    'pvalue' = sm$coefficients[1, 5]
  )
  # Format the results to 2 decimal places
  for (i in 1:3) {
    dt1[, i] <- sprintf("%0.2f", dt1[, i])
  }
  dt1 <- dt1 %>% mutate('HR' = paste0(hr, "(", lower, ",", upper, ")"))
  return(dt1)
}

# Set up parallel processing with 60 workers
plan(multisession, workers = 64)

# Set the working directory
setwd("~/data/PSO/cox/ca/")

# Process all files in the directory and summarize outcomes
outcome_case <- future_lapply(list.files(), FUN = function(x) {
  dt1 <- input_dat(x) # Preprocess each file
  dt2 <- data.frame(
    'disease' = x,
    'case_ca' = length(which(dt1$ca_outcome == 1)),
    'case_ca_male' = length(which(dt1$ca_outcome == 1 & dt1$sex == 1)),
    'case_ca_female' = length(which(dt1$ca_outcome == 1 & dt1$sex == 0))
  )
  return(dt2)
}) %>% ldply(data.frame) # Combine the results into a single data frame

# Write the outcome summary to a CSV file
fwrite(outcome_case, 
       "~/data/PSO/cox/outcome_case.csv", 
       row.names = F)


# Define lists of file names for different categories of diseases
colorectum <- paste0("C", 18:20, ".csv") # Colorectal cancer files
skin <- paste0("C", 43:44, ".csv") # Skin cancer files
female_genital <- paste0("C", 51:57, ".csv") # Female genital cancer files
male_genital <- paste0("C", 60:62, ".csv") # Male genital cancer files

# Calculate the list of outcome files to be processed
outcome1 <- setdiff(
  list.files(), # List all files in the directory
  c(
    outcome_case$disease[outcome_case$case_ca < 50], # Exclude diseases with less than 50 cases
    colorectum, 
    skin, 
    female_genital, 
    male_genital 
  )
)

# Helper function to process data and perform Cox regression
process_data <- function(file_list, formula, sex_filter = NULL) {
  future_lapply(file_list, FUN = function(x) {
    dat <- input_dat(x)
    if (!is.null(sex_filter)) {
      dat <- dat %>% filter(sex == sex_filter)
    }
    tb <- cox_reg(formula, dat) %>%
      mutate(
        disease = substring(x, 1, (nchar(x) - 4)),
        case_total = paste0(length(which(dat$ca_outcome == 1)), "/", nrow(dat))
      )
    return(tb)
  }) %>% ldply(data.frame)
}

# Define survival analysis formulas for different models using the Surv function from the survival package

# Formula 1: Includes pso_outcome, age, sex, Ethnic, smoking_status, alcohol_frequency, 
# physical_activity, BMI, GC, methotrexate, and cyclosporin.
fml1 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + sex + 
                    Ethnic + smoking_status + alcohol_frequency + 
                    physical_activity + BMI + GC + methotrexate + cyclosporin")

# Formula 2: Similar to fml1 but adds dietary variables red_meat and fiber.
fml2 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + sex + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + red_meat + fiber")

# Formula 3: Similar to fml1 but adds a variable for time spent outdoors.
fml3 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + sex + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + time_outdoors")

# Formula 4: Similar to fml1 but excludes the sex variable.
fml4 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin")

# Formula 5: Similar to fml4 but adds dietary variables red_meat and fiber.
fml5 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + red_meat + fiber")

# Formula 6: Similar to fml4 but adds a variable for time spent outdoors.
fml6 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + time_outdoors")

# Formula 7: Specifically tailored for female-related outcomes, includes menarche_age,
# menopause, and live_births_number.
fml7 <- as.formula("Surv(time, ca_outcome) ~ pso_outcome + age + Ethnic + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + menarche_age + 
                    menopause + live_births_number")

# Define a function to process and analyze data
res <- function(i) {
  # All cases
  res1 <- process_data(outcome1, fml1)
  res2 <- process_data(colorectum, fml2)
  res3 <- process_data(skin, fml3)
  
  # Female cases
  res4 <- process_data(outcome1, fml4, sex_filter = 0)
  res5 <- process_data(colorectum, fml5, sex_filter = 0)
  res6 <- process_data(skin, fml6, sex_filter = 0)
  res7 <- process_data(female_genital, fml7, sex_filter = 0)
  
  # Male cases
  res8 <- process_data(c(outcome1, male_genital), fml4, sex_filter = 1)
  res9 <- process_data(colorectum, fml5, sex_filter = 1)
  res10 <- process_data(skin, fml6, sex_filter = 1)
  
  # Combine and arrange results
  tb1 <- rbind(res1, res2, res3) %>% arrange(disease)
  tb2 <- rbind(res4, res5, res6, res7) %>% arrange(disease)
  tb3 <- rbind(res8, res9, res10) %>% arrange(disease)
  
  return(list(tb1, tb2, tb3))
}
cox_reg_res<-res()
fwrite(cox_reg_res[[1]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/pso_cox_all.csv",
       row.names = F)
fwrite(cox_reg_res[[2]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/pso_cox_female.csv",
       row.names = F)
fwrite(cox_reg_res[[3]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/pso_cox_male.csv",
       row.names = F)
##############################################################################


#########################Estimating_multicollinearity#########################
#VIF
library(car)
dat1<-input_dat('C34.csv')
df1<-vif(coxph(fml1,dat1)) %>%
  data.frame %>%
  mutate('var' = rownames(df1)) %>%
  dplyr::select(4,1)
fwrite(df1,"~/data/PSO/cox/res/vif.csv",row.names = T)
##############################################################################


###########################PRS_logistic_regression############################
rm(list=ls())
gc()
setwd("~/data/PSO/cox/ca")
covar<-fread('~/data/PSO/cox/covar.csv')
#prs logistic regression excluded baseline ca
input_dat<-function(file_name){
  dt1<-fread(file_name) %>%
    dplyr::rename(ca_date = time)
  prs<-fread('~/data/PSO/cox/prs.csv') %>%
    left_join(covar,by='ID') %>%
    filter(!is.na(attend.0)) %>%
    left_join(dt1,by='ID') %>%
    mutate('time' = round(as.numeric(
      difftime(ca_date,attend.0,units = 'days')/365),2)) %>%
    filter(time > 0 | is.na(time)) %>%
    mutate('ca_outcome' = ifelse(is.na(time),0,1)) %>%
    data.frame
  for (i in 19:29){prs[,i]<-factor(prs[,i])}
  return(prs)
}

plan(multisession,workers = 64)
outcome_case<-future_lapply(list.files(),FUN=function(x){
  dt1<-input_dat(x)
  dt2<-data.frame('disease' = x,
                  'case_ca' = length(which(dt1$ca_outcome == 1)),
                  'case_ca_male' = length(which(dt1$ca_outcome == 1 & dt1$sex == 1)),
                  'case_ca_female' = length(which(dt1$ca_outcome == 1 & dt1$sex == 0)))
  return(dt2)
}) %>% ldply(data.frame)
fwrite(outcome_case,
       "~/data/PSO/cox/outcome_case_prs_lr.csv",
       row.names = F)

colorectum<-paste0("C",18:20,".csv")
skin<-paste0("C",43:44,".csv")
female_genital<-paste0("C",51:57,".csv")
male_genital<-paste0("C",60:62,".csv")
outcome1<-setdiff(list.files(),
                  c(outcome_case$disease[outcome_case$case_ca < 50],
                    colorectum,
                    skin,
                    female_genital,
                    male_genital))

lr_reg<-function(fml,dat){
  fit<-glm(fml,family = 'binomial',data = dat)
  beta<-summary(fit)$coefficients[2,1]
  se<-summary(fit)$coefficients[2,2]
  p<-summary(fit)$coefficients[2,4]
  dt1<-data.frame('OR' = exp(beta),
                  'OR_lower' = exp(beta - 1.96*se),
                  'OR_upper' = exp(beta + 1.96*se),
                  'pvalue' = p)
  for (i in 1:3){dt1[,i]<-sprintf("%.2f",dt1[,i])}
  dt1<-dt1 %>% mutate('OR_CI' = paste0(OR,"(",OR_lower,",",OR_upper,")"))
  return(dt1)
}

process_data <- function(file_list, formula, sex_filter = NULL) {
  future_lapply(file_list, FUN = function(x) {
    dat <- input_dat(x) %>%
      filter(Ethnic == 1) %>% # Filter for specific ethnicity
      { if (!is.null(sex_filter)) filter(., sex == sex_filter) else . } # Optional sex filter
    tb <- lr_reg(formula, dat) %>%
      mutate(
        disease = substring(x, 1, (nchar(x) - 4)),
        case_total = paste0(length(which(dat$ca_outcome == 1)), "/", nrow(dat))
      )
    return(tb)
  }) %>% ldply(data.frame)
}

res <- function(i) {
  # All cases
  res1 <- process_data(outcome1, fml1)
  res2 <- process_data(colorectum, fml2)
  res3 <- process_data(skin, fml3)
  
  # Female cases
  res4 <- process_data(outcome1, fml4, sex_filter = 0)
  res5 <- process_data(colorectum, fml5, sex_filter = 0)
  res6 <- process_data(skin, fml6, sex_filter = 0)
  res7 <- process_data(female_genital, fml7, sex_filter = 0)
  
  # Male cases
  res8 <- process_data(c(outcome1, male_genital), fml4, sex_filter = 1)
  res9 <- process_data(colorectum, fml5, sex_filter = 1)
  res10 <- process_data(skin, fml6, sex_filter = 1)
  
  # Combine and arrange results
  tb1 <- rbind(res1, res2, res3) %>% arrange(disease)
  tb2 <- rbind(res4, res5, res6, res7) %>% arrange(disease)
  tb3 <- rbind(res8, res9, res10) %>% arrange(disease)
  
  return(list(tb1, tb2, tb3))
}

# Define logistic regression formulas for different models

# Formula 1: Includes prs, age, sex, smoking_status, alcohol_frequency, physical_activity,
# BMI, GC, methotrexate, cyclosporin, and principal components 1 to 10.
fml1 <- as.formula("ca_outcome ~ prs + age + sex + 
                    smoking_status + alcohol_frequency + 
                    physical_activity + BMI + GC + methotrexate + cyclosporin + 
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 2: Similar to fml1 but adds dietary variables red_meat and fiber.
fml2 <- as.formula("ca_outcome ~ prs + age + sex + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + red_meat + fiber + 
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 3: Similar to fml1 but adds a variable for time spent outdoors.
fml3 <- as.formula("ca_outcome ~ prs + age + sex + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + time_outdoors + 
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 4: Similar to fml1 but excludes the sex variable.
fml4 <- as.formula("ca_outcome ~ prs + age + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + pc1 + pc2 + pc3 + 
                    pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 5: Similar to fml4 but adds dietary variables red_meat and fiber.
fml5 <- as.formula("ca_outcome ~ prs + age + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + red_meat + fiber + 
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 6: Similar to fml4 but adds a variable for time spent outdoors.
fml6 <- as.formula("ca_outcome ~ prs + age + 
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + time_outdoors + 
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Formula 7: Specifically tailored for female-related outcomes, includes menarche_age,
# menopause, and live_births_number.
fml7 <- as.formula("ca_outcome ~ prs + age +  
                    smoking_status + alcohol_frequency + physical_activity + 
                    BMI + GC + methotrexate + cyclosporin + menarche_age + 
                    menopause + live_births_number + pc1 + pc2 + pc3 + pc4 + 
                    pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

lr_reg_res<-res()
fwrite(lr_reg_res[[1]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/prs_lr_all.csv",
       row.names = F)
fwrite(lr_reg_res[[2]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/prs_lr_female.csv",
       row.names = F)
fwrite(lr_reg_res[[3]][,c(6,7,5,4)],
       "~/data/PSO/cox/res/prs_lr_male.csv",
       row.names = F)

