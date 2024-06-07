################################One_sample_mr###############################
#This R script performed one-sample Mendelian Randomization (MR) analysis using logistic 
#regression models, focusing on genetic predispositions to diseases. It definesd a 
#function input_dat2 to preprocess input data files by merging them with covariate 
#data and calculating relevant variables for MR analysis. The main function for MR, one_sample_mr, 
#fited  two-step logistic regression models and calculated the causal effect estimate, 
#its standard error, confidence intervals, and p-value.The script set up formulas for 
#the logistic regression models, considering polygenic risk scores (PRS), age, sex, and 
#principal components to adjust for population stratification. It then identified relevant 
#files for analysis, excluding those with fewer than 50 cases. Using future_lapply, the 
#script processed each file in parallel, standardizing the PRS, fitting the logistic models, 
#and summarizing the results.The analysis is divided into three parts: overall cases, 
#female-specific cases, and male-specific cases. This approach efficiently handled 
#large-scale genetic data to assess the causal effects of genetic variants on disease 
#risk across different populations.

#by: Li ruolin
#date:2024-01-16
############################################################################


library(data.table)
library(dplyr)
library(future)
library(plyr)

# Define a function to preprocess input data
input_dat2 <- function(x) {
  # Load PRS data
  prs <- fread('~/Documents/pso_ca/prs_cov.csv', select = c(1, 2, 5:14)) 
  # Load PSO data
  pso <- fread("~/data/PSO/cox/pso.csv") 
  dt1 <- fread(x) %>%
    dplyr::rename(ca_date = time) %>%
    inner_join(pso, by = 'ID') %>%
    inner_join(prs, by = 'ID') %>%
    left_join(covar, by = 'ID') %>%
    dplyr::rename(pso_date = time) %>%
    mutate(
      pso = ifelse(is.na(pso_date), 0, 1),
      ca = ifelse(is.na(ca_date), 0, 1)
    ) %>%
    filter(!is.na(attend.0) & Ethnic == 1)
  return(dt1)
}

# Define a function to perform one-sample Mendelian Randomization
one_sample_mr <- function(fml1, fml2, data) {
  fit1 <- glm(fml1, family = 'binomial', data = data)
  fit2 <- glm(fml2, family = 'binomial', data = data)
  
  beta1 <- summary(fit1)$coefficients[2, 1]
  beta2 <- summary(fit2)$coefficients[2, 1]
  se1 <- summary(fit1)$coefficients[2, 2]
  se2 <- summary(fit2)$coefficients[2, 2]
  
  beta <- beta2 / beta1
  se <- sqrt(se2^2 / beta1^2 + beta2^2 * se1^2 / beta1^4)
  beta_lower <- beta - 1.96 * se
  beta_upper <- beta + 1.96 * se
  z <- beta / se
  p <- 2 * (1 - pnorm(abs(z)))
  
  dt1 <- data.frame(
    OR = exp(beta),
    OR_lower = exp(beta_lower),
    OR_upper = exp(beta_upper),
    pvalue = p
  )
  
  for (i in 1:3) {
    dt1[, i] <- sprintf("%.2f", dt1[, i])
  }
  
  dt1 <- dt1 %>%
    mutate(OR_CI = paste0(OR, "(", OR_lower, ",", OR_upper, ")"))
  
  return(dt1)
}

# Define formulas for logistic regression
fml1 <- as.formula("pso ~ prs + age + sex + pc1 + pc2 + pc3 + pc4 + 
                    pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
fml2 <- as.formula("ca ~ prs + age + sex + pc1 + pc2 + pc3 + pc4 + 
                    pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
fml3 <- as.formula("pso ~ prs + age + pc1 + pc2 + pc3 + pc4 + 
                    pc5 + pc6 + pc7 + pc8 + pc9 + pc10")
fml4 <- as.formula("ca ~ prs + age + pc1 + pc2 + pc3 + pc4 + 
                    pc5 + pc6 + pc7 + pc8 + pc9 + pc10")

# Load outcome cases and define relevant files
outcome_case <- fread('~/data/NC_PSO/cox/outcome_case.csv')
female_genital <- paste0("C", 51:57, ".csv")
male_genital <- paste0("C", 60:62, ".csv")
outcome2 <- setdiff(
  list.files(),
  c(outcome_case$disease[outcome_case$case_ca < 50], female_genital, male_genital)
)

# Process all relevant outcome files
res1 <- future_lapply(outcome2, FUN = function(x) {
  dat <- input_dat2(x)
  dat$prs <- (dat$prs - mean(dat$prs)) / sd(dat$prs)
  tb <- one_sample_mr(fml1, fml2, dat) %>%
    mutate(
      disease = substring(x, 1, (nchar(x) - 4)),
      case_total = paste0(length(which(dat$ca == 1)), "/", nrow(dat))
    )
  return(tb)
}) %>% ldply(data.frame)

# Process female cases
res2 <- future_lapply(c(outcome2, female_genital), FUN = function(x) {
  dat <- input_dat2(x) %>%
    filter(sex == 0)
  dat$prs <- (dat$prs - mean(dat$prs)) / sd(dat$prs)
  tb <- one_sample_mr(fml3, fml4, dat) %>%
    mutate(
      disease = substring(x, 1, (nchar(x) - 4)),
      case_total = paste0(length(which(dat$ca == 1)), "/", nrow(dat))
    )
  return(tb)
}) %>% ldply(data.frame)

# Process male cases
res3 <- future_lapply(c(outcome2, male_genital), FUN = function(x) {
  dat <- input_dat2(x) %>%
    filter(sex == 1)
  dat$prs <- (dat$prs - mean(dat$prs)) / sd(dat$prs)
  tb <- one_sample_mr(fml3, fml4, dat) %>%
    mutate(
      disease = substring(x, 1, (nchar(x) - 4)),
      case_total = paste0(length(which(dat$ca == 1)), "/", nrow(dat))
    )
  return(tb)
}) %>% ldply(data.frame)

fwrite(res1,
       "~/data/NC_PSO/cox/res/one_sample_mr_all.csv",
       row.names=F)
fwrite(res2,
       "~/data/NC_PSO/cox/res/one_sample_mr_female.csv",
       row.names=F)
fwrite(res3,
       "~/data/NC_PSO/cox/res/one_sample_mr_male.csv",
       row.names=F)
