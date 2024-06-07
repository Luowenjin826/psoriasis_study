#############################Two_sample_mr##################################
#This script conducted a two-sample Mendelian Randomization (MR) analysis to explore 
#potential causal relationships between gpsoriasis and various cancer outcomes. The R 
#script leveraged the TwoSampleMR and MRPRESSO packages to process and analyze the data. 

#by: Li ruolin
#date:2024-01-18
############################################################################


############################Shell_script####################################
#download anus and anal canal cancer GWAS summary data from GWAS catalog
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90043001-GCST90044000/GCST90043915/GCST90043915_buildGRCh37.tsv.gz
#download lung cancer GWAS summary data from GWAS catalog
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011812/GCST90011812_buildGRCh37.tsv.gz
#download kidney cancer GWAS summary data from GWAS catalog
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011818/GCST90011818_buildGRCh37.tsv.gz
#download non-Hodgkin's lymphoma GWAS summary data from GWAS catalog
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90042001-GCST90043000/GCST90042739/GCST90042739_buildGRCh37.tsv.gz

# Decompress the downloaded files
gunzip GCST90043915_buildGRCh37.tsv.gz
gunzip GCST90011812_buildGRCh37.tsv.gz
gunzip GCST90011818_buildGRCh37.tsv.gz
gunzip GCST90042739_buildGRCh37.tsv.gz
############################################################################


############################R_script########################################

library(TwoSampleMR)
library(MRPRESSO)

# Obtain exposure data
exposure_dat <- extract_instruments("finn-b-L12_PSORIASIS",
                                    p1 = 5e-08,
                                    r2 = 0.001,
                                    kb = 10000)

# Extract outcome data for various cancers
# Skin cancer
outcome_dat_C44 <- extract_outcome_data(snps = exposure_dat$SNP, 
                                        outcomes = "ieu-b-4959")

# Breast cancer
outcome_dat_C50 <- extract_outcome_data(snps = exposure_dat$SNP, 
                                        outcomes = "ebi-a-GCST004988")

# Secondary neoplasm of other sites
outcome_dat_C79 <- extract_outcome_data(snps = exposure_dat$SNP, 
                                        outcomes = "ukb-d-C79")

# Anus cancer
outcome_dat_C21 <- read_outcome_data(
  filename = "GCST90043915_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  chr_col = 'chromosome',
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

# Lung cancer
outcome_dat_C34 <- read_outcome_data(
  filename = "GCST90011812_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  chr_col = 'chromosome',
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

# Kidney cancer
outcome_dat_C64 <- read_outcome_data(
  filename = "GCST90011818_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  chr_col = 'chromosome',
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

# Non-Hodgkin's lymphoma
outcome_dat_C85 <- read_outcome_data(
  filename = "GCST90042739_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  chr_col = 'chromosome',
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

# Perform MR analysis for each outcome data
for (i in list(outcome_dat_C21, outcome_dat_C34, outcome_dat_C64, 
               outcome_dat_C85, outcome_dat_C44, outcome_dat_C50, 
               outcome_dat_C79)) {
  # Harmonise data
  dat <- harmonise_data(exposure_dat, i)
  # Perform MR analysis
  res <- mr(dat)
  print(res)
  # Single SNP analysis
  res_single <- mr_singlesnp(dat)
  mr_forest_plot(res_single)
  mr_scatter_plot(res, dat)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  mr_funnel_plot(res_single)
  
  # MR-PRESSO analysis
  mr_presso(
    BetaOutcome = "beta.outcome", 
    BetaExposure = "beta.exposure", 
    SdOutcome = "se.outcome", 
    SdExposure = "se.exposure", 
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE, 
    data = dat, 
    NbDistribution = 1000,  
    SignifThreshold = 0.05, 
    seed = 6
  )
}
