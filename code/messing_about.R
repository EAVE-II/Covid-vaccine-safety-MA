library(readxl)


# This is data for England
DaCVaPVaccineSafetyforCR <- read_excel('DaCVaPVaccineSafetyforCR.xlsx')

# haemorrhage data for all countries
haem <- read_excel('haem.xlsx')

# I believe this is the result of the meta analysis for haemorrhage events
ma_res <- readRDS('ma_res.RDS')

# This is data for Wales
Meta_Analysis_Wales_Coeficients <- read.csv('Meta-Analysis_Wales_Coeficients.csv')

# MetaAnalysis.zip contains files that are unzipped into the folder

scotland_ma_results <- read.csv('scotland_ma_results.csv')

