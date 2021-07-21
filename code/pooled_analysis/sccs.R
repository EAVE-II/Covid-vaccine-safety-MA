######################################################################
## Title: Meta-analysis of first dose ChAdOx1 and BNT162b2 COVID-19 vaccinations and thrombocytopenic, venous thromboembolic
##        and haemorrhagic events in the UK
##
## Short title: Vaccine safety meta-analysis
##
## Code author: Steven Kerr
##
## Description: This 'rehydrates' count data from England and Wales with synthetic
#               and combines with Scottish data, generating synthetic IDs
######################################################################

library(readxl)
library(tidyverse)
library(finalfit)
library(survival)

# Set working directory on PHS for SRK
setwd('/conf/EAVE/GPanalysis/progs/SRK/Covid-vaccine-safety-MA')

###################### FUNCTIONS ##################################

rehydrate <- function(df){
  
  output <- data.frame()
  
  for(row in 1:nrow(df)){
    
    #row <-1
    
    period_of_event <- df[row, 'period']
    
    N <- df[row, 'N']
    
    if(period_of_event == 'Pre-risk'){event_col <- c(1,0,0)} else 
      if (period_of_event == 'Clearance'){event_col <- c(0,1,0)}else 
      if (period_of_event == 'Risk'){event_col <- c(0,0,1)}
    
    template <- data.frame( period = c('Pre-risk', 'Clearance', 'Risk'),
                            interval = c(90, 14, 28),
                            event = event_col,
                            vaccine_type = rep(df[row, 'vaccine_type'], 3))
    
    block <- bind_rows(replicate(N, template, simplify = FALSE))
    
    output <- bind_rows(output, block)
  }
  return(output)
}


######################## LOAD DATA ############################

# Values here are taken from data/pooled/analysis/sccs.analysis.html
wales <- data.frame( 'period' = rep( c('Pre-risk', 'Clearance', 'Risk'), 2),
                      vaccine_type = c( rep('AZ', 3), rep('PB', 3)),
                      N =c(6,0,1,1,1,1))

# Values here are taken from data/pooled/analysis/RCGP_CVST_SCCS_eventcounts.xlsx
eng <- data.frame( 'period' = rep( c('Pre-risk', 'Clearance', 'Risk'), 2),
                     vaccine_type = c( rep('AZ', 3), rep('PB', 3)),
                     N =c(3,1,9,4,0,0))


scot <- readRDS('./data/pooled_analyses/scot_sccs_data.rds')

scot <- dplyr::rename(scot, period = expgr, vaccine_type = Vacc.Type) %>%
        select(period, interval, event, vaccine_type)

scot <- mutate(scot, period = gsub("Pre.Vacc","Pre-risk", period))


wales <- rehydrate(wales)

eng <- rehydrate(eng)

df <- rbind(eng, scot)
df<- rbind(df, wales)

# Artificial IDs
df$ID <- ceiling(as.numeric(rownames(df))/3)

event_count <- group_by(df, vaccine_type, period) %>% filter(event == 1) %>% tally %>% pull(n)
  
########################## RESULTS ############################

# Ensure Pre-risk is the baseline level 
df$period <- factor(df$period, levels = c("Pre-risk", "Clearance", "Risk"))

sccs_AZ <- clogit(event ~ period+ strata(ID) + offset(log(interval)),data=df, subset=vaccine_type=="AZ")

sccs_PB <- clogit(event ~ period + strata(ID) + offset(log(interval)),data=df, subset=vaccine_type=="PB")



sccs_results <- data.frame( Period= c('AZ', 'Pre-risk', 'Clearance', 'Risk', 'PB', 'Pre-risk', 'Clearance','Risk'),
                            `Number of events` = c('', event_count[1:3], '', event_count[4:6]),
                            OR = c('', 1,  round(exp(sccs_AZ$coef), 2), '', 1, round(exp(sccs_PB$coef), 2)), 
                            `CI` = '', stringsAsFactors = FALSE)

se <- c(0.6110, 0.3000, 1.05409, 0.79408)

upper <- sprintf('%.2f', exp( c(sccs_AZ$coef, sccs_PB$coef) + 1.96 * se) )
lower <- sprintf('%.2f', exp( c(sccs_AZ$coef, sccs_PB$coef) - 1.96 * se) )

CI <- paste0('[', lower, '-', upper, ']'     )

sccs_results[c(3,4,7,8), 'CI'] <- CI

names(sccs_results)[4] <- '95% CI'

write.csv(sccs_results, './output/pooled_analysis/sccs_cvst.csv', row.names = FALSE)

