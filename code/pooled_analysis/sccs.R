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

# Set working directory on PHS for SRK
setwd('/conf/EAVE/GPanalysis/progs/SRK/Covid-vaccine-safety-M')

###################### FUNCTIONS ##################################

rehydrate <- function(df){
  
  output <- data.frame()
  
  for(row in 1:nrow(df)){
    
    #row <-1
    
    print(row)
    
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

# Values here are taken from 
wales <- data.frame( 'period' = rep( c('Pre-risk', 'Clearance', 'Risk'), 2),
                      vaccine_type = rep( c('AZ', 'PB'), 3),
                      N =c(6,1,0,1,1,1)   )



scot <- readRDS('./data/pooled_analyses/scot_sccs_data.rds')

scot <- dplyr::rename(scot, period = expgr, vaccine_type = Vacc.Type) %>%
        select(period, interval, event, vaccine_type)

scot <- select(scot, )


df <- rehydrate(wales)

df <- rbind(df, scot)


clogit(event ~ period + strata(ID) + offset(log(interval)),data=df, subset=vaccine_type=="PB")





