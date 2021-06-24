######################################################################
## Title: Meta-analysis of first dose ChAdOx1 and BNT162b2 COVID-19 vaccinations and thrombocytopenic, venous thromboembolic
##        and haemorrhagic events in the UK
##
## Short title: Vaccine safety meta-analysis
##
## Code author: Chris Robertson, Steven Kerr
##
## Description: This does the meta analysis using estimated odds ratios for adverse events
##              following vaccination from each country
######################################################################

library(readxl)
library(plyr)
library(tidyverse)
library(meta)
library(stringr)

source('./code/01_MA_Input.R')

####################### FUNCTIONS #############################################

# This sets upper and lower confidence interval entries for treatment effects
# to NA if either of their exponentials is infinite.
setNA <- function(upper, lower){
  
  indices <- is.infinite(exp(upper)) | is.infinite(exp(lower))
  
  upper[indices] <- NA
  lower[indices] <- NA
  
  return(list(upper, lower))
}

# This sets upper and lower confidence interval entries for treatment effects
# to 0 if either of their exponentials is infinite.
set0 <- function(upper, lower){
  
  indices <- is.infinite(exp(upper)) | is.infinite(exp(lower))
  
  upper[indices] <- 0
  lower[indices] <- 0
  
  return(list(upper, lower))
}

# Replace time entries with something more readable.
time_replace <- function(vector){
  str_replace_all(vector, c("v1_0:6" = 'Day 0-6', 
                            "v1_7:13" = 'Day 7-13',
                            "v1_14:20" = 'Day 14-20',
                            "v1_21:27" = 'Day 21-27',
                            "v1_0:27" = 'Day 0-27',
                            "v1_28" = 'Day 28'))
}

# Test
# vacc <- 'AZ'
#event <- 'itp'

# Create individual table + meta-analysis for a given vaccine and event
create_table_ma_plot <- function(vacc, event){
  
  df$time <- time_replace(df$time)
  
  names(df) <- gsub('time', 'Time period', names(df))
  
  times <- c("Day 0-6", "Day 7-13", "Day 14-20", "Day 21-27", "Day 0-27","Day 28+") 
  
  df <- filter(df, `Time period` %in% times & vaccine==vacc & group == event)
  
  df <- df[order(df$country),]
  
  table <- select(df, `Time period`, country, N, R) %>% arrange(`Time period`, country) %>% 
    pivot_wider(id_cols= `Time period`, names_from=country, values_from=c(N,R)) %>%
    select(1, 5, 2, 6, 3, 7, 4)
  
  table[is.na(table)] <- 0

  table[, 2:ncol(table)][ table[, 2:ncol(table)] < 5 & table[, 2:ncol(table)] > 0 ] <- '<5'
  
  table <- table[match(times, pull(table, 'Time period')),]
  
  # Format numbers with commas per three decimal places
  table <- mutate_if(table, is.numeric, ~formatC(round(.), format = "f", big.mark = ",", drop0trailing = TRUE) )
  
  new_row <- c(endpoints[[event]], rep('', 8))
  
  table <- rbind(new_row, table)
  
  ma <- metagen(TE=df$log_rr, seTE=df$se_log_rr, studlab=df$country, byvar= pull(df, 'Time period'),
                  backtransf=TRUE, sm="OR", comb.fixed=TRUE, bylab = 'Time period')

  TE <- ma$TE
  TE[exp(TE) > 1000] <- NA
  ma$TE <- TE

  # Upper and lower limits for confidence intervals have to be set strategically in order to not get an 
  # error from forest(), but also to create a forest plot that isnt nonsensical.
  limits <- setNA(ma[['upper']], ma[['lower']])

  ma[['upper']] <- limits[[1]]
  ma[['lower']] <- limits[[2]]

  limits <- set0(ma[['upper.fixed']], ma[['lower.fixed']])

  ma[['upper.fixed']] <- limits[[1]]
  ma[['lower.fixed']] <- limits[[2]]

  limits <- set0(ma[['upper.fixed.w']], ma[['lower.fixed.w']])

  ma[['upper.fixed.w']] <- limits[[1]]
  ma[['lower.fixed.w']] <- limits[[2]]

  return(list(table, ma))
}

# Combine individual tables into one for publication
create_pub_table <- function(){
  for (i in 1:length(endpoints)) { 
    #i <- 2

    input <- readRDS(paste0(path, 'ma_res_', names(endpoints)[i],".rds") )
    
    if (i == 1){
      AZ_table <- input$az_tab
      PB_table <- input$pb_tab
    } else {
      AZ_table <- rbind(AZ_table, input$az_tab)
      PB_table <- rbind(PB_table, input$pb_tab)
    }
  }
  
  names(AZ_table) <- c(' ', 'England - RCGP', ' ', 'Scotland', ' ', 'Wales', ' ')
  names(PB_table) <- c(' ', 'England - RCGP', ' ', 'Scotland', ' ', 'Wales', ' ')
  
  new_row <- c('Time period', rep( c('Number of cases', 'Number of controls'), 3))
  
  AZ_table <- rbind(new_row, AZ_table)
  PB_table <- rbind(new_row, PB_table)
  
  write.csv(AZ_table, paste(path, 'AZ_table.csv', sep='' ), row.names = FALSE)
  write.csv(PB_table, paste(path, 'PB_table.csv', sep='' ), row.names = FALSE) 
}

###############################################################################

# Named vector of endpoints
endpoints <- c( "Arterial_thromb" = "Arterial Thrombosis",
                "any_haem" =  "Haemorrhagic events",
                "itp" = "Idiopathic Thrombocytopenic Purpura",
                "itp_gen" = "Thrombocytopenic events (excluding ITP)",
                "throm_cvst" = "Venous thromboembolic events" )

# Change this depending on whether main analysis, sensitivty analysis etc
#study <- 'SCCS'
study <- 'case-control-sensitivity'
#study <- 'case-control'

path <- paste0('./output/', study, '/')

############################################################################

# Create figures output lists etc in a loop
for (i in 1:length(endpoints) ) { 
  output_list <- list()
  i <- 3
  
  if(i==3){
    i <- 4
  }
  
  print(i)
  
  output_list$Endpoint = endpoints[[i]]
  output_list$group = names(endpoints)[i]

  analysis_objects <- create_table_ma_plot('AZ', output_list$group) 
  
  output_list$az_tab <- analysis_objects[[1]]
  output_list$az <- analysis_objects[[2]]
  
  png(paste(path, 'AZ_', output_list$group, '_fig.png', sep = ''), width = 1200, height = 750)

  forest(output_list$az, comb.random=FALSE, comb.fixed=TRUE, overall=FALSE, leftcols=c("studlab"), leftlabs=c("Country"),
       label.right = "Higher Risk", label.left="Lower Risk", main="log(OR)", plotwidth = unit(8, "cm"),
       colgap=unit(4, "cm"))

  dev.off()
  

  analysis_objects <- create_table_ma_plot('PB', output_list$group) 
  
  output_list$pb_tab <- analysis_objects[[1]]
  output_list$pb <- analysis_objects[[2]]
  
  png(paste(path, 'PB_', output_list$group, '_fig.png', sep = ''), width = 900, height = 750)
  
  forest(output_list$pb, comb.random=FALSE, comb.fixed=TRUE, overall=FALSE, leftcols=c("studlab"), leftlabs=c("Country"), 
         label.right = "Higher Risk", label.left="Lower Risk", main="log(OR)", plotwidth = unit(8, "cm"),
         colgap=unit(4, "cm"))
  
  dev.off()
  
 

saveRDS(output_list,paste0(path, "/ma_res_" , output_list$group, ".rds"))
}


create_pub_table()
