library(readxl)
library(plyr)
library(tidyverse)
library(meta)

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

#Read in the cohrt from 01a_Vaccination_Input.R in Vaccine folder
project_path <- paste0(Location,"EAVE/GPanalysis/progs/CR/Safety_Meta_Analysis")

fun_calc_0_27 <- function(df){
  df <- df %>% mutate(rlrr = R*log_rr, rvlrr = R*(se_log_rr^2) ) %>% 
    dplyr::summarise_at(vars("N","R","rlrr","rvlrr"), ~sum(.) )
  df <- df %>% mutate(log_rr=rlrr/R, se_log_rr = sqrt(rvlrr/R), 
                      log_lcl = log_rr - 1.96*se_log_rr, log_ucl = log_rr + 1.96*se_log_rr ) %>% 
    mutate(RR= exp(log_rr), LCL = exp(log_lcl), UCL = exp(log_ucl)  ) %>% 
    dplyr::select(N, R, RR, LCL, UCL, log_rr, log_lcl, log_ucl, se_log_rr)
  df
}


#read in scotland data
scot <- read.csv(paste0(project_path,"/data/scotland_ma_results.csv"))
scot <- scot %>%  dplyr::rename(RR=HR, Status=vs_type) %>% 
  mutate(log_rr = log(RR), log_lcl = log(LCL), log_ucl = log(UCL)) %>% 
  mutate(se_log_rr = (log_ucl - log_lcl)/4)


#read _ew data
ew_group_names <- c("throm_cvst","any_throm","any_haem","any_itp","itp_gen","itp")
rm(z)
for (i in 1:6) { 
#i <- 1
z_in <- read_xlsx(paste0(project_path,"/data/DaCVaPVaccineSafetyforCR.xlsx"), sheet=i)
names(z_in) <- c(names(z_in)[1:5],"RR","LCL","UCL", "log_rr","se_log_rr")
df <- filter(z_in, Status %in% paste0("AstraZeneca_v1_", c("0:6","7:13","14:20","21:27")))
df <- fun_calc_0_27(df)
df_az <- bind_cols(data.frame(Endpoint=unique(z_in$Endpoint), Country="England (RGCP)", Status = "AstraZeneca_v1_0:27"), df)
df <- filter(z_in, Status %in% paste0("Pfizer_v1_", c("0:6","7:13","14:20","21:27")))
df <- fun_calc_0_27(df)
df_pb <- bind_cols(data.frame(Endpoint=unique(z_in$Endpoint), Country="England (RGCP)", Status = "Pfizer_v1_0:27"), df)
z_in <- bind_rows(z_in, df_az, df_pb)
z_in$group <- ew_group_names[i]
if (exists("z")) z <- bind_rows(z, z_in) else z <- z_in
}

#Can't get dplyr::relocate to work for some reason. Have to do it manually
eng <- z %>% dplyr::relocate(group, .after=Endpoint) %>% dplyr::relocate(se_log_rr, .after=last_col())


eng <- eng %>% mutate(Status = gsub("AstraZeneca", "AZ",Status)) %>% 
  mutate(Status = gsub("Pfizer", "PB",Status))
eng <- eng %>% dplyr::rename(country = Country)
eng <- eng %>%  mutate(country="England_RGCP")


wales <- read_csv(paste0(project_path,"/data/Meta-Analysis_Wales_Coeficients.csv"))
wales <- wales %>% filter(model_type %in% c("fully_adjusted", "reference")) %>% dplyr::select(-p_event, -statistic, -p.value)
names(wales) <- c("Endpoint","model_type","Status","N","R","log_rr","se_log_rr","RR","LCL","UCL")
wales <- wales %>% dplyr::select(-model_type)
wales <- wales %>% mutate( group = c( rep("any",13), rep("throm_cvst", 13), rep("any_haem",13), rep("any_itp", 13), 
                                      rep("any_throm", 13), rep("itp",13),  rep("itp_gen",13)  ))
wales <- wales %>% dplyr::select(Endpoint, group, Status, N, R, RR:UCL, log_rr, se_log_rr) 
wales <- wales %>% mutate(Status = gsub("UV","uv",Status),
                      Status = gsub(" Dose 1 Day ","_v1_",Status),
                      Status = gsub("00-06","0:6",Status),
                      Status = gsub("07-13","7:13",Status),
                      Status = gsub("14-20","14:20",Status),
                      Status = gsub("21-27","21:27",Status),
                      Status = gsub("00-28","0:27",Status))
wales <- wales %>%  mutate(country = "Wales") %>% 
  dplyr::relocate(country, .after=group)
# This doesn't work, looks incomplete. I remove rows with NA values for RR instead
#wales <- wales %>% filter(!is.)

wales <- filter(wales, !is.na(RR))



df <- bind_rows(dplyr::select(eng, -log_lcl, -log_ucl),
                dplyr::select(scot, -log_lcl, -log_ucl) ,wales)
df <- df %>% mutate(Status = factor(Status, levels =c("uv", "AZ_v1_0:6","AZ_v1_7:13","AZ_v1_14:20", "AZ_v1_21:27", "AZ_v1_28+" ,"AZ_v1_0:27", "AZ_v2",
                                                      "PB_v1_0:6", "PB_v1_7:13",   "PB_v1_14:20", "PB_v1_21:27", "PB_v1_28+", "PB_v1_0:27","PB_v2"))) %>% 
  arrange(Endpoint, country, Status)
df <- df %>% mutate(vaccine = case_when(grepl("AZ", Status) ~ "AZ",
                                        grepl("PB", Status) ~ "PB",
                                        TRUE ~ "uv")) %>% 
  mutate(time = as.character(Status)) %>% 
  mutate(time = gsub("AZ_","",time)) %>% 
  mutate(time=gsub("PB_", "", time)) %>% 
  mutate(time=factor(time, levels =c("uv","v1_0:6","v1_7:13","v1_14:20","v1_21:27","v1_0:27","v1_28+" ,  "v2")))

