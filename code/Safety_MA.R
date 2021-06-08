library(readxl)
library(plyr)
library(tidyverse)
library(meta)

z <- read_excel("haem.xlsx")
z <- z %>% mutate(log_rr = log(RR), log_lcl = log(LCL), log_ucl = log(UCL)) %>% 
  mutate(se_log_rr = (log_ucl - log_lcl)/4)

z_ew <- filter(z, Country=="RGCP") %>% 
  dplyr::select(Status, N, R, log_rr, se_log_rr)


fun_calc_0_27 <- function(df){
  df <- df %>% mutate(rlrr = R*log_rr, rvlrr = R*(se_log_rr^2) ) %>% 
    dplyr::summarise_at(vars("N","R","rlrr","rvlrr"), ~sum(.) )
  df <- df %>% mutate(log_rr=rlrr/R, se_log_rr = sqrt(rvlrr/R), 
                      log_lcl = log_rr - 1.96*se_log_rr, log_ucl = log_rr + 1.96*se_log_rr ) %>% 
    mutate(RR= exp(log_rr), LCL = exp(log_lcl), UCL = exp(log_ucl)  ) %>% 
    dplyr::select(N, R, RR, LCL, UCL, log_rr, log_lcl, log_ucl, se_log_rr)
  df
}
df <- filter(z_ew, Status %in% paste0("AZ_v1_", c("0:6","7:13","14:20","21:27")))
df <- fun_calc_0_27(df)
df_az <- bind_cols(data.frame(Endpoint="Haemorrhage", Country="RGCP", Status = "AZ_v1_0:27"), df)
df <- filter(z_ew, Status %in% paste0("PB_v1_", c("0:6","7:13","14:20","21:27")))
df <- fun_calc_0_27(df)
df_pb <- bind_cols(data.frame(Endpoint="Haemorrhage", Country="RGCP", Status = "PB_v1_0:27"), df)
df <- bind_rows(z,df_az,df_pb)
df <- df %>% mutate(Status = factor(Status, levels =c("uv", "AZ_v1_0:6","AZ_v1_7:13","AZ_v1_14:20", "AZ_v1_21:27", "AZ_v1_28+" ,"AZ_v1_0:27", "AZ_v2",
                        "PB_v1_0:6", "PB_v1_7:13",   "PB_v1_14:20", "PB_v1_21:27", "PB_v1_28+", "PB_v1_0:27","PB_v2"))) %>% 
     arrange(Endpoint, Country, Status)
df <- df %>% mutate(vaccine = case_when(grepl("AZ", Status) ~ "AZ",
                                       grepl("PB", Status) ~ "PB",
                                       TRUE ~ "uv")) %>% 
  mutate(time = as.character(Status)) %>% 
  mutate(time = gsub("AZ_","",time)) %>% 
  mutate(time=gsub("PB_", "", time)) %>% 
  mutate(time=factor(time, levels =c("uv","v1_0:6","v1_7:13","v1_14:20","v1_21:27","v1_0:27","v1_28+" ,  "v2")))

z.df <- filter(df, Status=="AZ_v1_0:27")
z.ma <- metagen(TE=z.df$log_rr, seTE=z.df$se_log_rr, studlab=z.df$Country, backtransf=TRUE, sm="OR", comb.fixed=TRUE)
forest(z.ma, comb.random=FALSE, leftcols=c("studlab"), leftlabs=c("Country"), 
       label.right = "Higher Risk", label.left="Lower Risk", main="log(OR)")


output_list <- list()
output_list$Endpoint = "Haemorrhage"


z.df <- filter(df, time=="v1_0:27")
z.ma <- metagen(TE=z.df$log_rr, seTE=z.df$se_log_rr, studlab=z.df$Country, byvar=z.df$vaccine,
                backtransf=TRUE, sm="OR", comb.fixed=TRUE)
output_list$days_0_27 <- z.ma
forest(z.ma, comb.random=FALSE, comb.fixed=TRUE, overall=FALSE, leftcols=c("studlab"), leftlabs=c("Country"), 
       label.right = "Higher Risk", label.left="Lower Risk", main="log(OR)")


z.df <- filter(df, time %in% c("v1_0:6","v1_7:13","v1_14:20","v1_21:27","v1_28+") & vaccine=="AZ")
z.ma <- metagen(TE=z.df$log_rr, seTE=z.df$se_log_rr, studlab=z.df$Country, byvar=z.df$time,
                backtransf=TRUE, sm="OR", comb.fixed=TRUE)
forest(z.ma, comb.random=FALSE, comb.fixed=TRUE, overall=FALSE, leftcols=c("studlab"), leftlabs=c("Country"), 
       label.right = "Higher Risk", label.left="Lower Risk")
output_list$az <- z.ma
z.df <- filter(df, time %in% c("v1_0:6","v1_7:13","v1_14:20","v1_21:27","v1_28+") & vaccine=="PB")
z.ma <- metagen(TE=z.df$log_rr, seTE=z.df$se_log_rr, studlab=z.df$Country, byvar=z.df$time,
                backtransf=TRUE, sm="OR", comb.fixed=TRUE)
output_list$pb <- z.ma

saveRDS(output_list,"ma_res.RDS")
