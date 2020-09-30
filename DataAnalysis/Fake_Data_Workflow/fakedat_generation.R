rm(list=ls())
library(dplyr)
library(survival)
library(nnet)
gc_final <- read.csv("~/Desktop/Thesis_Code/Fully_Cont_MSM_MPP/DataAnalysis/Fake_Data_Workflow/GC_FINAL.csv")

get_gc <- function(gc_final){
  gc_final$servdate <- as.Date(as.character(gc_final$servdate), format = "%d-%b-%y")
  gc_final$enddate <- as.Date(as.character(gc_final$enddate), format = "%d-%b-%y")
  gc_final$followdate <- as.Date(as.character(gc_final$followdate), format = "%d-%b-%y")
  gc_final$admcensdate <- as.Date(as.character(gc_final$admcensdate), format = "%d-%b-%y")
  
  gc_final2 <- gc_final %>%
    mutate(Dtrt = case_when(daily <= 5 & daily > 0 ~ "0 < daily dose <= 5mg",
                            daily <= 10 & daily > 5 ~ "5 < daily dose <= 10mg",
                            daily < 30 & daily > 10 ~ "10 < daily dose < 30mg",
                            daily < 50 & daily >= 30 ~ "30 <= daily dose < 50mg",
                            daily <= 60 & daily >= 50  ~ "50mg <= daily dose <= 60mg"
          ),
    Start = as.numeric(servdate - followdate),
    End = as.numeric(enddate - followdate)) %>%
    mutate(Start = ifelse(Start != 0 & Start %% 5 == 0, Start-0.1, Start),
           Servdate = Start)
  
  gc_final2$End <- ifelse(gc_final2$End > 365, 365, gc_final2$End)
  
  gc_final2 <- gc_final2 %>%
    mutate(status = case_when(ikn == 1 ~ 0,
                              ikn == 2 ~ 1,
                              ikn == 3 ~ 1, 
                              ikn == 4 ~ 0,
                              ikn == 7 ~ 1), 
           yt = case_when(ikn == 1 ~ 365,
                          ikn == 2 ~ 173,
                          ikn == 3 ~ 199, 
                          ikn == 4 ~ 365,
                          ikn == 7 ~ 30))
  gc_final2 <- gc_final2[gc_final2$Start < gc_final2$yt, ]
  gc_final2$End <- ifelse(gc_final2$End > gc_final2$yt, gc_final2$yt, gc_final2$End)
  gc_final2$Duration <- gc_final2$End - gc_final2$Start
  
  gc_final3 <- gc_final2[, c("ikn",  "Start", "End", "Duration", "Dtrt", "status", "yt",
                             "cumdose", "daily", "Servdate")]

  gc_final3$Dtrt <- as.factor(gc_final3$Dtrt)
  levels(gc_final3$Dtrt) <- factor(c(1:nlevels(gc_final3$Dtrt)))
  
  gc_final3$Start <- as.numeric(gc_final3$Start)
  gc_final3$Servdate <- gc_final3$Start
  gc_final3$End <- as.numeric(gc_final3$End)
  gc_final3$Duration <- as.numeric(gc_final3$Duration)
  gc_final3$ikn <- as.factor(gc_final3$ikn)
  levels(gc_final3$ikn) <- as.factor(c(1:5))
  
  gc_final3$status <- as.numeric(gc_final3$status)
  gc_final3$drugtype <- 1
  return(gc_final3)
}

# get_gc_ikn_cum_servd <- function(gc_final){
#   gc_final$servdate <- as.Date(as.character(gc_final$servdate), format = "%d-%b-%y")
#   gc_final$enddate <- as.Date(as.character(gc_final$enddate), format = "%d-%b-%y")
#   gc_final$followdate <- as.Date(as.character(gc_final$followdate), format = "%d-%b-%y")
#   gc_final$admcensdate <- as.Date(as.character(gc_final$admcensdate), format = "%d-%b-%y")
#   
#   gc_final2 <- gc_final %>%
#     mutate(Dtrt = case_when(daily <= 5 & daily > 0 ~ "0 < daily dose <= 5mg",
#                             daily <= 10 & daily > 5 ~ "5 < daily dose <= 10mg",
#                             daily < 30 & daily > 10 ~ "10 < daily dose < 30mg",
#                             daily < 50 & daily >= 30 ~ "30 <= daily dose < 50mg",
#                             daily <= 60 & daily >= 50  ~ "50mg <= daily dose <= 60mg"
#     ),
#     Start = as.numeric(servdate - followdate),
#     End = as.numeric(enddate - followdate)) 
#   
#   gc_final2$End <- ifelse(gc_final2$End > 365, 365, gc_final2$End)
#   
#   gc_final2 <- gc_final2 %>%
#     mutate(status = case_when(ikn == 1 ~ 0,
#                               ikn == 2 ~ 1,
#                               ikn == 3 ~ 1, 
#                               ikn == 4 ~ 0,
#                               ikn == 7 ~ 1), 
#            yt = case_when(ikn == 1 ~ 365,
#                           ikn == 2 ~ 173,
#                           ikn == 3 ~ 199, 
#                           ikn == 4 ~ 365,
#                           ikn == 7 ~ 30))
#   gc_final2 <- gc_final2[gc_final2$Start < gc_final2$yt, ]
#   gc_final2$End <- ifelse(gc_final2$End > gc_final2$yt, gc_final2$yt, gc_final2$End)
#   gc_final2$Duration <- gc_final2$End - gc_final2$Start
#   
#   gc_final2$Servdate <- gc_final2$Start
#   gc_final2$Enddate <- gc_final2$End
#   gc_final2 <- gc_final2 %>%
#     group_by(ikn) %>%
#     mutate(cumdose_projection = lead(cumdose)) %>%
#     mutate(cumdose_projection = ifelse(
#       is.na(cumdose_projection), Duration * daily + cumdose, cumdose_projection
#   ))
# 
#   gc_final2$Dtrt <- as.factor(gc_final2$Dtrt)
#   levels(gc_final2$Dtrt) <- factor(c(1:nlevels(gc_final2$Dtrt)))
#   gc_final2$ikn <- as.factor(gc_final2$ikn)
#   levels(gc_final2$ikn) <- as.factor(c(1:5))
#   gc_ikn_cum_servd <- gc_final2[, c("ikn", "cumdose", "Servdate", "Enddate", "Dtrt",
#                                     "Start", "End", "cumdose_projection", "status", "daily")]
#   return(gc_ikn_cum_servd)
# }

get_baseline <- function(){
  baseline <- data.frame(
    ikn = as.factor(c(1,2,3,4,5)),
    diab =  c(1, 0, 1, 1, 0),
    hormone = c(0, 0, 1, 1, 1)
  )
  return(baseline)
}

get_ltc_final <- function(){
  ltc_final <- data.frame(matrix(vector(), 0, 4))
  ltc_final <- rbind(ltc_final, 
                    c(1, 60, 0, 365),
                    c(3, 20, 1, 199),
                    c(4, 228,0, 365),
                    c(5, 21, 1, 30))
  colnames(ltc_final) <- c("ikn",  "Ltcdate", "status", "eventdate")
  ltc_final$ikn <- as.factor(ltc_final$ikn)
  return(ltc_final)
}

get_bp_days_final <- function(){
  bp_final <- data.frame(matrix(vector(), 0, 4))
  bp_final <- rbind(bp_final, 
                     c(1, 22, 99,  0),
                     c(1, 122, 2, 0),
                     c(1, 124, 2, 0),
                     c(1, 127, 3, 0),
                     c(1, 226, 5, 0),
                     c(1, 310, 20,0),
                     c(2, 101, 17,1),
                     c(2, 137, 8, 1),
                     c(2, 151, 20, 1),
                     c(3, 35,  5, 1),
                     c(3, 84,  30, 1),
                     c(3, 129, 70, 1),
                     c(4, 321, 1, 0))
  colnames(bp_final) <- c("ikn",  "BPdate", "duration", "status")
  bp_final$ikn <- as.factor(bp_final$ikn)
  return(bp_final)
}

# #######################################################################
# ## Stacked Bar Plot - Sequential Ordinal Logistic
# # 15 Days Interval
# df$End <- df$Start+df$Duration
# df <- df[, c("ikn",  "Start", "End", "Daily_cat", "status")]
# library(survival)
# day_sep_list <- seq(from = 0, to = 365.25, length.out = 25)[-1] # 15.21875 as separation
# df_new = survSplit(Surv(Start, End, status) ~., df, cut=day_sep_list, episode ="half_month")
# df_new$half_month[which(!(df_new$End %in% day_sep_list))] = 0
# df_new <- df_new[df_new$half_month != 0, ]
# 
# tab <- table(df_new$Daily_cat, df_new$half_month)
# data <- rbind(tab[1,], tab[2,], tab[3,], tab[4,])
# 
# colnames(data) <- as.factor(sort(unique(df_new$half_month)))
# rownames(data) <- unique(df_new$Daily_cat)
# 
# data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
# 
# library(RColorBrewer)
# coul <- brewer.pal(4, "Pastel2")
# par(mar = c(4, 4, 2,9))
# barplot(data_percentage, beside = FALSE, legend=FALSE, xlab="Half Months", main="Half Monthly dosage proportion table",
#         col=coul , border="white")
# legend("right", xpd = TRUE, inset = c(-0.25, 0),
#        legend=c(#'Unexposed', 
#          'daily dose < 5mg', '5 <= daily dose < 20mg',
#          '20 <= daily dose < 30mg', '30 <= daily dose < 50mg'),
#        pch=c(20,20,20,20), col=coul, cex=0.7, bty='n')

