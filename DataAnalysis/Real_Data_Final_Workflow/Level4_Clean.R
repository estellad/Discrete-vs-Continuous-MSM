############# Load data ############
# Read in GC and impute 0 at the gc, need to do the same for BP
setwd("~/Files/projects/cdp/p0904.521.000/level4/for_edong_student")
# library(plyr) # For density comparison of death and fx times.
library(dplyr)
library(stringr)
#library(nnet)
library(survival)
library(splines)

################################################################################
#################### Exposure Data Discrete Time (5 days) #####################
################################################################################
# 5 days split
dayinterval <- 5
split <- 73 #for 5 days, length.out = 74
#split <- 365 #for 1 day, length.out = 366  # equivalent as continuous time
#split <- 52 #for 7 days, length.out = 53
#split <- 12 #for 30 days, length.out = 13

dtimes_V <- floor(seq(from=0, to=365, length.out=74)[-1]) # for 5 days
dtimes_Vs <- c(0, dtimes_V[1:length(dtimes_V)-1]) # for 5 days
#dtimes_V <- floor(seq(from=0, to=365, length.out=366)[-1]) # for 1 day  # equivalent as continuous time
#dtimes_V <- floor(seq(from=0, to=365, length.out=53)[-1]) # for 7 days
#dtimes_V <- floor(seq(from=0, to=365, length.out=13)[-1]) # for 30 days

gc <- read.csv("gc_final.csv", sep= "\t", header= T)
gc$ikn <- droplevels(factor(gc$study_id))
gc$study_id <- droplevels(factor(gc$study_id))
levels(gc$study_id) <- factor(c(1:nlevels(gc$study_id)))
gc$days_to_servdate <- ifelse(is.na(gc$days_to_servdate), Inf, gc$days_to_servdate)
gc$study_id <- as.numeric(gc$study_id)
gc$yt <- as.numeric(gc$yt)
gc <- gc %>%
  arrange(study_id, Start)

# Recalculate Dtrt in gc, cumdose and cumdose projection
gc <- gc %>%
  mutate(Dtrt = ifelse(daily <= 5 & daily > 0, 1, 
                       ifelse(daily <= 10 & daily > 5, 2, 
                              ifelse(daily < 30 & daily > 10, 3, 
                                     ifelse(daily < 50 & daily >= 30, 4, 5))
                       )))

sapply(as.character(c(1:5)), function(x){summary(gc[gc$Dtrt == x, 'daily'])})
#                   1         2        3        4         5
# Median          ____      ____     ____     ____      ____    # We take median to determine daily dose. 

gc <- gc %>%
  mutate(
    daily_exact = daily,
    daily = ifelse(Dtrt == '1', 5, 
                   ifelse(Dtrt == '2', 10, 
                          ifelse(Dtrt == '3', 20, 
                                 ifelse(Dtrt == '4', 34, 50
                                 )))
                   
    )) 

gc$status <- as.numeric(gc$status)


df_Dinit <- read.csv("df_Dinit.csv", sep= "\t", header= T)
df_Dinit$ikn <- droplevels(factor(df_Dinit$study_id))
df_Dinit$study_id <- droplevels(factor(df_Dinit$study_id))
levels(df_Dinit$study_id) <- factor(c(1:nlevels(df_Dinit$study_id)))
df_Dinit$study_id <- as.numeric(df_Dinit$study_id)
df_Dinit <- df_Dinit %>%
  arrange(study_id)


ltc_conf <- read.csv("ltc_conf.csv", sep= "\t", header= T)
ltc_conf$study_id <- droplevels(factor(ltc_conf$study_id))
levels(ltc_conf$study_id) <- factor(c(1:nlevels(ltc_conf$study_id)))
ltc_conf$days_to_ltcdate <- ifelse(is.na(ltc_conf$days_to_ltcdate), Inf, ltc_conf$days_to_ltcdate)
ltc_conf$study_id <- as.numeric(ltc_conf$study_id)
ltc_conf <- ltc_conf %>%
  arrange(study_id)


bp_final <- read.csv("bp_final.csv", sep= "\t", header= T)
bp_final$ikn <- droplevels(factor(bp_final$study_id))
bp_final <- subset(bp_final, select = -study_id)
bp_final <- left_join(bp_final, df_Dinit[, c('ikn', 'study_id')])
bp_final$study_id <- as.numeric(bp_final$study_id)
bp_final$status <- as.numeric(bp_final$status)
bp_final <- bp_final %>%
  arrange(study_id, Start)%>%
  filter(duration != 0)

bp_final <- bp_final %>%
  group_by(study_id) %>%
  mutate(nextstart = lead(Start)) %>%
  mutate(nextstart = as.numeric(ifelse(is.na(nextstart), 365, nextstart))) %>%
  mutate(overlapdays = ifelse(nextstart < End, End - nextstart, 0), 
         cumoverlapdays = cumsum(overlapdays)) %>%
  mutate(nextstart = nextstart + cumoverlapdays) %>%
  mutate(lagnextstart = lag(nextstart))%>%
  mutate(lagnextstart = as.numeric(ifelse(is.na(lagnextstart), Start, lagnextstart))) %>%
  mutate(Start = lagnextstart,
         End = Start + duration) %>%
  filter(Start < 365) %>%
  mutate(End = ifelse(End > 365, 365, End),
         duration = End - Start, 
         days_to_BPdate = Start) %>%
  dplyr::select(-nextstart, -lagnextstart, -overlapdays, -cumoverlapdays) 

nobs <- length(unique(gc$study_id))  # _____ pt


# Universal Parameters
stepsize = 1 # for sw over time plot in Diagnostics.R
tlim = 365
################# Exposure Data and Models ######################

################# Outcome Models ######################

################# Diagnostics ####################
library(survey)
library(tableone) 



