######################### Survsplit ######################
df_new = survSplit(Surv(Start, End, status) ~., gc, cut = dtimes_V, episode="cycle")
df_newsub <- df_new[, c("study_id", "Start", "End", "Dtrt", "cycle", "daily", "cumdose", "cumdose_projection", "days_to_servdate")]
df_newsub$Dtrt <- as.numeric(df_newsub$Dtrt)

# Manual impute 0 Daily_cat for non-exist months, no need for yt, status, days_to_servdate yet
gc_full <- data.frame(
  study_id=as.numeric(rep(1:nobs, times=rep(length(dtimes_V), nobs))),
  Start=rep(dtimes_Vs, nobs),
  End=rep(dtimes_V, nobs),
  Dtrt = 0,
  cycle = rep(c(1:split), nobs),
  daily = 0
)
gc_full$Dtrt <- as.numeric(gc_full$Dtrt)

# Merge these two together
dat <- full_join(df_newsub, gc_full)

gc_original <- gc %>%
  dplyr::select(study_id, yt)
dat <- dat %>%
  arrange(study_id, cycle, End, Dtrt) %>%
  full_join(gc_original)

# After filter to per study_id per cycle, exposure data is ready 
dat_expo <- dat %>%
  group_by(study_id, cycle) %>%
  filter(row_number() == n())

dat_expo <- dat_expo %>%
  mutate(Start = ifelse(!(Start %in% dtimes_Vs), dtimes_Vs[cycle], Start))
#,End = ifelse(End>yt, yt, End))     # No need to cut short End, for both simulation and the data, because no futime 

dat_expo <- dat_expo[dat_expo$Start < dat_expo$yt, ]

################ After exclude after death, we plot line plot, ###########################################
################ so that Dtrt=0 means strictly unexposed, rather than dead or stopped follow-up ##########
Dinit <- gc %>%
  group_by(study_id) %>%
  filter(row_number() == 1) %>%
  dplyr::rename(D0 = Dtrt)

if(doselineplot){
  source("/users/edong/Code/Dose_line_plot.R", echo=T)
}

# Initial D data set
df_Dinit <- Dinit[, c('study_id', 'D0')] %>%
  full_join(subset(df_Dinit, select = -D0))

####################################### Map in initial cumdose (take out init and map to all) ############
gc$cumdose <- as.numeric(ifelse(is.na(gc$cumdose), 0, as.character(gc$cumdose)))
gc_cum_all_expo_init <- gc %>%
  group_by(study_id)%>%
  filter(row_number() == 1) %>%
  select(study_id, cumdose)

dat_expo <- dat_expo %>%
  select(study_id, Start, End, Dtrt, cycle, daily, days_to_servdate, yt) 

# dat_expo$study_id <- as.factor(dat_expo$study_id)
# gc_cum_all_expo_init$study_id <- as.factor(gc_cum_all_expo_init$study_id)
dat_expo_cum <- dat_expo %>%
  full_join(gc_cum_all_expo_init)

########## Calculate GC cumdose here #########
gc_cum_cleaned <- dat_expo_cum %>%
  group_by(study_id) %>%
  mutate(
    Duration = End-Start,
    gccuminterval = cumsum(daily*Duration),
    gccumdose_projection = cumdose + gccuminterval,
    gccumdose_out = lag(gccumdose_projection),
    gccumdose = lag(gccumdose_out)) %>%
  mutate(gccumdose_out = ifelse(is.na(gccumdose_out), cumdose, gccumdose_out),
         gccumdose = ifelse(is.na(gccumdose), cumdose, gccumdose)) %>%
  select(-gccuminterval, -cumdose)

gc_cum_cleaned <- gc_cum_cleaned %>%
  mutate(days_to_servdate = ifelse(!is.na(days_to_servdate) & days_to_servdate < Start, Inf, days_to_servdate)) 

gc_cum_cleaned$days_to_servdate[which(gc_cum_cleaned$days_to_servdate == 0.0)] <- Inf
gc_cum_cleaned$days_to_servdate <- ifelse(is.na(gc_cum_cleaned$days_to_servdate), Inf, gc_cum_cleaned$days_to_servdate)

# Lag one for D_k-1, for label visit
dat_lagd_cum <- gc_cum_cleaned %>%
  group_by(study_id) %>%
  mutate(Dtrtp = lag(Dtrt)) 
dat_lagd_cum$Dtrtp = ifelse(is.na(dat_lagd_cum$Dtrtp), 
                            dat_lagd_cum$Dtrt, dat_lagd_cum$Dtrtp) # label the time0 Dtrt to as epi1 Dtrtp

# Label Visit
dat_lagd_cum$Visit <- ifelse(dat_lagd_cum$Dtrt != dat_lagd_cum$Dtrtp, 1, 0); table(dat_lagd_cum$Visit) 
# For those had Visit happen in an interval, but stayed on the same dose
dat_lagd_cum$Visit <- ifelse(dat_lagd_cum$Dtrt == dat_lagd_cum$Dtrtp & is.finite(dat_lagd_cum$days_to_servdate)&
                               dat_lagd_cum$days_to_servdate <= dat_lagd_cum$End & dat_lagd_cum$days_to_servdate > dat_lagd_cum$Start,
                             1, dat_lagd_cum$Visit)
table(dat_lagd_cum$Visit) 

# subset to those have Visit = 1 and check
# dat_visit_check <- dat_lagd_cum[dat_lagd_cum$Visit == 1, 
#                                 c("study_id", "Start", "End", "Dtrt", "Dtrtk_1", "Visit", "days_to_servdate")]


################## LAG LTC ###################
# TODO: Merge in LTC and code Ltc_conf -> TDconf delay one interval
long_V_slm <- full_join(dat_lagd_cum, ltc_conf, by="study_id")
long_V_slm <- long_V_slm %>%
  mutate(Ltcconf = ifelse(days_to_ltcdate <= End, 1, 0)) %>%
  group_by(study_id) %>%
  mutate(Ltc = lag(Ltcconf)) %>%
  mutate(Ltc = ifelse(is.na(Ltc), 0, Ltc)) %>%
  select(-days_to_ltcdate, -Ltcconf)


################# BP cum days ################
nobsbp <- length(unique(bp_final$study_id))

bp_cum_new <- survSplit(Surv(Start, End, status) ~., bp_final, cut = dtimes_V, episode="cycle")

bp_cum_new <- bp_cum_new %>%
  mutate(Duration = End - Start) %>%
  group_by(study_id, cycle) %>%
  dplyr::summarise(BPdays = sum(Duration))

bpstudy_idlist <- unique(bp_final$study_id)
bp_full <- data.frame(
  study_id= as.numeric(rep(bpstudy_idlist, times=rep(length(dtimes_V), nobsbp))),
  cycle = rep(c(1:split), nobsbp),
  BPdays = 0
)

bp_cum_all <- full_join(bp_cum_new, bp_full)
bp_cum_all <- bp_cum_all %>%
  arrange(study_id, cycle, BPdays) %>%
  group_by(study_id, cycle) %>%
  filter(row_number() == n()) 

bp_cum_all <- bp_cum_all %>%
  group_by(study_id) %>%
  mutate(cumBPdays = cumsum(BPdays), 
         cumBPdaysp = lag(cumBPdays)) %>%
  mutate(cumBPdaysp = ifelse(is.na(cumBPdaysp), 0, cumBPdaysp)) %>%
  select(study_id, cycle, cumBPdays, cumBPdaysp)


#### After cleaning those with BP days, merge in full long_V_slm and clean those with 0 BP days 
long_V_sub <- long_V_slm %>%
  select(study_id, cycle) 

bp_cum_all_cleaned <- right_join(bp_cum_all, long_V_sub)

bp_cum_all_cleaned <- bp_cum_all_cleaned %>%
  mutate(cumBPdays = ifelse(is.na(cumBPdays), 0, cumBPdays),
         cumBPdaysp = ifelse(is.na(cumBPdaysp), 0, cumBPdaysp)) 
bp_cum_all_cleaned <- bp_cum_all_cleaned[, c("study_id", "cycle", "cumBPdays", "cumBPdaysp")]

long_V_slm <- long_V_slm %>%
  full_join(bp_cum_all_cleaned, by=c("study_id", "cycle"))

# Merge in baseline conf
long_V_slm <- long_V_slm %>%
  full_join(df_Dinit)


########################## Factorize and relevel #######################
# Exposure D0 data
df_Dinit$D0 <- as.factor(df_Dinit$D0); df_Dinit$D0 <- relevel(df_Dinit$D0, ref="1")
facvector <- c('arthritis', 'pfx', 'inhalegc', 'sexhor', 'thiazide', 
               'falldr', 'fallcomorb', 'emerhosp','diabetes', 'pralcal')
df_Dinit[, facvector] <- lapply(df_Dinit[, facvector], as.factor)

# Exposure V data
long_V_slm$D0 <- as.factor(long_V_slm$D0); long_V_slm$D0 <- relevel(long_V_slm$D0, ref="1")
long_V_slm$Dtrt <- as.factor(long_V_slm$Dtrt); long_V_slm$Dtrt <- relevel(long_V_slm$Dtrt, ref="1")
long_V_slm$Dtrtp <- as.factor(long_V_slm$Dtrtp); long_V_slm$Dtrtp <- relevel(long_V_slm$Dtrtp, ref="1")

# ### TO Delete: Actually no need of cumGCdays for PO graph, as always treat means POdaily(End-0)
# long_V_slm <- long_V_slm %>%
#   mutate(gcexpostatus = ifelse(daily != 0, 1, 0)) %>%
#   group_by(study_id) %>%
#   mutate(GCdaysinterval = cumsum(gcexpostatus * Duration)) # %>%
# #   mutate(cumGCdays = lag(GCdaysinterval)) %>%
# #   mutate(cumGCdays = ifelse(is.na(cumGCdays), 0, cumGCdays))


# Exposure D data
Ddata <- long_V_slm[long_V_slm$Visit == 1,]
Ddata$D0 <- as.factor(Ddata$D0); Ddata$D0 <- relevel(Ddata$D0, ref="1")
Ddata$Dtrt <- as.factor(Ddata$Dtrt); Ddata$Dtrt <- relevel(Ddata$Dtrt, ref="1")
Ddata$Dtrtp <- as.factor(Ddata$Dtrtp); Ddata$Dtrtp <- relevel(Ddata$Dtrtp, ref="1")


########### Center all continuous vars in exposure data set ############
df_Dinit$agescale <- as.numeric(scale(df_Dinit$age, scale = FALSE))

long_V_slm$agescale <- as.numeric(scale(long_V_slm$age, scale = FALSE))
long_V_slm$gccumdosescale <- as.numeric(scale(long_V_slm$gccumdose, scale = FALSE))
long_V_slm$cumBPdayspscale <- as.numeric(scale(long_V_slm$cumBPdaysp, scale = FALSE))

Ddata$agescale <- as.numeric(scale(Ddata$age, scale = FALSE)) 
Ddata$gccumdosescale <- as.numeric(scale(Ddata$gccumdose, scale = FALSE)) 
Ddata$cumBPdayspscale <- as.numeric(scale(Ddata$cumBPdaysp, scale = FALSE)) 


##### Appendix #####
# gc_cum_cleaned <- NULL
# for(i in 1:nobs){
#   gc_cumi <- gc_cum_all_expo[gc_cum_all_expo$study_id == i, ]
#   if(nrow(gc_cumi)==1){
#     daygap = gc_cumi$End
#   }else{
#     daygap = dayinterval
#   }
#   gc_cumi$cumdose_projection[1] <- gc_cumi$cumdose[1] + gc_cumi$daily[1] * daygap
#   
#   if(nrow(gc_cumi)>2){
#     for(j in 2:nrow(gc_cumi)){
#       gc_cumi$cumdose[j] <- gc_cumi$cumdose_projection[j-1]
#       gc_cumi$cumdose_projection[j] <- gc_cumi$cumdose[j] + gc_cumi$daily[j] * daygap
#     }
#   }
#   gc_cum_cleaned <- rbind(gc_cum_cleaned, gc_cumi)
#   if(i %% 10 ==0 ) print(i)
# }
