############### For cb ##############
gc_init_all <- gc %>%
  group_by(study_id) %>%
  filter(row_number() == 1) 

################### Manual imput data for unexposed here #################
gc_new <- gc %>%
  group_by(study_id) %>%
  mutate(lagEnd = lag(End)) %>%
  mutate(lagEnd = ifelse(is.na(lagEnd), 0, lagEnd))

# gc_new_first2 <- gc_new %>%
#   filter(study_id %in% c(1,2))

gc_new_first2_try <- gc_new[, c('study_id', 'Start', 'End', 'Dtrt', 'yt', 'days_to_servdate', 'daily', 'lagEnd')] %>%
  mutate(study_idn = ifelse(Start != lagEnd, study_id, Inf),
         Startn = ifelse(Start != lagEnd, lagEnd, Inf),
         Endn = ifelse(Start != lagEnd, Start, Inf),
         Dtrtn = ifelse(Start != lagEnd, 0, Inf),
         dailyn = ifelse(Start != lagEnd, 0, Inf),
         days_to_servdaten = Inf,
         ytn = ifelse(Start != lagEnd, yt, Inf))

gc_new_first2_try_manually <- gc_new_first2_try[, c('study_idn', 'Startn', 'Endn', 'Dtrtn', 'dailyn', 'days_to_servdaten', 'ytn')] %>%
  filter(is.finite(study_idn)) %>%
  dplyr::rename(study_id = study_idn,
                Start = Startn,
                End = Endn,
                Dtrt = Dtrtn,
                daily = dailyn,
                days_to_servdate = days_to_servdaten,
                yt = ytn)

gc_new_first2_try <- gc_new_first2_try[, !(colnames(gc_new_first2_try) 
                                           %in% c("lagEnd", "study_idn", "Startn", "Endn", "Dtrtn", 
                                                  "dailyn", "days_to_servdaten", "ytn"))]
#gc_new_first2_try_manually$study_id <- as.factor(gc_new_first2_try_manually$study_id)
gc_first2 <- full_join(gc_new_first2_try_manually, gc_new_first2_try)

gc_first2 <- gc_first2 %>%
  arrange(study_id, Start)

################## Makeup last interval until event
gc_first2_last <- gc_first2 %>%
  group_by(study_id) %>%
  filter(row_number() == n()) 
gc_first2_last$study_idn = ifelse(gc_first2_last$Start != gc_first2_last$yt, gc_first2_last$study_id, Inf)
gc_first2_last$Startn = ifelse(gc_first2_last$Start != gc_first2_last$yt, gc_first2_last$End, Inf)
gc_first2_last$Endn = ifelse(gc_first2_last$Start != gc_first2_last$yt, gc_first2_last$yt, Inf)
gc_first2_last$Dtrtn = ifelse(gc_first2_last$Start != gc_first2_last$yt, 0, Inf)
gc_first2_last$dailyn = ifelse(gc_first2_last$Start != gc_first2_last$yt, 0, Inf)
gc_first2_last$days_to_servdaten = Inf
gc_first2_last$ytn = ifelse(gc_first2_last$Start != gc_first2_last$yt, gc_first2_last$yt, Inf)

gc_new_first2_last_manually <- gc_first2_last[, c('study_idn', 'Startn', 'Endn', 'Dtrtn', 'dailyn', 'days_to_servdaten', 'ytn')] %>%
  filter(is.finite(study_idn)) %>%
  dplyr::rename(study_id = study_idn,
                Start = Startn,
                End = Endn,
                Dtrt = Dtrtn,
                daily = dailyn,
                days_to_servdate = days_to_servdaten,
                yt = ytn)

# get rid of repetitive study_id col before merge
gc_full_final <- full_join(gc_first2, gc_new_first2_last_manually)
gc_done <- gc_full_final %>%
  arrange(study_id, Start)

# Merge in status 
gc_done <- full_join(gc_done, gc_init_all[, c("study_id", "status")])

# For CTM, I don't need to modify 20 (numbers that %% 5=0) to 19.9, so change back to original
gc_done$days_to_servdate = ifelse(gc_done$days_to_servdate %% 1 != 0, gc_done$days_to_servdate + 0.1, gc_done$days_to_servdate)
gc_done$Start = ifelse(gc_done$Start %% 1 != 0, gc_done$Start + 0.1, gc_done$Start)
gc_done$End = ifelse(gc_done$End %% 1 != 0, gc_done$End + 0.1, gc_done$End)

########################## Too huge data set, I decide to split female male here ##########################
# Use male as an example for smaller dataset (later need to split bp as well)
# merge in sex to split here 
base_sex <- df_Dinit[,c('study_id', 'sex')]

gc_done <- gc_done %>%
  inner_join(base_sex)
gc_done <- gc_done[gc_done$sex == gend, ]   # male n=____, actual sampled in outcome: ____ -> now ____ split before cb: ____ unselected
length(unique(gc_done$study_id))            # female n=____, actuallu sampled -> now ____, so ____ unselected
gc_init <- gc_init_all[gc_init_all$study_id %in% unique(gc_done$study_id), ]  # only sex==1 or 0

################################ Fracture times for survsplit #########################
tlim = 365
cs <- gc_init[gc_init$status==1, c("study_id", "yt")]
cs$Stop <- cs$yt
cs$Yevent <- 1                                      # male ___ events, total cb time point = ___*201 = ____

dtimes_Y <- sort(unique(cs$yt))
########### Base series data ##########
ftime = as.numeric(pmin(gc_init$yt, tlim))
ms <- nrow(cs)*200

set.seed(5)
persons <- rep(unique(gc_init$study_id), as.numeric(rmultinom(1, ms, ftime/sum(ftime))))
length(unique(persons)) # male: ____, female: ___
list <- union(unique(persons), unique(cs$study_id))
length(unique(list))    # male tot: ____; female tot: ____

ftimedf <- as.data.frame(cbind(gc_init$study_id, ftime))
ftimerow_list <- match(persons, ftimedf$V1)
moments <- rep(NA, ms)
for (i in 1:ms){
  moments[i] <- runif(1, min=0.0, max=ftime[ftimerow_list[i]])
}

bs <- gc_init_all[persons, c("study_id", "yt")]               # now back to all pt gc_init to take out corresponding row with index persons
bs$Stop <- moments
bs$Yevent <- 0
a1 <- length(unique(bs$study_id)); a2<- length(unique(cs$study_id)); 
cbid <- unique(union(bs$study_id, cs$study_id)); length(cbid)

casebase <- rbind(cs,bs)   ; length(unique(casebase$study_id))         # male = ____ pt; female = ____ pt; total ____, ____ unselected
casebase$study_id <- as.numeric(casebase$study_id)
casebase <- casebase %>%
  group_by(study_id)%>%
  arrange(study_id, Stop) 

#length(unique(casebase$study_id))
casebase$Yevent <- as.logical(casebase$Yevent)

casebase <- casebase %>%
  group_by(study_id) %>%
  mutate(Begin = lag(Stop)) %>%
  mutate(Begin = ifelse(is.na(Begin), 0, Begin))

gc_done <- gc_done[gc_done$study_id %in% unique(casebase$study_id), ] # subset to only those selected in casebase sampling
#length(unique(gc_done$study_id))
gc_done$status <- as.numeric(gc_done$status) # must need this for Survsplit


#### Survsplit to get GC cumdose so far, also now need the TD confs: BP days so far ####
cb <- NULL
cb_study_idlist <- unique(casebase$study_id)
for(i in 1:length(cb_study_idlist)){
  casebasei <- casebase[casebase$study_id == cb_study_idlist[i], ]
  cbi_dtimesY <- casebasei$Stop
  cbi_dtimesYs <- casebasei$Begin
  gc_donei <- as.data.frame(gc_done[gc_done$study_id == cb_study_idlist[i], ])
  cbi_new <- survSplit(Surv(Start, End, status)~., gc_donei, cut = cbi_dtimesY, episode = "cycle")
  cbi_new <- cbi_new %>%
    filter(!is.na(Start)) %>%
    filter(End <= cbi_dtimesY[length(cbi_dtimesY)]) %>%
    group_by(cycle) %>%
    filter(row_number() == n()) %>%
    mutate(Start = ifelse(!(Start %in% cbi_dtimesYs), cbi_dtimesYs[cycle], Start))
  
  cb <- rbind(cb, cbi_new)
  if(i %% 10 == 0) print(i)
}
cb <- cb[!is.na(cb$Start), ]

#View(cb[cb$End != cb$Endcurious, ])
# Can merge in baseline long, but let's wait 
cb$futime <- sum(ftime)/ms                              # universal futime
cb$fiveinterval <- findInterval(cb$End, c(0, dtimes_V)) # TODO: replace cycle with fiveinterval

# label Visit
cbv <- cb 
cbv$days_to_servdate = ifelse(cbv$days_to_servdate < cbv$Start, Inf, cbv$days_to_servdate)

cbv <- cbv %>%
  group_by(study_id) %>%
  mutate(
    Dtrtp = lag(Dtrt)) %>%
  mutate(
    days_to_servdate = ifelse(is.na(days_to_servdate), Inf, days_to_servdate),
    Dtrtp = ifelse(is.na(Dtrtp), Dtrt, Dtrtp),
    Visit = ifelse(Dtrtp != Dtrt, 1, 0)) 
table(cbv$Visit) #V=1 male: ___; female: 

cbv$Visit <- ifelse(cbv$Dtrt == cbv$Dtrtp & is.finite(cbv$days_to_servdate)&
                      cbv$days_to_servdate <= cbv$End & cbv$days_to_servdate > cbv$Start,
                    1, cbv$Visit) 
table(cbv$Visit) #V=1 male: ___; female:

# No need because these values will be matched from discrete exposure dataset
#source("/users/edong/Code/calc_cumGC_cumBP_outcome.R", echo=T) 

###(no lag for LTC, because CTM outcome model)
ltcsubgender <- ltc_conf[ltc_conf$study_id %in% unique(cbv$study_id), ]
long_Y_cb <- cbv %>% #bp_cum_all %>%
  left_join(ltcsubgender) %>%
  mutate(Ltc = ifelse(days_to_ltcdate <= End, 1, 0))

# Label the Yevent outcome by joining with casebase, which is male only already
casebase$End <- casebase$Stop
long_Y_cb <- left_join(long_Y_cb, casebase[, c('study_id', 'End', 'Yevent')], by=c('study_id', 'End'))

# Replace cycle as five interval
long_Y_cb$cycle <- long_Y_cb$fiveinterval
long_Y_cb$cycle <- ifelse(long_Y_cb$End %% 5 == 0, long_Y_cb$cycle - 1, long_Y_cb$cycle)

# merge in df_Dinit male
df_Dinit <- df_Dinit[ , colnames(df_Dinit)[!colnames(df_Dinit) %in% c('X', 'X.1')]]
df_Dinitsubgender <- df_Dinit[df_Dinit$study_id %in% unique(cbv$study_id), ]
long_Y_cb <-  long_Y_cb[, c("study_id", colnames(long_Y_cb)[!colnames(long_Y_cb) %in% colnames(df_Dinit)])]
long_Y_cb <- left_join(long_Y_cb, df_Dinitsubgender )
long_Y_cb$Dtrt <- as.factor(long_Y_cb$Dtrt)
long_Y_cb$Dtrtp <- as.factor(long_Y_cb$Dtrtp)
long_Y_cb$D0 <- as.factor(long_Y_cb$D0)

# casebasecheck <- casebase[casebase$study_id %in% unique(long_Y_cb$study_id), ]
# dim(casebasecheck); table(casebasecheck$Yevent)
# dim(long_Y_cb); table(long_Y_cb$Yevent)
# View(long_Y_cb[is.na(long_Y_cb$Yevent), ])
long_Y_cb <- subset(long_Y_cb, select=-c(fiveinterval))
long_Y_cb$timebasis_epi <- bs(long_Y_cb$End, Boundary.knots=c(0, tlim), degree = 2)
long_Y_cb$A <- long_Y_cb$Dtrt


########################## Factorize and relevel #######################
facvec <- c("Dtrt",      "Dtrtp",  "Visit",      "Ltc",      "Yevent",          "D0",     "sex",
            "diabetes", "arthritis", "pfx",    "antisexhor", "inhalegc", "bronchodilators", "sexhor", "thiazide", 
            "pralcal",  "pdmab",     "falldr", "fallcomorb", "emerhosp", "totprevBPdaysc",  "status", "ikn", "A") 

long_Y_cb[, facvec] <- lapply(long_Y_cb[, facvec], as.factor); 

long_Y_cb$A <- relevel(long_Y_cb$A, ref="1"); long_Y_cb$D0 <- relevel(long_Y_cb$D0, ref="1"); 
long_Y_cb$Dtrt <- relevel(long_Y_cb$Dtrt, ref="1"); long_Y_cb$Dtrtp <- relevel(long_Y_cb$Dtrtp, ref="1");

########### Center all continuous vars in exposure data set ############
long_Y_cb$agescale <- as.numeric(scale(long_Y_cb$age, scale = FALSE))


###################### Calculate weights ######################
# Initial D Multinom marginal and conditional models
d0mproball <- predict(d0modm, newdata = long_Y_cb, type = "prob")
d0cproball <- predict(d0modc, newdata = long_Y_cb, type = "prob")

long_Y_cb$d0mprob <- d0mproball[cbind(1:nrow(d0mproball), long_Y_cb$D0)]
long_Y_cb$d0cprob <- d0cproball[cbind(1:nrow(d0cproball), long_Y_cb$D0)]

# Initial D ratio
long_Y_cb$d0p <- long_Y_cb$d0mprob/long_Y_cb$d0cprob; summary(long_Y_cb$d0p)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # baseline in numerator
# 1       1       1       1       1       1 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # baseline not in numerator
# 0.2418  0.8892  0.9670  0.9988  1.0770  4.1834

## V model 
V_gender$vmprob <- predict(vmodm, newdata = V_gender, type = "response")
V_gender$vcprob <- predict(vmodc, newdata = V_gender, type = "response")
V_gender$vmprob <- ifelse(V_gender$Visit == 0, 1-V_gender$vmprob, V_gender$vmprob)
V_gender$vcprob <- ifelse(V_gender$Visit == 0, 1-V_gender$vcprob, V_gender$vcprob)

V_gender <- V_gender %>%
  group_by(ikn) %>%
  mutate(cumprodm = cumprod(vmprob),
         cumprodc = cumprod(vcprob)) %>%
  mutate(vp = cumprodm/cumprodc)


## D model: Dtrtp = 0
dp0mproball <- predict(dmodp0m, newdata = Dtrtp0dat, type = "prob")
dp0cproball <- predict(dmodp0c, newdata = Dtrtp0dat, type = "prob")

Dtrtp0dat$dp0mprob <- dp0mproball[cbind(1:nrow(dp0mproball), Dtrtp0dat$Dtrt)]
Dtrtp0dat$dp0cprob <- dp0cproball[cbind(1:nrow(dp0cproball), Dtrtp0dat$Dtrt)]
Dtrtp0dat$dp <- Dtrtp0dat$dp0mprob/Dtrtp0dat$dp0cprob
Dtrtp0dat$dp0p <- Dtrtp0dat$dp0mprob/Dtrtp0dat$dp0cprob


## D model: Dtrtp > 0
dpnot0mproball <- predict(dmodpnot0m, newdata = Dtrtpnot0dat, type = "prob")
dpnot0cproball <- predict(dmodpnot0c, newdata = Dtrtpnot0dat, type = "prob")

Dtrtpnot0dat$dpnot0mprob <- dpnot0mproball[cbind(1:nrow(dpnot0mproball), Dtrtpnot0dat$Dtrt)]
Dtrtpnot0dat$dpnot0cprob <- dpnot0cproball[cbind(1:nrow(dpnot0cproball), Dtrtpnot0dat$Dtrt)]
Dtrtpnot0dat$dp <-Dtrtpnot0dat$dpnot0mprob/Dtrtpnot0dat$dpnot0cprob
Dtrtpnot0dat$dpnot0p <- Dtrtpnot0dat$dpnot0mprob/Dtrtpnot0dat$dpnot0cprob


# Append these two D data together before map back to V
Dtrtpdat <- rbind(Dtrtp0dat[, c('study_id', 'cycle', 'dp')], 
                  Dtrtpnot0dat[, c('study_id', 'cycle', 'dp')])

# Map back two D weights to V_gender
V_gender <- V_gender %>%
  full_join(Dtrtpdat) %>%
  mutate(dp = ifelse(is.na(dp), 1, dp)) %>% # This is enough, but for diagnostic purpose, we also map each D weight
  full_join(Dtrtp0dat[, c('study_id', 'cycle', 'dp0p')]) %>%
  mutate(dp0p = ifelse(is.na(dp0p), 1, dp0p)) %>%
  full_join(Dtrtpnot0dat[, c('study_id', 'cycle', 'dpnot0p')]) %>%
  mutate(dpnot0p = ifelse(is.na(dpnot0p), 1, dpnot0p))


# Lag the V D weights in the exposure data
V_gender <- V_gender %>%
  group_by(study_id) %>%
  mutate(vp.lag = lag(vp),
         dp.lag = lag(dp),
         dp0p.lag = lag(dp0p),
         dpnot0p.lag = lag(dpnot0p)) %>%
  mutate(vp.lag = ifelse(is.na(vp.lag), 1, vp.lag),
         dp.lag = ifelse(is.na(dp.lag), 1, dp.lag),
         dp0p.lag = ifelse(is.na(dp0p.lag), 1, dp0p.lag),
         dpnot0p.lag = ifelse(is.na(dpnot0p.lag), 1, dpnot0p.lag)
  )

# Combine weights in exposure dataset
V_gender$sw <- V_gender$vp.lag * V_gender$dp.lag; summary(V_gender$sw)


########### Map V weight to outcome data, also map in GC cumdose and BP cumdays ######
excludeindex <- setdiff(unique(V_gender$study_id), unique(long_Y_cb$study_id));length(excludeindex) # number of not cb sampled pt
long_Y_cb <- long_Y_cb %>%
  left_join( V_gender[!(V_gender$study_id %in% excludeindex), 
                      c('study_id', 'cycle', 'vp.lag', 'dp.lag', 'dp0p.lag', 'dpnot0p.lag', 'sw',
                        'gccumdose_out', 'cumBPdaysp'#, 'cumGCdays'
                      )]) %>%
  dplyr::rename(vp = vp.lag,
                dp = dp.lag,
                dp0p = dp0p.lag,
                dpnot0p = dpnot0p.lag,
                gccumdose = gccumdose_out)

summary(long_Y_cb$vp)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000    0.9829    0.9978    1.0920    1.0050 2255.0000 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. # baseline not in numerator 
# 0.0      0.9      1.0      5.4      1.1 433858.2 

summary(long_Y_cb$dp0p)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6085  1.0000  1.0000  1.0000  1.0000  1.8120 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   # baseline not in numerator
# 0.2348  1.0000  1.0000  1.0001  1.0000  3.1850 

summary(long_Y_cb$dpnot0p)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7964  1.0000  1.0000  1.0000  1.0000  1.2720 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # baseline not in numerator
# 0.4026  1.0000  1.0000  0.9998  1.0000  2.0380 

summary(long_Y_cb$dp)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6085  1.0000  1.0000  1.0000  1.0000  1.8120 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   # baseline not in numerator
# 0.2348  1.0000  1.0000  1.0001  1.0000  3.1850 

long_Y_cb$sw = long_Y_cb$sw * long_Y_cb$d0p
summary(long_Y_cb$sw)
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000    0.9828    0.9978    1.0900    1.0060    2255.0000 
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. # baseline not in numerator
# 0.0       0.8       1.0       9.0       1.1     807105.7 

qsl <- quantile(long_Y_cb$sw, prob=0.005)
qsu <- quantile(long_Y_cb$sw, prob=0.995)
long_Y_cb$sw.trunc <- ifelse(long_Y_cb$sw > qsu, qsu, long_Y_cb$sw)
long_Y_cb$sw.trunc <- ifelse(long_Y_cb$sw.trunc < qsl, qsl, long_Y_cb$sw.trunc)
summary(long_Y_cb$sw.trunc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.       # baseline in female
# 0.6354  0.9810  0.9969  1.0085  1.0056  3.0127 
#   Min. 1st Qu.  Median    Mean  3rd Qu.    Max.    # baseline not in numerator
# 0.1177  0.8357  0.9667  1.0295  1.1338  4.2156 

## Center continuous vars that enter the model
long_Y_cb$gccumdosescale <- as.numeric(scale(long_Y_cb$gccumdose, scale = FALSE))
long_Y_cb$cumBPdayspscale <- as.numeric(scale(long_Y_cb$cumBPdaysp, scale = FALSE))

# (do not center cumGCdays, because that is for calculating potential outcome cumdose) Actually no need GC cum days









