#source("~/Desktop/Thesis_Code/Analysis/fakedat_gccum_bpdays_exposure.R") # This sources generation
source("~/Desktop/Thesis_Code/Fully_Cont_MSM_MPP/DataAnalysis/Fake_Data_Workflow/fakedat_generation.R")
# this is for stop of treatment. 
gc <- get_gc(gc_final)
nobs <- length(unique(gc$ikn))
gc_init <- gc %>%         #Initial data set ->  This is my df_Dinit
  group_by(ikn) %>%
  filter(row_number() == 1)
gc_ikn_yt <- gc_init[, c("ikn","yt")]

gc_ikn <- gc_init %>% # a gc data set with only ikn
  select(ikn)
gc_ikn$ikn <- as.factor(gc_ikn$ikn)

#gc_cum <- get_gc_ikn_cum_servd(gc_final)

ltc_final <- get_ltc_final()
ltc_final$ikn <- as.factor(ltc_final$ikn)
ltc_conf <- full_join(gc_ikn, ltc_final, by="ikn") 
ltc_conf$Ltcdate <- ifelse(is.na(ltc_conf$Ltcdate), Inf, ltc_conf$Ltcdate)
ltc_conf <- ltc_conf[, c("ikn", "Ltcdate")]

baseline <- get_baseline()

############################# Only GC survsplit ###############################
split <- 73
dayinterval <- 5
dtimes_V <- seq(from = 0, to = 365, length.out = 74)[-1] # 5 days interval
dtimes_Vs <- c(0, dtimes_V[1:length(dtimes_V)-1])
df_new <- survSplit(Surv(Start, End, status) ~., gc, cut=dtimes_V, episode = "episode")
df_newsub <- df_new[, c("ikn", "Start", "End", "Dtrt", "episode", "daily", "Servdate")]
df_newsub$Dtrt <- as.numeric(df_newsub$Dtrt)


gc_full<- data.frame(
  ikn= as.factor(rep(1:nobs, times=rep(length(dtimes_V), nobs))),
  Start = rep(dtimes_Vs, nobs),
  End = rep(dtimes_V, nobs),
  Dtrt= 0,
  episode = rep(c(1:split), nobs),
  daily = 0
)
gc_full$Dtrt <- as.numeric(gc_full$Dtrt)

dat <- full_join(df_newsub, gc_full)
dat <- dat %>%
  arrange(ikn, episode, End, Dtrt) %>%
  full_join(gc_ikn_yt) %>%
  mutate(yt = as.numeric(yt))

########## Now take filter the last per ikn per episode, and the expo data is ready ########### 
# TODO (except cumulative GC dose not merged in): Servsplit cumdose, and them merge with full and merge in here? 
# Can't advance because not know how to add cumdose to full
dat_expo <- dat %>%
  group_by(ikn, episode) %>%
  filter(row_number() == n()) %>%
  mutate(Start = ifelse(!(Start %in% dtimes_Vs), dtimes_Vs[episode], Start),
         End = ifelse(End > yt, yt, End))  # TODO: if discrete don't cut end point, then delete this

dat_expo <- dat_expo[dat_expo$Start < dat_expo$yt, ]

######### Map in initial cumdose  
# take out initial cumdose and map to all
gc_cum_all_expo_init <- gc %>%  # gc
  group_by(ikn) %>%
  filter(row_number() == 1) %>%
  select(ikn, cumdose)

dat_expo_cum <- dat_expo %>%
  full_join(gc_cum_all_expo_init)

##### calculate the cumdose and cumdose_projection
gc_cum_cleaned <- dat_expo_cum %>%
  group_by(ikn) %>%
  mutate(gccuminterval = cumsum(daily*dayinterval), 
         gccumdose_projection = cumdose + gccuminterval, 
         gccumdose = lag(gccumdose_projection)) %>%
  mutate(gccumdose = ifelse(is.na(gccumdose), cumdose, gccumdose))



########### create a lag of Dtrt
dat_lagd <- gc_cum_cleaned %>%
  group_by(ikn) %>%
  mutate(Dtrtp = lag(Dtrt))
dat_lagd$Dtrtp <- ifelse(is.na(dat_lagd$Dtrtp), dat_lagd$Dtrt, dat_lagd$Dtrtp) # lag the time0 Dtrt for episode 1 Dtrtp

# TODO: cumdose already merged, only need servdate Servdate for Visit labeling. 
dat_lagd_cum <- dat_lagd %>%
  #full_join(gc_cum_cleaned) %>%
  mutate(Servdate = ifelse(Servdate<Start | Servdate == 0, NA, Servdate))
#dat_lagd_cum <- dat_lagd_cum[dat_lagd_cum$Start < dat_lagd_cum$yt, ]

dat_lagd_cum$Visit <- ifelse(dat_lagd_cum$Dtrt != dat_lagd_cum$Dtrtp, 1, 0)
# For those changed treatment at times within an interval, but stayed on the same dose
dat_lagd_cum$Visit <- ifelse(dat_lagd_cum$Dtrt == dat_lagd_cum$Dtrtp & !is.na(dat_lagd_cum$Servdate) &
                           dat_lagd_cum$Servdate <= dat_lagd_cum$End & dat_lagd_cum$Servdate > dat_lagd_cum$Start, # TODO: question about endpoints.
                         1, dat_lagd_cum$Visit) 
# dat_lagd_cum <- dat_lagd_cum %>%                   # Lead one for inclusion in the End interval, might need to do this for other covars start date. 
#   group_by(ikn) %>%
#   mutate(leadServdate = lead(Servdate)) 
# 
# # TODO: Question about end points
# dat_lagd_cum$Visit[which(dat_lagd_cum$leadServdate == dat_lagd_cum$End & dat_lagd_cum$Visit==0)] <- 1

dat_lagd_cum <- dat_lagd_cum %>%
  select(-Servdate, -gccuminterval, -cumdose) 
# TODO: wait after LTC and BP Merge in baseline
long_V_slm <- dat_lagd_cum %>%
  full_join(baseline)

# Merge in LTC
ltc_conf$ikn <- as.factor(ltc_conf$ikn)  # Somehow as factor is as integer, 
long_V_slm <- long_V_slm %>%
  full_join(ltc_conf)
  
long_V_slm$Ltcconf <- (long_V_slm$Ltcdate <= long_V_slm$End) # TODO: push one 
long_V_slm$Ltcconf <- 

# Merge in BP days so far.
# get bp and servsplit, exclude beyond yt date.
# LTC set to Inf because we model use label
# BP set to 0, because we are calculating days so far.
# merge in yt, modify end interval and delete
long_V_sub <- long_V_slm %>%
  select(ikn, episode)

bp_cum_all_cleaned <- right_join(bp_cum_all, long_V_sub)

long_V_slm <- long_V_slm %>%
  inner_join(bp_cum_all_cleaned, by=c("ikn", "episode")) %>%
  mutate(cumBPdays = ifelse(is.na(cumBPdays), 0, cumBPdays),
         cumBPdaysp = ifelse(is.na(cumBPdaysp), 0, cumBPdaysp))





#TODO: Stratify by sex
###################### Initial D model fitting ###########################
model_init_m <- multinom(D0 ~ diabetic, data=Ddata)

model_init_c <- multinom(D0 ~ diabetic + hormone, data=Ddata)

#################### Pooled Logisitic Visit Model Fitting ################
# Marginal
model_logi_V_slm_marginal <- glm(Visit ~  D0 + diabetic + hormone 
                                 + offset(log(futime)), family=binomial(link=logit), data=long_V_slm)
# coef(model_logi_V_slm_marginal)

# Conditional
model_logi_V_slm_conditional <- glm(Visit ~ Ltcconf + Dtrtp + D0 + diabetic + hormone + cumdose + episode
                                    + offset(log(futime)), family=binomial(link=logit), data=long_V_slm)
# coef(model_logi_V_slm_conditional)


# D model: Directly multinomial, because there is 0-dose with 3 days cut-off.
# No need to subset data.
#################### D Model fitting ################
Ddata <- long_V_slm[long_V_slm$Visit == 1, ]
# Marginal
model_logi_D_slm_marginal <- multinom(Dtrt ~ Dtrtp + D0, data=Ddata)
# coef(model_logi_D_slm_marginal)

# Conditional
model_logi_D_slm_conditional <- multinom(Dtrt ~ Ltcconf+Dtrtp + D0, data=Ddata)
#coef(model_logi_D_slm_conditional)

