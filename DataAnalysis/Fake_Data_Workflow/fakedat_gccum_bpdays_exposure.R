# Servsplit Cumdose
source("~/Desktop/Thesis_Code/Fully_Cont_MSM_MPP/DataAnalysis/Fake_Data_Workflow/fakedat_generation.R")
nobs <- 5
gc_cum <- get_gc_ikn_cum_servd(gc_final)

gc <- get_gc(gc_final)
gc_init <- gc %>%         #Initial data set
  group_by(ikn) %>%
  filter(row_number() == 1)
gc_ikn_yt <- gc_init[, c("ikn","yt")]

split <- 73
dayinterval <- 5
dtimes_V <- seq(from = 0, to = 365, length.out = 74)[-1] # 5 days interval
dtimes_Vs <- c(0, dtimes_V[1:length(dtimes_V)-1])
gc_cum_new <- survSplit(Surv(Start, End, status) ~., gc_cum, cut=dtimes_V, episode = "episode")
gc_cum_new$ikn <- as.factor(gc_cum_new$ikn)
gc_cum_new$Dtrt <- as.factor(gc_cum_new$Dtrt)

gc_full<- data.frame(
  ikn= as.factor(rep(1:nobs, times=rep(length(dtimes_V), nobs))),
  Start = rep(dtimes_Vs, nobs),
  End = rep(dtimes_V, nobs),
  Dtrt= 0,
  episode = rep(c(1:split), nobs),
  daily = 0
)
gc_full$Dtrt <- as.factor(gc_full$Dtrt)

gc_cum_all <- full_join(gc_cum_new, gc_full)
gc_cum_all <- gc_cum_all %>%
  arrange(ikn, episode, End, Dtrt) %>%
  full_join(gc_ikn_yt) %>%
  mutate(yt = as.numeric(yt))

gc_cum_all_expo <- gc_cum_all %>%
  group_by(ikn, episode) %>%
  filter(row_number() == n()) %>%
  select(-status) %>%
  mutate(Start = ifelse(!(Start %in% dtimes_Vs), dtimes_Vs[episode], Start),
         End = ifelse(End > yt, yt, End))  # TODO: if discrete don't cut end point, then delete this

gc_cum_all_expo <- gc_cum_all_expo[gc_cum_all_expo$Start < gc_cum_all_expo$yt, ]

# TODO: here need to correct sub
# TODO: instead of take out from gc_cum_all_expo, take from gc_cum/gc (include cumdose in gc)
# take out initial cumdose and map to all
gc_cum_all_expo_init <- gc_cum %>%
  group_by(ikn) %>%
  filter(row_number() == 1) %>%
  select(ikn, cumdose)

gc_cum_all_expo <- gc_cum_all_expo %>%
  select(-cumdose, -cumdose_projection) %>%
  full_join(gc_cum_all_expo_init)

gc_cum_cleaned <- gc_cum_all_expo %>%
  group_by(ikn) %>%
  mutate(gccuminterval = cumsum(daily*dayinterval),
         gccumdose_projection = cumdose + gccuminterval,
         gccumdose = lag(gccumdose_projection)) %>%
  mutate(gccumdose = ifelse(is.na(gccumdose), cumdose, gccumdose))


# gc_cum_cleaned <- NULL
# for(i in 1:nobs){
#   gc_cumi <- gc_cum_all_expo[gc_cum_all_expo$ikn == i, ]
#   gc_cumi$cumdose_projection[1] <- gc_cumi$cumdose[1] + gc_cumi$daily[1] * dayinterval
#   for(j in 2:nrow(gc_cumi)){
#     gc_cumi$cumdose[j] <- gc_cumi$cumdose_projection[j-1]
#     gc_cumi$cumdose_projection[j] <- gc_cumi$cumdose[j] + gc_cumi$daily[j] * dayinterval
#   }
#   gc_cum_cleaned <- rbind(gc_cum_cleaned, gc_cumi)
# }
gc_cum_cleaned <- gc_cum_cleaned %>%
  select(ikn, Start, End, gccumdose, gccumdose_projection, episode, Dtrt, Servdate,
         daily)

# Now can merge back to gc_expose !!!!
gc_cum_cleaned$Dtrt <- as.numeric(gc_cum_cleaned$Dtrt)

#################################### BP ########################################
bp_conf <- get_bp_days_final()
bp_conf$Start <- bp_conf$BPdate
bp_conf$End <- bp_conf$Start + bp_conf$duration
nobsbp <- length(unique(bp_conf$ikn))

split <- 73
dayinterval <- 5
dtimes_V <- seq(from = 0, to = 365, length.out = 74)[-1] # 5 days interval
dtimes_Vs <- c(0, dtimes_V[1:length(dtimes_V)-1])
bp_cum_new <- survSplit(Surv(Start, End, status) ~., bp_conf, cut=dtimes_V, episode = "episode")
bp_cum_new$ikn <- as.factor(bp_cum_new$ikn)

bp_cum_new <- bp_cum_new %>%
  mutate(Duration = End - Start) %>%
  group_by(ikn, episode) %>% 
  summarise(BPdays = sum(Duration))

bp_full<- data.frame(
  ikn= as.factor(rep(1:nobsbp, times=rep(length(dtimes_V), nobsbp))),
  episode = rep(c(1:split), nobsbp),
  BPdays = 0
)

bp_cum_all <- full_join(bp_cum_new, bp_full)
bp_cum_all <- bp_cum_all %>%
  arrange(ikn, episode, BPdays) %>%
  group_by(ikn, episode) %>%
  filter(row_number() == n())%>%
  group_by(ikn) %>%
  mutate(cumBPdays = cumsum(BPdays), 
         cumBPdaysp = lag(cumBPdays)) %>% # cumBPdays at k-1, would use for modeling
  mutate(cumBPdaysp = ifelse(is.na(cumBPdaysp), 0, cumBPdaysp)) %>% # fix first as 0
  select(ikn, episode, cumBPdays, cumBPdaysp)
  

########## After cleaning existing, Merge in full long_V_slm ikn list, and clean those with 0 BP ###########


