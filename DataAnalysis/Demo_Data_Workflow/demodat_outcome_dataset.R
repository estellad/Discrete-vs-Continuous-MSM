###################################### Construct Outcome dataset #########################################
# Visit Conf Dtrt lagging? No because CTM (case-base)

# Case base
# Case series dataset (outcome model):
# length(dtimes_Y) = 2936, 
tlim = 365
cs <- gc_init[gc_init$status == 1, c("ikn",'yt')]
cs$stop <- cs$yt
cs$Yevent <- 1

dtimes_Y <- sort(unique(cs$yt))
################### # Base series dataset (outcome model):
ftime <- as.numeric(pmin(gc_init$yt, tlim))
ms <- nrow(cs) * 200
# TODO: persons = nrow(cs)*200
set.seed(2)
# based on how much you contribute to the event time, persons determines how many times this individual will appear in the base dataset.
persons <- rep(1:nobs, as.numeric(rmultinom(1, ms, ftime/sum(ftime))))
# based on the length of the event time of ith person, we randomly take a non-event time within this frame.
moments <- rep(NA, ms)
for (i in 1:ms) {
  moments[i] <- runif(1, min=0.0, max=ftime[persons[i]])
}
while(length(unique(moments)) < ms){
  for (i in 1:ms) {
    moments[i] <- runif(1, min=0.0, max=ftime[persons[i]])
  }
}

bs <- gc_init[persons, ]
bs$stop <- moments
bs$Yevent <- 0

casebase <- rbind(cs, bs)
casebase <- casebase[order(casebase$ikn, casebase$stop),]
casebase$Yevent <- as.logical(casebase$Yevent)

# Individual survsplit, for Dtrt and GC cum derivation (don't forget to lag one?)
cb <- NULL
for(i in 1:nobs){
  casebasei <- casebase[casebase$ikn == i, ]
  cbi_dtimesY <- casebasei$End
  cbi_dtimesYs <- casebasei$Start
  cbi_new <- survSplit(Surv(Begin, Stop, status) ~., gc_cum, cut=cbi_dtimesY, episode = "episode")
  cbi_new <- cbi_new %>%
    group_by(episode) %>%
    filter(row_number() == n()) %>%
    mutate(Begin = ifelse(!(Begin %in% Start), cbi_dtimesYs[episode], Begin))
  
  # label visit 
  
}

cb_bp <- NULL
for(i in 1:nobs){
  casebasei <- casebase[casebase$ikn == i, ]
  cbi_dtimesY <- casebasei$End
  cbi_dtimesYs <- casebasei$Start
  cbi_new <- survSplit(Surv(Begin, Stop, status) ~., gc_cum, cut=cbi_dtimesY, episode = "episode")
  cbi_new <- cbi_new %>%
    group_by(episode) %>%
    filter(row_number() == n()) %>%
    mutate(Begin = ifelse(!(Begin %in% Start), cbi_dtimesYs[episode], Begin))
  
  # label visit 
  
}


# Merge in baseline -> this can wait 
# for (i in 1:nobs) {
#   idx <- casebase$idx == i
#   if (sum(idx) > 0) {
#     casebase$start[idx] <- c(0,casebase$stop[idx][1:(length(casebase$stop[idx])-1)])
#     casebase$Dprev[idx] <- long$Dfrom[long$idx == i][1]
#   }
# }
