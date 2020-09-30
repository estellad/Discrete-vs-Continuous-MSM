# rm(list=ls())
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_X_01/Sim1/helper_sim1_slm.R")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_X_01/Sim1/estfun-glm.q")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_X_01/Sim1/infjack-glm.q")
# library(survival)
# library(splines)
# library(dplyr)
# nobs <- 1000
# seed = 2
# m <- 4

marginalnew_one_round <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs = nobs, states = states, trates = trates, tlim=tlim, seed = seed)
  times <- times_wide_generation(long, nobs)
  dtimes_A <- seq(from = 0, to = tlim, length.out = 61)[-1]  # ~ 29.94 days interval
  
  ## Exposure Data: A as outcome, X as cov (delay one interval)
  long_A_slm <- get_long_A_slm(times, dtimes_A)  
  
  ###################### Logistic Exposure Model Fitting ################## 
  # Marginal
  timebasis_mo <- bs(long_A_slm$month, Boundary.knots = c(0, 73), degree = 2)
  model_logi_A_slm_marginal <- glm(Aexpo ~ 1 , family=binomial(link=logit), data=long_A_slm) # Marginal
  # coef(model_logi_A_slm_marginal)
  
  # Conditional
  model_logi_A_slm_conditional <- glm(Aexpo ~ Xconf, family=binomial(link=logit), data=long_A_slm) # Conditional
  # coef(model_logi_A_slm_conditional)
  
  long_A_slm$m_prob <- predict(model_logi_A_slm_marginal, newdata = long_A_slm, type = "response")  # marginal 
  long_A_slm$c_prob <- predict(model_logi_A_slm_conditional, newdata = long_A_slm, type = "response")  # conditional
  long_A_slm$m_prob <- ifelse(long_A_slm$Aexpo == 0, 1-long_A_slm$m_prob, long_A_slm$m_prob)
  long_A_slm$c_prob <- ifelse(long_A_slm$Aexpo == 0, 1-long_A_slm$c_prob, long_A_slm$c_prob)
  long_A_slm <- long_A_slm %>%
    group_by(idx) %>%
    mutate(cumprodm = cumprod(m_prob),
           cumprodc = cumprod(c_prob))

  ##################### Outcome Data: Y as outcome, A as cov ##########################
  dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))

  long_Y_slm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times
  long_Y_slm$stop <- dtimes_Y                                          # automatic replicate dstops when finished.
  long_Y_slm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  long_Y_slm$futime <- long_Y_slm$stop - long_Y_slm$start

  long_Y_slm$month <- findInterval(dtimes_Y, c(0,dtimes_A))
  long_Y_slm <- long_Y_slm[order(long_Y_slm$idx, long_Y_slm$stop),]            # actually useless, already ordered
  long_Y_slm$Xconf <- (long_Y_slm$xt <= long_Y_slm$stop)
  long_Y_slm$Aexpo <- (long_Y_slm$at <= long_Y_slm$stop)
  long_Y_slm$tlim <- pmin(long_Y_slm$at, long_Y_slm$yt, tlim) # rep at (of idx 1) 4.44 100 times

  for (i in 1:nobs) {

    # TODO: need to map over the time interval in the outcome data set
    idx <- long_Y_slm$idx == i
    idxa <- long_A_slm$idx == i
    
    long_Y_slm$Yevent[idx] <- ifelse(times$yt[i] == long_Y_slm$stop[idx], TRUE, FALSE)
    t <- pmin(long_Y_slm$stop[idx], long_Y_slm$tlim[idx])

    ## Logistic Exposure Model
    long_Y_slm$m_prob_cum[idx] <- as.numeric(unlist(long_A_slm[idxa,][long_Y_slm[long_Y_slm$idx == i,]$month,"cumprodm"]))
    long_Y_slm$c_prob_cum[idx] <- as.numeric(unlist(long_A_slm[idxa,][long_Y_slm[long_Y_slm$idx == i,]$month,"cumprodc"]))
    
    long_Y_slm$m_prob_cum[idx][which(is.na(long_Y_slm$m_prob_cum[idx]))] <- as.numeric(unlist(long_A_slm[idxa,'cumprodm']))[nrow(long_A_slm[idxa,])]
    long_Y_slm$c_prob_cum[idx][which(is.na(long_Y_slm$c_prob_cum[idx]))] <- as.numeric(unlist(long_A_slm[idxa,'cumprodc']))[nrow(long_A_slm[idxa,])]
    
    j <- findInterval(times$at[i], long_Y_slm$start[idx])
    long_Y_slm$sw_logi[idx] <- ifelse(times$at[i] <= t, long_Y_slm$m_prob_cum[idx][j]/long_Y_slm$c_prob_cum[idx][j],
                                      long_Y_slm$m_prob_cum[idx]/long_Y_slm$c_prob_cum[idx])

    # if (i %% 100 == 0)
    #     print(i)
  }
  long_Y_slm <- subset(long_Y_slm, start < yt)
  long_Y_slm$timea <- long_Y_slm$Aexpo * ifelse(is.finite(long_Y_slm$at), long_Y_slm$stop - long_Y_slm$at, 0.0)
  # long_Y_slm$futime <- long_Y_slm$stop - long_Y_slm$start
  timebasisa_cox <- bs(long_Y_slm$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
  ### Case base
  # Case series dataset (outcome model)
  cs <- times[is.finite(times$yt),]
  cs$stop <- cs$yt
  cs$Yevent <- 1
  
  ################### # Base series dataset (outcome model):
  
  ftime <- pmin(times$yt, tlim)
  ms <- nrow(cs) * 200
  # based on how much you contribute to the event time, persons determines how many times this individual will appear in the base dataset.
  persons <- rep(1:nobs, as.numeric(rmultinom(1, ms, ftime/sum(ftime)))) # length=20600 
  # based on the length of the event time of ith person, we randomly take a non-event time within this frame. 
  moments <- rep(NA, ms)                                                 # length=20600
  for (i in 1:ms) {
    moments[i] <- runif(1, min=0.0, max=ftime[persons[i]])
  }
  while(length(unique(moments)) < ms){
    for (i in 1:ms) {
      moments[i] <- runif(1, min=0.0, max=ftime[persons[i]])
    }
  }
  bs <- times[persons,]
  bs$stop <- moments
  bs$Yevent <- 0
  
  casebase <- rbind(cs, bs)
  casebase <- casebase[order(casebase$idx, casebase$stop),]
  casebase$Yevent <- as.logical(casebase$Yevent)
  for (i in 1:nobs) {
    idx <- casebase$idx == i
    if (sum(idx) > 0) {
      casebase$start[idx] <- c(0,casebase$stop[idx][1:(length(casebase$stop[idx])-1)])
    }
  }
  casebase$futime <- casebase$stop - casebase$start   
  casebase$month <- findInterval(casebase$stop, c(0,dtimes_A))
  
  timebasis_cb <- bs(casebase$stop, Boundary.knots=c(0, tlim), degree=2)
  casebase$Xconf <- (casebase$xt <= casebase$stop)
  casebase$Aexpo <- (casebase$at <= casebase$stop)
  
  
  casebase$tlim <- pmin(casebase$at, casebase$yt, tlim)
  for (i in 1:nobs) {
    idx <- casebase$idx == i
    idxa <- long_A_slm$idx == i
    if (sum(idx) > 0) {
      t <- pmin(casebase$stop[idx], casebase$tlim[idx])
      
      # Logistic Exposure Model
      casebase$m_prob_cum[idx] <- as.numeric(unlist(long_A_slm[idxa,][casebase[casebase$idx == i,]$month,"cumprodm"]))
      casebase$c_prob_cum[idx] <- as.numeric(unlist(long_A_slm[idxa,][casebase[casebase$idx == i,]$month,"cumprodc"]))
      
      casebase$m_prob_cum[idx][which(is.na(casebase$m_prob_cum[idx]))] <- as.numeric(unlist(long_A_slm[idxa,'cumprodm']))[nrow(long_A_slm[idxa,])]
      casebase$c_prob_cum[idx][which(is.na(casebase$c_prob_cum[idx]))] <- as.numeric(unlist(long_A_slm[idxa,'cumprodc']))[nrow(long_A_slm[idxa,])]
      
      j <- which(casebase$Aexpo[idx] == 1)[1]
      if(is.na(j)){
        casebase$sw_logi[idx] <- casebase$m_prob_cum[idx]/casebase$c_prob_cum[idx]
      }else{
        casebase$sw_logi[idx][1:(j-1)] <- (casebase$m_prob_cum[idx]/casebase$c_prob_cum[idx])[1:(j-1)]
        casebase$sw_logi[idx][j:length(which(idx))] <- casebase$m_prob_cum[idx][j]/casebase$c_prob_cum[idx][j]
      }
      
    }
    # if (i %% 100 == 0)
    #     print(i)
  }

  casebase$timea <- casebase$Aexpo * ifelse(is.finite(casebase$at), casebase$stop - casebase$at, 0.0)
  timebasisa_cb <- bs(casebase$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
  ################## Outcome Models #################
  pointest <- matrix(rep(NA,length(m)), nrow = 1)
  varest <- matrix(rep(NA,length(m)), nrow = 1)
  
  ## Cox Outcome Models
  coxout_logi_ctm <- coxph(Surv(start, stop, Yevent) ~ Aexpo + timebasisa_cox + cluster(idx), weight=sw_logi, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[1] <- coef(coxout_logi_ctm)[1]
  varest[1] <- vcov(coxout_logi_ctm)[1,1]

  coxout_unweight_ctm <- coxph(Surv(start, stop, Yevent) ~ Aexpo + timebasisa_cox + cluster(idx), data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[2] <- coef(coxout_unweight_ctm)[1]
  varest[2] <- vcov(coxout_unweight_ctm)[1,1]
  
  
  ## Case-base Outcome Models
  cbout_poi_ctm <- glm(Yevent ~ timebasis_cb + Aexpo + timebasisa_cb + offset(log(futime)), family=binomial(link=logit), weights=sw_logi, data=casebase)
  pointest[3] <- coef(cbout_poi_ctm)[4]
  varest[3] <- diag(infjack.glm(cbout_poi_ctm, casebase$idx))[4]
  
  cbout_unweight_ctm <- glm(Yevent ~ timebasis_cb + Aexpo + timebasisa_cb + offset(log(futime)), family=binomial(link=logit), data=casebase)
  pointest[4] <- coef(cbout_unweight_ctm)[4]
  varest[4] <- diag(vcov(cbout_unweight_ctm))[4]
  
  
  # # A look at weights
  # summary(long_Y_slm$sw_logi)
  # summary(casebase$sw_logi)
  # #
  # par(mfrow = c(1,2))
  # plot(long_Y_slm$stop, long_Y_slm$sw_logi, main = "Logistic Exposure with Cox Outcome")
  # plot(casebase$stop, casebase$sw_logi, main = "Logistic Exposure with Case-base Outcome")

  
  marginal_result <- list(pointest, varest)
  return(marginal_result)
}

#system.time(marginal_result <- marginalnew_one_round(seed, m, nobs, states, trates, tlim))




