# rm(list=ls())
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim2/estfun-glm.q")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim2/infjack-glm.q")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim2/helper_sim2.R")
# library(survival)
# library(splines)
# library(dplyr)
# #
# # ### This part will be called in Marginal_sim-repeat.R; Only exist for testing this one round function ###
# nobs <- 1000
# seed = 2
# m <- 6

marginalnew_one_round <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs = nobs, states = states, trates = trates, tlim=tlim, seed = seed)
  times <- times_wide_generation(long, nobs)
  dtimes_A <- sort(unique(times$at[is.finite(times$at)]))  # Stop of treatment time
  
  ## Exposure Data: A as outcome, X as cov (delay one interval)
  long_A_ctm <- get_long_A_ctm(times, dtimes_A)
  
  #################### Poisson Exposure Model Fitting ################
  # Marginal
  model_poi_A_ctm_marginal <- glm((!Aexpo) ~ 1 + offset(log(futime)), family=poisson(link=log), data=long_A_ctm)
  coef(model_poi_A_ctm_marginal)
  
  # Conditional
  model_poi_A_ctm_conditional <- glm((!Aexpo) ~ Xconf + offset(log(futime)), family=poisson(link=log), data=long_A_ctm)
  coef(model_poi_A_ctm_conditional)
  
  ###################### Cox Exposure Model Fitting ##################
  # Marginal
  model_cox_A_ctm_marginal <- coxph(Surv(start, stop, (!Aexpo)) ~ 1, data=long_A_ctm, control = coxph.control(timefix=FALSE))
  coef(model_cox_A_ctm_marginal)

  # Conditional
  model_cox_A_ctm_conditional <- coxph(Surv(start, stop, (!Aexpo)) ~ Xconf, data=long_A_ctm, control = coxph.control(timefix=FALSE))
  coef(model_cox_A_ctm_conditional)
  
  #################### Poisson Exposure Model Weighting Functions ################
  hz <- function(x, xt=Inf) {exp(coef(model_poi_A_ctm_conditional)[1] + coef(model_poi_A_ctm_conditional)[2] * I(x >= xt))}
  shz <- function(x, xt=Inf) {exp(coef(model_poi_A_ctm_marginal)[1] + 0 * x)}

  ch <- function(x, xt=Inf) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
      xt_times <- sort(c(0, x[k], xt))
      for (i in 2:(2 + I(xt < x[k]))) {
        int[k] <- int[k] + integrate(hz, lower=xt_times[i-1], upper=xt_times[i], xt=xt)$value
      }
    }
    return(int)
  }
  sch <- function(x, xt=Inf) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
      xt_times <- sort(c(0, x[k], xt))
      for (i in 2:(2 + I(xt < x[k]))) {
        int[k] <- int[k] + integrate(shz, lower=xt_times[i-1], upper=xt_times[i], xt=xt)$value
      }
    }
    return(int)
  }
  
  #################### Cox Exposure Model Weights ###################
  # marginal
  Psi0s_cox <- basehaz(model_cox_A_ctm_marginal, centered = FALSE)
  Psi0s_cox <- Psi0s_cox[Psi0s_cox[, 2] %in% dtimes_A, ]
  Psi0s_cox <- rbind(c(0.0, 0.0, 0.0), Psi0s_cox)

  # conditional
  Psi0_cox <- basehaz(model_cox_A_ctm_conditional, centered = FALSE) # 1331 2
  Psi0_cox <- Psi0_cox[Psi0_cox[, 2] %in% dtimes_A, ]
  Psi0_cox <- rbind(c(0.0, 0.0, 0.0), Psi0_cox)
  
  
  
  ###################### Outcome Data: Y as outcome, X, A as cov (delay one interval) ##########################
  # dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))
  # 
  # long_Y_ctm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times 
  # long_Y_ctm$stop <- dtimes_Y                                          # automatic replicate dstops when finished. 
  # long_Y_ctm <- long_Y_ctm[order(long_Y_ctm$idx, long_Y_ctm$stop),]            # actually useless, already ordered
  # long_Y_ctm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  # long_Y_ctm$Xconf <- (long_Y_ctm$xt <= long_Y_ctm$stop)
  # long_Y_ctm$Aexpo <- (long_Y_ctm$at > long_Y_ctm$stop)        
  # long_Y_ctm$Yevent <- NA
  # long_Y_ctm$tlim <- pmin(long_Y_ctm$at, long_Y_ctm$yt, tlim) # rep at (of idx 1) 4.44 100 times
  # for (i in 1:nobs) {
  #   idx <- long_Y_ctm$idx == i
  #   long_Y_ctm$Yevent[idx] <- ifelse(times$yt[i] == long_Y_ctm$stop[idx], TRUE, FALSE)
  #   t <- pmin(long_Y_ctm$stop[idx], long_Y_ctm$tlim[idx])
  #   
  #   # Poisson Exposure Model 
  #   long_Y_ctm$zsurv[idx] <- exp(-ch(t, xt=times$xt[i]))
  #   long_Y_ctm$szsurv[idx] <- exp(-sch(t, xt=times$xt[i]))
  #   long_Y_ctm$zch[idx] <- ch(t, xt=times$xt[i])
  #   long_Y_ctm$szch[idx] <- sch(t, xt=times$xt[i])
  #   
  #   long_Y_ctm$zhz[idx] <- ifelse(times$at[i] <= t, hz(t, xt=times$xt[i]), 1.0)       
  #   long_Y_ctm$szhz[idx] <- ifelse(times$at[i] <= t, shz(t, xt=times$xt[i]), 1.0)
  #   long_Y_ctm$num[idx] <- long_Y_ctm$szsurv[idx] * long_Y_ctm$szhz[idx]
  #   long_Y_ctm$denom[idx] <- long_Y_ctm$zsurv[idx] * long_Y_ctm$zhz[idx]
  #   long_Y_ctm$iptw[idx] <- long_Y_ctm$num[idx]/long_Y_ctm$denom[idx]
  #   
  #   
  #   # Cox Exposure Model
  #   # Marginal weights
  #   long_Y_ctm$Psis_cox_cum[idx] <- Psi0s_cox$hazard[findInterval(t, Psi0s_cox$time)]
  #   j <- findInterval(times$at[i], Psi0s_cox$time)
  #   long_Y_ctm$Psis_cox_noncum[idx] <- ifelse(times$at[i] <= t, Psi0s_cox$hazard[j]-Psi0s_cox$hazard[max(1, j-1)], 1)
  #   long_Y_ctm$spc_cox[idx] <- exp(-long_Y_ctm$Psis_cox_cum[idx])
  #   
  #   # Conditional weights; as formula
  #   Psi_cox <- Psi0_cox
  #   Psi_cox$hazard <- Psi_cox$hazard - c(0, Psi0_cox$hazard[1:(nrow(Psi0_cox)-1)])
  #   Psi_cox$hazard <- ifelse(Psi_cox$time > times$xt[i], Psi_cox$hazard * exp(coef(model_cox_A_ctm_conditional)[1]), Psi_cox$hazard)
  #   Psi_cox$hazard <- cumsum(Psi_cox$hazard)
  #   
  #   long_Y_ctm$Psi_cox_cum[idx] <- Psi_cox$hazard[findInterval(t, Psi_cox$time)]
  #   long_Y_ctm$Psi_cox_noncum[idx] <- ifelse(times$at[i] <= t, Psi_cox$hazard[j]-Psi_cox$hazard[max(1, j-1)], 1)
  #   long_Y_ctm$pc_cox[idx] <- exp(-long_Y_ctm$Psi_cox_cum[idx])
  #   
  #   long_Y_ctm$sw_cox[idx] <- (long_Y_ctm$spc_cox[idx]*long_Y_ctm$Psis_cox_noncum[idx])/
  #     (long_Y_ctm$pc_cox[idx]*long_Y_ctm$Psi_cox_noncum[idx])
  #   
  #   
  #   # if (i %% 100 == 0)
  #   #     print(i)
  # }
  # long_Y_ctm <- subset(long_Y_ctm, stop <= yt)
  # long_Y_ctm$futime <- long_Y_ctm$stop - long_Y_ctm$start
  # # TODO: Time not on A
  # long_Y_ctm$timea <- as.numeric(!(long_Y_ctm$Aexpo)) * ifelse(is.finite(long_Y_ctm$at), long_Y_ctm$stop - long_Y_ctm$at, 0.0)
  # # TODO: Time on A
  # # long_Y_ctm$timea <- long_Y_ctm$Aexpo * ifelse(is.finite(long_Y_ctm$at), pmin(long_Y_ctm$stop, long_Y_ctm$at) - 0, long_Y_ctm$stop)
  # # TODO: Before A = 0 -> 1
  # # long_Y_ctm$timea <- long_Y_ctm$Aexpo * ifelse(is.finite(long_Y_ctm$at), long_Y_ctm$stop - long_Y_ctm$at, 0.0)
  # timebasisa_cox <- bs(long_Y_ctm$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
  # Case base
  # Case series dataset (outcome model):
  
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
  casebase$futime <- sum(ftime)/ms             # universal futime
  
  timebasis_cb <- bs(casebase$stop, Boundary.knots=c(0, tlim), degree=2)
  casebase$Xconf <- (casebase$xt <= casebase$stop)
  casebase$Aexpo <- (casebase$at > casebase$stop)
  
  casebase$tlim <- pmin(casebase$at, casebase$yt, tlim)
  for (i in 1:nobs) {
    idx <- casebase$id == i
    if (sum(idx) > 0) {
      t <- pmin(casebase$stop[idx], casebase$tlim[idx])
      
      casebase$zsurv[idx] <- exp(-ch(t, xt=times$xt[i]))
      casebase$szsurv[idx] <- exp(-sch(t, xt=times$xt[i]))
      casebase$zhz[idx] <- ifelse(times$at[i] <= t, hz(t, xt=times$xt[i]), 1.0)
      casebase$szhz[idx] <- ifelse(times$at[i] <= t, shz(t, xt=times$xt[i]), 1.0)
      casebase$num[idx] <- casebase$szsurv[idx] * casebase$szhz[idx]
      casebase$denom[idx] <- casebase$zsurv[idx] * casebase$zhz[idx]
      casebase$iptw[idx] <- casebase$num[idx]/casebase$denom[idx]
      
      
      # Cox Exposure Model
      # Marginal weights
      casebase$Psis_cox_cum[idx] <- Psi0s_cox$hazard[findInterval(t, Psi0s_cox$time)]
      j <- findInterval(times$at[i], Psi0s_cox$time)                                                    
      casebase$Psis_cox_noncum[idx] <- ifelse(times$at[i] <= t, Psi0s_cox$hazard[j]-Psi0s_cox$hazard[max(1, j-1)], 1)
      casebase$spc_cox[idx] <- exp(-casebase$Psis_cox_cum[idx])
      
      # Conditional weights; as formula
      Psi_cox <- Psi0_cox
      Psi_cox$hazard <- Psi_cox$hazard - c(0, Psi0_cox$hazard[1:(nrow(Psi0_cox)-1)])
      Psi_cox$hazard <- ifelse(Psi_cox$time > times$xt[i], Psi_cox$hazard * exp(coef(model_cox_A_ctm_conditional)[1]), Psi_cox$hazard)
      Psi_cox$hazard <- cumsum(Psi_cox$hazard)
      
      casebase$Psi_cox_cum[idx] <- Psi_cox$hazard[findInterval(t, Psi_cox$time)]
      casebase$Psi_cox_noncum[idx] <- ifelse(times$at[i] <= t, Psi_cox$hazard[j]-Psi_cox$hazard[max(1, j-1)], 1)
      casebase$pc_cox[idx] <- exp(-casebase$Psi_cox_cum[idx])
      
      casebase$sw_cox[idx] <- (casebase$spc_cox[idx]*casebase$Psis_cox_noncum[idx])/
        (casebase$pc_cox[idx]*casebase$Psi_cox_noncum[idx])
      
    }
    #print(i)
  }
  # TODO: Time of not on A
  casebase$timea <- as.numeric(!(casebase$Aexpo)) * ifelse(is.finite(casebase$at), casebase$stop - casebase$at , 0.0)
  # TODO: Time of on A
  # casebase$timea <- casebase$Aexpo * ifelse(is.finite(casebase$at), pmin(casebase$stop, casebase$at) - 0, casebase$stop)
  timebasisa_cb <- bs(casebase$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
  ################## Outcome Models #################
  pointest <- matrix(rep(NA,length(m)), nrow = 1)
  varest <- matrix(rep(NA,length(m)), nrow = 1)
  
  # ## Cox Outcome Models
  # long_Y_ctm$notA <- !(long_Y_ctm$Aexpo)
  # #coxout_poi_ctm <- coxph(Surv(start, stop, Yevent) ~ Aexpo + timebasisa_cox + cluster(idx), weight=iptw, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # coxout_poi_ctm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), weight=iptw, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # pointest[1] <- coef(coxout_poi_ctm)[1]
  # varest[1] <- vcov(coxout_poi_ctm)[1,1]
  # 
  # #coxout_cox_ctm <- coxph(Surv(start, stop, Yevent) ~ Aexpo + timebasisa_cox + cluster(idx), weight=sw_cox, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # coxout_cox_ctm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), weight=sw_cox, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # pointest[2] <- coef(coxout_cox_ctm)[1]
  # varest[2] <- vcov(coxout_cox_ctm)[1,1]
  # 
  # #coxout_unweight_ctm <- coxph(Surv(start, stop, Yevent) ~ Aexpo + timebasisa_cox + cluster(idx), data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # coxout_unweight_ctm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # pointest[3] <- coef(coxout_unweight_ctm)[1]
  # varest[3] <- vcov(coxout_unweight_ctm)[1,1]
  
  
  ## Case-base Outcome Models
  casebase$notA <- !(casebase$Aexpo)
  cbout_poi_ctm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime)), family=binomial(link=logit), weights=iptw, data=casebase)
  pointest[4] <- coef(cbout_poi_ctm)[4]
  varest[4] <- diag(infjack.glm(cbout_poi_ctm, casebase$idx))[4]
  
  cbout_cox_ctm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime)), family=binomial(link=logit), weights=sw_cox, data=casebase)
  pointest[5] <- coef(cbout_cox_ctm)[4]
  varest[5] <- diag(infjack.glm(cbout_cox_ctm, casebase$idx))[4]
  
  cbout_unweight_ctm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime)), family=binomial(link=logit), data=casebase)
  pointest[6] <- coef(cbout_unweight_ctm)[4]
  varest[6] <- diag(vcov(cbout_unweight_ctm))[4]
  
  
  # # A look at weights
  # summary(na.omit(long_Y_ctm$iptw))
  # summary(na.omit(long_Y_ctm$sw_cox))
  # summary(na.omit(casebase$iptw))
  # summary(na.omit(casebase$sw_cox))
  # #
  # par(mfrow = c(2,2))
  # plot(long_Y_ctm$stop, long_Y_ctm$iptw, main = "Poisson Exposure with Cox Outcome")
  # plot(long_Y_ctm$stop, long_Y_ctm$sw_cox, main = "Cox Exposure with Cox Outcome")
  # plot(casebase$stop, casebase$iptw, main = "Poisson Exposure with Case-base Outcome")
  # plot(casebase$stop, casebase$sw_cox, main = "Cox Exposure with Case-base Outcome")

  
  marginal_result <- list(pointest, varest)
  return(marginal_result)
}

#system.time(marginal_result <- marginalnew_one_round(seed, m, nobs, states, trates, tlim))




