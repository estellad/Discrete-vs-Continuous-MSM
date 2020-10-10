# rm(list=ls())
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/helper_sim4_new_slm.R")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/estfun-glm.q")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/infjack-glm.q")
# library(survival)
# library(splines)
# library(nnet)
# library(dplyr)
# logit <- function(p) {log(p)-log(1-p)}
# # ### This part will be called in Marginal_sim-repeat.R; Only exist for testing this one round function ###
# nobs <- 1000
# seed = 2
# m <- 8

marginalnew_one_round <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1 = trates1, tlim, seed)
  times <- times_wide_generation(long, nobs)
  dtimes_V <- seq(from = 0, to = tlim, length.out = 61)[-1] # monthly interval
  
  ## Exposure Data: V as outcome, X as cov (delay one interval)
  long_V_slm <- get_long_V_slm(times, dtimes_V, long)
  
  #################### Pooled Logisitic Visit Model Fitting ################
  # Marginal
  model_logi_V_slm_marginal <- glm(Visit ~ Dprev , family=binomial(link=logit), data=long_V_slm)
  # coef(model_logi_V_slm_marginal)
  
  # Conditional
  model_logi_V_slm_conditional <- glm(Visit ~ Xconf + Dprev, family=binomial(link=logit), data=long_V_slm)
  # coef(model_logi_V_slm_conditional)
  
  long_V_slm$m_prob <- predict(model_logi_V_slm_marginal, newdata = long_V_slm, type = "response")  # marginal 
  long_V_slm$c_prob <- predict(model_logi_V_slm_conditional, newdata = long_V_slm, type = "response")  # conditional
  long_V_slm$m_prob <- ifelse(long_V_slm$Visit == 0, 1-long_V_slm$m_prob, long_V_slm$m_prob)
  long_V_slm$c_prob <- ifelse(long_V_slm$Visit == 0, 1-long_V_slm$c_prob, long_V_slm$c_prob)
  long_V_slm <- long_V_slm %>%
    group_by(idx) %>%
    mutate(cumprodm = cumprod(m_prob),
           cumprodc = cumprod(c_prob))
  
  ################## Multinomial Logistic D Treatment Model Fitting: Changed Dose #############
  # Treatment Model after visit
  D_expo_data <- long_V_slm[long_V_slm$Visit == 1, ]
  
  # Marginal
  model_logi_D_slm_marginal <- multinom(Dtrt ~ Dprev, data=D_expo_data)
  # coef(model_logi_D_slm_marginal)
  
  # Conditional
  model_logi_D_slm_conditional <- multinom(Dtrt ~ Xconf+Dprev, data=D_expo_data)
  #coef(model_logi_D_slm_conditional)
  
  ###################### Outcome Data: Y as outcome, X, V as cov (delay one interval) ##########################
  dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))

  long_Y_slm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times
  long_Y_slm$stop <- dtimes_Y                                          # automatic replicate dstops when finished.
  long_Y_slm <- long_Y_slm[order(long_Y_slm$idx, long_Y_slm$stop),]            # actually useless, already ordered
  long_Y_slm$month <- findInterval(dtimes_Y, c(0,dtimes_V))
  long_Y_slm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  long_Y_slm$futime <- long_Y_slm$stop - long_Y_slm$start

  long_Y_slm$Xconf <- (long_Y_slm$xt <= long_Y_slm$stop)
  long_Y_slm$Visit <- (long_Y_slm$vt <= long_Y_slm$stop)   #keep V=0->1
  long_Y_slm$Yevent <- NA
  long_Y_slm$tlim <- pmin(long_Y_slm$vt, long_Y_slm$yt, tlim) # rep vt (of idx 1) 4.44 100 times

  long_Y_slm$Dprev <- NA
  for (i in 1:nobs){
    long_Y_slm$Dprev[long_Y_slm$idx == i] <- long$Dfrom[long$idx == i][1]
  }

  m_prob_d_all <- predict(model_logi_D_slm_marginal, long_Y_slm, type = "prob");
  c_prob_d_all <- predict(model_logi_D_slm_conditional, long_Y_slm, type = "prob");

  for (i in 1:nobs) {
    idx <- long_Y_slm$idx == i
    idxa <- long_V_slm$idx == i
    long_Y_slm$Yevent[idx] <- ifelse(times$yt[i] == long_Y_slm$stop[idx], TRUE, FALSE)
    t <- pmin(long_Y_slm$stop[idx], long_Y_slm$tlim[idx])

    ## Logistic Exposure Model
    long_Y_slm$m_prob_cum[idx] <- as.numeric(unlist(long_V_slm[idxa,][long_Y_slm[long_Y_slm$idx == i,]$month,"cumprodm"]))
    long_Y_slm$c_prob_cum[idx] <- as.numeric(unlist(long_V_slm[idxa,][long_Y_slm[long_Y_slm$idx == i,]$month,"cumprodc"]))
    
    long_Y_slm$m_prob_cum[idx][which(is.na(long_Y_slm$m_prob_cum[idx]))] <- as.numeric(unlist(long_V_slm[idxa,'cumprodm']))[nrow(long_V_slm[idxa,])]
    long_Y_slm$c_prob_cum[idx][which(is.na(long_Y_slm$c_prob_cum[idx]))] <- as.numeric(unlist(long_V_slm[idxa,'cumprodc']))[nrow(long_V_slm[idxa,])]
    j <- findInterval(times$vt[i], long_Y_slm$start[idx])
    long_Y_slm$sw_logi[idx] <- ifelse(times$vt[i] <= t, long_Y_slm$m_prob_cum[idx][j]/long_Y_slm$c_prob_cum[idx][j],
                                      long_Y_slm$m_prob_cum[idx]/long_Y_slm$c_prob_cum[idx])

    ##############################
    # Treatment D Model
    long_Y_slm$Dtrt[idx] <- ifelse(times$dt[i] <= t,
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]],
                                   long$Dfrom[long$idx == i][1])

    long_Y_slm$m_prob_d[idx] <- ifelse(times$dt[i] <= t,
                                       m_prob_d_all[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
                                       m_prob_d_all[idx, (1+long$Dfrom[long$idx == i][1])])
    long_Y_slm$c_prob_d[idx] <- ifelse(times$dt[i] <= t,
                                       c_prob_d_all[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
                                       c_prob_d_all[idx, (1+long$Dfrom[long$idx == i][1])])

    ## Logistic D Model
    j_d <- which(long_Y_slm$stop[idx] > times$dt[i])[1]
    if(!is.na(j_d)){
      long_Y_slm$sw_logi_d[idx] <- ifelse(times$dt[i] <= t, long_Y_slm$m_prob_d[idx][j_d]/long_Y_slm$c_prob_d[idx][j_d], # on and after visit
                                          #long_Y_slm$m_prob_d[idx][1]/long_Y_slm$c_prob_d[idx][1]
                                          1) # before visit
    }else{
      long_Y_slm$sw_logi_d[idx] <- 1 #long_Y_slm$m_prob_d[idx][1]/long_Y_slm$c_prob_d[idx][1] # before visit
    }
    # if (i %% 100 == 0)
    #     print(i)
    #
  }
  long_Y_slm <- subset(long_Y_slm, stop <= yt)
  long_Y_slm$Dtrt <- as.factor(long_Y_slm$Dtrt)
  long_Y_slm$Dtrt <- relevel(long_Y_slm$Dtrt, ref = "1")
  long_Y_slm$A <- long_Y_slm$Dtrt

  # TODO:
  # Combined weights
  long_Y_slm$sw_logilogi_a <- long_Y_slm$sw_logi * long_Y_slm$sw_logi_d

  # Time not A = 1, A=0
  long_Y_slm$timea0 <- as.numeric(long_Y_slm$A == 0) * ifelse(is.finite(long_Y_slm$vt), long_Y_slm$stop - long_Y_slm$vt, 0.0)
  timebasisa0_cox <- bs(long_Y_slm$timea0, Boundary.knots=c(0, tlim), degree=1)

  # Time on A=2
  # long_Y_slm$timea2 <- as.numeric(long_Y_slm$A == 2) * ifelse(is.finite(long_Y_slm$vt), long_Y_slm$stop - long_Y_slm$vt, 0.0)
  # long_Y_slm$timea2 <- ifelse(long_Y_slm$timea2 < 0, 0.0, long_Y_slm$timea2)
  # timebasisa2_cox <- bs(long_Y_slm$timea2, Boundary.knots=c(0, tlim), degree=1)
  
  # Case base
  # Case series dataset (outcome model):
  cs <- times[is.finite(times$yt),]
  cs$stop <- cs$yt
  cs$Yevent <- 1
  
  ################### # Base series dataset (outcome model):
  
  ftime <- pmin(times$yt, tlim)
  ms <- nrow(cs) * 200
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
      casebase$Dprev[idx] <- long$Dfrom[long$idx == i][1]
    }
  }
  casebase$futime <- casebase$stop - casebase$start                   # non-universal futime
  casebase$month <- findInterval(casebase$stop, c(0,dtimes_V))

  timebasis_cb <- bs(casebase$stop, Boundary.knots=c(0, tlim), degree=2)
  casebase$Xconf <- (casebase$xt <= casebase$stop)
  casebase$Visit <- (casebase$vt <= casebase$stop)  # V 0->1
  casebase$tlim <- pmin(casebase$vt, casebase$yt, tlim)
  
  # TODO: D model
  m_prob_d_all_cb <- predict(model_logi_D_slm_marginal, newdata = casebase, type = "prob")  # marginal
  c_prob_d_all_cb <- predict(model_logi_D_slm_conditional, newdata = casebase, type = "prob")  # conditional
  
  for (i in 1:nobs) {
    idx <- casebase$id == i
    idxa <- long_V_slm$idx == i
    if (sum(idx) > 0) {
      t <- pmin(casebase$stop[idx], casebase$tlim[idx])
      # Logistic Exposure Model
      casebase$m_prob_cum[idx] <- as.numeric(unlist(long_V_slm[idxa,][casebase[casebase$idx == i,]$month,"cumprodm"]))
      casebase$c_prob_cum[idx] <- as.numeric(unlist(long_V_slm[idxa,][casebase[casebase$idx == i,]$month,"cumprodc"]))
      
      casebase$m_prob_cum[idx][which(is.na(casebase$m_prob_cum[idx]))] <- as.numeric(unlist(long_V_slm[idxa,'cumprodm']))[nrow(long_V_slm[idxa,])]
      casebase$c_prob_cum[idx][which(is.na(casebase$c_prob_cum[idx]))] <- as.numeric(unlist(long_V_slm[idxa,'cumprodc']))[nrow(long_V_slm[idxa,])]
      j <- which(casebase$stop[idx] > times$vt[i])[1]
      if(is.na(j)){
        casebase$sw_logi[idx] <- casebase$m_prob_cum[idx]/casebase$c_prob_cum[idx]
      }else{
        casebase$sw_logi[idx][1:(j-1)] <- (casebase$m_prob_cum[idx]/casebase$c_prob_cum[idx])[1:(j-1)]
        casebase$sw_logi[idx][j:length(which(idx))] <- casebase$m_prob_cum[idx][j]/casebase$c_prob_cum[idx][j]
      }
      
      ##############################
      # Treatment D Model
      casebase$Dtrt[idx] <- ifelse(times$dt[i] <= t,
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]],
                                   long$Dfrom[long$idx == i][1])
      
      casebase$m_prob_d[idx] <- ifelse(times$dt[i] <= t,
                                       m_prob_d_all_cb[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
                                       m_prob_d_all_cb[idx, (1+long$Dfrom[long$idx == i][1])])
      casebase$c_prob_d[idx] <- ifelse(times$dt[i] <= t,
                                       c_prob_d_all_cb[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
                                       c_prob_d_all_cb[idx, (1+long$Dfrom[long$idx == i][1])])
      
      ## Logistic D Model
      j_d <- which(casebase$stop[idx] > times$dt[i])[1]
      if(!(is.na(j_d))){
        casebase$sw_logi_d[idx] <- ifelse(times$dt[i] <= t, casebase$m_prob_d[idx][j_d]/casebase$c_prob_d[idx][j_d],
                                          1)
      }else{
        casebase$sw_logi_d[idx] <- 1
      }
    }
    
    #print(i)
  }
  
  # TODO:
  # Combined Treatment A
  casebase$Dtrt <- as.factor(casebase$Dtrt)
  casebase$Dtrt <- relevel(casebase$Dtrt, ref = "1")
  casebase$A <- casebase$Dtrt
  
  
  # Combined weights sw_poi_a & sw_cox_a
  casebase$sw_logilogi_a <- casebase$sw_logi * casebase$sw_logi_d
  
  # Time on A=0
  casebase$timea0 <- as.numeric(casebase$A == 0) * ifelse(is.finite(casebase$vt), casebase$stop - casebase$vt, 0.0)
  timebasisa0_cb <- bs(casebase$timea0, Boundary.knots=c(0, tlim), degree=1)
  
  # Time on A=2
  # casebase$timea2 <- as.numeric(casebase$A == 2) * ifelse(is.finite(casebase$vt), casebase$stop - casebase$vt, 0.0)
  # casebase$timea2 <- ifelse(casebase$timea2 < 0, 0.0, casebase$timea2)
  # timebasisa2_cb <- bs(casebase$timea2, Boundary.knots=c(0, tlim), degree=1)
  
  
  ################## Outcome Models #################
  pointest <- matrix(rep(NA,length(m)), nrow = 1)
  varest <- matrix(rep(NA,length(m)), nrow = 1)
  pointest2 <- matrix(rep(NA,length(m)), nrow = 1)
  varest2 <- matrix(rep(NA,length(m)), nrow = 1)
  
  ## Cox Outcome on A Model
  # Only V Cox model
  coxout_Vlogi_slm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox + Dprev #+ timebasisa2_cox
                           + cluster(idx), weight=sw_logi, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  # Only D Logistic model
  coxout_Dlogi_slm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox + Dprev #+ timebasisa2_cox
                            + cluster(idx), weight=sw_logi_d, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  # V D Combined Poisson * Logistic Model
  coxout_VlogiDlogi_slm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox + Dprev #+ timebasisa2_cox
                                + cluster(idx), weight=sw_logilogi_a, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  # Unweighted
  coxout_unweight_slm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox + Dprev #+ timebasisa2_cox
                               + cluster(idx), data=long_Y_slm, control = coxph.control(timefix=FALSE))

  # For parameter A=0
  pointest[1] <- coef(coxout_Vlogi_slm)[1]
  varest[1] <- vcov(coxout_Vlogi_slm)[1,1]
  pointest[2] <- coef(coxout_Dlogi_slm)[1]
  varest[2] <- vcov(coxout_Dlogi_slm)[1,1]
  pointest[3] <- coef(coxout_VlogiDlogi_slm)[1]
  varest[3] <- vcov(coxout_VlogiDlogi_slm)[1,1]
  pointest[4] <- coef(coxout_unweight_slm)[1]
  varest[4] <- vcov(coxout_unweight_slm)[1,1]

  # For parameter A=2
  pointest2[1] <- coef(coxout_Vlogi_slm)[2]
  varest2[1] <- vcov(coxout_Vlogi_slm)[2,2]
  pointest2[2] <- coef(coxout_Dlogi_slm)[2]
  varest2[2] <- vcov(coxout_Dlogi_slm)[2,2]
  pointest2[3] <- coef(coxout_VlogiDlogi_slm)[2]
  varest2[3] <- vcov(coxout_VlogiDlogi_slm)[2,2]
  pointest2[4] <- coef(coxout_unweight_slm)[2]
  varest2[4] <- vcov(coxout_unweight_slm)[2,2]
  
  
  # Case-base Outcome Models
  # Only V cox
  cbout_Vlogi_slm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb + Dprev #+ timebasisa2_cb 
                        + offset(log(futime)), family=binomial(link=logit), weights=sw_logi, data=casebase)
  # Only D logistic
  cbout_Dlogi_slm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb + Dprev #+ timebasisa2_cb 
                         + offset(log(futime)), family=binomial(link=logit), weights=sw_logi_d, data=casebase)
  # V * D combined: poisson * logistic
  cbout_VlogiDlogi_slm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb + Dprev #+ timebasisa2_cb 
                             + offset(log(futime)), family=binomial(link=logit), weights=sw_logilogi_a, data=casebase)
  # Unweighted
  cbout_unweight_slm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb + Dprev #+ timebasisa2_cb 
                            + offset(log(futime)), family=binomial(link=logit), data=casebase)
  
  # For parameter A=0
  pointest[5] <- coef(cbout_Vlogi_slm)[4]
  varest[5] <- diag(infjack.glm(cbout_Vlogi_slm, casebase$idx))[4]
  pointest[6] <- coef(cbout_Dlogi_slm)[4]
  varest[6] <- diag(infjack.glm(cbout_Dlogi_slm, casebase$idx))[4]
  pointest[7] <- coef(cbout_VlogiDlogi_slm)[4]
  varest[7] <- diag(infjack.glm(cbout_VlogiDlogi_slm, casebase$idx))[4]
  pointest[8] <- coef(cbout_unweight_slm)[4]
  varest[8] <- diag(vcov(cbout_unweight_slm))[4]
  
  # For parameter A=2
  pointest2[5] <- coef(cbout_Vlogi_slm)[5]
  varest2[5] <- diag(infjack.glm(cbout_Vlogi_slm, casebase$idx))[5]
  pointest2[6] <- coef(cbout_Dlogi_slm)[5]
  varest2[6] <- diag(infjack.glm(cbout_Dlogi_slm, casebase$idx))[5]
  pointest2[7] <- coef(cbout_VlogiDlogi_slm)[5]
  varest2[7] <- diag(infjack.glm(cbout_VlogiDlogi_slm, casebase$idx))[5]
  pointest2[8] <- coef(cbout_unweight_slm)[5]
  varest2[8] <- diag(vcov(cbout_unweight_slm))[5]
  
  # # A look at weights
  # summary(long_Y_slm$sw_logi)
  # summary(long_Y_slm$sw_logi_d)
  # summary(long_Y_slm$sw_logilogi_a)
  # 
  # summary(casebase$sw_logi)
  # summary(casebase$sw_logi_d)
  # summary(casebase$sw_logilogi_a)
  # #
  # par(mfrow = c(2,2))
  # plot(long_Y_slm$stop, long_Y_slm$sw_logi, main = "Logistic Exposure with Cox Outcome")
  # plot(casebase$stop, casebase$sw_logi, main = "Logistic Exposure with Case-base Outcome")
  # plot(long_Y_slm$stop, long_Y_slm$sw_logilogi_a, main = "Logistic Exposure with Cox Outcome")
  # plot(casebase$stop, casebase$sw_logilogi_a, main = "Logistic Exposure with Case-base Outcome")
  
  marginal_result <- list(pointest, varest, pointest2, varest2)
  return(marginal_result)
}
#system.time(marginal_result <- marginalnew_one_round(seed, m, nobs, states, trates, tlim))
# # 
# marginal_result[1]
# marginal_result[3]


