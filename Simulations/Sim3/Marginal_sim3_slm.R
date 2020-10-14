# rm(list=ls())
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/helper_sim3_slm.R")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/estfun-glm.q")
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/infjack-glm.q")
# library(survival)
# library(splines)
# logit <- function(p) {log(p)-log(1-p)}
# # ### This part will be called in Marginal_sim-repeat.R; Only exist for testing this one round function ###
# nobs <- 1000
# seed = 2
# m <- 8

marginalnew_one_round <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1=trates1, tlim, seed)
  times <- times_wide_generation(long, nobs)
  dtimes_V <- seq(from = 0, to = tlim, length.out = 61)[-1] # monthly interval
  
  ## Exposure Data: V as outcome, X as cov (delay one interval)
  long_V_slm <- get_long_V_slm(times, dtimes_V, long)
  
  #################### Logistic Exposure Model Fitting ################
  # Marginal
  # Exposure model should invert and in outcome model also invert V
  model_logi_V_slm_marginal <- glm(Visit ~ 1 , family=binomial(link=logit), data=long_V_slm) # Marginal
  coef(model_logi_V_slm_marginal)
  
  # Conditional
  model_logi_V_slm_conditional <- glm(Visit ~ Xconf, family=binomial(link=logit), data=long_V_slm) # Conditional
  coef(model_logi_V_slm_conditional)
  
  long_V_slm$m_prob <- predict(model_logi_V_slm_marginal, newdata = long_V_slm, type = "response")  # marginal 
  long_V_slm$c_prob <- predict(model_logi_V_slm_conditional, newdata = long_V_slm, type = "response")  # conditional
  long_V_slm$m_prob <- ifelse(long_V_slm$Visit == 0, 1-long_V_slm$m_prob, long_V_slm$m_prob)
  long_V_slm$c_prob <- ifelse(long_V_slm$Visit == 0, 1-long_V_slm$c_prob, long_V_slm$c_prob)
  long_V_slm <- long_V_slm %>%
    group_by(idx) %>%
    mutate(cumprodm = cumprod(m_prob),
           cumprodc = cumprod(c_prob))
  
  ################## Logistic D Treatment Model Fitting #############
  # Treatment Model
  D_expo_data <- long_V_slm[long_V_slm$Visit == 1, ]
  
  # Marginal
  model_logi_D_slm_marginal <- glm(Dtrt ~ 1, family=binomial(link=logit), data=D_expo_data)
  # coef(model_logi_D_slm_marginal)
  
  # Conditional                        
  model_logi_D_slm_conditional <- glm(Dtrt ~ Xconf, family=binomial(link=logit), data=D_expo_data)
  # coef(model_logi_D_slm_conditional)
  
  # ###################### Outcome Data: Y as outcome, X, V as cov (delay one interval) ##########################
  dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))

  long_Y_slm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times
  long_Y_slm$stop <- dtimes_Y                                          # automatic replicate dstops when finished.
  long_Y_slm <- long_Y_slm[order(long_Y_slm$idx, long_Y_slm$stop),]            # actually useless, already ordered
  long_Y_slm$month <- findInterval(dtimes_Y, c(0,dtimes_V))
  long_Y_slm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  long_Y_slm$futime <- long_Y_slm$stop - long_Y_slm$start

  long_Y_slm$Xconf <- (long_Y_slm$xt <= long_Y_slm$stop)
  long_Y_slm$Visit <- (long_Y_slm$vt <= long_Y_slm$stop)
  long_Y_slm$Yevent <- NA
  long_Y_slm$tlim <- pmin(long_Y_slm$vt, long_Y_slm$yt, tlim) # rep vt (of idx 1) 4.44 100 times

  # TODO: D model
  long_Y_slm$m_prob_d <- predict(model_logi_D_slm_marginal, newdata = long_Y_slm, type = "response")  # marginal
  long_Y_slm$c_prob_d <- predict(model_logi_D_slm_conditional, newdata = long_Y_slm, type = "response")  # conditional


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
    long_Y_slm$Dtrt[idx] <- ifelse(times$dt[i] <= t, long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 1)
    # TODO: Here ifelse Dtrt == 1 for the probs, then parameter near 1.
    long_Y_slm$m_prob_d[idx] <- ifelse(long_Y_slm$Dtrt[idx] == 0, 1-long_Y_slm$m_prob_d[idx], long_Y_slm$m_prob_d[idx])
    long_Y_slm$c_prob_d[idx] <- ifelse(long_Y_slm$Dtrt[idx] == 0, 1-long_Y_slm$c_prob_d[idx], long_Y_slm$c_prob_d[idx])

    ## Logistic D Model
    j_d <- which(long_Y_slm$stop[idx] > times$dt[i])[1]
    if(!is.na(j_d)){
      long_Y_slm$sw_logi_d[idx] <- ifelse(times$dt[i] <= t, long_Y_slm$m_prob_d[idx][j_d]/long_Y_slm$c_prob_d[idx][j_d],
                                          1)
    }else{
      long_Y_slm$sw_logi_d[idx] <- 1
    }
    # if (i %% 100 == 0)
    #     print(i)
    #
  }
  long_Y_slm <- subset(long_Y_slm, stop <= yt)
  long_Y_slm$A <- ifelse(long_Y_slm$Dtrt == 1, 1, 0)

  # Combined weights two logistic models
  long_Y_slm$sw_logilogi_a <- long_Y_slm$sw_logi * long_Y_slm$sw_logi_d

  # Time not on A
  long_Y_slm$timea <- as.numeric(!(long_Y_slm$A)) * ifelse(is.finite(long_Y_slm$vt), long_Y_slm$stop - long_Y_slm$vt, 0.0)
  timebasisa_cox <- bs(long_Y_slm$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
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
    }
  }
  casebase$futime <- sum(ftime)/ms             # universal futime
  casebase$month <- findInterval(casebase$stop, c(0,dtimes_V))
  
  timebasis_cb <- bs(casebase$stop, Boundary.knots=c(0, tlim), degree=2)
  casebase$Xconf <- (casebase$xt <= casebase$stop)
  casebase$Visit <- (casebase$vt <= casebase$stop)  
  casebase$tlim <- pmin(casebase$vt, casebase$yt, tlim)
  
  # TODO: D model
  casebase$m_prob_d <- predict(model_logi_D_slm_marginal, newdata = casebase, type = "response")  # marginal
  casebase$c_prob_d <- predict(model_logi_D_slm_conditional, newdata = casebase, type = "response")  # conditional
  
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
      casebase$Dtrt[idx] <- ifelse(times$dt[i] <= t, long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 1)
      # TODO: Same here ifelse(Dtrt == 1, 1-prob, prob) works. Currently ifelse Dtrt == 0
      casebase$m_prob_d[idx] <- ifelse(casebase$Dtrt[idx] == 0, 1-casebase$m_prob_d[idx], casebase$m_prob_d[idx])
      casebase$c_prob_d[idx] <- ifelse(casebase$Dtrt[idx] == 0, 1-casebase$c_prob_d[idx], casebase$c_prob_d[idx])   # now flip the probability?
      
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
  casebase$A <- ifelse(casebase$Dtrt == 1, 1, 0)
  # Combined weights two logistic models 
  casebase$sw_logilogi_a <- casebase$sw_logi * casebase$sw_logi_d
  
  # Time not on A
  casebase$timea <- as.numeric(!(casebase$A)) * ifelse(is.finite(casebase$vt), casebase$stop - casebase$vt, 0.0)
  timebasisa_cb <- bs(casebase$timea, Boundary.knots=c(0, tlim), degree=1)
  
  
  ################## Outcome Models #################
  pointest <- matrix(rep(NA,length(m)), nrow = 1)
  varest <- matrix(rep(NA,length(m)), nrow = 1)
  
  ## Cox Outcome Models
  # Flip the outcome model:
  long_Y_slm$notA <- as.numeric(!(long_Y_slm$A))
  # V only
  coxout_Vlogi_slm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), weight=sw_logi, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[1] <- coef(coxout_Vlogi_slm)[1]
  varest[1] <- vcov(coxout_Vlogi_slm)[1,1]

  # D only
  coxout_Dlogi_slm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), weight=sw_logi_d, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[2] <- coef(coxout_Dlogi_slm)[1]
  varest[2] <- vcov(coxout_Dlogi_slm)[1,1]

  # V*D combined
  coxout_VlogiDlogi_slm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), weight=sw_logilogi_a, data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[3] <- coef(coxout_VlogiDlogi_slm)[1]
  varest[3] <- vcov(coxout_VlogiDlogi_slm)[1,1]

  # unweighted
  coxout_unweight_slm <- coxph(Surv(start, stop, Yevent) ~ notA + timebasisa_cox + cluster(idx), data=long_Y_slm, control = coxph.control(timefix=FALSE))
  pointest[4] <- coef(coxout_unweight_slm)[1]
  varest[4] <- vcov(coxout_unweight_slm)[1,1]
  
  
  ## Case-base Outcome Models
  casebase$notA <- as.numeric(!(casebase$A))
  # V only
  cbout_Vlogi_slm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime))
                         , family=binomial(link=logit), weights=sw_logi, data=casebase)
  pointest[5] <- coef(cbout_Vlogi_slm)[4]
  varest[5] <- diag(infjack.glm(cbout_Vlogi_slm, casebase$idx))[4]
  
  # D only
  cbout_Dlogi_slm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime))
                         , family=binomial(link=logit), weights=sw_logi_d, data=casebase)
  pointest[6] <- coef(cbout_Dlogi_slm)[4]
  varest[6] <- diag(infjack.glm(cbout_Dlogi_slm, casebase$idx))[4]
  
  # V*D combined
  cbout_VlogiDlogi_slm <- glm(Yevent ~ timebasis_cb + notA + timebasisa_cb + offset(log(futime))
                              , family=binomial(link=logit), weights=sw_logilogi_a, data=casebase)
  pointest[7] <- coef(cbout_VlogiDlogi_slm)[4]
  varest[7] <- diag(infjack.glm(cbout_VlogiDlogi_slm, casebase$idx))[4]
  
  # unweighted
  cbout_unweight_slm <- glm(Yevent ~ timebasis_cb + notA+ timebasisa_cb + offset(log(futime))
                            , family=binomial(link=logit), data=casebase)
  pointest[8] <- coef(cbout_unweight_slm)[4]
  varest[8] <- diag(vcov(cbout_unweight_slm))[4]
  
  
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
  
  marginal_result <- list(pointest, varest)
  return(marginal_result)
}

#system.time(marginal_result <- marginalnew_one_round(seed, m, nobs, states, trates, tlim))
# 
# marginal_result[1]


