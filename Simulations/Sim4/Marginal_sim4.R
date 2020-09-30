# rm(list=ls())
# source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/helper_sim4_new.R")
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
# m <- 12

marginalnew_one_round <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1 = trates1, tlim, seed)
  times <- times_wide_generation(long, nobs)
  dtimes_V <- sort(unique(times$vt[is.finite(times$vt)]))  # Stop of treatment time
  
  ## Exposure Data: V as outcome, X as cov (delay one interval)
  long_V_ctm <- get_long_V_ctm(times, dtimes_V, long)
  
  #################### Poisson Exposure Model Fitting ################
  # Marginal
  model_poi_V_ctm_marginal <- glm(Visit ~ Dprev + offset(log(futime)), family=poisson(link=log), data=long_V_ctm)
  # coef(model_poi_V_ctm_marginal)

  # Conditional
  model_poi_V_ctm_conditional <- glm(Visit ~ Xconf + Dprev + offset(log(futime)), family=poisson(link=log), data=long_V_ctm)
  # coef(model_poi_V_ctm_conditional)

  # ###################### Cox Exposure Model Fitting ##################
  # # Marginal
  model_cox_V_ctm_marginal <- coxph(Surv(start, stop, Visit) ~ Dprev, data=long_V_ctm, control = coxph.control(timefix=FALSE))
  # coef(model_cox_V_ctm_marginal)
  #
  # # Conditional
  model_cox_V_ctm_conditional <- coxph(Surv(start, stop, Visit) ~ Xconf + Dprev, data=long_V_ctm, control = coxph.control(timefix=FALSE))
  # coef(model_cox_V_ctm_conditional)


  ################## Multinomial Logistic D Treatment Model Fitting: Changed Dose #############
  # Treatment Model after visit
  D_expo_data <- long_V_ctm[long_V_ctm$stop == long_V_ctm$dt, ]

  # Marginal
  model_logi_V_ctm_marginal_d <- multinom(Dtrt ~ Dprev, data=D_expo_data)
  # coef(model_logi_V_ctm_marginal_d)

  # Conditional
  model_logi_V_ctm_conditional_d <- multinom(Dtrt ~ Xconf+Dprev, data=D_expo_data)
  #coef(model_logi_V_ctm_conditional_d)


  #################### Poisson Exposure Model Weighting Functions ################
  hz <- function(x, xt=Inf) {exp(coef(model_poi_V_ctm_conditional)[1] + coef(model_poi_V_ctm_conditional)[2] * I(x >= xt)
                                 + coef(model_poi_V_ctm_conditional)[3])}
  shz <- function(x, xt=Inf) {exp(coef(model_poi_V_ctm_marginal)[1] + coef(model_poi_V_ctm_marginal)[2] + 0 * x)}

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
  Psi0s_cox <- basehaz(model_cox_V_ctm_marginal, centered = FALSE)
  Psi0s_cox <- Psi0s_cox[Psi0s_cox[, 2] %in% dtimes_V, ]
  Psi0s_cox <- rbind(c(0.0, 0.0, 0.0), Psi0s_cox)

  # conditional
  Psi0_cox <- basehaz(model_cox_V_ctm_conditional, centered = FALSE) # 1331 2
  Psi0_cox <- Psi0_cox[Psi0_cox[, 2] %in% dtimes_V, ]
  Psi0_cox <- rbind(c(0.0, 0.0, 0.0), Psi0_cox)

  ###################### Outcome Data: Y as outcome, X, V as cov (delay one interval) ##########################
  # dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))
  # 
  # long_Y_ctm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times
  # long_Y_ctm$stop <- dtimes_Y                                          # automatic replicate dstops when finished.
  # long_Y_ctm <- long_Y_ctm[order(long_Y_ctm$idx, long_Y_ctm$stop),]            # actually useless, already ordered
  # long_Y_ctm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  # long_Y_ctm$futime <- long_Y_ctm$stop - long_Y_ctm$start
  # long_Y_ctm$Xconf <- (long_Y_ctm$xt <= long_Y_ctm$stop)
  # long_Y_ctm$Visit <- (long_Y_ctm$vt <= long_Y_ctm$stop)   #keep V=0->1
  # long_Y_ctm$Yevent <- NA
  # long_Y_ctm$tlim <- pmin(long_Y_ctm$vt, long_Y_ctm$yt, tlim) # rep vt (of idx 1) 4.44 100 times
  # 
  # long_Y_ctm$Dprev <- NA
  # for (i in 1:nobs){
  #   long_Y_ctm$Dprev[long_Y_ctm$idx == i] <- long$Dfrom[long$idx == i][1]
  # }
  # 
  # m_prob_d_all <- predict(model_logi_V_ctm_marginal_d, long_Y_ctm, type = "prob");
  # c_prob_d_all <- predict(model_logi_V_ctm_conditional_d, long_Y_ctm, type = "prob");
  # 
  # 
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
  #   long_Y_ctm$zhz[idx] <- ifelse(times$vt[i] <= t, hz(t, xt=times$xt[i]), 1.0)
  #   long_Y_ctm$szhz[idx] <- ifelse(times$vt[i] <= t, shz(t, xt=times$xt[i]), 1.0)
  #   long_Y_ctm$num[idx] <- long_Y_ctm$szsurv[idx] * long_Y_ctm$szhz[idx]
  #   long_Y_ctm$denom[idx] <- long_Y_ctm$zsurv[idx] * long_Y_ctm$zhz[idx]
  #   long_Y_ctm$iptw[idx] <- long_Y_ctm$num[idx]/long_Y_ctm$denom[idx]
  # 
  # 
  #   # Cox Exposure Model
  #   # Marginal weights
  #   # TODO: change to below style
  #   j <- findInterval(times$vt[i], Psi0s_cox$time)
  # 
  #   Psis_cox <- Psi0s_cox
  #   Psis_cox$hazard <- Psis_cox$hazard - c(0, Psi0s_cox$hazard[1:(nrow(Psi0s_cox)-1)])
  #   Psis_cox$hazard <- ifelse(Psis_cox$time > times$xt[i], Psis_cox$hazard * exp(coef(model_cox_V_ctm_marginal)[1]), Psis_cox$hazard)
  #   Psis_cox$hazard <- cumsum(Psis_cox$hazard)
  # 
  #   long_Y_ctm$Psis_cox_cum[idx] <- Psis_cox$hazard[findInterval(t, Psis_cox$time)]
  #   long_Y_ctm$Psis_cox_noncum[idx] <- ifelse(times$vt[i] <= t, Psis_cox$hazard[j]-Psis_cox$hazard[max(1, j-1)], 1)
  #   long_Y_ctm$spc_cox[idx] <- exp(-long_Y_ctm$Psis_cox_cum[idx])
  # 
  # 
  #   # Conditional weights; as formula
  #   Psi_cox <- Psi0_cox
  #   Psi_cox$hazard <- Psi_cox$hazard - c(0, Psi0_cox$hazard[1:(nrow(Psi0_cox)-1)])
  #   Psi_cox$hazard <- ifelse(Psi_cox$time > times$xt[i], Psi_cox$hazard * exp(coef(model_cox_V_ctm_conditional)[1]
  #                                                                             + coef(model_cox_V_ctm_conditional)[2]), Psi_cox$hazard)
  #   Psi_cox$hazard <- cumsum(Psi_cox$hazard)
  # 
  #   long_Y_ctm$Psi_cox_cum[idx] <- Psi_cox$hazard[findInterval(t, Psi_cox$time)]
  #   long_Y_ctm$Psi_cox_noncum[idx] <- ifelse(times$vt[i] <= t, Psi_cox$hazard[j]-Psi_cox$hazard[max(1, j-1)], 1)
  #   long_Y_ctm$pc_cox[idx] <- exp(-long_Y_ctm$Psi_cox_cum[idx])
  # 
  #   long_Y_ctm$sw_cox[idx] <- (long_Y_ctm$spc_cox[idx]*long_Y_ctm$Psis_cox_noncum[idx])/
  #     (long_Y_ctm$pc_cox[idx]*long_Y_ctm$Psi_cox_noncum[idx])
  # 
  # 
  #   ##############################
  #   # Treatment D Model
  #   long_Y_ctm$Dtrt[idx] <- ifelse(times$dt[i] <= t,
  #                                  long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]],
  #                                  long$Dfrom[long$idx == i][1])
  # 
  #   long_Y_ctm$m_prob_d[idx] <- ifelse(times$dt[i] <= t,
  #                                      m_prob_d_all[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
  #                                      m_prob_d_all[idx, (1+long$Dfrom[long$idx == i][1])])
  #   long_Y_ctm$c_prob_d[idx] <- ifelse(times$dt[i] <= t,
  #                                      c_prob_d_all[idx, (1+long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]])],
  #                                      c_prob_d_all[idx, (1+long$Dfrom[long$idx == i][1])])
  # 
  #   ## Logistic D Model
  #   j_d <- which(long_Y_ctm$stop[idx] > times$dt[i])[1]
  #   if(!is.na(j_d)){
  #     long_Y_ctm$sw_logi_d[idx] <- ifelse(times$dt[i] <= t, long_Y_ctm$m_prob_d[idx][j_d]/long_Y_ctm$c_prob_d[idx][j_d], # on and after visit
  #                                         #long_Y_ctm$m_prob_d[idx][1]/long_Y_ctm$c_prob_d[idx][1]
  #                                         1) # before visit
  #   }else{
  #     long_Y_ctm$sw_logi_d[idx] <- 1 #long_Y_ctm$m_prob_d[idx][1]/long_Y_ctm$c_prob_d[idx][1] # before visit
  #   }
  #   # if (i %% 100 == 0)
  #   #     print(i)
  #   #
  # }
  # long_Y_ctm <- subset(long_Y_ctm, stop <= yt)
  # long_Y_ctm$Dtrt <- as.factor(long_Y_ctm$Dtrt)
  # long_Y_ctm$Dtrt <- relevel(long_Y_ctm$Dtrt, ref = "1")
  # long_Y_ctm$A <- long_Y_ctm$Dtrt
  # 
  # # TODO:
  # # Combined weights sw_poi_a & sw_cox_a
  # long_Y_ctm$sw_poi_a <- long_Y_ctm$sw_logi_d * long_Y_ctm$iptw
  # long_Y_ctm$sw_cox_a <- long_Y_ctm$sw_logi_d * long_Y_ctm$sw_cox
  # 
  # # Time not A = 1, A=0
  # long_Y_ctm$timea0 <- as.numeric(long_Y_ctm$A == 0) * ifelse(is.finite(long_Y_ctm$vt), long_Y_ctm$stop - long_Y_ctm$vt, 0.0)
  # timebasisa0_cox <- bs(long_Y_ctm$timea0, Boundary.knots=c(0, tlim), degree=1)
  # 
  # # Time on A=2
  # # long_Y_ctm$timea2 <- as.numeric(long_Y_ctm$A == 2) * ifelse(is.finite(long_Y_ctm$vt), long_Y_ctm$stop - long_Y_ctm$vt, 0.0)
  # # long_Y_ctm$timea2 <- ifelse(long_Y_ctm$timea2 < 0, 0.0, long_Y_ctm$timea2)
  # # timebasisa2_cox <- bs(long_Y_ctm$timea2, Boundary.knots=c(0, tlim), degree=1)
  
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
  casebase$futime <- sum(ftime)/ms             # universal futime

  timebasis_cb <- bs(casebase$stop, Boundary.knots=c(0, tlim), degree=2)
  casebase$Xconf <- (casebase$xt <= casebase$stop)
  casebase$Visit <- (casebase$vt <= casebase$stop)  # V 0->1
  casebase$tlim <- pmin(casebase$vt, casebase$yt, tlim)


  # TODO: D model
  m_prob_d_all_cb <- predict(model_logi_V_ctm_marginal_d, newdata = casebase, type = "prob")  # marginal
  c_prob_d_all_cb <- predict(model_logi_V_ctm_conditional_d, newdata = casebase, type = "prob")  # conditional

  for (i in 1:nobs) {
    idx <- casebase$id == i
    if (sum(idx) > 0) {
      t <- pmin(casebase$stop[idx], casebase$tlim[idx])

      # Poisson V Exposure Model
      casebase$zsurv[idx] <- exp(-ch(t, xt=times$xt[i]))
      casebase$szsurv[idx] <- exp(-sch(t, xt=times$xt[i]))
      casebase$zhz[idx] <- ifelse(times$vt[i] <= t, hz(t, xt=times$xt[i]), 1.0)
      casebase$szhz[idx] <- ifelse(times$vt[i] <= t, shz(t, xt=times$xt[i]), 1.0)
      casebase$num[idx] <- casebase$szsurv[idx] * casebase$szhz[idx]
      casebase$denom[idx] <- casebase$zsurv[idx] * casebase$zhz[idx]
      casebase$iptw[idx] <- casebase$num[idx]/casebase$denom[idx]


      # Cox V Exposure Model
      # Marginal weights
      j <- findInterval(times$vt[i], Psi0s_cox$time)
      
      Psis_cox <- Psi0s_cox
      Psis_cox$hazard <- Psis_cox$hazard - c(0, Psi0s_cox$hazard[1:(nrow(Psi0s_cox)-1)])
      Psis_cox$hazard <- ifelse(Psis_cox$time > times$xt[i], Psis_cox$hazard * exp(coef(model_cox_V_ctm_marginal)[1]), Psis_cox$hazard)
      Psis_cox$hazard <- cumsum(Psis_cox$hazard)
      
      casebase$Psis_cox_cum[idx] <- Psis_cox$hazard[findInterval(t, Psis_cox$time)]
      casebase$Psis_cox_noncum[idx] <- ifelse(times$vt[i] <= t, Psis_cox$hazard[j]-Psis_cox$hazard[max(1, j-1)], 1)
      casebase$spc_cox[idx] <- exp(-casebase$Psis_cox_cum[idx])

      
      # Conditional weights
      Psi_cox <- Psi0_cox
      Psi_cox$hazard <- Psi_cox$hazard - c(0, Psi0_cox$hazard[1:(nrow(Psi0_cox)-1)])
      Psi_cox$hazard <- ifelse(Psi_cox$time > times$xt[i], Psi_cox$hazard * exp(coef(model_cox_V_ctm_conditional)[1]
                                                                                + coef(model_cox_V_ctm_conditional)[2]), Psi_cox$hazard)
      Psi_cox$hazard <- cumsum(Psi_cox$hazard)

      casebase$Psi_cox_cum[idx] <- Psi_cox$hazard[findInterval(t, Psi_cox$time)]
      casebase$Psi_cox_noncum[idx] <- ifelse(times$vt[i] <= t, Psi_cox$hazard[j]-Psi_cox$hazard[max(1, j-1)], 1)
      casebase$pc_cox[idx] <- exp(-casebase$Psi_cox_cum[idx])

      casebase$sw_cox[idx] <- (casebase$spc_cox[idx]*casebase$Psis_cox_noncum[idx])/
        (casebase$pc_cox[idx]*casebase$Psi_cox_noncum[idx])


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
  casebase$sw_poi_a <- casebase$sw_logi_d * casebase$iptw
  casebase$sw_cox_a <- casebase$sw_logi_d * casebase$sw_cox

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

  # ## Cox Outcome on A Model
  # # Only V Poisson model
  # coxout_Vpoi_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                          + cluster(idx), weight=iptw, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # Only V Cox model
  # coxout_Vcox_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                          + cluster(idx), weight=sw_cox, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # Only D Logistic model
  # coxout_Dlogi_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                           + cluster(idx), weight=sw_logi_d, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # V D Combined Poisson * Logistic Model
  # coxout_VpoiDlogi_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                               + cluster(idx), weight=sw_poi_a, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # V D Combined Cox * Logistic Model
  # coxout_VcoxDlogi_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                               + cluster(idx), weight=sw_cox_a, data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # Unweighted
  # coxout_unweight_ctm <- coxph(Surv(start, stop, Yevent) ~ A + timebasisa0_cox #+ timebasisa2_cox
  #                              + cluster(idx), data=long_Y_ctm, control = coxph.control(timefix=FALSE))
  # 
  # # For parameter A=0
  # pointest[1] <- coef(coxout_Vpoi_ctm)[1]
  # varest[1] <- vcov(coxout_Vpoi_ctm)[1,1]
  # pointest[2] <- coef(coxout_Vcox_ctm)[1]
  # varest[2] <- vcov(coxout_Vcox_ctm)[1,1]
  # pointest[3] <- coef(coxout_Dlogi_ctm)[1]
  # varest[3] <- vcov(coxout_Dlogi_ctm)[1,1]
  # pointest[4] <- coef(coxout_VpoiDlogi_ctm)[1]
  # varest[4] <- vcov(coxout_VpoiDlogi_ctm)[1,1]
  # pointest[5] <- coef(coxout_VcoxDlogi_ctm)[1]
  # varest[5] <- vcov(coxout_VcoxDlogi_ctm)[1,1]
  # pointest[6] <- coef(coxout_unweight_ctm)[1]
  # varest[6] <- vcov(coxout_unweight_ctm)[1,1]
  # 
  # # For parameter A=2
  # pointest2[1] <- coef(coxout_Vpoi_ctm)[2]
  # varest2[1] <- vcov(coxout_Vpoi_ctm)[2,2]
  # pointest2[2] <- coef(coxout_Vcox_ctm)[2]
  # varest2[2] <- vcov(coxout_Vcox_ctm)[2,2]
  # pointest2[3] <- coef(coxout_Dlogi_ctm)[2]
  # varest2[3] <- vcov(coxout_Dlogi_ctm)[2,2]
  # pointest2[4] <- coef(coxout_VpoiDlogi_ctm)[2]
  # varest2[4] <- vcov(coxout_VpoiDlogi_ctm)[2,2]
  # pointest2[5] <- coef(coxout_VcoxDlogi_ctm)[2]
  # varest2[5] <- vcov(coxout_VcoxDlogi_ctm)[2,2]
  # pointest2[6] <- coef(coxout_unweight_ctm)[2]
  # varest2[6] <- vcov(coxout_unweight_ctm)[2,2]


  # Case-base Outcome Models
  # Only V poisson
  cbout_Vpoi_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                        + offset(log(futime)), family=binomial(link=logit), weights=iptw, data=casebase)

  # Only V cox
  cbout_Vcox_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                        + offset(log(futime)), family=binomial(link=logit), weights=sw_cox, data=casebase)

  # Only D logistic
  cbout_Dlogi_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                         + offset(log(futime)), family=binomial(link=logit), weights=sw_logi_d, data=casebase)

  # V * D combined: poisson * logistic
  cbout_VpoiDlogi_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                             + offset(log(futime)), family=binomial(link=logit), weights=sw_poi_a, data=casebase)

  # V * D combined: cox * logistic
  cbout_VcoxDlogi_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                             + offset(log(futime)), family=binomial(link=logit), weights=sw_cox_a, data=casebase)

  # Unweighted
  cbout_unweight_ctm <- glm(Yevent ~ timebasis_cb + A + timebasisa0_cb #+ timebasisa2_cb 
                            + offset(log(futime)), family=binomial(link=logit), data=casebase)

  # For parameter A=0
  pointest[7] <- coef(cbout_Vpoi_ctm)[4]
  varest[7] <- diag(infjack.glm(cbout_Vpoi_ctm, casebase$idx))[4]
  pointest[8] <- coef(cbout_Vcox_ctm)[4]
  varest[8] <- diag(infjack.glm(cbout_Vcox_ctm, casebase$idx))[4]
  pointest[9] <- coef(cbout_Dlogi_ctm)[4]
  varest[9] <- diag(infjack.glm(cbout_Dlogi_ctm, casebase$idx))[4]
  pointest[10] <- coef(cbout_VpoiDlogi_ctm)[4]
  varest[10] <- diag(infjack.glm(cbout_VpoiDlogi_ctm, casebase$idx))[4]
  pointest[11] <- coef(cbout_VcoxDlogi_ctm)[4]
  varest[11] <- diag(infjack.glm(cbout_VcoxDlogi_ctm, casebase$idx))[4]
  pointest[12] <- coef(cbout_unweight_ctm)[4]
  varest[12] <- diag(vcov(cbout_unweight_ctm))[4]

  # For parameter A=2
  pointest2[7] <- coef(cbout_Vpoi_ctm)[5]
  varest2[7] <- diag(infjack.glm(cbout_Vpoi_ctm, casebase$idx))[5]
  pointest2[8] <- coef(cbout_Vcox_ctm)[5]
  varest2[8] <- diag(infjack.glm(cbout_Vcox_ctm, casebase$idx))[5]
  pointest2[9] <- coef(cbout_Dlogi_ctm)[5]
  varest2[9] <- diag(infjack.glm(cbout_Dlogi_ctm, casebase$idx))[5]
  pointest2[10] <- coef(cbout_VpoiDlogi_ctm)[5]
  varest2[10] <- diag(infjack.glm(cbout_VpoiDlogi_ctm, casebase$idx))[5]
  pointest2[11] <- coef(cbout_VcoxDlogi_ctm)[5]
  varest2[11] <- diag(infjack.glm(cbout_VcoxDlogi_ctm, casebase$idx))[5]
  pointest2[12] <- coef(cbout_unweight_ctm)[5]
  varest2[12] <- diag(vcov(cbout_unweight_ctm))[5]
  
  # 
  # 
  # # # # A look at weights
  # # TODO: in simulation scenario 2 of not separating V and D for exposure, max weights is 5 for V exposure 
  # # V exposure models weights
  # summary(long_Y_ctm$iptw)
  # summary(long_Y_ctm$sw_cox)
  # # summary(casebase$iptw)
  # # summary(casebase$sw_cox)
  # 
  # # par(mfrow = c(2,2))
  # # plot(long_Y_ctm$stop, long_Y_ctm$iptw, main = "Poisson Exposure with Cox Outcome")
  # # plot(long_Y_ctm$stop, long_Y_ctm$sw_cox, main = "Cox Exposure with Cox Outcome")
  # # plot(casebase$stop, casebase$iptw, main = "Poisson Exposure with Case-base Outcome")
  # # plot(casebase$stop, casebase$sw_cox, main = "Cox Exposure with Case-base Outcome")
  # 
  # # D models weights
  # summary(long_Y_ctm$sw_logi_d)
  # #summary(casebase$sw_logi_d)
  # 
  # # # V D combined weights
  #  summary(long_Y_ctm$sw_poi_a)
  #  summary(long_Y_ctm$sw_cox_a)
  # # summary(casebase$sw_poi_a)
  # # summary(casebase$sw_cox_a)
  # 
  # # par(mfrow = c(2,2))
  # # plot(long_Y_ctm$stop, long_Y_ctm$sw_poi_a, main = "Poisson Exposure with Cox Outcome")
  # # plot(long_Y_ctm$stop, long_Y_ctm$sw_cox_a, main = "Cox Exposure with Cox Outcome")
  # # plot(casebase$stop, casebase$sw_poi_a, main = "Poisson Exposure with Case-base Outcome")
  # # plot(casebase$stop, casebase$sw_cox_a, main = "Cox Exposure with Case-base Outcome")
  # 
  marginal_result <- list(pointest, varest, pointest2, varest2)
  return(marginal_result)
}
#system.time(marginal_result <- marginalnew_one_round(seed, m, nobs, states, trates, tlim))
# # 
# marginal_result[1]
# marginal_result[3]


