rm(list=ls())
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/helper_sim3.R")
library(survival)
library(splines)
logit <- function(p) {log(p)-log(1-p)}
# ### This part will be called in Marginal_sim-repeat.R; Only exist for testing this one round function ###
nobs <- 1000
seed = 80
m <- 12


long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1 = trates1, tlim, seed)
times <- times_wide_generation(long, nobs)
dtimes_V <- sort(unique(times$vt[is.finite(times$vt)]))  # Stop of treatment time
dtimes_X <- sort(unique(times$xt[is.finite(times$xt)]))
dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))


# Check conditional X model generated correctly
############### Exposure Data: X as outcome, D as cov #########################

get_long_X_ctm <- function(times, dtimes_X, long){
  long_X_ctm <- times[rep(1:nobs, times=rep(length(dtimes_X), nobs)),] 
  long_X_ctm$stop <- dtimes_X                                           
  long_X_ctm <- long_X_ctm[order(long_X_ctm$id, long_X_ctm$stop),]  
  long_X_ctm$start <- NA
  long_X_ctm$Xconf <- NA
  long_X_ctm$Dtrt <- NA
  
  for (i in 1:nobs) {
    idx <- long_X_ctm$id == i
    long_X_ctm$start[idx] <- c(0, dtimes_X[1:(length(dtimes_X)-1)])
    long_X_ctm$Xconf[idx] <- ifelse(times$xt[i] == long_X_ctm$stop[idx], 1, 0)   # X 1->0 
    long_X_ctm$Dtrt[idx] <- ifelse(times$dt[i] < long_X_ctm$stop[idx], 
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 
                                   long$Dfrom[long$idx == i][1])  
    # because start with V=0 & D = Dinit
    
    # if (i %% 100 == 0)
    #     print(i)
  }
  long_X_ctm <- subset(long_X_ctm, stop <= xt)
  long_X_ctm <- subset(long_X_ctm, start <= yt)
  long_X_ctm$stop <- ifelse(long_X_ctm$stop > long_X_ctm$yt, long_X_ctm$yt, long_X_ctm$stop)
  long_X_ctm$futime <- long_X_ctm$stop - long_X_ctm$start
  
  return(long_X_ctm)
}

long_X_ctm <- get_long_X_ctm(times, dtimes_X, long)
long_V_ctm <- get_long_V_ctm(times, dtimes_V, long)


# TODO: check X Y generation 
model_X <- glm(Xconf ~ Dtrt+ offset(log(futime)), family=poisson(link=log), data = long_X_ctm)
coef(model_X)
# (Intercept)        Dtrt 
# -2.895251       1.159817  
# = log(0.05528515)
# alpha_x = 0.04 beta_x = 1.5




# Check conditional Y model generated correctly
###################### Outcome Data: Y as outcome, X, V as cov (delay one interval) ##########################

get_long_Y_ctm <- function(times, dtimes_Y, long){
  long_Y_ctm <- times[rep(1:nobs, times=rep(length(dtimes_Y), nobs)),] # rep 1:1000, each 100 (dtimes_Y) times
  long_Y_ctm$stop <- dtimes_Y                                          # automatic replicate dstops when finished.
  long_Y_ctm <- long_Y_ctm[order(long_Y_ctm$idx, long_Y_ctm$stop),]            # actually useless, already ordered
  long_Y_ctm$start <- c(0, dtimes_Y[1:(length(dtimes_Y)-1)])
  long_Y_ctm$futime <- long_Y_ctm$stop - long_Y_ctm$start
  long_Y_ctm$Xconf <- (long_Y_ctm$xt <= long_Y_ctm$stop)
  long_Y_ctm$Visit <- (long_Y_ctm$vt <= long_Y_ctm$stop)   #keep V=0->1
  long_Y_ctm$Yevent <- NA
  long_Y_ctm$tlim <- pmin(long_Y_ctm$vt, long_Y_ctm$yt, tlim) # rep vt (of idx 1) 4.44 100 times
  
  
  for (i in 1:nobs) {
    idx <- long_Y_ctm$idx == i
    long_Y_ctm$Yevent[idx] <- ifelse(times$yt[i] == long_Y_ctm$stop[idx], TRUE, FALSE)
    t <- pmin(long_Y_ctm$stop[idx], long_Y_ctm$tlim[idx])
    
    ##############################
    # Treatment D Model
    long_Y_ctm$Dtrt[idx] <- ifelse(times$dt[i] <= t, 
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 
                                   long$Dfrom[long$idx == i][1])
    
  }
  
  long_Y_ctm <- subset(long_Y_ctm, stop <= yt)
  long_Y_ctm$Dtrt <- as.factor(long_Y_ctm$Dtrt)
  long_Y_ctm$Dtrt <- relevel(long_Y_ctm$Dtrt, ref = "0")
  long_Y_ctm$A <- long_Y_ctm$Dtrt
  
  return(long_Y_ctm)
}
long_Y_ctm <- get_long_Y_ctm(times, dtimes_Y, long)

model_Y <- glm(Yevent ~ A + Xconf+ offset(log(futime)), family=poisson(link=log), data = long_Y_ctm)
coef(model_Y)
# >   coef(model_Y)
# (Intercept)          A1   XconfTRUE 
# -3.8606708  -0.9881394   1.4509973 
# = log(0.02)
# alpha_y = 0.02  beta_y = -1  gamma_y = 1.5
