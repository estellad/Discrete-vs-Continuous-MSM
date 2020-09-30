rm(list=ls())
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/helper_sim4_new.R")
library(parallel)
library(survival)
library(splines)
library(dplyr)
library(nnet)
expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}

nobs <- 1000
outpath <- 'C:/Users/dongy30/Desktop'

ncores <- max(1, detectCores())
simn <- 1000# num simulations, acutally less
numsim <- floor(simn/ncores) # 6 sim per core
mtot <- ncores * numsim # actual total num sim 48
m <- 3
#m <- 7
#m <- 4
#m <- 7    # how many columns

setwd("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4")
###### Per simulation specification of st, ss ######

sim <- function(no) { 
  source("helper_sim4_new.R")
  
  from <- seq(1,mtot,by=numsim)[no] # seq = c(1, 26, 51, 76) # here "no" decides which core      #seq(1+100,mtot+100,by=numsim)
  
  paras <- matrix(NA, numsim, m)
  
  for(k in 1:numsim){
    iter <- from + (k - 1)
    set.seed(iter)
    
    long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1 = trates1, tlim, seed=iter)
    times <- times_wide_generation(long, nobs)
    dtimes_V <- sort(unique(times$vt[is.finite(times$vt)]))

    long_V_ctm <- get_long_V_ctm(times, dtimes_V, long)
    
    
    D_expo_data <- long_V_ctm[long_V_ctm$stop == long_V_ctm$dt, ]
    # Multinom + Dprev
    model_D_D0 <- multinom(Dtrt ~ Xconf + Dprev, data=D_expo_data)
    # A=1 intercept
    paras[k, 1] <- coef(model_D_D0)[1,1]
    paras[k, 2] <- coef(model_D_D0)[1,2]
    # A=2 intercept
    paras[k, 3] <- coef(model_D_D0)[2,1]
    paras[k, 4] <- coef(model_D_D0)[2,2]

    # Effect of D0=2 on A=1
    paras[k, 5] <- coef(model_D_D0)[1,3]
    # Effect of D0=2 on A=2
    paras[k, 6] <- coef(model_D_D0)[2,3]


    D0_expo_data <- long_V_ctm %>%
      group_by(idx)%>%
      filter(row_number() == 1)
    D0_expo_data$Dprev <- as.factor(D0_expo_data$Dprev)
    D0_expo_data$Dprev <- relevel(D0_expo_data$Dprev, ref = "1")
    model_D_0 <- glm(Dprev ~ 1, family=binomial(link=logit), data=D0_expo_data)
    paras[k, 7] <- coef(model_D_0)
    
    
    #################################################################################
    # # Multinomial
    # model_D <- multinom(Dtrt ~ Xconf, data=D_expo_data)
    # # coef(model_D)
    # # A=1 intercept
    # paras[k, 1] <- coef(model_D)[1,1]
    # paras[k, 2] <- coef(model_D)[1,2]
    # # A=2 intercept
    # paras[k, 3] <- coef(model_D)[2,1]
    # paras[k, 4] <- coef(model_D)[2,2]
    
    #################################################################################
    #dtimes_X <- sort(unique(times$xt[is.finite(times$xt)]))
    #dtimes_Y <- sort(unique(times$yt[is.finite(times$yt)]))


    # Check conditional X model generated correctly
    ############### Exposure Data: X as outcome, D as cov #########################
    # long_X_ctm <- get_long_X_ctm(times, dtimes_X, long)
    # model_X <- glm(Xconf ~ Dtrt+ offset(log(futime)), family=poisson(link=log), data = long_X_ctm)
    # #coef(model_X)
    # paras[k, 1] <- coef(model_X)[1]
    # paras[k, 2] <- coef(model_X)[2]
    # paras[k, 3] <- coef(model_X)[3]
    # 
    # 
    # # Check conditional Y model generated correctly
    # ###################### Outcome Data: Y as outcome, X, V as cov (delay one interval) ##########################
    # long_Y_ctm <- get_long_Y_ctm(times, dtimes_Y, long)
    # model_Y <- glm(Yevent ~ A + Xconf+ offset(log(futime)), family=poisson(link=log), data = long_Y_ctm)
    # #coef(model_Y)
    # paras[k, 4] <- coef(model_Y)[1]
    # paras[k, 5] <- coef(model_Y)[2]
    # paras[k, 6] <- coef(model_Y)[3]
    # paras[k, 7] <- coef(model_Y)[4]

  }
  write.table(paras, file=file.path(outpath, paste0('paras', no)))
  
  return(NULL)
}# sim for each core


cl <- makeCluster(ncores)
clusterExport(cl, c('expit','logit','outpath','m','numsim','mtot'))
parallel::clusterEvalQ(cl, c(library(survival), library(splines), library(nnet), library(dplyr)))
system.time(
  parLapply(cl, 1:ncores, sim)
)
stopCluster(cl)


##### Organize Results SLM: Run all at once read from desktop ######
## Combine old and new _all 1&2 results as _all1
parasall <- data.frame(matrix(vector(), 0, m))
for (i in 1:ncores) {
  # Result Folder
  paras <- read.table(file=file.path(outpath, paste0('paras', i)))
  parasall <- rbind(parasall, paras)
}
parasall <- na.omit(parasall); nrow(parasall)

########### End of SLM result organizaiton #############
results <- colMeans(parasall)
###### D conditional model ######
est_result <- rbind(
  c(0.1, 1, 0.3, 1.5, 
    0.5, 0.5, -0.3),
  results
)

colnames(est_result) <- c("intercept_wD0_A=1", "coef_wD0_A=1", "intercept_wD0_A=2", "coef_wD0_A=2", 
                          "D0=2_on_A=1", "D0=2_on_A=2", "coef_D0")
rownames(est_result) <- c("trueval", "estimates")
round(est_result,3)
write.table(round(est_result, 3), file=file.path(outpath, 'CTM_200_sim4_checDgeneration'))


###### X Y conditional models #######
# For A=1 and A=2, where ref = 0
# est_result <- rbind(
#   c(alpha_x,     beta_x,      beta_x, alpha_y, beta_y1, beta_y2, gamma_y),
#   c(exp(results[1]), results[c(2,3)], exp(results[4]), results[c(5:7)])
# )
# 
# colnames(est_result) <- c("alpha_x",   "beta_x_A=1",  "beta_x_A=2",  "alpha_y", "beta_y1_A=1",  "beta_y2_A=2",  "gamma_y_Xconf")
# rownames(est_result) <- c("trueval", "estimates")
# round(est_result,3)
# write.table(round(est_result, 3), file=file.path(outpath, 'CTM_200_sim4_checkXYgeneration'))

