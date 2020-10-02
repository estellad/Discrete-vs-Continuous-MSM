# This is over 50 rounds of simulations.
# each round is documented in Marginal.R
#rm(list=ls())
library(parallel)

expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}

# #### Helper Functions ####
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/Marginal_sim3_slm.R") 
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/helper_sim3_slm.R")
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/estfun-glm.q")
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3/infjack-glm.q")



#### Default Specification ####
outpath <- 'C:/Users/dongy30/Desktop/Results/Scenario3/SLM'
#outpath <- 'C:/Users/dongy30/Desktop'
ncores <- max(1, detectCores())
simn <- 1000  # num simulations, acutally less
numsim <- floor(simn/ncores) # 25 sim per core
mtot <- ncores * numsim # actual total num sim 100
#numoutcome <- 2  # case-base and cox
numoutcome <- 2
numexposure <- 4
numdatasplit <- 1 # slm only
m <- numoutcome * numexposure * numdatasplit    # how many columns

# My effects of coef, lambda, transition matrix
input_coef = 1

setwd("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim3")
###### Per simulation specification of st, ss ######

sim <- function(no) { # sim for each core
  
  #### Helper Functions ####
  source("Marginal_sim3_slm.R") 
  source("helper_sim3_slm.R")
  source("estfun-glm.q")
  source("infjack-glm.q")
  
  from <- seq(1, mtot, by=numsim)[no] # seq = c(1, 26, 51, 76) # here "no" decides which core      #seq(1+100,mtot+100,by=numsim)
  
  pointest <- matrix(NA, numsim, m) 
  varest <- matrix(NA, numsim,  m)
  coverage <- matrix(NA, numsim,  m)
  power <- matrix(NA, numsim,  m) 
  
  
  for (j in 1:numsim) {
    iter <- from + (j - 1)
    
    set.seed(iter)
    
    marginal_result <- marginalnew_one_round(iter, m, nobs, states, trates, tlim)
    
    pointest[j, 1:m] <- marginal_result[[1]]
    varest[j, 1:m] <- marginal_result[[2]]
    
    cil <- marginal_result[[1]] + qnorm(0.025) * sqrt(marginal_result[[2]])
    ciu <- marginal_result[[1]] + qnorm(0.975) * sqrt(marginal_result[[2]])
    coverage[j,1:m] <- (input_coef >= cil) & (input_coef <= ciu)
    power[j,1:m] <- (cil > 0.0) | (ciu < 0.0)
  }
  
  write.table(pointest, file=file.path(outpath, paste0('pointest', no)))
  write.table(varest, file=file.path(outpath, paste0('varest', no)))
  write.table(coverage, file=file.path(outpath, paste0('coverage', no)))
  write.table(power, file=file.path(outpath, paste0('power', no)))
  return(NULL)
}

# # This will work on Mac/Linux:
# 
# system.time(
#   mclapply(1:ncores, sim, mc.cores=ncores, mc.silent=FALSE)
# )

# This should work on Windows:

cl <- makeCluster(ncores)
clusterExport(cl, c('expit','logit','outpath','m','input_coef','numsim','mtot'))
parallel::clusterEvalQ(cl, c(library(survival), library(splines)))
system.time(
  parLapply(cl, 1:ncores, sim)
)
stopCluster(cl)


##### Organize Results SLM: Run all at once read from desktop ######
## Combine old and new _all 1&2 results as _all1
pointestall <- data.frame(matrix(vector(), 0, m))
varestall <- data.frame(matrix(vector(), 0, m))
coverageall <- data.frame(matrix(vector(), 0, m))
powerall <- data.frame(matrix(vector(), 0, m))
for (i in 1:ncores) {
  # Result Folder
  pointest <- read.table(file=file.path(outpath, paste0('pointest', i)))
  varest <- read.table(file=file.path(outpath, paste0('varest', i)))
  coverage <- read.table(file=file.path(outpath, paste0('coverage', i)))
  power <- read.table(file=file.path(outpath, paste0('power', i)))
  
  pointestall <- rbind(pointestall, pointest)
  varestall <- rbind(varestall, varest)
  coverageall <- rbind(coverageall, coverage)
  powerall <- rbind(powerall, power)
}

pointestall <- na.omit(pointestall); nrow(pointestall)
varestall <- na.omit(varestall)
coverageall <- na.omit(coverageall)
powerall <- na.omit(powerall)

#### Only CB #######
# pointestall <- pointestall[, 5:8]
# varestall <- varestall[, 5:8]
# coverageall <- coverageall[, 5:8]
# powerall <- powerall[, 5:8]
# m<-4

########### End of SLM result organizaiton #############

results <- cbind(colMeans(pointestall),
                 colMeans(pointestall) - rep(input_coef,m), # Bias # Repeat true value 1, 3 times cause currently ctm 8 cols
                 apply(pointestall, 2, sd),
                 sqrt(colMeans(varestall)),
                 100*(apply(pointestall, 2, var) + (colMeans(pointestall) - rep(input_coef,m))^2), # 100*MSE
                 100*sqrt(colMeans(varestall)/nrow(pointestall)),    # 100*MCE
                 colMeans(coverageall),
                 colMeans(powerall))

vec <- c('_Vlogi_slm', '_Dlogi_slm', '_VlogiDlogi_slm', '_unweighted_slm')
#rownames(results) <- c(paste0("Case_base_Outcome", vec))
rownames(results) <- c(paste0("Cox_Outcome", vec), paste0("Case_base_Outcome", vec))
colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', '100xMCE', 'Coverage', 'Power')
round(results, 3)
write.table(round(results, 3), file=file.path(outpath, 'SLM_1000_results_sim3'))

