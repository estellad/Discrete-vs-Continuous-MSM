# This is over 50 rounds of simulations.
# each round is documented in Marginal.R
#rm(list=ls())
library(parallel)

expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}

# #### Helper Functions ####
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/Marginal_sim4.R") 
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/helper_sim4_new.R")
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/estfun-glm.q")
source("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4/infjack-glm.q")


#### Default Specification ####
outpath <- 'C:/Users/dongy30/Desktop/Results/Scenario4/CTM'
#outpath <- 'C:/Users/dongy30/Desktop/Results/Scenario4/SLM'
ncores <- max(1, detectCores())
simn <- 1000# num simulations, acutally less
numsim <- floor(simn/ncores) # 25 sim per core
mtot <- ncores * numsim # actual total num sim 100
#numoutcome <- 2  # case-base and cox
numoutcome <- 2
numexposure <- 6
#numexposure <- 3 # poi_weighted, unweighted
numdatasplit <- 1 # slm only
m <- numoutcome * numexposure * numdatasplit    # how many columns

# My effects of coef, lambda, transition matrix
input_coef = 1 # A=0
input_coef2 = -1 # A=2

setwd("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4")
###### Per simulation specification of st, ss ######

sim <- function(no) { # sim for each core
  
  #### Helper Functions ####
  source("Marginal_sim4.R") 
  source("helper_sim4_new.R")
  source("estfun-glm.q")
  source("infjack-glm.q")
  
  from <- seq(1,mtot,by=numsim)[no] # seq = c(1, 26, 51, 76) # here "no" decides which core      #seq(1+100,mtot+100,by=numsim)
  
  pointest <- matrix(NA, numsim, m) 
  varest <- matrix(NA, numsim,  m)
  coverage <- matrix(NA, numsim,  m)
  power <- matrix(NA, numsim,  m) 
  
  pointest2 <- matrix(NA, numsim, m) 
  varest2 <- matrix(NA, numsim,  m)
  coverage2 <- matrix(NA, numsim,  m)
  power2 <- matrix(NA, numsim,  m) 
  
  
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
    
    pointest2[j, 1:m] <- marginal_result[[3]]
    varest2[j, 1:m] <- marginal_result[[4]]
    
    cil <- marginal_result[[3]] + qnorm(0.025) * sqrt(marginal_result[[4]])
    ciu <- marginal_result[[3]] + qnorm(0.975) * sqrt(marginal_result[[4]])
    coverage2[j,1:m] <- (input_coef2 >= cil) & (input_coef2 <= ciu)
    power2[j,1:m] <- (cil > 0.0) | (ciu < 0.0)
  }
  
  write.table(pointest, file=file.path(outpath, paste0('pointest', no)))
  write.table(varest, file=file.path(outpath, paste0('varest', no)))
  write.table(coverage, file=file.path(outpath, paste0('coverage', no)))
  write.table(power, file=file.path(outpath, paste0('power', no)))
  
  write.table(pointest2, file=file.path(outpath, paste0('pointest2', no)))
  write.table(varest2, file=file.path(outpath, paste0('varest2', no)))
  write.table(coverage2, file=file.path(outpath, paste0('coverage2', no)))
  write.table(power2, file=file.path(outpath, paste0('power2', no)))
  return(NULL)
}

# # This will work on Mac/Linux:
# 
# system.time(
#   mclapply(1:ncores, sim, mc.cores=ncores, mc.silent=FALSE)
# )

# This should work on Windows:

cl <- makeCluster(ncores)
clusterExport(cl, c('expit','logit','outpath','m','input_coef','input_coef2','numsim','mtot'))
parallel::clusterEvalQ(cl, c(library(survival), library(splines), library(nnet)))
#parLapply(cl, 1:ncores, sim)
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

pointestall2 <- data.frame(matrix(vector(), 0, m))
varestall2 <- data.frame(matrix(vector(), 0, m))
coverageall2 <- data.frame(matrix(vector(), 0, m))
powerall2 <- data.frame(matrix(vector(), 0, m))
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
  
  
  pointest2 <- read.table(file=file.path(outpath, paste0('pointest2', i)))
  varest2 <- read.table(file=file.path(outpath, paste0('varest2', i)))
  coverage2 <- read.table(file=file.path(outpath, paste0('coverage2', i)))
  power2 <- read.table(file=file.path(outpath, paste0('power2', i)))
  
  pointestall2 <- rbind(pointestall2, pointest2)
  varestall2 <- rbind(varestall2, varest2)
  coverageall2 <- rbind(coverageall2, coverage2)
  powerall2 <- rbind(powerall2, power2)
}

# pointestall <- na.omit(pointestall); nrow(pointestall)
# varestall <- na.omit(varestall)
# coverageall <- na.omit(coverageall)
# powerall <- na.omit(powerall)
# 
# pointestall2 <- na.omit(pointestall2); nrow(pointestall2)
# varestall2 <- na.omit(varestall2)
# coverageall2 <- na.omit(coverageall2)
# powerall2 <- na.omit(powerall2)

#### Only CB #######
pointestall <- pointestall[, 7:12]
varestall <- varestall[, 7:12]
coverageall <- coverageall[, 7:12]
powerall <- powerall[, 7:12]

pointestall2 <- pointestall2[, 7:12]
varestall2 <- varestall2[, 7:12]
coverageall2 <- coverageall2[, 7:12]
powerall2 <- powerall2[, 7:12]
m<-6

########### End of SLM result organizaiton #############
# A=0
results <- cbind(colMeans(pointestall),
                 colMeans(pointestall) - rep(input_coef,m), # Bias # Repeat true value 1, 3 times cause currently ctm 12 cols
                 apply(pointestall, 2, sd),
                 sqrt(colMeans(varestall)),
                 sqrt(apply(pointestall, 2, var) + (colMeans(pointestall) - rep(input_coef,m))^2), # RMSE
                 100*sqrt(colMeans(varestall)/nrow(pointestall)),    # 100*MCE
                 colMeans(coverageall),
                 colMeans(powerall))

#vec <- c('_poi_ctm', '_cox_ctm', '_unweighted_ctm')
vec <- c('_Vpoi_ctm', '_Vcox_ctm', '_Dlogi_ctm', '_VDpoilogi_ctm', '_VDcoxlogi_ctm', '_unweighted_ctm')
#vec <- c('_poi_ctm', '_unweighted_ctm')
#vec <- c('_logi_slm', '_cox_slm', '_unweighted_slm')
#rownames(results) <- c(paste0("Cox_Outcome", vec))
rownames(results) <- c(paste0("Case_base_Outcome", vec))
#rownames(results) <- c(paste0("Cox_Outcome", vec), paste0("Case_base_Outcome", vec))
#colnames(results) <- c('Mean', 'SD', 'Mean SE', '100xMCE', 'Power')
colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', 'RMSE', '100xMCE', 'Coverage', 'Power')
round(results, 3)
write.table(round(results, 3), file=file.path(outpath, 'CTM_200_sim4_A0'))


# A=2
results2 <- cbind(colMeans(pointestall2),
                 colMeans(pointestall2) - rep(input_coef2,m), # Bias # Repeat true value 1, 3 times cause currently ctm 12 cols
                 apply(pointestall2, 2, sd),
                 sqrt(colMeans(varestall2)),
                 sqrt(apply(pointestall2, 2, var) + (colMeans(pointestall2) - rep(input_coef2,m))^2), # RMSE
                 100*sqrt(colMeans(varestall2)/nrow(pointestall2)),    # 100*MCE
                 colMeans(coverageall2),
                 colMeans(powerall2))

#vec <- c('_poi_ctm', '_cox_ctm', '_unweighted_ctm')
vec <- c('_Vpoi_ctm', '_Vcox_ctm', '_Dlogi_ctm', '_VDpoilogi_ctm', '_VDcoxlogi_ctm', '_unweighted_ctm')
#vec <- c('_poi_ctm', '_unweighted_ctm')
#vec <- c('_logi_slm', '_cox_slm', '_unweighted_slm')
#rownames(results2) <- c(paste0("Cox_Outcome", vec))
rownames(results2) <- c(paste0("Case_base_Outcome", vec))
#rownames(results2) <- c(paste0("Cox_Outcome", vec), paste0("Case_base_Outcome", vec))
#colnames(results2) <- c('Mean', 'SD', 'Mean SE', '100xMCE', 'Power')
colnames(results2) <- c('Mean', 'Bias', 'SD', 'Mean SE', 'RMSE', '100xMCE', 'Coverage', 'Power')
round(results2, 3)
write.table(round(results2, 3), file=file.path(outpath, 'CTM_200_sim4_A2'))
