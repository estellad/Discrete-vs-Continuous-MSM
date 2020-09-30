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
simn <- 1000 # num simulations, acutally less
numsim <- floor(simn/ncores) # 6 sim per core
mtot <- ncores * numsim # actual total num sim 48
m <- 6      # how many columns

setwd("C:/Users/dongy30/Desktop/Marginal_Sim-master/Z_10_X_01/Sim4")
###### Per simulation specification of st, ss ######

sim <- function(no) { 
  source("helper_sim4_new.R")
  
  from <- seq(1,mtot,by=numsim)[no] 
  
  paras <- matrix(NA, numsim, m)
  
  for(k in 1:numsim){
    iter <- from + (k - 1)
    set.seed(iter)
    marginal_result <- checkAYtrue(iter, m, nobs, states, trates, tlim)
    paras[k, 1:m] <- marginal_result[[1]]
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
parasall <- as.data.frame(matrix(as.numeric(unlist(parasall)), simn, m))

########### End of SLM result organizaiton #############
results <- rbind(apply(parasall,2,min),
                 colMeans(parasall),
                 apply(parasall,2,max))

###### tab AY and D ######
colnames(results) <- c("YtrueA=0", "YtrueA=1", "YtrueA=2", 
                       "DtrtD=0",  "DtrtD=1",  "DtrtD=2")
rownames(results) <- c("estimates_min", "estimates_mean", "estimates_max")
round(results,3) 
# 1. c(-2, -1.5)                    | 1000, min at 61 134 323 409 935.     # 1. 7. 
#                YtrueA=0 YtrueA=1 YtrueA=2 DtrtD=0 DtrtD=1 DtrtD=2  
# estimates_min    15.000   13.000    5.000 217.000  52.000 139.000
# estimates_mean   28.925   28.224   13.355 256.796  79.705 183.278
# estimates_max    47.000   48.000   26.000 296.000 111.000 225.000
## 480 A=0 dose only worse, A=2 good.
#
# 2. c(-1.5, -2)                    | 1000, min at 429, YtrueA=2 is 2.     # 2. 4.
#                YtrueA=0 YtrueA=1 YtrueA=2 DtrtD=0 DtrtD=1 DtrtD=2
# estimates_min    15.000   15.000    2.000 216.000  95.000  86.000
# estimates_mean   30.206   31.667   11.621 265.203 135.492 119.142
# estimates_max    48.000   52.000   24.000 308.000 183.000 154.000
## 480 A=0 dose largest bias, A=2 weird.
#
# 3. c(-1.5, -2) + phi1s_x2->0 on D=2| , m1000in at 429, YtrueA=2 is 2.   # 3. 5.
#                YtrueA=0 YtrueA=1 YtrueA=2 DtrtD=0 DtrtD=1 DtrtD=2
# estimates_min    19.000   17.000    2.000 253.000 119.000  37.000
# estimates_mean   36.396   33.445   10.065 302.488 162.381  54.968
# estimates_max    55.000   51.000   22.000 352.000 203.000  85.000
## TODOing
# 4. c(-1.5, -2) + A=2 interaction (should be the same as 2.)
## 
#
# 5. c(-1.5, -2) + phi1s_x2->0 on D=2 + A=2 interaction (should be the same as 3.)
## 
# 
# 6. c(-2, -1.5) + phi1s_x2->0 on D=2| 1000, min at 429, YtrueA=1 is 1.  # 6. 8.
# #                YtrueA=0 YtrueA=1 YtrueA=2 DtrtD=0 DtrtD=1 DtrtD=2
# estimates_min    23.000    16.00    1.000 273.000  80.000   71.00
# estimates_mean   39.462    29.96   10.885 318.363 108.056   94.36
# estimates_max    58.000    54.00   22.000 357.000 145.000  122.00
## TODO
# 7. c(-2, -1.5) + A=2 interaction (should be the same as 1.)
## 
#
# 8. c(-2, -1.5) + phi1s_x2->0 on D=2 + A=2 interaction (should be the same as 6.)
##   

write.table(round(results, 3), file=file.path(outpath, 'CTM_200_sim4_checkAY'))

