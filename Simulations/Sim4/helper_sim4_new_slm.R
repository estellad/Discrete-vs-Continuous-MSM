# This is everything needed to run the scripts: multistate_3_outcome_cov.R and multistate_subset_each_lambda.R
# and sim_repeat.R

#####  Default Parameters  #####
library(survival)
library(dplyr)
expit <- function(x) {1/(1+exp(-x))}
nobs <- 1000

tlim <- 5     # 1   

# # My new trying parameters (for more X before V so no positivity violation): 
alpha_v <- 0.08
alpha_x <- 0.04
alpha_y <- 0.02

beta_v <- 1.5
beta_x0 <- 0; beta_x <- 1.5; beta_x2 <- 1.5;
beta_y0 <- 0; beta_y1 <- -1; beta_y2 <- -2; 
gamma_y <- 1.5


# in D=0 matrix
# lambdaV                       
lambda210 <- alpha_v; lambda210
# lambdaX
lambda130 <- alpha_x*exp(beta_x0); lambda130
# lambdaY
lambda140 <- alpha_y*exp(beta_y0); lambda140
# lambdaX
lambda250 <- alpha_x*exp(beta_x0); lambda250
# lambdaY
lambda260 <- alpha_y*exp(beta_y0); lambda260
# lambdaV
lambda530 <- alpha_v*exp(beta_v); lambda530
# lambdaY
lambda370 <- alpha_y*exp(beta_y0 + gamma_y); lambda370
# lambdaY
lambda580 <- alpha_y*exp(beta_y0 + gamma_y); lambda580

# D=0 matrix
rates0 <- matrix(c(0,        0, lambda130, lambda140,         0,         0,         0,         0,
                   lambda210,0,         0,         0, lambda250, lambda260,         0,         0,
                   0,        0,         0,         0,         0,         0, lambda370,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0, lambda530,         0,         0,         0,         0, lambda580,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0), 
                 8, 8, byrow=TRUE)
rownames(rates0) <- colnames(rates0) <- c('1', '2', '3', '4', '5', '6', '7', '8')
trates0 <- rowSums(rates0)


# in D=1 matrix
# lambdaV                       
lambda211 <- alpha_v; lambda211
# lambdaX
lambda131 <- alpha_x*exp(beta_x); lambda131
# lambdaY
lambda141 <- alpha_y*exp(beta_y1); lambda141
# # lambdaX
lambda251 <- alpha_x*exp(beta_x); lambda251
# # lambdaY
lambda261 <- alpha_y*exp(beta_y1); lambda261 
# lambdaV
lambda531 <- alpha_v*exp(beta_v); lambda531
# lambdaY
lambda371 <- alpha_y*exp(beta_y1 + gamma_y); lambda371
# lambdaY
lambda581 <- alpha_y*exp(beta_y1 + gamma_y); lambda581

# D=1 matrix in 
rates1 <- matrix(c(0,        0, lambda131, lambda141,         0,         0,         0,         0,
                   lambda211,0,         0,         0, lambda251, lambda261,         0,         0,
                   0,        0,         0,         0,         0,         0, lambda371,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0, lambda531,         0,         0,         0,         0, lambda581,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0), 
                 8, 8, byrow=TRUE)
rownames(rates1) <- colnames(rates1) <- c('1', '2', '3', '4', '5', '6', '7', '8')
trates1 <- rowSums(rates1)


# in D=2 matrix
# lambdaV                       
lambda212 <- alpha_v; lambda212
# lambdaX
lambda132 <- alpha_x*exp(beta_x2); lambda132
# lambdaY
lambda142 <- alpha_y*exp(beta_y2); lambda142
# lambdaX
lambda252 <- alpha_x*exp(beta_x2); lambda252
# lambdaY
lambda262 <- alpha_y*exp(beta_y2); lambda262 
# lambdaV
lambda532 <- alpha_v*exp(beta_v); lambda532
# lambdaY
lambda372 <- alpha_y*exp(beta_y2 + gamma_y); lambda372
# lambdaY
lambda582 <- alpha_y*exp(beta_y2 + gamma_y); lambda582

# D=2 matrix in 
rates2 <- matrix(c(0,        0, lambda132, lambda142,         0,         0,         0,         0,
                   lambda212,0,         0,         0, lambda252, lambda262,         0,         0,
                   0,        0,         0,         0,         0,         0, lambda372,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0, lambda532,         0,         0,         0,         0, lambda582,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0,
                   0,        0,         0,         0,         0,         0,         0,         0), 
                 8, 8, byrow=TRUE)
rownames(rates2) <- colnames(rates2) <- c('1', '2', '3', '4', '5', '6', '7', '8')
trates2 <- rowSums(rates2)


# rates[which(rates == 0)] <- NA
# library(mstate)
# paths(rates, start = 2)
#       [,1] [,2] [,3] [,4]
# [1,]    2   NA   NA   NA
# [2,]    2    1   NA   NA
# [3,]    2    1    3   NA
# [4,]    2    1    3    7
# [5,]    2    1    4   NA
# [6,]    2    5   NA   NA
# [7,]    2    5    3   NA
# [8,]    2    5    3    7
# [9,]    2    5    8   NA
# [10,]    2    6   NA   NA

states <- as.numeric(colnames(rates0)) # c(1, 2, 3, 4, 5, 6, 7, 8) # Here rates1 and rates0 are the same
numstates <- length(states)
tstates <- c(4, 6, 7, 8) 
cols <- c('gray','green','red','black', 'orange', 'blue', 'coral','cyan')

# input_coef <- matrix(c(alpha_v, alpha_x, alpha_y, alpha_x+beta_x, alpha_y+beta_y, alpha_v+beta_v, alpha_y+gamma_y, alpha_y+gamma_y+beta_y), nrow=1); input_coef
# colnames(input_coef) <- c('alpha_v', 'alpha_x', 'alpha_y', 'alpha_x+beta_x', 'alpha_y+beta_y', 'alpha_v+beta_v', 'alpha_y+gamma_y', 
#                        'alpha_y+gamma_y+beta_y')
# rownames(input_coef) <- "input_coef"

mat_haz0 <- matrix(c(lambda210, lambda130, lambda140, lambda250, lambda260, lambda530, lambda370, lambda580), nrow = 1)
colnames(mat_haz0) <- c('lambda210', 'lambda130', 'lambda140', 'lambda250', 'lambda260', 'lambda530', 'lambda370', 'lambda580')
rownames(mat_haz0) <- "sim_hazard"

mat_haz1 <- matrix(c(lambda211, lambda131, lambda141, lambda251, lambda261, lambda531, lambda371, lambda581), nrow = 1)
colnames(mat_haz1) <- c('lambda211', 'lambda131', 'lambda141', 'lambda251', 'lambda261', 'lambda531', 'lambda371', 'lambda581')
rownames(mat_haz1) <- "sim_hazard"

mat_haz2 <- matrix(c(lambda212, lambda132, lambda142, lambda252, lambda262, lambda532, lambda372, lambda582), nrow = 1)
colnames(mat_haz2) <- c('lambda212', 'lambda132', 'lambda142', 'lambda252', 'lambda262', 'lambda532', 'lambda372', 'lambda582')
rownames(mat_haz2) <- "sim_hazard"

############ Helper Functions ##########

# A=1,V=0,D=1 -> A=1,V=1,D=1
#             -> A=0,V=1,D=0
long_dat_one_sim <- function(nobs, states, trates0, trates1, tlim, seed){
  
  set.seed(seed)
  
  # Generate long dataset:
  st <- matrix(0, nobs, 10)
  ss <- matrix(NA, nobs, 10)
  is <- rep('2', nobs)
  # TODO: 
  Ds <- matrix(NA, nobs, 10)
  Dinit <- matrix(NA, nobs, 1)
  
  # D_0 binary intercept, for level 2 (ref is 1, check D0=2 param)
  phis <- c(0) 
  
  # D multiordinal parameters, for level 1 and 2 (ref is 0, check D=1 & 2 params)
  intercepts <- c(-1.5, -2) 
  #intercepts <- c(-2, -1.5) 
  #phi1s <- c(1, 1.5)
  phi1s <- c(1, 0)
  num.cat <- 3
  
  sample_d <- function(X, phis, intercepts, phi1s, num.cat, init, Dinit_i){
    if (init==TRUE) {  # c(1,2)
      X <- 0
      nu <- phis
      p <- expit(nu)
      Dval <- sample(c(1,2), 1, prob = c(1-p, p)) #c(0.5, 0.5)
      #Dval <- sample(c(1:(num.cat-1)), 1, prob=prob.d[2:num.cat]) # This is for starting with multinom
    }else{             
      # D_0 effect on D
      beta_d0 <- c(0.5)#, 0.5)
      
      multireg <- intercepts + phi1s * X + beta_d0 * c(Dinit_i==1, Dinit_i==2)
      prob.d <- exp(multireg)
      prob.d <- cbind(1/(1+ sum(prob.d)), prob.d/(1+ matrix(sum(prob.d), 1, num.cat-1)))
      Dval <- sample(c(0,1,2), 1, prob=prob.d)                 
    }
    return(Dval)
  }
  
  for (i in 1:nobs) {
    # TODO: the initial dose to start with
    Dinit[i] <- sample_d(X=0, phis, intercepts, phi1s, num.cat, init = TRUE, Dinit_i = NA)
    # To start with D=Dinit matrix, might change later if a V happens
    if(Dinit[i] == 1){           
      trates = trates1; rates = rates1
    }else if(Dinit[i] == 2){
      trates = trates2; rates = rates2
    }
    
    current <- match(is[i], states)
    old <- match(is[i], states)
    counter <- 1
    st[i,counter] <- rexp(1, trates[current])
    ss[i,counter] <- sample(states, 1, prob=rates[current,]/trates[current])
    new <- match(ss[i,counter], states);new
    
    if (st[i,1:counter] > tlim) {
      ss[i,counter] <- 0
      st[i,counter] <- tlim
    }
    
    # TODO:
    if(ss[i,counter] == 1){  #1 state after 2->1, no X
      X <- 0
      Ds[i,counter] <- sample_d(X, phis, intercepts, phi1s, num.cat, init = FALSE, Dinit_i = Dinit[i])
    }else{
      Ds[i,counter] <- Dinit[i]
    }
    
    if(ss[i,counter] == 1 & Ds[i,counter] == 0){
      trates = trates0; rates = rates0
    }else if(ss[i,counter] == 1 & Ds[i,counter] == 1){
      trates = trates1; rates = rates1
    }else if(ss[i,counter] == 1 & Ds[i,counter] == 2){
      trates = trates2; rates = rates2
    }
    
    # All other transitions keep using the D=Dinit matrix
    ctime <- st[i,counter]; ctime
    while (ctime < tlim & !(ss[i,counter] %in% tstates)) {
      current <- match(ss[i,counter], states) 
      old <- match(ss[i,counter], states) 
      counter <- counter + 1
      st[i,counter] <- rexp(1, trates[current]) 
      ss[i,counter] <- sample(states, 1, prob=rates[current,]/trates[current]) 
      new <- match(ss[i,counter], states)
      ctime <- min(ctime+st[i,counter],tlim)  # == sum(st[i,1:counter]) 
      if (sum(st[i,1:counter]) > tlim) {
        ss[i,counter] <- 0
        st[i,counter] <- tlim - sum(st[i,1:(counter-1)])  
      }
      
      # TODO:
      if(ss[i,counter] == 3 & ss[i,counter-1] == 5){  #1 state after 5->3, where X has happened
        X <- 1
        Ds[i,counter] <- sample_d(X, phis, intercepts, phi1s, num.cat, FALSE, Dinit_i = Dinit[i])
      }else{
        Ds[i,counter] <- Ds[i,counter-1] 
      }
      
      if(ss[i,counter] == 3 & ss[i,counter-1] == 5 & Ds[i,counter] == 0){
        trates = trates0; rates = rates0
      }else if(ss[i,counter] == 3 & ss[i,counter-1] == 5 & Ds[i,counter] == 1){
        trates = trates1; rates = rates1
      }else if(ss[i,counter] == 3 & ss[i,counter-1] == 5 & Ds[i,counter] == 2){
        trates = trates2; rates = rates2
      }#else{} keep using the D=Ds[i, counter-1] matrix
      
    }
  }
  
  # Long-format dataset:
  idx <- matrix(rep(1:nobs, 10), nobs, 10, byrow=FALSE)
  
  long <- data.frame('idx'=as.numeric(t(idx)), 'futime'=as.numeric(t(st)), 'to'=as.integer(t(ss)), 'Dto'=as.integer(t(Ds)))
  long <- subset(long, !is.na(to)) 
  
  stc <- t(apply(st, 1, FUN=cumsum))
  long$from <- t(cbind(2,ss[,1:(ncol(st)-1)]))[t(!is.na(ss))]
  long$Dfrom <- t(cbind(Dinit, Ds[,1:(ncol(st)-1)]))[t(!is.na(Ds))]  
  long$start <- t(cbind(0,stc[,1:(ncol(stc)-1)]))[t(!is.na(ss))]
  long$stop <- t(stc)[t(!is.na(ss))]
  
  long <- long[, c('idx', 'futime', 'start', 'stop', 'from', 'to', 'Dfrom', 'Dto')]
  
  return(long)
}
#long <- long_dat_one_sim(nobs, states, trates0, trates1, tlim, seed)

## Helpers for Generating V exposure model and Y outcome model dataset
# Here at the same time generated an incomplete long_v for times wide dataset and then complete long_v
times_wide_generation <- function(long, nobs){
  long$vt <- ifelse((long$from == 2 & long$to == 1) | (long$from == 5 & long$to == 3), long$stop, Inf)
  long$xt <- ifelse((long$from == 2 & long$to == 5) | (long$from == 1 & long$to == 3), long$stop, Inf)
  long$yt <- ifelse(long$to %in% c(4,6,7,8), long$stop, Inf)
  long$dt <- long$vt
  
  times_wide <- data.frame('idx'=as.numeric(c(1:nobs)))
  for(i in 1:nobs){
    idx <- long$idx == i
    times_wide$vt[i] <- sort(long$vt[idx])[1]  # later on if backward error allowed, then vt2 will be [2]. Add if else condition for only one vt.
    times_wide$xt[i] <- sort(long$xt[idx])[1]
    times_wide$yt[i] <- sort(long$yt[idx])[1]
    times_wide$dt[i] <- times_wide$vt[i]
  }
  
  return(times_wide)
}

# times <- times_wide_generation(long, nobs)
# dtimes_V <- seq(from = 0, to = tlim, length.out = 61)[-1] # monthly interval
get_long_V_slm <- function(times, dtimes_V, long){
  long_V_slm <- times[rep(1:nobs, times=rep(length(dtimes_V), nobs)),] # rep 1:1000, each 100 (dtimes) times 
  long_V_slm$stop <- dtimes_V                                          # automatic replicate dtimes when finished. 
  long_V_slm <- long_V_slm[order(long_V_slm$id, long_V_slm$stop),]            # actually useless, already ordered
  long_V_slm$Xconf <- (long_V_slm$xt <= long_V_slm$stop)
  long_V_slm$start <- NA
  long_V_slm$Visit <- NA
  long_V_slm$Dtrt <- NA
  
  for (i in 1:nobs) {
    idx <- long_V_slm$id == i
    long_V_slm$start[idx] <- c(0, dtimes_V[1:(length(dtimes_V)-1)])
    long_V_slm$Visit[idx] <- ifelse(times$vt[i] <= long_V_slm$stop[idx], 1, 0)   # Keep V 0->1
    long_V_slm$month[idx] <- c(1:length(dtimes_V))
    long_V_slm$Dtrt[idx] <- ifelse(times$dt[i] <= long_V_slm$stop[idx], 
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 
                                   long$Dfrom[long$idx == i][1])  
    long_V_slm$Dprev[idx] <- long[long$idx == i, ]$Dfrom[1]
    # because start with V=0 & D = Dinit
    
    # if (i %% 100 == 0)
    #     print(i)
  }
  long_V_slm <- subset(long_V_slm, start <= vt)
  long_V_slm <- subset(long_V_slm, start <= yt)
  #long_V_slm$stop <- ifelse(long_V_slm$stop > long_V_slm$yt, long_V_slm$yt, long_V_slm$stop)
  long_V_slm$futime <- long_V_slm$stop - long_V_slm$start
  
  return(long_V_slm)
}

#long_V_slm <- get_long_V_slm(times, dtimes_V, long)