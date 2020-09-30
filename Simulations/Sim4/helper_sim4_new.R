# This is everything needed to run the scripts: multistate_3_outcome_cov.R and multistate_subset_each_lambda.R
# and sim_repeat.R

#####  Default Parameters  #####
library(survival)
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

### No positivity violation

# num_event <- vector()
# num_x <- vector()
# num_v <- vector()
# x_v_both <- vector()
# x_before_v <- vector()
# x_before_d0 <- vector()
# x_before_d1 <- vector()
# x_before_d2 <- vector()
# x_d_both <- vector()
# x_before_dnum <- vector()
# length_x_before_d <- vector()
# 
# for(i in 1:nobs){
#   long_i <- long_dat_one_sim(nobs, states, trates0, trates1, tlim, seed = i)
#   times_i <- times_wide_generation(long_i, nobs)
#   y_events_i <- nrow(times_i[is.finite(times_i$yt),])
#   num_event[i] <- y_events_i
# 
#   x_events_i <- nrow(times_i[is.finite(times_i$xt),])
#   num_x[i] <- x_events_i
# 
#   v_events_i <- nrow(times_i[is.finite(times_i$vt),])
#   num_v[i] <- v_events_i
# 
#   x_v_times_i <- times_i[is.finite(times_i$xt) & is.finite(times_i$vt), ]
#   x_d_times_i <- times_i[is.finite(times_i$xt) & is.finite(times_i$dt), ]
#   x_v_both[i] <- nrow(x_v_times_i)
#   x_d_both[i] <- nrow(x_d_times_i)
# 
#   x_before_v[i] <- nrow(x_v_times_i[x_v_times_i$xt<x_v_times_i$vt, ])
#   x_before_dnum[i] <- nrow(x_d_times_i[x_d_times_i$xt<x_d_times_i$dt, ])
# 
#   x_before_d <- long_i[long_i$idx %in% x_v_times_i[x_v_times_i$xt<x_v_times_i$vt, ]$idx, ]
#   length_x_before_d[i] <- length(unique(x_before_d$idx))
#   x_before_d0[i] <- nrow(x_before_d[x_before_d$from == 5 & x_before_d$to == 3 & x_before_d$Dto == 0, ])
#   x_before_d1[i] <- nrow(x_before_d[x_before_d$from == 5 & x_before_d$to == 3 & x_before_d$Dto == 1, ])
#   x_before_d2[i] <- nrow(x_before_d[x_before_d$from == 5 & x_before_d$to == 3 & x_before_d$Dto == 2, ])
# }
# summary(num_event)
# summary(num_x)
# summary(num_v)
# summary(x_v_both)
# summary(x_d_both)
# 
# summary(x_before_v)
# summary(x_before_dnum)
# summary(length_x_before_d)
# 
# summary(x_before_d0)
# summary(x_before_d1)
# summary(x_before_d2)
#
#
# > summary(num_event)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 42.00   62.00   68.00   67.77   73.00   92.00 
# > summary(num_x)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501.0   538.0   549.0   548.5   559.0   588.0 
# > summary(num_v)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 464.0   508.0   519.0   518.9   530.0   564.0 
# > summary(x_v_both)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 291.0   338.0   348.0   348.2   359.0   388.0 
# > summary(x_d_both)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 291.0   338.0   348.0   348.2   359.0   388.0 
# > 
#   > summary(x_before_v)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 246.0   286.8   297.0   296.6   307.0   343.0 
# > summary(x_before_dnum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 246.0   286.8   297.0   296.6   307.0   343.0 
# > summary(length_x_before_d)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 246.0   286.8   297.0   296.6   307.0   343.0 
# > 
#   > summary(x_before_d0)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 53.00   73.00   79.00   78.46   84.00  106.00 
# > summary(x_before_d1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 32.00   54.00   59.00   59.13   64.00   92.00 
# > summary(x_before_d2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 127     151     159     159     166     203 


#times <- times_wide_generation(long, nobs)
#dtimes_V <- sort(unique(times$vt[is.finite(times$vt)]))
get_long_V_ctm <- function(times, dtimes_V, long){
  long_V_ctm <- times[rep(1:nobs, times=rep(length(dtimes_V), nobs)),] # rep 1:1000, each 100 (dtimes) times 
  long_V_ctm$stop <- dtimes_V                                          # automatic replicate dtimes when finished. 
  long_V_ctm <- long_V_ctm[order(long_V_ctm$id, long_V_ctm$stop),]            # actually useless, already ordered
  long_V_ctm$Xconf <- (long_V_ctm$xt <= long_V_ctm$stop)
  long_V_ctm$start <- NA
  long_V_ctm$Visit <- NA
  long_V_ctm$Dtrt <- NA
  
  for (i in 1:nobs) {
    idx <- long_V_ctm$id == i
    long_V_ctm$start[idx] <- c(0, dtimes_V[1:(length(dtimes_V)-1)])
    long_V_ctm$Visit[idx] <- ifelse(times$vt[i] == long_V_ctm$stop[idx], 1, 0)   # V 1->0 
    long_V_ctm$Dtrt[idx] <- ifelse(times$dt[i] == long_V_ctm$stop[idx], 
                                   long$Dto[long$idx == i][long$stop[long$idx == i] == times$dt[i]], 
                                   long$Dfrom[long$idx == i][1])  
    long_V_ctm$Dprev[idx] <- long[long$idx == i, ]$Dfrom[1]
    # because start with V=0 & D = Dinit
    
    # if (i %% 100 == 0)
    #     print(i)
  }
  long_V_ctm <- subset(long_V_ctm, stop <= vt)
  long_V_ctm <- subset(long_V_ctm, start <= yt)
  long_V_ctm$stop <- ifelse(long_V_ctm$stop > long_V_ctm$yt, long_V_ctm$yt, long_V_ctm$stop)
  long_V_ctm$futime <- long_V_ctm$stop - long_V_ctm$start
  
  return(long_V_ctm)
}

#long_V_ctm <- get_long_V_ctm(times, dtimes_V, long)


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
  
  long_X_ctm$Dtrt <- as.factor(long_X_ctm$Dtrt)
  long_X_ctm$Dtrt <- relevel(long_X_ctm$Dtrt, ref = "0")
  
  return(long_X_ctm)
}

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

checkAYtrue <- function(seed, m, nobs, states, trates, tlim){
  long <- long_dat_one_sim(nobs, states, trates0=trates0, trates1 = trates1, tlim, seed)
  times <- times_wide_generation(long, nobs)
  dtimes_V <- sort(unique(times$vt[is.finite(times$vt)]))  # Stop of treatment time
  
  ## Exposure Data: V as outcome, X as cov (delay one interval)
  long_V_ctm <- get_long_V_ctm(times, dtimes_V, long)
  
  ################## Multinomial Logistic D Treatment Model Fitting: Changed Dose #############
  # Treatment Model after visit
  D_expo_data <- long_V_ctm[long_V_ctm$stop == long_V_ctm$dt, ]
  
  # Marginal
  model_logi_V_ctm_marginal_d <- multinom(Dtrt ~ 1, data=D_expo_data)
  # coef(model_logi_V_ctm_marginal_d)
  
  # Conditional
  model_logi_V_ctm_conditional_d <- multinom(Dtrt ~ Xconf+Dprev, data=D_expo_data)
  #coef(model_logi_V_ctm_conditional_d)
  
  ####################################### Outcome model ######################################
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
  
  ytruedtab <- table(casebase$A, casebase$Yevent)[,2] # 1,0,2
  ytruedtab <- ytruedtab[c(2,1,3)] # 0, 1, 2
  dtab <- table(D_expo_data$Dtrt)  # 0, 1, 2
  tabs <- as.numeric(c(ytruedtab, dtab)) # c(20,  25,  13, 225, 187, 107)

  marginal_result <- list(tabs)
  return(marginal_result)
}