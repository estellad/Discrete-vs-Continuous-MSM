# This is everything needed to run the scripts: multistate_3_outcome_cov.R and multistate_subset_each_lambda.R
# and sim_repeat.R

#####  Default Parameters  #####
library(survival)
nobs <- 1000
# tlim <- 1
# 
# alpha_a <- 0.8 
# alpha_x <- 0.4 
# alpha_y <- 0.2
# 
# beta_a <- 0.6
# beta_x <- 0.1
# beta_y <- 0.7
# gamma_y <- 0.9 
# 
# # lambdaA
# lambda12 <- exp(alpha_a); lambda12
# # lambdaX
# lambda13 <- exp(alpha_x); lambda13
# # lambdaY
# lambda14 <- exp(alpha_y); lambda14
# # lambdaX
# lambda25 <- exp(alpha_x + beta_x); lambda25
# # lambdaY
# lambda26 <- exp(alpha_y + beta_y); lambda26
# # lambdaA
# lambda35 <- exp(alpha_a + beta_a); lambda35
# # lambdaY
# lambda37 <- exp(alpha_y + gamma_y); lambda37
# # lambdaY
# lambda58 <- exp(alpha_y + beta_y + gamma_y); lambda58

tlim <- 5       # 1

# Original A=0 -> A=1 parameters:
# alpha_a <- 0.1  # 0.8
# alpha_x <- 0.05 # 0.4
# alpha_y <- 0.02  # 0.2
# 
# beta_a <- 1.5   # 0.6
# beta_x <- -1.5  # 0.1
# beta_y <- -1    # 0.7
# gamma_y <- 1.5  # 0.9

# Me trying: to increase Xt before at events
# alpha_a <- 0.08  # 0.8
# alpha_x <- 0.05 # 0.4
# alpha_y <- 0.04 # 0.2
# 
# beta_a <- -1.5   # 0.6
# beta_x <- 1  # 0.1
# beta_y <- -1    # 0.7        # True  Effect
# gamma_y <- 1.5  # 0.9

# Olli suggestion: 
alpha_a <- 0.2
alpha_x <- 0.05
alpha_y <- 0.02

beta_a <- 1.5
beta_x <- 1.5
beta_y <- -1
gamma_y <- 1.5

# # lambdaA
# lambda12 <- alpha_a; lambda12
# # lambdaX
# lambda13 <- alpha_x; lambda13
# # lambdaY
# lambda14 <- alpha_y; lambda14
# # lambdaX
# lambda25 <- alpha_x*exp(beta_x); lambda25
# # lambdaY
# lambda26 <- alpha_y*exp(beta_y); lambda26
# # lambdaA
# lambda35 <- alpha_a*exp(beta_a); lambda35
# # lambdaY
# lambda37 <- alpha_y*exp(gamma_y); lambda37
# # lambdaY
# lambda58 <- alpha_y*exp(beta_y + gamma_y); lambda58

# lambdaA                          # This is the graph without changing state name, starting = 2
lambda21 <- alpha_a; lambda21
# lambdaX
lambda13 <- alpha_x; lambda13
# lambdaY
lambda14 <- alpha_y; lambda14
# lambdaX
lambda25 <- alpha_x*exp(beta_x); lambda25
# lambdaY
lambda26 <- alpha_y*exp(beta_y); lambda26
# lambdaA
lambda53 <- alpha_a*exp(beta_a); lambda53
# lambdaY
lambda37 <- alpha_y*exp(gamma_y); lambda37
# lambdaY
lambda58 <- alpha_y*exp(beta_y + gamma_y); lambda58

rates <- matrix(c(0,        0, lambda13, lambda14,        0,        0,        0,        0,
                  lambda21, 0,        0,        0, lambda25, lambda26,        0,        0,
                  0,        0,        0,        0,        0,        0, lambda37,        0,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0, lambda53,        0,        0,        0,        0, lambda58,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0,        0,        0,        0,        0,        0,        0), 
                8, 8, byrow=TRUE)

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

rownames(rates) <- colnames(rates) <- c('1', '2', '3', '4', '5', '6', '7', '8')
trates <- rowSums(rates)

states <- as.numeric(colnames(rates)) # c(1, 2, 3, 4, 5, 6, 7, 8)
numstates <- length(states)
tstates <- c(4, 6, 7, 8) 
cols <- c('gray','green','red','black', 'orange', 'blue', 'coral','cyan')

# input_coef <- matrix(c(alpha_a, alpha_x, alpha_y, alpha_x+beta_x, alpha_y+beta_y, alpha_a+beta_a, alpha_y+gamma_y, alpha_y+gamma_y+beta_y), nrow=1); input_coef
# colnames(input_coef) <- c('alpha_a', 'alpha_x', 'alpha_y', 'alpha_x+beta_x', 'alpha_y+beta_y', 'alpha_a+beta_a', 'alpha_y+gamma_y', 
#                        'alpha_y+gamma_y+beta_y')
# rownames(input_coef) <- "input_coef"

# Did not change state name: starting = 2
mat_haz <- matrix(c(lambda21, lambda13, lambda14, lambda25, lambda26, lambda53, lambda37, lambda58), nrow = 1)
colnames(mat_haz) <- c('lambda21', 'lambda13', 'lambda14', 'lambda25', 'lambda26', 'lambda53', 'lambda37', 'lambda58')
rownames(mat_haz) <- "sim_hazard"

############ Helper Functions ##########

## Common 1 simulation long dataset at 1 seed. for two methods.
long_dat_one_sim <- function(nobs, states, trates, tlim, seed){
  
  set.seed(seed)
  
  # Generate long dataset:
  st <- matrix(0, nobs, 10)
  ss <- matrix(NA, nobs, 10)
  is <- rep('2', nobs)
  
  for (i in 1:nobs) {
    current <- match(is[i], states)
    old <- match(is[i], states)
    counter <- 1
    st[i,counter] <- rexp(1, trates[current])
    ss[i,counter] <- sample(states, 1, prob=rates[current,]/trates[current])
    new <- match(ss[i,counter], states)
    
    if (st[i,1:counter] > tlim) {
      ss[i,counter] <- 0
      st[i,counter] <- tlim
    }
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
    }
  }
  
  # Long-format dataset:
  idx <- matrix(rep(1:nobs, 10), nobs, 10, byrow=FALSE)
  
  long <- data.frame('idx'=as.numeric(t(idx)), 'futime'=as.numeric(t(st)), 'to'=as.integer(t(ss)))
  long <- subset(long, !is.na(to)) 
  
  stc <- t(apply(st, 1, FUN=cumsum))
  long$from <- t(cbind(2,ss[,1:(ncol(st)-1)]))[t(!is.na(ss))]
  long$start <- t(cbind(0,stc[,1:(ncol(stc)-1)]))[t(!is.na(ss))]
  long$stop <- t(stc)[t(!is.na(ss))]
  
  return(long)
}

# long <- long_dat_one_sim(nobs, states, trates, tlim, seed = 2)


## Helpers for Generating A exposure model and Y outcome model dataset
# Here at the same time generated an incomplete long_A for times wide dataset and then complete long_A
times_wide_generation <- function(long, nobs){
  long$at <- ifelse((long$from == 2 & long$to == 1) | (long$from == 5 & long$to == 3), long$stop, Inf)
  long$xt <- ifelse((long$from == 2 & long$to == 5) | (long$from == 1 & long$to == 3), long$stop, Inf)
  long$yt <- ifelse(long$to %in% c(4,6,7,8), long$stop, Inf)
  
  times_wide <- data.frame('idx'=as.numeric(c(1:nobs)))
  for(i in 1:nobs){
    idx <- long$idx == i
    times_wide$at[i] <- sort(long$at[idx])[1]  # later on if backward error allowed, then at2 will be [2]. Add if else condition for only one at.
    times_wide$xt[i] <- sort(long$xt[idx])[1]
    times_wide$yt[i] <- sort(long$yt[idx])[1]
  }
  return(times_wide)
}

# num_event <- NULL
# num_x <- NULL
# num_a <- NULL
# x_a_both <- NULL
# x_before_a <- NULL
# for(i in 1:100){
#   long_i <- long_dat_one_sim(nobs, states, trates, tlim, seed = i)
#   times_i <- times_wide_generation(long_i, nobs)
#   y_events_i <- nrow(times_i[is.finite(times_i$yt),])
#   num_event[i] <- y_events_i
# 
#   x_events_i <- nrow(times_i[is.finite(times_i$xt),])
#   num_x[i] <- x_events_i
# 
#   a_events_i <- nrow(times_i[is.finite(times_i$at),])
#   num_a[i] <- a_events_i
# 
#   x_a_times_i <- times_i[is.finite(times_i$xt) & is.finite(times_i$at), ]
#   x_a_both[i] <- nrow(x_a_times_i)
#   # 100
#   x_before_a[i] <- nrow(x_a_times_i[x_a_times_i$xt<x_a_times_i$at, ])
# }
# summary(num_event)
# summary(num_x)
# summary(num_a)
# summary(x_a_both)
# summary(x_before_a)

# With parameter setting in A=0->1
# > summary(num_event)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 80.0    95.0   102.0   101.6   108.0   124.0
# > summary(num_x)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 157.0   171.8   181.0   180.7   189.0   204.0 
# > summary(num_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 395     426     436     437     448     485
# > summary(x_a_both)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 91.0   104.0   110.5   110.2   117.0   137.0 
# > summary(x_before_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 74.00   93.00   99.00   99.68  106.25  129.00 

# My new parameter setting
# > summary(num_event)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 120.0   142.0   150.5   150.3   160.0   178.0 
# > summary(num_x)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 396.0   423.0   432.5   433.8   445.0   471.0 
# > summary(num_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 225.0   247.0   255.0   256.6   266.0   288.0 
# > summary(x_a_both)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.00   43.75   49.00   48.87   53.25   72.00 
# > summary(x_before_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.00   15.00   18.00   18.32   21.00   29.00 

# Olli's suggested parameter setting
# > summary(num_event)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 110.0   137.0   146.0   145.4   152.0   175.0 
# > summary(num_x)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 480.0   510.0   520.0   519.9   530.2   555.0 
# > summary(num_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 785.0   799.0   808.0   808.0   816.2   835.0 
# > summary(x_a_both)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 422.0   447.0   460.0   457.7   466.0   493.0 
# > summary(x_before_a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 361.0   387.0   397.5   396.9   406.0   430.0 

 
# times <- times_wide_generation(long, nobs)
# 
# nrow(times[is.finite(times$yt),])
# # 165
# nrow(times[is.finite(times$at),])
# # 811
# nrow(times[is.finite(times$xt),])
# # 518
# x_a_times <- times[is.finite(times$xt) & is.finite(times$at), ]
# nrow(x_a_times) # X A both happen
# # 454
# nrow(x_a_times[x_a_times$xt<x_a_times$at, ])  # X happen before A
# # 398
# 
# dtimes_A <- sort(unique(times$at[is.finite(times$at)]))  # Stop of treatment time

get_long_A_ctm <- function(times, dtimes_A){
  long_A_ctm <- times[rep(1:nobs, times=rep(length(dtimes_A), nobs)),] # rep 1:1000, each 100 (dtimes) times
  long_A_ctm$stop <- dtimes_A                                          # automatic replicate dtimes when finished.
  long_A_ctm <- long_A_ctm[order(long_A_ctm$id, long_A_ctm$stop),]            # actually useless, already ordered
  long_A_ctm$Xconf <- (long_A_ctm$xt <= long_A_ctm$stop)

  for (i in 1:nobs) {
    idx <- long_A_ctm$id == i
    long_A_ctm$start[idx] <- c(0, dtimes_A[1:(length(dtimes_A)-1)])
    long_A_ctm$Aexpo[idx] <- ifelse(times$at[i] == long_A_ctm$stop[idx], 0, 1)
    # if (i %% 100 == 0)
    #     print(i)
  }
  long_A_ctm <- subset(long_A_ctm, stop <= at)
  long_A_ctm <- subset(long_A_ctm, start <= yt)
  long_A_ctm$stop <- ifelse(long_A_ctm$stop > long_A_ctm$yt, long_A_ctm$yt, long_A_ctm$stop)
  long_A_ctm$futime <- long_A_ctm$stop - long_A_ctm$start

  return(long_A_ctm)
}

# long_A_ctm <- get_long_A_ctm(times, dtimes_A)
# 
# # Marginal
# # Exposure model should invert, but in outcome model do not invert A
# model_poi_A_ctm_marginal <- glm(!(Aexpo) ~ 1 + offset(log(futime)), family=poisson(link=log), data=long_A_ctm)
# coef(model_poi_A_ctm_marginal)
# # Interpretation: hazard of stop the treatment, survival prob of being on treatment. Before, it was hazard of start the treatment.
# 
# # Conditional
# model_poi_A_ctm_conditional <- glm(!(Aexpo) ~ Xconf + offset(log(futime)), family=poisson(link=log), data=long_A_ctm)
# coef(model_poi_A_ctm_conditional)

