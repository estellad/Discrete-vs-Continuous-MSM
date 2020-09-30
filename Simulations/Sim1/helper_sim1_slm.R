# This is everything needed to run the scripts: multistate_3_outcome_cov.R and multistate_subset_each_lambda.R
# and sim_repeat.R

#####  Default Parameters  #####
library(survival)
nobs <- 1000
tlim <- 5       # 1

alpha_a <- 0.1  # 0.8
alpha_x <- 0.05 # 0.4
alpha_y <- 0.02  # 0.2

beta_a <- 1.5   # 0.6
beta_x <- -1.5  # 0.1
beta_y <- -1    # 0.7
gamma_y <- 1.5  # 0.9

# lambdaA
lambda12 <- alpha_a; lambda12
# lambdaX
lambda13 <- alpha_x; lambda13
# lambdaY
lambda14 <- alpha_y; lambda14
# lambdaX
lambda25 <- alpha_x*exp(beta_x); lambda25
# lambdaY
lambda26 <- alpha_y*exp(beta_y); lambda26
# lambdaA
lambda35 <- alpha_a*exp(beta_a); lambda35
# lambdaY
lambda37 <- alpha_y*exp(gamma_y); lambda37
# lambdaY
lambda58 <- alpha_y*exp(beta_y + gamma_y); lambda58


rates <- matrix(c(0, lambda12, lambda13, lambda14,        0,        0,        0,        0,
                  0,        0,        0,        0, lambda25, lambda26,        0,        0,
                  0,        0,        0,        0, lambda35,        0, lambda37,        0,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0,        0,        0,        0,        0,        0, lambda58,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0,        0,        0,        0,        0,        0,        0,
                  0,        0,        0,        0,        0,        0,        0,        0), 
                8, 8, byrow=TRUE)
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

mat_haz <- matrix(c(lambda12, lambda13, lambda14, lambda25, lambda26, lambda35, lambda37, lambda58), nrow = 1)
colnames(mat_haz) <- c('lambda12', 'lambda13', 'lambda14', 'lambda25', 'lambda26', 'lambda35', 'lambda37', 'lambda58')
rownames(mat_haz) <- "sim_hazard"

############ Helper Functions ##########

## Common 1 simulation long dataset at 1 seed. for two methods.
long_dat_one_sim <- function(nobs, states, trates, tlim, seed){
  
  set.seed(seed)
  
  # Generate long dataset:
  st <- matrix(0, nobs, 10)
  ss <- matrix(NA, nobs, 10)
  is <- rep('1', nobs)
  
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
  long$from <- t(cbind(1,ss[,1:(ncol(st)-1)]))[t(!is.na(ss))]
  long$start <- t(cbind(0,stc[,1:(ncol(stc)-1)]))[t(!is.na(ss))]
  long$stop <- t(stc)[t(!is.na(ss))]
  
  return(long)
}
#long <- long_dat_one_sim(nobs, states, trates, tlim, seed)

## Helpers for Generating A exposure model and Y outcome model dataset
# Here at the same time generated an incomplete long_A for times wide dataset and then complete long_A
times_wide_generation <- function(long, nobs){
  long$at <- ifelse((long$from == 1 & long$to == 2) | (long$from == 3 & long$to == 5), long$stop, Inf)
  long$xt <- ifelse((long$from == 1 & long$to == 3) | (long$from == 2 & long$to == 5), long$stop, Inf)
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

# times <- times_wide_generation(long, nobs)
# dtimes_A <- seq(from = 0, to = tlim, length.out = 61)[-1]
get_long_A_slm <- function(times, dtimes_A){
  long_A_slm <- times[rep(1:nobs, times=rep(length(dtimes_A), nobs)),] # rep 1:1000, each 100 (dtimes) times 
  long_A_slm$stop <- dtimes_A                                          # automatic replicate dtimes when finished. 
  long_A_slm <- long_A_slm[order(long_A_slm$id, long_A_slm$stop),]            # actually useless, already ordered
  long_A_slm$Xconf <- (long_A_slm$xt <= long_A_slm$stop)
  
  for (i in 1:nobs) {
    idx <- long_A_slm$id == i
    long_A_slm$start[idx] <- c(0, dtimes_A[1:(length(dtimes_A)-1)])
    long_A_slm$Aexpo[idx] <- ifelse(times$at[i] <= long_A_slm$stop[idx], 1, 0)
    long_A_slm$month[idx] <- c(1:length(dtimes_A))
    # if (i %% 100 == 0)
    #     print(i)
  }
  #long_A_slm <- subset(long_A_slm, stop <= at)
  long_A_slm <- subset(long_A_slm, start < at)
  long_A_slm <- subset(long_A_slm, start < yt)
  # We only discretize at pre-specified time points, so we don't shorten the time interval, but add keep and extend the current one.
  # long_A_slm$stop <- ifelse(long_A_slm$stop > long_A_slm$yt, long_A_slm$yt, long_A_slm$stop)
  long_A_slm$futime <- long_A_slm$stop - long_A_slm$start
  
  return(long_A_slm)
}

# long_A_slm <- get_long_A_slm(times, dtimes_A)
# 
# ###################### Logistic Exposure Model Fitting ################## 
# # Marginal
# model_logi_A_slm_marginal <- glm(Aexpo ~ 1 + offset(log(futime)), family=binomial(link=logit), data=long_A_slm) # Marginal
# coef(model_logi_A_slm_marginal)
# 
# # Conditional
# model_logi_A_slm_conditional <- glm(Aexpo ~ Xconf+month + offset(log(futime)), family=binomial(link=logit), data=long_A_slm) # Conditional
# coef(model_logi_A_slm_conditional)

