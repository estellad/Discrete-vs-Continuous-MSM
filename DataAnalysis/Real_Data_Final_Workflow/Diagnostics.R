if(female){
  pdf('/users/edong/Output/Female/Diagnostics/densityfemale_diffscale_qs.pdf', width=13, height=10)
}else{
  pdf('/users/edong/Output/Male/Diagnostics/densitymale_diffscale_qs.pdf', width=13, height=10) 
}
pltm <- rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5)
              #, c(6,6,6,7,7,7)
)
layout(pltm)
par(mar = c(3,4,4,1))

plot(density(long_Y_cb$sw), main = paste('Density Plot of Combined Stabilized Weight \n for Exposure A -', gender.name), 
     xlab="Inverse Probability Weights", col='darkred')
rug(long_Y_cb$sw, col='darkred')

plot(density(long_Y_cb$sw.trunc), main = paste('Density Plot of Combined Stabilized Weight \n Truncated at 0.5% and 99.5% for Exposure A -', gender.name), 
     xlab="Inverse Probability Weights", col='darkred', xlim=c(floor(qsl), ceil(qsu)))
rug(long_Y_cb$sw, col='darkred')
abline(v=qsl, lty='dotted', col='darkred', lwd=2)
abline(v=qsu, lty='dotted', col='darkred', lwd=2)

plot(density(long_Y_cb$d0p), main = paste('Density Plot of Initial Dose sw \n for Exposure D0 -', gender.name), 
     xlab="Inverse Probability Weights")#, xlim=xrange)
rug(long_Y_cb$d0p)

plot(density(long_Y_cb$vp), main = paste('Density Plot of Subsequent Visiting sw \n for Exposure V -', gender.name), 
     xlab="Inverse Probability Weights")#, xlim=xrange)
rug(long_Y_cb$vp)

plot(density(long_Y_cb$dp), main = paste('Density Plot of Subsequent Dose sw \n for Exposure D -', gender.name), 
     xlab="Inverse Probability Weights")#, xlim=xrange)
rug(long_Y_cb$dp)

# ## Two D models based on Dtrtp=0 or >0 status
# plot(density(long_Y_cb$dp0p), main = paste('Density Plot of Subsequent Dose sw \n where Prev. Dose is 0 -', gender.name), 
#      xlab="Inverse Probability Weights")#, xlim=xrange)
# rug(long_Y_cb$dp0p)
# 
# plot(density(long_Y_cb$dpnot0p), main = paste('Density Plot of Subsequent Dose sw \n where Prev. Dose is not 0 -', gender.name), 
#      xlab="Inverse Probability Weights")#, xlim=xrange)
# rug(long_Y_cb$dpnot0p)

dev.off()


# library(lattice)
# densityplot(~weight | weighttype, data= superlongswdenscompare)

########################## Density percentage interval over time plot ###################
swinterval <- subset(long_Y_cb,select = c(sw, sw.trunc, End)) %>%
  arrange(End) %>% 
  mutate(index = c(1:nrow(long_Y_cb)))

cutoff <- head(seq(0,365, by=stepsize), -1)
lab <- NULL
for(i in 1:length(cutoff)){
  lab[i] <- which(swinterval$End > cutoff[i])[1]
}
length(lab); head(lab); vec <- lab-lag(lab); min(vec[-1])

swinterval$tp <- ifelse(swinterval$index %in% lab, 'a', 1)
swinterval$tp <- cumsum(swinterval$tp == 'a'); swinterval$tp <- as.factor(swinterval$tp)
swinterval <- subset(swinterval, select = -index)


# time 0 weight 
tp0sw <- data.frame(
  tp=0,
  minsw=min(long_Y_cb$d0p),
  maxsw=max(long_Y_cb$d0p),
  q1=quantile(long_Y_cb$d0p, probs = 0.25),
  q3=quantile(long_Y_cb$d0p, probs = 0.75),
  q10=quantile(long_Y_cb$d0p, probs = 0.05),
  q90=quantile(long_Y_cb$d0p, probs = 0.95),
  median= quantile(long_Y_cb$d0p, probs = 0.50)
)

logtp0sw <- log(tp0sw); logtp0sw[1,1] <- 0


if(trunc){
  # trunc
  swintervalsumlog <-  swinterval %>% 
    group_by(tp) %>%
    summarise(
      minsw = log(min(sw.trunc)),
      maxsw = log(max(sw.trunc)),
      q1 = log(quantile(sw.trunc, probs = 0.25)),
      q3 = log(quantile(sw.trunc, probs = 0.75)),
      q10 = log(quantile(sw.trunc, probs = 0.05)),
      q90 = log(quantile(sw.trunc, probs = 0.95)),
      median = log(quantile(sw.trunc, probs = 0.50)))
}else{
  # no truncation
  swintervalsumlog <-  swinterval %>% 
    group_by(tp) %>%
    summarise(
      minsw = log(min(sw)),
      maxsw = log(max(sw)),
      q1 = log(quantile(sw, probs = 0.25)),
      q3 = log(quantile(sw, probs = 0.75)),
      q10 = log(quantile(sw, probs = 0.05)),
      q90 = log(quantile(sw, probs = 0.95)),
      median = log(quantile(sw, probs = 0.50)))
}
swintervalsumlog <- rbind(logtp0sw, swintervalsumlog) %>%
  mutate(step = stepsize) %>%
  mutate(step = cumsum(step)) %>%
  mutate(step = lag(step)) %>%
  mutate(step = ifelse(is.na(step), 0, step))

# swintervalsum <-  swinterval %>% 
#   group_by(tp) %>%
#   summarise(
#     minsw = min(sw.trunc),
#     maxsw = max(sw.trunc),
#     q1 = quantile(sw.trunc, probs = 0.25),
#     q3 = quantile(sw.trunc, probs = 0.75),
#     q10 = quantile(sw.trunc, probs = 0.05),
#     q90 = quantile(sw.trunc, probs = 0.95),
#     median = quantile(sw.trunc, probs = 0.50)) 
# swintervalsum <- rbind(tp0sw, swintervalsum)

# log graph 0.05
# log graph 0.1
# log graph 0.5
# log graph 1
if(female & trunc){
  pdf('/users/edong/Output/Female/Diagnostics/swtimefemale_trunc.pdf', width=6.5, height=6.5)
  titleplt <- 'Estimated Combined Stabilized Weights \n Truncated at 0.5% and 99.5% Over Time - '
}else if(female & !trunc){
  pdf('/users/edong/Output/Female/Diagnostics/swtimefemale.pdf', width=6.5, height=6.5)
  titleplt <- 'Estimated Combined Stabilized Weights Over Time - '
}else if(!female & trunc){
  pdf('/users/edong/Output/Male/Diagnostics/swtimemale_trunc.pdf', width=6.5, height=6.5)
  titleplt <- 'Estimated Combined Stabilized Weights \n Truncated at 0.5% and 99.5% Over Time - '
}else{
  pdf('/users/edong/Output/Male/Diagnostics/swtimemale.pdf', width=6.5, height=6.5)
  titleplt <- 'Estimated Combined Stabilized Weights Over Time - '
}
par(mfrow=c(1,1))
# yrange <- range(c(swintervalsumlog$minsw, swintervalsumlog$maxsw)); yrange
# ymaxrange <- max(abs(floor(yrange[1]*2)/2), ceil(yrange[2]*2)/2)
# yrange <- c(-ymaxrange, ymaxrange)
yrange <- c(-1.5, 1.5)
plot(median ~ step, lty=1, main = paste0(titleplt, gender.name),
     xlab='Days since start of follow-up', ylab='IPTW (log-scale)', lwd=2, type='s',
     xlim = c(0, max(swintervalsumlog$step)), ylim = yrange, yaxt= 'n',
     data = swintervalsumlog)
lines(minsw ~ step, data = swintervalsumlog, type = 's', lty=2, col='red', lwd=1.2)
lines(maxsw ~ step, data = swintervalsumlog, type = 's', lty=2, col='red', lwd=1.2)
lines(q1 ~ step, data = swintervalsumlog, type = 's', lty=3, col='blue', lwd=1.2)
lines(q3 ~ step, data = swintervalsumlog, type = 's', lty=3, col='blue', lwd=1.2)
lines(q10 ~ step, data = swintervalsumlog, type = 's', lty=4, col='purple', lwd=1.2)
lines(q90 ~ step, data = swintervalsumlog, type = 's', lty=4, col='purple', lwd=1.2)

vec <- seq(yrange[1], yrange[2], length.out= 5)
axis(side = 2, at = vec, 
     labels = round(exp(vec),3))
legend('bottomleft', c('Median', '50% Interval', '90% Interval', 'Min./Max.'), cex = 0.7, 
       col=c('black', 'blue', 'purple', 'red'), lty=c(1,3,4,2), lwd=1.5)
dev.off()



######## log graph without step, just lines to linear interpolate the points.
# pdf('swtimemale_interpolate005.pdf', width=6.5, height=6.5)
# yrange = c(-2,2)
# plot(x=c(0,cutoff), swintervalsumlog$median, lty=1, main = 'Estimated stabilized weights as a function of time - Male',
#      xlab='Days since start of follow-up', ylab='IPTW (log-scale)', lwd=2, ylim = c(yrange[1], yrange[2]), type = 'n', 
#      yaxt = 'n'
#      )
# lines(swintervalsumlog$median, lty=1, col='black', lwd=2)
# lines(swintervalsumlog$minsw, lty=2, col='red', lwd=2)
# lines(swintervalsumlog$maxsw, lty=2, col='red', lwd=2)
# lines(swintervalsumlog$q1, lty=3, col='blue', lwd=2)
# lines(swintervalsumlog$q3, lty=3, col='blue', lwd=2)
# lines(swintervalsumlog$q10, lty=4, col='purple', lwd=2)
# lines(swintervalsumlog$q90, lty=4, col='purple', lwd=2)
# axis(side = 2, at = seq(yrange[1], yrange[2], length.out=5), labels = round(exp(seq(yrange[1], yrange[2], length.out=5)),1))
# legend('topright', c('Median', '50% Interval', '90% Interval', 'Min./Max.'), cex = 0.8, 
#        col=c('black', 'blue', 'purple', 'red'), lty=c(1,3,4,2), lwd=2)
# dev.off()
# 
#### graph
# plot(x=c(0,cutoff), y=swintervalsum$median, lty='longdash', main = 'Estimated stabilized weights as a function of time',
#      xlab='Days since start of follow-up', ylab='IPTW', ylim = c(0,5.5), lwd=2, type = 'n')
# lines(x=c(0,cutoff), swintervalsum$median, lty=1, col='black', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$minsw, lty=2, col='red', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$maxsw, lty=2, col='red', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$q1, lty=3, col='blue', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$q3, lty=3, col='blue', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$q10, lty=4, col='purple', lwd=2)
# lines(x=c(0,cutoff), swintervalsum$q90, lty=4, col='purple', lwd=2)





###################### Patients died without outcome ##
# #nrow(gc_init[gc_init$yt == 365 & gc_init$status == 0, ])
# bp_init <- bp_final %>%
#   group_by(study_id) %>%
#   filter(row_number() == n())
# 
# bpstats <- left_join(bp_init, base_sex)
# table(bpstats$sex)
# 
# 
# ltcstats <- left_join(ltc_conf, base_sex)
# ltcstats <- ltcstats[is.finite(ltcstats$days_to_ltcdate),]
# table(ltcstats$sex)

########## Useless #########
# Diagnostics
# Parameters to specify before source:
# stepsize = 2

##### Baseline variables association with outcome #####
#library(nnet)

# gc <- read.csv("gc_final.csv", sep= ",", header= T)
# df_Dinit <- read.csv("df_Dinit.csv", sep= ",", header= T)
# long_Y_cbg <- long_Y_cb     # continue from Outcome_models.R
# #long_Y_cb <- read.csv("long_Y_cb.csv", sep= ",", header= T)
# female = FALSE
# 
# # for male              
# #long_Y_cbg <- long_Y_cb[long_Y_cb$sex == 1, ]  # n=_____, actual sampled in outcome: _______
# summary(long_Y_cbg$futime) # Get the average cut-off and compare with dayssupl for GC and BP to determine 
#                            # the underestimation due to discretization -> mean cut: 45 days
# 
# df0 <- df_Dinit[df_Dinit$study_id %in% unique(long_Y_cbg$study_id), ]


############################ Prescreening #################
# gc_yt <- gc[ , c('study_id', 'yt')]
# basecheck <- left_join(df_Dinit, gc_yt, by='study_id')
# # "diabetes", "arthritis",  "pfx",     "antisexhor",   "inhalegc",  "bronchodilators", 
# # "sexhor",   "thiazide",   "pralcal", "pdmab",        "falldr",    "fallcomorb", 
# # "emerhosp", "totprevBPdaysc"
# # Stratify by Sex
# basecheckf <- basecheck[basecheck$sex == 0, ]
# basecheckm <- basecheck[basecheck$sex == 1, ]
# 
# # With outcome
# survdiff(Surv(yt, status) ~ as.numeric(age), data = basecheckf)
# survdiff(Surv(yt, status) ~ as.numeric(age), data = basecheckm)
# 
# # With initial dose
# chisq.test(basecheckm$bronchodilators,basecheckm$D0)


#################### Density Rugs plot ####################
# par(mfrow=c(2,2))
# longycb.densw <- long_Y_cb %>%
#   dplyr::select(sw, d0p, vp, dp)
# 
# swdf <- data.frame(weight = longycb.densw$sw, weighttype = rep('Combined sw', nrow(longycb.densw)))
# d0pdf <- data.frame(weight = longycb.densw$d0p, weighttype = rep('Initial dose sw', nrow(longycb.densw)))
# vpdf <- data.frame(weight = longycb.densw$vp, weighttype = rep('Subsequent visiting sw', nrow(longycb.densw)))
# dpdf <- data.frame(weight = longycb.densw$dp, weighttype = rep('Subsequent dose sw', nrow(longycb.densw)))
# 
# superlongswdenscompare <- rbind(swdf, d0pdf, vpdf, dpdf)

# par(mfrow=c(1,1))
# sm.density.compare(superlongswdenscompare$weight, superlongswdenscompare$weighttype, 
#                    xlab = 'Treatment duration/ days-supplied values', xlim = c(0,6), lty=11:14, 
#                    col=c('red', 'dark gray', 'deeppink4', 'darkgreen'), lwd=2, h=0.3)
# title(main='Density comparison of different weighting methods')
# rug(long_Y_cb$sw)


# pdf('/users/edong/Output/densitymale.pdf', width=9.5, height=8)
## Different scale 
# par(mfrow=c(2,2))
# plot(density(long_Y_cb$sw), main = paste('Density plot of combined stabilized weight (sw) \n for exposure A -', gender.name), 
#      xlab="Inverse Probability Weights", col='darkred')
# rug(long_Y_cb$sw, col='darkred')
# 
# plot(density(long_Y_cb$d0p), main = paste('Density plot of initial dose sw \n for exposure D0 -', gender.name), 
#      xlab="Inverse Probability Weights")
# rug(long_Y_cb$d0p)
# 
# plot(density(long_Y_cb$vp), main = paste('Density plot of subsequent visiting sw \n for exposure V -', gender.name), 
#      xlab="Inverse Probability Weights")
# rug(long_Y_cb$vp)
# 
# plot(density(long_Y_cb$dp), main = paste('Density plot of subsequent dose sw \n for exposure D -', gender.name), 
#      xlab="Inverse Probability Weights")
# rug(long_Y_cb$dp)
#dev.off()

## All same scale
#par(mfrow=c(2,2))
#xrange = range(long_Y_cb$sw)
