library(UsingR)
library(sm)  # Only in R studio environment

female = FALSE
base_sex <- df_Dinit[, c('ikn', 'sex')]

gc_init <- gc %>%
  group_by(study_id) %>%
  filter(row_number() == 1) %>%
  full_join(base_sex)

deathyt <- gc_init[gc_init$status == 0 & gc_init$yt<365, c('yt', 'sex')] # use color 6
deathyt <- deathyt %>%
  group_by(sex) %>%
  arrange(sex)
deathyt$sex <- factor(deathyt$sex, labels = c('Female', 'Male'))

densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}

deathyt <- ddply(deathyt, 'sex', transform, 
                 val_peak = signif(densMode(yt)$x, 3)
)


fxyt <- gc_init[gc_init$status == 1, c('study_id', 'yt', 'sex')] # use color 10
fxyt <-fxyt %>%
  group_by(sex) %>%
  arrange(sex, study_id) 
#fxyt_male <- fxyt[fxyt$sex == 1, ]
fxyt$sex <- factor(fxyt$sex, labels = c('Female', 'Male'))
fxyt <- ddply(fxyt, 'sex', transform, 
              val_peak = signif(densMode(yt)$x, 3)
)

############################ Density comparison on death and fracture by sex #################
pdf('/users/edong/Output/Descriptive/deathfxdays.pdf', width=9.5, height=4)
par(mfrow=c(1,2), xaxs = 'i')
sm.density.compare(deathyt$yt, deathyt$sex, xlab = 'Death date', xlim = c(0,365), lty=1:2, col=c(6,6), 
                   lwd=2, ylim=c(0,0.005))
title(main='Density Comparison of Death Time by Sex')
L=legend(210, 0.005, legend=c('Female', 'Male'), col=c(6,6), lty=1:2, ncol=1, bty='n', lwd=2, cex=0.8,
         x.intersp=0.5, inset=0.02, title=NA)
legend(x=210, y=0.005, legend=c('____', '____'), cex=0.8,
       col=rep(NA,2), lty=c(1,1), ncol=1, x.intersp=5, bg=NA, title="Sex/ Num. Death")
#axis(side = 1, at = seq(0, 365, by=30))
abline(v=deathyt[1, 'val_peak'], col = 6, lty=1)
abline(v=deathyt[nrow(deathyt),'val_peak'], col = 6, lty = 2)


sm.density.compare(fxyt$yt, fxyt$sex, xlab = 'Fracture date', xlim = c(0,365), lty=1:2, col=c('purple','purple'), 
                   lwd=2, ylim=c(0,0.005))
title(main='Density Comparison of Fracture Time by Sex')
# legend('topright', legend=c('Female', 'Male'), cex = 0.8, col=c(10,10), lty=1:2, title="Sex")
L=legend(210, 0.005, legend=c('Female', 'Male'), col=c('purple','purple'), lty=1:2, ncol=1, bty='n', lwd=2, cex=0.8,
         x.intersp=0.5, inset=0.02, title=NA)
legend(x=210, y=0.005, legend=c('___', '   ____'), cex=0.8,
       col=rep(NA,2), lty=c(1,1), ncol=1, x.intersp=5, bg=NA, title="Sex/ Num. Fracture")
#axis(side = 1, at = seq(0, 365, by=30))
abline(v=fxyt[1, 'val_peak'], col = 'purple', lty=1)
abline(v=fxyt[nrow(fxyt),'val_peak'], col = 'purple', lty = 2)

dev.off()



# hist(nonadmcensyt, breaks=100, xlab='days', xlim=c(0,365), 
#      main='Event time for patients without')
# 
# hist(deathyt, breaks=100, xlab='days', xlim=c(0,365), 
#      main='Event time for patients without')
# 
# hist(fxyt, breaks=100, xlab='days', xlim=c(0,365), 
#      main='Event time for patients without')
# 
# 
# # nonadmdf <- data.frame(
# #   day = nonadmcensyt, 
# #   group=rep("Non-admin. censor", length(nonadmcensyt)))


# mode(density(deathyt))
# 
# library(sm)
# sm.density.compare(dftot$day, dftot$group, xlab = 'days', xlim = c(-50,380), lty=c(2,1), col=c(6,10))
# title(main='Event time distributions for patients not administratively censored')
# legend('topright', levels(dftot$group), cex = 0.8, col=c(6,10), lty=c(2,1), 
#        title="Event Types")
# 
# getmode <- function(v){
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v,uniqv)))]
# }


# deathytdf <- data.frame(
#   day = deathyt$yt, 
#   group=rep("Death", nrow(deathyt)))
# fxytdf <- data.frame(
#   day = fxyt$yt, 
#   group=rep("Fracture", nrow(fxyt)))
# 
# dftot <- rbind(deathytdf, fxytdf)
# Superimpose a density on histogram
# par(mfrow=c(1,2))
# hist(deathyt, freq=FALSE, breaks = 25, main='Death time', xlab='days', xlim=c(0,365), xaxt='n')
# #lines(density(deathyt), col=6, lwd=2) # density estimates with default
# lines(density(deathyt, adjust=2), col=6, lwd=2, lty=2) # add another smoother density
# axis(side = 1, at = seq(0,360, by=20))
# 
# 
# hist(fxyt, freq=FALSE, breaks = 25, main='Fracture time', xlab='days', xlim=c(0,365), xaxt='n')
# #lines(density(fxyt), col=10, lwd=2, lty=1) # density estimates with default
# lines(density(fxyt, adjust=2), col=10, lwd=2, lty=2) # add another smoother density
# #legend('topright', c('density', 'smoother density'), cex = 0.8, lty = c(1,2), col=10)
# axis(side = 1, at = seq(0,360, by=20))



table(gc$Dtrt)
table(long_V_slm$Dtrt)

gc$Start = ifelse(gc$Start %% 1 != 0, gc$Start + 0.1, gc$Start)
gc$End = ifelse(gc$End %% 1 != 0, gc$End + 0.1, gc$End)

long.gc.by.dtrt <- gc %>%
  mutate(dayssupl = End-Start) %>%
  select(Dtrt, dayssupl) %>%
  arrange(Dtrt) %>%
  mutate(Dtrt = as.factor(Dtrt))

# # table won't use, will choose density plot compare
# dtrt.dur.summary <- long.gc.by.dtrt %>%
#   group_by(Dtrt) %>%
#   summarise(
#     minduration = min(dayssupl),
#     maxduration = max(dayssupl),
#     meantrtduration = mean(dayssupl),
#     sdduration = sd(dayssupl),
#     q1 = quantile(dayssupl, probs = 0.25),
#     median = quantile(dayssupl, probs = 0.50),
#     q3 = quantile(dayssupl, probs = 0.75),
#     n=n()) 



pdf('/users/edong/Output/Descriptive/density_compare_doses.pdf', width=9, height=6)
par(mfrow=c(1,1))
sm.density.compare(long.gc.by.dtrt$dayssupl, long.gc.by.dtrt$Dtrt, xlab = 'Treatment duration/ days-supplied values', 
                   xlim = c(-10,110), lty=1:5, col=1:5, lwd=2)
title(main='Density Comparison of Treatment Duration by Daily Dose Categories')

# legend('topright', legend=c("0<daily dose <= 5mg: n=____ ", 
#                             "5 < daily dose <= 10mg: n=____", 
#                             "10 < daily dose < 30mg: n=____",
#                             "30 <= daily dose < 50mg: n=____", 
#                             "50mg <= daily dose: n=____"), 
#        cex = 0.8, col=1:5, lty=1:5, title="Daily dose level/ Num. dispensations")

L=legend(64, 0.125, legend=c("0<daily dose <= 5mg", 
                             "5 < daily dose <= 10mg", 
                             "10 < daily dose < 30mg",
                             "30 <= daily dose < 50mg", 
                             "50mg <= daily dose"), col=1:5, lty=1:5, ncol=1, bty='n', lwd=2, cex=0.9,
         x.intersp=0.5, inset=0.02, title=NA)
legend(x=65, y=0.125, legend=c('____', '____', '____', '____', '____'), cex=0.9,
       col=rep(NA,5), lty=c(1,1), ncol=1, x.intersp=13, bg=NA, title="Daily dose level/ Num. dispensations")
dev.off()
