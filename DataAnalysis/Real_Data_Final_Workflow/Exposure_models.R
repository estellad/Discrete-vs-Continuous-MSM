################################## Initial D model: D0_gender ###########################
d0modm <- polr(D0 ~  
                 agescale + arthritis + pfx #+ antisexhor + pralcal + pdmab + bronchodilators
               + inhalegc + sexhor + thiazide + falldr + fallcomorb
               + emerhosp + totprevBPdaysc, data=D0_gender, 
               control = list(maxit = 2000, trace = 1, REPORT=1))
summary(d0modm)

## Male Conditional:
d0modc <- polr(D0 ~ #diabetes 
                 agescale + arthritis + pfx #+ antisexhor + pralcal + pdmab + bronchodilators
               + inhalegc + sexhor + thiazide + falldr + fallcomorb
               + emerhosp + totprevBPdaysc, data=D0_gender, 
               control = list(maxit = 2000, trace = 1, REPORT=1))
summary(d0modc)

######################################### V model: V_gender ############################
## Male Marginal
V_gender$timebasis_epi <- bs(V_gender$cycle, Boundary.knots=c(0, 73), degree = 2) 
vmodm <- glm(Visit ~  Dtrtp #+ offset(log(futime)) + diabetes 
             + timebasis_epi + gccumdosescale #+ as.numeric(cycle) 
             #+ Ltc + cumBPdaysp 
             + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
             + inhalegc + sexhor + thiazide + falldr + fallcomorb 
             + emerhosp + totprevBPdaysc,
             family = binomial(link=logit), data=V_gender)
summary(vmodm)

# Male Conditional
vmodc <- glm(Visit ~ Dtrtp #+ diabetes + offset(log(futime)) 
             + timebasis_epi + gccumdosescale #+ as.numeric(cycle)
             + Ltc + cumBPdayspscale 
             + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
             + inhalegc + sexhor + thiazide + falldr + fallcomorb 
             + emerhosp + totprevBPdaysc,
             family = binomial(link=logit), data=V_gender)
summary(vmodc)



######################## Subset D_gender data ###########################
Dtrtp0dat = D_gender[D_gender$Dtrtp == 0, ] # Dtrtp=0; Dtrtp=1,2,3,4,5 (0 is impossible)
Dtrtp0dat$Dtrtp <- droplevels(factor(Dtrtp0dat$Dtrtp))
Dtrtp0dat$Dtrt <- droplevels(factor(Dtrtp0dat$Dtrt))
Dtrtp0dat$Dtrt <- relevel(Dtrtp0dat$Dtrt, ref="1")

Dtrtpnot0dat = D_gender[D_gender$Dtrtp != 0, ] # Dtrtp=1,2,3,4,5; Dtrt=0,1,2,3,4,5
Dtrtpnot0dat$Dtrtp <- droplevels(factor(Dtrtpnot0dat$Dtrtp))
Dtrtpnot0dat$Dtrtp <- relevel(Dtrtpnot0dat$Dtrtp, ref="1")
Dtrtpnot0dat$Dtrt <- relevel(Dtrtpnot0dat$Dtrt, ref="0")

########################### D model: Only those with Dtrtp = 0 ############################
# Male Marginal
Dtrtp0dat$timebasis_epi <- bs(Dtrtp0dat$cycle, Boundary.knots=c(0, 73), degree = 2)
dmodp0m <- polr(Dtrt ~ #Dtrtp #+ diabetes 
                  + timebasis_epi + gccumdosescale #+ Ltc + cumBPdaysp                      
                + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
                + inhalegc + sexhor + thiazide + falldr + fallcomorb 
                + emerhosp + totprevBPdaysc, data=Dtrtp0dat, 
                control = list(maxit = 2000, trace = 1, REPORT=1))
summary(dmodp0m)                  

# Male Conditional
dmodp0c <- polr(Dtrt ~ #Dtrtp #+ diabetes 
                  + timebasis_epi + Ltc + cumBPdayspscale + gccumdosescale 
                + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
                + inhalegc + sexhor + thiazide + falldr + fallcomorb 
                + emerhosp + totprevBPdaysc, data=Dtrtp0dat, 
                control = list(maxit = 2000, trace = 1, REPORT=1))
summary(dmodp0c)

##################### D model: Only those with Dtrtp > 0 ############################
# Male Marginal
Dtrtpnot0dat$timebasis_epi <- bs(Dtrtpnot0dat$cycle, Boundary.knots=c(0, 73), degree = 2)
dmodpnot0m <- polr(Dtrt ~ Dtrtp #+ diabetes 
                   + timebasis_epi + gccumdosescale + Ltc + cumBPdaysp                      
                   + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
                   + inhalegc + sexhor + thiazide + falldr + fallcomorb 
                   + emerhosp + totprevBPdaysc, data=Dtrtpnot0dat, 
                   control = list(maxit = 2000, trace = 1, REPORT=1))                
summary(dmodpnot0m)

# Male Conditional
dmodpnot0c <- polr(Dtrt ~ Dtrtp #+ diabetes 
                   + timebasis_epi + Ltc + cumBPdayspscale + gccumdosescale 
                   + agescale + arthritis + pfx #+ pralcal + antisexhor + pdmab + bronchodilators
                   + inhalegc + sexhor + thiazide + falldr + fallcomorb 
                   + emerhosp + totprevBPdaysc, data=Dtrtpnot0dat, 
                   control = list(maxit = 2000, trace = 1, REPORT=1))
summary(dmodpnot0c)




