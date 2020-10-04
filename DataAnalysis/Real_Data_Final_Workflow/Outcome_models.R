rangecumd <- range(long_Y_cb$gccumdosescale)
# Fit cumulative dose as a spline
long_Y_cb$timebasis_cumdosescale <- bs(long_Y_cb$gccumdosescale, 
                                       Boundary.knots = c(rangecumd[1], rangecumd[2]), degree = 2)
if(female){
  ################ A - Female
  # Combined weights
  cbout.sw.slm.a <- glm(Yevent ~ timebasis_epi + A + offset(log(futime)) 
                        + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                        + fallcomorb + emerhosp + pfx
                        + pralcal + diabetes, family = binomial(link=logit),
                        weights = sw.trunc, 
                        data = long_Y_cb)
  # Unweighted
  cbout.un.slm.a <- glm(Yevent ~ timebasis_epi + A + offset(log(futime)) 
                        + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                        + fallcomorb + emerhosp + pfx
                        + pralcal + diabetes, family = binomial(link=logit),
                        data = long_Y_cb)
  # Weighted null model for likelihood ratio test
  cbout.sw.null.a <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                         + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                         + fallcomorb + emerhosp + pfx
                         + pralcal + diabetes, family = binomial(link=logit),
                         weights = sw.trunc, 
                         data = long_Y_cb)
  ################ Continuous cumdose (centered) - Female
  # Combined weights
  cbout.sw.slm.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                           + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                           + fallcomorb + emerhosp + pfx + gccumdosescale
                           + pralcal + diabetes, family = binomial(link=logit),
                           weights = sw.trunc, 
                           data = long_Y_cb)
  # Unweighted
  cbout.un.slm.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                           + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                           + fallcomorb + emerhosp + pfx + gccumdosescale
                           + pralcal + diabetes, family = binomial(link=logit),
                           data = long_Y_cb)
  # Weighted null model for likelihood ratio test
  cbout.sw.null.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                            + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                            + fallcomorb + emerhosp + pfx
                            + pralcal + diabetes, family = binomial(link=logit),
                            weights = sw.trunc, 
                            data = long_Y_cb)
  ################ Spline cumdose (centered) - Female
  # Combined weights
  cbout.sw.slm.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                             + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                             + fallcomorb + emerhosp + pfx + timebasis_cumdosescale
                             + pralcal + diabetes, family = binomial(link=logit),
                             weights = sw.trunc, 
                             data = long_Y_cb)
  # Unweighted
  cbout.un.slm.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                             + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                             + fallcomorb + emerhosp + pfx + timebasis_cumdosescale
                             # Female additional
                             + pralcal + diabetes, family = binomial(link=logit),
                             data = long_Y_cb)
  # Weighted null model for likelihood ratio test
  cbout.sw.null.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                              + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                              + fallcomorb + emerhosp + pfx 
                              + pralcal + diabetes, family = binomial(link=logit),
                              weights = sw.trunc, 
                              data = long_Y_cb)
}else{
  ############### A - Male
  # Combined
  cbout.sw.slm.a <- glm(Yevent ~ timebasis_epi + A + offset(log(futime)) 
                        + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                        + fallcomorb + emerhosp + pfx, family = binomial(link=logit), 
                        weights = sw.trunc,
                        data = long_Y_cb)
  # Unweighted
  cbout.un.slm.a <- glm(Yevent ~ timebasis_epi + A + offset(log(futime)) 
                        + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                        + fallcomorb + emerhosp + pfx, family = binomial(link=logit),
                        data = long_Y_cb)
  # Weighted null model for likelihood ratio test
  cbout.sw.null.a <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                         + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                         + fallcomorb + emerhosp + pfx, family = binomial(link=logit), 
                         weights = sw.trunc,
                         data = long_Y_cb)
  
  ############### Continuous cumdose (centered) - Male
  # Combined
  cbout.sw.slm.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                           + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                           + fallcomorb + emerhosp + pfx + gccumdosescale, family = binomial(link=logit), 
                           weights = sw.trunc,
                           data = long_Y_cb)
  # Unweighted
  cbout.un.slm.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                           + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                           + fallcomorb + emerhosp + pfx + gccumdosescale, family = binomial(link=logit),
                           data = long_Y_cb)
  
  # Weighted null model for likelihood ratio test
  cbout.sw.null.cumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                            + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                            + fallcomorb + emerhosp + pfx, family = binomial(link=logit), 
                            weights = sw.trunc,
                            data = long_Y_cb)
  
  ############### Spline cumdose (centered) - Male
  # Combined
  cbout.sw.slm.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                             + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                             + fallcomorb + emerhosp + pfx + timebasis_cumdosescale, family = binomial(link=logit), 
                             weights = sw.trunc,
                             data = long_Y_cb)
  # Unweighted
  cbout.un.slm.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                             + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                             + fallcomorb + emerhosp + pfx + timebasis_cumdosescale, family = binomial(link=logit),
                             data = long_Y_cb)
  
  # Weighted null model for likelihood ratio test
  cbout.sw.null.spcumd <- glm(Yevent ~ timebasis_epi + offset(log(futime)) 
                              + agescale + falldr + sexhor + totprevBPdaysc + inhalegc + thiazide + arthritis
                              + fallcomorb + emerhosp + pfx, family = binomial(link=logit), 
                              weights = sw.trunc,
                              data = long_Y_cb)
  
}
