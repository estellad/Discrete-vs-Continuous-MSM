################################## Exposure models ##################################
if(expo & !outcome){ 
  get_multinom_exposure_tables <- function(model, numintercept, order, nota, label){
    ctable <- coef(summary(model))
    ctable <- ctable[order, ]
    p <- format(round(pnorm(abs(ctable[,'t value']), lower.tail = F)*2,3), nsmall=2)
    pval <- ifelse(p<0.001, '<0.001 *', 
                   ifelse(p<0.05, paste0(as.character(p), ' *'), paste0(as.character(p), '  ')))
    ctablesub <- ctable[1:(nrow(ctable)-numintercept),]
    ci2 <- cbind(ctablesub[,1]-1.96*ctablesub[,2], ctablesub[,1]+1.96*ctablesub[,2])
    or <- coef(model)[order][1:(length(order) - numintercept)]
    ors <- format(round(exp(cbind(OR=or, ci2)), digits = 2), nsmall=2)
    orsdone <- cbind(ors[,1], paste0('(', ors[,2],', ', ors[,3], ')'))
    orsdone <- rbind(orsdone, matrix("  ", nrow=numintercept, ncol=2, byrow=T))
    
    # OR table 
    resulttab <- as.data.frame(cbind(nota, format(round(ctable[,-3],2), nsmall=2), orsdone, pval))
    colnames(resulttab) <- c('Notation', 'Coef. Estimate', 'Std.Error', 'OR', '95% CI', 'p-value')
    rownames(resulttab) <- label; resulttab
    
    return(list(ors, resulttab))
  }
  
  get_V_exposure_table <- function(vmodc, vorder, vnota, vlabel){
    ctable <- coef(summary(vmodc))
    ctable <- ctable[vorder, ]
    p <- format(round(ctable[,4], 3), nsmall=2)
    pval <- ifelse(p<0.001, '<0.001 *', 
                   ifelse(p<0.05, paste0(as.character(p), ' *'), paste0(as.character(p), '  ')))
    ci2 <- cbind(ctable[,1]-1.96*ctable[,2], ctable[,1]+1.96*ctable[,2])
    hr <- coef(vmodc)[vorder]
    hrsv <- format(round(exp(cbind(HR=hr, ci2)), digits = 2), nsmall=2)
    hrsdone <- cbind(hrsv[,1], paste0('(', hrsv[,2],', ', hrsv[,3], ')'))
    
    # HR table 
    resultvtab <- as.data.frame(cbind(vnota, format(round(ctable[,-c(3,4)],2), nsmall=2), hrsdone, pval))
    colnames(resultvtab) <- c('Notation', 'Coef. Estimate', 'Std.Error', 'HR', '95% CI', 'p-value')
    rownames(resultvtab) <- vlabel; resultvtab
    
    return(list(hrsv, resultvtab))
  }
  
  # ########################### D0 conditional
  d0tables <- get_multinom_exposure_tables(d0modc, numintercept=4, d0order, d0nota, d0label)
  orsd0 <- d0tables[[1]]
  resultd0tab <- d0tables[[2]]
  
  # ############################# V conditional
  vtables <- get_V_exposure_table(vmodc, vorder, vnota, vlabel)
  hrsv <- vtables[[1]]
  resultvtab <- vtables[[2]]
  
  # ############################### D conditional -> Dtrtp=0
  dp0tables <- get_multinom_exposure_tables(dmodp0c, numintercept=4, dp0order, dp0nota, dp0label)
  orsdp0 <- dp0tables[[1]]
  resultdp0tab <- dp0tables[[2]]
  
  # ############################### D conditional -> Dtrtp>0
  dpnot0tables <- get_multinom_exposure_tables(dmodpnot0c, numintercept=5, dpnot0order, dpnot0nota, dpnot0label)
  orsdpnot0 <- dpnot0tables[[1]]
  resultdpnot0tab <- dpnot0tables[[2]]
  
  
  ## Output D0 latex and csv table
  if(female){
    print(xtable(resultd0tab, type = 'latex'), file = '/users/edong/Output/Female/Exposure/femaled0.tex', sanitize.text.function = function(x){x})
    write.csv(resultd0tab, '/users/edong/Output/Female/Exposure/femaled0.csv')
  }else{
    print(xtable(resultd0tab, type = 'latex'), file = '/users/edong/Output/Male/Exposure/maled0.tex', sanitize.text.function = function(x){x})
    write.csv(resultd0tab, '/users/edong/Output/Male/Exposure/maled0.csv')
  }
  
  ## Output V latex and csv table
  if(female){
    print(xtable(resultvtab, type = 'latex'), file = '/users/edong/Output/Female/Exposure/femalev.tex', sanitize.text.function = function(x){x})
    write.csv(resultvtab, '/users/edong/Output/Female/Exposure/femalev.csv')
  }else{
    print(xtable(resultvtab, type = 'latex'), file = '/users/edong/Output/Male/Exposure/malev.tex', sanitize.text.function = function(x){x})
    write.csv(resultvtab, '/users/edong/Output/Male/Exposure/malev.csv')
  }
  
  ## Output D latex and csv table, previous dose = 0
  if(female){
    print(xtable(resultdp0tab, type = 'latex'), file = '/users/edong/Output/Female/Exposure/femaledp0.tex', sanitize.text.function = function(x){x})
    write.csv(resultdp0tab, '/users/edong/Output/Female/Exposure/femaledp0.csv')
  }else{
    print(xtable(resultdp0tab, type = 'latex'), file = '/users/edong/Output/Male/Exposure/maledp0.tex', sanitize.text.function = function(x){x})
    write.csv(resultdp0tab, '/users/edong/Output/Male/Exposure/maledp0.csv')
  }
  
  ## Output D latex and csv table, previous dose > 0
  if(female){
    print(xtable(resultdpnot0tab, type = 'latex'), file = '/users/edong/Output/Female/Exposure/femaledpnot0.tex', sanitize.text.function = function(x){x})
    write.csv(resultdpnot0tab, '/users/edong/Output/Female/Exposure/femaledpnot0.csv')
  }else{
    print(xtable(resultdpnot0tab, type = 'latex'), file = '/users/edong/Output/Male/Exposure/maledpnot0.tex', sanitize.text.function = function(x){x})
    write.csv(resultdpnot0tab, '/users/edong/Output/Male/Exposure/maledpnot0.csv')
  }
  
  
  ################################## Forest plot ###################################
  plot_forest <- function(resulttab, ratios, riskmeasure, numintercept, modname, female, titlename, width, height){
    rownames <- rownames(resulttab) 
    rownames <- gsub('$', '', rownames, fixed = TRUE)
    rownames <- gsub('\\geq', '>=', rownames, fixed = TRUE)
    rownames <- gsub('d_0', 'd0', rownames, fixed = TRUE)
    rownames <- gsub('\\leq', '<=', rownames, fixed = TRUE)
    cleanedrownames <- gsub('\\mid', '|', rownames, fixed = TRUE)
    
    if(modname == 'V'){
      ratiolist <- c(riskmeasure, format(as.numeric(ratios[, 1]), nsmall = 2))
    }else{
      ratiolist <- c(riskmeasure, format(as.numeric(ratios[, 1]), nsmall = 2), rep(NA, numintercept))
    }
    
    tabletext <- cbind(c('Covariate', cleanedrownames), 
                       c('p-value', as.character(resulttab[,6])),
                       ratiolist
    )
    
    if(modname == 'V'){
      m <- c(NA, as.numeric(ratios[, 1]))
      l <- c(NA, as.numeric(ratios[, 2]))
      u <- c(NA, as.numeric(ratios[, 3]))
    }else{
      m <- c(NA, as.numeric(ratios[, 1]), rep(NA, numintercept))
      l <- c(NA, as.numeric(ratios[, 2]), rep(NA, numintercept))
      u <- c(NA, as.numeric(ratios[, 3]), rep(NA, numintercept))
    }
    
    if(female){
      pdf(paste0('/users/edong/Output/Female/Exposure/', modname, 'modelforest.pdf'), width = width, height = height, onefile = FALSE)
    }else{
      pdf(paste0('/users/edong/Output/Male/Exposure/', modname, 'modelforest.pdf'), width = width, height = height, onefile = FALSE)
    }
    
    par(mar = c(2,3,2,2))
    forestplot(tabletext, m, l, u, zero=1, 
               is.summary=c(TRUE, rep(FALSE, nrow(resulttab))),
               xticks = seq(range(as.numeric(ratios[, 1]))[1], 
                            range(as.numeric(ratios[, 1]))[2], length.out = 5), 
               xlog = TRUE, boxsize = 0.2,
               col = fpColors(box = 'royalblue', line = 'darkblue'),
               title = titlename
    )
    dev.off()
  }
  
  
  plot_forest(resultd0tab, orsd0, 'OR', 4, 'D0', female, paste0('Initial Dose Model - ', gender.name), width=8, height=8)
  plot_forest(resultvtab, hrsv, 'HR', NA, 'V', female, paste0('Subsequent Visiting Model - ', gender.name), width=8, height=10)
  plot_forest(resultdp0tab, orsdp0, 'OR', 4, 'Dp=0', female, paste0('Subsequent Dosage Model For Previous Dose Is 0 - ', gender.name), width=8, height=10)
  plot_forest(resultdpnot0tab, orsdpnot0, 'OR', 5, 'Dpnot0', female, paste0('Subsequent Dosage Model For Previous Dose Is Not 0 - ', gender.name), width=8, height=12)
  
  
  ################## Check Cross Tabulation of Large ORs in the Forest Plot ###################
  crosstab_Dtrtpnot0dat <- table(Dtrtpnot0dat$Dtrt, Dtrtpnot0dat$Dtrtp) # subseted in Exposure_model.R
  if(female){
    write.csv(crosstab_Dtrtpnot0dat, '/users/edong/Output/Female/Exposure/crosstabfemale_Dtrtpnot0dat.csv')
  }else{
    write.csv(crosstab_Dtrtpnot0dat, '/users/edong/Output/Male/Exposure/crosstabmale_Dtrtpnot0dat.csv')
  }
  
  
}else if(!expo & outcome){
  ################ Outcome Model output cleaning ##############
  ## Combined weights table
  ResultSwTab <- function(sw.model, order, nota, label){
    ctable <- coef(summary(sw.model))
    ctable <- ctable[order, ]
    p <- format(round(ctable[,4], 3), nsmall=2)
    pval <- ifelse(p<0.001, '<0.001 *', 
                   ifelse(p<0.05, paste0(as.character(p), ' *'), paste0(as.character(p), '  ')))
    ci2 <- cbind(ctable[,1]-1.96*ctable[,2], ctable[,1]+1.96*ctable[,2])
    hr <- coef(sw.model)[order]
    hrs <- format(round(exp(cbind(HR=hr, ci2)), digits = 2), nsmall=2)
    hrsdone <- cbind(hrs[,1], paste0('(', hrs[,2],', ', hrs[,3], ')'))
    
    # HR table 
    resultswtab <- as.data.frame(cbind(nota, format(round(ctable[,-c(3,4)],2), nsmall=2), hrsdone, pval))
    colnames(resultswtab) <- c('Notation', 'Coef. Estimate', 'Std.Error', 'HR', '95% CI', 'p-value')
    rownames(resultswtab) <- label; resultswtab
    
    return(resultswtab)
  }
  
  
  ## Unweighted
  ResultUnTab <- function(un.model, order, nota, label){
    ctable <- coef(summary(un.model))
    ctable <- ctable[order, ]
    p <- format(round(ctable[,4], 3), nsmall=2)
    pval <- ifelse(p<0.001, '<0.001 *', 
                   ifelse(p<0.05, paste0(as.character(p), ' *'), paste0(as.character(p), '  ')))
    ci2 <- cbind(ctable[,1]-1.96*ctable[,2], ctable[,1]+1.96*ctable[,2])
    hr <- coef(un.model)[order]
    hrs <- format(round(exp(cbind(HR=hr, ci2)), digits = 2), nsmall=2)
    hrsdone <- cbind(hrs[,1], paste0('(', hrs[,2],', ', hrs[,3], ')'))
    
    # HR table 
    resultuntab <- as.data.frame(cbind(nota, format(round(ctable[,-c(3,4)],2), nsmall=2), hrsdone, pval))
    colnames(resultuntab) <- c('Notation', 'Coef. Estimate', 'Std.Error', 'HR', '95% CI', 'p-value')
    rownames(resultuntab) <- label; resultuntab
    
    return(resultuntab)
  }
  
  # Amodels
  resultswtab.a <- ResultSwTab(cbout.sw.slm.a, aorder, anota, alabel)
  resultuntab.a <- ResultUnTab(cbout.un.slm.a, aorder, anota, alabel)
  # cumd models
  resultswtab.cumd <- ResultSwTab(cbout.sw.slm.cumd, cumdorder, cumdnota, cumdlabel); 
  resultuntab.cumd <- ResultUnTab(cbout.un.slm.cumd, cumdorder, cumdnota, cumdlabel)
  # spcumd models
  resultswtab.spcumd <- ResultSwTab(cbout.sw.slm.spcumd, spcumdorder, spcumdnota, spcumdlabel); 
  resultuntab.spcumd <- ResultUnTab(cbout.un.slm.spcumd, spcumdorder, spcumdnota, spcumdlabel)
  
  
  ## Output Amodel weighted & unweighted latex table
  if(female){
    #SW
    print(xtable(resultswtab.a, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleswa.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.a, '/users/edong/Output/Female/Outcome/femaleswa.csv')
    #Un
    print(xtable(resultuntab.a, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleuna.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.a, '/users/edong/Output/Female/Outcome/femaleuna.csv')
  }else{
    print(xtable(resultswtab.a, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleswa.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.a, '/users/edong/Output/Male/Outcome/maleswa.csv')
    print(xtable(resultuntab.a, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleuna.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.a, '/users/edong/Output/Male/Outcome/maleuna.csv')
  }
  
  ## Output cumd model latex table
  if(female){
    print(xtable(resultswtab.cumd, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleswcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.cumd, '/users/edong/Output/Female/Outcome/femaleswcumd.csv')
    print(xtable(resultuntab.cumd, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleuncumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.cumd, '/users/edong/Output/Female/Outcome/femaleuncumd.csv')
  }else{
    print(xtable(resultswtab.cumd, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleswcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.cumd, '/users/edong/Output/Male/Outcome/maleswcumd.csv')
    print(xtable(resultuntab.cumd, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleuncumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.cumd, '/users/edong/Output/Male/Outcome/maleuncumd.csv')
  }
  
  ## Output spcumd model latex table
  if(female){
    print(xtable(resultswtab.spcumd, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleswspcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.spcumd, '/users/edong/Output/Female/Outcome/femaleswspcumd.csv')
    print(xtable(resultuntab.spcumd, type = 'latex'), file = '/users/edong/Output/Female/Outcome/femaleunspcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.spcumd, '/users/edong/Output/Female/Outcome/femaleunspcumd.csv')
  }else{
    print(xtable(resultswtab.spcumd, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleswspcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultswtab.spcumd, '/users/edong/Output/Male/Outcome/maleswspcumd.csv')
    print(xtable(resultuntab.spcumd, type = 'latex'), file = '/users/edong/Output/Male/Outcome/maleunspcumd.tex', sanitize.text.function = function(x){x})
    write.csv(resultuntab.spcumd, '/users/edong/Output/Male/Outcome/maleunspcumd.csv')
  }
  
}








