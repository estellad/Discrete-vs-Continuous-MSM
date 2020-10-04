####### Plot SMD Before and After Weighting #######
omnismd <- function(data, weight, vars, factorvars, smdplottitle, strata, col, smdplotname, labelxaxis){
  ################## Unweighted baseline #################
  (smd.table <- CreateTableOne(
    vars = vars,
    data = data, strata = strata, 
    factorVars = factorvars))
  
  smd.values <- as.data.frame(cbind(`Comparison` = rownames(t(ExtractSmd(smd.table))[-1, ]),
                                    t(ExtractSmd(smd.table))[-1, ])) %>%
    mutate(Comparison = as.character(Comparison)) %>%
    mutate(`Drug1` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[1])),
           `Drug2` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[3]))) %>%
    gather(Variable, SMD, agescale:arthritis, factor_key = T) %>%
    rename(`Covariates` = Variable) %>%
    mutate(`SMD` = as.numeric(SMD))
  
  ################## Weighted baseline ###################
  wgt.data <- svydesign(ids=~1, data = data, weights = ~weight)
  
  (smd.wgt.table <- svyCreateTableOne(
    vars = vars,
    data = wgt.data, strata = strata, 
    factorVars = factorvars))
  
  wgt.smd.values <- as.data.frame(cbind(`Comparison` = rownames(t(ExtractSmd(smd.wgt.table))[-1, ]),
                                        t(ExtractSmd(smd.wgt.table))[-1, ])) %>%
    mutate(Comparison = as.character(Comparison)) %>%
    mutate(`Drug1` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[1])),
           `Drug2` = sapply(Comparison, function(x) as.integer(unlist(strsplit(x, " "))[3]))) %>%
    gather(Variable, SMD, agescale:arthritis, factor_key = T) %>%
    rename(`Covariates` = Variable) %>%
    mutate(`SMD` = as.numeric(SMD))
  
  
  # Combined SMD Results
  smd.all.a <- smd.values %>% mutate(`State` = 'Pre',
                                     `Category` = paste0('Pre-', `Covariates`))
  smd.all.b <- wgt.smd.values %>% mutate(`State` = 'Post',
                                         `Category` = paste0('Pre-', `Covariates`))
  smd.all <- as.data.frame(rbind(smd.all.a, smd.all.b))
  smd.all$State <- factor(smd.all$State)
  smd.all$State <- relevel(smd.all$State, ref = 'Pre')
  smd.all$Category <- factor(smd.all$Category)
  
  
  ### Function plot before and after SMD box plot
  smd.both.vplot <- function(data, title= NULL, col, labelxaxis){
    w = 9.75
    h = 9.75
    
    v.plot <- ggplot(data, aes(x=`Covariates`, y=SMD)) + 
      geom_boxplot(aes(fill = `State`)) + 
      
      coord_flip() + theme_bw() + 
      scale_x_discrete(limits = rev(levels(data$`Covariates`)), 
                       labels = labelxaxis) +
      labs(x = 'Covariates', y = 'Pairwise Standardized Mean Difference') + 
      scale_fill_manual('Weighting', labels = c('Unweighted', 'IPW'), values = col) + 
      theme(axis.text.x = element_text(size=10, margin=margin(t=10)),
            axis.text.y = element_text(size=10, margin=margin(t=10)),
            axis.title.x = element_text(size=12, margin=margin(t=24)),
            axis.title.y = element_text(size=12, margin=margin(t=24)),
            plot.title = element_text(hjust = 0.5, face = 'bold', size = 14, margin = margin(b=24)),
            legend.title = element_text(size = 10, face = 'bold'),
            legend.text = element_text(size = 10),
            
            # Specifies Top Right Corner (1,1 is coordinate for proportionate location on graph)
            legend.position=c(1,1),
            legend.justification=c(1,1),
            
            #White Fill, Black Border
            legend.background = element_rect(color = 'black')) +
      ggtitle(title) + 
      geom_hline(linetype = 2, col = 'black', yintercept=0.2) +
      geom_hline(linetype = 2, col = 'black', yintercept=0.5) +
      geom_hline(linetype = 2, col = 'black', yintercept=0.8)
    print(v.plot)
  }
  
  if(female){
    path <- "/users/edong/Output/Female/Diagnostics/"
  }else{
    path <- "/users/edong/Output/Male/Diagnostics/"
  }
  pdf(paste0(path, smdplotname), width=6.5, height=6.5)
  smd.both.vplot(data = smd.all, 
                 title = smdplottitle,
                 col =  col,
                 labelxaxis = labelxaxis)
  dev.off()
}

######### Correct column types #########
### Common params ###
# TP 0
male_vars <- c( "agescale", "falldr", "sexhor", "totprevBPdaysc",
                "inhalegc", "thiazide", 
                "fallcomorb", "emerhosp", "pfx", "arthritis")

male_facvars <- c("falldr", "sexhor", "totprevBPdaysc",
                  "inhalegc", "thiazide", 
                  "fallcomorb", "emerhosp", "pfx", "arthritis")

male_labelxaxis <- c('agescale' = 'Age', 'falldr' = 'Fall related drugs use', 
                     'sexhor'= 'Sex hormone', 'totprevBPdaysc' = 'Prev. BP duration', 
                     'inhalegc' = 'Inhale GC', 'thiazide' = 'Thiazide', 
                     'fallcomorb' = 'Fall related conditions', 
                     'emerhosp' = 'Emergency or hospitalization', 'pfx' = 'Prev. fracture', 
                     'arthritis' = 'Arthritis')

# TP 123
male_varsTD <- c(male_vars, 'cumBPdayspscale', 'Ltc')

male_facvarsTD <- c(male_facvars, 'Ltc')

male_labelxaxisTD <- c(male_labelxaxis,
                       'cumBPdayspscale' = 'Cumulative BP days', 
                       'Ltc' = 'LTC')

#################
# TP 0
female_vars <- c(male_vars, "pralcal",  "diabetes")

female_facvars <- c(male_facvars, "pralcal", "diabetes")

female_labelxaxis <- c(male_labelxaxis, 
                       "pralcal" = 'Prev. raloxifene or calcitonin',
                       'diabetes' = 'Diabetes')
# TP 123
female_varsTD <- c(male_varsTD, "pralcal", "diabetes")

female_facvarsTD <- c(male_facvarsTD, "pralcal", "diabetes")

female_labelxaxisTD <- c(male_labelxaxisTD,
                         "pralcal" = 'Prev. raloxifene or calcitonin',
                         'diabetes' = 'Diabetes')


if(female){ 
  vars = female_vars
  factorvars = female_facvars
  labelxaxis = female_labelxaxis
  varsTD = female_varsTD
  factorvarsTD = female_facvarsTD
  labelxaxisTD = female_labelxaxisTD
}else{
  vars = male_vars
  factorvars = male_facvars
  labelxaxis = male_labelxaxis
  varsTD = male_varsTD
  factorvarsTD = male_facvarsTD
  labelxaxisTD = male_labelxaxisTD
}

# TP 0
cb_sub_fac0 <- long_Y_cb[, c('study_id', factorvars, 'D0')]
cb_sub_fac0 <- apply(cb_sub_fac0, 2, as.factor)
cb_sub_num0 <- long_Y_cb[, c(setdiff(vars, factorvars), 'End','d0p')]

cb_sub0 <- cbind(cb_sub_num0, cb_sub_fac0)

# TP 123
cb_sub_fac <- long_Y_cb[, c('study_id', factorvarsTD, 'Dtrt')]
cb_sub_fac <- apply(cb_sub_fac, 2, as.factor)
cb_sub_num <- long_Y_cb[, c(setdiff(varsTD, factorvarsTD), 'End','sw.trunc')]

cb_sub <- cbind(cb_sub_num, cb_sub_fac)

########################## Range-wise 0-100, 101-183, 184-365 ######################
tp0 = 'Baseline'
tp1 = '0< Time Points <= 100 Days'
tp2 = '100 < Time Points <= 183 Days'
tp3 = '183 < Time Points <= 365 Days'

long_Y_cb_baseline_sw <- cb_sub0 %>% # combined weights truncated, d0 weight not truncated.
  group_by(study_id) %>%
  filter(row_number() == 1)


long_Y_cb_init_sw <- cb_sub %>%
  group_by(study_id) %>%
  filter(End <= 100 & End > 0)
summary(long_Y_cb_init_sw$End) 


long_Y_cb_med_sw <- cb_sub %>%
  group_by(study_id) %>%
  filter(End <= 183 & End >100)
summary(long_Y_cb_med_sw$End) 


long_Y_cb_last_sw <- cb_sub %>%
  group_by(study_id) %>%
  filter(End <= 365 & End >183)
summary(long_Y_cb_last_sw$End)  

##### TP 0
omnismd(data = long_Y_cb_baseline_sw, weight = long_Y_cb_baseline_sw$d0p, vars = vars, factorvars = factorvars, 
        smdplottitle = paste0('Pairwise SMD of Baseline and TD Covariates \n at ', tp0, ' \n Before and After IP Weighting - ', gender.name), 
        strata = 'D0', col = c('orangered', 'midnightblue'), smdplotname = paste0('smd_tp0_sw_',gender.name,'r.pdf'), labelxaxis = labelxaxis)
##### TP 1
omnismd(data = long_Y_cb_init_sw, weight = long_Y_cb_init_sw$sw.trunc, vars = varsTD, factorvars = factorvarsTD, 
        smdplottitle = paste0('Pairwise SMD of Baseline and TD Covariates \n at ', tp1, ' \n Before and After IP Weighting - ', gender.name), 
        strata = 'Dtrt', col = c('orange', 'navy'), smdplotname = paste0('smd_tp1_sw_',gender.name,'r.pdf'), labelxaxis = labelxaxisTD)
##### TP 2
omnismd(data = long_Y_cb_med_sw, weight = long_Y_cb_med_sw$sw.trunc, vars = varsTD, factorvars = factorvarsTD, 
        smdplottitle = paste0('Pairwise SMD of Baseline and TD Covariates \n at ', tp2, ' \n Before and After IP Weighting - ', gender.name), 
        strata = 'Dtrt', col = c('gold', 'dodgerblue'), smdplotname = paste0('smd_tp2_sw_',gender.name,'r.pdf'), labelxaxis = labelxaxisTD)
##### TP 3
omnismd(data = long_Y_cb_last_sw, weight = long_Y_cb_last_sw$sw.trunc, vars = varsTD, factorvars = factorvarsTD, 
        smdplottitle = paste0('Pairwise SMD of Baseline and TD Covariates \n at ', tp3, ' \n Before and After IP Weighting - ', gender.name), 
        strata = 'Dtrt', col = c('yellow', 'skyblue1'), smdplotname = paste0('smd_tp3_sw_',gender.name,'r.pdf'), labelxaxis = labelxaxisTD)




