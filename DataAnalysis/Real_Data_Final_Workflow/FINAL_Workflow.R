# Load libraries
source("/users/edong/Code/Level4_Clean.R", echo=T)
doselineplot = FALSE
nonadmcensplt = FALSE
source("/users/edong/Code/Exposure_dataset.R", echo=T)

female = FALSE
gend <- ifelse(female, 0, 1)
gender.name <- ifelse(female, 'Female', 'Male')
#source("/users/edong/Code/Death_Fx_Time_plot.R", echo=T) # need to load plyr density comparison of death/fx(female/male); /doses(all)

# subset to gender
D0_gender <- df_Dinit[df_Dinit$sex == gend, ] # male n=____, actual sampled in casebase: ____
V_gender <- long_V_slm[long_V_slm$sex == gend, ]
D_gender <- Ddata[Ddata$sex == gend, ] 

library(MASS)
if(female){
  source("/users/edong/Code/Exposure_models_female.R", echo=T); requireoutput = TRUE; weightplot = TRUE
  #source("/users/edong/Code/Exposure_models_female_nonum.R", echo=T); requireoutput = FALSE; weightplot = FALSE # for SMD diag plot only
}else{
  source("/users/edong/Code/Exposure_models.R", echo=T); requireoutput = TRUE; weightplot = TRUE
  #source("/users/edong/Code/Exposure_models_nonum.R", echo=T); requireoutput = FALSE; weightplot = FALSE # for SMD diag plot only
}

################################### Can save exposure model output now #################################
###################### Exposure Model outputs (could comment out from entire workflow) #################
library(xtable)
source("/users/edong/Code/Model_Outputs_names.R", echo=T)
#requireoutput = TRUE
if(requireoutput){
  library(forestplot)
  expo = TRUE; outcome = FALSE
  source("/users/edong/Code/Model_Outputs.R", echo=T) # forest plot has to run in RStudio, due to package
}

################## Outcome weighting/GCcumdose/BPcumdays mapping ##############
source("/users/edong/Code/Outcome_dataset.R", echo=T)

################# Weight Diagnostics ###############
library(Hmisc)
library(dplyr)
#weightplot = TRUE
if(weightplot){
  trunc <- TRUE         # for sw over time plot, should check both
  source("/users/edong/Code/Diagnostics.R", echo=T)  
}
library(tidyverse)
library(ggplot2)
source("/users/edong/Code/DiagnosticSMD.R", echo=T) # SMD has to run in Command, due to package

##################### Outcome models ###################
source("/users/edong/Code/Outcome_models.R", echo=T)

library(lmtest) 
########## Likelihood Ratio Test #########
# Likelihood ratio test between each one of the three weighted outcome models 
# and their own null model
lrtest(cbout.sw.null.a, cbout.sw.slm.a)
lrtest(cbout.sw.null.cumd, cbout.sw.slm.cumd)
lrtest(cbout.sw.null.spcumd, cbout.sw.slm.spcumd)
###################################
requireoutput = TRUE
if(requireoutput){
  expo = FALSE; outcome = TRUE
  source("/users/edong/Code/Model_Outputs.R", echo=T)
}

########### Plot risk over time & Plot potential outcomes hazards ############
modelnames <- c('Dose Level Model', 'Cumulative Dose Model', 'Flexible Cumulative Dose Model', NA)
#whichmodel <- 'Dose Level Model' 
whichmodel <- NA
source("/users/edong/Code/Hazard_over_time_plt.R", echo=T) 





