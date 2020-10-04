##################### Table 1 #####################
# Read in baseline, df
source("/users/edong/Code/Cleandat.R")
############# Take out the first dispense daily dose category of each person ############
library(dplyr)
df <- gc_final[, c("ikn", "Dtrt")]
df_initial <- df %>%
  group_by(ikn) %>%
  filter(row_number() == 1)

######################### Merge together #########################
dat <- inner_join(baseline, df_initial, by="ikn")
#View(dat)

dat_female <- dat[dat$sex == 0,]
dat_male <- dat[dat$sex == 1,]

############################# Table One, two sex tables stratified by initial dose ###########################
library(tableone)
myVars <- c( "age",     "falldr", "sexhor", "antisexhor", "totprevBPdaysc",  "pralcal", "pdmab",  "inhalegc",   "bronchodilators", 
             "thiazide",  "arthritis", "fallcomorb", "emerhosp", "diabetes",  "pfx", "status")
catVars <- c("totprevBPdaysc")

# Female table
tb1_female <- CreateTableOne(data=dat_female, test = FALSE,
                             vars = myVars, strata = "Dtrt", factorVars = catVars)
tabwrt_female <- print(tb1_female, smd=F, noSpaces=T, showAllLevels=T)
write.csv(tabwrt_female, file = "/users/edong/Output/Descriptive/tab1_female.csv")

# Male table
tb1_male <- CreateTableOne(data=dat_male, test = FALSE,
                           vars = myVars, strata = "Dtrt", factorVars = catVars)
tabwrt_male <- print(tb1_male, smd=F, noSpaces=T, showAllLevels=T)
write.csv(tabwrt_male, file = "/users/edong/Output/Descriptive/tab1_male.csv")

# # All table
# tb1 <- CreateTableOne(data=dat, test = FALSE,
#                              vars = myVars, strata = "Dtrt", factorVars = catVars)
# tabwrt <- print(tb1, smd=F, noSpaces=T, showAllLevels=T)
# write.csv(tabwrt, file = "/users/edong/Output/tab1.csv")
#
# ### Stratify by sex instead of initial dose, which is now a cat var, only one table
# myVars_bysex <- c( "age",      "diabetes", "arthritis",  "pfx",      "antisexhor",  "inhalegc",   "bronchodilators", "sexhor",         "Dtrt",
#                    "thiazide", "pralcal",  "pdmab",      "falldr",   "fallcomorb",  "emerhosp",   "prevlastBPdtcat", "totprevBPdaysc", "status")
# catVars_bysex <- c("prevlastBPdtcat", "totprevBPdaysc", "Dtrt")
# 
# ## Female table
# tb1_bysex <- CreateTableOne(data=dat, test = FALSE,
#                              vars = myVars_bysex, strata = "sex", factorVars = catVars_bysex)
# tabwrt_bysex <- print(tb1_bysex, smd=F, noSpaces=T, showAllLevels=T)
# write.csv(tabwrt_bysex, file = "/users/edong/Output/tab1_bysex.csv")



