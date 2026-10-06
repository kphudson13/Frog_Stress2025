# =========================================================
# Model selection for the 2026 study comparing
# metabolic rate and corticosterone in frogs
#
# This script should be ran from control script, if not 
# you'll need to manually load in data
# 
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================



# Build models ------------------------------------------------------------

# Simple models 
MSMRWeightTemp_Model <- lm(data = Data, log(MSMR) ~ log(Weight) + Temperature)
summary(MSMRWeightTemp_Model)
CortWeightTemp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature)
summary(CortWeightTemp_Model)
CortMSMR_Model <- lm(data = Data, log(Cort) ~ log(MSMR)) 
summary(CortMSMR_Model)

# Reproduction only models
MSMRWeightTemp_Rep_Model <- lm(data = Data, log(MSMR) ~ log(Weight) + Temperature + Reproductive)
summary(MSMRWeightTemp_Rep_Model)
CortWeightTemp_Rep_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Reproductive)
summary(CortWeightTemp_Rep_Model)
CortMSMR_Rep_Model <- lm(data = Data, log(Cort) ~ log(MSMR) + Reproductive) 
summary(CortMSMR_Rep_Model)

# Sex only models 
MSMRWeightTemp_Sex_Model <- lm(data = Data, log(MSMR) ~ log(Weight) + Temperature + Sex)
summary(MSMRWeightTemp_Sex_Model)
CortWeightTemp_Sex_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Sex)
summary(CortWeightTemp_Sex_Model)
CortMSMR_Sex_Model <- lm(data = Data, log(Cort) ~ log(MSMR) + Sex) 
summary(CortMSMR_Sex_Model)

#ReproSex models 
MSMRWeightTemp_RepSex_Model <- lm(data = Data, log(MSMR) ~ log(Weight) + Temperature + ReproSex)
summary(MSMRWeightTemp_RepSex_Model)
CortWeightTemp_RepSex_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + ReproSex)
summary(CortWeightTemp_RepSex_Model)
CortMSMR_RepSex_Model <- lm(data = Data, log(Cort) ~ log(MSMR) + ReproSex) 
summary(CortMSMR_RepSex_Model)

# Species models 
MSMRWeightTemp_Spp_Model <- lm(data = Data, log(MSMR) ~ log(Weight) + Temperature + Species)
summary(MSMRWeightTemp_Spp_Model)
CortWeightTemp_Spp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Species)
summary(CortWeightTemp_Spp_Model)
CortMSMR_Spp_Model <- lm(data = Data, log(Cort) ~ log(MSMR) + Species) 
summary(CortMSMR_Spp_Model)

# Model selection ---------------------------------------------------------

# MSMR models
MSMRWeightTemp_AIC <- AIC(MSMRWeightTemp_Model, 
                          MSMRWeightTemp_Rep_Model, 
                          MSMRWeightTemp_Sex_Model, 
                          MSMRWeightTemp_RepSex_Model,
                          MSMRWeightTemp_Spp_Model)

# MSMRWeightTemp_BIC <- BIC(MSMRWeightTemp_Model, 
#                           MSMRWeightTemp_Rep_Model, 
#                           MSMRWeightTemp_Sex_Model, 
#                           MSMRWeightTemp_RepSex_Model,
#                           MSMRWeightTemp_Spp_Model)

# Cort models 
CortWeightTemp_AIC <- AIC(CortWeightTemp_Model, 
                          CortWeightTemp_Rep_Model, 
                          CortWeightTemp_Sex_Model,
                          CortWeightTemp_RepSex_Model,
                          CortWeightTemp_Spp_Model)

# CortWeightTemp_BIC <- BIC(CortWeightTemp_Model, 
#                           CortWeightTemp_Rep_Model, 
#                           CortWeightTemp_Sex_Model,
#                           CortWeightTemp_RepSex_Model,
#                           CortWeightTemp_Spp_Model)

# Cort MSMR models
CortMSMR_AIC <- AIC(CortMSMR_Model, 
                    CortMSMR_Rep_Model, 
                    CortMSMR_Sex_Model,
                    CortMSMR_RepSex_Model,
                    CortMSMR_Spp_Model)

# CortMSMR_BIC <- BIC(CortMSMR_Model, 
#                     CortMSMR_Rep_Model, 
#                     CortMSMR_Sex_Model,
#                     CortMSMR_RepSex_Model,
#                     CortMSMR_Spp_Model)




# AIC Table ---------------------------------------------------------------

AIC_table <- as.data.frame(rbind(CortMSMR_AIC, CortWeightTemp_AIC, MSMRWeightTemp_AIC)) %>%
  mutate(across(2, \(x) round(x, digits = 2))) %>%
  add_row(.before = 1) %>% add_row(.before = 7) %>% add_row(.before = 13) %>% #add blank rows to divide cortisol and crtstn and blank model
  mutate(Model = c(NA,
                   "Corticosterone ~ MSMR", 
                   "Corticosterone ~ MSMR + Reproductive", 
                   "Corticosterone ~ MSMR + Sex", 
                   "Corticosterone ~ MSMR + ReproSex",
                   "Corticosterone ~ MSMR + Species",
                   NA,
                   "Corticosterone ~ Temp + Mass",
                   "Corticosterone ~ Temp + Mass + Reproductive",
                   "Corticosterone ~ Temp + Mass + Sex",
                   "Corticosterone ~ Temp + Mass + Reprosex",
                   "Corticosterone ~ Temp + Mass + Species",
                   NA,
                   "MSMR ~ Temp + Mass",
                   "MSMR ~ Temp + Mass + Reproductive",
                   "MSMR ~ Temp + Mass + Sex",
                   "MSMR ~ Temp + Mass + Reprosex",
                   "MSMR ~ Temp + Mass + Species"
                   )) %>%
  `rownames<-`(NULL) %>%
  select(3,1,2)

AIC_table[is.na(AIC_table)] <- " "  # remove NAs

write.csv(AIC_table, file = paste("Figures/", directory, "/AIC_Table.csv", sep = ""))



