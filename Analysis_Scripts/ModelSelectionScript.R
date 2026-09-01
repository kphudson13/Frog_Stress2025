# =========================================================
# Model selection for the 2026 study comparing
# metabolic rate and corticosterone in frogs
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================

library(tidyverse)

rm(list=ls()) # Clear environment

Data <- read.csv("Data_Spreadsheet.csv") %>%
  mutate(VCO2 = ifelse(VCO2 <= 0, NA, VCO2),
         mW = VCO2 * 21.1 / 60 * 100, # convert vco2 to watts
         Reproductive = ifelse(Reproductive == "Gravid", Reproductive, "Other"),
         Sex = ifelse(is.na(Sex) | Sex == "", "Other", Sex), 
         ReproSex = case_when(
           Reproductive == "Gravid" ~ "Gravid",
           !is.na(Sex) & Sex != "" ~ Sex,
           TRUE ~ "Other"),
         JulianDate = as.numeric(format(as.Date(CaptureDate, format = "%B %d %Y"), "%j")),
         MSMR = mW/Weight)

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

MSMRWeightTemp_BIC <- BIC(MSMRWeightTemp_Model, 
                          MSMRWeightTemp_Rep_Model, 
                          MSMRWeightTemp_Sex_Model, 
                          MSMRWeightTemp_RepSex_Model,
                          MSMRWeightTemp_Spp_Model)

# Cort models 
CortWeightTemp_AIC <- AIC(CortWeightTemp_Model, 
                          CortWeightTemp_Rep_Model, 
                          CortWeightTemp_Sex_Model,
                          CortWeightTemp_RepSex_Model,
                          CortWeightTemp_Spp_Model)

CortWeightTemp_BIC <- BIC(CortWeightTemp_Model, 
                          CortWeightTemp_Rep_Model, 
                          CortWeightTemp_Sex_Model,
                          CortWeightTemp_RepSex_Model,
                          CortWeightTemp_Spp_Model)

# Cort MSMR models
CortMSMR_AIC <- AIC(CortMSMR_Model, 
                    CortMSMR_Rep_Model, 
                    CortMSMR_Sex_Model,
                    CortMSMR_RepSex_Model,
                    CortMSMR_Spp_Model)

CortMSMR_BIC <- BIC(CortMSMR_Model, 
                    CortMSMR_Rep_Model, 
                    CortMSMR_Sex_Model,
                    CortMSMR_RepSex_Model,
                    CortMSMR_Spp_Model)
