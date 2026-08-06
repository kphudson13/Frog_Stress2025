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
  mutate(VCO2 = ifelse(VCO2 <= 0, NA, VCO2)) %>%
  mutate(mW = VCO2 * 21.1) %>% # convert vco2 to watts
  mutate(Reproductive = ifelse(Reproductive == "Gravid", Reproductive, "Other"))

# Build models ------------------------------------------------------------

# Simple models 
MSMRWeightTemp_Model <- lm(data = Data, log(mW/Weight) ~ log(Weight) + Temperature)
summary(MSMRWeightTemp_Model)
CortMSMR_Model <- lm(data = Data, log(Cort) ~ log(mW/Weight)) 
summary(CortMSMR_Model)
CortWeightTemp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature)
summary(CortWeightTemp_Model)

# Reproduction models
MSMRWeightTemp_Rep_Model <- lm(data = Data, log(mW/Weight) ~ log(Weight) + Temperature + Reproductive)
summary(MSMRWeightTemp_Rep_Model)
CortMSMR_Rep_Model <- lm(data = Data, log(Cort) ~ log(mW/Weight) + Reproductive) 
summary(CortMSMR_Rep_Model)
CortWeightTemp_Rep_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Reproductive)
summary(CortWeightTemp_Rep_Model)

# Species models 
MSMRWeightTemp_Spp_Model <- lm(data = Data, log(mW/Weight) ~ log(Weight) + Temperature + Species)
summary(MSMRWeightTemp_Spp_Model)
CortMSMR_Spp_Model <- lm(data = Data, log(Cort) ~ log(mW/Weight) + Species) 
summary(CortMSMR_Spp_Model)
CortWeightTemp_Spp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Species)
summary(CortWeightTemp_Spp_Model)

# Model selection ---------------------------------------------------------

# MSMR models
MSMRWeightTemp_AIC <- AIC(MSMRWeightTemp_Model, MSMRWeightTemp_Rep_Model, MSMRWeightTemp_Spp_Model)
MSMRWeightTemp_BIC <- BIC(MSMRWeightTemp_Model, MSMRWeightTemp_Rep_Model, MSMRWeightTemp_Spp_Model)

# Cort models 
CortWeightTemp_AIC <- AIC(CortWeightTemp_Model, CortWeightTemp_Rep_Model, CortWeightTemp_Spp_Model)
CortWeightTemp_BIC <- BIC(CortWeightTemp_Model, CortWeightTemp_Rep_Model, CortWeightTemp_Spp_Model)

# Cort MSMR models
CortMSMR_AIC <- AIC(CortMSMR_Model, CortMSMR_Rep_Model, CortMSMR_Spp_Model)
CortMSMR_BIC <- BIC(CortMSMR_Model, CortMSMR_Rep_Model, CortMSMR_Spp_Model)


