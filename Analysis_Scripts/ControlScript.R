# =========================================================
# Control script to test different data combinations
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================


library(tidyverse)
library(grid)
library(gridExtra) # for stats table theme
library(rsq) # for partial r2

rm(list=ls()) # Clear environment
source("Analysis_Scripts/DirectoryFunction.R")
source("Analysis_Scripts/CooksDistFunction.R")

Whole_data <- read.csv("Data_Spreadsheet.csv") %>%
  mutate(VCO2 = ifelse(VCO2 <= 0, NA, VCO2),
         mW = VCO2 * 21.1 / 60 * 1000, # energetic conversion(jouls per ml), min to sec, and Watts to mW
         Reproductive = ifelse(Reproductive == "Gravid", Reproductive, "Other"),
         Sex = ifelse(is.na(Sex) | Sex == "", "Other", Sex),
         MSMR = mW/Weight) 

# All Data ----------------------------------------------------------------

Data <- Whole_data
directory <- "AllData" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")

# No high CV --------------------------------------------------------------

Data = Whole_data %>%
  filter(Notes != "High CV")

directory <- "FilterCV" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")

# No high CV or single well failure ---------------------------------------

Data = Whole_data %>%
  filter(Notes != "High CV" & Notes != "Single well failiure")

directory <- "FilterCVSWF" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")
