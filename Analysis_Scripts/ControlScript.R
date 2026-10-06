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
library(emmeans)
library(cowplot) #to combine plots 
library(janitor) # for extra plots script

rm(list=ls()) # Clear environment
source("Analysis_Scripts/DirectoryFunction.R")
source("Analysis_Scripts/CooksDistFunction.R")
source("Analysis_Scripts/AssignLettersFunction.R")

Whole_data <- read.csv("Data_Spreadsheet.csv") %>%
  mutate(VCO2 = ifelse(VCO2 <= 0, NA, VCO2),
         mW = VCO2 * 21.1 / 60 * 1000, # energetic conversion(jouls per ml), min to sec, and Watts to mW
         MSMR = VCO2/Weight,
         Reproductive = ifelse(Reproductive == "Gravid", Reproductive, "Other"),
         Sex = ifelse(is.na(Sex) | Sex == "", "Other", Sex),
         JulianDate = as.numeric(format(as.Date(CaptureDate, format = "%B %d %Y"), "%j")),
         ReproSex = case_when( Reproductive == "Gravid" ~ "Gravid",
                               !is.na(Sex) & Sex != "" ~ Sex,
                               TRUE ~ "Other"))

theme1 <- theme(
  legend.background = element_blank(), # make clear 
  legend.title = element_blank(), # remove title
  legend.key.height = unit(0.6, "lines"), # fit legend lines closer together
  legend.position=c(.75, 0.1)
) # set standard theme for plots

# All Data ----------------------------------------------------------------

Data <- Whole_data
directory <- "AllData" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")
source("Analysis_Scripts/ModelSelectionScript.R")
source("Analysis_Scripts/SummaryStatsScript.R")

# No high CV --------------------------------------------------------------

Data = Whole_data %>%
  filter(Notes != "High CV")

directory <- "FilterCV" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")
source("Analysis_Scripts/ModelSelectionScript.R")
source("Analysis_Scripts/SummaryStatsScript.R")

# No high CV or single well failure ---------------------------------------

Data = Whole_data %>%
  filter(Notes != "High CV" & Notes != "Single well failure")

directory <- "FilterCVSWF" 
CreateDR(directory)

source("Analysis_Scripts/AnalysisScript.R")
source("Analysis_Scripts/ModelSelectionScript.R")
source("Analysis_Scripts/SummaryStatsScript.R")

