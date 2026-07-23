# =========================================================
# Core statistics and plots for the 2026 study comparing
# metabolic rate and corticosterone in frogs
#
# Author: Kyle Hudson
# Circa 2026
# =========================================================

library(tidyverse)
library(grid)
library(gridExtra) # for stats table theme
library(rsq)

rm(list=ls()) #clear environment

Data <- read.csv("Data_Spreadsheet.csv") %>%
  mutate(VCO2 = ifelse(VCO2 <= 0, NA, VCO2)) %>%
  mutate(mW = VCO2 * 21.1) #convert vco2 to watts

dir.create("Figures", showWarnings = FALSE) # create directory for figures

# VCO2 Models ------------------------------------------------------------------

MSMRWeightTemp_Model <- lm(data = Data, log(mW/Weight) ~ log(Weight) + Temperature)
summary(MSMRWeightTemp_Model)
CortMSMR_Model <- lm(data = Data, log(Cort) ~ log(mW/Weight)) 
summary(CortMSMR_Model)
CortWeightTemp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature)
summary(CortWeightTemp_Model)

# Plots -------------------------------------------------------------------

theme1 <- theme(legend.background = element_blank(),
                legend.title = element_blank(),
                legend.key.height = unit(0.6, "lines"),
                legend.position=c(.8,0.1)) # set standard theme for plots

# normalized points for plotting
Data$MSMR_NormTemp <- log(Data$mW/Data$Weight) - MSMRWeightTemp_Model[["coefficients"]][["Temperature"]]*Data$Temperature

(MSMRWeight_Plot <- ggplot(Data, aes(x = log(Weight), y = MSMR_NormTemp)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
                slope = coefficients(summary(MSMRWeightTemp_Model))[2,1]) + # pull lines from model
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 1, y = -10,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(MSMRWeightTemp_Model))[2,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(MSMRWeightTemp_Model)$partial.rsq[1], 2))))),
             parse = TRUE) +
    scale_x_continuous(limits = c(-1.5, 4)) +
    scale_y_continuous(limits = c(-12, -6))
)

ggsave(filename = "Figures/MSMRWeight_Plot.png", width=90, height=90, units="mm") #save a picture

# normalized points for plotting
Data$MSMR_NormWeight <- log(Data$mW/Data$Weight) - MSMRWeightTemp_Model[["coefficients"]][["log(Weight)"]]*log(Data$Weight)

(MSMRTemp_Plot <- ggplot(data = Data, aes(x = Temperature, y = MSMR_NormWeight)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
                slope = coefficients(summary(MSMRWeightTemp_Model))[3,1]) +
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 22, y = -10,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(MSMRWeightTemp_Model))[3,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(MSMRWeightTemp_Model)$partial.rsq[2], 2))))),
             parse = TRUE)+
    scale_x_continuous(limits = c(12, 36)) +
    scale_y_continuous(limits = c(-12, -6))
)

ggsave(filename = "Figures/MSMRTemp_Plot.png", width=90, height=90, units="mm") #save a picture

(CortMSMR_Plot <- ggplot(data = Data, aes(x=log(mW/Weight), y = log(Cort))) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(CortMSMR_Model))[1,1],
                slope = coefficients(summary(CortMSMR_Model))[2,1]) + 
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = -10, y = 6,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(CortMSMR_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(CortMSMR_Model))[2,1], 2)),
                                      ~R^2 ==~ .(round(summary(CortMSMR_Model)$r.squared, 2))))),
             parse = TRUE) +
    scale_x_continuous(limits = c(-11, -6)) +
    scale_y_continuous(limits = c(-4, 8))
)

ggsave(filename = "Figures/CortMSMR_Plot.png", width=90, height=90, units="mm") #save a picture

# normalized points for plotting
Data$Cort_NormTemp <- log(Data$Cort) - CortWeightTemp_Model[["coefficients"]][["Temperature"]]*Data$Temperature

(CortWeight_Plot <- ggplot(data = Data, aes(x=log(Weight), y = Cort_NormTemp)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(CortWeightTemp_Model))[1,1],
                slope = coefficients(summary(CortWeightTemp_Model))[2,1]) +
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 3, y = 5,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(CortWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(CortWeightTemp_Model))[2,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(CortWeightTemp_Model)$partial.rsq[1], 2))))),
             parse = TRUE) +
    scale_x_continuous(limits = c(-1, 4)) +
    scale_y_continuous(limits = c(-2, 6))
)

ggsave(filename = "Figures/CortWeight_Plot.png", width=90, height=90, units="mm") #save a picture

# normalized points for plotting
Data$Cort_NormWeight <- log(Data$Cort) - CortWeightTemp_Model[["coefficients"]][["log(Weight)"]]*log(Data$Weight)

(CortTemp_Plot <- ggplot(data = Data, aes(x=Temperature, y = Cort_NormWeight)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(CortWeightTemp_Model))[1,1],
                slope = coefficients(summary(CortWeightTemp_Model))[3,1]) +
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 23, y = 5,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(CortWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(CortWeightTemp_Model))[3,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(CortWeightTemp_Model)$partial.rsq[2], 2))))),
             parse = TRUE) +
    scale_x_continuous(limits = c(12, 36)) +
    scale_y_continuous(limits = c(-2, 6))
)

ggsave(filename = "Figures/CortTemp_Plot.png", width=90, height=90, units="mm") #save a picture

# Stats Table -------------------------------------------------------------


Stats_Tab <- rbind(coefficients(summary(MSMRWeightTemp_Model)),
                   coefficients(summary(CortWeightTemp_Model)),
                   coefficients(summary(CortMSMR_Model))) %>%
  as.data.frame(.) %>%
  slice(-c(1,4,7)) %>%  #cut out all the rows of intercept stats
  select(., -"t value") %>%
  cbind(., c(rsq.partial(MSMRWeightTemp_Model)$partial.rsq, 
             rsq.partial(CortWeightTemp_Model)$partial.rsq, 
             summary(CortMSMR_Model)$r.squared)) %>% # bind in r sq
  mutate(across(c(1:4), \(x) round(x, digits = 2))) %>% #new way to round w/ anonymous function
  `colnames<-`(c("Estimate", "SE (Slope)", "p value", "R2")) %>%
  `rownames<-`(c("MSMR ~ Weight", "MSMR ~ Temp", "Cort ~ Weight", "Cort ~ Temp", "Cort ~ MSMR"))

tt1 <- ttheme_default(rowhead=list(fg_params=list(fontface = "bold"),
                                   bg_params=list(fill="grey80"))) # theme for stats table

write.csv(Stats_Tab, file = "Figures/StatsTab.csv", row.names = TRUE)

#export stats table 
png("Figures/StatsTab.png",
    height = 180*nrow(Stats_Tab), 
    width = 500*ncol(Stats_Tab),
    res = 300)
grid.newpage()
grid.table(Stats_Tab, theme = tt1)
grid.text("Stats Table", x = 0.2, y = 0.9, gp = gpar(fontface = "bold"))
dev.off()



