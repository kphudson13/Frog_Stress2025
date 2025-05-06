
library(tidyverse) #for everything 
library(gridExtra) #to export the stats tables
library(grid)
library(cowplot) #to arrange ggplots in a grid 


# Simulate data -----------------------------------------------------------


rm(list = ls()) #clear workspace
set.seed(123)

nsim <- 1000 #number of simulations
MSMRMods <- list() #parking lots for models 
CortMods <- list()
CortMSMR <- list()

#MSMR determined from Goetz 2018 figure 3 control 
#female average 0.170281, standard error 0.0107
#Male average 0.213707571809, standard error 0.0169689
MSMRSD <- 0.015*sqrt(25) #standard deviation = standard error * sqrt(n)

#Cort release rate from figure 4B Gabor 2013 
#SD is estimated by 3/4 x IQR
CortSD <- 2.85
pb <- progress_bar$new(total = length(nsim)) #for progress bar

##for gamma, from lab 5 
#scale parameter = variance/mean = coef.var. 
#shape parameter = mean^2/variance

N <- 72 #sample size
SimulatedData <- data.frame(ID = as.character(seq(1,N,1)), 
                            Temp = as.numeric(rep(c(14,18,22,26,30,34), length = N)), # six temp treatments 
                            sex = NA,
                            Mass = NA) 

#Mass estimates by sex from Goetz 2018 Table 1
for (i in 1:nsim){
  SimulatedData$sex = as.character(sample(c("M","F"), N, replace = TRUE)) #assume equal numbers
  for (j in 1:N) { 
    if (SimulatedData$sex[j] == "M") {
      SimulatedData$Mass[j] <- rgamma(1, shape = 6.8^2/0.4^2, scale = 0.4^2/6.8)
    } else {
      SimulatedData$Mass[j] <- rgamma(1, shape = 20^2/2^2, scale = 2^2/20)
    }
  }
  MSMR_det <- SimulatedData$Mass^-0.25 * exp(SimulatedData$Temp/10) * (0.2/6.04) #6.04 is the average MSMR before it is scaled to Goetz 2013
  Cort_det <- SimulatedData$Mass^-0.25 * exp(SimulatedData$Temp/10) * (7.8/6.04) #7.8 is the average MSMR before it is scaled to Gabor 2018 
  SimulatedData$MSMR <- rgamma(N, shape = (MSMR_det^2)/(MSMRSD^2), scale = (MSMRSD^2)/MSMR_det) 
  SimulatedData$Cort <- rgamma(N, shape = (Cort_det^2)/(CortSD^2), scale = (CortSD^2)/Cort_det) 
  MSMRModel <- lm(log(MSMR) ~ log(Mass) + Temp, data = SimulatedData)
  CortModel <- lm(log(Cort) ~ log(Mass) + Temp, data = SimulatedData)
  CortMSMRModel <- lm(log(Cort) ~ log(MSMR), data = SimulatedData)
  MSMRMods[[i]] <- MSMRModel
  CortMods[[i]] <- CortModel
  CortMSMR[[i]] <- CortMSMRModel
}

#store the final example to visualize one dataframe 
write.csv(SimulatedData, "Power_Analysis/SimulatedData.csv", row.names = FALSE)


# Stats table -------------------------------------------------------------


#blank dataframes for averages
MSMRaverage <- matrix(NA, nrow = 3, ncol = 4)
Cortaverage <- matrix(NA, nrow = 3, ncol = 4)
CortMSMRaverage <- matrix(NA, nrow = 2, ncol = 4)

#get averages across the 1000 simulated models
for (i in 1:3){
  for (j in 1:4){
    MSMRaverage[i, j] <- mean(sapply(1:nsim, function(n) coef(summary(MSMRMods[[n]]))[i, j]))
    MSMRrsq <- mean(sapply(1:nsim, function(n) summary(MSMRMods[[n]])$r.squared))
    Cortaverage[i, j] <- mean(sapply(1:nsim, function(n) coef(summary(CortMods[[n]]))[i, j]))
    Cortrsq <- mean(sapply(1:nsim, function(n) summary(CortMods[[n]])$r.squared))
  }
}

#confidence intervals 
MSMRConfint <- matrix(NA, nrow = 3, ncol = 2)
CortConfint <- matrix(NA, nrow = 3, ncol = 2)

for (i in 1:3){
  for (j in 1:2) {
    MSMRConfint[i, j] <- mean(sapply(1:nsim, function(n) confint(MSMRMods[[n]])[i, j]))
    CortConfint[i, j] <- mean(sapply(1:nsim, function(n) confint(CortMods[[n]])[i, j]))
  }
}

#Cort ~ MSMR model averages
for (i in 1:2){
  for (j in 1:4){
    CortMSMRaverage[i, j] <- mean(sapply(1:nsim, function(n) coef(summary(CortMSMR[[n]]))[i, j]))
        CortMSMRrsq <- mean(sapply(1:nsim, function(n) summary(CortMSMR[[n]])$r.squared))
  }
}

CortMSMRConfint <- matrix(NA, nrow = 2, ncol = 2)

for (i in 1:2){
  for (j in 1:2) {
    CortMSMRConfint[i, j] <- mean(sapply(1:nsim, function(n) confint(CortMSMR[[n]])[i, j]))
  }
}

#organize statistics in a way I like 
StatTab <- rbind(MSMRaverage, Cortaverage, CortMSMRaverage) %>%
  as.data.frame(.) %>%
  slice(-c(1,4,7)) %>% #remove intercept rows
  cbind(., c(MSMRConfint[2,1],
             MSMRConfint[3,1],
             CortConfint[2,1],
             CortConfint[3,1],
             CortMSMRConfint[2,1]), 
        c(MSMRConfint[2,2],
          MSMRConfint[3,2],
          CortConfint[2,2],
          CortConfint[3,2],
          CortMSMRConfint[2,2]), #add back in confidence intervals as column
        c(MSMRaverage[1,1],
          MSMRaverage[1,1],
          Cortaverage[1,1],
          Cortaverage[1,1],
          CortMSMRaverage[1,1]), #add back in intercepts as column
        c(MSMRrsq,
          MSMRrsq,
          Cortrsq,
          Cortrsq,
          CortMSMRrsq)) %>% #and r squared values 
  `colnames<-`(c("Estimate", "SE Est.", "T value",  "p value", "Lower 95% CI", "Upper 95% CI", "Intercept", "R squared")) %>%
  `rownames<-`(c("MSMR ~ mass",
                 "MSMR ~ Temp",
                 "Cort ~ mass", 
                 "Cort ~ Temp",
                 "Cort ~ MSMR")) %>%
  mutate(Intercept = as.numeric(Intercept),
         `R squared` = as.numeric(`R squared`)) %>% #make them numeric for rounding 
  mutate(Estimate = round(Estimate, 3),
         `SE Est.` = round(`SE Est.`, 3),
         `T value` = round(`T value`, 2),
         `p value` = round(`p value`, 3),
         Intercept = round(Intercept, 2),
         `R squared` = round(`R squared`, 3)) %>%
  mutate(`p value` = ifelse(`p value` < 0.001, "< 0.001", `p value`)) #change very small p values to < 0.001
  
  
StatTab$Intercept[c(2,4)] <- " " #repeated values 
StatTab$'R squared'[c(2,4)] <- " "
  

#set the theme for the table png
tt1 <- ttheme_minimal(rowhead=list(fg_params=list(fontface = "bold")))

#export stats table as a csv
write.csv(StatTab, "Power_Analysis/StatsTab.csv", row.names = TRUE)

#export stats table asa pretty png 
png("Power_Analysis/StatsTab.png", 
    height = 190*nrow(StatTab), 
    width = 430*ncol(StatTab),
    res = 300)
grid.newpage()
grid.table(StatTab, theme = tt1)
dev.off()  


# Plots -------------------------------------------------------------------

### BE WARNED ###
#this is not computationally efficient, but it works
#visualizing with the last plot_grid will take a while

#Start with an emptly plot to put expected lines on 
MSMRMassPlot <- ggplot() +
  geom_point(aes(x = log(SimulatedData$Mass), y = log(SimulatedData$MSMR)), color = NA) + #just need this to set the scale for the axis 
  labs(x = "ln Body Mass (g)",
       y = "ln MSMR") +
  scale_y_continuous(limits = c(-7, 0)) +
  theme_classic () 

#plot the 1000 lines for each model
i <- 1
while (i <= length(MSMRMods)){
  MSMRMassPlot <- MSMRMassPlot +
    geom_abline(intercept = summary(MSMRMods[[i]])$coefficients[1], 
                slope = summary(MSMRMods[[i]])$coefficients[2], #2nd row is mass
                color = "blue",
                alpha = 0.02) #fairly clear so we can see overlap
  i <- i + 1
}

MSMRTempPlot <- ggplot() +
  geom_point(aes(x = SimulatedData$Temp, y = log(SimulatedData$MSMR)), color = NA) + #just need this to set the scale for the axis 
  labs(x = "Temperature (C)",
       y = "ln MSMR") +
  theme_classic () 

i <- 1
while (i <= length(MSMRMods)){
  MSMRTempPlot <- MSMRTempPlot +
    geom_abline(intercept = summary(MSMRMods[[i]])$coefficients[1], 
                slope = summary(MSMRMods[[i]])$coefficients[3], #3rd row is temp
                color = "blue",
                alpha = 0.02) #fairly clear so we can see overlap
  i <- i + 1
}

CortMassPlot <- ggplot() +
  geom_point(aes(x = log(SimulatedData$Mass), y = log(SimulatedData$Cort)), color = NA) + #just need this to set the scale for the axis 
  labs(x = "ln Body Mass (g)",
       y = "ln Cort") +
  scale_y_continuous(limits = c(-5, 2)) +
  theme_classic ()

i <- 1
while (i <= length(CortMods)){
  CortMassPlot <- CortMassPlot +
    geom_abline(intercept = summary(CortMods[[i]])$coefficients[1], 
                slope = summary(CortMods[[i]])$coefficients[2], #2nd row is mass
                color = "blue",
                alpha = 0.02) #fairly clear so we can see overlap
  i <- i + 1
}

CortTempPlot <- ggplot() +
  geom_point(aes(x = SimulatedData$Temp, y = log(SimulatedData$Cort)), color = NA) + #just need this to set the scale for the axis 
  labs(x = "Temperature (C)",
       y = "ln Cort") +
  theme_classic ()

i <- 1
while (i <= length(CortMods)){
  CortTempPlot <- CortTempPlot +
    geom_abline(intercept = summary(CortMods[[i]])$coefficients[1], 
                slope = summary(CortMods[[i]])$coefficients[3], #3rd row is temp
                color = "blue",
                alpha = 0.02) #fairly clear so we can see overlap
  i <- i + 1
}

CortMSMRPlot <- ggplot() +
  geom_point(aes(x = log(SimulatedData$MSMR), y = log(SimulatedData$Cort)), color = NA) + #just need this to set the scale for the axis 
  labs(x = "ln MSMR",
       y = "ln Cort") +
  theme_classic ()

i <- 1
while (i <= length(CortMSMR)){
  CortMSMRPlot <- CortMSMRPlot +
    geom_abline(intercept = summary(CortMSMR[[i]])$coefficients[1], 
                slope = summary(CortMSMR[[i]])$coefficients[2], #2nd row is msmr
                color = "blue",
                alpha = 0.02) #fairly clear so we can see overlap
  i <- i + 1
}

#put all plots in a grid 
plot_grid(MSMRMassPlot, MSMRTempPlot, CortMassPlot, CortTempPlot, CortMSMRPlot, 
          labels = c('A', 'B', 'C', 'D', 'E'), label_size = 12)










