
library(tidyverse)
library(progress)

#For sources of a priori data see script 'Multi_sim.R'

rm(list = ls()) #clear workspace
set.seed(123)

#vary sample size 
nsim <- 1000 #number of simulations
numvec <- seq(6, 250, by=3) #sample size 
MSMRMassPowerVec <- numeric(length(numvec)) #A bunch of parking lots
MSMRTempPowerVec <- numeric(length(numvec))
MSMRMasspval <- numeric(nsim)
MSMRTemppval <- numeric(nsim)
CortMassPowerVec <- numeric(length(numvec))
CortTempPowerVec <- numeric(length(numvec))
CortMasspval <- numeric(nsim)
CortTemppval <- numeric(nsim)
CortMSMRPowerVec <- numeric(length(numvec))
CortMSMRpval <- numeric(nsim)

MSMRSD <- 0.015*sqrt(25)
CortSD <- 2.85
pb <- progress_bar$new(total = length(numvec)) #for progress bar

for (j in 1:length(numvec)) {
  pb$tick() #for progress bar
  Sys.sleep(1 / nsim)
  N <- numvec[j] #sample size
  SimulatedData <- data.frame(ID = as.character(seq(1,N,1)), 
                              Temp = as.numeric(rep(c(14,18,22,26,30,34), length = N)), # six temp treatments 
                              sex = NA,
                              Mass = NA) 
  for (i in 1:nsim){
    SimulatedData$sex <- as.character(sample(c("M","F"), N, replace = TRUE)) #random sexes assigned to each temp
    for (k in 1:N) { 
      if (SimulatedData$sex[k] == "M") {
        SimulatedData$Mass[k] <- rgamma(1, shape = 6.8^2/0.4^2, scale = 0.4^2/6.8) #mass for males
      } else {
        SimulatedData$Mass[k] <- rgamma(1, shape = 20^2/2^2, scale = 2^2/20) #mass for females
      }
    }
    MSMR_det <- SimulatedData$Mass^-0.25 * exp(SimulatedData$Temp/10) * (0.2/6.04) 
    Cort_det <- SimulatedData$Mass^-0.25 * exp(SimulatedData$Temp/10) * (7.8/6.04)
    SimulatedData$MSMR <- rgamma(N, shape = (MSMR_det^2)/(MSMRSD^2), scale = (MSMRSD^2)/MSMR_det) #shape = mean^2/variance, scale = variance/mean
    SimulatedData$Cort <- rgamma(N, shape = (Cort_det^2)/(CortSD^2), scale = (CortSD^2)/Cort_det) 
    MSMRModel <- lm(log(MSMR) ~ log(Mass) + Temp, data = SimulatedData)
    CortModel <- lm(log(Cort) ~ log(Mass) + Temp, data = SimulatedData)
    CortMSMRModel <- lm(log(Cort) ~ log(MSMR), data = SimulatedData)
    MSMRMasspval[i] <- coef(summary(MSMRModel))["log(Mass)", "Pr(>|t|)"]
    MSMRTemppval[i] <- coef(summary(MSMRModel))["Temp", "Pr(>|t|)"]
    CortMasspval[i] <- coef(summary(CortModel))["log(Mass)", "Pr(>|t|)"]
    CortTemppval[i] <- coef(summary(CortModel))["Temp", "Pr(>|t|)"]
    CortMSMRpval[i] <- coef(summary(CortMSMRModel))["log(MSMR)", "Pr(>|t|)"]
  }
  MSMRMassPowerVec[j] <- sum(MSMRMasspval < 0.05)/nsim
  MSMRTempPowerVec[j] <- sum(MSMRTemppval < 0.05)/nsim
  CortMassPowerVec[j] <- sum(CortMasspval < 0.05)/nsim
  CortTempPowerVec[j] <- sum(CortTemppval < 0.05)/nsim
  CortMSMRPowerVec[j] <- sum(CortMSMRpval < 0.05)/nsim
}

par(mfrow=c(2,3)) #set up 2x2 plot layout)

#MSMR mass plot
plot(numvec, MSMRMassPowerVec, xlab="Sample size", ylab="Statistical power", main = "MSMR ~ Mass", type="p", pch=16, bty = "l")
lines(lowess(numvec, MSMRMassPowerVec), col="red", lwd = 2) #add a lowess line to the plot
grid() #add gridlines
abline(h=0.8, col="blue") #add a line at 0.8 power 
mtext(text="A", side = 3, adj = 0, font = 2)#add a letter to the plot

#MSMR temp plots
plot(numvec[1:10], MSMRTempPowerVec[1:10], xlab="Sample size", ylab="Statistical power", main = "MSMR ~ Temperature", type="p", pch=16, bty = "l")
lines(lowess(numvec[1:10], MSMRTempPowerVec[1:10], f = 0.1), col="red", lwd = 2) #add a lowess line to the plot
grid() #add gridlines
abline(h=0.8, col="blue") #add a line at 0.8 power
mtext(text="B", side = 3, adj = 0, font = 2)

#Cort mass plot
plot(numvec, CortMassPowerVec, xlab="Sample size", ylab="Statistical power", main = "Cort ~ Mass", type="p", pch=16, bty = "l")
lines(lowess(numvec, CortMassPowerVec), col="red", lwd = 2) #add a lowess line to the plot
grid() #add gridlines
abline(h=0.8, col="blue") #add a line at 0.8 power
mtext(text="C", side = 3, adj = 0, font = 2)

#Cort temp plots
plot(numvec[1:10], CortTempPowerVec[1:10], xlab="Sample size", ylab="Statistical power", main = "Cort ~ Temperature", type="p", pch=16, bty = "l")
lines(lowess(numvec[1:10], CortTempPowerVec[1:10], f = 0.1), col="red", lwd = 2) #add a lowess line to the plot
grid() #add gridlines
abline(h=0.8, col="blue") #add a line at 0.8 power
mtext(text="D", side = 3, adj = 0, font = 2)

#Cort MSMR plot
plot(numvec[1:20], CortMSMRPowerVec[1:20], xlab="Sample size", ylab="Statistical power", main = "Cort ~ MSMR", type="p", pch=16, bty = "l")
lines(lowess(numvec[1:20], CortMSMRPowerVec[1:20], f = 0.2), col="red", lwd = 2) #add a lowess line to the plot
grid() #add gridlines
abline(h=0.8, col="blue") #add a line at 0.8 power
mtext(text="E", side = 3, adj = 0, font = 2)


# dev.off() #close the plot device







