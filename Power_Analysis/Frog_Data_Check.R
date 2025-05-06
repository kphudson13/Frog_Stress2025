
library(tidyverse)

rm(list = ls()) #clear workspace

#one simulated data frame
SimulatedData <- read.csv("Power_Analysis/SimulatedData.csv")

summary(SimulatedData)

# Exploratory Plots -------------------------------------------------------

#histogram of mass for each sex 
par(mfrow = c(1,2))
hist(subset(SimulatedData, sex == "M")$Mass, main = "Male mass", xlab = "Mass")
hist(subset(SimulatedData, sex == "F")$Mass, main = "Female mass", xlab = "Mass")
dev.off()

#histogram of MSMR for each sex
par(mfrow = c(3,2))
hist(subset(SimulatedData, Temp == 14)$MSMR, breaks = seq(0,0.8,by=0.02), main = "14c", xlab = "MSMR")
hist(subset(SimulatedData, Temp == 18)$MSMR, breaks = seq(0,0.8,by=0.02), main = "18c", xlab = "MSMR")
hist(subset(SimulatedData, Temp == 22)$MSMR, breaks = seq(0,0.8,by=0.02), main = "22c", xlab = "MSMR")
hist(subset(SimulatedData, Temp == 26)$MSMR, breaks = seq(0,0.8,by=0.02), main = "26c", xlab = "MSMR")
hist(subset(SimulatedData, Temp == 30)$MSMR, breaks = seq(0,0.8,by=0.02), main = "30c", xlab = "MSMR")
hist(subset(SimulatedData, Temp == 34)$MSMR, breaks = seq(0,0.8,by=0.02), main = "34c", xlab = "MSMR")
dev.off()

#histogram of cort for each sex 
par(mfrow = c(3,2))
hist(subset(SimulatedData, Temp == 14)$Cort, breaks = seq(0,28,by=2), main = "14c", xlab = "Cort")
hist(subset(SimulatedData, Temp == 18)$Cort, breaks = seq(0,28,by=2), main = "18c", xlab = "Cort")
hist(subset(SimulatedData, Temp == 22)$Cort, breaks = seq(0,28,by=2), main = "22c", xlab = "Cort")
hist(subset(SimulatedData, Temp == 26)$Cort, breaks = seq(0,28,by=2), main = "26c", xlab = "Cort")
hist(subset(SimulatedData, Temp == 30)$Cort, breaks = seq(0,28,by=2), main = "30c", xlab = "Cort")
hist(subset(SimulatedData, Temp == 34)$Cort, breaks = seq(0,28,by=2), main = "34c", xlab = "Cort")
dev.off()

#MSMR ~ Mass
ggplot(SimulatedData, aes(x = log(Mass), y = log(MSMR), color = as.factor(Temp))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() 

#MSMR ~ Temp
ggplot(SimulatedData, aes(x = Temp, y = log(MSMR))) +
  geom_jitter(width = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

#Cort ~ Mass
ggplot(SimulatedData, aes(x = log(Mass), y = log(Cort), color = as.factor(Temp))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

#Cort ~ Temp
ggplot(SimulatedData, aes(x = Temp, y = Cort)) +
  geom_jitter(width = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

#Cort ~ MSMR log
ggplot(SimulatedData, aes(x = log(MSMR), y = log(Cort))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()

#Cort ~ MSMR not logged 
ggplot(SimulatedData, aes(x = MSMR, y = Cort)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic()








