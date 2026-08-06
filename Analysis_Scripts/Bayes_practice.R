

library(tidyverse)
library(rjags)

data <- read.csv("Power_Analysis/SimulatedData.csv")

plot(data)
N=nrow(data)

model_string <- "model{
# Prior for mass
beta1 ~ dnorm(-0.25,0.0001)
# prior for temp
beta2 ~ dnorm(0.1,0.0001)
#Prior for alpha
alpha ~ dnorm(0,0.0001)
# Likelihood
for(i in 1:n){
Y[i] ~ dnorm(mu[i],inv.var)
mu[i] <- alpha+(mass[i]^beta1)*exp(beta2*temp[i])*(0.2/6.04)
}
# Prior for the inverse variance (precision)
inv.var ~ dunif(0, 100)
sigma <- 1/sqrt(inv.var)
}"



model <- jags.model(textConnection(model_string),
                    data = list(Y = data$MSMR, mass=data$Mass, temp = data$Temp, n=N))
update(model, 10000, progress.bar="none") # Burnin for 10000 samples

samp <- coda.samples(model,
                     variable.names=c("alpha","beta1","beta2","sigma"),
                     n.iter=20000, progress.bar="none")
summary(samp)
plot(samp)
