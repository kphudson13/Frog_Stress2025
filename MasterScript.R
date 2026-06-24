
library(tidyverse)

Data <- read.csv("Data_Spreadsheet.csv")


# Models ------------------------------------------------------------------

MSMRWeightTemp_Model <- lm(data = Data, log(VCO2) ~ log(Weight) + Temperature)
summary(MSMRWeightTemp_Model)
CortMSMR_Model <- lm(data = Data, log(Cort) ~ log(VCO2/Weight))
summary(CortMSMR_Model)
CortWeightTemp_Model <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature)
summary(CortWeightTemp_Model)


# Plots -------------------------------------------------------------------


(MSMRWeight_Plot <- ggplot(Data, aes(x = log(Weight), y = log(VCO2/Weight))) +
   geom_point(aes(colour = Species)) +
   geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
               slope = coefficients(summary(MSMRWeightTemp_Model))[2,1]) + 
   theme_classic() +
   annotate("text", size = 3.5, x = 2, y = -12,
            label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                     ~x^.(round(coefficients(summary(MSMRWeightTemp_Model))[2,1], 2)),
                                     ~R^2 ==~ .(round(summary(MSMRWeightTemp_Model)$r.squared, 2))))),
            parse = TRUE)
)

(MSMRTemp_Plot <- ggplot(data = Data, aes(x = Temperature, y = log(VCO2/Weight))) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
                slope = coefficients(summary(MSMRWeightTemp_Model))[3,1]) +
    theme_classic() +
    annotate("text", size = 3.5, x = 15, y = -11,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(MSMRWeightTemp_Model))[3,1], 2)),
                                      ~R^2 ==~ .(round(summary(MSMRWeightTemp_Model)$r.squared, 2))))),
             parse = TRUE)
)

(CortMSMR_Plot <- ggplot(data = Data, aes(x=log(VCO2/Weight), y = log(Cort), shape = Notes == "High CV")) +
    geom_point(aes(colour = Species)) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_shape_manual(
      values = c(`TRUE` = 15,   # square
                 `FALSE` = 16)) +# circle
    theme_classic()
)

(CortWeight_Plot <- ggplot(data = Data, aes(x=log(Weight), y = log(Cort))) +
    geom_point(aes(colour = Species)) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic()
)

(CortTemp_Plot <- ggplot(data = Data, aes(x=log(Temperature), y = log(Cort))) +
    geom_point(aes(colour = Species)) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic()
)




