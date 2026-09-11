# =========================================================
# Core statistics and plots for the 2026 study comparing
# metabolic rate and corticosterone in frogs
# 
# This script should be ran from control script, if not 
# you'll need to manually load in data
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================


# Raw Models ------------------------------------------------------------------

MSMRWeightTemp_Unfiltered <- lm(data = Whole_data, log(MSMR) ~ log(Weight) + Temperature + Species) 
CortMSMR_Unfiltered <- lm(data = Data, log(Cort) ~ log(MSMR) + Species) 
CortWeightTemp_Unfiltered <- lm(data = Data, log(Cort) ~ log(Weight) + Temperature + Species)


# Run Cook's distance filtering -----------------------------------------------

CDist_fun(MSMRWeightTemp_Unfiltered,
          log(MSMR) ~ log(Weight) + Temperature + Species,
          Whole_data)

summary(MSMRWeightTemp_Model)

CDist_fun(CortWeightTemp_Unfiltered,
          log(Cort) ~ log(Weight) + Temperature + Species,
          Data)

summary(CortWeightTemp_Model)

CDist_fun(CortMSMR_Unfiltered,
          log(Cort) ~ log(MSMR) + Species,
          Data)

summary(CortMSMR_Model)

# Plots -------------------------------------------------------------------

theme1 <- theme(legend.background = element_blank(), # make clear 
                legend.title = element_blank(), # remove title
                legend.key.height = unit(0.6, "lines"), # fit legend lines closer together
                legend.position=c(.8, 0.1)
) # set standard theme for plots

# normalized points for plotting
MSMRWeightTemp_Data$MSMR_NormTemp <- log(MSMRWeightTemp_Data$MSMR) - MSMRWeightTemp_Model[["coefficients"]][["Temperature"]]*MSMRWeightTemp_Data$Temperature

(MSMRWeight_Plot <- ggplot(MSMRWeightTemp_Data, aes(x = log(Weight), y = MSMR_NormTemp)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
                slope = coefficients(summary(MSMRWeightTemp_Model))[2,1]) + # pull lines from model
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 1, y = -6,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(MSMRWeightTemp_Model))[2,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(MSMRWeightTemp_Model)$partial.rsq[1], 2))))),
             parse = TRUE) +
    labs(x = "Weight(ln(g))",
         y = expression("Normalized MSMR (ln(ml CO2 min"^-1*" g"^-1*"))")) +
    scale_x_continuous(limits = c(-1, 4)) +
    scale_y_continuous(limits = c(-12, -4))
) 

ggsave(filename = paste("Figures/", directory, "/MSMRWeight_Plot.png", sep = ""), 
       width=90, height=90, units="mm") #save a picture

# normalized points for plotting
MSMRWeightTemp_Data$MSMR_NormWeight <- log(MSMRWeightTemp_Data$MSMR) - MSMRWeightTemp_Model[["coefficients"]][["log(Weight)"]]*log(MSMRWeightTemp_Data$Weight)

(MSMRTemp_Plot <- ggplot(MSMRWeightTemp_Data, aes(x = Temperature, y = MSMR_NormWeight)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(MSMRWeightTemp_Model))[1,1],
                slope = coefficients(summary(MSMRWeightTemp_Model))[3,1]) +
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 20, y = -4,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(MSMRWeightTemp_Model))[1,1], 2))
                                      ~e^{.(round(coefficients(summary(MSMRWeightTemp_Model))[3,1], 2))*x},
                                      Partial ~R^2 ==~ .(round(rsq.partial(MSMRWeightTemp_Model)$partial.rsq[2], 2))))),
             parse = TRUE) +
    labs(x = "Temperature (c)",
         y = expression("Normalized MSMR (ln(ml CO2 min"^-1*" g"^-1*"))")) +
    scale_x_continuous(limits = c(12, 40)) +
    scale_y_continuous(limits = c(-10, -3))
)

ggsave(filename = paste("Figures/", directory, "/MSMRTemp_Plot.png", sep = ""), 
       width=90, height=90, units="mm") #save a picture

# normalized points for plotting
CortWeightTemp_Data$Cort_NormTemp <- log(CortWeightTemp_Data$Cort) - CortWeightTemp_Model[["coefficients"]][["Temperature"]]*CortWeightTemp_Data$Temperature

(CortWeight_Plot <- ggplot(CortWeightTemp_Data, aes(x=log(Weight), y = Cort_NormTemp)) +
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
    labs(x = "Weight(ln(g))",
         y = "Normalized Cort (ln(ng/ml))") +
    scale_x_continuous(limits = c(-1, 4)) +
    scale_y_continuous(limits = c(-2, 6))
)

ggsave(filename = paste("Figures/", directory, "/CortWeight_Plot.png", sep = ""), 
       width=90, height=90, units="mm") #save a picture

# normalized points for plotting
CortWeightTemp_Data$Cort_NormWeight <- log(CortWeightTemp_Data$Cort) - CortWeightTemp_Model[["coefficients"]][["log(Weight)"]]*log(CortWeightTemp_Data$Weight)

(CortTemp_Plot <- ggplot(CortWeightTemp_Data, aes(x=Temperature, y = Cort_NormWeight)) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(CortWeightTemp_Model))[1,1],
                slope = coefficients(summary(CortWeightTemp_Model))[3,1]) +
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = 20, y = 5.5,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(CortWeightTemp_Model))[1,1], 2))
                                      ~e^{.(round(coefficients(summary(CortWeightTemp_Model))[3,1], 2))*x},
                                      Partial ~R^2 ==~ .(round(rsq.partial(CortWeightTemp_Model)$partial.rsq[2], 2))))),
             parse = TRUE) +
    labs(x = "Temperature (c)",
         y = "Normalized Cort (ln(ng/ml))") +
    scale_x_continuous(limits = c(12, 40)) +
    scale_y_continuous(limits = c(-2, 7))
)

ggsave(filename = paste("Figures/", directory, "/CortTemp_Plot.png", sep = ""), 
       width=90, height=90, units="mm") #save a picture

(CortMSMR_Plot <- ggplot(CortMSMR_Data, aes(x=log(MSMR), y = log(Cort))) +
    geom_point(aes(colour = Species)) +
    geom_abline(intercept = coefficients(summary(CortMSMR_Model))[1,1],
                slope = coefficients(summary(CortMSMR_Model))[2,1]) + 
    theme_classic() +
    theme1 +
    annotate("text", size = 3.5, x = -9, y = 6,
             label = list(bquote(atop(y==~ .(round(coefficients(summary(CortMSMR_Model))[1,1], 2))
                                      ~x^.(round(coefficients(summary(CortMSMR_Model))[2,1], 2)),
                                      Partial ~R^2 ==~ .(round(rsq.partial(CortMSMR_Model)$partial.rsq[1], 2))))),
             parse = TRUE) +
    labs(x = expression("MSMR (ln(ml CO2 min"^-1*" g"^-1*"))"),
         y = "Cort (ln(ng/ml))") +
    scale_x_continuous(limits = c(-10, -4)) +
    scale_y_continuous(limits = c(-1, 7))
)

ggsave(filename = paste("Figures/", directory, "/CortMSMR_Plot.png", sep = ""), 
       width=90, height=90, units="mm") #save a picture

# Stats Table -------------------------------------------------------------

# Extract CIs beforehand
ci_list <- rbind(confint(MSMRWeightTemp_Model),
                 confint(CortWeightTemp_Model),
                 confint(CortMSMR_Model)) %>%
  as.data.frame() %>%
  slice(-c(1,4,5,6,7,10,11,12,13,15,16,17))

# Extract n for each model, only shown in first row per model
n_list <- c(nobs(MSMRWeightTemp_Model), NA,
            nobs(CortWeightTemp_Model), NA,
            nobs(CortMSMR_Model))

Stats_Tab <- rbind(coefficients(summary(MSMRWeightTemp_Model)),
                   coefficients(summary(CortWeightTemp_Model)),
                   coefficients(summary(CortMSMR_Model))) %>%
  as.data.frame(.) %>%
  slice(-c(1,4,5,6,7,10,11,12,13,15,16,17)) %>%
  select(., -"t value") %>%
  cbind(., partial_r2 = c(rsq.partial(MSMRWeightTemp_Model)$partial.rsq[1:2],
                          rsq.partial(CortWeightTemp_Model)$partial.rsq[1:2],
                          rsq.partial(CortMSMR_Model)$partial.rsq[1]),
        model_r2 = c(summary(MSMRWeightTemp_Model)$r.squared, NA,
                     summary(CortWeightTemp_Model)$r.squared, NA,
                     summary(CortMSMR_Model)$r.squared)) %>%
  mutate(across(c(1:5), \(x) round(x, digits = 2))) %>%
  `colnames<-`(c("Est. [95% CI]", "S.E.", "p val.", "Partial R2", "Model R2")) %>%
  mutate(`Est. [95% CI]` = paste0(`Est. [95% CI]`, " [", round(ci_list[,1], 2), ", ", round(ci_list[,2], 2), "]"),
         `p val.` = ifelse(`p val.` < 0.001, "< 0.001", `p val.`),
         `Model R2` = ifelse(is.na(`Model R2`), "", as.character(`Model R2`)),
         n = ifelse(is.na(n_list), "", as.character(n_list))) %>%
  select(`Est. [95% CI]`, n, `S.E.`, `p val.`, `Partial R2`, `Model R2`) %>%
  `rownames<-`(c("MSMR ~ Weight", "MSMR ~ Temp", "Cort ~ Weight", "Cort ~ Temp", "Cort ~ MSMR"))


tt1 <- ttheme_minimal(rowhead=list(fg_params=list(fontface = "bold"))
) #set the theme for the table png)

write.csv(Stats_Tab, file = paste("Figures/", directory, "/StatsTab.csv", sep = ""), row.names = TRUE)

#export stats table
png(paste("Figures/", directory, "/StatsTab.png", sep = ""),
    height = 180*nrow(Stats_Tab),
    width = 500*ncol(Stats_Tab),
    res = 300)
grid.newpage()
grid.table(Stats_Tab, theme = tt1)
# grid.text("Stats Table", x = 0.2, y = 0.9, gp = gpar(fontface = "bold"))
dev.off()






