# =========================================================
# Some stuff to check how the collection is coming.
# More or less for IACUC info, this one won't filter out
# frogs that did't work
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================


Data <- read.csv("Data_Spreadsheet.csv")

library(dplyr)
library(tidyr)
library(janitor)

SummaryTable <- Data %>%
  mutate(Temperature = if_else(is.na(Temperature), "NA", as.character(Temperature))) %>%
  group_by(Species, Temperature) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Temperature, values_from = N, values_fill = 0) %>%
  adorn_totals(where = c("row", "col"))

SummaryTable

# Plots -------------------------------------------------------------------

SizeCount <- Data %>%
  group_by(Temperature) %>%
  summarise(n = n())

SppCount <- Data %>%
  group_by(Species) %>%
  summarise(n = n())

ggplot(Data, aes(x = Temperature, y = Weight)) +
  geom_jitter(width = 0.1) +
  theme_classic()

ggplot(Data, aes(x = Species, y = Weight)) +
  geom_jitter(width = 0.1) +
  theme_classic()