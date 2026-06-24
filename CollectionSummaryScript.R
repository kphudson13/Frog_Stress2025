


Data <- read.csv("Data_Spreadsheet.csv")

library(dplyr)
library(tidyr)
library(janitor)

SummaryTable <- Data %>%
  filter(!is.na(VCO2)) %>% 
  group_by(Species, Temperature) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Temperature,
    values_from = N,
    values_fill = 0
  ) %>%
  adorn_totals(where = c("row", "col"))

SummaryTable

