# =========================================================
# Some stuff to check how the collection is coming.
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================

SummaryTable <- Data %>%
  mutate(Temperature = if_else(is.na(Temperature), "NA", as.character(Temperature))) %>%
  group_by(Species, Temperature) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Temperature, values_from = N, values_fill = 0) %>%
  adorn_totals(where = c("row", "col"))

SummaryTable

# Plots -------------------------------------------------------------------


# Calculate counts per Temperature (regardless of species)
counts_temp_cort <- CortWeightTemp_Data %>%
  group_by(Temperature) %>%
  summarise(n = n(), y_pos = max(log(Weight)), .groups = "drop")

# Calculate counts per Species (for legend)
counts_species_cort <- CortWeightTemp_Data %>%
  group_by(Species) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Species, " (n = ", n, ")"))

# Create named vector for legend labels
legend_labels_Cort <- setNames(counts_species_cort$label, counts_species_cort$Species)

(WeightTemp_PlotCORT <- ggplot(CortWeightTemp_Data, aes(x = Temperature, y = log(Weight), colour = Species)) +
    geom_jitter(width = 0.1) +
    geom_text(data = counts_temp_cort,
              aes(x = Temperature, y = y_pos + 0.2, label = paste0("n = ", n)),
              colour = "black",
              size = 3) +
    theme_classic() +
    theme1 +
    theme(legend.position=c(0.7, 0.1)) +
    labs(x = "Temperature (°C)",
         y = " ") +
    scale_x_continuous(limits = c(9, 40), breaks = c(14, 24, 34)) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
    scale_colour_discrete(labels = legend_labels_Cort))



# Calculate counts per Temperature (regardless of species)
counts_temp_MSMR <- MSMRWeightTemp_Data %>%
  group_by(Temperature) %>%
  summarise(n = n(), y_pos = max(log(Weight)), .groups = "drop")

# Calculate counts per Species (for legend)
counts_species_MSMR <- MSMRWeightTemp_Data %>%
  group_by(Species) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Species, " (n = ", n, ")"))

# Create named vector for legend labels
legend_labels_MSMR <- setNames(counts_species_MSMR$label, counts_species_MSMR$Species)

(WeightTemp_PlotMSMR <- ggplot(MSMRWeightTemp_Data, aes(x = Temperature, y = log(Weight), colour = Species)) +
    geom_jitter(width = 0.15) +
    geom_text(data = counts_temp_MSMR,
              aes(x = Temperature, y = y_pos + 0.2, label = paste0("n = ", n)),
              colour = "black",
              size = 3) +
    theme_classic() +
    theme1 +
    theme(legend.position=c(0.7, 0.1)) +
    labs(x = "Temperature (°C)",
         y = "Weight(ln(g))") +
    scale_x_continuous(limits = c(9, 40), breaks = c(14, 24, 34)) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
    scale_colour_discrete(labels = legend_labels_MSMR))


save_plot(plot_grid(WeightTemp_PlotMSMR, WeightTemp_PlotCORT, 
                    labels = c("A", "B"), label_fontface = "plain"),
          filename = paste("Figures/", directory, "/WeightTemp_Plot.png", sep = ""),
          base_height = 3.3, base_width = 6.6)











