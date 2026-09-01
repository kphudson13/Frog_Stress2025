# =========================================================
# Flow-through respirometry CO2 extraction with
# asymptotic fitting for each animal segment
#
# Author: Kyle Hudson
# Circa 2026
# Live laugh love
# =========================================================

# Setup -------------------------------------------------------------------

# devtools::install_github("daniel1noble/metabR",
#                          dependencies = TRUE,
#                          force = TRUE)

library(metabR)
library(tidyverse)
library(minpack.lm)

rm(list = ls()) # Clear environment

dir.create("Data_Extraction/Raw_Figures",
           showWarnings = FALSE,
           recursive = TRUE) # Create directory for figures

files <- list.files(path = "Data_Extraction",
                    pattern = "\\.exp$",
                    full.names = TRUE) # List all .exp files 

filetemp <- read.exp("Data_Extraction/14cNS-Mar13_0001.exp")

# Experimental settings ---------------------------------------------------

flow <- 50                 # flow rate (mL/min)
discard <- 40              # seconds discarded after switching
cycle_length <- 361 + 601  # total cycle duration (sec)

# FUNCTION: Fit asymptotic exponential model ------------------------------
#
# Model:
#
# CO2(t) = A + (y0 - A) * exp(-k*t)
#
# A  = asymptote (steady-state equilibrium CO2)
# y0 = starting CO2
# k  = wash-in rate constant

fit_asymptote <- function(df){
  
  df <- df %>%
    mutate(t = row_number() - 1)
  
  # Starting parameter estimates
  start_A  <- mean(tail(df$CO2, 30), na.rm = TRUE) # last 30 points
  start_y0 <- first(df$CO2) # first point after discard
  start_k  <- 0.01 
  
  # Nonlinear fit
  fit <- try(nlsLM(
    CO2 ~ A + (y0 - A) * exp(-k * t),
    data = df,
    start = list(A = start_A,
                 y0 = start_y0,
                 k = start_k),
    control = nls.lm.control(maxiter = 200)),
    silent = TRUE)
  
  if(inherits(fit, "try-error")){
    return(NULL)
  } # If fit fails
  
  df$fitted <- predict(fit) # Predicted fitted values
  
  list(asymptote = coef(fit)["A"],
       fit = fit,
       data = df) # Return results
}

# FUNCTION: Process one .exp file -----------------------------------------

process_exp <- function(file_path){
  
  base_name <- tools::file_path_sans_ext(basename(file_path)) # save the title but not full path

  meta <- str_match(base_name, "^(\\d+c)(NS|S)-([A-Za-z0-9]+)_") # save title as metadata
  
  temp    <- str_remove(meta[,2], "c") # remove the c
  stress  <- meta[,3]
  date    <- meta[,4]
  run_num <- str_extract(base_name, "\\d+$") # save individual parts from title
  
  fig_name <- paste(temp, stress, date, run_num, sep = "_") # name figure based off metadata
  
  # Read Sable Systems file
  sableDat <- read.exp(file_path) %>%
    as.data.frame(.) %>%
    slice(-1) %>% # Cut out first row
    mutate(CO2_Raw = CO2, 
           CO2 = CO2*(BP-WVP)/BP, # correct for WVP
           time = row_number(),
           cycle = floor((time - 1) / cycle_length),
           time_in_cycle = (time - 1) %% cycle_length,
           phase = ifelse(time_in_cycle < 361,
                          "control",
                          "animal"), # second part is the actual animal
           valid = (phase == "control" & time_in_cycle >= discard) |
             (phase == "animal" & time_in_cycle >= (361 + discard)))
  
  # PLOT RAW CO2 TRACE + FITTED ASYMPTOTES
  fig_file <- file.path("Data_Extraction/Raw_Figures", 
                        paste0(fig_name, "_RawCO2.png"))
  
  # Open PNG device
  grDevices::png(
    filename = fig_file,
    width = 12,
    height = 5,
    units = "in",
    res = 300)
  
  # Always close graphics device if function exits unexpectedly
  on.exit(try(grDevices::dev.off(),
              silent = TRUE),
          add = TRUE)
  
  # Base plot
  plot(sableDat$CO2, type = "l", col = "black", lwd = 1, 
       xlab = "Time (sec)", ylab = "CO2",
       ylim = c(min(-0.01, min(sableDat$CO2, na.rm = TRUE)),
                max(sableDat$CO2, na.rm = TRUE)),
       main = paste(fig_name, "\nRaw CO2 Trace with Asymptotic Fits"))
  
  # Add vertical lines for chamber switching
  num_cycles <- max(sableDat$cycle, na.rm = TRUE)
  
  for(i in 0:num_cycles){
    
    control_start <- which(sableDat$cycle == i &
                             sableDat$time_in_cycle == 0)[1] # Find control start
    
    if(!is.na(control_start)){
      abline(
        v = control_start,
        col = "forestgreen",
        lty = 2,
        lwd = 1.5)
    } # Plot control start
    
    animal_start <- which( sableDat$cycle == i &
                             sableDat$time_in_cycle == 361)[1] # Find animal starts
    
    if(!is.na(animal_start)){
      
      abline(
        v = animal_start,
        col = "purple",
        lty = 2,
        lwd = 1.5)
    } # plot animal starts
  } # close chamber switch line loop
  
  
  # Fit each cycle and plot asymptotes
  for(i in 0:num_cycles){
    
    # CONTROL / BLANK
    # Use average of final 30 seconds
    control_df <- sableDat %>%
      filter(cycle == i, phase == "control", valid) # use only control cycles
    
    if(nrow(control_df) >= 30){
      control_asym <- mean(
        tail(control_df$CO2, 30),
        na.rm = TRUE) # Average CO2 from the final 30 seconds
      
      # Locate control indices in full trace
      control_idx <- which(sableDat$cycle == i &
                             sableDat$phase == "control" &
                             sableDat$valid)
      
      lines(x = c(max(control_idx) - 180,
                  max(control_idx)),
            y = c(control_asym,
                  control_asym),
            col = "red",
            lty = 2,
            lwd = 2 ) # Plot control average
    }
    
    # ANIMAL ASYMPTOTE
    animal_df <- sableDat %>%
      filter(cycle == i, phase == "animal", valid)
    
    if(nrow(animal_df) < 180){
      next
    } # Skip tiny datasets
    
    fit_result <- fit_asymptote(animal_df) # Fit animal asymptote
    
    if(is.null(fit_result)){
      next
    } # Skip failed animal fits
    
    fitted_df <- fit_result$data
    
    asym <- fit_result$asymptote
    
    # Locate animal indices in full trace
    idx <- which(
      sableDat$cycle == i &
        sableDat$phase == "animal" &
        sableDat$valid
    )
    
    lines(idx,
          fitted_df$fitted,
          col = "blue",
          lwd = 2) # Add fitted exponential curve
    
    lines(x = c(min(idx), max(idx)),
          y = c(asym, asym),
          col = "red",
          lty = 2,
          lwd = 2) # Add animal asymptote
    
    text(x = max(idx),
         y = asym,
         labels = paste0( "A = ", round(asym, 5)),
         pos = 3,
         cex = 0.7,
         col = "red") # Label animal asymptote
  }
  
  grDevices::dev.off()  # Explicitly close PNG device
  
  # CALCULATE CONTROL VALUES
  # Average of final 30 seconds of each control chamber
  control_vals <- sableDat %>%
    filter(phase == "control", valid) %>%
    group_by(cycle) %>%
    summarise(control_co2 = mean(tail(CO2, 30),na.rm = TRUE), 
              .groups = "drop")
  
  # CALCULATE ANIMAL ASYMPTOTES
  animal_summary <- sableDat %>%
    filter(phase == "animal", valid) %>%
    group_by(cycle) %>%
    group_modify(~{
      fit_result <- fit_asymptote(.x)
      tibble(total_co2 = ifelse(is.null(fit_result), NA, fit_result$asymptote))}) %>%
    ungroup()
  
  # COMBINE SUMMARIES
  co2_summary <- animal_summary %>%
    left_join(
      control_vals,
      by = "cycle") %>%
    mutate(delta_co2 = total_co2 - control_co2, # Baseline-corrected CO2
           VCO2_ml_min = flow * delta_co2 / 100, # Convert to mL/min (sable unit is percent)
           file = base_name, # Metadata
           temp = temp,
           stress = stress,
           date = date,
           line = case_when(cycle == 0 ~ 2, # match order to lines from machine
                            cycle == 1 ~ 3, # this only matters for another column in the printout
                            cycle == 2 ~ 5,
                            cycle == 3 ~ 6,
                            cycle == 4 ~ 7,
                            cycle == 5 ~ 8, 
                            TRUE       ~ NA_real_ )) # safety net – should never hit
  
  return(co2_summary) # Return summary
}

all_co2 <- map_dfr(files, process_exp) # PROCESS ALL FILES

head(all_co2)

write.csv(
  all_co2,
  "Data_Extraction/CO2_Asymptote_Summary.csv",
  row.names = FALSE
) # EXPORT


