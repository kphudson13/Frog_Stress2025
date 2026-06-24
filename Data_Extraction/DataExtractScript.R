# =========================================================
# Flow-through respirometry CO2 extraction with
# asymptotic fitting for each animal segment
#
# Author: Kyle Hudson
# =========================================================


# Setup -------------------------------------------------------------------

# devtools::install_github("daniel1noble/metabR",
#                          dependencies = TRUE,
#                          force = TRUE)

library(metabR)
library(tidyverse)
library(minpack.lm)

dir.create("Data_Extraction/Raw_Figures",
           showWarnings = FALSE) # create directory for plots


# Experimental settings
flow <- 50                 # flow rate (mL/min)
discard <- 20              # seconds discarded after switching
cycle_length <- 361 + 601  # total cycle duration (sec)

files <- list.files(
  path = "Data_Extraction",
  pattern = "\\.exp$",
  full.names = TRUE
) # List all .exp files in directory

# FUNCTION: Fit asymptotic exponential model
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
    mutate(t = row_number() - 1) # create relative time starting at zero
  
  # starting parameter estimates
  start_A  <- mean(tail(df$CO2, 30), na.rm = TRUE)
  start_y0 <- first(df$CO2)
  start_k  <- 0.01
  
  # nonlinear fit
  fit <- try(
    
    nlsLM(
      CO2 ~ A + (y0 - A) * exp(-k * t),
      data = df,
      start = list(
        A = start_A,
        y0 = start_y0,
        k = start_k
      ),
      control = nls.lm.control(
        maxiter = 200
      )
    ),
    silent = TRUE
  )
  
  if(inherits(fit, "try-error")){
    return(NULL)
  } # if fit fails
  
  df$fitted <- predict(fit) # predicted fitted values
  
  # return results
  list(
    asymptote = coef(fit)["A"],
    fit = fit,
    data = df
  )
}

# FUNCTION: Process one .exp file
process_exp <- function(file_path){
  
  # File metadata
  file_name <- basename(file_path)
  base_name <- tools::file_path_sans_ext(file_name)
  # example filename: 24cNS-Feb18_0001
  meta <- str_match(
    base_name, "^(\\d+c)(NS|S)-([A-Za-z0-9]+)_")
  
  temp   <- meta[,2]
  stress <- meta[,3]
  date   <- meta[,4]
  run_num <- str_extract(base_name, "\\d+$")
  
  fig_name <- paste(temp, stress, date, run_num, sep = "_")
  sableDat <- read.exp(file_path) # Read Sable Systems file
  
  # Convert to dataframe and define cycles
  sableDat <- as.data.frame(sableDat) %>%
    slice(-1) %>%
    mutate(
      time = row_number(), # overall time index
      cycle = floor((time - 1) / cycle_length), # cycle number
      time_in_cycle = (time - 1) %% cycle_length, # time within cycle
      phase = ifelse(
        time_in_cycle < 361,
        "control",
        "animal"), # define chamber phase
      valid =
        (phase == "control" &
           time_in_cycle >= discard) |
        (phase == "animal" &
           time_in_cycle >= (361 + discard)) # remove switching artifact period
    )
  
  # PLOT RAW CO2 TRACE + FITTED ASYMPTOTES
  fig_file <- file.path(
    "Data_Extraction/Raw_Figures",
    paste0(fig_name, "_RawCO2.png"))
  
  png(
    filename = fig_file,
    width = 12,
    height = 5,
    units = "in",
    res = 300)
  
  # Raw trace
  plot(
    sableDat$CO2,
    type = "l",
    col = "black",
    lwd = 1,
    
    xlab = "Time (sec)",
    ylab = "CO2",
    
    main = paste(
      fig_name,
      "\nRaw CO2 Trace with Asymptotic Fits"
    )
  )
  
  # Add vertical lines for chamber switching
  num_cycles <- max(sableDat$cycle)
  
  for(i in 0:num_cycles){
    control_start <- which(
      sableDat$cycle == i &
        sableDat$time_in_cycle == 0
    )[1]
    
    abline(
      v = control_start,
      col = "forestgreen",
      lty = 2,
      lwd = 1.5) # line for control starts 
    
    animal_start <- which(
      sableDat$cycle == i &
        sableDat$time_in_cycle == 361
    )[1]
    
    abline(
      v = animal_start,
      col = "purple",
      lty = 2,
      lwd = 1.5
    ) #line for animal chamber starts
  }
  
  # =======================================================
  # Fit each animal cycle
  # =======================================================
  
  for(i in 0:num_cycles){
    
    # -----------------------------------------------------
    # Extract animal phase
    # -----------------------------------------------------
    
    animal_df <- sableDat %>%
      
      filter(
        cycle == i,
        phase == "animal",
        valid
      )
    
    if(nrow(animal_df) < 20) next # skip tiny datasets
    
    # -----------------------------------------------------
    # Fit asymptotic model
    # -----------------------------------------------------
    
    fit_result <- fit_asymptote(animal_df)
    
    # skip failed fits
    if(is.null(fit_result)) next
    
    fitted_df <- fit_result$data
    
    asym <- fit_result$asymptote
    
    # -----------------------------------------------------
    # Locate indices in full trace
    # -----------------------------------------------------
    
    idx <- which(
      
      sableDat$cycle == i &
        sableDat$phase == "animal" &
        sableDat$valid
    )
    
    # -----------------------------------------------------
    # Add fitted curve
    # -----------------------------------------------------
    
    lines(
      idx,
      fitted_df$fitted,
      
      col = "blue",
      lwd = 2
    )
    
    # -----------------------------------------------------
    # Add asymptote line ONLY over animal segment
    # -----------------------------------------------------
    
    lines(
      x = c(min(idx), max(idx)),
      y = c(asym, asym),
      
      col = "red",
      lty = 2,
      lwd = 2
    )
    
    # -----------------------------------------------------
    # Label asymptote
    # -----------------------------------------------------
    
    text(
      x = max(idx),
      y = asym,
      
      labels = paste0(
        "A = ",
        round(asym, 5)
      ),
      
      pos = 3,
      cex = 0.7,
      col = "red"
    )
  }
  
  dev.off()
  
  # =======================================================
  # Calculate CONTROL asymptotes
  # =======================================================
  
  control_vals <- sableDat %>%
    
    filter(
      phase == "control",
      valid
    ) %>%
    
    group_by(cycle) %>%
    
    group_modify(~{
      
      fit_result <- fit_asymptote(.x)
      
      tibble(
        
        control_co2 =
          
          ifelse(
            is.null(fit_result),
            NA,
            fit_result$asymptote
          )
      )
    }) %>%
    
    ungroup()
  
  # =======================================================
  # Calculate ANIMAL asymptotes
  # =======================================================
  
  animal_summary <- sableDat %>%
    
    filter(
      phase == "animal",
      valid
    ) %>%
    
    group_by(cycle) %>%
    
    group_modify(~{
      
      fit_result <- fit_asymptote(.x)
      
      tibble(
        
        total_co2 =
          
          ifelse(
            is.null(fit_result),
            NA,
            fit_result$asymptote
          )
      )
    }) %>%
    
    ungroup()
  
  # =======================================================
  # Combine summaries
  # =======================================================
  
  co2_summary <- animal_summary %>%
    
    left_join(
      control_vals,
      by = "cycle"
    ) %>%
    
    mutate(
      
      # baseline-corrected CO2
      delta_co2 =
        total_co2 - control_co2,
      
      # convert to mL/min
      #
      # assumes CO2 is in ppm
      VCO2_ml_min =
        (flow * delta_co2 / 1e6) * 60,
      
      # metadata
      file = base_name,
      temp = temp,
      stress = stress,
      date = date,
      
      run = row_number()
    )
  
  return(co2_summary)
}

# =========================================================
# Process ALL files
# =========================================================

all_co2 <- map_dfr(
  files,
  process_exp
)

# =========================================================
# Final output
# =========================================================

print(all_co2)

# optional export
write.csv(
  all_co2,
  "Data_Extraction/CO2_Asymptote_Summary.csv",
  row.names = FALSE
)

