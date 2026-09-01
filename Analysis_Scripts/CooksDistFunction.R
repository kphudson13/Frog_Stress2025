
# Cook's Distance Function ---------------------------------------------------

CDist_fun <- function(mod, formula, data) {
  CD <- cooks.distance(mod)
  keep <- CD <= 4 / nobs(mod) 
  model_rows <- as.integer(names(residuals(mod)))
  model_data <- data[model_rows, , drop = FALSE]
  filtered_data <- model_data[keep, , drop = FALSE]
  
  model_name <- deparse(substitute(mod))
  obj_name <- sub("_Unfiltered$", "", model_name)
  data_name <- paste0(obj_name, "_Data")
  
  assign(data_name, filtered_data, envir = .GlobalEnv) # Save filtered data
  
  # do.call forces the data argument to resolve to the actual named global object
  filtered_model <- do.call(lm, list(formula = formula, 
                                     data = as.name(data_name)),
                            envir = .GlobalEnv)
  
  assign(paste0(obj_name, "_Model"), filtered_model, envir = .GlobalEnv)
}