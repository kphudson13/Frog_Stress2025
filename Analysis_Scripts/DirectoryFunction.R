


#make a figure directory for initial pull users

CreateDR <- function(DR) {
  if (!file.exists("Figures")) {
    dir.create("Figures")
  }
  if (!file.exists(paste("Figures/", DR, sep = ""))) {
    dir.create(paste("Figures/", DR, sep = ""))
  }
}




