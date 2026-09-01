
# Directory function ------------------------------------------------------

#make a figure directory for initial pull users
CreateDR <- function(DR) {
  if (!file.exists("Figures")) {
    dir.create("Figures")
  }
  if (!file.exists(paste("Figures/", DR, sep = ""))) {
    dir.create(paste("Figures/", DR, sep = ""))
  }
}


# #Pull model 
# 
# PullModel <- function(mod1) {
#   if(file.exists(paste(DR, mod1, sep = ""))) {
#     load(paste(DR, mod1, sep = ""))
#   } else {
#     stop("Model not found. Run 'MasterScript' first")
#   }
# }
# 
# 
# }


