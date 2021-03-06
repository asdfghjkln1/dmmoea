#!bin/usr/env Rstudio
run_normalization <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 2){
    print(paste0("Not enough parameters (", argnum, "/1)"))
    return(-1)
  }
  path <- args[1] 
  test <- args[2]
  setwd(path)
  source("dmmoea_functions.R")
  library(ggplot2)
  #source("dmmoea_parameters.R")
  #source("dmmoea_libraries.R")
  #source("dmmoea_distances.R")
  #source("dmmoea_irace_conf.R")
  
  test.path <<- file.path(path, test) #"Tests", paste0("tuning_", test))
  
  if(!file.exists(file.path(test.path, "limits.csv"))){
    print("limits.csv not exists")
    if(length(list.dirs(test.path) > 0)){
      print("Updating normalizaton limits...")
      get_normalization_limits(test.path) 
      print("Ready.")  
    }
  }else{
    print("limits.csv exists")
    limits <- read.table(file.path(test.path, "limits.csv"), sep=",", row.names = NULL, header = TRUE)
    print(limits)
    if(length(list.dirs(test.path) > 0)){
      print("Updating anyway...")
      get_normalization_limits(test.path)
      print("Ready.")
    }
  }
  
  return(0)
}

run_normalization()
