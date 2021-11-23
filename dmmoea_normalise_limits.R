#!bin/usr/env Rstudio
run_normalization <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 1){
    print(paste0("Not enough parameters (", argnum, "/1)"))
    return(-1)
  }
  path <- args[1] 
  test <- args[2]
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  test.path <<- file.path(path, "Tests", paste0("tuning_", test))
  
  if(!file.exists(file.path(test.path, "limits.csv"))){
    print("Limits not exists")
    if(length(list.dirs(test.path) > 0)){
      print("Updating normalizaton limits...")
      get_normalization_limits(test.path) 
      print("Ready.")  
    }
  }else{
    print("Limits exists")
    limits <- read.table(file.path(test.path, "limits.csv"), sep=",", row.names = NULL, header = TRUE)
    print(limits)
  }
  
  return(0)
}

run_normalization()
