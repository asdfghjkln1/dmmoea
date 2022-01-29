calculate_hypervolume_manual <- function(){
  #library(eaf)
  source("dmmoea_functions.R")
  source("dmmoea_libraries.R")
  args <- commandArgs(trailingOnly = TRUE)
  
  path <- args[1]
  results.path <- args[2]

  setwd(path)
  source("dmmoea_functions.R")
  #source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  #source("dmmoea_irace_conf.R")
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  #path <- "~/dmmoea/Tests/runs/BallHall"
  
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  print(algorithms)
  algorithms <- c("dmnsga2", "dnsga2", "nsga2")

  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    exp.path <- file.path(path, algorithm, dataset)
    limits <- read.table(file.path(exp.path, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
    scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
    scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
    print(datasets)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        print(data)
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv[k] <- calculate_hypervolume(data, c(1,1), maximize=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
        print(paste(algorithm, dataset, experiment, "HV =", round(hv[k], 3)))
      }
      print(paste(algorithm,dataset, "Mean HV=", mean(hv)))
    }
  }
  
}


calculate_hypervolume_manual()
