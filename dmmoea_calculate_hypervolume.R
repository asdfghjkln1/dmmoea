calculate_hypervolume_manual <- function(){
  #library(eaf)
  source("dmmoea_functions.R")
  source("dmmoea_libraries.R")
  args <- commandArgs(trailingOnly = TRUE)
  
  path <- args[1]
  results.path <- args[2]
  dataset <- args[3]

  setwd(path)
  source("dmmoea_functions.R")
  #source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  #source("dmmoea_irace_conf.R")
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  #path <- "~/dmmoea/Tests/runs/BallHall"
  
  #algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  #print(algorithms)
  algorithms <- c("dmnsga2", "dnsga2", "nsga2")

  #limits <- read.table(file.path(path, "dnsga2", dataset, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
  #print(limits)
  #data <- read.table(file.path(path, "data_pareto.csv"), sep=",", header=TRUE, row.names=NULL)
  #scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
  #scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
  
  ideal <- data[data$Algorithm == "Ideal Pareto"]
  ideal <- ideal[!duplicated(ideal[, 1:2]) , ]
  print(paste0("Pareto front size: ", nrow(ideal)))
  
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    exp.path <- file.path(path, algorithm, dataset)
    limits <- read.table(file.path(exp.path, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
    scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
    scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
    contribution <- as.data.frame(matrix(nrow = 0, ncol = 7))
    colnames(contribution) <- c("id", "Algorithm", "Dataset", "hypervolume", "contributed", "total", "ratio")
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      ideal.dataset <- ideal[ideal$Dataset == dataset, ]
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        print(data)
        
        nrow.pareto <- nrow(data)
        contrib <- do.call(paste0, data) %in% do.call(paste0, ideal)
        contrib <- sum(contrib)
        print(paste0("Exp ", experiment, ": ", contrib, " / ", nrow.pareto, " -> ", round(contrib/nrow.pareto, 2)))
        #print(paste0("Experiment ", experiment, "has contributed ", contrib, " solutions out of ", nrow.pareto))
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv <- calculate_hypervolume(data, c(1,1), maximize=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
        contribution <- rbind(contribution, c(experiment, algorithm, dataset, hv, contrib, nrow.pareto, round(contrib/nrow.pareto, 2)))
        print(paste(algorithm, dataset, experiment, "HV =", round(hv[k], 3)))
      }
      print(paste(algorithm,dataset, "Mean HV=", mean(hv)))
    }
  }
  
  write.table(contribution, file=file.path(path, "contribution.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
}

calculate_hypervolume <- function(pareto, point, maximise=FALSE){
  pareto <- pareto[order(pareto[, 1], decreasing=TRUE), ]
  if(maximise){
    hv <- (pareto[1, 1] - point[1])*(pareto[1, 2] - point[2])
    if(nrow(pareto) == 1){
      return(hv)
    }
    for(i in 1:(nrow(pareto)-1)){
      h <- (pareto[i+1, 1] - pareto[i, 1])
      w <- (point[2] - pareto[i+1, 2])
      hv <- hv + w*h
    }
  }else{
    hv <- (point[1] - pareto[1, 1])*(point[2] - pareto[1, 2])
    if(nrow(pareto) == 1){
      return(hv)
    }
    for(i in 1:(nrow(pareto)-1)){
      w <- (pareto[i, 1] - pareto[i+1, 1])
      h <- (point[2] - pareto[i+1, 2])
      hv <- hv + w*h
    }
  }
  return(hv)
}


calculate_hypervolume_manual_2 <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  path <- args[1]
  results.path <- args[2]
  dataset <- args[3]
  
  setwd(path)
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  
  #library(eaf)
  #source("dmmoea_functions.R")
  #source("dmmoea_libraries.R")
  
  #path <- "~/dmmoea/Tests/runs/BallHall"
  
  algorithms <- c("dmnsga2", "dnsga2", "nsga2")
  #dataset <- "arabidopsis"
  
  limits <- read.table(file.path(path, "limits", dataset, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
  print(limits)
  print(file.path(path, "data_pareto.csv"))
  data <- read.table(file.path(path, "data_pareto.csv"), sep=",", header=TRUE, row.names=NULL)
  scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
  scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
  
  data <- data[data$Dataset == dataset, ]
  print("Length: ")
  print(nrow(data))
  hv <- c()
  #ideal <- data[data$Algorithm == "Ideal Pareto"]
  #ideal <- ideal[!duplicated(ideal[, 1:2]) , ]
  #print(paste0("Pareto front size: ", nrow(ideal)))
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    data.alg <- data[data$Algorithm == algorithm, ]
    print(paste0("Length algorithm ", algorithm, ": ", nrow(data.alg)))
    data.alg <- data.alg[!duplicated(data.alg[, 1:2]), ]
    print(paste0("Length algorithm ", algorithm, " (after): ", nrow(data.alg)))
    #print(data.alg)
    
    data.norm <- data.frame("f1"=scaler.f1(data.alg[, 1]), "f2"=scaler.f2(data.alg[, 2]))
    hv[i] <- eaf::hypervolume(data.norm, c(1,1), maximise=FALSE)  #calculate_hypervolume(data.norm, c(1,1), maximise=FALSE) 
    print(paste(algorithm, dataset, "HV =", hv[i]))
  }
  print(algorithms)
  print(hv)
}


calculate_hypervolume_manual()
