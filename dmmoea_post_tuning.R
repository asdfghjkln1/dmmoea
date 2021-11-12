execute_tests <- function(params, path, output.folder, algorithm, dataset, n.times=1){
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  distances <- load.gene.distance(dataset, params$alpha)
  for(i in 1:n.times){
    output.exp <- file.path(output.folder, i)#file.path(basename(params$test.path), "Debug", "test")
    if(dir.exists(output.exp)){
      next
    }
    #output.exp <- get_new_dirname(output.exp)
    print(paste0("Starting ", algorithm, " in ", dataset, " run: ", i))
    dir.create(output.folder, showWarnings=FALSE, recursive=TRUE)
    exp.id <- basename(output.exp)
    if(algorithm == "dmnsga2"){
      #print("DMNSGA2")
      res <- diverse_memetic_nsga2(distances, params, output.exp, debug=TRUE, plot=TRUE)
    }else if(algorithm == "dnsga2"){
      #print("DNSGA2")
      res <- dnsga2(distances, params, output.exp, debug=TRUE, plot=TRUE)
    }else if(algorithm == "nsga2"){
      #print("NSGA2")
      res <- nsga2(distances, params, output.exp, debug=TRUE, plot=TRUE)
    }else{
      print("Algorithm not supported!!")
    }
    
    evaluate_solutions(res$population, res$clustering, distances, params$K, 
                       params$objDim, params$obj_maximize, dirname(output.exp), exp.id, algorithm, dataset, plot=TRUE)
  }
}

get_evaluation_limits <- function(path){
  max.f1 <- 0
  min.f1 <- Inf
  max.f2 <- 0
  min.f2 <- Inf
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){
      next
    }
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      exp.path <- file.path(path, algorithm, dataset)
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE)
        max.values <- apply(data, 2, max)
        min.values <- apply(data, 2, min)
        if(max.values[1] > max.f1){
          max.f1 <- unname(max.values[1])
        }
        if(max.values[2] > max.f2){
          max.f2 <- unname(max.values[2])
        }
        if(min.values[1] < min.f1){
          min.f1 <- unname(min.values[1])
        }
        if(min.values[2] < min.f2){
          min.f2 <- unname(min.values[2])
        }
      }
    }
  }
  limits <- data.frame("min.f1"=min.f1, "max.f1"=max.f1, "min.f2"=min.f2, "max.f2"=max.f2)
  write.table(limits, file=file.path(path, "limits.csv"), sep=",", append=FALSE, row.names = FALSE, quote = FALSE)
}

evaluate_run_results <- function(path){
  limits <- read.csv(file.path(path, "limits.csv"), header = TRUE)
  scaler.f1 <- function(x){ (x-limits$min.f1)/(limits$max.f1-limits$min.f1) }
  scaler.f2 <- function(x){ (x-limits$min.f2)/(limits$max.f2-limits$min.f2) }
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  plot.data <- as.data.frame(matrix(nrow=0, ncol=4))
  colnames(plot.data) <- c("id", "Algorithm", "Dataset", "Hypervolume")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      exp.path <- file.path(path, algorithm, dataset)
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE)
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv <- eaf::hypervolume(data=data, reference=c(2,2), maximise = TRUE)
        values <- data.frame("id"=experiment, "Algorithm"=algorithm, "Dataset"=dataset,"Hypervolume"=hv)
        plot.data <- rbind(plot.data, values)
      }
    }
  }
  write.table(plot.data, file=file.path(path, "plot_data.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
}



test_best_configurations <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 2){
    print(paste0("Not enough parameters (", argnum, "/2)"))
    return(-1)
  }
  path <- args[1] 
  trials <- args[2]
  
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  tune.path <- file.path(path, "Tests", "tuning")
  test.path <- file.path(path, "Tests", "runs")
  
  algorithms <- list.dirs(path=tune.path, full.names = FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    best_params <- read.table(file.path(tune.path, algorithm, "best_configurations.csv"), sep=",", header=TRUE)
    # Initialize params
    params <- init_parameters()
    params$K <- best_params$K
    params$objectives <- best_params$objectives
    params$evaluations <- best_params$evaluations
    params$popSize <- best_params$popSize
    params$mating_rate <- best_params$mating_rate
    params$mutation_rate <- best_params$mutation_rate
    params$tourSize <- best_params$tourSize
    params$alpha <- best_params$alpha
    params$is_random_population <- best_params$is_random_population
    params$auto_adjust_initial_params <- best_params$auto_adjust_initial_params
    if(!is.na(params$auto_adjust_initial_params)){
      params$min_density_radius <- best_params$min_density_radius
      params$max_density_radius <- best_params$max_density_radius
      params$density_tol <- best_params$density_tol
    }
    params$phases <- best_params$phases
    params$agents <- best_params$agents
    params$sync_off <- best_params$sync_off
    params$convergence_tol <- best_params$convergence_tol
    params$mutation_radius <- best_params$mutation_radius
    params$seed <- runif(1, 0, 1)*1235
    datasets <- list.dirs(path=file.path(tune.path, algorithm),full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      output.folder <- file.path(test.path, algorithm, dataset)
      execute_tests(params, path, output.folder, algorithm, dataset, n.times=trials) 
      
    }
  }
  
  get_evaluation_limits(test.path)
  evaluate_run_results(test.path)
  
  plot.data <- read.table(file.path(test.path, "plot_data.csv"), sep=",", header=TRUE)
  plot_algorithm_comparison(test.path, plot.data)
  plot_algorithm_comparison_pareto(test.path)
}


test_best_configurations()
