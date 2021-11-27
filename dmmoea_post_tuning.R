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
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
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

evaluate_run_results <- function(path, maximize=FALSE, alpha=0.5){
  limits <- read.csv(file.path(path, "limits.csv"), header = TRUE)
  scaler.f1 <- function(x){ (x-limits$min.f1)/(limits$max.f1-limits$min.f1) }
  scaler.f2 <- function(x){ (x-limits$min.f2)/(limits$max.f2-limits$min.f2) }
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  plot.data <- as.data.frame(matrix(nrow=0, ncol=6))
  plot.data.diveristy <- as.data.frame(matrix(nrow=0, ncol=6))
  colnames(plot.data) <- c("id", "Algorithm", "Dataset", "Hypervolume", "Silhouette", "Delta")
  colnames(plot.data.diversity) <- c("id", "Algorithm", "Dataset", "Metric", "Diversity", "Cluster_Ratio")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      distances <- load.gene.distance(dataset, alpha=alpha)
      exp.path <- file.path(path, algorithm, dataset)
      evaluation.file <- read.table(file.path(exp.path, "evaluations.csv"), header=TRUE, sep=",", row.names=NULL)
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        data.pareto <- read.table(file.path(exp.path, experiment, "population.csv"), sep=",", header=FALSE, row.names=NULL)
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv <- calculate_hypervolume(data, c(1,1), maximize)
        sil <- evaluation.file[k, "avg_sil"]
        delta <- evaluation.file[k, "delta"]
        diversity.jaccard <- diversity_analysis(data.pareto, distances, metric="jaccard", 
                                        path=file.path(exp.path, experiment), alpha=alpha, plot=TRUE)
        diversity.NMI <- diversity_analysis(data.pareto, distances, metric="NMI", 
                                        path=file.path(exp.path, experiment), alpha=alpha, plot=TRUE)
        values.diversity <- data.frame("id"=rep(experiment, 2), "Algorithm"=rep(algorithm, 2), 
                                       "Dataset"=rep(dataset,2), "Metric"=c("jaccard", "NMI"), 
                                       "Diversity"=c(diversity.jaccard$avg.dist, diversity.NMI$avg.dist),
                                       "Cluster_Ratio"=c(diversity.jaccard$k.best, diversity.NMI$k.best))
        values <- data.frame("id"=experiment, "Algorithm"=algorithm, "Dataset"=dataset,
                             "Hypervolume"=hv, "Silhouette"=sil, "Delta"=delta)
        plot.data <- rbind(plot.data, values)
        plot.data.diversity <- rbind(plot.data.diversity, values.diversity)
      }
    }
  }
  write.table(plot.data, file=file.path(path, "plot_data.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(plot.data.diveristy, file=file.path(path, "plot_data_diversity.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
}

test_best_configurations <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 3){
    print(paste0("Not enough parameters (", argnum, "/3)"))
    return(-1)
  }
  path <- args[1] 
  trials <- args[2]
  obj_fun <- args[3]
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  tune.path <- file.path(path, "Tests", paste0("tuning_", obj_fun))
  test.path <- file.path(path, "Tests", "runs", obj_fun)
  
  algorithms <- list.dirs(path=tune.path, full.names = FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    best_params <- read.table(file.path(tune.path, algorithm, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
    # Initialize params
    params <- init_parameters(objectives=best_params$objectives)
    params$K <- best_params$K
    #params$objectives <- best_params$objectives
    params$evaluations <- best_params$evaluations
    params$popSize <- best_params$popSize
    params$mating_rate <- best_params$mating_rate
    params$mutation_rate <- best_params$mutation_rate
    params$alpha <- best_params$alpha
    params$is_random_population <- best_params$is_random_population
    params$auto_adjust_initial_params <- best_params$auto_adjust_initial_params
    if(!is.na(params$auto_adjust_initial_params)){
      params$min_density_radius <- best_params$min_density_radius
      params$max_density_radius <- best_params$max_density_radius
      params$density_tol <- best_params$density_tol
    }
    params$diversity_metric <- best_params$diversity_metric
    params$phases <- best_params$phases
    params$agents <- best_params$agents
    params$sync_off <- best_params$sync_off
    params$convergence_tol <- best_params$convergence_tol
    params$mutation_radius <- best_params$mutation_radius
    params$seed <- runif(1, 0, 1)*1235
    datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      print("Starting dataset:")
      print(dataset)
      output.folder <- file.path(test.path, algorithm, dataset)
      execute_tests(params, path, output.folder, algorithm, dataset, n.times=trials) 
      
    }
  }
  
  get_evaluation_limits(test.path)
  evaluate_run_results(test.path, params$obj_maximize)
  
  plot.data <- read.table(file.path(test.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data.diversity <- read.table(file.path(test.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  plot_algorithm_comparison(test.path, plot.data)
  plot_algorithm_comparison_diversity(test.path, plot.data.diversity)
  plot_algorithm_comparison_pareto(test.path)
}

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

test_best_configurations()
