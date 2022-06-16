test_best_configurations <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 6){
    print(paste0("Not enough parameters (", argnum, "/6)"))
    return(-1)
  }
  path <- args[1] 
  obj_fun <- args[2]
  dataset <- args[3]
  evaluations <- as.numeric(args[4])
  trial.start <- as.numeric(args[5])
  trial.stop <- as.numeric(args[6])
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  tune.path <- file.path(path, "Tests", paste0("tuning_", obj_fun))
  test.path <- file.path(path, "Tests", "runs", obj_fun)
  #print("Searching limits in:")
  #print(file.path(tune.path, "limits.csv"))
  #limits <- read.table(file.path(tune.path, "limits.csv"), sep=",", header = TRUE)
  #print("Limits are: ")
  #print(limits)
  algorithms <- list.dirs(path=tune.path, full.names = FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    #algorithm <- strsplit(algorithms[i], "_")[[1]][1]
    if(algorithm == "figures"){ next }
    best_params <- read.table(file.path(tune.path, algorithm, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
    # Initialize params
    params <- init_parameters(objectives=best_params$objectives)
    params$K <- best_params$K
    #params$objectives <- best_params$objectives
    params$evaluations <- evaluations #best_params$evaluations
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
    params$diversity_level <- best_params$diversity_level
    params$phases <- best_params$phases
    params$agents <- best_params$agents
    params$sync_off <- ifelse(is.na(as.numeric(best_params$sync_off)), 0, as.numeric(best_params$sync_off))
    params$convergence_tol <- best_params$convergence_tol
    params$mutation_radius <- best_params$mutation_radius
    params$seed <- runif(1, 0, 1)*1235
    #datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
    #for(j in 1:length(datasets)){
    #dataset <- datasets[j]
    print("Starting dataset:")
    print(dataset)
    output.folder <- file.path(test.path, algorithm, dataset)
    #limits <- read.table(file.path(tune.path, algorithm, dataset, "limits.csv"), sep=",", row.names=NULL, header=TRUE)
    execute_tests(params, path, output.folder, algorithm, dataset, trial.start=trial.start, trial.stop = trial.stop) 
    #}
  }
  
  #get_evaluation_limits(test.path)
  #evaluate_run_results(test.path, params$obj_maximize)
  
  #plot.data <- read.table(file.path(test.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  #plot.data.diversity <- read.table(file.path(test.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  #plot_algorithm_comparison(test.path, plot.data)
  #plot_algorithm_comparison_diversity(test.path, plot.data.diversity)
  #plot_algorithm_comparison_pareto(test.path)
}

test_best_configurations_paired <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 7){
    print(paste0("Not enough parameters (", argnum, "/7)"))
    return(-1)
  }
  path <- args[1] 
  trials <- args[2]
  obj_fun <- args[3]
  evaluations <- as.numeric(args[4])
  dataset <- args[5]
  pop.size <- as.numeric(args[6])
  K <- as.numeric(args[7])
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  tune.path <- file.path(path, "Tests", paste0("tuning_", obj_fun))
  test.path <- file.path(path, "Tests", "runs", "random") #obj_fun)
  algorithms <- "nsga2" #list.dirs(path=tune.path, full.names = FALSE, recursive = FALSE)
  params <- init_parameters(dataset)
  params$K <- K
  params$popSize <- pop.size
  alpha <- 0.5
  distances <- load.gene.distance(dataset, alpha)
  ## Super hard-coded... 
  best_params <- read.table(file.path(tune.path, "nsga2", "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
  params$is_random_population <- best_params$is_random_population
  params$auto_adjust_initial_params <- best_params$auto_adjust_initial_params
  if(!is.na(params$auto_adjust_initial_params)){
    params$min_density_radius <- best_params$min_density_radius
    params$max_density_radius <- best_params$max_density_radius
    params$density_tol <- best_params$density_tol
  }
  ##
  seeds <- data.frame("id"=rep(1:trials), "dataset"=rep(dataset, trials), "seed"=rep(NA, trials))
  for(i in 1:trials){
    seed <- as.numeric(Sys.time())
    seeds[i, "seed"] <- seed
    #if(dataset=="arabidopsis" || dataset=="cell_cycle"){
    #  P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE, seed=seed)
    #}else{
      P <- generate_initial_pop(pop.size, K, distances$n.genes, seed) 
    #}
    print(dataset)
    print(algorithms)
    for(j in 1:length(algorithms)){
      algorithm <- algorithms[j]
      #algorithm <- strsplit(algorithms[i], "_")[[1]][1]
      print(algorithm)
      if(algorithm == "figures"){ next }
      print(file.path(test.path, algorithm, dataset, i))
      if(dir.exists(file.path(test.path, algorithm, dataset, i))){
        next
      }
      ### Initialize parameters ###
      best_params <- read.table(file.path(tune.path, algorithm, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
      #params$objectives <- best_params$objectives
      #params$K <- best_params$K
      params$objectives <- best_params$objectives
      params$evaluations <- evaluations - nrow(P) #best_params$evaluations
      #params$popSize <- pop.size
      params$mating_rate <- best_params$mating_rate
      params$mutation_rate <- best_params$mutation_rate
      #params$alpha <- best_params$alpha
      params$diversity_metric <- best_params$diversity_metric
      params$diversity_level <- best_params$diversity_level
      params$phases <- best_params$phases
      params$agents <- best_params$agents
      params$sync_off <- ifelse(is.na(as.numeric(best_params$sync_off)), 0, as.numeric(best_params$sync_off))
      params$convergence_tol <- best_params$convergence_tol
      params$mutation_radius <- best_params$mutation_radius
      params$seed <- seed
      output.exp <- file.path(test.path, algorithm, dataset, i)
      dir.create(output.exp, showWarnings = FALSE, recursive = TRUE)
      print(paste("Algorithm", algorithm, "dataset", dataset, "run", i, "..."))
      if(algorithm == "dmnsga2"){
        res <- diverse_memetic_nsga2(distances, params, output.exp, initial_population=P, limits=NULL, debug=FALSE, plot=FALSE)
      }else if(algorithm == "dnsga2"){
        res <- dnsga2(distances, params, output.exp, initial_population=P, limits=NULL, debug=FALSE, plot=FALSE)
      }else if(algorithm == "nsga2"){
        res <- nsga2(distances, params, output.exp, initial_population=P, limits=NULL, debug=FALSE, plot=FALSE)
      }else{
        print("Algorithm not supported!!")
      }
      print("Finished. Evaluating solutions...")
      print(dirname(output.exp))
      print(i)
      evaluate_solutions(res$population, res$clustering, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp), i, algorithm, dataset, plot=FALSE)
    }
  }
}

execute_tests <- function(params, path, output.folder, algorithm, dataset, trial.start=1, trial.stop=31){
  #setwd(path)
  #source("dmmoea_functions.R")
  #source("dmmoea_parameters.R")
  #source("dmmoea_libraries.R")
  #source("dmmoea_distances.R")
  #source("dmmoea_irace_conf.R")
  algorithm.name <- strsplit(algorithm, "_")[[1]][1]

  distances <- load.gene.distance(dataset, params$alpha)
  for(i in trial.start:trial.stop){
    output.exp <- file.path(output.folder, i)#file.path(basename(params$test.path), "Debug", "test")
    if(dir.exists(output.exp)){
      next
    }
    #output.exp <- get_new_dirname(output.exp)
    print(paste0("Starting ", algorithm, " in ", dataset, " run: ", i))
    dir.create(output.folder, showWarnings=FALSE, recursive=TRUE)
    exp.id <- basename(output.exp)
    dir.create(file.path(output.exp), recursive = TRUE, showWarnings = FALSE) 
    
    save_timestamps(status=0,output.path = output.exp)
    if(algorithm.name == "dmnsga2"){
      #print("DMNSGA2")
      res <- diverse_memetic_nsga2(distances, params, output.exp, debug=TRUE, plot=FALSE)
    }else if(algorithm.name == "dnsga2"){
      #print("DNSGA2")
      res <- dnsga2(distances, params, output.exp, debug=TRUE, plot=FALSE)
    }else if(algorithm.name == "nsga2"){
      #print("NSGA2")
      res <- nsga2(distances, params, output.exp, debug=TRUE, plot=FALSE)
    }else{
      print("Algorithm not supported!!")
    }
    save_timestamps(status=1,output.path = output.exp)
    
    evaluate_solutions(res$population, res$clustering, distances, params$K, 
                       params$objDim, params$obj_maximize, dirname(output.exp), exp.id, algorithm, dataset, plot=FALSE)
  }
}

test_best_configurations()
#test_best_configurations_paired()
