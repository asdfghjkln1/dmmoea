literature_comparison_experiments <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 6){
    print(paste0("Not enough parameters (", argnum, "/6)"))
    return(-1)
  }
  path <- args[1]
  params.path <- args[2]
  algorithms.param <- args[3]
  ref.algorithm <- args[4]
  trials <- args[5]
  evaluations <- args[6]
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  source("moc.gapbk/R/libraries.R")
  source("moc.gapbk/R/main.R")
  
  tune.path <- file.path(path, params.path)
  test.path <- file.path(path, "Tests", "experiments")
  
  #Load best params
  best_params <- read.table(file.path(tune.path, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
  # Initialize params
  params <- init_parameters(objectives=best_params$objectives)
  params$K <- best_params$K
  #params$objectives <- best_params$objectives
  params$evaluations <- evaluations
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
  
  algorithms <- strsplit(algorithms.param, ",")[[1]]
  #algorithms <- c(literature.algorithm) # list.dirs(path=tune.path, full.names = FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    print("Starting algorithm:")
    print(algorithm)
    if(algorithm == "figures"){ next }
    datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      print("Starting dataset:")
      print(dataset)
      output.folder <- file.path(test.path, algorithm, dataset)
      limits <- read.table(file.path(tune.path, ref.algorithm, dataset, "limits.csv"), sep=",", row.names=NULL, header=TRUE)
      execute_tests(params, path, output.folder, algorithm, dataset, limits, n.times=trials) 
    }
  }
  
  #get_evaluation_limits(test.path)
  #evaluate_run_results(test.path, params$obj_maximize)
  
  #plot.data <- read.table(file.path(test.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  #plot.data.diversity <- read.table(file.path(test.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  #plot_algorithm_comparison(test.path, plot.data)
  #plot_algorithm_comparison_diversity(test.path, plot.data.diversity)
  #plot_algorithm_comparison_pareto(test.path)
}

execute_tests <- function(params, path, output.folder, algorithm, dataset, limits, n.times=1){
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
      res <- diverse_memetic_nsga2(distances, params, output.exp, limits, debug=TRUE, plot=TRUE)
    }else if(algorithm == "dnsga2"){
      #print("DNSGA2")
      res <- dnsga2(distances, params, output.exp, limits, debug=TRUE, plot=TRUE)
    }else if(algorithm == "nsga2"){
      #print("NSGA2")
      res <- nsga2(distances, params, output.exp, limits, debug=TRUE, plot=TRUE)
    }else if(algorithm == "moc.gapbk"){
      res <- run_moc_gapbk(distances, params, output.exp, limits)
    }else{
      print("Algorithm not supported!!")
    }
    
    evaluate_solutions(res$population, res$clustering, distances, params$K, 
                       params$objDim, params$obj_maximize, dirname(output.exp), exp.id, algorithm, dataset, plot=TRUE)
  }
}

run_moc_gapbk <- function(distances, params, output.exp, limits){
  
  dmatrix1 <- distances$exp.dist
  dmatrix2 <- distances$bio.dist
  num_k <- params$K
  
  # Call MOC_GaPBK with default values.
  # Set local_search=TRUE because its reccomended
  # Set generation as a very high number and set stop criteria by evaluations used.
  res <- moc.gabk(dmatrix1, dmatrix2, num_k, generation=999999, pop_size=10, rat_cross=0.80, 
           rat_muta=0.01, tour_size=2, neighborhood=0.10, local_search=TRUE, 
           cores=params$cores, evaluations=params$evaluations, output.path=output.exp, debug=TRUE)
  colnames(res$population) <- c(paste0("V", 1:params$K), "f1", "f2", "rnkIndex", "density")
  return(res)
}

literature_comparison_experiments()
