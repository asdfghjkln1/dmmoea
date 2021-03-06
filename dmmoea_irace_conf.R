target.runner <- function(experiment, scenario){
  instance <- experiment$instance
  if(!is.na(scenario$parallel)){
    source("dmmoea_libraries.R")
  }
  configuration <- experiment$configuration
  actions <- strsplit(instance, " ")[[1]]
  dataset.name <- strsplit(instance, " ")[[1]][1]
  debug <- FALSE
  plot <- FALSE
  if(length(actions) > 1){
    for(i in 2:length(actions)){
      tag <- strsplit(actions[i], "=")[[1]][1]
      if(tag == "--plot"){
        plot <- as.logical(strsplit(actions[i], "=")[[1]][2])
      }else if(tag == "--print"){
        debug <- as.logical(strsplit(actions[i], "=")[[1]][2])
      }
    } 
  }
  
  #dataset <- configuration[['dataset']]
  experiment$configuration[['dataset']] <- dataset.name
  
  params <- init_parameters(dataset.name = dataset.name, objectives=configuration[['objectives']])

  algorithm <- configuration[['algorithm']]
  params$K <- as.numeric(configuration[['K']])
  #params$objectives <- configuration[['objectives']]
  params$evaluations <- as.numeric(configuration[['evaluations']])
  params$popSize <- as.numeric(configuration[['popSize']])
  params$mating_rate <- as.numeric(configuration[['mating_rate']])
  params$mutation_rate <- as.numeric(configuration[['mutation_rate']])
  params$alpha <- as.numeric(configuration[['alpha']])
  params$is_random_population <- as.numeric(configuration[['is_random_population']])
  params$auto_adjust_initial_params <- as.numeric(configuration[['auto_adjust_initial_params']])
  if(!is.na(params$auto_adjust_initial_params)){
    params$min_density_radius <- as.numeric(configuration[['min_density_radius']])
    params$max_density_radius <- as.numeric(configuration[['max_density_radius']])
    params$density_tol <- as.numeric(configuration[['density_tol']]) 
  }
  params$diversity_metric <- configuration[['diversity_metric']]
  params$diversity_level <- configuration[['diversity_level']]
  params$phases <- as.numeric(configuration[['phases']])
  params$agents <- as.numeric(configuration[['agents']])
  params$sync_off <- as.numeric(configuration[['sync_off']])
  if(is.na(params$sync_off)){ params$sync_off <- 0}
  params$convergence_tol <- as.numeric(configuration[['convergence_tol']])
  params$mutation_radius <- as.numeric(configuration[['mutation_radius']])
  params$seed <- as.numeric(experiment$seed)
  #Load distances
  distances <- load.gene.distance(dataset.name, params$alpha)
  
  output <- file.path(test.path, algorithm, params$dataset)
  exp.id <- paste0(dataset.name, "_", experiment$id.configuration)
  output.exp <- file.path(output, exp.id)
  output.exp <- get_new_dirname(output.exp)
  dir.create(output.exp, recursive = TRUE, showWarnings = FALSE)
  exp.id <- basename(output.exp)
  if (plot){ 
    limits <- read.csv(test.path, "limits.csv")
  }else {
    limits <- NULL
  }
  tic(quiet = TRUE)
  if(algorithm == "dmnsga2"){
    res <- diverse_memetic_nsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else if(algorithm == "dnsga2"){
    res <- dnsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else if(algorithm == "nsga2"){
    res <- nsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else{
    print("Algorithm not supported!!")
    return(list(cost=Inf))
  }
  t <- toc(quiet=TRUE)
  t <- unname(t[["toc"]] - t[["tic"]])
  
  
  eval <- evaluate_solutions(res$population, res$clustering, distances, params$K, 
                             params$objDim, params$obj_maximize, output, exp.id, algorithm, dataset.name, t, plot=plot)
  #end_time <- Sys.time()
  #time <- round(end_time - start_time, 2)
  #print(paste("Instance",exp.id,"ended in", time, "seconds"))
  #cost <- -1*eval$results[1, "avg_sil"]
  #print(paste0("Cost is: ", cost))
  
  return(list()) #"time"=as.numeric(time)))
}

target.runner.diversity <- function(experiment, scenario){
  instance <- experiment$instance
  if(!is.na(scenario$parallel)){
    source("dmmoea_libraries.R")
  }
  configuration <- experiment$configuration
  actions <- strsplit(instance, " ")[[1]]
  dataset.name <- strsplit(instance, " ")[[1]][1]
  debug <- FALSE
  plot <- FALSE
  if(length(actions) > 1){
    for(i in 2:length(actions)){
      tag <- strsplit(actions[i], "=")[[1]][1]
      if(tag == "--plot"){
        plot <- as.logical(strsplit(actions[i], "=")[[1]][2])
      }else if(tag == "--print"){
        debug <- as.logical(strsplit(actions[i], "=")[[1]][2])
      }
    } 
  }
  
  #dataset <- configuration[['dataset']]
  experiment$configuration[['dataset']] <- dataset.name
  
  params <- init_parameters(dataset.name = dataset.name, objectives=configuration[['objectives']])
  
  algorithm <- configuration[['algorithm']]
  params$K <- as.numeric(configuration[['K']])
  #params$objectives <- configuration[['objectives']]
  params$evaluations <- as.numeric(configuration[['evaluations']])
  params$popSize <- as.numeric(configuration[['popSize']])
  params$mating_rate <- as.numeric(configuration[['mating_rate']])
  params$mutation_rate <- as.numeric(configuration[['mutation_rate']])
  params$alpha <- as.numeric(configuration[['alpha']])
  params$is_random_population <- as.numeric(configuration[['is_random_population']])
  params$auto_adjust_initial_params <- as.numeric(configuration[['auto_adjust_initial_params']])
  if(!is.na(params$auto_adjust_initial_params)){
    params$min_density_radius <- as.numeric(configuration[['min_density_radius']])
    params$max_density_radius <- as.numeric(configuration[['max_density_radius']])
    params$density_tol <- as.numeric(configuration[['density_tol']]) 
  }
  params$diversity_metric <- configuration[['diversity_metric']]
  params$diversity_level <- configuration[['diversity_level']]
  params$phases <- as.numeric(configuration[['phases']])
  params$agents <- as.numeric(configuration[['agents']])
  params$sync_off <- as.numeric(configuration[['sync_off']])
  if(is.na(params$sync_off)){ params$sync_off <- 0}
  params$convergence_tol <- as.numeric(configuration[['convergence_tol']])
  params$mutation_radius <- as.numeric(configuration[['mutation_radius']])
  params$seed <- as.numeric(experiment$seed)
  #Load distances
  distances <- load.gene.distance(dataset.name, params$alpha)
  
  output <- file.path(test.path, algorithm, params$dataset)
  exp.id <- paste0(dataset.name, "_", experiment$id.configuration)
  output.exp <- file.path(output, exp.id)
  output.exp <- get_new_dirname(output.exp)
  dir.create(output.exp, recursive = TRUE, showWarnings = FALSE)
  exp.id <- basename(output.exp)
  if (plot){ 
    limits <- read.csv(test.path, "limits.csv")
  }else {
    limits <- NULL
  }
  tic(quiet = TRUE)
  if(algorithm == "dmnsga2"){
    res <- diverse_memetic_nsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else if(algorithm == "dnsga2"){
    res <- dnsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else if(algorithm == "nsga2"){
    res <- nsga2(distances, params, output.exp, limits=limits, debug=debug, plot=plot)
  }else{
    print("Algorithm not supported!!")
    return(list(cost=Inf))
  }
  t <- toc(quiet=TRUE)
  t <- unname(t[["toc"]] - t[["tic"]])
  
  rows <- nrow(res$population)
  print(paste0("Evaluating ", rows, "solution; Clustering number", length(res$clustering)))
  pareto.clustering <- res$clustering
  #pareto.clustering <- res$clustering[res$population[, "rnkIndex"] == min(res$population[, "rnkIndex"])]
  #if(rows < 2){
  #  pareto.clustering <- res$clustering[res$population[, "rnkIndex"] <= 2]
  #}
  diversity <- calculate_diversity_matrix(pareto.clustering, "jaccard")
  avg.dist <- mean(diversity[lower.tri(diversity, diag = FALSE)])
  res <- data.frame("id"=exp.id, "Algorithm"=algorithm, "Dataset"=dataset.name, "Diversity"=avg.dist, "time"=t)
  write.table(res, file=file.path(test.path, algorithm, "diversity.csv"), append=TRUE, sep=",", row.names = FALSE, col.names = FALSE)
  return(list("cost"=avg.dist)) #"time"=as.numeric(time)))
}

target.evaluator <- function(experiment, num.configurations, all.conf.id,scenario, target.runner.call){
  
  instance <- experiment$instance
  configuration <- experiment$configuration
  dataset.name <- strsplit(instance, " ")[[1]][1]
  
  #print(paste0("Starting evaluator for ", instance.name, ":"))
  
  exp.id <- paste0(dataset.name, "_", experiment$id.configuration)
  base.path <- file.path(test.path, configuration[["algorithm"]]) #instance.name)
  instance.path <- file.path(base.path, dataset.name, exp.id, paste0(exp.id, ".csv"))
  
  ##** CAREFUL WHEN ADDING MORE OBJECTIVES !! **## 
  if(configuration[['objectives']] == "XieBeni" || configuration[['objectives']] == "Dunn" || configuration[['objectives']] == "BallHall"){
    if(!file.exists(file.path(test.path, "limits.csv"))){
      limits <- get_normalization_limits(test.path) #normalise_results(test.path)
    }else{
      limits <- read.csv(file.path(test.path, "limits.csv"), header = TRUE) 
      #read.table(file.path(test.path, "limits.csv"), sep=",", header = TRUE, row.names = NULL)
    }
    
    max.f1 <- as.numeric(limits["max.f1"])
    max.f2 <- as.numeric(limits["max.f2"])
    min.f1 <- as.numeric(limits["min.f1"])
    min.f2 <- as.numeric(limits["min.f2"])
    
    scaler.f1 <- function(x){ max(1, (x-min.f1)/(max.f1-min.f1)) }
    scaler.f2 <- function(x){ max(1, (x-min.f2)/(max.f2-min.f2)) }
    pareto <- read.csv(instance.path, header=FALSE, sep=",")
    pareto <- data.frame("f1"=scaler.f1(pareto[, 1]), "f2"=scaler.f2(pareto[, 2]))
    pareto[1] <- lapply(pareto[1], function(x) ifelse(x < max.f1, x, max.f1))
    pareto[2] <- lapply(pareto[2], function(x) ifelse(x < max.f2, x, max.f2))
    #write.table(norm.pareto, file = instance.path, append=FALSE, sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
    maximize <- FALSE
    ref.point <- c(1,1)
  }else if(configuration[['objectives']] == "Silhouette"){
    ref.point <- c(0,0)
    pareto <- read.csv(instance.path, header=FALSE, sep=",")
    maximize <- TRUE
  }
  hv <- calculate_hypervolume(pareto, ref.point, maximize)

  res <- data.frame("id"=exp.id, "Algorithm"=configuration[['algorithm']], "Dataset"=dataset.name,"Hypervolume"=hv)
  if(file.exists(file.path(base.path, "plot_data.csv"))){
    write.table(res, file=file.path(base.path, "plot_data.csv"), sep=",", append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
  }else{
    write.table(res, file=file.path(base.path, "plot_data.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
  hv <- ifelse(maximize, hv, -hv)
  return(list("cost"=hv))
}

#target.evaluator <- function(experiment, num.configurations, all.conf.id, scenario, target.runner.call){}