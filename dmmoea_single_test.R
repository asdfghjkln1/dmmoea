single_test <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 4){
    print(paste0("Not enough parameters (", argnum, "/4)"))
    return(-1)
  }
  path <- args[1] #"X:\\Universidad\\dmmoea" #args[1] # "X:\\Universidad\\dmmoea" #args[1] #
  results.path <- args[2] #"Tests\\operators" # args[2] #"Tests\\runs\\XieBeni" # args[2]
  param.path <- args[3]
  runs <- args[4]
  
  setwd(path)
  library(ggpubr)
  library(rstatix)
  library(ggplot2)
  library(gridExtra)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  source("moc.gapbk/R/libraries.R")
  source("moc.gapbk/R/main.R")
  
  results.path <- file.path(path, results.path)
  dir.create(results.path, recursive = TRUE, showWarnings = FALSE)
  print("Path in:")
  print(results.path)
  #params <- init_parameters()
  
  best_params <- read.table(file.path(path, param.path, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
  params <- init_parameters(objectives=best_params$objectives)
  params$K <- 6 #as.numeric(best_params$K)
  params$objectives <- best_params$objectives
  params$popSize <- as.numeric(best_params$popSize)
  params$mating_rate <- best_params$mating_rate
  params$mutation_rate <- best_params$mutation_rate
  params$alpha <- best_params$alpha
  params$is_random_population <- as.numeric(best_params$is_random_population)
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
  params$convergence_tol <- -1 #best_params$convergence_tol
  params$mutation_radius <- best_params$mutation_radius
  params$evaluations <- 4000 #params$popSize*(generations+1)
  
  time <- as.data.frame(matrix(ncol=4, nrow=0))
  colnames(time) <- c("id", "Algorithm", "Dataset", "elapsed")
  if(!file.exists(file.path(results.path, "time.csv"))){
    write.table(time, file=file.path(results.path, "time.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
  algorithms <- c("mfuzz", "moc.gapbk", "tmix", "dmnsga2", "dnsga2")
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    params$dataset <- dataset
    distances <- load.gene.distance(dataset, params$alpha)
    test_algorithm(algorithms, params, distances, dataset, results.path, runs)
  }
  #}
  print("Tests finished!...")
}

test_algorithm <- function(algorithms, params, distances, dataset, results.path, runs){
  
  time <- read.table(file.path(results.path, "time.csv") , sep=",", header = TRUE, row.names=NULL)
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    #P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed)
    for(j in 1:length(algorithms)){
      algorithm <- algorithms[j]
      output.path <- file.path(results.path, algorithm, dataset)
      print(paste0("algorithm ", algorithm, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      start <- Sys.time()
      if(algorithm == "dmnsga2"){
        res <- diverse_memetic_nsga2(distances, params, output.exp, limits=limits, debug=TRUE, plot=FALSE)
      }else if(algorithm == "dnsga2"){
        res <- dnsga2(distances, params, output.exp,  limits=limits, debug=TRUE, plot=FALSE)
      }else if(algorithm == "nsga2"){
        res <- nsga2(distances, params, output.exp,  limits=limits, debug=TRUE, plot=FALSE)
      }else if(algorithm == "moc.gapbk"){
        res <- run_moc_gapbk(distances, params, output.exp, limits, debug=TRUE)
      }else if(algorithm == "tmix"){
        res <- run_tmix_clust(distances, params, output.exp, limits, debug=TRUE)
      }else if(algorithm == "mfuzz"){
        res <- run_mfuzz_clust(distances, params, output.exp, limits)
      }else if(algorithm == "random"){
        res <- run_random(distances, params, output.exp, limits)
      }else {
        warning("Algorithm not supported!!")
        return(NULL)
      }
      end <- Sys.time()
      elapsed <- end - start
      time.elapsed <- data.frame(i, algorithm, dataset, as.numeric(elapsed))
      write.table(time.elapsed, file=file.path(results.path, "time.csv"), sep=",", append = TRUE,  quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      evaluate_solutions(res$population, res$clustering, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp), as.character(i), algorithm, dataset, plot=TRUE)
    }
  } 
}

select_pareto_solutions <- function(exp.path, limit.run=Inf){
  folder.path <- file.path(exp.path)
  algorithms <- list.dirs(path=folder.path, full.names=FALSE, recursive = FALSE)
  algorithms <- algorithms[!(algorithms %in% c("figures", "nsga2", "random"))]
  plot.data <- as.data.frame(matrix(nrow=0, ncol=5))
  colnames(plot.data) <- c("f1", "f2", "rnkIndex", "Dataset", "Algorithm")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    datasets <- list.dirs(path=file.path(folder.path, algorithm), recursive = FALSE, full.names=FALSE)
    for(j in 1:length(datasets)){
      pareto.dataset <- as.data.frame(matrix(nrow=0, ncol=2)) #** N of objectives hardcoded! **
      dataset <- datasets[j]
      print(paste0("Starting dataset ", dataset, " in ", algorithm))
      dataset.path <- file.path(folder.path, algorithm, dataset)
      experiments <- list.dirs(path=dataset.path, recursive = FALSE, full.names=FALSE)
      experiments <- experiments[as.numeric(experiments) <= limit.run]
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        if(file.exists(file.path(dataset.path, experiment, paste0(experiment, ".csv")))){
          obj.values <- read.table(file=file.path(dataset.path, experiment, paste0(experiment, ".csv")), sep=",", header = FALSE, row.names=NULL)
          pareto.dataset <- rbind(pareto.dataset, obj.values)
        }
      }
      ranking <- nsga2R::fastNonDominatedSorting(as.matrix(pareto.dataset))
      rnkIndex <- integer(nrow(pareto.dataset))
      i <- 1
      while (i <= length(ranking)) {
        rnkIndex[ranking[[i]]] <- i
        i <- i + 1
      } 
      pareto.dataset[, "rnkIndex"] <- rnkIndex
      pareto.dataset <- pareto.dataset[order(rnkIndex), ]
      pareto.dataset <- pareto.dataset[pareto.dataset$rnkIndex == 1, ]
      colnames(pareto.dataset) <- c("f1", "f2", "rnkIndex")
      pareto.dataset[, "Dataset"] <- rep(dataset, nrow(pareto.dataset))
      pareto.dataset[, "Algorithm"] <- rep(algorithm, nrow(pareto.dataset))
      plot.data <- rbind(plot.data, pareto.dataset)
    }
  }
  plot.data.norm <- as.data.frame(matrix(ncol=6, nrow=0))
  colnames(plot.data.norm) <- c("f1", "f2", "rnkIndex", "Dataset", "Algorithm")
  
  # Generate ideal pareto for each dataset
  for(j in 1:length(datasets)){
    data <- plot.data[plot.data$Dataset == datasets[j], ]
    ranking <- nsga2R::fastNonDominatedSorting(as.matrix(data[, c("f1", "f2")]))
    rnkIndex <- integer(nrow(data))
    i <- 1
    while (i <= length(ranking)) {
      rnkIndex[ranking[[i]]] <- i
      i <- i + 1
    }
    data[, "rnkIndex"] <- rnkIndex
    ideal.pareto <- data[data$rnkIndex == 1, ]
    ideal.pareto[, "Algorithm"] <- rep("Ideal pareto", nrow(ideal.pareto))
    plot.data <- rbind(ideal.pareto, plot.data)
  }
    
  
  datasets <- unique(plot.data$Dataset)
  #plot.data$Algorithm <- factor(plot.data$Algorithm, levels=c("nsga2", "dnsga2", "dmnsga2", "Ideal pareto"))
  plot.data$Algorithm <- factor(plot.data$Algorithm, levels=c("dnsga2", "dmnsga2", "moc.gapbk", "tmix", "mfuzz", "Ideal pareto"))
  
  for(i in 1:length(datasets)){
    dataset <- datasets[i]
    data <- plot.data[plot.data$Dataset == dataset, ]
    dataset.path <- file.path(folder.path, algorithms[1], dataset) # A little hardcoded, but it should work
    limits <- read.table(file.path(dataset.path, "limits.csv"), sep=",", header = TRUE, row.names=NULL)
    data.norm <- data
    data.norm[, c("f1", "f2")] <- normalise_pareto(as.matrix(data.norm[, c("f1", "f2")]), limits=limits)
    plot.data.norm <- rbind(plot.data.norm, data.norm)
  }
  write.table(plot.data, file=file.path(folder.path, "data_pareto.csv"), sep=",", row.names = FALSE, col.names = TRUE)
  write.table(plot.data.norm, file=file.path(folder.path, "data_pareto_norm.csv"), sep=",", row.names = FALSE, col.names = TRUE)
}



run_moc_gapbk <- function(distances, params, output.exp, limits, debug=TRUE){
  
  dmatrix1 <- distances$exp.dist
  dmatrix2 <- distances$bio.dist
  num_k <- params$K
  
  # Call MOC_GaPBK with default values.
  # Set local_search=TRUE because its reccomended
  # Set generation as a very high number and set stop criteria by evaluations used.
  res <- moc.gabk(dmatrix1, dmatrix2, num_k, generation=999999, pop_size=10, rat_cross=0.80, 
                  rat_muta=0.01, tour_size=2, neighborhood=0.10, local_search=TRUE, 
                  cores=params$cores, evaluations=params$evaluations, output.path=output.exp, debug=debug)
  colnames(res$population) <- c(paste0("V", 1:params$K), "f1", "f2", "rnkIndex", "density")
  return(res)
}

run_tmix_clust <- function(distances, params, output.exp, limits, debug=TRUE){
  dir.create(file.path(output.exp), recursive = TRUE, showWarnings = FALSE) 
  if(debug){ 
    output.log.file <- file.path(output.exp, "log.txt")
    sink(output.log.file, append=FALSE)
    writeLines(c(""), output.log.file)
    print("Initiating TMixClust...")
  }
  D <- as.data.frame(distances$data.matrix)
  D.exp <- distances$exp.dist
  #params$popSize <- round(params$popSize/5)
  population <- as.data.frame(matrix(nrow=params$popSize, ncol=params$K))
  clustering.groups <- list()
  i <- 1
  while(i <= params$popSize){
    print(paste0("Solution ", i, "/", params$popSize))
    tmix.res <- TMixClust(D, nb_clusters = params$K, em_iter_max = 1000/params$popSize, mc_em_iter_max=2)
    groups <- tmix.res$em_cluster_assignment
    medoids <- get.medoid.diss.matrix(groups, D.exp)
    if(nrow(unique(t(medoids))) == params$K){
      #print("New solution:")
      print(medoids)
      population[i, ] <- medoids
      clustering.groups[[i]] <- groups
      i <- i + 1
    }
    #print("Population status:")
    #print(population)
  }
  row.names(population) <- 1:nrow(population)
  P.rows <- row.names(population)
  population <- remove_duplicated(population, params$K)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  P.rows <- row.names(population)
  population <- evaluate_population(population, distances, clustering.groups, params)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  print("Popualtion Results:")
  print(population)
  sink(type="output")
  return(list("population"=population, "clustering"=clustering.groups))
}

run_mfuzz_clust <- function(distances, params, output.exp, limits){
  D <- as.matrix(distances$data.matrix)
  D.exp <- distances$exp.dist
  eset <- new("ExpressionSet",exprs=D)
  eset <- standardise(eset)
  mest <- mestimate(eset)
  
  population <- as.data.frame(matrix(nrow=params$popSize, ncol=params$K))
  clustering.groups <- list()
  i <- 1
  set.seed(Sys.time())
  while(i <= params$popSize){
    mfuzz.res <-mfuzz(eset, c=params$K, m=mest)
    groups <- mfuzz.res$cluster
    medoids <- get.medoid.diss.matrix(groups, D.exp)
    if(nrow(unique(t(medoids))) == params$K){
      population[i, ] <- medoids
      clustering.groups[[i]] <- groups
      i <- i + 1
    }
  }
  row.names(population) <- 1:nrow(population)
  P.rows <- row.names(population)
  population <- remove_duplicated(population, params$K)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  P.rows <- row.names(population)
  population <- evaluate_population(population, distances, clustering.groups, params)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  return(list("population"=population, "clustering"=clustering.groups))
}

run_random <- function(distances, params, output.exp, limits){
  
  P <- params$Pop.size * 10
  population <- generate_initial_pop(P, params$K, distances$n.genes) 
  
  row.names(population) <- 1:nrow(population)
  P.rows <- row.names(population)
  population <- remove_duplicated(population, params$K)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  P.rows <- row.names(population)
  population <- evaluate_population(population, distances, clustering.groups, params)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  
  return(list("population"=population, "clustering"=clustering.groups))
}


get.medoid.diss.matrix <- function(groups, dist){
  k <- length(unique(groups))
  medoids <- as.data.frame(matrix(ncol=k, nrow=1))
  for(i in 1:k){
    elem.group <- dist[groups == i, groups == i]
    dist.sum <- apply(elem.group,2, sum)
    medoids[1, i] <- which.min(dist.sum)
  }
  return(medoids)
}

single_test()