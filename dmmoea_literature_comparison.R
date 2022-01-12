literature_comparison_experiments <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 6){
    print(paste0("Not enough parameters (", argnum, "/6)"))
    return(-1)
  }
  path <- args[1] # "X:\\Universidad\\dmmoea"
  params.path <- args[2] #"Tests\\a"
  algorithms.param <- args[3] # "tmix,mfuzz,dnsga2"
  ref.algorithm <- args[4] # "dnsga2"
  trials <- as.numeric(args[5]) # "1"
  evaluations <- as.numeric(args[6]) # 2000
  
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  source("moc.gapbk/R/libraries.R")
  source("moc.gapbk/R/main.R")

  tune.path <- file.path(path, params.path)
  
  #Load best params
  best_params <- read.table(file.path(tune.path, ref.algorithm, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
  
  test.path <- file.path(path, "Tests", "experiments", best_params$objectives)
  # Initialize params
  params <- init_parameters(objectives=best_params$objectives)
  params$K <- as.numeric(best_params$K)
  params$objectives <- best_params$objectives
  params$evaluations <- evaluations
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
  params$convergence_tol <- best_params$convergence_tol
  params$mutation_radius <- best_params$mutation_radius
  params$seed <- runif(1, 0, 1)*1235
  
  algorithms <- strsplit(algorithms.param, ",")[[1]]
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
      res <- diverse_memetic_nsga2(distances, params, output.exp, limits, debug=TRUE, plot=FALSE)
    }else if(algorithm == "dnsga2"){
      res <- dnsga2(distances, params, output.exp, limits, debug=TRUE, plot=FALSE)
    }else if(algorithm == "nsga2"){
      res <- nsga2(distances, params, output.exp, limits, debug=TRUE, plot=FALSE)
    }else if(algorithm == "moc.gapbk"){
      res <- run_moc_gapbk(distances, params, output.exp, limits)
    }else if(algorithm == "tmix"){
      res <- run_tmix_clust(distances, params, output.exp, limits, debug=TRUE)
    }else if(algorithm == "mfuzz"){
      res <- run_mfuzz_clust(distances, params, output.exp, limits)
    }else{
      warning("Algorithm not supported!!")
      return(NULL)
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
  params$popSize <- round(params$popSize/10)
  population <- as.data.frame(matrix(nrow=params$popSize, ncol=params$K))
  clustering.groups <- list()
  i <- 1
  while(i <= params$popSize){
    print(paste0("Solution ", i, "/", params$popSize))
    tmix.res <- TMixClust(D, nb_clusters = params$K, em_iter_max = 100, mc_em_iter_max=50)
    groups <- tmix.res$em_cluster_assignment
    medoids <- get.medoid.diss.matrix(groups, D.exp)
    if(nrow(unique(t(medoids))) == params$K){
      print("New solution:")
      print(medoids)
      population[i, ] <- medoids
      clustering.groups[[i]] <- groups
      i <- i + 1
    }
    print("Population status:")
    print(population)
  }
  row.names(population) <- 1:nrow(population)
  P.rows <- row.names(population)
  population <- remove_duplicated(population, params$K)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  P.rows <- row.names(population)
  population <- evaluate_population(population, distances, clustering.groups, params)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
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

literature_comparison_experiments()