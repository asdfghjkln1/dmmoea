

init_parameters <- function(dataset.name="arabidopsis", objectives="XieBeni", diversity_metric="jaccard"){
  
  #test.path <- "X:\\Universidad\\dmmoea\\Tests\\tuning"
  
  #Default gene expression data, folder name of distances matrices
  dataset <- dataset.name
  
  #Number of cores.
  cores <- 4
  
  #Number of objectives functions. Always 2.
  objDim <-2 
  
  #Objective functions.
  obj_names<- objectives
  maximize <- rep(FALSE, objDim)
  for(i in 1:objDim){
    if(obj_names == "Deviation" || obj_names == "Silhouette" || obj_names == "Separation"){
      maximize[i] = TRUE
    }
  }
  #maximize <- ifelse(obj_names == "Deviation" || obj_names == "Silhouette" || obj_names == "Separation", TRUE, FALSE)
  
  # Ponderation of expression matrix to biologic matrix
  alpha <- 0.5
  
  #Range of number of cluster k.
  K <- 4
  #serie_k <- 4
  
  ##Initial population params
  min_density_radius <- 0.01 # Minimum radius to calculate cluster density
  max_density_radius <- 0.2 # Maximum radius to calculate cluster density
  density_tol <- 0.01 # Percentage of dataset points as minimum acceptable density
  auto_adjust_initial_params <- TRUE # If true ignore all 3 previous parameters and auto adjust
  
  ##Memetic parameters
  agents <- 4
  phases <- 4
  sync_off <- FALSE
  sync_method <- "hc"
  diversity_level <- 4
  diversity_metric <- diversity_metric
  
  ##Evolutionary parameters
  poblacion_inicial_aleatorio<-FALSE
  mutation_radius <- 0.1
  #generation<- 20
  evaluations <- 400 # Number of maximum solutions an algorithm can evaluate until stopped
  convergencia <- 15 # Maximum number of generations without change (not used, set to 3/4 of generations)
  popSize <- 40 #total number of population
  ratCruz<-0.5
  ratMuta<-0.1
  mutType<-"all" # "mut.only", "selective", "all"
  tourSize <- 2 #ceiling(popSize*0.20)
  tolerancia <- 0.05 # Minimum convergence index to consider "no changes"
  seed <- 123
  
  return(list( dataset=dataset, cores=cores, objectives=obj_names, objDim=objDim, obj_maximize=maximize,
               alpha=alpha, K=K, min_density_radius=min_density_radius, max_density_radius=max_density_radius, density_tol=density_tol, 
               auto_adjust_initial_params=auto_adjust_initial_params, is_random_population=poblacion_inicial_aleatorio,
               popSize=popSize, mating_rate=ratCruz, mutation_rate=ratMuta, mutation_type=mutType, mutation_radius=mutation_radius,
               tourSize=tourSize, convergence_tol=tolerancia, evaluations=evaluations, agents=agents, phases=phases, 
               diversity_level=diversity_level, diversity_metric=diversity_metric, sync_off=sync_off, sync_method=sync_method, seed=seed))
              #obj1_criteria=obj1_criteria, obj2_criteria=obj2_criteria, 
}