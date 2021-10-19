## Parameters based in Jorge Párraga's work

parameters=function(){
  
  #Path: Default gene expression data, folder name of distances matrices
  path="arabidopsis"
  
  #Number of cores.
  nucleos<- 4
  #Objective functions.
  nombre_objetivos<-c("XieBeni","XieBeni")
  maximize <- c(FALSE, FALSE)
  #"XB-exp", "XB-bio"; "Dev-exp", "Dev-bio"; "Sep-exp", "Sep-bio"
  
  #Number of objectives functions. Always 2.
  objDim <-2 
  #Real. Commonly 0.10
  vecinos<-0.05    
  #Always 1.
  serie_alfa<-seq(1,1,-0.10) # Only 1 (uses expression data only )
  
  #Range of number of cluster k.
  serie_k<-c(4:6) 
  #serie_k <- 4
  
  ##Initial population params
  min_radius <- 0.1 # Minimum radius to calculate cluster density
  max_radius <- 0.5 # Maximum radius to calculate cluster density
  density_tol <- 0.05 # Percentage of data points as minimum acceptable density
  auto_adjust_initial_params <- TRUE # If true ignore all 3 previous parameters and auto adjust
  
  ##Memetic parameters
  agents <- 4
  phases <- 3
  
  ##Evolutionary parameters
  poblacion_inicial_aleatorio<-FALSE
  generation<- 20
  convergencia <- 15 # Maximum number of generations without change (not used, set to 3/4 of generations)
  popSize <- 20 #total number of population
  ratCruz<-0.90
  ratMuta<-0.2
  tourSize <- 2 #ceiling(popSize*0.20)
  tolerancia <- 0.05 # Minimum convergence index to consider "no changes"
  evaluations <- 5e6 # Number of maximum solutions an algorithm can search until stopped
  seed <- 123
  
  return(list( path=path, cores=nucleos, objectives=nombre_objetivos, objDim=objDim, obj_maximize=maximize, neighbors=vecinos,
               alpha_series=serie_alfa, k_series=serie_k, min_radius=min_radius, max_radius=max_radius, density_tol=density_tol, 
               auto_adjust_initial_params=auto_adjust_initial_params, is_random_population=poblacion_inicial_aleatorio,
               generations=generation, popSize=popSize, mating.rate=ratCruz, mutation.rate=ratMuta, tourSize=tourSize,
               convergence.tol=tolerancia, evaluations=evaluations, agents=agents, phases=phases, seed=seed))
              #obj1_criteria=obj1_criteria, obj2_criteria=obj2_criteria, 
}
