## Parameters based in Jorge Párraga's work

parameters=function(){
  
  #Path: Default gene expression data, folder name of distances matrices
  path="arabidopsis"
  
  #Number of cores.
  nucleos<- 4
  #Objective functions.
  nombre_objetivos<-c("XieBeni","XieBeni") 
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
  
  ##Evolutionary parameters
  poblacion_inicial_aleatorio<-TRUE
  generation<- 100
  convergencia <- 10 #Integer. Commonly 10.
  popSize <- 40 #total number of population
  ratCruz<-0.80
  ratMuta<-0.01
  tourSize <- 2 #ceiling(popSize*0.20)
  tolerancia <- 0.05 # ADDED
  evaluations <- 5e6 # Number of maximum solutions an algorithm can search until stopped
  seed <- 123
  
  return(list( path=path, cores=nucleos, objectives=nombre_objetivos, objDim=objDim, neighbors=vecinos,
               alpha_series=serie_alfa, k_series=serie_k, is_random_population=poblacion_inicial_aleatorio,
               generations=generation, popSize=popSize, ratCruz=ratCruz, ratMuta=ratMuta, tourSize=tourSize,
               tol=tolerancia, evaluations=evaluations, convergency.limit=convergencia, seed=seed))
              #obj1_criteria=obj1_criteria, obj2_criteria=obj2_criteria, 
}
