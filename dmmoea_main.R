########################################################################################
#           Vicente Rivera Rebolledo                                                   #
#           vicente.rivera.r@usach.cl                                                  #
# DMMOEA: Diverse Multimodal Multiobjective Evolutionary Algorithm for Gene Clustering #
########################################################################################

# Set working directory
dir <- rstudioapi::getActiveDocumentContext()$path #*** Doesnt work when sourcing ***
dir <- paste0(dirname(dir), "/")
setwd( dir )

source("dmmoea_functions.R")

source("dmmoea_parameters.R")

source("dmmoea_libraries.R")

source("dmmoea_distances.R")

params <- parameters()
distances <- load.gene.distance("arabidopsis")

diversity.metric <- "jaccard"
diversity.level <- 3

output.base <- file.path(dir, "/Tests/arabidopsis")
  
path.nsga2 <- "nsga2/test"
path.dnsga2 <- "dnsga2/test"
path.dmnsga2 <- "dmnsga2/test"

output.path.nsga2 <- file.path(output.base, path.nsga2)
output.path.dnsga2 <- file.path(output.base, path.dnsga2)
output.path.dmnsga2 <- file.path(output.base, path.dmnsga2)
K <- 4

start_time <- Sys.time()
# Diverse memetic NSGA-2 tests
source("dmmoea_functions.R")
pareto.results.1 <- diverse_memetic_nsga2(distances, K, diversity.metric, diversity.level,  params, output.path.dmnsga2)

end_time <- Sys.time()

# Extracting results
results <- evaluate_solutions(pareto.results.1$population, pareto.results.1$clustering, K, dist, alpha, params, output.dnsga2)

print(paste("DMNSGA-2 ended!. Execution finished in", round(end_time - start_time, 2), "minutes"))



start_time <- Sys.time()

pareto.results.2 <- dnsga2(distances, K, diversity.metric, diversity.level, params, output.path.dnsga2)

end_time <- Sys.time()

# Extracting results
#source("dmmoea_functions.R")
results <- evaluate_solutions(pareto.results.2$population, pareto.results.1$clustering, K, distances, params, output.path.dnsga2)

print(paste("DNSGA-2 ended!. Execution finished in", round(end_time - start_time, 2), "minutes"))






#projected.pareto <- hipervolume_projection(pareto[, (4+1):(4+params$objDim)], output)

# Original NSGA-II Tests

# Seeds used for arabidopsis: 4321 for test A, 1234 for test B

pareto.results.2 <- nsga2(dist, k, params, output.nsga2)

results <- evaluate_solutions(pareto.results.2$population, pareto.results.2$clustering, k, dist, alpha, output.nsga2)
#pareto <- NSGA2(initial_pareto, 1, param$popSize, param$evaluations, "B_main_log_moc.txt", test_output, experiment)
#write.table(pareto, path=paste0(test_path, , sep = ",",col.names = TRUE, quote = FALSE)

# Test with CRAN package no funciona...
#install.packages("moc.gapbk", dependencies = TRUE)
#library("amap")
#library("moc.gapbk")













# Only this function is needed in level 1 agent's session, this is to save loading packages for some agents
library("nsga2R")

NonDomSorting <- nsga2R::fastNonDominatedSorting



#Parameters are inside param
pareto <- DeepMemeticNSGA2(param, "B_main_log_memetic.txt", test_path, seed)


pareto_nsga2 <- read.csv(file=paste0(test_path, "B_nsga2.csv"), header=TRUE, sep=",")
pareto_memetic <- read.csv(file=paste0(test_path, "B_memetic.csv"), header=TRUE, sep=",")
pareto_moc <- read.csv(file=paste0(test_path, "B_moc.csv"), header=TRUE, sep=",")

plotPareto("Arabidopsis", "B", pareto_memetic, pareto_moc, pareto_nsga2)
plot_convergence()

DeepMemeticNSGA2 <- function(param, outputfile, output_result_file, seed){
  
  setwd("<root directory here>\\NSGA2\\MOC-GaPBK")
  
  library("doParallel")
  source("MOC-GaPBK_nuevas_funciones.R")
  
  carga_distancias=dget("MOC-GaPBK_load_distancias.R")
  carga_poblacion_inicial=dget("MOC-GaPBK_genera_poblacion_inicial.R")
  carga_guarda_pareto=dget("MOC-GaPBK_obtiene_pareto_soluciones.R")
  carga_visualiza_silueta=dget("MOC-GaPBK_obtiene_siluetas.R")
  distancias=carga_distancias(param$path)
  
  pop_size <- param$popSize
  agentslv1 <- param$agentslv1
  agentslv2 <- param$agentslv2
  phaseslv1 <- param$phaseslv1
  phaseslv2 <- param$phaseslv2
  num_k <- param$serie_k
  alfa <- param$serie_alfa

  cl <- makeCluster(agentslv1)
  registerDoParallel(cl)
  
  phase = 1
  Agents <- list()
  
  exported_func <- c("NSGA2", "MemeticNSGA2", "Log", "NonDomSorting")
  clusterExport(cl, exported_func)
  
  sink(outputfile, append=TRUE)
  writeLines(c(""), outputfile)
  
  Log("Starting Deep Memetic NSGA2!\n")

  initial_pareto <- carga_poblacion_inicial(param$poblacion_inicial_aleatorio, num_k, pop_size, distancias$numgenes, seed)
  
  pop_per_agent <- floor(pop_size / agentslv1)
  
  # Assign each agent with a random set of solutions from the solution pool
  #if(nrow(poblacion_pareto) > starting_population){
  solutions <- rep(1:agentslv1, each=pop_per_agent)
  solutions <- sample(solutions)  
  solutions <- split(1:pop_size, solutions)
  
  
  for(i in 1:agentslv1){
    Agents <- c(Agents, list( initial_pareto[ unlist( solutions[i] ) , ] ) )  
  }
  
  #}#}else{
  #  solutions <- rep(1:agents, each=floor(nrow(poblacion_pareto)/agents))
  #  solutions <- sample(solutions)  
  #}
  
  ev <- param$evaluations
  sync_evaluations <- phaseslv1*pop_per_agent^2*agentslv1^2
  Log("Sync evaluation overhead is: %d", sync_evaluations)
  ev_per_agent = floor(ev/(phaseslv1*agentslv1)) - sync_evaluations
  
  Log("Number of evaluations per agent per lv1 phase: %d", ev_per_agent)

  while (phase <= phaseslv1) {

    Agents <- foreach(i=1:agentslv1, .combine = c, .export=c("pop_per_agent", "ev_per_agent", "phase"), .inorder = FALSE) %dopar% {
      
      pareto_agent <- MemeticNSGA2(Agents[[i]], i, pop_per_agent, ev_per_agent, param, phase)
      
      return( list( pareto_agent ) )
    }
    
    Log("Level 1 Agents finished!")
    
    # Synchronize agents
    for(i in 1:agentslv1){
      for(j in 1:agentslv1){
        if(i == j) {
          next
        }
        #solutions_taken <- population / param$num_agents
        #threshold <- population - solutions_taken
        #temp <- tournamentSelection(rbind(Agents[[i]], Agents[[j]]), pop_agent_i + pop_agent_j, param$tourSize)
        Agents[[i]] <- fitness(Agents[[i]], Agents[[j]], "min", "min", pop_per_agent)
      }
    }
    Log("Phase %d synchronization ended!\n", phase)
    Log("Phase %d solutions: \n", phase)
    print(Agents)
    
    #We need to discard the other information before using it for the new phase population
    if(phase != phaseslv1){
      for(i in 1:agentslv1){
        Agents[[i]] <- Agents[[i]][ , 1:num_k]
      }
    }
    
    phase <- phase + 1
    
    
  }
    
    #best_solutions_target <- Agents[[j]][1:solutions_taken, ]
    #worst_solutions_agent <- Agents[[i]][threshold: .N, ]
    #if(best_solutions_target$obj1[1,] < worst_solutions_agent$obj1[threshold, ] &&
    #   best_solutions_target$obj2[1,] < worst_solutions_agent$obj2[threshold, ]) { #POR AHORA CON MINIMIZACION
    #  Agent[[i]][threshold:.N, ] <- best_solutions_target
    #  cat(paste(Sys.time()),"Agent",i,"takes the",solutions_taken,"best solutions of agent",j)
    #}
  
  pareto <- data.frame()
  for(i in 1:agentslv1){
    pareto <- rbind(pareto, Agents[[i]])
  }
  pareto<- pareto[!duplicated(pareto),] # Delete duplicates
  
  Log("Deep Memetic NSGA-2 Algorithm finished!!")
  
  
  # Pruning non pareto-front solutions
  objetives <- as.matrix(pareto[, c("obj1", "obj2")]) 
  ranks <- NonDomSorting(objetives) 
  pareto <- pareto[ unlist(ranks[[1]]) , ] #Select the rank 1 indices
  # No need to limit number of population now
  #if(nrow(pareto) > pop_size){
  #  pareto <- pareto[1:pop_size, ]
  #}
  
  #Obtain results
  guarda_paretos=carga_guarda_pareto(alfa, num_k,  distancias$matriz_exp,  distancias$distancia_exp, distancias$distancia_bio, pop_per_agent, pareto, FALSE)
  
  #File result
  exp.name <- paste0(output_result_file, "\\B_memetic")
  
  #save results
  write.table(guarda_paretos$pareto_clustering, file = paste0(exp.name,"_clustering.csv"),sep = ",",col.names = TRUE, quote = FALSE)
  write.table(guarda_paretos$pareto_objectives, file = paste0(exp.name,".csv"),sep = ",",col.names = TRUE, row.names = TRUE,quote = FALSE)
  write.table(guarda_paretos$pareto_silhouette, file = paste0(exp.name,"_silhouettes.csv"),sep = ",",col.names = FALSE, quote = FALSE)
  
  #show results
  Log("Showing results: \n")
  print(pareto)
  #print(guarda_paretos)

  return(pareto)
}

# experiment: name of experiment. "memetic" "nsga2" "moc"
NSGA2 <- function(populationP, agente, popSize, evaluations, outfile, test_name, experiment){
  
  sink(outfile, append=TRUE)

  source("MOC-GaPBK_librerias.R")
  source("MOC-GaPBK_funciones.R")
  #Load functions
  carga_parametro=dget("MOC-GaPBK_parametros.R")
  carga_distancias=dget("MOC-GaPBK_load_distancias.R")
  carga_guarda_pareto=dget("MOC-GaPBK_obtiene_pareto_soluciones.R")
  carga_visualiza_silueta=dget("MOC-GaPBK_obtiene_siluetas.R")
  #Functions as variables
  param=carga_parametro()
  distancias=carga_distancias(param$path)
  #Constant value
  #for(num_k in param$serie_k) {
    #varNo=num_k
  #}
  varNo <- param$serie_k
  num_k <- varNo
  alfa = param$serie_alfa
  
  #convergencia <- rbind( data.frame(), rep(0, num_k))
  #colnames(convergencia) <- paste("V", 1:num_k, sep="")
  
  # If theres only one solution, next phase will fail, so at least add one random solution
  if(nrow(populationP) == 1){
    population_random <- as.data.frame(t(array(sample(1:distancias$numgenes, num_k, replace=F))))
    #Verify singletons
    reparar<-repararSingletons(as.matrix(population_random), distancias$distancia_exp, distancias$distancia_bio, distancias$matriz_exp, distancias$numgenes, alfa)
    population_random=reparar[[2]] #2 Cause, first are the clusters
    populationP <- rbind(populationP, population_random)
  }
  
  conv_count <- 0 #Only used for non-memetic algorithms
  g <- 1
  
  # Number of evaluations in function of p 
  if(experiment == "moc"){ # Memetic has always lower population, and not scales very well
    p <- abs(popSize*0.05) # Need to scale down or evaluation magnitude will be too big
  }else{
    p <- abs(popSize*0.1)
  }
  ev <- 0 # Current evaluations
  
  while (g<=param$generation && ev < evaluations) {
    
    ######################## Poblacion P ########################
    verficaSingletonsP<-repararSingletons(populationP, distancias$distancia_exp, distancias$distancia_bio, distancias$matriz_exp, distancias$numgenes, alfa)
    tablaGruposP<-verficaSingletonsP[[1]] #1 first argument returned. It is partition (groups)
    populationP<- as.matrix(verficaSingletonsP[[2]]) #2 second argument returned. It is poblacion without singletons
    populationP<-calculaJerarquiasDensidad(popSize, populationP, tablaGruposP, param$objDim, FALSE,distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, num_k,alfa,varNo, param$nombre_objetivos)
    if(g == 1){
      last_pareto <- populationP
    }
    
    Log("[Agent %d] Entering selection, crossover, mutation for generation %d\n", agente, g)
    ######################## Selection, Crossover, Mutation  ########################
    #Selection
    matingPool <- tournamentSelection(populationP, popSize, param$tourSize)
    populationQ <- t(sapply(1:popSize, function(u) array(rep(0,num_k))))
    #Crossover and Mutation
    cruza <-cruzamiento_k_puntos(popSize,num_k,matingPool, populationQ, param$ratCruz)
    populationQ<-controlaFactibilidad(cruza,populationQ,distancias$numgenes)
    muta <-mutacion_controller_random(popSize,num_k,populationQ, param$ratMuta, distancias$numgenes)
    populationQ<-controlaFactibilidad(muta,populationQ,distancias$numgenes)
    ######################## Population Q  ########################
    Log("[Agent %d] Generating Population Q for generation %d\n", agente, g)
    verficaSingletonsQ<-repararSingletons(populationQ, distancias$distancia_exp, distancias$distancia_bio, distancias$matriz_exp, distancias$numgenes, alfa)
    tablaGruposQ<-verficaSingletonsQ[[1]]
    populationQ<- as.matrix(verficaSingletonsQ[[2]])
    populationQ<-calculaJerarquiasDensidad(popSize, populationQ, tablaGruposQ, param$objDim, FALSE, distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, num_k,alfa,varNo, param$nombre_objetivos)
    ######################## Population R   ########################
    Log("[Agent %d] Generating Population R for generation %d\n", agente, g)
    populationR <- rbind(populationP,populationQ)
    rownames(populationR)<-1:nrow(populationR)
    
    Log("[Agent %d] Computing pareto ranking for population R of generation %d\n", agente, g)
    populationR<-populationR[, -((num_k+1):(varNo+param$objDim+2))] #Mantain only chromosomes, I mean, only integer without others columns
    tablaGruposR<- generaGrupos(nrow(populationR), populationR, distancias$distancia_exp, distancias$distancia_bio, alfa)
    #Re-compute Pareto ranking and crowding in population R
    populationR<-calculaJerarquiasDensidad(popSize*2, populationR, tablaGruposR, param$objDim, FALSE, distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, num_k, alfa, varNo, param$nombre_objetivos)
    populationR<-as.data.frame(populationR)
    poblacion_pareto<-subset(populationR, rnkIndex=='1')
    poblacion_pareto<-poblacion_pareto[!duplicated(poblacion_pareto[,1:(num_k)]),]
    #Delete solutions with singletons or groups with no gene.
    tablaGruposPareto<-generaGrupos(nrow(poblacion_pareto), as.matrix(poblacion_pareto[,1:num_k]), distancias$distancia_exp, distancias$distancia_bio, alfa)
    arreglar<-quitarSolucionesConSingletons(tablaGruposPareto, poblacion_pareto, num_k)
    poblacion_pareto <- arreglar$poblacion # Se cambio = por <-
    populationR<-poblacion_pareto
    
    #Path relinking and Pareto Local Search wont be used
    
    if(experiment == "moc"){
      Log("Performing path relinking for generation %d\n",g)
      #Perform Path relinking
      populationR<-generaPathRelinking(poblacion_pareto, alfa,num_k, distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, param$objDim, varNo, param$nombre_objetivos) #Return population R but only Pareto solutions without duplicated solutions
      Log("Performing pareto local search for generation  %d\n",g)
      #Perform Pareto Local Search")
      populationR<-generaParetoLocalSearch(populationR, param$vecinos, alfa, num_k, distancias$numgenes, popSize, distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, param$objDim, varNo, param$nombre_objetivos) #Input only Pareto (I mean, solutions in Population R after PR procedure) it includes Rnk y Densidad
      
      local_evaluations <- floor(96*distancias$numgenes*p^3 + p^4)
      ev <- ev + local_evaluations
      Log("Adding %d evaluations by using local search", local_evaluations)
      Log("Current evaluation count: %d out of %d (%f)", ev, evaluations, round(ev/evaluations, 2))
    }
    
    ######################################################################
    #                     Verify convergence
    ######################################################################
    
    Log("[Agent %d] Verifying convergence for generation %d\n", agente, g)
    
    #combined <- rbind(pareto, pareto)
    #new <- rep("new", nrow(pareto))
    #old <- rep("old", nrow(pareto))
    #gen <- data.frame(gen=c(new, old))
    #combined <- cbind(combined, gen)
    dominated_new_old <- dominated_solutions(poblacion_pareto, last_pareto)
    Log("[Agent %d] New generation has %d dominated solutions by previous generation", agente, dominated_new_old)
    dominated_old_new <- dominated_solutions(last_pareto, poblacion_pareto)
    Log("[Agent %d] Previous generation has %d dominated solutions by the new generation", agente, dominated_old_new)
    
    convergence_index <- dominated_old_new/nrow(last_pareto) - dominated_new_old/nrow(poblacion_pareto)
    Log("[Agent %d] Convergence index is %f", agente, round(convergence_index, 2))
    
    # A break is added because a refill with random solutions is not wanted (will worsen quality)
    
    if(convergence_index < param$tol && convergence_index >= 0){
        conv_count = conv_count + 1
        Log("[Agent %d] Convergence count is now %d", agente, conv_count)
        if(conv_count == param$convergencia){
          Log("[Agent %d] Reached convergence in generation %d", agente, g)
          break
        }
    }else{
      conv_count = 0 #Reset convergence if a change happened
    }
    
    # Save the last population P to compare with a new population P+1 for convergence 
    last_pareto <- poblacion_pareto
    
    if(g == param$generation){
      Log("[Agent %d] Maximum generations reached, returning...", agente)
      break
    }
    ######################################################################
    Log("[Agent %d] Population p+1 for generation %d\n", agente, g)
    ######################## Population P+1  ########################
    if(nrow(poblacion_pareto) < popSize) {
      
      populationP<- as.matrix(poblacion_pareto[,1:num_k])
      
      fill_num <- popSize - nrow(poblacion_pareto)
      Log("[Agent %d] Refilling %d populations for generation %d:\n",agente, fill_num , g)
      random_solutions <- fill_num + ceiling(popSize/3)   # Adding a surplus, just in case of duplication
      
      #Refill Population with PR
      population_random<-as.data.frame(t(sapply(1:random_solutions, function(u) array(sample(1:distancias$numgenes, num_k, replace=F)))))
      #Verify singletons
      reparar<-repararSingletons(as.matrix(population_random), distancias$distancia_exp, distancias$distancia_bio, distancias$matriz_exp, distancias$numgenes, alfa)
      population_random=reparar[[2]] #2 Cause, first are the clusters
      
      #Only a chromosome, I mean no repeated.
      if(nrow(populationP)>1){
        populationP<-populationP[!duplicated(populationP[,1:(num_k)]),]
        populationP<-as.data.frame(populationP)
      }
      
      populationP<-rbind(populationP[,1:(num_k)], population_random[1:fill_num, ])
      rownames(populationP)<-c(1:nrow(populationP))
      populationP<-as.matrix(populationP)
      
    }else{
      populationP<- as.matrix(poblacion_pareto[1:popSize,1:num_k])
    }
    
    Log("[Agent %d] Generation %d finished\n", agente, g)
    g <- g + 1
    
    local_evaluations <- floor(36*distancias$numgenes*p + 4*p^2 + 5*p)
    ev <- ev + local_evaluations
    Log("[Agent %d] Adding %d evaluations in this generation", agente, local_evaluations)
    Log("[Agent %d] Current evaluation count: %d out of %d (%f)", agente, ev, evaluations, round(ev/evaluations, 2))
 
  }
  
  
  ################# Convergence Analisis
  #Log("[Phase %d] Verifying convergence for every agent\n", phase)
  #  
  #  siluetas = carga_visualiza_silueta(Agents[[i]], distancias$matriz_exp, distancias$distancia_exp, distancias$distancia_bio, alfa, num_k)
  #  siluetas = siluetas$estadistica[, 1:num_k] # POR ALGUNA RAZON ENTREGA UNA COLUMNA DE MAS
  #  convergencia<-rbind(convergencia[[i]], siluetas)
  #  convergencia_silueta <- sum(convergencia[(g+1), ] - convergencia[g, ])
    
  #  Log("[Agent %d] Silhouette variation is %f at generation %d\n", agente, convergencia_silueta, g)
    
   # print(convergencia)
    
    #convergence_table <- data.frame(Generation=c(1:nrow(convergencia)), Silhouette.Mean=rowMeans(convergencia))
    
    #convergence_plot <- ggplot(convergence_table, aes(x=Generation, y=Silhouette.Mean), xlab="Generation", ylab="Avg. Silhouette index") +
    #  ggtitle(paste("Convergence for agent", agente)) +
    #  geom_line() +
    #  geom_point()
    
    #plot(convergence_plot)
  
  #################
  #Obtain results
  if(experiment != "memetic"){
    guarda_paretos=carga_guarda_pareto(alfa, num_k,  distancias$matriz_exp,  distancias$distancia_exp, distancias$distancia_bio, popSize, poblacion_pareto, FALSE)
    
    #File result
    exp.name <- paste0("<root directory here>\\NSGA2\\MOC-GaPBK\\Tests\\", test_name ,experiment)
    
    #save results
    write.table(guarda_paretos$pareto_clustering, file = paste0(exp.name,"_clustering.csv"),sep = ",",col.names = TRUE, quote = FALSE)
    write.table(guarda_paretos$pareto_objectives, file = paste0(exp.name,".csv"),sep = ",",col.names = TRUE, row.names = TRUE,quote = FALSE)
    write.table(guarda_paretos$pareto_silhouette, file = paste0(exp.name,"_silhouettes.csv"),sep = ",",col.names = FALSE, quote = FALSE)
    Log("The pareto solutions are:\n")
    print(poblacion_pareto)
  }
  
  return(poblacion_pareto)
}

