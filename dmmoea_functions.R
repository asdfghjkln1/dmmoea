generate.initial.pop <- function(params, k, num.genes){
  
  if (params$is_random_population==TRUE) {
    
    # Generate random population
    set.seed(params$seed) # set seed
    population <- t(sapply(1:params$popSize, function(x)  sample(1:num.genes, k ,replace=F)  ))
    #rm(.Random.seed, envir=globalenv()) ?
    return(as.matrix(population))
  }else{
    # Generate population using single objective solutions
    #Add code to read file with single objective solutions
  }
  
}

cluster_data <- function(distances, population, alpha){
  K <- ncol(population)
  pop.size <- nrow(population)
  bio.dist <- distances$bio.dist
  exp.dist <- distances$exp.dist
  gene.len <- distances$n.genes
  genes <- 1:gene.len
  
  clustering.results <- list(1:pop.size)
  grouping<-rep(NA,gene.len)
  
  p <- 1
  while(p < pop.size){
    for(gene in 1:gene.len){
      medoid.distances <- alpha*exp.dist[gene, population[p, 1:K]] + (1-alpha)*bio.dist[gene, population[p, 1:K]]
      gene.clust <- unname(which.min(medoid.distances))
      grouping[gene] <- gene.clust
    }
    clustering.results[[p]] <- grouping
    names(clustering.results[[p]])<-rownames(exp.dist)
    
    col <- 1
    # For every cluster, check for singletons. 
    # Replace medoid until all clusters in the solution has more than 2 elements
    while(col < K){
      if(length(which(unlist(grouping) == col)) < 3){
        print("A singleton was found!!")
        medoides.solution <- clustering.results[p, 1:K]
        reemplazo <- setdiff(genes, intersect(genes, medoides.solution))
        clustering.results[p, col] <- sample(reemplazo, 1)
        col <- 1 # Check for singletons in previous medoids, until every medoid has no singletons
      }else{
        col <- col + 1
      }
    }
    
    p <- p + 1
  }
  
  return(clustering.results)
}

evaluate_population <- function(population, cluster_results, distances, params){
  
  objectives <- params$objectives
  obj.dim <- length(objectives)
  obj1 <- objectives[1]
  obj2 <- objectives[2]
  
  if(obj1 == "XieBeni"){
    f1 <- evaluate_xie_beni(population, cluster_results, distances$exp.dist)
  }else if(obj1 != "XieBeni"){
    #TO DO implement
    print("Objective function not available!")
    return(NULL)
  }
  
  if(obj2 == "XieBeni"){
    f2 <- evaluate_xie_beni(population, cluster_results, distances$bio.dist)
  }else if(obj1 != "XieBeni"){
    #TO DO implement
    print("Objective function not available!")
    return(NULL)
  }
  print("Population before sorting:")
  print(population)
  obj.values <- cbind(f1,f2)
  population <- cbind(population, obj.values)
  
  ranking <- nsga2R::fastNonDominatedSorting(obj.values)
  rnkIndex <- integer(params$popSize)
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  } 
  population <- cbind(population, rnkIndex)
  
  # Calculate Crowding Distance
  objRange <- apply(obj.values, 2, max) - apply(obj.values, 2, min)
  cd <- nsga2R::crowdingDist4frnt(population, ranking, objRange) # Calculate CD
  cd.density <- apply(cd,1,sum)
  population <- cbind(population, cd.density) # Add cd density column
  population <- population[order(rnkIndex, -cd.density), ] # Order by cd
  
  print("Population after sorting:")
  print(population)
  return(as.data.frame(population))
}

# Xie-Beni clustering coefficient evaluation
# Evaluates a single objective
evaluate_xie_beni <- function(population, cluster_results, distances){
  
  obj.function <-rep(0, nrow(population))
  
  for (p in 1:nrow(population)) {
    
    # Initialize matrix
    n <- nrow(distances) # number of elements
    K <- ncol(population) # number of medoids
    D <- matrix(0, K, n) #Matriz de distancia centros x elementos
    
    # Calculate square distances of each element to its medoid
    for (k in 1:K) {
      for (i in 1:n) {
        cluster.medoid=population[p, k]
        D[k, i] = distances[cluster.medoid, i]^2
      }
    }
    
    XB.numerator = sum(D)
  
    combs <- t(combn(population[p, ], 2)) # Generate pairs of medoids
    s <- vector()
    for (pair in 1:nrow(combs)) {
      s[pair]= distances[combs[pair, 1], combs[pair, 2]]^2
    }
    
    XB.denominator <- min(s)

    obj.function[[p]]=XB.numerator/(n*XB.denominator)
  }
  return(unlist(obj.function))
}

nsga2 <- function(distances, initial_pop, K, alpha, params){
  P <- initial_pop
  g <- 0 # Current generation
  while( g <= params$generation ){
    
    clustering.groups <- cluster_data(distances, P, alpha)
    print("Groups set!")
    P <- evaluate_population(P, clustering.groups, distances, params)
    
    print("Population evaluated!")
    
    #if(nrow(P) < 2){
    #  print("ERROR: Population has length 1")
    #  return(NULL)
    #  #TO DO: Implementar Reparar singletons
    #}
    
    g <- g + 1
    return(P)
  }
  
}



#####
MemeticNSGA2 <- function(pareto_population, agent_id, pop_size, evaluations, param, current_phase){
  
  library("doParallel")
  source("MOC-GaPBK_nuevas_funciones.R")
  
  # (Number of) agents and pop_size refers of the second level depth
  phaseslv2 <- param$phaseslv2
  agents <- param$agentslv2
  
  cl <- makeCluster(agents)
  
  registerDoParallel(cl)
  
  phase=1
  num_k <- param$serie_k
  
  #This will contain a list of each agent population
  Agents <- list()
  
  exported_func <- c("NSGA2", "Log")
  clusterExport(cl, exported_func)
  
  agent_output_file <- paste("agent_", agent_id, ".txt", sep="")
  
  
  sink(agent_output_file, append=TRUE)
  if(current_phase == 1){ # Not doing this would erase the previous phase logs when called again...
    writeLines(c(""), agent_output_file)
  }
  
  # If theres only one solution, next phase will fail (one or more agents without solutions), so random refill is needed
  if(nrow(pareto_population) < agents){
    
    source("MOC-GaPBK_funciones.R")
    carga_distancias=dget("MOC-GaPBK_load_distancias.R")
    distancias=carga_distancias(param$path)
    
    solutions_needed= agents-nrow(pareto_population)
    
    population_random<-as.data.frame(t(sapply(1:solutions_needed, function(u) array(sample(1:distancias$numgenes, num_k, replace=F)))))
    #Verify singletons
    reparar<-repararSingletons(as.matrix(population_random), distancias$distancia_exp, distancias$distancia_bio, distancias$matriz_exp, distancias$numgenes, alfa)
    population_random=reparar[[2]] #2 Cause, first are the clusters
    pareto_population <- rbind(pareto_population, population_random)
    #rownames(pareto_population)<-c(1:nrow(pareto_population))
  }
  
  lv1population <- pop_size
  lv2population <- floor(pop_size / agents)
  
  #Splitting solutions to asign an agent. It shouldnt be random otherwise agent loses explotation
  if(nrow(pareto_population) < lv1population){
    num_per_agent <- floor(nrow(pareto_population) / agents)
    solutions <- rep(1:agents, each=num_per_agent)  
  }else{
    solutions <- rep(1:agents, each=lv2population)
  }
  
  solutions <- split(1:nrow(pareto_population), solutions)
  
  for(i in 1:agents){
    Agents <- c(Agents, list( pareto_population[ unlist(solutions[i]) , ] ) )
  }
  
  #Number of maximum evaluations per phase per agent
  sync_evaluations <- phaseslv2*agents^2*lv2population^2
  Log("Sync evaluation overhead is %d", sync_evaluations)
  ev_per_agent <- floor(evaluations/(phaseslv2*agents)) - sync_evaluations
  Log("Number of evaluations per agent per lvl2 phase: %d", ev_per_agent)
  
  exported_vars <- c("agent_output_file", "Agents", "lv2population", "ev_per_agent")
  
  current_evaluations = 0
  
  while (phase <= phaseslv2) {
    # Por cada agente i, generar soluciones
    Log("[Phase %d] Population for level 2 phase %d:\n", current_phase, phase)
    print(Agents[[1]])
    print(Agents[[2]])
    
    Agents <- foreach(i=1:agents, .combine = c , .export=exported_vars, .inorder = FALSE) %dopar% {
      
      Log("[Agent Level 1] Initializing agent level 3\n")
      poblacion_pareto <- as.matrix(Agents[[i]])
      poblacion_pareto <- NSGA2(poblacion_pareto, i, lv2population, ev_per_agent, agent_output_file, "NA", "memetic") #parameter "NA" will not be used with memetic
      
      return( list(poblacion_pareto) ) 
    }
    
    Log("[Agent Level 1] Agents finished\n")
    Log("Printing results:\n")
    print(Agents)
    
    # Synchronize agents
    for(i in 1:agents){
      for(j in 1:agents){
        if(i == j) {
          next
        }
        Agents[[i]] <- fitness(Agents[[i]], Agents[[j]], "min", "min", lv2population)
      }
    }
    Log("[Phase %d] Level 2 phase %d synchronization ended!\n", current_phase, phase)
    
    #We need to discard the other information before using it for the new phase population
    #After the final phase, return the complete results
    if(phase < phaseslv2){
      for(i in 1:agents){
        Agents[[i]] <- Agents[[i]][ , 1:num_k]
      }
    }
    
    phase <- phase + 1
    
  }
  Log("[Agent Level 1] Memetic algorithm finished!\n")
  stopCluster(cl)
  
  Log("These are the raw results:\n")
  print(Agents)
  
  pareto <- data.frame()
  for(i in 1:agents){
    pareto <- rbind(pareto, Agents[[i]])
  }
  
  pareto <- pareto[!duplicated(pareto),] # Delete duplicates
  
  # Pruning non pareto-front solutions
  objetives <- as.matrix(pareto[, c("obj1", "obj2")]) 
  ranks <- NonDomSorting(objetives) 
  pareto <- pareto[ unlist(ranks[[1]]) , ] # Select the rank 1 indices
  
  if(nrow(pareto) > lv1population){ # Cut if theres too much population
    pareto <- pareto[1:lv1population, ]
  }
  
  Log("These are the final results:\n")
  print(pareto)
  
  return(pareto)
}
#####

#####
plotPareto <- function(dataset_name, exp_name ,ParetoMemetic, ParetoMOCGaPBK, ParetoNGSA2){
  
  A <- as.data.frame(ParetoMemetic)
  B <- as.data.frame(ParetoMOCGaPBK)
  C <- as.data.frame(ParetoNGSA2)
  A <- A[, c("obj1", "obj2")]
  A <- cbind(A, exp=rep("A", nrow(A)))
  B <- B[, c("obj1", "obj2")]
  B <- cbind(B, exp=rep("B", nrow(B)))
  C <- C[, c("obj1", "obj2")]
  C <- cbind(C, exp=rep("C", nrow(C)))
  S <- rbind(A,B,C)
  S[,1] <- (S[,1]-min(S[,1]))/(max(S[,1])-min(S[,1]))
  S[,2] <- (S[,2]-min(S[,2]))/(max(S[,2])-min(S[,2]))
  colnames(S) <- c("obj1", "obj2", "method")
  
  pareto_plot <- ggplot(S, aes(x=obj1, y=obj2, colour=method)) +
    labs(title=dataset_name, color="Algorithm", tag = exp_name, x="Expression Index", y="Biological Index") +
    geom_line(aes(color = method)) +
    geom_point(size=2)+
    scale_color_manual(name = "Algorithm", values=c("red", "green", "blue"), labels = c("DMNGSA-2", "MOC-GaPBK", "NGSA-2"))
  
  plot(pareto_plot)
}
#####

#####
plot_convergence <- function(){
  
  #install.packages("gridExtra")
  library(gridExtra)
  
  # Test A
  label <- rep( c("1", "2"), each=6)
  phases <- rep( 1:6, 2)
  
  # Test B
  #label <- rep( c("1", "2"), each=10)
  #phases <- rep(1:10, 2)
  
  # Test A
  agent1 <- c(25, 12, 10, 11, 10, 10, 22, 12, 10, 11, 10, 10)
  agent2 <- c(15, 14, 16, 11, 10, 10, 19, 20, 21, 16, 11, 10)
  agent3 <- c(22, 12, 10, 11, 10, 10, 28, 18, 12, 20, 11, 10)
  agent4 <- c(36, 11, 10, 10, 10, 10, 21, 11, 31, 15, 10, 10)
  
  a1 <- data.frame(phases,agent1, label)
  a2 <- data.frame(phases,agent2, label)
  a3 <- data.frame(phases,agent3, label)
  a4 <- data.frame(phases,agent4, label)
  # Test B
  #a5 <- data.frame(phases,agent5, label)
  
  # Test A
  phase <- 1:6
  phase.2 <- 3
  rows<-2
  # Test B
  #phase <- c(1:10)
  #phase.2 <- 5
  #rows<-3
  
  conv_plot1 <- ggplot(data=a1, aes(x=phases, y=agent1, col=label, title="Agent 1 (lv1)"), title="Agent 1 (lv1)") +
    labs(x="Phase", y="Convergence") +
    geom_line() +
    scale_x_continuous(breaks = phase)+
    ggtitle("Agent 1 (lv1)") + 
    theme(legend.position = "none") +
    geom_vline(xintercept = phase.2, linetype=2)
  #geom_line(aes(color = label)) +
  conv_plot2 <- ggplot(data=a2, aes(x=phases, y=agent2, col=label)) +
    labs(x="Phase", y="Convergence") +
    geom_line() + 
    scale_x_continuous(breaks = phase)+
    ggtitle("Agent 2 (lv1)") + 
    theme(legend.position = "none") +
    geom_vline(xintercept = phase.2, linetype=2)
  
  conv_plot3 <- ggplot(data=a3, aes(x=phases, y=agent3, col=label)) +
    labs(x="Phase", y="Convergence") +
    geom_line() + 
    ggtitle("Agent 3 (lv1)") + 
    scale_x_continuous(breaks = phase) +
    theme(legend.position = "none") +
    geom_vline(xintercept = phase.2, linetype=2)
  
  conv_plot4 <- ggplot(data=a4, aes(x=phases, y=agent4, col=label)) +
    labs(x="Phase", y="Convergence") +
    geom_line() + 
    ggtitle("Agent 4 (lv1)") + 
    scale_x_continuous(breaks = phase) +
    theme(legend.position = "none") +
    geom_vline(xintercept = phase.2, linetype=2)
  # Test B
  conv_plot5 <- ggplot(data=a5, aes(x=phases, y=agent5, col=label)) +
    labs(x="Phase", y="Convergence") +
    geom_line() + 
    scale_x_continuous(breaks = phase) +
    ggtitle("Agent 5 (lv1)") + 
    theme(legend.position = "none") +
    geom_vline(xintercept = phase.2, linetype=2)
  
  grid.arrange(conv_plot1, conv_plot2, conv_plot3, conv_plot4, nrow=rows)
  #grid.arrange(conv_plot1, conv_plot2, conv_plot3, conv_plot4, conv_plot5, nrow=rows)
}
#####

