generate_initial_pop <- function(P.size, K, num.genes, seed, is_random_population){
  
  if (is_random_population) {
    
    # Generate random population
    set.seed(seed) # set seed
    population <- t(sapply(1:P.size, function(x)  sample(1:num.genes, K ,replace=F)  ))
    #rm(.Random.seed, envir=globalenv()) ?
    return(as.matrix(population))
  }else{
    # Generate population using single objective solutions
    #Add code to read file with single objective solutions
  }
  
}

generate_diverse_initial_pop <- function(distances, params, K, p.size=NULL, diverse_population=FALSE){
  
  n.genes <- distances$n.genes
  if(missing(p.size)){
    p.size <- params$popSize
  }
  # Normalize distance matrix
  D <- distances$exp.dist # params$alpha*distances$bio.dist + (1-params$alpha)*distances$exp.dist # Doesn't work? Weird numbers
  D <- as.data.frame(eaf::normalise(D, c(0,1)))
  D.size <- nrow(D)
  
  population <- matrix(nrow=p.size, ncol=K)
  
  # Create a chromosome (single element of population) for every radius of density
  # Set density tolerance T (number of elem. within density radius)
  if(params$auto_adjust_initial_params){
    T <- n.genes/K
    radius <- T/n.genes
    min.r <- max(0.1, radius*0.5)
    max.r <- min(0.5, radius*1.5)
    r.series <- seq(min.r, max.r, by=(max.r - min.r)/p.size)
  }else{
    T <- n.genes * params$density_tol # *** Might need to be adjusted ***
    r.series <- seq(params$min_radius,params$max_radius, by=(params$max_radius-params$min_radius)/p.size)
  }
  used.in.population <- list()
  for(p in 1:p.size){
    used.in.chromosome <- rep(FALSE, D.size) # Keeps track of discarded (used in cluster) elements
    candidates <- sample(1:D.size, D.size, replace = FALSE)
    candidate.search.index <- 1
    for(k in 1:K){
      gene.selected <- FALSE
      densities <- rep(0, D.size) 
      for(i in candidate.search.index:length(candidates)){
        g <- candidates[i]
        candidate.search.index <- candidate.search.index + 1 # Avoids repeating search
        if(used.in.chromosome[g] == TRUE){
          next
        }
        in.cluster <- D[g, ] <= r.series[p]
        ro <- sum(in.cluster) - 1
        densities[g] <- ro
        # IF density is enough, select as medoid and stop searching
        if(ro > T){
          if(diverse_population){
            # Check first if >80% of population is marked as "used" in population, reset restrictions for diversity
            # Otherwise, not enough elements available to use as medoids. Only used when diverse_population=TRUE
            if(length(used.in.population)/D.size > 0.8){
              used.in.population <- unlist(population[p, 1:max(k-1, 1)]) # Clean population register of elements used
            }
            if(g %in% used.in.population){
              next # Discard candidate if it was already used in previous chromosomes 
            }
          }
          #print(paste("Dense medoid found at trial", candidate.search.index, "using", g))
          population[p, k] <- g
          gene.selected <- TRUE
          break
        }
      }
      # Since no dense enough medoid was found (search was exhausted), fill the chromosome with best options
      if(!gene.selected){
        # *** Change this if using dynamic K ****
        # Fill remaining medoids with the highest density candidates, because no genes are dense enough
        candidates.left <- order(densities, decreasing = TRUE)
        index <- 1
        for(m in k:K){
          dense.medoid <- candidates.left[index]
          used <- used.in.chromosome[dense.medoid]
          while(used){
            dense.medoid <- sample(1:D.size, 1)
            index <- index + 1
            used <- used.in.chromosome[dense.medoid]
            if(diverse_population){ 
              # Check again if elements used in population are too many, if so, reset
              if(length(used.in.population)/D.size > 0.8){
                used.in.population <- unlist(population[p, 1:max(k-1, 1)]) # Clean population register of elements used
              }
              used <- dense.medoid %in% used.in.population
            }
          }
          population[p, m] <- dense.medoid
          # Update used genes
          in.cluster <- D[dense.medoid, ] <= r.series[p]
          used.in.chromosome <- used.in.chromosome | in.cluster # Add used elements
          # Update used genes in population
          used.in.population <- unlist(append(used.in.population, which(in.cluster)))
          
          index <- index + 1
        }
        break # Exit for since all medoids have been filled
      }
      # Records used elements from this chromosome
      used.in.chromosome <- used.in.chromosome | in.cluster # Add used elements
      # Record used elements in population to increase selection diversity
      used.in.population <- unlist(append(used.in.population, which(in.cluster)))
    }
  }
  return(population)
}

cluster_data <- function(distances, population, alpha){
  K <- ncol(population)
  pop.size <- nrow(population)
  bio.dist <- distances$bio.dist
  exp.dist <- distances$exp.dist
  #dist.compounded <- alpha*bio.dist + (1-alpha)*exp.dist
  #dist.compounded <- as.data.frame(dist.compounded)
  gene.len <- distances$n.genes
  genes <- 1:gene.len
  clustering.results <- list(1:pop.size)
  grouping<-rep(NA,gene.len)
  p <- 1
  while(p <= pop.size){
    for(gene in 1:gene.len){
      medoid.distances <- alpha*exp.dist[gene, population[p, 1:K]] + (1-alpha)*bio.dist[gene, population[p, 1:K]]
      #medoid.distances <- dist.compounded[gene, population[p, 1:K]] # Doesn't work this way, not sure why
      # Note: sometimes all distances are 1, it chooses always index 1
      gene.clust <- unname(which.min(medoid.distances))
      grouping[gene] <- gene.clust
      
    }
    clustering.results[[p]] <- grouping
    names(clustering.results[[p]])<-rownames(exp.dist)
    
    col <- 1
    # For every cluster, check for singletons. 
    # Replace medoid until all clusters in the solution has more than 2 elements
    while(col <= K){
      if(length(which(unlist(grouping) == col)) < 3){
        replacement <- setdiff(genes, intersect(genes, population[p, 1:K]))
        population[p, col] <- sample(replacement, 1)
        new.medoids <- population[p, 1:K]
        #** Update groups for new medoids
        for(gene in 1:gene.len){
          medoid.distances <- alpha*exp.dist[gene, new.medoids ] + (1-alpha)*bio.dist[gene, new.medoids]
          gene.clust <- unname(which.min(medoid.distances))
          grouping[gene] <- gene.clust
        }
        clustering.results[[p]] <- grouping
        names(clustering.results[[p]])<-rownames(exp.dist)
        col <- 1 # Check for singletons in previous medoids, until every medoid has no singletons
      }else{
        col <- col + 1
      }
    }
    p <- p + 1
  }
  return(list("population"=population, "clustering.results"=clustering.results))
}

evaluate_population <- function(population, cluster_results, distances, params){
  
  objectives <- params$objectives
  obj.dim <- length(objectives)
  obj1 <- objectives[1]
  obj2 <- objectives[2]
  
  p <- length(cluster_results) # not used atm
  
  if(obj1 == "XieBeni"){
    f1 <- evaluate_xie_beni(population, distances$exp.dist)
  }else if(obj1 != "XieBeni"){
    #TO DO implement
    print("Objective function not available!")
    return(NULL)
  }
  
  if(obj2 == "XieBeni"){
    f2 <- evaluate_xie_beni(population, distances$bio.dist)
  }else if(obj1 != "XieBeni"){
    #TO DO implement
    print("Objective function not available!")
    return(NULL)
  }
  
  obj.values <- cbind(f1,f2)
  population <- dominance_ranking_sorting(population, obj.values)
  
  return(population)
}

dominance_ranking_sorting <- function(population, obj.values){
  
  ranking <- nsga2R::fastNonDominatedSorting(obj.values)
  rnkIndex <- integer(nrow(population))
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  } 
  population <- cbind(population, obj.values, rnkIndex)
  
  # Calculate Crowding Distance
  objRange <- apply(obj.values, 2, max) - apply(obj.values, 2, min)
  cd <- nsga2R::crowdingDist4frnt(population, ranking, objRange) # Calculate CD
  cd.density <- apply(cd,1,sum)
  population <- cbind(population, cd.density) # Add cd density column
  population <- population[order(rnkIndex, -cd.density), ] # Order by cd
  return(as.data.frame(population))
}

# Xie-Beni clustering coefficient evaluation
# Evaluates a single objective
evaluate_xie_beni <- function(population, distances){
  
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

#### NSGA-2 ####
nsga2 <- function(distances, K, params, output.path){
  
  # Initialize Random population
  P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed, params$is_random_population) 
  P.size <- nrow(P) # Population size
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, P.clustering.groups, distances, params)
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  
  generations.no.changes <- 0
  has.converged <- FALSE
  
  ## Main generation loop
  while(!has.converged){
    
    print(paste("Entering selection, crossover, mutation for generation ", g))
    
    ###### Selection, Crossover, Mutation  ######
    
    ## Selection
    mating_pool <- tournamentSelection(P, P.size, params$tourSize)
    
    ##Crossover and Mutation
    Q <- population_mating_and_mutation(mating_pool, K, num.genes, params)
    
    ## Evaluate solutions in population
    Q.data <- cluster_data(distances, Q, params$alpha)
    Q <- Q.data$population
    Q.clustering.groups <- Q.data$clustering.results
    Q <- evaluate_population(Q, Q.clustering.groups, distances, params)
    Q.clustering.groups <- Q.clustering.groups[order(as.numeric(rownames(Q)))]
    rownames(Q) <- 1:nrow(Q)
    
    ## Create population R as the combination of P and Q
    R <- rbind(P,Q)
    R.clustering.groups <- c(P.clustering.groups, Q.clustering.groups)
    dup.index <- duplicated(R[, 1:K])
    R <- as.matrix(R[!dup.index, 1:(K+params$objDim)])
    R.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(R)))] # Update clustering
    
    obj.values <- R[, (K+1):ncol(R)] # Select objective values
    R <- dominance_ranking_sorting(R[, 1:K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(R)))] # Update clustering
    rownames(R)<-1:nrow(R)
    
    ## Population fitness selection
    P_next_generation <- fitness_selection_crowding_distance(R, P.size, K) 
    P.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(P_next_generation)))] # Update clustering
    rownames(P_next_generation) <- 1:nrow(P_next_generation)
    new.pareto.front <- P_next_generation[P_next_generation$rnkIndex == 1, ]
    pareto.clustering <- P.clustering.groups[ order(as.numeric(rownames(new.pareto.front))) ]
    
    ## Output pareto front plot
    plot_pareto(current.pareto.front, new.pareto.front, g, output.path) # Output pareto front
    
    ## Measure convergence of pareto front
    convergence.index <- convergence_coefficient(current.pareto.front, new.pareto.front, g)
    
    ## Check how different is the new pareto front, count generations with no changes
    if(convergence.index <= params$convergence.tol){
      generations.no.changes <- generations.no.changes + 1 
    }else{
      generations.no.changes <- 0 # Reset counter if changes detected
    }
    ## Check for stop criteria
    if(generations.no.changes > params$convergence.limit || g >= params$generation){
      has.converged <- TRUE
    }
    ## Continue to next generation
    g <- g + 1
    P <- P_next_generation
    current.pareto.front <- new.pareto.front
  }
  
  return(list("population"=P_next_generation, "clustering"=P.clustering.groups))
}

#### Diverse NSGA-2 ####
dnsga2 <- function(distances, K, diversity.metric, diversity.level, params, output.path, generations=NULL, P.size=NULL, agent=NULL, initial_population=NULL){
  
  if(missing(agent)){
    print(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"..."))
  }else{
    log.file <- file.path(output.path, "log.txt")
    sink(log.file, append=TRUE) # Register in log output file
    source("dmmoea_libraries.R")
    source("dmmoea_functions.R") # Load function workspace
    Log(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"..."), agent=agent)
  }
  
  # Initialize population
  if(!missing(initial_population) && !is.null(initial_population)){
    P <- as.matrix(initial_population)[, 1:K]
    #*** TODO: Check initial population validity as input **** #
  }else if(params$is_random_population){
    # Random population
    P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed, params$is_random_population) 
  }else if(diversity.level >= 1){
    P <- generate_diverse_initial_pop(distances, params, K, diverse_population = TRUE)
  }else{
    P <- generate_diverse_initial_pop(distances, params, K, diverse_population = FALSE)
  }
  
  if(missing(P.size)){
    P.size <- params$popSize # Population size
  } 
  if(missing(generations)){
    generations <- params$generation # Maximum generations
  }
  
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  
  Log("A", agent=agent)
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, P.clustering.groups, distances, params)
  P.clustering.groups <- P.clustering.groups[order(as.numeric(rownames(P)))]
  rownames(P) <- 1:nrow(P)

  Log("B", agent=agent)
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  
  generations.no.changes <- 0
  has.converged <- FALSE
  
  ## Main generation loop
  while(!has.converged){
    
    if(missing(agent)){
      print(paste("Entering selection, crossover, mutation for generation ", g))
    }else{
      Log(paste("Entering selection, crossover, mutation for generation ", g), agent=agent)
    }
    
    ###### Selection, Crossover, Mutation  ######
    
    ## Selection
    mating_pool <- tournamentSelection(P, P.size, params$tourSize)
    
    ##Crossover and Mutation
    if(diversity.level >= 3){
      Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, K, params, P.size=P.size) 
    }else{
      Q <- population_mating_and_mutation(mating_pool, K, num.genes, params, P.size=P.size)
    }
    
    ## Evaluate solutions in population
    Q.data <- cluster_data(distances, Q, params$alpha)
    Q <- Q.data$population
    Q.clustering.groups <- Q.data$clustering.results
    
    Q <- evaluate_population(Q, Q.clustering.groups, distances, params)
    Q.clustering.groups <- Q.clustering.groups[order(as.numeric(rownames(Q)))]
    rownames(Q) <- 1:nrow(Q)
    
    ## Create population R as the combination of P and Q
    R <- rbind(P,Q)
    R.clustering.groups <- c(P.clustering.groups, Q.clustering.groups)
    dup.index <- duplicated(R[, 1:K])
    R <- as.matrix(R[!dup.index, 1:(K+params$objDim)])
    R.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(R)))] # Update clustering

    obj.values <- R[, (K+1):ncol(R)] # Select objective values
    R <- dominance_ranking_sorting(R[, 1:K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(R)))] # Update clustering
    rownames(R)<-1:nrow(R)
  
    ## Population fitness selection
    if(diversity.level >= 2){
      P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, P.size, K, diversity.metric)
    }else{
      P_next_generation <- fitness_selection_crowding_distance(R, P.size, K) 
    }
    P.clustering.groups <- R.clustering.groups[order(as.numeric(rownames(P_next_generation)))] # Update clustering
    rownames(P_next_generation) <- 1:nrow(P_next_generation)
    new.pareto.front <- P_next_generation[P_next_generation$rnkIndex == 1, ]
    pareto.clustering <- P.clustering.groups[ order(as.numeric(rownames(new.pareto.front))) ]
    
    if(missing(agent)){
      plot_pareto(current.pareto.front, new.pareto.front, g, output.path) # Output pareto front
    }else{
      plot_pareto(current.pareto.front, new.pareto.front, g, output.path, agent=agent) # Output pareto front
    }
    
    ## Measure convergence of pareto front
    convergence.index <- convergence_coefficient(current.pareto.front, new.pareto.front, g)
    
    ## Check how different is the new pareto front, count generations with no changes
    if(convergence.index <= params$convergence.tol){
      generations.no.changes <- generations.no.changes + 1 
    }else{
      generations.no.changes <- 0 # Reset counter if changes detected
    }
    ## Check for stop criteria
    if(generations.no.changes > params$convergence.limit || g >= generations){
      has.converged <- TRUE
    }
    
    ## Continue to next generation
    g <- g + 1
    P <- P_next_generation
    current.pareto.front <- new.pareto.front
    
    Log("C", agent=agent)
  }
  
  Log("Finished dnsga2!", agent)
  print(P_next_generation)
  
  return(list("population"=P_next_generation, "clustering"=P.clustering.groups))
}

convergence_coefficient <- function(current.pareto, new.pareto, generation){
  if(generation == 0){
    return(FALSE) # Cant converge at first generation
  }
  old.dominated <- 0
  new.dominated <- 0
  current.pareto <- as.data.frame(current.pareto)
  new.pareto <- as.data.frame(new.pareto)
  for(i in 1:nrow(new.pareto)){
    for(j in 1:nrow(current.pareto)){
      if((new.pareto[i,"f1"] < current.pareto[j,"f1"]) && (new.pareto[i,"f2"] < current.pareto[j,"f2"])){
        old.dominated <- old.dominated + 1
      }else if((new.pareto[i,"f1"] > current.pareto[j,"f1"]) && (new.pareto[i,"f2"] > current.pareto[j,"f2"])){
        new.dominated <- new.dominated + 1
      }
    } 
  }
  
  convergence.index <- old.dominated/nrow(current.pareto) - new.dominated/nrow(new.pareto)
  print(paste0("Convergence index: ", convergence.index))
  return(convergence.index)
}

fitness_selection_crowding_distance <- function(R, P.size, K){
  R <- as.data.frame(R)
  last.rank <- R[P.size, "rnkIndex"]
  rank <- last.rank
  i <- P.size
  while(rank == last.rank){
    rank <- R[i, "rnkIndex"]
    i <- i - 1
  }
  preselected <- R[1:(i+1), ]
  remaining <- P.size - nrow(preselected)
  last.ranking.solutions <- R[R$rnkIndex == last.rank, ]
  last.ranking.solutions <- last.ranking.solutions[order(last.ranking.solutions$cd.density, decreasing=TRUE), ]
  last.ranking.solutions <- last.ranking.solutions[1:remaining, ]
  return(rbind(preselected, last.ranking.solutions))
}

fitness_selection_diversity_metric <- function(R, groups, P.size, K, metric){
  if(nrow(R) <= P.size){
    Log("Population size is lower than original population, returning...")
    return(R)
  }
  
  R <- as.data.frame(R)
  last.rank <- R[P.size, "rnkIndex"]
  if(R[P.size + 1, "rnkIndex"] > last.rank){ # No need to check if there's dominance
    return(R[1:P.size, ])
  }
  rank <- last.rank
  i <- P.size
  while(rank == last.rank){
    rank <- R[i, "rnkIndex"]
    i <- i - 1
  }
  preselected <- R[1:(i+1), ]
  remaining <- P.size - nrow(preselected)
  last.ranking.solutions.index <- which(R$rnkIndex == last.rank)
  last.ranking.solutions <- R[last.ranking.solutions.index, ]
  groups <- groups[last.ranking.solutions.index]
  diversity.matrix <- calculate_diversity_matrix(groups, metric)
  
  # Calculate mean distance from a solution to every other, and order it by decreasing
  mean.solution.distance <- apply(diversity.matrix, 1, function(x) mean(x))
  last.ranking.solutions <- last.ranking.solutions[order(mean.solution.distance, decreasing=TRUE), ]
  last.ranking.solutions <- last.ranking.solutions[1:remaining, ]
  return(rbind(preselected, last.ranking.solutions))
}

population_mating_and_mutation <- function(mating_pool, K, num.genes, params){
  Q <- matrix(0, nrow = params$popSize, ncol = K)
  mat.rate <- params$mating.rate
  mut.rate <- params$mutation.rate
  pop.size <- params$popSize
  genes <- 1:num.genes
  for(p in 1:pop.size){
    pair <- sample(1:nrow(mating_pool), 2, replace = FALSE)
    
    for(k in 1:K){
      genome.prob <- runif(1,0,1)
      mutation.prob <- runif(1,0,1)
      
      if(genome.prob < mat.rate){
        gene <- mating_pool[pair[1], k]
        if(is.element(gene, Q[p, ])){
          gene <- mating_pool[pair[2], k]
        }
      }else{
        gene <- mating_pool[pair[2], k]
        if(is.element(gene, Q[p, ])){
          gene <- mating_pool[pair[1], k]
        }
      }
      if(mutation.prob < mut.rate){
        gene <- sample(genes, 1)
      }
      while(is.element(gene, Q[p, ])){ 
        gene <- sample(genes, 1)
      }
      Q[p, k] <- gene
    }
  }
  rownames(Q) <- (nrow(mating_pool)+1):(nrow(mating_pool)*2)
  return(Q)
}

diverse_population_mating_and_mutation <- function(mating_pool, distances, groups, K, params, P.size=NULL){
  
  if(missing(P.size)){
    P.size <- params$popSize
  }
  D <- distances$exp.dist #******** Using expression distance temporally *******
  Q <- matrix(0, nrow = P.size, ncol = K)
  mat.rate <- params$mating.rate
  mut.rate <- params$mutation.rate
  
  genes <- 1:nrow(D)
  for(p in 1:P.size){
    pair <- sample(1:nrow(mating_pool), 2, replace = FALSE)
    chr1 <- mating_pool[pair[1], 1:K] # Chromosome 1
    chr2 <- mating_pool[pair[2], 1:K] # Chromosome 2
    for(k in 1:K){
      genome.prob <- runif(1,0,1)
      mutation.prob <- runif(1,0,1)
      gene.chr.1 <- chr1[1,k] # Cluster of chr1 k-th gene
      clust.dist <- D[gene.chr.1, unlist(chr2)] # Distances of chr1 to chr2's genes 
      target <- which.min(clust.dist)
      
      if(genome.prob < mat.rate){
        gene <- chr2[1, target]
        #If gene was already added, use the other chromosome. Low chance of occurrence
        if(is.element(gene, Q[p, ]) && !is.element(gene.chr.1, Q[p, ])){
          gene <- gene.chr.1
        }
      }else{
        gene <- gene.chr.1
        #If gene was already added, use the other chromosome. Low chance of occurrence
        if(is.element(gene, Q[p, ]) && !is.element(chr2[1, target], Q[p, ])){
          gene <- chr2[1, target]
        }
      }
      if(mutation.prob < mut.rate){
        #*** More diversity criteria can be added here ***
        group.chr1 <- groups[[p]][gene.chr.1]
        elems.group <- which(groups[[p]] == group.chr1)
        gene <- sample(elems.group, 1)
        #print(paste("Gene ", k, "mutated"))
      }
      # This only occurs in rare cases when random value was already in offspring chrmosome
      # Adding a while makes sure no error can occur
      while(is.element(gene, Q[p, ])){ 
        #*** More diversity criteria can be added here ***
        group.chr1 <- groups[[p]][gene.chr.1]
        elems.group <- which(groups[[p]] == group.chr1)
        gene <- sample(elems.group, 1)
      }
      Q[p, k] <- gene
    }
  }
  rownames(Q) <- (nrow(mating_pool)+1):(nrow(mating_pool)*2)
  return(Q)
}

##### Diversity metrics functions #####

calculate_diversity_matrix <- function(groups, metric){

  size <- length(groups)
  diversity.matrix <- matrix(0, ncol=size, nrow=size)
  for(i in 1:(size-1)){
    for(j in (i+1):size){
      diversity.matrix[i, j] = diversity_metric(unname(groups[[i]]), unname(groups[[j]]), metric)
    }
  }
  return(diversity.matrix)
}

diversity_metric <- function(group1, group2, metric){
  n.elem <- length(group1)
  pairs <- t(combn(1:n.elem, 2))
  yy <- 0 # True positive
  nn <- 0 # True negative
  yn <- 0 # False positive
  ny <- 0 # False negative
  for(pair in 1:nrow(pairs)){
    p1 <- pairs[pair, 1] # Gene 1 
    p2 <- pairs[pair, 2] # Gene 2
    g1.p1 <- group1[p1] # Group of gene 1 in group1
    g1.p2 <- group1[p2] # Group of gene 2 in group1
    g2.p1 <- group2[p1] # Group of gene 1 in group2
    g2.p2 <- group2[p2] # Group of gene 2 in group2
    if((g1.p1 == g1.p2) && (g2.p1 == g2.p2)){ # genes are in same cluster, in both solutions
      yy <- yy + 1
    }else if(g1.p1 == g1.p2){ # genes are in same cluster in group1, but not in group2
      yn <- yn + 1
    }else if(g2.p1 == g2.p2){ # genes are in same cluster in group2, but not in group1
      ny <- ny + 1
    }else{ # genes are not in the same cluster according to both solutions
      nn <- nn + 1
    }
  }
  if(metric == "jaccard"){
    return(jaccard_index(yy,ny,yn,nn))
  }else if(metric == "rand"){
    return(rand_index(yy,ny,yn,nn))
  #...
  }else{
    print(paste("ERROR: Metric", metric, "not supported!"))
    return(FALSE)
  } 
  #conf.matrix <- matrix(c(yy, ny, yn, nn), nrow=2, ncol=2)
  #return(conf.matrix)
  
}

jaccard_index <- function(yy, ny, yn, nn){
  #yy <- matrix[1,1]
  #ny <- matrix[2,1]
  #yn <- matrix[1,2]
  #nn <- matrix[2,2]
  return( yy / (yy + yn + ny) )
}

rand_index <- function(yy, ny, yn, nn){
  return ( (yy + nn) / sum(matrix) )
}

#### Evaluation and result plotting functions ####

# Plot partial pareto front
# **** Use normalization or absolute values??? To be reviewed *****
plot_pareto <- function(old.pareto, new.pareto, generation, output.path, agent=NULL){
  #print(paste("Pareto front for generation", generation))
  #print(new.pareto)
  output.path <- file.path(output.path, "pareto")
  new.pareto <- new.pareto[, c("f1", "f2")]
  old.pareto <- old.pareto[, c("f1", "f2")]
  #new.pareto <- as.data.frame(eaf::normalise(new.pareto[, c("f1", "f2")], c(0,1)))
  #old.pareto <- as.data.frame(eaf::normalise(old.pareto[, c("f1", "f2")], c(0,1)))
  
  # Define an scale function and normalize
  #scaler <- function(x){ (x-min(x))/(max(x)-min(x)) }
  #if(nrow(new.pareto) == 1 && nrow(old.pareto) == 1){ # If new pareto has 1 solution, normalize in reference
  #  temp <- rbind(old.pareto, new.pareto)
  #  temp <- as.data.frame(eaf::normalise(temp, c(0,1)))
  #  old.pareto <- temp[1, ]
  #  new.pareto <- temp[2, ]
  #}else if(nrow(old.pareto) == 1){
  #  temp <- rbind(old.pareto, new.pareto)
  #  temp <- as.data.frame(eaf::normalise(temp, c(0,1)))
  #  old.pareto <- temp[1, ]
  #}else if(nrow(new.pareto) == 1){
  #  temp <- rbind(new.pareto, old.pareto)
  #  temp <- as.data.frame(eaf::normalise(temp, c(0,1)))
  #  new.pareto <- temp[1, ]
  #}
  #}else{
  #new.pareto <- (new.pareto-min(new.pareto))/(max(new.pareto)-min(new.pareto))
  #old.pareto <- (old.pareto-min(old.pareto))/(max(old.pareto)-min(old.pareto))
  #}
  
  label <- c(rep("current generation", nrow(new.pareto)), rep("last generation", nrow(old.pareto)))
  pareto <- rbind(new.pareto, old.pareto)
  pareto <- cbind(pareto, label)
  colnames(pareto) <- c("f1", "f2", "Pareto")
  alpha <- ifelse(label == "last generation", 0.5, 1)
  #pareto$color <- rep("grey", nrow(pareto))
  suppressWarnings(ggplot(pareto, aes(x=f1, y=f2, colour=Pareto)) +
    labs(title=paste0("Generation ", generation), x="Expression Index", y="Biological Index") +
    geom_point(aes(color = Pareto), size=2) +
    geom_line(aes(group = Pareto, alpha=alpha)) +
    scale_alpha_continuous(guide=FALSE) +
    xlim(0, max(pareto$f1)) +
    ylim(0, max(pareto$f2))  
  )
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE)
  if(missing(agent)){
    filename <- file.path(output.path, paste0("g_", generation, "agent_", agent,".png")) 
  }else{
    filename <- file.path(output.path, paste0("g_", generation, ".png"))
  }
  suppressWarnings(ggsave(filename, height=5, width=7))
}

evaluate_solutions <- function(population, clustering, K, distances, params, output.path){
  pareto <- population[population$rnkIndex == 1, ]
  obj.index <- (K+1):(K+params$objDim)
  N <- nrow(pareto)
  #dist.compounded <- params$alpha*distances$bio.dist + (1-params$alpha)*distances$exp.dist
  output <- file.path(output.path, "results")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  
  # Create silhouette results table
  labels <- rep(NA, K)
  for(i in 1:K){ labels[i] <- paste("Sil. cluster", i) }
  labels[K+1] <- "Average sil"
  silhouette.res <- data.frame(matrix(ncol=K+1, nrow=N))
  colnames(silhouette.res) <- labels
  
  for(i in 1:N){
    # Calculate Silhouette metric
    sil <- cluster::silhouette(clustering[[i]], distances$exp.dist)#dist.compounded)
    pam.res <- pam(distances$exp.dist, metric = "euclidean", k = K, medoids = pareto[i, 1:K], do.swap = FALSE)
    plot.pam <- factoextra::fviz_cluster(pam.res, geom = "point")
    plot.sil <- factoextra::fviz_silhouette(sil)
    
    out.file <- file.path(output, paste0("p", i, "_pam.png"))
    png(out.file)
    print(plot.pam)
    dev.off()
    out.file.2 <- file.path(output, paste0("p", i, "_sil.png"))
    png(out.file.2)
    print(plot.sil)
    dev.off()
    #plot(sil, main=paste("Silhouette plot: Solution", i), border=1, col=1:N, nmax.lab=10)
    res <- summary(sil)
    
    silhouette.res[i, 1:K] <- unname(res$clus.avg.widths)
    silhouette.res[i, K+1] <- res$avg.width
    #results[i, 1] <- mean(sil[,"sil_width"])
  }
  # Calculate Hypervolume metric
  # Calculate nadir point for hypervolume
  nadir.point <- rep(NA, params$objDim)
  for(obj in 1:params$objDim){
    if(params$obj_maximize[obj]){
      nadir.point[obj] <- min(population[, K+obj])*1.2 # worse sol. + delta distance  
    }else{
      nadir.point[obj] <- max(population[, K+obj])*1.2 # worse sol. + delta distance  
    }
  }
  hv.res <- eaf::hypervolume(population[ , obj.index], nadir.point, params$obj_maximize)
  rownames(silhouette.res) <- 1:N
  
  results <- list("silhouette"=silhouette.res, "hypervolume"=hv.res)
  write.csv(silhouette.res, file.path(output.path, "silhouette.csv"), row.names=TRUE)
  write(hv.res, file.path(output.path, "hypervolume.txt"))
  
  return(results)
}

#NOT USED
hipervolume_projection <- function(pareto, ref.vector, output.path){
  projected.pareto <- pareto
  size <- nrow(pareto)
  obj.dim <- ncol(pareto)
  for(i in 1:size){
    for(obj in 1:obj.dim){
      projected.pareto[i, obj] <- pareto[i, obj] + (1 - sum(pareto[i, ]))/obj.dim
    }
  }
  labels <- c( rep("Original", size), rep("Projected", size) )
  points <- rbind(pareto, projected.pareto)
  points <- cbind( points, labels)
  colnames(points) <- c("f1", "f2", "Hypervolume")
  print(points)
  # Scatterplot of HV
  ggplot(points, aes(x=f1, y=f2, color=Hypervolume, shape=Hypervolume)) +
    geom_point() + 
    labs(main="Hypervolume projection", xlab="Genetic expression", ylab="Biological function")
  # Save the figure
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE)
  filename <- file.path(output.path, "pareto.png")
  ggsave(filename, height=7, width=7)
  
  return(projected.pareto)
}


#### Memetic algorithm functions #####

#https://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
#From hendalst
Log <- function(text, agent, ...) {
  if(missing(agent)){
    msg <- sprintf(paste0(as.character(Sys.time()), ": " , text, "\n"), ...)
  }else{
    msg <- sprintf(paste0(as.character(Sys.time()), ": [Agent ", agent, "] ", text, "\n"), ...)
  }
  cat(msg)
}

diverse_memetic_nsga2 <- function(distances, K, diversity.metric, diversity.level, params, output.path){
  
  agents <- params$agents
  P.size <- params$popSize
  phases <- params$phases
  phase <- 1
  Agents <- list()
  
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE)
  output.log.file <- file.path(output.path, "log.txt")
  
  sink(output.log.file, append=TRUE)
  writeLines(c(""), output.log.file)
  
  Log(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"..."))
  # Initialize population
  if(params$is_random_population){
    P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed, params$is_random_population) # Random population
  }else if(diversity.level >= 1){
    P <- generate_diverse_initial_pop(distances, params, K, p.size=P.size, diverse_population = TRUE)
  }else{
    P <- generate_diverse_initial_pop(distances, params, K, p.size=P.size, diverse_population = FALSE)
  }
  P.size <- nrow(P) # Population size
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, P.clustering.groups, distances, params)
  P.clustering.groups <- P.clustering.groups[order(as.numeric(rownames(P)))]
  rownames(P) <- 1:nrow(P)
  
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  
  generations.no.changes <- 0
  has.converged <- FALSE
  
  pop.per.agent <- floor(P.size/agents)
  gen.per.phase <- params$generations/phases
  # Assign each agent with a random set of solutions from the solution pool
  solutions <- rep(1:agents, each=pop.per.agent)
  solutions <- sample(solutions)
  solutions <- split(1:P.size, solutions)
  
  for(i in 1:agents){
    Agents <- c(Agents, list( P[ unlist( solutions[i] ) , ] ) )  
  }
  
  cl <<- makeCluster(agents)
  registerDoParallel(cl)
  
  source("dmmoea_functions.R")
  exported_func <- c("dnsga2")
  #exported_func <- c("dnsga2", "cluster_data", "evaluate_population", "evaluate_xie_beni", "dominance_ranking_sorting", "Log")
  #clusterExport(cl, c("agents", "phases", exported_func))
  clusterExport(cl, exported_func)
  
  while (phase <= phases) {
    
    Agents <- foreach(i=1:params$agents, .combine = c, .export=c("pop.per.agent", "gen.per.phase", "phase"), .inorder = FALSE) %dopar% {
      
      #pareto_agent <- dnsga2(Agents[[i]], i, pop_per_agent, ev_per_agent, param, phase)
      pareto_agent <- dnsga2(initial_population=Agents[[i]], agent=i, P.size=pop.per.agent, generations=gen.per.phase, distances, K, diversity.metric, diversity.level, params, output.path)
      return( list( pareto_agent ) )
    }
    print(Agents)
    
    Log(paste("Phase", phase, "Agents finished!"))
    
    # Synchronize agents
    print(agents)
    for(i in 1:agents){
      for(j in 1:agents){
        if(i == j) {
          next
        }
        #solutions_taken <- population / param$num_agents
        #threshold <- population - solutions_taken
        #temp <- tournamentSelection(rbind(Agents[[i]], Agents[[j]]), pop_agent_i + pop_agent_j, param$tourSize)
        
        #Agents[[i]] <- fitness(Agents[[i]], Agents[[j]], "min", "min", pop_per_agent)
        print(paste("Testing agent", i, ",",j))
      }
    }
    Log("Phase %d synchronization ended!\n", phase)
    Log("Phase %d solutions: \n", phase)
    print(Agents)
    
    #We need to discard the other information before using it for the new phase population
    if(phase != phases){
      for(i in 1:agents){
        Agents[[i]] <- Agents[[i]][ , 1:K]
      }
    }
    
    print("Next phase")
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
  for(i in 1:agents){
    pareto <- rbind(pareto, Agents[[i]])
  }
  pareto<- pareto[!duplicated(pareto),] # Delete duplicates
  
  Log("Deep Memetic NSGA-2 Algorithm finished!!")
  
  
  # Pruning non pareto-front solutions
  #objetives <- as.matrix(pareto[, c("f1", "f2")]) 
  #ranks <- NonDomSorting(objetives) 
  #pareto <- pareto[ unlist(ranks[[1]]) , ] #Select the rank 1 indices
  # No need to limit number of population now
  #if(nrow(pareto) > pop_size){
  #  pareto <- pareto[1:pop_size, ]
  #}
  
  #Obtain results
  #guarda_paretos=carga_guarda_pareto(alfa, num_k,  distancias$matriz_exp,  distancias$distancia_exp, distancias$distancia_bio, pop_per_agent, pareto, FALSE)
  
  #File result
  #exp.name <- paste0(output_result_file, "\\B_memetic")
  
  #save results
  #write.table(guarda_paretos$pareto_clustering, file = paste0(exp.name,"_clustering.csv"),sep = ",",col.names = TRUE, quote = FALSE)
  #write.table(guarda_paretos$pareto_objectives, file = paste0(exp.name,".csv"),sep = ",",col.names = TRUE, row.names = TRUE,quote = FALSE)
  #write.table(guarda_paretos$pareto_silhouette, file = paste0(exp.name,"_silhouettes.csv"),sep = ",",col.names = FALSE, quote = FALSE)
  
  #show results
  #Log("Showing results: \n")
  #print(pareto)
  #print(guarda_paretos)
  
  return(pareto)
}
