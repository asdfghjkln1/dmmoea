generate_initial_pop <- function(P.size, K, num.genes, seed){
  # Generate random population
  #set.seed(seed) # set seed
  population <- t(sapply(1:P.size, function(x)  sample(1:num.genes, K ,replace=F)  ))
  #rm(.Random.seed, envir=globalenv()) ?
  return(as.matrix(population))
  
}

generate_diverse_initial_pop <- function(distances, params, p.size=NULL, diverse_population=FALSE){
  K <- params$K
  n.genes <- distances$n.genes
  if(missing(p.size)){
    p.size <- params$popSize
  }
  # Normalise distance matrix
  #D <- distances$comp.dist # params$alpha*distances$bio.dist + (1-params$alpha)*distances$exp.dist # Doesn't work? Weird numbers
  #D <- distances$exp.dist
  D <- distances$comp.dist
  #D <- as.data.frame(eaf::normalise(D, c(0,1)))
  D.size <- nrow(D)
  population <- matrix(nrow=p.size, ncol=K)
  
  # Create a chromosome (single element of population) for every radius of density
  # Set density tolerance T (number of elem. within density radius)
  auto.adjust <- params$auto_adjust_initial_params
  
  if(!is.na(auto.adjust)){
    if(auto.adjust){
      T <- n.genes/K
      radius <- T/n.genes
      min.r <- max(0.1, radius*0.5)
      max.r <- min(0.5, radius*1.5)
      r.series <- seq(min.r, max.r, by=(max.r - min.r)/p.size) 
    }else{
      T <- n.genes * params$density_tol # *** Might need to be adjusted ***
      r.series <- seq(params$min_density_radius,params$max_density_radius, by=(params$max_density_radius-params$min_density_radius)/p.size)
    }
  }else{
    T <- n.genes * params$density_tol # *** Might need to be adjusted ***
    r.series <- seq(params$min_density_radius,params$max_density_radius, by=(params$max_density_radius-params$min_density_radius)/p.size)
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
          population[p, k] <- g
          gene.selected <- TRUE
          candidate.search.index <- min(i + 1, length(candidates)) # Avoids repeating search
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
  #bio.dist <- distances$bio.dist
  #exp.dist <- distances$exp.dist
  dist.compounded <- distances$comp.dist
  gene.len <- distances$n.genes
  clustering.results <- list(1:pop.size)
  grouping<-rep(NA,gene.len)
  p <- 1
  while(p <= pop.size){
    for(gene in 1:gene.len){
      #medoid.distances <- alpha*exp.dist[gene, population[p, 1:K]] + (1-alpha)*bio.dist[gene, population[p, 1:K]]
      medoid.distances <- dist.compounded[gene, population[p, 1:K]] # Doesn't work this way, not sure why
      gene.clust <- unname(which.min(medoid.distances))
      grouping[gene] <- gene.clust
      
    }
    clustering.results[[p]] <- grouping
    names(clustering.results[[p]])<-rownames(dist.compounded)
    col <- 1
    # For every cluster, check for singletons. 
    # Replace medoid until all clusters in the solution has more than 2 elements
    while(col <= K){
      if(length(which(unlist(grouping) == col)) < 3){
        replacement <- setdiff(1:gene.len, intersect(1:gene.len, population[p, 1:K]))
        population[p, col] <- sample(replacement, 1)
        new.medoids <- population[p, 1:K]
        #** Update groups for new medoids
        for(gene in 1:gene.len){
          #medoid.distances <- alpha*exp.dist[gene, new.medoids ] + (1-alpha)*bio.dist[gene, new.medoids]
          medoid.distances <- dist.compounded[gene, population[p, 1:K]]
          gene.clust <- unname(which.min(medoid.distances))
          grouping[gene] <- gene.clust
        }
        clustering.results[[p]] <- grouping
        names(clustering.results[[p]])<-rownames(dist.compounded)
        col <- 1 # Check for singletons in previous medoids, until every medoid has no singletons
      }else{
        col <- col + 1
      }
    }
    p <- p + 1
  }
  return(list("population"=population, "clustering.results"=clustering.results))
}

evaluate_silhouette <- function(gene.dist, groups, N){
  # Create silhouette results table
  f <- rep(NA, N)
  for(i in 1:N){
    sil <- cluster::silhouette(groups[[i]], gene.dist)
    res <- summary(sil)
    f[i] <- res$avg.width
    
  }
  return(f)
}

evaluate_population <- function(population, distances, groups, params){
  K <- params$K
  objectives <- rep(params$objectives, 2)
  obj.dim <- length(objectives)
  obj1 <- objectives[1]
  obj2 <- objectives[2]
  population <- as.matrix(population)
  #p <- length(cluster_results) # not used atm
  if(obj1 == "XieBeni"){
    f1 <- evaluate_xie_beni(population, distances$exp.dist)
  }else if(obj1 == "BallHall"){
    f1 <- evaluate_ball_hall(population, distances$exp.dist, groups)
  }else if(obj1 == "Dunn"){
    f1 <- evaluate_dunn_index(distances$exp.dist, groups, nrow(population))
  }else if(obj1 == "Silhouette"){
    f1 <- evaluate_silhouette(distances$exp.dist, groups, nrow(population))
  }else{
    #TO DO implement
    warning("Objective function not available!")
    return(NULL)
  }
  
  if(obj2 == "XieBeni"){
    f2 <- evaluate_xie_beni(population, distances$bio.dist)
  }else if(obj2 == "BallHall"){
    f2 <- evaluate_ball_hall(population, distances$bio.dist, groups)
  }else if(obj2 == "Dunn"){
    f2 <- evaluate_dunn_index(distances$bio.dist, groups, nrow(population))
  }else if(obj2 == "Silhouette"){
    f2 <- evaluate_silhouette(distances$bio.dist, groups, nrow(population))
  }else{
    #TO DO implement
    warning("Objective function not available!")
    return(NULL)
  }
  
  # Fixing possible outliers from evaluation metrics, replace bugged values (NA) with worst values
  custom.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), Inf)
  custom.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), -Inf)
  worst.f1 <- ifelse(params$obj_maximize[1], custom.min(f1), custom.max(f1))
  worst.f2 <- ifelse(params$obj_maximize[2], custom.min(f2), custom.max(f2))
  for(i in 1:nrow(population)){
    f1[i] <- ifelse(is.na(f1[i]), worst.f1, f1[i])
    f2[i] <- ifelse(is.na(f2[i]), worst.f2, f2[i])
  }
  #apply(f1, 1, function(x) ifelse(is.na(x), worst.f1, x))
  #apply(f2, 1, function(x) ifelse(is.na(x), worst.f2, x))
  
  obj.values <- cbind(f1,f2)
  population <- dominance_ranking_sorting(population, obj.values)
  return(population)
}

# Xie-Beni clustering coefficient evaluation
# Evaluates a single objective
evaluate_xie_beni <- function(population, distances){
  
  obj.function <-rep(0, nrow(population))
  n <- nrow(distances) # number of elements
  K <- ncol(population) # number of medoids
  for (p in 1:nrow(population)) {
    
    # Initialize matrix
    D <- matrix(0, K, n) #Matriz de distancia centros x elementos
    
    # Calculate square distances of each element to its medoid
    for (k in 1:K) {
      cluster.medoid <- population[p, k]
      for (i in 1:n) {
        D[k, i] <- distances[cluster.medoid, i]^2
      }
    }
    
    XB.numerator = sum(D)
  
    combs <- t(combn(population[p, ], 2)) # Generate pairs of medoids
    s <- vector()
    for (pair in 1:nrow(combs)) {
      s[pair] = distances[combs[pair, 1], combs[pair, 2]]^2
    }
    s <- s[s > 0]
    XB.denominator <- min(s)

    if(is.finite(XB.denominator)){
      obj.function[[p]]=XB.numerator/(n*XB.denominator) 
    }else{
      obj.function[[p]]=NA
    }
  }
  return(unlist(obj.function))
}

evaluate_ball_hall <- function(population, distances, groups){
  
  obj.function <-rep(0, nrow(population))
  
  for (p in 1:nrow(population)) {
    
    # Initialize matrix
    n <- nrow(distances) # number of elements
    K <- ncol(population) # number of medoids
    
    # Calculate square distances of each element to its medoid
    sum.dist <- rep(0, K)
    for (k in 1:K) {
      cluster.medoid <- population[p, k]
      group.medoid <- groups[[p]][cluster.medoid]
      in.cluster <- which(groups[[p]] == group.medoid)
      
      dist.to.medoid <- distances[cluster.medoid, in.cluster]
      sum.dist[k] <- sum(dist.to.medoid^2)
    }
    
    obj.function[p] <- sum(sum.dist)/K
  }
  return(obj.function)
}

# implement PMB -> evaluate_PMB_index ?
# implement Inverted generational distance -> evaluate_IGD_index ?

evaluate_dunn_index <- function(distances, groups, N){
  
  obj.function <-rep(0, N)
  
  for (p in 1:N) {
    
    obj.function[p] <- clValid::dunn(distance=distances, clusters=groups[[p]])
  }
  
  return(obj.function)
}

get_new_dirname <- function(path){
  dir.name <- basename(path)
  dir.base <- dirname(path)
  #exp.id <- strsplit(file.name, "_")[[1]][2]
  #files <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
  if(dir.exists(path)){
    new.name <- paste0(dir.name, "(2)")
    counter <- 3
    while(dir.exists(file.path(dir.base, new.name))){
      new.name <- paste0(dir.name, "(", counter, ")")
      counter <- counter + 1
    }
    return(file.path(dir.base, new.name))
    #match <- grep(paste0("_", exp.id), files, value=TRUE)
  }else{
    return(file.path(path))
  }
}

#### NSGA-2 ####
nsga2 <- function(distances, params, output.path, debug=FALSE, plot=FALSE){
  evaluation.count <- 0
  evaluations <- params$evaluations
  K <- params$K
  
  #output.path <- get_new_dirname(output.path)
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE) 
  
  if(debug){
    output.log.file <- file.path(output.path, "log.txt")
    sink(output.log.file, append=FALSE)
    writeLines(c(""), output.log.file)
    Log(paste("Initiating NSGA-2...")) 
  }
  
  # Initialize Random population
  P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed) 
  P.size <- nrow(P) # Population size
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, distances, P.clustering.groups, params)
  evaluation.count <- evaluation.count + nrow(P)
  P.clustering.groups <- P.clustering.groups[as.numeric(rownames(P))]
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  #Set condition for early convergence if pareto front does not change
  generations.no.changes <- 0
  estimated.generations <- params$evaluations/P.size
  generations.no.change.limit <- round(estimated.generations*3/4) 
  
  has.converged <- FALSE
  ## Main generation loop
  while(!has.converged){
    
    if(debug){
      print(paste("Entering selection, crossover, mutation for generation ", g)) 
      print(paste0("Evaluations used (", evaluation.count, "/", evaluations, ")"))
    }
    
    if(nrow(P) < P.size){ # Fill solutions if population size is not enough. Does not happen commonly
      to.fill <- P.size-nrow(P)
      print(paste("Not enough solutions in population. Filling", to.fill, "solutions"))
      P.fill.solutions <- fill_population(P, K, num.genes, fill=to.fill)
      P.fill.data <- cluster_data(distances, P.fill.solutions, params$alpha)
      P.fill <- as.data.frame(P.fill.data$population)
      rownames(P.fill) <- as.character(round(runif(to.fill, 0, 1)*1e4))
      P.fill.clustering <- P.fill.data$clustering.results
      P.fill.rows <- row.names(P.fill)
      P.fill <- evaluate_population(P.fill, distances, P.fill.clustering, params)
      evaluation.count <- evaluation.count + nrow(P.fill)
      P.fill.clustering <- P.fill.clustering[match((rownames(P.fill)), P.fill.rows)]
      P <- rbind(P, P.fill)
      P.clustering.groups <- c(P.clustering.groups, P.fill.clustering)
      obj.values <- P[, (K+1):(K+params$objDim)] # Select objective values
      P.rows <- rownames(P)
      P <- dominance_ranking_sorting(P[, 1:K], obj.values) # Recalculate ranking
      P.clustering.groups <- P.clustering.groups[match((rownames(P)), P.rows)]  # Update clustering
      rownames(P)<-1:nrow(P)
    }
    
    ###### Selection, Crossover, Mutation  ######
    
    ## Selection
    mating_pool <- nsga2R::tournamentSelection(P, P.size, params$tourSize)
    
    ##Crossover and Mutation
    Q <- population_mating_and_mutation(mating_pool, num.genes, params)
    
    ## Evaluate solutions in population
    Q.data <- cluster_data(distances, Q, params$alpha)
    Q <- Q.data$population
    Q.clustering.groups <- Q.data$clustering.results
    Q.rows <- rownames(Q) 
    Q <- evaluate_population(Q, distances, Q.clustering.groups, params)
    evaluation.count <- evaluation.count + nrow(Q)
    Q.clustering.groups <- Q.clustering.groups[match((rownames(Q)), Q.rows) ]
    rownames(Q) <- 1:nrow(Q)
    
    ## Create population R as the combination of P and Q
    R <- rbind(P,Q)
    R.clustering.groups <- c(P.clustering.groups, Q.clustering.groups)
    R.rows <- rownames(R)
    R <- remove_duplicated(R, K)
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering
    
    obj.values <- R[, (K+1):(K+params$objDim)] # Select objective values
    R.rows <- rownames(R)
    R <- dominance_ranking_sorting(R[, 1:K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering
    rownames(R)<-1:nrow(R)
    
    ## Population fitness selection
    P_next_generation <- fitness_selection_crowding_distance(R, P.size, K) 
    P.clustering.groups <- R.clustering.groups[as.numeric(rownames(P_next_generation))] # Update clustering
    rownames(P_next_generation) <- 1:nrow(P_next_generation)
    new.pareto.front <- P_next_generation[P_next_generation$rnkIndex == 1, ]
    pareto.clustering <- P.clustering.groups[as.numeric(rownames(new.pareto.front))]
    
    ## Output pareto front plot
    if(plot){
      plot_pareto(current.pareto.front, new.pareto.front, g, output.path) # Output pareto front 
    }
    
    ## Measure convergence of pareto front
    convergence.index <- convergence_coefficient(current.pareto.front, new.pareto.front, g, params$obj_maximize)
    
    ## Check how different is the new pareto front, count generations with no changes
    if(convergence.index <= params$convergence_tol){
      generations.no.changes <- generations.no.changes + 1 
    }else{
      generations.no.changes <- 0 # Reset counter if changes detected
    }
    ## Check for stop criteria
    if(generations.no.changes >= generations.no.change.limit || evaluation.count > evaluations){
      has.converged <- TRUE
    }
    ## Continue to next generation
    g <- g + 1
    P <- P_next_generation
    current.pareto.front <- new.pareto.front
    if(debug){
      print(P) 
    }
  }
  if(debug){
    sink(type="output") 
  }
  
  return(list("population"=P_next_generation, "clustering"=P.clustering.groups))
}


#### Diverse NSGA-2 ####
dnsga2 <- function(distances, params, output.path, debug=FALSE, plot=FALSE){
  evaluation.count <- 0
  evaluations <- params$evaluations
  diversity.metric <- params$diversity_metric
  diversity.level <- params$diversity_level
  K <- params$K
  
  #output.path <- get_new_dirname(output.path)
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE) 
  if(debug){ 
    output.log.file <- file.path(output.path, "log.txt")
    sink(output.log.file, append=FALSE)
    writeLines(c(""), output.log.file)
    Log(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"...")) 
  }
  
  # Initialize population
  if(params$is_random_population){
    # Random population
    P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed) 
  }else if(diversity.level >= 1){
    P <- generate_diverse_initial_pop(distances, params, diverse_population = TRUE)
  }else{
    P <- generate_diverse_initial_pop(distances, params, diverse_population = FALSE)
  }
  
  P.size <- params$popSize # Population size
  #Set condition for early convergence if pareto front does not change
  estimated.generations <- params$evaluations/P.size
  generations.no.change.limit <- round(estimated.generations*3/4) 
  
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, distances, P.clustering.groups, params)
  evaluation.count <- evaluation.count + nrow(P)
  P.clustering.groups <- P.clustering.groups[as.numeric(rownames(P))]
  rownames(P) <- 1:P.size
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  generations.no.changes <- 0
  has.converged <- FALSE
  
  ## Main generation loop
  while(!has.converged){
    if(debug){
      Log(paste("Entering selection, crossover, mutation for generation ", g))
      Log(paste0("Evaluations used (", evaluation.count, "/", evaluations, ")"))
    }
    if(nrow(P) < P.size){ # Fill solutions if population size is not enough. Does not happen commonly
      to.fill <- P.size-nrow(P)
      Log(paste("Not enough solutions in population. Filling", to.fill, "solutions"))
      P.fill.solutions <- fill_population(P, K, num.genes, fill=to.fill)
      P.fill.data <- cluster_data(distances, P.fill.solutions, params$alpha)
      P.fill <- as.data.frame(P.fill.data$population)
      rownames(P.fill) <- as.character(round(runif(to.fill, 0, 1)*1e4))
      P.fill.clustering <- P.fill.data$clustering.results
      P.fill.rows <- row.names(P.fill)
      P.fill <- evaluate_population(P.fill, distances, P.fill.clustering, params)
      evaluation.count <- evaluation.count + nrow(P.fill)
      P.fill.clustering <- P.fill.clustering[match((rownames(P.fill)), P.fill.rows)]
      P <- rbind(P, P.fill)
      P.clustering.groups <- c(P.clustering.groups, P.fill.clustering)
      obj.values <- P[, (K+1):(K+params$objDim)] # Select objective values
      P.rows <- rownames(P)
      P <- dominance_ranking_sorting(P[, 1:K], obj.values) # Recalculate ranking
      P.clustering.groups <- P.clustering.groups[match((rownames(P)), P.rows)]  # Update clustering
      rownames(P)<-1:nrow(P)
    }
    
    ###### Selection, Crossover, Mutation  ######
    
    ## Selection
    mating_pool <- nsga2R::tournamentSelection(P, P.size, params$tourSize)
    
    ##Crossover and Mutation
    if(diversity.level >= 3){
      Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, params, P.size=P.size) 
    }else{
      Q <- population_mating_and_mutation(mating_pool, num.genes, params)
    }
    
    ## Evaluate solutions in population
    Q.data <- cluster_data(distances, Q, params$alpha)
    Q <- Q.data$population
    Q.clustering.groups <- Q.data$clustering.results
    Q.rows <- rownames(Q) 
    Q <- evaluate_population(Q, distances, Q.clustering.groups, params)
    evaluation.count <- evaluation.count + nrow(Q)
    Q.clustering.groups <- Q.clustering.groups[match((rownames(Q)), Q.rows) ]
    rownames(Q) <- 1:nrow(Q)
    
    ## Create population R as the combination of P and Q
    R <- rbind(P,Q)
    R.clustering.groups <- c(P.clustering.groups, Q.clustering.groups)
    R.rows <- rownames(R)
    R <- remove_duplicated(R, K)
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering
    
    obj.values <- R[, (K+1):(K+params$objDim)] # Select objective values
    R.rows <- rownames(R)
    R <- dominance_ranking_sorting(R[, 1:K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering
    rownames(R)<-1:nrow(R)
    
    ## Population fitness selection
    if(diversity.level >= 2){
      P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, P.size, K, diversity.metric)
    }else{
      P_next_generation <- fitness_selection_crowding_distance(R, P.size, K) 
    }
    P.clustering.groups <- R.clustering.groups[as.numeric(rownames(P_next_generation))] # Update clustering
    rownames(P_next_generation) <- 1:nrow(P_next_generation)
    new.pareto.front <- P_next_generation[P_next_generation$rnkIndex == 1, ]
    pareto.clustering <- P.clustering.groups[as.numeric(rownames(new.pareto.front)) ]
    
    if(plot){
      plot_pareto(current.pareto.front, new.pareto.front, g, output.path) # Output pareto front 
    }
    
    ## Measure convergence of pareto front
    convergence.index <- convergence_coefficient(current.pareto.front, new.pareto.front, g, params$obj_maximize)
    
    ## Check how different is the new pareto front, count generations with no changes
    if(convergence.index <= params$convergence_tol){
      generations.no.changes <- generations.no.changes + 1 
    }else{
      generations.no.changes <- 0 # Reset counter if changes detected
    }
    ## Check for stop criteria
    if(generations.no.changes > generations.no.change.limit ||  evaluation.count > evaluations){
      has.converged <- TRUE
    }
    
    ## Continue to next generation
    g <- g + 1
    P <- P_next_generation
    current.pareto.front <- new.pareto.front
    if(debug){
      print(P) 
    }
  }
  if(debug){
    sink(type="output") 
  }
  return(list("population"=P_next_generation, "clustering"=P.clustering.groups))
}

#### Diverse NSGA-2 ####
dnsga2_agent <- function(distances, params, output.path, P.size, agent, phase, evaluations, initial_population, debug=FALSE, plot=FALSE){
  evaluation.count <- 0
  diversity.metric <- params$diversity_metric
  diversity.level <- params$diversity_level
  K <- params$K
  # Load function workspace
  source("dmmoea_libraries.R")
  source("dmmoea_functions.R")
  
  if(debug){
    log.file <- file.path(output.path, paste0("log_",agent,".txt"))
    if(file.exists(log.file)){
      sink(log.file, append=TRUE) # Register in log output file
    }else{
      sink(log.file, append=FALSE)
    }
    Log(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"..."), agent=agent) 
  }
  
  # Initialize population
  if(!is.null(initial_population)){
    P <- as.data.frame(initial_population$population)
    P.clustering.groups <- initial_population$clustering
    current.pareto.front <- P[P$rnkIndex == max(P$rnkIndex), ] # Current pareto front
    #*** TODO: Check initial population validity as input **** #
  }else{
    warning("Error: No initial population provided.")
    return(NULL)
  }
  
  #Set condition for early convergence if pareto front does not change
  estimated.generations <- params$evaluations/P.size
  generations.no.change.limit <- round(estimated.generations*3/4) 
  num.genes <- distances$n.genes
  g <- 1 # Current generation
  generations.no.changes <- 0
  has.converged <- FALSE
  
  ## Main generation loop
  while(!has.converged){
    if(debug){
      Log(paste("Entering selection, crossover, mutation for generation ", g), agent=agent)
      Log(paste0("Evaluations used (", evaluation.count, "/", evaluations, ")"))
    }
    if(nrow(P) < P.size){ # Fill solutions if population size is not enough. Does not happen commonly
      to.fill <- P.size-nrow(P)
      Log(paste("Not enough solutions in population. Filling", to.fill, "solutions", agent=agent))
      P.fill.solutions <- fill_population(P, K, num.genes, fill=to.fill)
      P.fill.data <- cluster_data(distances, P.fill.solutions, params$alpha)
      P.fill <- as.data.frame(P.fill.data$population)
      rownames(P.fill) <- as.character(round(runif(to.fill, 0, 1)*1e4))
      P.fill.clustering <- P.fill.data$clustering.results
      P.fill.rows <- row.names(P.fill)
      P.fill <- evaluate_population(P.fill, distances, P.fill.clustering, params)
      evaluation.count <- evaluation.count + nrow(P.fill)
      P.fill.clustering <- P.fill.clustering[match((rownames(P.fill)), P.fill.rows)]
      P <- rbind(P, P.fill)
      P.clustering.groups <- c(P.clustering.groups, P.fill.clustering)
      obj.values <- P[, (K+1):(K+params$objDim)] # Select objective values
      P.rows <- rownames(P)
      P <- dominance_ranking_sorting(P[, 1:K], obj.values) # Recalculate ranking
      P.clustering.groups <- P.clustering.groups[match((rownames(P)), P.rows)]  # Update clustering
      rownames(P)<-1:nrow(P)
    }
    
    ###### Selection, Crossover, Mutation  ######
    
    ## Selection
    mating_pool <- nsga2R::tournamentSelection(P, P.size, params$tourSize)
    ##Crossover and Mutation
    if(diversity.level >= 3){
      Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, params, P.size=P.size) 
    }else{
      Q <- population_mating_and_mutation(mating_pool, num.genes, params, P.size=P.size)
    }
    ## Evaluate solutions in population
    Q.data <- cluster_data(distances, Q, params$alpha)
    Q <- Q.data$population
    Q.clustering.groups <- Q.data$clustering.results
    Q.rows <- rownames(Q) 
    Q <- evaluate_population(Q, distances, Q.clustering.groups, params)
    evaluation.count <- evaluation.count + nrow(Q)
    Q.clustering.groups <- Q.clustering.groups[match((rownames(Q)), Q.rows) ]
    rownames(Q) <- 1:nrow(Q)
    
    ## Create population R as the combination of P and Q
    R <- rbind(P,Q)
    R.clustering.groups <- c(P.clustering.groups, Q.clustering.groups)
    R.rows <- rownames(R)
    R <- remove_duplicated(R, K)
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering

    obj.values <- R[, (K+1):(K+params$objDim)] # Select objective values
    R.rows <- rownames(R)
    R <- dominance_ranking_sorting(R[, 1:K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[match((rownames(R)), R.rows)] # Update clustering
    rownames(R)<-1:nrow(R)
  
    ## Population fitness selection
    if(diversity.level >= 2){
      P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, P.size, K, diversity.metric)
    }else{
      P_next_generation <- fitness_selection_crowding_distance(R, P.size, K) 
    }
    P.clustering.groups <- R.clustering.groups[as.numeric(rownames(P_next_generation))] # Update clustering
    rownames(P_next_generation) <- 1:nrow(P_next_generation)
    new.pareto.front <- P_next_generation[P_next_generation$rnkIndex == 1, ]
    pareto.clustering <- P.clustering.groups[as.numeric(rownames(new.pareto.front)) ]
    if(plot){
      plot_pareto(current.pareto.front, new.pareto.front, g, output.path, agent=agent, phase=phase) # Output pareto front
    }
    
    ## Measure convergence of pareto front
    convergence.index <- convergence_coefficient(current.pareto.front, new.pareto.front, g, params$obj_maximize)
    
    ## Check how different is the new pareto front, count generations with no changes
    if(convergence.index <= params$convergence_tol){
      generations.no.changes <- generations.no.changes + 1 
    }else{
      generations.no.changes <- 0 # Reset counter if changes detected
    }
    ## Check for stop criteria
    if(generations.no.changes > generations.no.change.limit || evaluation.count > evaluations){
      has.converged <- TRUE
    }
    
    ## Continue to next generation
    g <- g + 1
    P <- P_next_generation
    if(debug){
      print(P) 
    }
    current.pareto.front <- new.pareto.front
  }
  if(debug){
    sink(type="output") 
  }
  
  return(list("population"=P_next_generation, "clustering"=P.clustering.groups))
}

convergence_coefficient <- function(current.pareto, new.pareto, generation, maximize){
  if(generation == 0){
    return(FALSE) # Cant converge at first generation
  }
  old.dominated <- 0
  new.dominated <- 0
  current.pareto <- as.data.frame(current.pareto)
  new.pareto <- as.data.frame(new.pareto)
  
  obj <- c( ifelse(maximize[1], 1, -1), ifelse(maximize[2], 1, -1) )
  
  for(i in 1:nrow(new.pareto)){
    for(j in 1:nrow(current.pareto)){
      if((obj[1]*new.pareto[i,"f1"] < obj[1]*current.pareto[j,"f1"]) && (obj[2]*new.pareto[i,"f2"] < obj[2]*current.pareto[j,"f2"])){
        old.dominated <- old.dominated + 1
      }else if((obj[1]*new.pareto[i,"f1"] > obj[1]*current.pareto[j,"f1"]) && (obj[1]*new.pareto[i,"f2"] > obj[1]*current.pareto[j,"f2"])){
        new.dominated <- new.dominated + 1
      }
    } 
  }
  
  convergence.index <- old.dominated/nrow(current.pareto) - new.dominated/nrow(new.pareto)
  #print(paste0("Convergence index: ", convergence.index))
  return(convergence.index)
}

dominance_ranking_sorting <- function(population, obj.values){
  if(nrow(population) == 1){
    population <- cbind(population, obj.values, "rnkIndex"=1, "cd.density"=1)
    return(population)
  }
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

fitness_selection_crowding_distance <- function(R, P.size, K){
  if(nrow(R) <= P.size){
    #print("Population size is lower than original population, returning...")
    return(R)
  }
  R <- as.data.frame(R)
  last.rank <- R[P.size, "rnkIndex"]
  if(R[P.size + 1, "rnkIndex"] > last.rank){ # No need to check if there's dominance
    return(R[1:P.size, ])
  }
  rank <- last.rank
  if(last.rank == 1){
    last.ranking.solutions <- R[R$rnkIndex == 1, ]
    last.ranking.solutions <- last.ranking.solutions[order(last.ranking.solutions$cd.density, decreasing=TRUE), ]
    last.ranking.solutions <- last.ranking.solutions[1:P.size, ]
    return(last.ranking.solutions)
  }else{
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
}

fitness_selection_diversity_metric <- function(R, groups, P.size, K, metric){
  if(nrow(R) <= P.size){
    #print("Population size is lower than original population, returning...")
    return(R)
  }
  R <- as.data.frame(R)
  last.rank <- R[P.size, "rnkIndex"]
  if(R[P.size + 1, "rnkIndex"] > last.rank){ # No need to check if there's dominance
    return(R[1:P.size, ])
  }
  rank <- last.rank
  if(last.rank == 1){ # If there are only rank 1 solutions in pareto, do this.
    last.ranking.solutions.index <- which(R$rnkIndex == 1)
    last.ranking.solutions <- R[last.ranking.solutions.index, ]
    groups <- groups[last.ranking.solutions.index]
    diversity.matrix <- calculate_diversity_matrix(groups, metric)
    # Calculate mean distance from a solution to every other, and order it by decreasing
    mean.solution.distance <- apply(diversity.matrix, 1, function(x) mean(x))
    last.ranking.solutions <- last.ranking.solutions[order(mean.solution.distance, decreasing=TRUE), ]
    last.ranking.solutions <- last.ranking.solutions[1:P.size, ]
    return (last.ranking.solutions)
  }else{ # Else, search for last pareto frontier and only compare those.
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
}

population_mating_and_mutation <- function(mating_pool, num.genes, params, P.size=NULL){
  K <- params$K
  if(missing(P.size)){
    P.size <- params$popSize
  }
  
  Q <- matrix(0, nrow = P.size, ncol = K)
  mat.rate <- params$mating_rate
  mut.rate <- params$mutation_rate
  genes <- 1:num.genes
  for(p in 1:P.size){
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

diverse_population_mating_and_mutation <- function(mating_pool, distances, groups, params, P.size=NULL){
  K <- params$K
  if(missing(P.size)){
    P.size <- params$popSize
  }
  
  #D <- distances$exp.dist 
  D <- distances$comp.dist
  Q <- matrix(0, nrow = P.size, ncol = K)
  mat.rate <- params$mating_rate
  mut.rate <- params$mutation_rate
  
  in.density.radius <- list()
  density.radius <- params$mutation_radius
  
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
        selected <- FALSE
        # Define a radius and randomly select a gene to mutate.
        density.radius <- params$mutation_radius
        # Grow radius is no gene is found nearby.
        while(!selected){
          in.radius <- which(distances$comp.dist[gene.chr.1, ] <= density.radius)
          in.radius <- which(!(in.radius %in% in.density.radius))
          if(length(in.radius) > 0){
            gene <- sample(in.radius, 1) 
            selected <- TRUE
          }else{
            density.radius <- density.radius*1.1
          }
        }
        #group.chr1 <- groups[[p]][gene.chr.1]
        #elems.group <- which(groups[[p]] == group.chr1)
        #gene <- sample(elems.group, 1)
        #print(paste("Gene ", k, "mutated"))
      }
      # This only occurs in rare cases when random value was already in offspring chromosome
      # Adding a while makes sure no error can occur
      while(is.element(gene, Q[p, ])){ 
        #*** More diversity criteria can be added here ***
        #group.chr1 <- groups[[p]][gene.chr.1]
        #elems.group <- which(groups[[p]] == group.chr1)
        #gene <- sample(elems.group, 1)
        gene <- sample(genes, 1)
      }
      Q[p, k] <- gene
      # Remember the selected gene neighbors, so a mutation cant be assigned those
      in.radius <- which(D[gene, ] <= density.radius)
      in.density.radius <- c(in.density.radius, in.radius)
      # If almost all of the genes are marked, reset.
      if(length(in.density.radius)/nrow(D) > 0.8){
        in.density.radius <- list()
      }
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
      diversity.matrix[j, i] = diversity.matrix[i, j]
    }
  }
  return(diversity.matrix)
}

diversity_metric <- function(group1, group2, metric){
  if(metric == "NMI"){
    return(normalized_mutual_information(group1, group2))
  }else{
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
      warning(paste("ERROR: Metric", metric, "not supported!"))
      return(FALSE) 
    }
  } 
  #conf.matrix <- matrix(c(yy, ny, yn, nn), nrow=2, ncol=2)
  #return(conf.matrix)
}

jaccard_index <- function(yy, ny, yn, nn){
  #yy <- matrix[1,1]
  #ny <- matrix[2,1]
  #yn <- matrix[1,2]
  #nn <- matrix[2,2]
  return( round(1 - (yy / (yy + yn + ny)), 3) )
}

# Do not use rand index
rand_index <- function(yy, ny, yn, nn){
  return ( round(1 - ((yy + nn) / (yy + yn + ny + nn)), 3) )
}

# Not working
normalized_mutual_information <- function(A.groups, B.groups){
  N <- length(A.groups)
  tab.a <- table(A.groups)
  P.a <- unname(apply(tab.a, 1, function(x) x/N))
  tab.b <- table(B.groups)
  P.b <- unname(apply(tab.b, 1, function(x) x/N))
  H_A <- -1*sum(P.a*log(P.a))
  H_B <- -1*sum(P.a*log(P.a))
  
  NMI_den <- sqrt(H_A*H_B)
  
  
  #K <- max(A.groups)
  NMI_num <- infotheo::mutinformation(A.groups, B.groups)
  
  #for(j in 1:K){
  #  for(k in 1:K){
  #    P.ab <- sum(A.groups == j & B.groups == k)/N
  #    print(paste0(j, ",", k))
  #    print(P.ab)
  #    NMI_num <- NMI_num + P.ab*log( P.ab / ( P.a[j]*P.b[k] ) ) 
  #    print("NMI")
  #    print(NMI_num)
  #  }
  #}
  return( round(1 - (NMI_num/NMI_den), 3) )
}

#### Evaluation and result plotting functions ####

# Plot partial pareto front
# **** Use normalization of all solutions *****
plot_pareto <- function(old.pareto, new.pareto, generation, output.path, agent=NULL, phase=NULL){
  #print(paste("Pareto front for generation", generation))
  #print(new.pareto)
  output.path <- file.path(output.path, "pareto")
  #new.pareto <- new.pareto[, c("f1", "f2")]
  #old.pareto <- old.pareto[, c("f1", "f2")]
  #new.pareto <- as.data.frame(eaf::normalise(new.pareto[, c("f1", "f2")], c(0,1)))
  #old.pareto <- as.data.frame(eaf::normalise(old.pareto[, c("f1", "f2")], c(0,1)))
  new.pareto <- normalise_pareto(new.pareto[, c("f1", "f2")])
  old.pareto <- normalise_pareto(old.pareto[, c("f1", "f2")])
  
  # Define an scale function and normalise
  #scaler <- function(x){ (x-min(x))/(max(x)-min(x)) }
  if(nrow(new.pareto) == 1 && nrow(old.pareto) == 1){ # If new pareto has 1 solution, normalise in reference
    temp <- rbind(old.pareto, new.pareto)
    #temp <- as.data.frame(eaf::normalise(temp, c(0,1)))
    temp <- normalise_pareto(temp)
    old.pareto <- temp[1, ]
    new.pareto <- temp[2, ]
  }else if(nrow(old.pareto) == 1){
    temp <- rbind(old.pareto, new.pareto)
    temp <- normalise_pareto(temp)
    old.pareto <- temp[1, ]
  }else if(nrow(new.pareto) == 1){
    temp <- rbind(new.pareto, old.pareto)
    temp <- normalise_pareto(temp)
    new.pareto <- temp[1, ]
  }
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
    xlim(0, 1) + #max(pareto$f1)) +
    ylim(0, 1)#max(pareto$f2))  
  )
  if(missing(agent)){
    dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output.path, paste0("g_", generation, ".png"))
  }else if(missing(phase)){
    folder <- file.path(output.path, paste0("agent_", agent))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(folder, paste0("g_", generation,".png")) 
  }else{
    folder <- file.path(output.path, paste0("agent_", agent))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(folder, paste0("_phase_", phase, "_g_", generation, ".png")) 
  }
  suppressWarnings(ggsave(filename, height=5, width=7))
}

plot_phase_population <- function(solutions, phase, output.path, is_final=FALSE, show.all=FALSE){
  output.path <- file.path(output.path, "pareto")
  if(missing(is_final) || is_final==FALSE){
    plot.title <- paste0("Phase ", phase)
    folder <- file.path(output.path, "phases")
    filename <- file.path(folder, paste0("phase_", phase, ".png"))
  }else{
    plot.title <- "Diverse Memetic NSGA-II results"
    folder <- file.path(output.path)
    filename <- file.path(output.path, "pareto_results.png") 
  }
  ranks <- solutions[, "rnkIndex"]
  #population <- as.data.frame(eaf::normalise(solutions[ , c("f1", "f2")], c(0,1)))
  population <- normalise_pareto(solutions[ , c("f1", "f2")])
  population <- cbind(population, ranks)
  colnames(population) <- c("f1", "f2", "Solutions")
  #pareto <- solutions[ranks == 1, ]
  population[ranks == 1, "Solutions"] <- "Pareto front"
  population[ranks != 1, "Solutions"] <- "Dominated"
  population$Solutions <- as.factor(population$Solutions) 
  if(show.all){
    color.values <- c("#6E4F59", "#F24949")
  }else{
    population <- population[population$Solutions == "Pareto front", ]
    color.values <- c("#F24949")
  }
  alpha <- ifelse(population$Solutions == "Pareto front", 1, 0.5)
  #pareto$color <- rep("grey", nrow(pareto))
  #** TODO: Add evaluation metric as subtitle in labels x and y 
  suppressWarnings(ggplot(population, aes(x=f1, y=f2, colour=Solutions)) +
                     labs(title=plot.title, x="Expression Index", y="Biological Index") +
                     geom_point(aes(color = Solutions), size=2, alpha=alpha) +
                     scale_color_manual(values=color.values) +
                     scale_alpha_continuous(guide=FALSE) +
                     xlim(0, 1) + #max(pareto$f1)) +
                     ylim(0, 1)#max(pareto$f2))  
  )

  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(ggsave(filename, height=5, width=7))
}

# exp.path: path until database name
plot_experiment_results <- function(exp.path){
  
  folder.path <- file.path(exp.path)
  plot.data <- read.table(file.path(folder.path,"plot_data.csv"), sep=",", header = TRUE, row.names=NULL)
  
  ggplot(plot.data, aes(x=Dataset, y=Hypervolume, fill=Dataset)) +
     labs(title="Hypervolume comparison") +
     geom_boxplot()
  
  ggsave(file.path(folder.path, "hv_results.png"), width = 6, height = 4)
  
  datasets <- list.dirs(path=folder.path, full.names=FALSE, recursive = FALSE)
  #datasets <- datasets[datasets %in% c("arabidopsis", "cell_cycle", "serum", "sporulation")]
  pareto <- list()
  for(i in 1:length(datasets)){
    dataset <- datasets[i]
    eval.file <- file.path(folder.path, dataset, "evaluations.csv")
    evaluations <- read.table(eval.file, sep=",", header = TRUE, row.names=NULL)
    evaluations <- cbind(evaluations, "Algorithm"=rep(dataset, nrow(evaluations)))
    if(i == 1){
      plot.data.metrics <- evaluations
    }else{
      plot.data.metrics <- rbind(plot.data.metrics, evaluations)
    }
  }
  ggplot(plot.data.metrics, aes(x=Dataset, y=avg_sil, fill=Dataset)) +
    labs(title="Average cluster silhouette comparison", y="Average sil.") +
    geom_boxplot()
  
  ggsave(file.path(folder.path, "sil_results.png"), width = 6, height = 4)
  
  ggplot(plot.data.metrics, aes(x=Dataset, y=delta, fill=Dataset)) +
    labs(title="Delta metric comparison", y="Delta spread index") +
    geom_boxplot()
  
  ggsave(file.path(folder.path, "delta_results.png"), width = 6, height = 4)
  
  ggplot(plot.data.metrics, aes(x=Dataset, y=time, fill=Dataset)) +
    labs(title="Execution time comparison", y="Time (seconds)") +
    geom_boxplot()
  
  ggsave(file.path(folder.path, "time_results.png"), width = 6, height = 4)
}

plot_algorithm_comparison <- function(exp.path, plot.data=NULL){
  if(!missing(plot.data)){
    dir.create(file.path(exp.path, "figures"), recursive=TRUE, showWarnings = FALSE)
    ggplot(plot.data, aes(x=Dataset, y=Hypervolume, fill=Algorithm)) +
      labs(title="Hypervolume values by algorithm", y="Hypervolume") +
      geom_boxplot() +
      facet_wrap(~Dataset, scale="free")
    ggsave(file.path(exp.path, "figures", "hv_results.png"), height=7, width=7)
    
    ggplot(plot.data, aes(x=Dataset, y=Silhouette, fill=Algorithm)) +
      labs(title="Silhouette index values by algorithm", y="Average population silhouette") +
      geom_boxplot() +
      facet_wrap(~Dataset, scale="free")
    ggsave(file.path(exp.path, "figures", "sil_results.png"), height=7, width=7)
    
    ggplot(plot.data, aes(x=Dataset, y=Delta, fill=Algorithm)) +
      labs(title="Delta index values by algorithm", y="Delta index") +
      geom_boxplot() +
      facet_wrap(~Dataset, scale="free")
    ggsave(file.path(exp.path, "figures", "delta_results.png"), height=7, width=7)
  }else{
    folder.path <- file.path(exp.path)
    algorithms <- list.dirs(path=folder.path, full.names=FALSE, recursive = FALSE)
    for(i in 1:length(algorithms)){
      algorithm <- algorithms[i]
      data.file <- file.path(folder.path, algorithm, "plot_data.csv")
      evaluations <- read.table(data.file, sep=",", header = TRUE, row.names=NULL)
      if(i == 1){
        plot.data.metrics <- evaluations
      }else{
        plot.data.metrics <- rbind(plot.data.metrics, evaluations)
      }
    }
    ggplot(plot.data.metrics, aes(x=Dataset, y=Hypervolume, fill=Algorithm)) +
      labs(title="Hypervolume values by algorithm", y="Hypervolume") +
      geom_boxplot() +
      facet_wrap(~Dataset, scale="free")
    
    dir.create(file.path(folder.path, "figures"), recursive=TRUE, showWarnings = FALSE)
    ggsave(file.path(folder.path, "figures", "hv_results.png"), height=5, width=7)
  }
}

plot_algorithm_comparison_diversity <- function(exp.path, plot.data){
  dir.create(file.path(exp.path, "figures"), recursive=TRUE, showWarnings = FALSE)
  ggplot(plot.data[plot.data$Metric=="jaccard", ], aes(x=Dataset, y=Diversity, fill=Algorithm)) +
    labs(title="Diversity in pareto front per algorithm", subtitle="Jaccard dissimilarity", y="Average separation in pareto front") +
    geom_boxplot() +
    facet_wrap(~Dataset, scale="free")
  ggsave(file.path(exp.path, "figures", "diversity_results_jaccard.png"), height=7, width=7)
  
  ggplot(plot.data[plot.data$Metric=="NMI", ], aes(x=Dataset, y=Diversity, fill=Algorithm)) +
    labs(title="Diversity in pareto front per algorithm", subtitle="NMI dissimilarity", y="Average separation in pareto front") +
    geom_boxplot() +
    facet_wrap(~Dataset, scale="free")
  ggsave(file.path(exp.path, "figures", "diversity_results_NMI.png"), height=7, width=7)
  
  ggplot(plot.data, aes(x=Dataset, y=Cluster_Ratio, fill=Algorithm)) +
    labs(title="Average cluster ratio in pareto front per algorithm", y="Ratio of clusters to number of solutions") +
    geom_boxplot() +
    facet_wrap(~Dataset, scale="free")
  ggsave(file.path(exp.path, "figures", "clust_ratio_results.png"), height=7, width=7)
  
  ggplot(plot.data[plot.data$Metric=="NMI", ], aes(x=Dataset, y=Cluster_Ratio, fill=Algorithm)) +
    labs(title="Average cluster ratio in pareto front per algorithm", y="Ratio of clusters to number of solutions") +
    geom_boxplot() +
    facet_wrap(~Dataset, scale="free")
  ggsave(file.path(exp.path, "figures", "clust_ratio_results_NMI.png"), height=7, width=7)
}

plot_algorithm_comparison_pareto <- function(exp.path){
  folder.path <- file.path(exp.path)
  algorithms <- list.dirs(path=folder.path, full.names=FALSE, recursive = FALSE)
  plot.data <- as.data.frame(matrix(nrow=0, ncol=5))
  colnames(plot.data) <- c("f1", "f2", "rnkIndex", "Algorithm", "Dataset")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){
      next
    }
    datasets <- list.dirs(path=file.path(folder.path, algorithm), recursive = FALSE, full.names=FALSE)
    for(j in 1:length(datasets)){
      pareto.dataset <- as.data.frame(matrix(nrow=0, ncol=2)) #** N of objectives hardcoded! **
      dataset <- datasets[j]
      dataset.path <- file.path(folder.path, algorithm, dataset)
      experiments <- list.dirs(path=dataset.path, recursive = FALSE, full.names=FALSE)
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
  
  dir.create(file.path(folder.path, "figures"), recursive=TRUE, showWarnings = FALSE)
  for(i in 1:length(datasets)){
    dataset <- datasets[i]
    data <- plot.data[plot.data$Dataset == dataset, ]
    data.norm <- data
    #data.norm[, c("f1", "f2")] <- eaf::normalise(data.norm[, c("f1", "f2")], c(0,1))
    data.norm[, c("f1", "f2")] <- normalise_pareto(data.norm[, c("f1", "f2")])
    
    ggplot(data, aes(x=f1, y=f2, color=Algorithm)) +
      labs(title=paste0("Pareto front for dataset: ",dataset), x="Genetic expression", y="Biological function") +
      geom_point() +
      geom_line() +
      theme(legend.position="top")
    
    ggsave(file.path(folder.path, "figures", paste0("pareto_comparison_",dataset ,".png")), height=7, width=7) 
    
    ggplot(data.norm, aes(x=f1, y=f2, color=Algorithm)) +
      labs(title=paste0("Normalised pareto front for dataset: ",dataset), x="Genetic expression", y="Biological function") +
      geom_point() +
      geom_line() +
      theme(legend.position="top") +
      xlim(0, 1) +
      ylim(0, 1)
    #facet_wrap(~Dataset, scale="free")
    
    ggsave(file.path(folder.path, "figures", paste0("pareto_comparison_",dataset ,"_norm.png")), height=7, width=7) 
  }
}

## Evaluation Metrics

evaluate_solutions <- function(population, clustering, distances, K, objDim, obj_maximize, output.base, exp.id, algorithm.name, dataset.name, time=-1, plot=FALSE){
  output <- file.path(output.base, exp.id) # Path of overall results from instance (like normalization limits)
  #if(plot){
  #  output.plots <- file.path(output, exp.id) # Path only used when plotting an individual instance's results plots
    #dir.create(file.path(output, exp.id, "clustering"), recursive = TRUE, showWarnings = FALSE)
  #}

  pareto <- population[population$rnkIndex == 1, ]
  obj.index <- (K+1):(K+objDim)
  N <- nrow(pareto)
  #dist.compounded <- distances$comp.dist #params$alpha*distances$exp.dist + (1-params$alpha)*distances$bio.dist
  #gene.dist <- distances$exp.dist 
  gene.dist <- distances$comp.dist 
  #gene.dist2 <- distances$comp.dist
  
  # Create silhouette results table
  labels <- rep(NA, K)
  for(i in 1:K){ labels[i] <- paste("Sil. cluster", i) }
  labels[K+1] <- "Average sil"
  silhouette.res <- data.frame(matrix(ncol=K+1, nrow=N))
  colnames(silhouette.res) <- labels
  if(plot){
    dir.create(file.path(output, "clustering"), recursive = TRUE, showWarnings = FALSE)
  }
  for(i in 1:N){
    # Calculate Silhouette metric
    sil <- cluster::silhouette(clustering[[i]], gene.dist)
    pam.res <- cluster::pam(gene.dist, metric = "euclidean", k = K, medoids = pareto[i, 1:K], do.swap = FALSE)
    plot.pam <- factoextra::fviz_cluster(pam.res, geom = "point")
    plot.sil <- factoextra::fviz_silhouette(sil, print.summary=FALSE)
    
    if(plot){
      out.file <- file.path(output, "clustering", paste0("p", i, "_pam.png"))
      png(out.file)
      print(plot.pam)
      dev.off()
      out.file.2 <- file.path(output, "clustering", paste0("p", i, "_sil.png"))
      png(out.file.2)
      print(plot.sil)
      dev.off() 
    }
    res <- summary(sil)
    
    silhouette.res[i, 1:K] <- unname(res$clus.avg.widths)
    silhouette.res[i, K+1] <- res$avg.width
  }
  rownames(silhouette.res) <- 1:N
  
  # Calculate Delta spread
  delta <- delta_spread(pareto, obj.index)
  
  results <- data.frame("silhouette"=silhouette.res, "delta"=delta)#, "hypervolume"=hv.res)
  # Save silhouette results
  if(plot){
    write.csv(silhouette.res, file.path(output, "clustering", "silhouette.csv"), row.names=FALSE) 
  }
  # Append this experiment evaluations (row) to instance's file (data.frame)
  eval.file <- file.path(output.base, "evaluations.csv")
  #print("At evaluate_solutions")
  #print(paste(exp.id, dataset.name, mean(silhouette.res[, "Average sil"]), delta, time))
  res <- data.frame("exp_id"=exp.id, "Dataset"=dataset.name, "avg_sil"=mean(silhouette.res[, "Average sil"]), "delta"=delta, "time"=time)#, "hypervolume"=hv.res)
  if(file.exists(eval.file)){
    write.table(res, file = eval.file, sep = ",", append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
  }else{
    #dir.create(file.path(output), recursive = TRUE, showWarnings = FALSE)
    write.table(res, file = eval.file, sep = ",", append = FALSE, quote = FALSE,
                col.names = TRUE, row.names = FALSE)
  }
  
  id <- file.path(paste0(exp.id,".csv"))
  #if(!dir.exists(output.base)){
  #  dir.create(file.path(output.base), recursive = TRUE, showWarnings = FALSE)
  #}
  write.table(pareto[, obj.index], file = file.path(output, id), sep=",", append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(pareto[, 1:K], file = file.path(output, "population.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  #write.csv(delta, file.path(output.path, "delta.csv"), row.names=FALSE)
  #write(hv.res, file.path(output.path, "hypervolume.txt"))
  return(list("results"=res, "pareto"=pareto[, obj.index]))
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

#Dev's Delta spread indicator
delta_spread <- function(pareto, obj.index){
  pareto <- pareto[, obj.index]
  N <- nrow(pareto)
  if(N < 2){
    warning("Delta spread can't be calculated from less than 2 solutions.")
    return(0)
  }
  order.f1 <- pareto[order(pareto$f1, decreasing = TRUE), ]
  order.f2 <- pareto[order(pareto$f2, decreasing = TRUE), ]
  
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  
  extreme.f1 <- order.f1[1, ]
  boundary.f1 <- order.f1[2, ]
  extreme.f2 <- order.f2[1, ]
  boundary.f2 <- order.f2[2, ]
  
  df <- euclidean(extreme.f1, boundary.f1)
  dl <- euclidean(extreme.f2, boundary.f2)
  
  D <- 0
  d <- array(dim=N-1)
  for(i in 1:(N-1)){
    d[i] <- euclidean(pareto[i+1, ], pareto[i, ])
  }
  d.mean <- mean(d)
  for(i in 1:(N-1)){
    D <- D + abs(d[i] - d.mean)
  }
  
  delta <- (df + dl + D) / (df + dl + (N-1)*d.mean)
  return(delta)
  
}


#### Memetic algorithm functions #####

diverse_fitness_sync <- function(Agent.A, Agent.B, diverse.metric, obj_indexes, pop_limit){
  Pop.A <- Agent.A$population[FALSE, ]
  Pop.B <- Agent.B$population[FALSE, ]
  Clust.A <- Agent.A$clustering
  Clust.B <- Agent.B$clustering
  p <- length(Clust.A)
  q <- length(Clust.B)
  # Create a distance matrix between chromosomes of the union of both populations
  clust <- c(Clust.A, Clust.B)
  rownames(Agent.A$population) <- 1:p
  rownames(Agent.B$population) <- (p+1):(p+q)
  Pop <- rbind(Agent.A$population, Agent.B$population)
  distance.matrix <- calculate_diversity_matrix(clust, diverse.metric)
  distance.matrix <- as.dist(distance.matrix)
  hc <- hclust(distance.matrix, method="average")
  pop.membership <- cutree(hc, k=2)

  Pop.A <- Pop[pop.membership == 1, ]
  Pop.B <- Pop[pop.membership == 2, ]
  Clust.A <- clust[pop.membership == 1]
  Clust.B <- clust[pop.membership == 2]
  A.rows <- rownames(Pop.A)
  B.rows <- rownames(Pop.B)
  # Population A non-dominated sorting
  objectives <- as.matrix(Pop.A[, obj_indexes])
  ranking <- nsga2R::fastNonDominatedSorting(objectives)
  rankIndex <-integer(nrow(Pop.A))
  i <- 1
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  Pop.A[, "rnkIndex"] <- rankIndex
  Pop.A <- Pop.A[order(rankIndex), ]
  if(nrow(Pop.A) > pop_limit){
    Pop.A <- Pop.A[1:pop_limit, ]
  }
  # Population B non-dominated sorting
  objectives <- as.matrix(Pop.B[, obj_indexes])
  ranking <- nsga2R::fastNonDominatedSorting(objectives)
  rankIndex <-integer(nrow(Pop.B))
  i <- 1
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  Pop.B[, "rnkIndex"] <- rankIndex
  Pop.B <- Pop.B[order(rankIndex), ]
  if(nrow(Pop.B) > pop_limit){
    Pop.B <- Pop.B[1:pop_limit, ] 
  }
  Agent.A$population <- Pop.A
  Agent.A$clustering <- Clust.A[match(rownames(Pop.A), A.rows)]
  Agent.B$population <- Pop.B
  Agent.B$clustering <- Clust.B[match(rownames(Pop.B), B.rows)]
  return(list("Agent.A"=Agent.A, "Agent.B"=Agent.B))
  
}

#Checks if a chromosome is already in a population
check_in_population <- function(P, chr){
  #matches <- apply(P, 1, function(x) all(diff(match(chr, x)) == 1))
  matches <- apply(P, 1, function(x) chr %in% x)
  return ( any(apply(matches, 2, function(x) Reduce("&", x))) )
}

fitness_sync <- function(Agent.A, Agent.B, obj_maximize, obj_indexes, pop_limit){
  Pop.A <- Agent.A$population
  Pop.B <- Agent.B$population
  Clust.A <- Agent.A$clustering
  Clust.B <- Agent.B$clustering
  pop_agent_a <- nrow(Pop.A)
  pop_agent_b <- nrow(Pop.B)
  rownames(Pop.A) <- 1:pop_agent_a
  rownames(Pop.B) <- (pop_agent_a+1):(pop_agent_a+pop_agent_b)
  
  # Create a factor 1 from maximization or -1 for minimization
  obj <- c( ifelse(obj_maximize[1], 1, -1), ifelse(obj_maximize[2], 1, -1) )
  
  # Check Agent B giving solutions to Agent A
  worst.A <- Pop.A[pop_agent_a, ]
  max.taken <- round(pop_agent_a/2) # Limit half of population to preserve solutions
  taken <- 0
  for(i in 1:pop_agent_b){
    if(obj[1]*Pop.B[i, "f1"] > obj[1]*worst.A$f1 && obj[2]*Pop.B[i, "f2"] > obj[2]*worst.A$f2 ){
      if(check_in_population(Pop.A[, 1:K], Pop.B[i, 1:K]) == TRUE){
        next
      }
      Pop.A <- rbind(Pop.A, Pop.B[i,])
      Clust.A <- c(Clust.A, list(Clust.B[[i]]))
      taken <- taken + 1
      if(taken == max.taken){ break }
    }
  }
  worst.B <- Pop.B[pop_agent_b, ]
  max.taken <- round(pop_agent_b/2)
  taken <- 0
  # Check Agent A giving solutions to Agent B
  for(i in 1:pop_agent_b){
    if(obj[1]*Pop.A[i, "f1"] > obj[1]*worst.B$f1 && obj[2]*Pop.A[i, "f2"] > obj[2]*worst.B$f2 ){
      if(check_in_population(Pop.B[, 1:K], Pop.A[i, 1:K]) == TRUE){
        next
      }
      Pop.B <- rbind(Pop.B, Pop.A[i,])
      Clust.B <- c(Clust.B, list(Clust.A[[i]]))
      taken <- taken + 1
      if(taken == max.taken){ break }
    }
  }
  A.rows <- rownames(Pop.A)
  B.rows <- rownames(Pop.B)
  # Population A non-dominated sorting
  objectives <- as.matrix(Pop.A[, obj_indexes])
  ranking <- nsga2R::fastNonDominatedSorting(objectives)
  rankIndex <-integer(nrow(Pop.A))
  i <- 1
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  Pop.A[, "rnkIndex"] <- rankIndex
  Pop.A <- Pop.A[order(rankIndex), ]
  if(nrow(Pop.A) > pop_limit){
    Pop.A <- Pop.A[1:pop_limit, ]
  }
  # Population B non-dominated sorting
  objectives <- as.matrix(Pop.B[, obj_indexes])
  ranking <- nsga2R::fastNonDominatedSorting(objectives)
  rankIndex <-integer(nrow(Pop.B))
  i <- 1
  while (i <= length(ranking)) {
    rankIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  Pop.B[, "rnkIndex"] <- rankIndex
  Pop.B <- Pop.B[order(rankIndex), ]
  if(nrow(Pop.B) > pop_limit){
    Pop.B <- Pop.B[1:pop_limit, ] 
  }

  Agent.A$population <- Pop.A
  Agent.A$clustering <- Clust.A[match(rownames(Pop.A), A.rows)]
  Agent.B$population <- Pop.B
  Agent.B$clustering <- Clust.B[match(rownames(Pop.B), B.rows)]
  
  return(list("Agent.A"=Agent.A, "Agent.B"=Agent.B))
}


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

diverse_memetic_nsga2 <- function(distances, params, output.path, debug=FALSE, plot=FALSE){
  evaluations <- params$evaluations
  evaluation.count <- 0
  diversity.metric <- params$diversity_metric
  diversity.level <- params$diversity_level
  K <- params$K 
  agents <- params$agents
  P.size <- params$popSize
  num.genes <- distances$n.genes
  phases <- params$phases + 1
  phase <- 1
  Agents <- list()
  
  #output.path <- get_new_dirname(output.path)
  dir.create(file.path(output.path), recursive = TRUE, showWarnings = FALSE)
  
  if(debug){
    output.log.file <- file.path(output.path, "log.txt")
    sink(output.log.file, append=FALSE)
    writeLines(c(""), output.log.file)
    Log(paste("Initiating DNSGA-2 with diversity level",diversity.level, "and metric", diversity.metric,"...")) 
  }
  # Initialize population
  if(params$is_random_population){
    P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed) # Random population
  }else if(diversity.level >= 1){
    P <- generate_diverse_initial_pop(distances, params, p.size=P.size, diverse_population = TRUE)
  }else{
    P <- generate_initial_pop(params$popSize, K, distances$n.genes, params$seed) # Random population
    #P <- generate_diverse_initial_pop(distances, params, p.size=P.size, diverse_population = FALSE)
  }
  
  g <- 1 # Current generation
  ## Initialize and evaluate population P
  P.data <- cluster_data(distances, P, params$alpha)
  P <- P.data$population
  P.clustering.groups <- P.data$clustering.results
  P <- evaluate_population(P, distances, P.clustering.groups, params)
  evaluation.count <- evaluation.count + nrow(P)
  P.clustering.groups <- P.clustering.groups[as.numeric(rownames(P))]
  rownames(P) <- 1:P.size
  current.pareto.front <- P[P$rnkIndex == 1, ] # Current pareto front
  
  if(plot){
    plot_phase_population(P, 0, output.path, show.all=TRUE)
  }
  generations.no.changes <- 0
  has.converged <- FALSE
  
  pop.per.agent <- floor(P.size/agents)
  evaluations.per.phase <- floor(evaluations/phases)
  evaluations.per.agent <- floor(evaluations.per.phase/agents)
  # Assign each agent with a random set of solutions from the solution pool
  solutions <- rep(1:agents, each=pop.per.agent)
  solutions <- sample(solutions)
  solutions <- split(1:P.size, solutions)
  for(i in 1:agents){
    Agents[[i]] <- list("population"=P[unlist(solutions[i]), ], "clustering"=P.clustering.groups[unlist(solutions[i])])
  }
        
  cores <- parallel::detectCores()
  if(cores < agents){
    warning(paste0("[Warning] Number of cores (", cores ,") is lower than agent number (",agents,"). Lower agent number is suggested"))
    cl <<- parallel::makeCluster(cores)
  }else{
    cl <<- parallel::makeCluster(agents)
  }
  doParallel::registerDoParallel(cl)
  
  source("dmmoea_functions.R")
  source("dmmoea_libraries.R")
  exported_func <- c("dnsga2_agent")
  #exported_func <- c("dnsga2", "cluster_data", "evaluate_population", "evaluate_xie_beni", "Log")
  #clusterExport(cl, c("agents", "phases", exported_func))
  parallel::clusterExport(cl, exported_func)
  while (phase <= phases) {
    if(debug){
      Log(paste0("Evaluations used (", evaluation.count, "/", evaluations, ")"))
    }
    Agents <- foreach(i=1:params$agents, .combine = c, .export=c("pop.per.agent", "evaluations.per.agent","phase"), .inorder = FALSE) %dopar% {
      # For every agent, call an NSGA-II procedure
      pareto_agent <- dnsga2_agent(distances, params, output.path, P.size=pop.per.agent, agent=i, phase=phase, 
                                   evaluations=evaluations.per.agent, 
                                   initial_population=Agents[[i]], debug=debug, plot=plot)
      #parameters: distances, diversity.metric, diversity.level, params, output.path, generations, P.size, agent, phase, initial_population
      return( list( pareto_agent ) )
    }
    evaluation.count <- evaluation.count + evaluations.per.phase
    if(phase != phases){
      # Aggregate solutions across agents
      print("Before aggregate")
      res <- aggregate_agents(Agents, agents, K, params$objDim)
      
      #evaluation.count <- evaluation.count + nrow(res$population)
      res$population <- remove_duplicated(res$population, K)
      print("After aggregate")
      
      # Plot current solutions for this phase
      if(plot){
        plot_phase_population(res$population, phase, output.path, show.all=TRUE)
      }
      # If syncronization is disable, keep agents as is and go to next phase
      if(params$sync_off){
        phase <- phase + 1
        next
      }
      # Remove all duplicates
      #Log("Testing population size:")
      #for(i in 1:agents){
      #  print(nrow(Agents[[i]]$population))
      #  #Agents[[i]]$population <- remove_duplicated(Agents[[i]]$population, K)
      #}
      
    }else{
      # If all phases ended, algorithm finished!
      if(debug){
        Log("Agent finished last phase!, exiting...") 
      }
      break
    }
    
    if(debug){
      Log(paste("Phase", phase, "agents finished!. Starting sync stage...")) 
    }
    # Synchronize every pair of agents
    obj_indexes <- (K+1):(K+params$objDim)
    for(i in 1:(agents-1)){
      for(j in (i+1):agents){
        if(diversity.level >= 4){
          res <- diverse_fitness_sync(Agents[[i]], Agents[[j]], diversity.metric, obj_indexes, pop.per.agent)
        }else{
          res <- fitness_sync(Agents[[i]], Agents[[j]], params$obj_maximize, obj_indexes, pop.per.agent)
        }
        Agents[[i]] <- res$Agent.A # Population agent A
        Agents[[j]] <- res$Agent.B # Population agent B
      }
    }
    if(debug){
      for(i in 1:agents){
        print(paste0("Agent ", i, " recieves:"))
        print(Agents[[i]]$population)
      }
      Log(paste("Phase", phase, "synchronization ended!"))
    }
    phase <- phase + 1
  }
  
  parallel::stopCluster(cl) # Clean workers
  closeAllConnections() # Close output log files
  
  # Aggregate solutions across agents and return
  solutions<- aggregate_agents(Agents, agents, K, params$objDim)
  solutions$population <- remove_duplicated(solutions$population, K)
  if(plot){
    plot_phase_population(solutions$population, phases, output.path, is_final=TRUE, show.all=TRUE) 
  }
  #print("Resulting population:")
  #print(solutions$population)
  return(solutions)
}

aggregate_agents <- function(Agents, n.agents, K, obj.dim){
  population <- Agents[[1]]$population
  clustering <- Agents[[1]]$clustering
  row.count <- nrow(population)
  row.names(population) <- 1:row.count
  for(i in 2:n.agents){
    row.names(Agents[[i]]$population) <- (row.count+1):(row.count+nrow(Agents[[i]]$population))
    population <- rbind(population, Agents[[i]]$population)
    clustering <- c(clustering, Agents[[i]]$clustering)
    row.count <- nrow(population)
  }
  #n.pop <- nrow(population)
  #rownames(population) <- 1:n.pop
  obj.values <- population[ , (K+1):(K+obj.dim)]
  population <- dominance_ranking_sorting(population[, 1:K], obj.values)
  clustering <- clustering[as.numeric(rownames(population))]
  return(list("population"=population, "clustering"=clustering))
}

fill_population <- function(P, K, n.genes, fill){
  filled <- 1
  filled.pop <- matrix(nrow = fill, ncol = K)
  while(filled <= fill){
    sample.chr <- sample(1:nrow(P), 1)
    chr <- P[sample.chr, 1:K]
    sample.medoid <- sample(1:K, 1)
    sample.gene <- sample(1:n.genes, 1)
    if(!is.element(sample.gene, P[sample.chr, ])){
      chr[1, sample.medoid] <- sample.gene
      filled.pop[filled,1:K] <- unlist(chr)
      filled <- filled + 1
    }
  }
  return(filled.pop)
}

remove_duplicated <- function(population, K){
  pop <- population[, 1:K]
  pop <- t(apply(pop, 1, function(x) sort(x)))
  return(population[!duplicated(pop), ])
}

normalise_pareto <- function(data, dims=2){
  max <- apply(data, 2, max)
  min <- apply(data, 2, min)
  for(i in 1:dims){
    scaler <- function(x){ (x-min[i])/(max[i]-min[i]) }
    data[, i] <- scaler(data[, i])
  }
  return(as.data.frame(data))
}

normalise_results <- function(base.path, algorithm){
  max.f1 <- 0
  min.f1 <- Inf
  max.f2 <- 0
  min.f2 <- Inf
  datasets <- list.dirs(path=file.path(base.path, algorithm), full.names=FALSE, recursive = FALSE)
  #pareto <- list()
  for(i in 1:length(datasets)){
    dataset <- datasets[i]
    experiments <- list.dirs(path=file.path(base.path, dataset), full.names=TRUE, recursive = FALSE)
    for(f in 1:length(experiments)){
      exp.name <- experiments[f] #strsplit(basename(experiments[f]), "\\(")[[1]][1]
      #pareto[[f]] <- read.csv(file.path(experiments[f], paste0(exp.name, ".csv")), header = FALSE)
      pareto <- read.table(file.path(exp.name, paste0(basename(exp.name), ".csv")), sep=",", header = FALSE, row.names=NULL)
      max.values <- apply(pareto, 2, max)
      min.values <- apply(pareto, 2, min)
      if(max.values[1] > max.f1){
        max.f1 <- unname(max.values[1])
      }
      if(max.values[2] > max.f2){
        max.f2 <- unname(max.values[2])
      }
      if(min.values[1] < min.f1){
        min.f1 <- unname(min.values[1])
      }
      if(min.values[2] < min.f2){
        min.f2 <- unname(min.values[2])
      }
    }
    limits <- data.frame("min.f1"=min.f1, "max.f1"=max.f1, "min.f2"=min.f2, "max.f2"=max.f2)
    write.table(limits, file=file.path(base.path, "limits.csv"), sep=",", append=FALSE, row.names = FALSE, quote = FALSE)
    #scaler.f1 <- function(x){ (x-min.f1)/(max.f1-min.f1) }
    #scaler.f2 <- function(x){ (x-min.f2)/(max.f2-min.f2) }
    #for(f in 1:length(files)){
    #  norm <- data.frame("f1"=scaler.f1(pareto[[f]][, 1]), "f2"=scaler.f2(pareto[[f]][, 2]))
    #  write.csv(norm, file = files[f], append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
    #}
  }
  return(limits)
}

get_normalization_limits <- function(base.path){
  #max.f1 <- 0
  min.f1 <- 0
  #max.f2 <- 0
  min.f2 <- 0
  data <- data.frame("max.f1"=NA, "max.f2"=NA)
  row <- 1
  algorithms <- list.dirs(path=file.path(base.path), full.names=FALSE, recursive = FALSE)
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    datasets <- list.dirs(path=file.path(base.path, algorithm), full.names=FALSE, recursive = FALSE)
    if(length(datasets) == 0) { next }
    #pareto <- list()
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      experiments <- list.dirs(path=file.path(base.path, algorithm, dataset), full.names=TRUE, recursive = FALSE)
      if(length(experiments) == 0){ next }
      for(k in 1:length(experiments)){
        exp.name <- experiments[k] #strsplit(basename(experiments[f]), "\\(")[[1]][1]
        if(!file.exists(exp.name, paste0(basename(exp.name), ".csv"))){ next }
        #pareto[[f]] <- read.csv(file.path(experiments[f], paste0(exp.name, ".csv")), header = FALSE)
        pareto <- read.table(file.path(exp.name, paste0(basename(exp.name), ".csv")), sep=",", header = FALSE, row.names=NULL)
        max.values <- apply(pareto, 2, max)
        data[row, ] <- round(max.values, 3)
        row <- row + 1
        #min.values <- apply(pareto, 2, min)
      }
    }
  }
  cut <- 0.9
  ggplot(data, aes(x=max.f1)) +
    geom_histogram() +
    labs(title="Distribution of expression objective function limits", x="Frecuency") +
    geom_vline(xintercept = quantile(data$max.f1, cut, na.rm=TRUE))
  
  ggsave(file.path(base.path, "limit_f1.png"), width = 6, height = 4)
  
  ggplot(data, aes(x=max.f2)) +
    geom_histogram() +
    labs(title="Distribution of biological objective function limits", x="Frecuency") +
    geom_vline(xintercept = quantile(data$max.f2, cut, na.rm=TRUE))
  
  ggsave(file.path(base.path, "limit_f2.png"), width = 6, height = 4)
  
  limits <- data.frame("min.f1"=min.f1, "max.f1"=quantile(data$max.f1, cut, na.rm=TRUE), "min.f2"=min.f2, "max.f2"=quantile(data$max.f2, cut, na.rm=TRUE))
  write.table(limits, file=file.path(base.path, "limits.csv"), sep=",", append=FALSE, row.names = FALSE, quote = FALSE)
  return(limits)
}

calculate_hypervolume <- function(pareto, point, maximize=FALSE){
  pareto <- pareto[order(pareto[, 1], decreasing=TRUE), ]
  if(maximize){
    hv <- (pareto[1, 1] - point[1])*(pareto[1, 2] - point[2])
    if(nrow(pareto) == 1){
      return(hv)
    }
    for(i in 1:(nrow(pareto)-1)){
      h <- (pareto[i+1, 1] - pareto[i, 1])
      w <- (point[2] - pareto[i+1, 2])
      hv <- hv + w*h
    }
  }else{
    hv <- (point[1] - pareto[1, 1])*(point[2] - pareto[1, 2])
    if(nrow(pareto) == 1){
      return(hv)
    }
    for(i in 1:(nrow(pareto)-1)){
      w <- (pareto[i, 1] - pareto[i+1, 1])
      h <- (point[2] - pareto[i+1, 2])
      hv <- hv + w*h
    }
  }
  return(hv)
}

# Comparation Metrics for post-tuning

inverse_generational_distance_plus <- function(S, R){
  # s: objective values pair (f1, f2) for solution of S
  # r: objective values pair (f1, f2) for solution of R
  # m: objective dimensions
  d.plus <- function(s,r) { 
    d <- sum(apply( s - r, 2, function(x) ifelse(x < 0, 0, x) ))^2 
    return(d)
  }
  IGD.plus <- 0
  for(i in 1:nrow(R)){
    r <- R[i, ]
    res <- apply(S, 1, function(x) d.plus(x, r))
    IGD.plus <- IGD.plus + min(res)
  }
  return(IGD.plus)
}

epsilon_multiplicative <- function(S, R){
  
  #for(i in 1:nrow(R)){
  #  r <- R[i, ]
  #  res <- apply(S, 1, function(x) min(max(x/r)))
  #}
  #apply(R, 1, function(x) max( apply(x, 1, function(y) min(max(y/x)))))
  inner.epsilon <- function(S, r) { apply(S, 1, function(x) max(x/r)) } 
  res <- apply(R, 1, function(x) min(inner.epsilon(S ,x)))
  return(max(res))
}

diversity_analysis <- function(P, distances, metric, exp.path=NULL, alpha=0.5, plot=FALSE){
  P <- as.matrix(P)
  res.P <- cluster_data(distances, P, alpha)
  P <- res.P$population
  P.groups <- res.P$clustering.results 
  d.matrix.P <- calculate_diversity_matrix(P.groups, metric)
  avg.dist <- mean(d.matrix.P[lower.tri(d.matrix.P, diag = FALSE)])
  d.matrix.P <- as.dist(d.matrix.P)
  min.k <- 2
  max.k <- max(nrow(P) - 2, min.k + 2) #max(round(sqrt(nrow(P)), 1), 4)
  print(paste0("(", min.k, " - " , max.k, ")"))
  if(nrow(P) < 4){
    #warning("This pareto front has too few solutions, diversity may not be accurate!!")
    return(list("diss"=d.matrix.P, "avg.dist"=avg.dist, "k.ratio"=mean(c(max.k, min.k))/nrow(P)))
  }
  
  #indexes <- c("frey", "dunn", "cindex", "silhouette", "mcclain")
  values <- list()
  for(i in 1:length(indexes)){
    best <- NbClust::NbClust(diss = d.matrix.P, distance = NULL, index=indexes[i],
                             min.nc = min.k, max.nc = max.k, method = "single")
    values[i] <- as.integer(best$Best.nc["Number_clusters"])
  }
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  best.k <- max(getmode(unlist(values)), min.k + 2)
  if(plot && !is.null(exp.path)){
    pam.res <- cluster::pam(x=d.matrix.P, diss=FALSE, k = best.k, do.swap = FALSE)
    plot.pam <- factoextra::fviz_cluster(pam.res, geom = "point", main=paste0("Pareto similarity clustering (", metric, ")"))
    dir.create(file.path(exp.path, "diversity"), recursive = TRUE, showWarnings = FALSE)
    out.file <- file.path(exp.path, "diversity", paste0("p", i, "_", metric, ".png"))
    png(out.file)
    print(plot.pam)
    dev.off()
  }
  return(list("diss"=d.matrix.P, "avg.dist"=avg.dist, "k.ratio"=best.k/nrow(P)))
}




