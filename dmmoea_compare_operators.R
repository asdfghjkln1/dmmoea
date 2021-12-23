#!bin/usr/env Rstudio
compare_operators <- function(){
  
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 3){
    print(paste0("Not enough parameters (", argnum, "/3)"))
    return(-1)
  }
  path <- args[1] #"X:\\Universidad\\dmmoea" #args[1] # "X:\\Universidad\\dmmoea" #args[1] #
  results.path <- args[2] #"Tests\\operators" # args[2] #"Tests\\runs\\XieBeni" # args[2]
  operator <- args[3] #"lv3"
  
  setwd(path)
  library(ggpubr)
  library(rstatix)
  library(ggplot2)
  library(gridExtra)
  source("dmmoea_distances.R")
  source("dmmoea_parameters.R")
  source("dmmoea_functions.R")

  runs <- 11
  results.path <- file.path(path, results.path, operator)
  dir.create(results.path, recursive = TRUE, showWarnings = FALSE)
  print("Path in:")
  print(results.path)
  params <- init_parameters()
  #params$popSize <- 20
  
  datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    distances <- load.gene.distance(dataset, params$alpha)
    if(operator == "lv1"){
      test_operator_lv1(params, distances, dataset, results.path, runs) 
    }else if(operator == "lv2"){
      test_operator_lv2(params, distances, dataset, results.path, runs)
    }else if(operator == "lv3"){
      test_operator_lv3(params, distances, dataset, results.path, runs)
    }else if(operator == "lv4"){
      test_operator_lv4(params, distances, dataset, results.path, runs)
    }
  }
  datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
  print("Finished running tests...")
  print(file.path(results.path, "plot_data.csv"))
  if(!file.exists(file.path(results.path, "plot_data.csv"))){
    warning("Limits not found, please run \"dmmoea_evaluate_results.R\" first.")
    return(-1)
  }
  if(!file.exists(file.path(results.path, "plot_data_diversity.csv"))){
    warning("Limits not found, please run \"dmmoea_evaluate_results.R\" first.")
    return(-1)
  }
  plot.pareto.norm <- read.table(file.path(results.path, "data_pareto_norm.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data <- read.table(file.path(results.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data.diversity <- read.table(file.path(results.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data$Algorithm <- as.factor(plot.data$Algorithm)
  plot.data$id <- as.factor(plot.data$id)
  figure.path <- file.path(results.path, "figures", "comparison")
  
  metrics <- c("Hypervolume", "Silhouette", "Delta")
  metrics.div <- c("Diversity", "Cluster_Ratio")
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    data <- plot.data[plot.data$Dataset == dataset, ]
    data.diversity <- plot.data.diversity[plot.data.diversity$Dataset == dataset, ]
    data.diversity.NMI <- data.diversity[data.diversity$Metric == "NMI", ]
    data.diversity.jaccard <- data.diversity[data.diversity$Metric == "jaccard", ]
    if(operator == "lv2" || operator == "lv3"){
      lapply(metrics, function(x) multi.variable.tests(data, x, "Algorithm", dataset, figure.path, paired=TRUE))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.NMI, x, "Algorithm", dataset, figure.path, paired=TRUE))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.jaccard, x, "Algorithm", dataset, figure.path, paired=TRUE))
    }else{
      lapply(metrics, function(x) multi.variable.tests(data, x, "Algorithm", dataset, figure.path))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.NMI, x, "Algorithm", dataset, figure.path))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.jaccard, x, "Algorithm", dataset, figure.path))
    }
    #res <- friedman.test.with.post.hoc(Hypervolume ~ Algorithm | id, data=data, to.plot.parallel = F, dataset=dataset)
    #print(res)
  }
}

test_operator_lv1 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Random", "Selective", "Selective_Diverse")
  for(j in 1:length(operators)){
    op <- operators[j]
    output.path <- file.path(output.folder, op, dataset)
    for(i in 1:runs){
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      if(op == "Random"){
        seed <- as.numeric(Sys.time())
        P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
      }else if(op == "Selective"){
        P  <- generate_diverse_initial_pop(distances, params, diverse_population=FALSE)
      }else if(op == "Selective_Diverse"){
        P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE)
      }
      P.data <- cluster_data(distances, P, params$alpha)
      P <- P.data$population
      P.clustering.groups <- P.data$clustering.results
      P <- evaluate_population(P, distances, P.clustering.groups, params)
      
      evaluate_solutions(P, P.clustering.groups, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp), 
                         basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    } 
  }
}

test_operator_lv2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Crowding_Distance", "Jaccard", "NMI")
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    R <- generate_initial_pop(2*params$popSize, params$K, distances$n.genes, seed) 
    R.data <- cluster_data(distances, R, params$alpha)
    R <- R.data$population
    R.clustering.groups <- R.data$clustering.results
    R <- evaluate_population(R, distances, R.clustering.groups, params)
    obj.values <- R[, (params$K+1):(params$K+params$objDim)] # Select objective values
    R.rows <- row.names(R)
    R <- dominance_ranking_sorting(R[, 1:params$K], obj.values) # Recalculate ranking
    R.clustering.groups <- R.clustering.groups[match((row.names(R)), R.rows)] # Update clustering
    row.names(R)<-1:nrow(R)
    for(j in 1:length(operators)){
      op <- operators[j]
      output.path <- file.path(output.folder, op, dataset)
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      
      ## Population fitness selection
      if(op=="Crowding_Distance"){
        P_next_generation <- fitness_selection_crowding_distance(R, params$popSize, params$K) 
      }else if(op=="Jaccard"){
        P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, params$popSize, params$K, "jaccard")
      }else if(op=="NMI"){
        P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, params$popSize, params$K, "NMI")
      }
      P.clustering.groups <- R.clustering.groups[as.numeric(row.names(P_next_generation))] # Update clustering
      row.names(P_next_generation) <- 1:nrow(P_next_generation)
      
      evaluate_solutions(P_next_generation, P.clustering.groups, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp), 
                         basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    }
  } 
}

test_operator_lv3 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Original", "Direct_Crossover", "Selective_Crossover", "Diverse_Mutation", "Combined")
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
    P.data <- cluster_data(distances, P, params$alpha)
    P <- P.data$population
    P.clustering.groups <- P.data$clustering.results
    P <- evaluate_population(P, distances, P.clustering.groups, params)
    for(j in 1:length(operators)){
      op <- operators[j]
      output.path <- file.path(output.folder, op, dataset)
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      
      mating_pool <- nsga2R::tournamentSelection(P, params$popSize, params$tourSize)
      if(op == "Original"){
        evaluate_solutions(P, P.clustering.groups, distances, params$K, 
                           params$objDim, params$obj_maximize, dirname(output.exp), 
                           basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
        next
      }else if(op == "Direct_Crossover"){
        Q <- population_mating_and_mutation(mating_pool, distances$n.genes, params)
      }else if(op == "Selective_Crossover"){
        Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, params, P.size=params$popSize, type="selective")
      }else if(op == "Diverse_Mutation"){
        Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, params, P.size=params$popSize, type="mut.only")
      }else if(op == "Combined"){
        Q <- diverse_population_mating_and_mutation(mating_pool, distances, P.clustering.groups, params, P.size=params$popSize, type="all")
      }else{
        print("Operator name not found!")
        print(op)
        return()
      }
      
      Q.data <- cluster_data(distances, Q, params$alpha)
      Q <- Q.data$population
      Q.clustering.groups <- Q.data$clustering.results
      Q <- evaluate_population(Q, distances, Q.clustering.groups, params)
      
      evaluate_solutions(Q, Q.clustering.groups, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp), 
                         basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    }
  } 
}

test_operator_lv4 <- function(params, distances, dataset, output.folder, runs=31){
  #operators <- c("Original", "Sync_Agent_1", "Sync_Agent_2", "Diverse_Agent_1", "Diverse_Agent_2")
  operators <- c("Original", "Sync_Agent", "Diverse_Agent")
  obj_indexes <- (params$K+1):(params$K+params$objDim)
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
    P.data <- cluster_data(distances, P, params$alpha)
    P <- P.data$population
    P.clustering.groups <- P.data$clustering.results
    P <- evaluate_population(P, distances, P.clustering.groups, params)
    cut <- round(params$popSize/2)
    Agent.A <- list("population"=P[1:cut, ], "clustering"=P.clustering.groups[1:cut])
    Agent.B <- list("population"=P[(cut+1):nrow(P), ], "clustering"=P.clustering.groups[(cut+1):nrow(P)])
    for(j in 1:length(operators)){
      op <- operators[j]
      output.path <- file.path(output.folder, op, dataset)
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      
      mating_pool <- nsga2R::tournamentSelection(P, params$popSize, params$tourSize)
      if(op == "Original"){
        evaluate_solutions(P, P.clustering.groups, distances, params$K, 
                           params$objDim, params$obj_maximize, dirname(output.exp), 
                           basename(output.exp), op, dataset, pareto.only=FALSE,  plot=FALSE)
        next
      }else if(op == "Sync_Agent"){
        res <- fitness_sync(Agent.A, Agent.B, params$obj_maximize, obj_indexes, cut)
        Agent.A <- res$Agent.A
        Agent.B <- res$Agent.B
        op <- paste0("Sync_Agent_", 1:2)
      }else if(op == "Diverse_Agent"){
        res <- diverse_fitness_sync(Agent.A, Agent.B, "jaccard", obj_indexes, cut)
        Agent.A <- res$Agent.A
        Agent.B <- res$Agent.B
        op <- paste0("Diverse_Agent_", 1:2)
      }else{
        print("Operator name not found!")
        print(op)
        return()
      }
      
      output.path.1 <- file.path(output.folder, op[1], dataset)
      output.path.2 <- file.path(output.folder, op[2], dataset)
      output.exp.1 <- file.path(output.path.1, i)
      output.exp.2 <- file.path(output.path.2, i)
      dir.create(output.exp.1, recursive=TRUE, showWarnings = FALSE)
      dir.create(output.exp.2, recursive=TRUE, showWarnings = FALSE)
      evaluate_solutions(Agent.A$population, Agent.A$clustering, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp.1), 
                         basename(output.exp.1), op, dataset, pareto.only=FALSE, plot=FALSE)
      evaluate_solutions(Agent.B$population, Agent.B$clustering, distances, params$K, 
                         params$objDim, params$obj_maximize, dirname(output.exp.2), 
                         basename(output.exp.2), op, dataset, pareto.only=FALSE, plot=FALSE)
    }
  }
}


multi.variable.tests <- function(data, metric, exp.group, dataset.name, output.path, paired=FALSE){
  dir.create(output.path, recursive=TRUE, showWarnings = FALSE)
  form <- as.formula(call("~", as.symbol(metric), as.symbol(exp.group)))
  if(paired == FALSE){
    res <- kruskal_test(data, formula=form) 
  }else{
    data$id <- as.factor(data$id)
    data$Algorithm <- as.factor(data$Algorithm)
    #print(data[1:20, ])
    form.friedman <- as.formula(paste0(metric, " ~ ", exp.group, " | id"))
    print(form.friedman)
    res <- friedman_test(data=data, formula=form.friedman)
    print("A")
  }
  pwc <- wilcox_test(data, formula=form, paired=paired, p.adjust.method="bonferroni")
  print("B")
  pwc <- pwc %>% add_xy_position(x = exp.group)
  w <- 4 + 0.6*(length(unique(data[, exp.group])) - 3)#ifelse(exp.group>3, 6,5)
  Y <- pwc$y.position
  gap.data <- max(data[, metric]) - min(data[, metric])
  gap <- max(Y[2] - Y[1], gap.data*0.07)
  pwc$y.position <- Y - gap*seq(from=1, to=0, length.out=length(unique(data[, exp.group])))
  
  ggplot(data, aes_string(x=exp.group, y=metric)) +
    geom_boxplot(aes_string(fill=exp.group)) +
    labs(subtitle = get_test_label(res, detailed = TRUE), 
         caption = get_pwc_label(pwc),
         title=paste0(metric, " comparison: ", dataset.name, " dataset")) +
    theme_minimal() +
    theme(strip.text.x = element_blank(),
          axis.text.x = element_text(angle=25)) +
    stat_pvalue_manual(pwc, label = "p = {p.adj}", hide.ns = TRUE)
  ggsave(file.path(output.path, paste0(metric, "_results_", dataset.name, ".png")), width = w, height = 6)
  print(paste(metric, dataset.name, "... Done."))
}

compare_pareto_front <- function(data, dataset.name, output.path){
  library(eaf)
  algorithms <- unique(data$Algorithm)
  algorithms <- algorithms[!(algorithms %in% "Ideal pareto")]
  epsilon.matrix <- matrix(0, ncol=length(algorithms),nrow = length(algorithms))
  epsilon.matrix[lower.tri(epsilon.matrix, diag=TRUE)] <- NA
  #IGD.matrix <- matrix(0, ncol=length(algorithms),nrow = length(algorithms))
  for(i in 1:(length(algorithms)-1)){
    for(j in (i+1):length(algorithms)){
      #print(paste(i,j))
      alg.1 <- algorithms[i]
      alg.2 <- algorithms[j]
      pareto.1 <- data[data$Algorithm == alg.1, 1:2]
      pareto.2 <- data[data$Algorithm == alg.2, 1:2]
      epsilon.matrix[i,j] <- eaf::epsilon_mult(pareto.1, pareto.2, maximise=FALSE)#epsilon_multiplicative(pareto.1, pareto.2)
      epsilon.matrix[j,i] <- eaf::epsilon_mult(pareto.2, pareto.1, maximise=FALSE)#epsilon_multiplicative(pareto.2, pareto.1)
      #IGD.matrix[i,j] <- inverse_generational_distance_plus(pareto.2, ideal.pareto)
    }
  }
  #print("Epsilon matrix:")
  #print(epsilon.matrix)
  epsilon.matrix <- reshape::melt(epsilon.matrix)
  ep.value <- round(epsilon.matrix$value,3)
  for(i in 1:length(algorithms)){
    epsilon.matrix[epsilon.matrix == i] <- algorithms[i]
  #  IGD.matrix[IGD.matrix == i] <- algorithms[i]
  }
  # in case that a value got accidentally replaced
  epsilon.matrix$value <- ep.value
  
  
  
  # --------- # 
  # Calculate IGD+ distance as implemented in: https://mlopez-ibanez.github.io/eaf/reference/igd.html
  IGD.table <- data.frame("Algoritmo"=algorithms, "IGD"=rep(0, length(algorithms)))
  colnames(IGD.table) <- c("Algoritmo", "IGD+")
  IGD.table$Algoritmo <- algorithms
  ideal.pareto <- data[data$Algorithm == "Ideal pareto", 1:2]
  for(i in 1:length(algorithms)){
    pareto <- data[data$Algorithm == algorithms[i], 1:2]
    IGD.table[i, "IGD+"] <- as.numeric(eaf::igd_plus(pareto, ideal.pareto, maximise=FALSE))
    #IGD.table[i, "IGD+"] <- inverse_generational_distance_plus(pareto, ideal.pareto)
  }
  IGD.table[, "IGD+"] <- as.numeric(formatC(IGD.table[,"IGD+"], format = "e", digits = 2))
  #IGD.table <- IGD.table[order(IGD.table[, "IGD+"]), ]
  order <- order(IGD.table[, "IGD+"], decreasing=FALSE)
  rnk <- integer(length(order))
  i <- 1
  while (i <= length(order)) {
    rnk[order[[i]]] <- i
    i <- i + 1
  } 
  IGD.table$Algoritmo <- paste0(algorithms, " (", rnk, ")")
  print("IGD table:")
  print(IGD.table)
  
  #par(mfrow = c(1, 2))
  p1 <- ggplot(epsilon.matrix, aes(x=X1, y=X2, fill=value)) + 
    geom_tile()+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient(low = "white", high = "red") +
    labs(title="Multiplicative Epsilon") +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(
      plot.title = element_text(size=8),
      axis.text.x = element_text(size=6, angle=25),
      axis.text.y = element_text(size=6)
    )
  p2 <- tableGrob(IGD.table, rows = NULL)
  #p2 <- ggplot(IGD.matrix, aes(x=X1, y=X2, fill=value)) + 
  #  geom_tile()+
  #  geom_text(aes(label = round(value, 3))) +
  #  scale_fill_gradient(low = "white", high = "red") +
  #  labs(title="Modified Inverted Generational Distance (IGD+)") +
  #  theme_minimal() +
  #  theme(axis.title.x = element_blank()) +
  #  theme(axis.title.y = element_blank()) +
  #  theme(
  #    plot.title = element_text(size=8),
  #    axis.text.x = element_text(size=6),
  #    axis.text.y = element_text(size=6)
  #  )
  #grid.arrange(p1, p2, ncol=2, top="Comparación entre fronteras de pareto")
  #tg <- textGrob('Title', gp = gpar(fontsize = 13, fontface = 'bold'))
  sg <- grid::textGrob(paste0("Dataset ", dataset.name), gp = grid::gpar(fontsize = 9))
  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output.path, paste0("pareto_comparison_", dataset.name, ".png")), 
         arrangeGrob(p1, p2, widths = c(5, 4), nrow=1, 
                     top="Comparación entre fronteras de pareto", 
                     bottom=sg), 
         height=2.5, width=6)
  #par(mfrow = c(1, 1))
  #grid.arrange(p1, p2, nrow = 1)
}


compare_operators()
