#!bin/usr/env Rstudio
compare_operators <- function(){
  Sys.setlocale("LC_CTYPE","spanish")
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 5){
    print(paste0("Not enough parameters (", argnum, "/5)"))
    return(-1)
  }
  path <- args[1] #"X:\\Universidad\\dmmoea" #args[1] # "X:\\Universidad\\dmmoea" #args[1] #
  results.path <- args[2] #"Tests\\operators" # args[2] #"Tests\\runs\\XieBeni" # args[2]
  operator <- args[3] #"lv3"
  runs <- as.numeric(args[4])
  generations <- as.numeric(args[5])
  
  setwd(path)
  library(ggpubr)
  library(rstatix)
  library(ggplot2)
  library(gridExtra)
  source("dmmoea_distances.R")
  source("dmmoea_parameters.R")
  source("dmmoea_functions.R")

  results.path <- file.path(path, results.path, operator)
  dir.create(results.path, recursive = TRUE, showWarnings = FALSE)
  print("Path in:")
  print(results.path)
  params <- init_parameters()
  params$evaluations <- params$popSize*(generations+1)
  
  datasets <- c("arabidopsis")#, "cell_cycle", "serum", "sporulation")
  
  #if(!file.exists(file.path(results.path, "diversity_evolution.csv"))){
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    params$dataset <- dataset
    distances <- load.gene.distance(dataset, params$alpha)
    if(operator == "lv1_v2"){
      test_operator_lv1_v2(params, distances, dataset, results.path, runs) 
    }else if(operator == "lv2_v2"){
      test_operator_lv2_v2(params, distances, dataset, results.path, runs)
    }else if(operator == "lv3_v2"){
      test_operator_lv3_v2(params, distances, dataset, results.path, runs)
    }else if(operator == "lv4_v2"){
      test_operator_lv4_v2(params, distances, dataset, results.path, runs)
    }
  }
  #}
  print("Tests finished! Plotting results...")
  data.diversity <- read.table(file.path(results.path, "diversity_evolution.csv"), sep=",", header=TRUE, row.names=NULL)
  data.diversity$Algorithm <- as.factor(data.diversity$Algorithm)
  data.diversity$id <- as.factor(data.diversity$id)
  
  plot_diversity_evolution(data.diversity, operator, results.path,paired=FALSE)
  return(0)
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
    if(operator == "lv2" || operator == "lv3" || operator == "lv4"){
      lapply(metrics, function(x) multi.variable.tests(data, x, "Algorithm", dataset, figure.path, paired=TRUE))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.NMI, x, "Algorithm", dataset, figure.path, paired=TRUE))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.jaccard, x, "Algorithm", dataset, figure.path, paired=TRUE))
    }else{
      lapply(metrics, function(x) multi.variable.tests(data, x, "Algorithm", dataset, figure.path))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.NMI, x, "Algorithm", dataset, figure.path))
      lapply(metrics.div, function(x) multi.variable.tests(data.diversity.jaccard, x, "Algorithm", dataset, figure.path))
    }
  }
}

plot_diversity_evolution <- function(data, exp.name, output.folder, paired=TRUE){
  form <- as.formula("diversity ~ Algorithm")
  data$generation <- as.factor(data$generation)
  #if(!paired){
  #  res <- kruskal_test(data, formula=form) 
  #}else{
  #  data$id <- as.factor(data$id)
  #  data$Algorithm <- as.factor(data$Algorithm)
  #  form.friedman <- as.formula("diversity ~ Algorithm | id")
  #  res <- friedman_test(data=data, formula=form.friedman)
  #}
  #pwc <- wilcox_test(data, formula=form, paired=paired, p.adjust.method="bonferroni")
  #pwc <- pwc %>% add_xy_position(x = "Algorithm")
  #w <- 4 + 0.6*(length(unique(data$Algorithm)) - 3)#ifelse(exp.group>3, 6,5)
  #Y <- pwc$y.position
  #gap.data <- max(data$diversity) - min(data$diversity)
  #gap <- max(Y[2] - Y[1], gap.data*0.07)
  #pwc$y.position <- Y - gap*seq(from=1, to=0, length.out=length(unique(data$Algorithm)))
  
  xlab <- "Generaci\U00F3n"
  if(exp.name == "lv4_v2"){
    xlab <- "Fase de comunicaci\U00F3n"
    title <- "Diversidad en DMNSGA-II seg\U00FAn operador de comunicaci\U00F3n"
  }else if(exp.name == "lv1_v2"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de poblaci\U00F3n inicial"
  }else if(exp.name == "lv2_v2"){
    title <- "Diversidad en DNSGA-II seg\U00FAn operador de selecci\U00F3n"
  }else if(exp.name == "lv3_v2"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de cruzamiento y mutaci\U00F3n"
  }
  #U+00E9
  data.se <- summarySE(data, measurevar="diversity", groupvars=c("Algorithm","generation"))
  print(head(data.se))
  
  pd <- position_dodge(0.1)
  ggplot(data.se, aes(x=generation, y=diversity, colour=Algorithm, group=Algorithm)) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=diversity-se, ymax=diversity+se), width=.1, position=pd) +
    labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", title=title) +
         #subtitle = get_test_label(res, detailed = TRUE), 
         #caption = get_pwc_label(pwc)) +
    theme_minimal()

  #ggplot(data, aes(x=generation, y=diversity)) +
  #  geom_boxplot(aes(fill=Algorithm)) +
  #  labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", 
  #       title=title) +
  #       #subtitle = get_test_label(res, detailed = TRUE), 
  #       #caption = get_pwc_label(pwc)) +
  #  theme_minimal()
  ggsave(file.path(output.folder, "diversity_evolution.png"))
}

test_operator_lv1_v2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Random", "Selective", "Selective_Diverse")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  
  for(j in 1:length(operators)){
    op <- operators[j]
    output.path <- file.path(output.folder, op, dataset)
    for(i in 1:runs){
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      if(op == "Random"){
        seed <- as.numeric(Sys.time())
        P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
      }else if(op == "Selective"){
        P  <- generate_diverse_initial_pop(distances, params, diverse_population=FALSE)
      }else if(op == "Selective_Diverse"){
        P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE)
      }
      
      nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder)
      
      #evaluate_solutions(P, P.clustering.groups, distances, params$K, 
      #                   params$objDim, params$obj_maximize, dirname(output.exp), 
      #                   basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    } 
  }
}

test_operator_lv2_v2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Crowding_Distance", "Jaccard", "NMI")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed)
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
        params$diversity_level <- 0
        P_next_generation <- dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="Jaccard"){
        params$diversity_level <- 2
        params$diversity_metric <- "jaccard"
        P_next_generation <- dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="NMI"){
        params$diversity_level <- 2
        params$diversity_metric <- "NMI"
        P_next_generation <- dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }

      
      #evaluate_solutions(P_next_generation, P.clustering.groups, distances, params$K, 
      #                   params$objDim, params$obj_maximize, dirname(output.exp), 
      #                   basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    }
  } 
}

test_operator_lv3_v2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Direct_Crossover", "Selective_Crossover", "Diverse_Mutation", "Combined")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed)
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
      if(op=="Direct_Crossover"){
        params$diversity_level <- 0
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder)
      }else if(op=="Selective_Crossover"){
        params$diversity_level <- 3
        params$mutation_type <- "selective"
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="Diverse_Mutation"){
        params$diversity_level <- 3
        params$mutation_type <- "mut.only"
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="Combined"){
        params$diversity_level <- 3
        params$mutation_type <- "all"
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }
    }
  } 
}

test_operator_lv4_v2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Elitist_sync", "Diverse_sync")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=5))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed)
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
      if(op=="Elitist_sync"){
        params$diversity_level <- 0
        diverse_memetic_nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder)
      }else if(op=="Diverse_sync"){
        params$diversity_level <- 0
        params$sync_method <- "pam"
        diverse_memetic_nsga2(distances, params, output.exp, initial_population=P, limits=NULL, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }
    }
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
        P_next_generation <- fitness_selection_crowding_distance(R, params$popSize, params$K, front.only=TRUE) 
      }else if(op=="Jaccard"){
        P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, params$popSize, params$K, "jaccard", front.only=TRUE)
      }else if(op=="NMI"){
        P_next_generation <- fitness_selection_diversity_metric(R, R.clustering.groups, params$popSize, params$K, "NMI", front.only=TRUE)
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
  operators <- c("Direct_Crossover", "Selective_Crossover", "Diverse_Mutation", "Combined")
  for(i in 1:runs){
    seed <- as.numeric(Sys.time())
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
    P.data <- cluster_data(distances, P, params$alpha)
    P <- P.data$population
    P.clustering.groups <- P.data$clustering.results
    P <- evaluate_population(P, distances, P.clustering.groups, params)
    mating_pool <- nsga2R::tournamentSelection(P, params$popSize, params$tourSize)
    for(j in 1:length(operators)){
      op <- operators[j]
      output.path <- file.path(output.folder, op, dataset)
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      if(dir.exists(output.exp)){
        next
      }
      dir.create(output.exp, showWarnings = FALSE, recursive=TRUE)
      
      
      #if(op == "Original"){
      #  evaluate_solutions(P, P.clustering.groups, distances, params$K, 
      #                     params$objDim, params$obj_maximize, dirname(output.exp), 
      #                     basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
      #  next
      #}
      if(op == "Direct_Crossover"){
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
  operators <- c("Sync_Agent", "Diverse_Agent")
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
      if(op == "Sync_Agent"){
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
  #grid.arrange(p1, p2, ncol=2, top="Comparaci\U00F3n entre fronteras de pareto")
  #tg <- textGrob('Title', gp = gpar(fontsize = 13, fontface = 'bold'))
  sg <- grid::textGrob(paste0("Dataset ", dataset.name), gp = grid::gpar(fontsize = 9))
  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output.path, paste0("pareto_comparison_", dataset.name, ".png")), 
         arrangeGrob(p1, p2, widths = c(5, 4), nrow=1, 
                     top="Comparaci\U00F3n entre fronteras de pareto", 
                     bottom=sg), 
         height=2.5, width=6)
  #par(mfrow = c(1, 1))
  #grid.arrange(p1, p2, nrow = 1)
}



## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                   conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


compare_operators()
