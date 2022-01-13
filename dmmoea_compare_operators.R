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
  #}
  print("Tests finished! Plotting results...")
  data.diversity <- read.table(file.path(results.path, "diversity_evolution.csv"), sep=",", header=TRUE, row.names=NULL)
  data.diversity$Algorithm <- as.factor(data.diversity$Algorithm)
  data.diversity$id <- as.factor(data.diversity$id)
  
  plot_diversity_evolution(data.diversity, operator, results.path)
  #plot_diversity_evolution(data.diversity, operator, metric="NMI", results.path)
}

#https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
plot_diversity_evolution_2 <- function(data, exp.name, output.folder){
  data <- data[as.numeric(data$generation) <= 20, ]
  data$generation <- as.factor(data$generation)
  #form <- as.formula("diversity ~ Algorithm")
  #if(metric=="NMI"){
  #  data$diversity <- data$diversity.2
  #}
  #data$diversity.2 <- NULL
  xlab <- "Generaci\U00F3n"
  if(exp.name == "lv4"){
    xlab <- "Fase de comunicaci\U00F3n"
    title <- "Diversidad en DMNSGA-II seg\U00FAn operador de comunicaci\U00F3n"
  }else if(exp.name == "lv1"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de poblaci\U00F3n inicial"
  }else if(exp.name == "lv2"){
    title <- "Diversidad en DNSGA-II seg\U00FAn operador de selecci\U00F3n"
  }else if(exp.name == "lv3"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de cruzamiento y mutaci\U00F3n"
  }
  
  stat.test <- compare_means(
    diversity ~ Algorithm, data = data, group.by = "generation",
    method = "t.test", ref.group = "Crowding_Distance"
  )
  print(stat.test)
  
  bp <- ggline(data, x = "generation", y = "diversity",
                  fill = "Algorithm", palette = "jco",
                  add = "mean_sd", add.params = list(group = "Algorithm"),
                  position = position_dodge(0.2)) + 
        labs(title=title)
  bp + stat_pvalue_manual(
    stat.test, x = "generation", y.position = 33, step.increase=-0.2,
    label = "p.signif",
    position = position_dodge(0.3)
  )
  print(bp)
  ggsave(file.path(output.folder, paste0("diversity_evolution.png")))
  print("plot finished!")
  return()
  
  if(exp.name != "lv4"){
    pd <- position_dodge(0.3)
    ggplot(data.se, aes(x=generation, y=diversity, colour=Algorithm, group=Algorithm)) +
      geom_line(position=pd) +
      geom_point(position=pd) +
      geom_errorbar(aes(ymin=diversity-se, ymax=diversity+se), width=.1, position=pd) +
      labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", title=title) +
      #subtitle = get_test_label(stat.test, detailed = TRUE),
      #caption = get_pwc_label(stat.test)) +
      stat_pvalue_manual(pvalues, label = "{p.adj.signif}", hide.ns = TRUE, xmin=NULL, xmax=NULL, x=1:20, y.position = 0.15, step.increase = -0.2) +
      theme_minimal()
  }else{
    ggplot(data, aes(x=generation, y=diversity)) +
      geom_boxplot(aes(fill=Algorithm)) +
      labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", 
           title=title,
           #subtitle = get_test_label(stat.test, detailed = TRUE), 
           caption = get_pwc_label(stat.test)) +
      theme_minimal()
  }
  print(output.folder)
  ggsave(file.path(output.folder, paste0("diversity_evolution.png")))
  print("plot finished!")
}

plot_diversity_evolution <- function(data, exp.name, output.folder){
  library(ggplot2)
  #install.packages("ggtext")
  library(ggtext)
  data <- data[as.numeric(data$generation) <= 20, ]
  last.gen <- max(as.numeric(data$generation))
  #data$generation <- as.factor(data$generation)
  #form <- as.formula("diversity ~ Algorithm")
  #if(metric=="NMI"){
  #  data$diversity <- data$diversity.2
  #}
  #data$diversity.2 <- NULL
  xlab <- "Generaci\U00F3n"
  if(exp.name == "lv1"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de poblaci\U00F3n inicial"
    title_2 <- "Calidad en NSGA-II seg\U00FAn operador de poblaci\U00F3n inicial"
    labels <- c("Random", "Densidad", "Densidad diverso", "Densidad diverso modif.")
    data$Algorithm <- factor(data$Algorithm, levels = c("Random", "Selective", "Selective_Diverse", "Selective_Diverse+"))
  }else if(exp.name == "lv2"){
    #print(head(data))
    title <- "Diversidad en NSGA-II seg\U00FAn operador de selecci\U00F3n"
    labels <- c("Crowding distance", "Jaccard", "NMI")
    title_2 <- "Calidad en NSGA-II seg\U00FAn operador de selecci\U00F3n"
    data$Algorithm <- factor(data$Algorithm, levels = c("Crowding_Distance", "Jaccard", "NMI"))
  }else if(exp.name == "lv3"){
    title <- "Diversidad en NSGA-II seg\U00FAn operador de cruzamiento y mutaci\U00F3n"
    title_2 <- "Calidad en NSGA-II seg\U00FAn operador de cruzamiento y mutaci\U00F3n"
    labels <- c("Uniforme", "Crossover diverso", "Mutaci\U00F3n diversa", "Combinado")
    data$Algorithm <- factor(data$Algorithm, levels = c("Direct_Crossover", "Selective_Crossover", "Diverse_Mutation", "Combined"))
  }else if(exp.name == "lv4"){
    xlab <- "Fase de comunicaci\U00F3n"
    title <- "Diversidad en DMNSGA-II seg\U00FAn operador de comunicaci\U00F3n"
    title_2 <- "Calidad en DMNSGA-II seg\U00FAn operador de comunicaci\U00F3n"
    labels <- c("Sincr. elitista (intra-agente)", "Sincr. elitista (inter-agente)", "Sincr. diversa (intra-agente)", "Sincr. diversa (inter-agente)")
    data$Algorithm <- factor(data$Algorithm, levels = c("Elitist_sync_within","Elitist_sync_between", "Diverse_sync_within", "Diverse_sync_between"))
    data.sil <- data[!(data$Algorithm %in% c("Elitist_sync_between", "Diverse_sync_between")), ]
  }
  if(exp.name == "lv2"){
    data.2 <-  data[data$generation <= 10, ]
    data.2$generation <- as.factor(data.2$generation)
    data.2 <- gather(data.2, Resultado, Soluciones, survived:discarded, factor_key=TRUE)
    data.2$Algorithm <- factor(data.2$Resultado, levels=c("survived", "discarded"))
    ggplot(data.2, aes(x=generation, y=Soluciones, fill=Resultado, group=Resultado)) +
      geom_bar(stat="summary", fun="median", position="stack") +
      #geom_errorbar(aes(ymin=Soluciones-se, ymax=Soluciones+se), width=.1, position=pd) +
      labs(x=xlab, y="N\U00FAmero de soluciones en selecci\U00E9n de supervivencia", title=title) +
      scale_fill_discrete(labels=c("Sobreviviente", "Descartado")) +
      facet_wrap(~Algorithm)
    theme_minimal()
    ggsave(file.path(output.folder, "survival_selection.png"))
  }
  
  if(exp.name == "lv4"){
    data.sil$Algorithm <- factor(data.sil$Algorithm, levels = c("Elitist_sync_within","Diverse_sync_within"))
    data.se.sil <- summarySE(data.sil, measurevar="silhouette", groupvars=c("Algorithm","generation"))
    labels.sil <- c("Sincronizaci\U00F3n elitista", "Sincronizaci\U00F3n diversa")
  }else{
    data.sil <- data
    data.se.sil <- summarySE(data, measurevar="silhouette", groupvars=c("Algorithm","generation"))
    labels.sil <- labels
  }
  pvalues <- pairwise_t_test(data.sil[data.sil$generation == last.gen, ], silhouette ~ Algorithm, paired = TRUE, p.adjust.method = "bonferroni")
  print("P-VALUES QUALITY:")
  print(pvalues)
  labels_test <- paste0(pvalues$group1, " - ", pvalues$group2, " (P-value: ", pvalues$p.adj, ")\n")#, fill=is.signif.color)
  print(labels_test)
  
  pd <- position_dodge(0.3)
  ggplot(data.se.sil, aes(x=generation, y=silhouette, colour=Algorithm, group=Algorithm)) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=silhouette-ci, ymax=silhouette+ci), width=.1, position=pd) +
    labs(x=xlab, y="Promedio coeficiente silueta", 
         title=title_2, 
         colour="Operador") +
    scale_colour_discrete(labels=labels.sil) +
    #subtitle = get_test_label(res, detailed = TRUE), 
    #caption = get_pwc_label(pwc)) +
    theme_minimal()
  ggsave(file.path(output.folder, paste0("quality_evolution.png")))
  
  #data.se$generation <- as.integer(data.se$generation)
  #data.se$Algorithm <- as.factor(data.se$Algorithm)
  #stat.test <- data %>%
  #  group_by(generation) %>%
  #  #do(tidy(pairwise_t_test(diversity ~ Algorithm, data = ., paired = TRUE, p.adjust.method = "bonferroni"))) %>%
  #  #print()
  #  pairwise_t_test(
  #    diversity ~ Algorithm, paired = TRUE, 
  #    p.adjust.method = "bonferroni"
  #  ) %>%
  #  select(-df, -statistic, -p) # Remove details
  #stat.test <- stat.test %>% add_xy_position(x = "Algorithm")
  #print("C")
  #stat.test <- c()
  #n.alg <- length(unique(data$Algorithm))
  #pwc <- data.frame("generation"=1:20, "p.adj.signif"=rep(NA,20))
  #pvalues <- pairwise_t_test(data[data$generation == 1, ], diversity ~ Algorithm, paired = TRUE, p.adjust.method = "bonferroni")
  #pvalues <- pvalues[ which.min(pvalues$p.adj), ]
  #pvalues$generation  <- rep(20, 3)
  #print(pvalues, n=40)
  #for(i in 2:max(data$generation)){
  #  data.gen <- data[data$generation == i, ]
  #  #print(data.gen)
  #  pvalue <- pairwise_t_test(data.gen, diversity ~ Algorithm, paired = TRUE, p.adjust.method = "bonferroni")
  #  max.signif <- pvalue[ which.min(pvalue$p.adj), ]
  #  #pvalue <- pvalue %>% select(p.adj, p.adj.signif)
  #  pvalues <- rbind(pvalues, max.signif)
  #}
  #pvalues <- pvalues %>% add_xy_position(x = Algorithm)
  #pvalues <- as.data.frame(pvalues)
  #pvalues$y.position <- rep(50, 20)
  data.se <- summarySE(data, measurevar="diversity", groupvars=c("Algorithm","generation"))
  last.gen <- max(as.numeric(data$generation))
  pvalues <- pairwise_t_test(data[data$generation == last.gen, ], diversity ~ Algorithm, paired = TRUE, p.adjust.method = "bonferroni")
  print("P-VALUES DIVERSITY:")
  print(pvalues)
  labels_test <- paste0(pvalues$group1, " - ", pvalues$group2, " (P-value: ", pvalues$p.adj, ")\n")#, fill=is.signif.color)
  print(labels_test)
  #box.df <- data.frame( label = labels_test, x = 20, y = 1, hjust = 0,
  #  vjust = 1, orientation = "upright", color = "black", fill = "cornsilk")
  if(exp.name != "lv4"){
    pd <- position_dodge(0.3)
    #print(data.se)
    ggplot(data.se, aes(x=generation, y=diversity, colour=Algorithm, group=Algorithm)) +
      geom_line(position=pd) +
      geom_point(position=pd) +
      geom_errorbar(aes(ymin=diversity-ci, ymax=diversity+ci), width=.1, position=pd) +
      labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", colour="Operador", 
           title=title
           #get_pwc_label(pvalues),
           #caption = labels_test#paste("Generaci\U00E9n 20: ", get_test_label(pvalues, detailed = TRUE,  type = "text"))
           ) +
      scale_colour_discrete(labels=labels) +
      #geom_textbox(data=box.df, aes(
      #  x=x, y=y, label = label, color = color, fill = fill,
      #  hjust = hjust, vjust = vjust, orientation = orientation))+#, width = unit(0.4, "npc")) +
      #get_pwc_label(stat.test)) +
      #stat_pvalue_manual(pvalues, label = "{p.adj.signif}", hide.ns = TRUE, y.position = 0.15, step.increase = -0.2) +
      theme_minimal()
  }else{
    data$generation <- as.factor(data$generation)
    ggplot(data, aes(x=generation, y=diversity)) +
      geom_boxplot( aes(fill=Algorithm)) +
      labs(x=xlab, y="Distancia promedio entre soluciones (Jaccard)", 
           title=title, fill="Operador") +
           #subtitle = get_test_label(stat.test, detailed = TRUE), 
           #caption = get_pwc_label(stat.test)) +
      scale_fill_discrete(labels=labels) +
      theme_minimal()
  }
  ggsave(file.path(output.folder, paste0("diversity_evolution.png")))
  print("plot finished!")
}

test_operator_lv1_v2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Random", "Selective", "Selective_Diverse", "Selective_Diverse+")
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=8))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation", "silhouette", "survived", "discarded")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  seeds <- c()
  count <- 1
  for(j in 1:length(operators)){
    op <- operators[j]
    output.path <- file.path(output.folder, op, dataset)
    for(i in 1:runs){
      params$auto_adjust_initial_params <- FALSE
      print(paste0("operator ", op, " dataset ", dataset, " run ", i))
      output.exp <- file.path(output.path, i)
      seed <- as.numeric(Sys.time())
      seeds[count] <- seed
      if(dir.exists(output.exp)){
        next
      }
      if(op == "Random"){
        P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
      }else if(op == "Selective"){
        P  <- generate_diverse_initial_pop(distances, params, diverse_population=FALSE, seed=seed)
      }else if(op == "Selective_Diverse"){
        params$auto_adjust_initial_params <- FALSE
        P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE, seed=seed)
      }else if(op == "Selective_Diverse+"){
        params$auto_adjust_initial_params <- TRUE
        P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE, seed=seed)
      }
      
      
      nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder)
      count <- count + 1
      #evaluate_solutions(P, P.clustering.groups, distances, params$K, 
      #                   params$objDim, params$obj_maximize, dirname(output.exp), 
      #                   basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    } 
  }
  write.table(as.data.frame(seeds), file=file.path(output.folder, "seeds.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
}

test_operator_lv1 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Random", "Selective", "Selective_Diverse", "Selective_Diverse+")
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=8))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation", "silhouette", "survived", "discarded")
    write.table(res, file=file.path(output.folder, "diversity_evolution.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
  }
  seeds <- c()
  count <- 1
  for(i in 1:runs){
    params$auto_adjust_initial_params <- FALSE
    print(paste0("dataset ", dataset, " run ", i))
    seed <- as.numeric(Sys.time())
    seeds[count:(count+4)] <- seed
  
    output.exp <- file.path(output.folder, "Random", dataset, i)
    if(dir.exists(output.exp)){ next }
    P <- generate_initial_pop(params$popSize, params$K, distances$n.genes, seed) 
    nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name="Random", base.path=output.folder)
    
    output.exp <- file.path(output.folder, "Selective", dataset, i)
    if(dir.exists(output.exp)){ next }
    params$auto_adjust_initial_params <- FALSE
    P  <- generate_diverse_initial_pop(distances, params, diverse_population=FALSE, seed=seed)
    nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name="Selective", base.path=output.folder)
    
    output.exp <- file.path(output.folder, "Selective_Diverse", dataset, i)
    if(dir.exists(output.exp)){ next }
    #params$auto_adjust_initial_params <- FALSE
    P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE, seed=seed)
    nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name="Selective_Diverse", base.path=output.folder)
    
    output.exp <- file.path(output.folder, "Selective_Diverse+", dataset, i)
    if(dir.exists(output.exp)){ next }
    params$auto_adjust_initial_params <- TRUE
    P <- generate_diverse_initial_pop(distances, params, diverse_population=TRUE, seed=seed)
    nsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name="Selective_Diverse+", base.path=output.folder)
    
    count <- count + 4
  } 
  write.table(as.data.frame(seeds), file=file.path(output.folder, "seeds.csv"), append=FALSE, sep=",", row.names = FALSE, col.names = TRUE)
}



test_operator_lv2 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Crowding_Distance", "Jaccard", "NMI")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=8))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation", "silhouette", "survived", "discarded")
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
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="Jaccard"){
        params$diversity_level <- 2
        params$diversity_metric <- "jaccard"
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }else if(op=="NMI"){
        params$diversity_level <- 2
        params$diversity_metric <- "NMI"
        dnsga2(distances, params, output.exp, initial_population=P, debug=FALSE, plot=FALSE, calculate.diversity=TRUE, experiment.name=op, base.path=output.folder) 
      }
      
      #evaluate_solutions(P_next_generation, P.clustering.groups, distances, params$K, 
      #                   params$objDim, params$obj_maximize, dirname(output.exp), 
      #                   basename(output.exp), op, dataset, pareto.only=FALSE, plot=FALSE)
    }
  } 
}

test_operator_lv3 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Direct_Crossover", "Selective_Crossover", "Diverse_Mutation", "Combined")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=8))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation", "silhouette", "survived", "discarded")
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

test_operator_lv4 <- function(params, distances, dataset, output.folder, runs=31){
  operators <- c("Elitist_sync", "Diverse_sync")
  
  if(!file.exists(file.path(output.folder, "diversity_evolution.csv"))){
    res <- as.data.frame(matrix(nrow=1, ncol=8))
    colnames(res) <- c("id", "Algorithm", "Dataset", "diversity", "generation", "silhouette", "survived", "discarded")
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

summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

compare_operators()
