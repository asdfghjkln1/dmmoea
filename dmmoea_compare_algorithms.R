#!bin/usr/env Rstudio
compare_algorithms <- function(){
  
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 2){
    print(paste0("Not enough parameters (", argnum, "/2)"))
    return(-1)
  }
  path <- args[1] # "X:\\Universidad\\dmmoea" #args[1] #
  results.path <- args[2] #"Tests\\runs\\XieBeni" # args[2]
  
  setwd(path)
  library(ggpubr)
  library(rstatix)
  library(ggplot2)
  library(gridExtra)
  source("dmmoea_functions.R")

  results.path <- file.path(path, results.path)
  figure.path <- file.path(results.path, "figures", "comparison")
  print("Path in:")
  print(results.path)
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
  
  datasets <- unique(plot.data$Dataset)
  metrics <- c("Hypervolume", "Silhouette", "Delta")
  metrics.div <- c("Diversity", "Cluster_Ratio")
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    data.pareto <- plot.pareto.norm[ plot.pareto.norm$Dataset == dataset, ]
    compare_pareto_front(data.pareto, dataset, figure.path)
    
    data <- plot.data[plot.data$Dataset == dataset, ]
    data.diversity <- plot.data.diversity[plot.data.diversity$Dataset == dataset, ]
    data.diversity.NMI <- data.diversity[data.diversity$Metric == "NMI", ]
    data.diversity.jaccard <- data.diversity[data.diversity$Metric == "jaccard", ]
    #formula <- as.formula(as.symbol(metric) ~ as.symbol("Algorithm") | as.symbol("Dataset"))
    lapply(metrics, function(x) kruskal.multi.variable.tests(data, x, "Algorithm", dataset, figure.path))
    lapply(metrics.div, function(x) kruskal.multi.variable.tests(data.diversity.NMI, x, "Algorithm", dataset, figure.path))
    lapply(metrics.div, function(x) kruskal.multi.variable.tests(data.diversity.jaccard, x, "Algorithm", dataset, figure.path))
    #res <- friedman.test.with.post.hoc(Hypervolume ~ Algorithm | id, data=data, to.plot.parallel = F, dataset=dataset)
    #print(res)
  }
}

kruskal.multi.variable.tests <- function(data, metric, exp.group, dataset.name, output.path){
  dir.create(output.path, recursive=TRUE, showWarnings = FALSE)
  form <- as.formula(call("~", as.symbol(metric), as.symbol(exp.group)))
  kruskal.res <- kruskal_test(data, formula=form)
  pwc <- wilcox_test(data, formula=form, p.adjust.method="bonferroni")
  pwc <- pwc %>% add_xy_position(x = exp.group)
  w <- 4 + 0.6*(length(unique(data[, exp.group])) - 3)#ifelse(exp.group>3, 6,5)
  Y <- pwc$y.position
  gap.data <- max(data[, metric]) - min(data[, metric])
  gap <- max(Y[2] - Y[1], gap.data*0.07)
  pwc$y.position <- Y - gap*seq(from=1, to=0, length.out=length(unique(data[, exp.group])))
  
  ggplot(data, aes_string(x=exp.group, y=metric)) +
    geom_boxplot(aes_string(fill=exp.group)) +
    labs(subtitle = get_test_label(kruskal.res, detailed = TRUE), 
         caption = get_pwc_label(pwc),
         title=paste0(metric, " comparison: ", dataset.name, " dataset")) +
    theme_minimal() +
    theme(strip.text.x = element_blank(), axis.text.x = element_text(angle=25)) +
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


compare_algorithms()
