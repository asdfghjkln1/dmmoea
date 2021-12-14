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
  plot.data <- read.table(file.path(results.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data.diversity <- read.table(file.path(results.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data$Algorithm <- as.factor(plot.data$Algorithm)
  plot.data$id <- as.factor(plot.data$id)
  
  datasets <- unique(plot.data$Dataset)
  metrics <- c("Hypervolume", "Silhouette", "Delta")
  metrics.div <- c("Diversity", "Cluster_Ratio")
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
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
  w <- 4 + 0.5*(length(unique(data[, exp.group])) - 3)#ifelse(exp.group>3, 6,5)
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
    theme(strip.text.x = element_blank()) +
    stat_pvalue_manual(pwc, label = "p = {p.adj}", hide.ns = TRUE)
  ggsave(file.path(output.path, paste0(metric, "_results_", dataset.name, ".png")), width = w, height = 6)
  print(paste(metric, dataset.name, "... Done."))
}

compare_algorithms()
