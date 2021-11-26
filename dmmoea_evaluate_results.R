#!bin/usr/env Rstudio
evaluate_results <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 3){
    print(paste0("Not enough parameters (", argnum, "/3)"))
    return(-1)
  }
  path <- args[1]
  results.path <- args[2]
  test <- args[3]
  
  setwd(path)
  source("dmmoea_functions.R")
  #source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  #source("dmmoea_irace_conf.R")

  results.path <- file.path(results.path)
  if(!file.exists(file.path(results.path, "limits.csv"))){
    #warning("Limits not found, please run \"dmmoea_normalize_limits.R\" first.")
    #return(-1)
    print("Limits not found. Getting normalization limits...")
    get_normalization_limits(results.path)
    print("Finished.")
  }
  
  evaluate_run_results(results.path)
  if(!file.exists(file.path(results.path, "plot_data.csv"))){
    warning("Quality plot data not found!. Aborting...")
    return(-1)
  }
  if(!file.exists(file.path(results.path, "plot_data_diversity.csv"))){
    warning("Diversity plot data not found!. Aborting...")
    return(-1)
  }
  plot.data <- read.table(file.path(results.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data.diversity <- read.table(file.path(results.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  plot_algorithm_comparison(results.path, plot.data)
  plot_algorithm_comparison_diversity(results.path, plot.data.diversity)
  plot_algorithm_comparison_pareto(results.path)
}

evaluate_run_results <- function(path, maximize=FALSE, alpha=0.5){
  limits <- read.csv(file.path(path, "limits.csv"), header = TRUE)
  scaler.f1 <- function(x){ (x-limits$min.f1)/(limits$max.f1-limits$min.f1) }
  scaler.f2 <- function(x){ (x-limits$min.f2)/(limits$max.f2-limits$min.f2) }
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  plot.data <- as.data.frame(matrix(nrow=0, ncol=6))
  plot.data.diveristy <- as.data.frame(matrix(nrow=0, ncol=6))
  colnames(plot.data) <- c("id", "Algorithm", "Dataset", "Hypervolume", "Silhouette", "Delta")
  colnames(plot.data.diversity) <- c("id", "Algorithm", "Dataset", "Metric", "Diversity", "Cluster_Ratio")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    if(algorithm == "figures"){ next }
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      distances <- load.gene.distance(dataset, alpha=alpha)
      exp.path <- file.path(path, algorithm, dataset)
      evaluation.file <- read.table(file.path(exp.path, "evaluations.csv"), header=TRUE, sep=",", row.names=NULL)
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        data.pareto <- read.table(file.path(exp.path, experiment, "population.csv"), sep=",", header=FALSE, row.names=NULL)
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv <- calculate_hypervolume(data, c(1,1), maximize)
        sil <- evaluation.file[k, "avg_sil"]
        delta <- evaluation.file[k, "delta"]
        diversity.jaccard <- diversity_analysis(data.pareto, distances, metric="jaccard", 
                                                path=file.path(exp.path, experiment), alpha=alpha, plot=TRUE)
        diversity.NMI <- diversity_analysis(data.pareto, distances, metric="NMI", 
                                            path=file.path(exp.path, experiment), alpha=alpha, plot=TRUE)
        values.diversity <- data.frame("id"=rep(experiment, 2), "Algorithm"=rep(algorithm, 2), 
                                       "Dataset"=rep(dataset,2), "Metric"=c("jaccard", "NMI"), 
                                       "Diversity"=c(diversity.jaccard$avg.dist, diversity.NMI$avg.dist),
                                       "Cluster_Ratio"=c(diversity.jaccard$k.best, diversity.NMI$k.best))
        values <- data.frame("id"=experiment, "Algorithm"=algorithm, "Dataset"=dataset,
                             "Hypervolume"=hv, "Silhouette"=sil, "Delta"=delta)
        plot.data <- rbind(plot.data, values)
        plot.data.diversity <- rbind(plot.data.diversity, values.diversity)
      }
    }
  }
  write.table(plot.data, file=file.path(path, "plot_data.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(plot.data.diveristy, file=file.path(path, "plot_data_diversity.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
}

evaluate_results()
