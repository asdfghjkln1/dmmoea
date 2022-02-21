#!bin/usr/env Rstudio
evaluate_results <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum < 2){
    print(paste0("Not enough parameters (", argnum, "/3)"))
    return(-1)
  }
  path <- args[1]
  results.path <- args[2]
  limit.run <- as.numeric(args[3])
  
  setwd(path)
  source("dmmoea_functions.R")
  #source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  #source("dmmoea_irace_conf.R")

  results.path <- file.path(path, results.path)
  print("Path in:")
  print(results.path)
  #if(!file.exists(file.path(results.path, "limits.csv"))){
    #warning("Limits not found, please run \"dmmoea_normalize_limits.R\" first.")
    #return(-1)
  if(length(list.dirs(results.path, recursive = FALSE)) > 1){
    print("Getting normalization limits...")
    limits <- get_normalization_limits(results.path) 
  }
  #}
  #limits <- read.csv(file.path(results.path, "limits.csv"), header = TRUE)
  print("Finished.")
  
  if(!file.exists(file.path(results.path, "plot_data.csv"))){
    print("Quality plot data not found!. Initiating evaluation...")
    evaluate_run_results(results.path, limits, limit.run = limit.run)
  }else if(!file.exists(file.path(results.path, "plot_data_diversity.csv"))){
    print("Diversity plot data not found!. Initiating evaluation...")
    evaluate_run_results(results.path, limits, limit.run = limit.run)
  }
  if(file.exists(file.path(results.path, "data_pareto.csv"))){
    print("Experiment evaluation done. Plotting results...")
    plot_algorithm_comparison_pareto(results.path, load.data=TRUE, limit.run = limit.run)
  }else{
    print("Getting pareto front and plotting...")
    plot_algorithm_comparison_pareto(results.path, load.data=FALSE, limit.run = limit.run)
  }
  
  plot.data <- read.table(file.path(results.path, "plot_data.csv"), sep=",", header=TRUE, row.names=NULL)
  plot.data.diversity <- read.table(file.path(results.path, "plot_data_diversity.csv"), sep=",", header=TRUE, row.names=NULL)
  print("Pareto comparison done. Plotting other comparison metrics...")
  plot_algorithm_comparison(results.path, plot.data)
  plot_algorithm_comparison_diversity(results.path, plot.data.diversity)
  print("Done.")
}

evaluate_run_results <- function(path, limits, maximize=FALSE, alpha=0.5, limit.run=Inf){
  algorithms <- list.dirs(path=path, full.names = FALSE, recursive = FALSE)
  algorithms <- algorithms[!(algorithms %in% "figures")]
  plot.data <- as.data.frame(matrix(nrow=0, ncol=6))
  plot.data.diversity <- as.data.frame(matrix(nrow=0, ncol=6))
  colnames(plot.data) <- c("id", "Algorithm", "Dataset", "Hypervolume", "Silhouette", "Delta")
  colnames(plot.data.diversity) <- c("id", "Algorithm", "Dataset", "Metric", "Diversity", "Cluster_Ratio")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      distances <- load.gene.distance(dataset, alpha=alpha)
      exp.path <- file.path(path, algorithm, dataset)
      limits <- read.table(file.path(exp.path, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
      scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
      scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
      evaluation.file <- read.table(file.path(exp.path, "evaluations.csv"), header=TRUE, sep=",", row.names=NULL)
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      experiments <- experiments[as.numeric(experiments) < limit.run]
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        data.pareto <- read.table(file.path(exp.path, experiment, "population.csv"), sep=",", header=TRUE, row.names=NULL)
        data <- data.frame("f1"=scaler.f1(data[, 1]), "f2"=scaler.f2(data[, 2]))
        hv <- calculate_hypervolume(data, c(1,1), maximize) #eaf::hypervolume(data, c(1,1), maximize=FALSE) #
        sil <- evaluation.file[k, "avg_sil"]
        delta <- evaluation.file[k, "delta"]
        diversity.jaccard <- diversity_analysis(data.pareto, distances, metric="jaccard", 
                                                exp.path=file.path(exp.path, experiment), alpha=alpha, plot=FALSE)
        diversity.NMI <- diversity_analysis(data.pareto, distances, metric="NMI", 
                                            exp.path=file.path(exp.path, experiment), alpha=alpha, plot=FALSE)
        values <- data.frame("id"=experiment, "Algorithm"=algorithm, "Dataset"=dataset,
                             "Hypervolume"=hv, "Silhouette"=sil, "Delta"=delta)
        plot.data <- rbind(plot.data, values)
        if(!is.na(diversity.jaccard)){
          values.diversity <- data.frame("id"=rep(experiment, 2), "Algorithm"=rep(algorithm, 2), 
                                         "Dataset"=rep(dataset,2), "Metric"=c("jaccard", "NMI"), 
                                         "Diversity"=c(diversity.jaccard$avg.dist, diversity.NMI$avg.dist),
                                         "Cluster_Ratio"=c(diversity.jaccard$k.ratio, diversity.NMI$k.ratio))
          plot.data.diversity <- rbind(plot.data.diversity, values.diversity)
        }
        print(paste(algorithm, dataset, experiment, "evaluated!"))
      }
    }
  }
  write.table(plot.data, file=file.path(path, "plot_data.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(plot.data.diversity, file=file.path(path, "plot_data_diversity.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
}

diversity_analysis <- function(P, distances, metric, exp.path=NULL, alpha=0.5, plot=FALSE){
  P <- as.matrix(P)
  if(nrow(P) < 2){
    return(list("diss"=NA, "avg.dist"=NA, "k.ratio"=NA))
  }
  res.P <- cluster_data(distances, P, alpha)
  P <- res.P$population
  P.groups <- res.P$clustering.results 
  d.matrix.P <- calculate_diversity_matrix(P.groups, metric)
  avg.dist <- mean(d.matrix.P[lower.tri(d.matrix.P, diag = FALSE)])
  d.matrix.P <- as.dist(d.matrix.P)
  if(nrow(P) < 4){
    return(list("diss"=d.matrix.P, "avg.dist"=avg.dist, "k.ratio"=NA))
  }
  
  
  hc <- hclust(d=d.matrix.P, method="single")
  sil <- c(-1)
  for(i in 2:(nrow(P)-1)){
    sil[i] <- mean(silhouette(cutree(hc, k=i), dist=d.matrix.P)[,"sil_width"])
  } 
  best.k <- which.max(sil)
  k.ratio <- best.k/(nrow(P)-1)
  #print("Silouettes:")
  #print(sil)
  #min.k <- 2
  #max.k <- max(nrow(P) - 2, min.k + 2) #max(round(sqrt(nrow(P)), 1), 4)
  #print(paste0("(", min.k, " - " , max.k, ")"))
  #if(nrow(P) < 5){
  #  #warning("This pareto front has too few solutions, diversity may not be accurate!!")
  #  return(list("diss"=d.matrix.P, "avg.dist"=avg.dist, "k.ratio"=mean(c(max.k, min.k))/nrow(P)))
  #}
  
  #indexes <- c("frey", "dunn", "cindex", "silhouette", "mcclain")
  #values <- list()
  #for(i in 1:length(indexes)){
  #  best <- NbClust::NbClust(diss = d.matrix.P, distance = NULL, index=indexes[i],
  #                           min.nc = min.k, max.nc = max.k, method = "single")
  #  values[i] <- as.integer(best$Best.nc["Number_clusters"])
  #}
  #getmode <- function(v) {
  #  uniqv <- unique(v)
  #  uniqv[which.max(tabulate(match(v, uniqv)))]
  #}
  #best.k <- max(getmode(unlist(values)), min.k + 2)
  if(plot){# && !dir.exists(file.path(exp.path, "diversity"))){
    pam.res <- cluster::pam(x=d.matrix.P, diss=FALSE, k = best.k, do.swap = FALSE)
    plot.pam <- factoextra::fviz_cluster(pam.res, geom = "point", main=paste0("Pareto similarity clustering (", metric, ")"))
    dir.create(file.path(exp.path, "diversity"), recursive = TRUE, showWarnings = FALSE)
    out.file <- file.path(exp.path, "diversity", paste0("diversity_", metric, ".png"))
    png(out.file)
    print(plot.pam)
    dev.off()
  }
  return(list("diss"=d.matrix.P, "avg.dist"=avg.dist, "k.ratio"=k.ratio))
}

evaluate_results()
