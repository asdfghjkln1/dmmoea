spaces_comparison() <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum < 7){
    print(paste0("Not enough parameters (", argnum, "/7)"))
    return(-1)
  }
  path <- args[1] # "X:\\Universidad\\dmmoea"
  params.path <- args[2] #"Tests\\a"
  algorithms.param <- args[3] # "tmix,mfuzz,dnsga2"
  ref.algorithm <- args[4] # "dnsga2"
  evaluations <- as.numeric(args[5]) # 2000
  trial.start <- as.numeric(args[6]) # 1
  trial.stop <- as.numeric(args[7]) # 31

  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  source("moc.gapbk/R/libraries.R")
  source("moc.gapbk/R/main.R")

  tune.path <- file.path(path, params.path)
  
  #Load best params
  best_params <- read.table(file.path(tune.path, ref.algorithm, "best_configurations.csv"), sep=",", header=TRUE, row.names=NULL)
  
  test.path <- file.path(path, "Tests", "experiments", paste0(best_params$objectives, "_2"))
  # Initialize params
  params <- init_parameters(objectives=best_params$objectives)
  params$K <- as.numeric(best_params$K)
  params$objectives <- best_params$objectives
  params$evaluations <- evaluations
  params$popSize <- as.numeric(best_params$popSize)
  params$mating_rate <- best_params$mating_rate
  params$mutation_rate <- best_params$mutation_rate
  params$alpha <- best_params$alpha
  params$is_random_population <- as.numeric(best_params$is_random_population)
  params$auto_adjust_initial_params <- best_params$auto_adjust_initial_params
  if(!is.na(params$auto_adjust_initial_params)){
    params$min_density_radius <- best_params$min_density_radius
    params$max_density_radius <- best_params$max_density_radius
    params$density_tol <- best_params$density_tol
  }
  params$diversity_metric <- best_params$diversity_metric
  params$diversity_level <- best_params$diversity_level
  params$phases <- best_params$phases
  params$agents <- best_params$agents
  params$sync_off <- ifelse(is.na(as.numeric(best_params$sync_off)), 0, as.numeric(best_params$sync_off))
  params$convergence_tol <- best_params$convergence_tol
  params$mutation_radius <- best_params$mutation_radius
  params$seed <- runif(1, 0, 1)*1235
  
  algorithms <- strsplit(algorithms.param, ",")[[1]]
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    print("Starting algorithm:")
    print(algorithm)
    if(algorithm == "figures"){ next }
    datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      print("Starting dataset:")
      print(dataset)
      output.folder <- file.path(test.path, algorithm, dataset)
      limits <- read.table(file.path(tune.path, ref.algorithm, dataset, "limits.csv"), sep=",", row.names=NULL, header=TRUE)
      execute_tests(params, path, output.folder, algorithm, dataset, limits, trial.start=trial.start, trial.stop = trial.stop) 
    }
  }
}

execute_tests <- function(params, path, output.folder, algorithm, dataset, limits, trial.start=1, trial.stop=31){
  distances <- load.gene.distance(dataset, params$alpha)
  output.exp <- file.path(output.folder, i)#file.path(basename(params$test.path), "Debug", "test")
  dir.create(output.folder, showWarnings=FALSE, recursive=TRUE)
  
  res <- run_spaces(distances, params, output.exp, limits, dataset)
  res2 <- run_cluster_n(distances, output.exp, dataset, trials=10)
}


run_cluster_n <- function(distances,output.exp, dataset, trials){
  best.nclusts <- data.frame()
  indexes <- c("frey", "mcclain", "cindex","silhouette","dunn")
  for(trial in 1:trials){
    nclusts <- integer(5)
    print(paste("Starting number of cluster analysis in dataset", dataset))
    require(NbClust)
    for(i in 1:length(indexes)){
      nb <- NbClust(data=NULL, diss = as.dist(distances$comp.dist), 
                  distance = NULL,
                  min.nc = 2, 
                  max.nc = 10, 
                  index = indexes[i],
                  method = "ward.D2")
      nclusts[i] <- as.numeric(nb$Best.nc)
    }
    #names(nclusts) <- indexes
    #print(nclusts)
    #print(median(nclusts))
    best.nclusts <- rbind(best.nclusts, as.data.frame(t(nclusts)))
  }
  colnames(best.nclusts) <- indexes
  print("NbClust results:")
  print(best.nclusts)
  best <- median(apply(best.nclusts, 2, median))
  print(paste0("Best cluster number for dataset ", dataset, ": ", best))
  return(best)
}

run_spaces <- function(distances, params, output.exp, limits, dataset){
  
  size <- 5000
  K <- params$K
  population <- generate_initial_pop(size, params$K, distances$n.genes) 
  P.data <- cluster_data(distances, population, params$alpha)
  population <- P.data$population
  clustering.groups <- P.data$clustering
  row.names(population) <- 1:nrow(population)
  P.rows <- row.names(population)
  population <- remove_duplicated(population, params$K)
  clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
  P.rows <- row.names(population)
  
  P <- evaluate_population(population, distances, clustering.groups, params)
  P.clustering.groups <- clustering.groups[match((row.names(P)), P.rows)]
  
  front <- P[P$rnkIndex == 1, ]
  print(front)
  clustering.front <- P.clustering.groups[as.numeric(row.names(front))]
  print(nrow(front))
  diversity.jaccard <- calculate_diversity_matrix(clustering.front, "jaccard")
  diversity.nmi <- calculate_diversity_matrix(clustering.front, "NMI")
  front.n <- nrow(front)
  dist.decision.jaccard <- rep(NA, front.n)
  obj.index <- (K+1):(params$objDim)
  dist.ratio.jaccard <- list()
  dist.ratio.NMI <- list()
  #comb <- 1
  #for(i in 1:(front.n-2)){
  #  for(j in (i+1):(front.n-1)){
  #    for(k in (j+1):front.n){
  #      decision.jaccard.a <- diversity.jaccard[i,j]
  #      decision.jaccard.b <- diversity.jaccard[i,k]
  #      decision.jaccard.c <- diversity.jaccard[j,k]
  #      decision.NMI.a <- diversity.NMI[i,j]
  #      decision.NMI.b <- diversity.NMI[i,k]
  #      decision.NMI.c <- diversity.NMI[j,k]
  #      dist.obj.a <- abs(front[i, obj.index] - front[j, obj.index])
  #      dist.obj.b <- abs(front[i, obj.index] - front[k, obj.index])
  #      dist.obj.c <- abs(front[j, obj.index] - front[k, obj.index])
  #      
  #      dist.ratio.jaccard[[comb]] <- c(dist.obj.a/decision.jaccard.a, dist.obj.b/decision.jaccard.b, dist.obj.c/decision.jaccard.c)
  #      dist.ratio.NMI[[comb]] <- c(dist.obj.a/decision.NMI.a, dist.obj.b/decision.NMI.b, dist.obj.c/decision.NMI.c)
  #      comb <- comb + 1
  #    }
  #  }
  #}
  if(FALSE){

  library(ggplot2)
  library(ggdendro)
  library(cowplot)
  
  hc.jaccard <- hclust(as.dist(diversity.jaccard))
  #png("hc_jaccard.png")
  
  plot.jaccard <- ggdendrogram(hc.jaccard, rotate = FALSE) +
                  labs(title=paste0("Soluciones pareto"), subtitle="Espacio de decisi\U00F3n (Jaccard)",
                       x="Id. Soluci\U00F3n", y="Distancia (Jaccard)")
  #plot.jaccard <- plot(hc.jaccard, 
  #                      main=paste0("Soluciones pareto \n Espacio de decisión (Jaccard)"),
  #                      xlab="Id. Solución", ylab="Distancia (Jaccard)")
  hc.nmi <- hclust(as.dist(diversity.nmi))
  
  plot.nmi <- ggdendrogram(hc.nmi, rotate = FALSE) +
    labs(title=paste0("Soluciones pareto"), subtitle="Espacio de decisi\U00F3n (NMI)",
         x="Id. Soluci\U00F3n", y="Distancia (NMI)")
  #plot.nmi <- plot(hc.nmi, 
  #                     main=paste0("Soluciones pareto \n Espacio de decisión (NMI)"),
  #                     xlab="Id. Solución", ylab="Distancia (NMI)")
  
  first.col <- plot_grid(plot.jaccard, plot.nmi, labels = c('A', 'B'), ncol = 1, align = 'v')
  
  ###
  pareto <- normalise_pareto(front[, c("f1", "f2")], limits)
  print(pareto)
  pareto <- cbind(pareto, 1:front.n)
  alpha <- ifelse(front[, "rnkIndex"] == 1, 1, 0.5)
  pareto <- cbind(pareto, alpha)
  colnames(pareto) <- c("f1", "f2", "id", "alpha")
  print(pareto)
  
  second.col <- ggplot(pareto, aes(x=f1, y=f2, label=id)) +
                     labs(title=paste0("Soluciones pareto"), subtitle="Espacio objetivo", x="Expresi\U00F3n g\U00E9nica (normalizado)", y="Funci\U00F3n biol\U00F3gica (normalizado)") +
                     geom_point(size=2) +
                     geom_text(hjust=-0.6, vjust=-0.9) +
                     geom_line(alpha=alpha) +
                     scale_alpha_continuous(guide=FALSE) +
                     #xlim(0, 1) + #max(pareto$f1)) +
                     #ylim(0, 1) + #max(pareto$f2))  
                     theme_minimal()
  
  plots <- plot_grid(first.col, second.col, labels = c('', 'C'), nrow = 1, rel_widths = c(1, 1.3))
  ggsave(paste0("spaces_", dataset,".png"), plots, height=5, width=8)
  #suppressWarnings(ggsave("spaces.png", plots, height=5, width=7))
  return(NULL)
  #ggdendrogram(hc.jaccard)
  #suppressWarnings(ggsave("hc_jaccard.png", height=5, width=7))
  
  ###
  
  
  
  
  
  suppressWarnings(ggsave("hc_nmi.png", height=5, width=7))
  
  return(NULL)
  
  print("Finished:")
  
  
  library("mstknnclust")
  library("igraph")
  
  g  <- graph_from_adjacency_matrix(as.matrix(diversity.jaccard), weighted=TRUE, mode="undirected")
  print(E(g)$weight)
  w <- E(g)$weight
  w.norm <- 9*(w-min(w))/(max(w)-min(w)) + 1
  print(w.norm)
  
  g_clust <- cluster_fast_greedy(g, weights = w.norm)
  #g_mst <- mst(g)
  #g_mst <- as.undirected(g)
  print(g_clust)
  print(g_clust$groups)
  png("mst_jaccard.png")
  plot(g, vertex.size=30, edge.width=w.norm,
       layout=igraph::layout.fruchterman.reingold(g, niter=10000),
       vertex.color=igraph::clusters(g)$membership, 
       main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisión (Jaccard)"))
  dev.off()
  
  
  g  <- graph_from_adjacency_matrix(as.matrix(diversity.nmi), weighted=TRUE)
  #print(E(g)$weight)
  g_mst <- mst(g)
  g_mst <- as.undirected(g_mst)
  w <- E(g_mst)$weight
  w.norm <- 9*(w-min(w))/(max(w)-min(w)) + 1
  print(w.norm)
  print(igraph::clusters(g_mst))
  png("mst_nmi.png")
  plot(g_mst, vertex.color=NA, vertex.size=30, edge.width=w.norm,
       layout=igraph::layout.fruchterman.reingold(g_mst, niter=10000),
       vertex.color=igraph::clusters(g_mst)$membership, 
       main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisión (NMI)"))
  dev.off()
  #print("c")

  
  return(NULL)
  }
  
  library(cowplot)
  library(ggplot2)
  library(ggdendro)
  
  hc.jaccard <- hclust(as.dist(diversity.jaccard))
  #png("hc_jaccard.png")
  
  plot.jaccard <- ggdendrogram(hc.jaccard, rotate = FALSE) +
    labs(title=paste0("Soluciones pareto"), subtitle="Espacio de decisi\U00F3n (Jaccard)",
         x="Id. Soluci\U00F3n", y="Distancia (Jaccard)")
  #plot.jaccard <- plot(hc.jaccard, 
  #                      main=paste0("Soluciones pareto \n Espacio de decisión (Jaccard)"),
  #                      xlab="Id. Solución", ylab="Distancia (Jaccard)")
  hc.nmi <- hclust(as.dist(diversity.nmi))
  
  plot.nmi <- ggdendrogram(hc.nmi, rotate = FALSE) +
    labs(title=paste0("Soluciones pareto"), subtitle="Espacio de decisi\U00F3n (NMI)",
         x="Id. Soluci\U00F3n", y="Distancia (NMI)")
  #plot.nmi <- plot(hc.nmi, 
  #                     main=paste0("Soluciones pareto \n Espacio de decisión (NMI)"),
  #                     xlab="Id. Solución", ylab="Distancia (NMI)")
  
  first.col <- plot_grid(plot.jaccard, plot.nmi, labels = c('A', 'B'), ncol = 1, align = 'v')

  knn.jaccard <- mstknnclust::mst.knn(diversity.jaccard)
  knn.nmi <- mstknnclust::mst.knn(diversity.nmi)
  
  print("Jaccard network info")
  print(knn.jaccard)
  print("NMI network info")
  print(knn.nmi)
  
  igraph::V(knn.jaccard$network)$label.cex <- seq(1.2,1.2,length.out=2)
  
  png(paste0("knn_jaccard_", dataset, ".png"))
  
  plot(knn.jaccard$network, vertex.size=30,
       vertex.color=igraph::clusters(knn.jaccard$network)$membership, 
       layout=igraph::layout.fruchterman.reingold(knn.jaccard$network, niter=10000),
       main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisi/U00F3n (Jaccard)"))
  
  #ggsave(paste0("jaccard_", dataset,".png"), height=5, width=8)
  
  #p1.jaccard <- recordPlot()
  dev.off()
  #plot.new()  
  
  
  igraph::V(knn.nmi$network)$label.cex <- seq(1.2,1.2,length.out=2)
  
  png(paste0("knn_nmi_", dataset, ".png"))
  
  plot(knn.nmi$network, vertex.size=30, 
       vertex.color=igraph::clusters(knn.nmi$network)$membership, 
       layout=igraph::layout.fruchterman.reingold(knn.nmi$network, niter=10000),
       main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisi/U00F3n (NMI)"))
  
  #ggsave(paste0("nmi_", dataset,".png"), height=5, width=8)
  
  #p1.nmi <- recordPlot()
  dev.off()
  #plot.new()  
  
  
  #first.col <- plot_grid("", "", labels = c('A', 'B'), ncol = 1, align = 'v')
  
  pareto <- normalise_pareto(front[, c("f1", "f2")], limits)
  print(pareto)
  pareto <- cbind(pareto, 1:front.n)
  alpha <- ifelse(front[, "rnkIndex"] == 1, 1, 0.5)
  pareto <- cbind(pareto, alpha)
  colnames(pareto) <- c("f1", "f2", "id", "alpha")
  print(pareto)
  
  second.col <- ggplot(pareto, aes(x=f1, y=f2, label=id)) +
    labs(title=paste0("Soluciones pareto"), subtitle="Espacio objetivo", x="Expresi\U00F3n g\U00E9nica (normalizado)", y="Funci\U00F3n biol\U00F3gica (normalizado)") +
    geom_point(size=2) +
    geom_text(hjust=-0.6, vjust=-0.9) +
    geom_line(alpha=alpha) +
    scale_alpha_continuous(guide=FALSE) +
    #xlim(0, 1) + #max(pareto$f1)) +
    #ylim(0, 1) + #max(pareto$f2))  
    theme_minimal()
  
  plots <- plot_grid(first.col, second.col, labels = c('', 'C'), nrow = 1, rel_widths = c(1, 1.3))
  ggsave(paste0("spaces_", dataset,".png"), plots, height=5, width=8)
  
  
  #print("Jaccard triangle semblance:")
  #print(dist.ratio.jaccard)
  #for(i in 1:comb){
  #  print(paste0("Pair ", i, ":  "))
  #  round(dist.ratio.jaccard[[i]], 3)
  #}
  #print("NMI triangle semblance:")
  #print(dist.ratio.NMI)
  #for(i in 1:comb){
  #  print(paste0("Pair ", i, ":  "))
  #  round(dist.ratio.NMI[[i]], 3)
  #}
  #front.2 <- P[P$rnkIndex == 2, ]
  #print(front.2)
  
  return(list("population"=front, "clustering"=clustering.front))#, "jaccard"=dist.ratio.jaccard, "NMI"=dist.ratio.NMI))
}


get.medoid.diss.matrix <- function(groups, dist){
  k <- length(unique(groups))
  medoids <- as.data.frame(matrix(ncol=k, nrow=1))
  for(i in 1:k){
    elem.group <- dist[groups == i, groups == i]
    dist.sum <- apply(elem.group,2, sum)
    medoids[1, i] <- which.min(dist.sum)
  }
  return(medoids)
}

spaces_comparison()