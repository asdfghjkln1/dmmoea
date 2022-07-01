spaces_comparison <- function(){
  #args <- commandArgs(trailingOnly = TRUE)
  #argnum <- length(args)
  #if(argnum < 2){
  #  print(paste0("Not enough parameters (", argnum, "/2)"))
  #  return(-1)
  #}
  path <- "X:\\Universidad\\dmmoea" #args[1] #
  trials <- 1 # as.numeric(args[2])
  size <- 50
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  
  test.path <- file.path(path, "spaces")
  datasets <- c("arabidopsis", "cell_cycle", "serum", "sporulation")
  
  #for(j in 1:length(datasets)){
  #  dataset <- datasets[j]
  #  print("Starting dataset:")
  #  print(dataset)
  #res <- run_spaces(test.path, datasets, size=size, K=5, trials=trials)
  #spaces_mstknn_plot(test.path, comp.measure="VI")
  #spaces_dist_ratio_plot(test.path, plot="scaled")
  #spaces_dist_ratio_plot(test.path, plot="non_scaled")
  spaces_dist_ratio_plot(test.path, plot="percentile")
}




run_cluster_n <- function(distances, output.folder, dataset, trials){
  best.nclusts <- data.frame()
  indexes <- c("frey", "mcclain", "cindex","silhouette","dunn")
  nclusts <- integer(5)
  print(paste("Starting number of cluster analysis in dataset", dataset))
  require(NbClust)
  
  for(j in 1:trials){
    output.exp <- file.path(output.folder, j)
    for(i in 1:length(indexes)){
      nb <- NbClust(data=NULL, diss = as.dist(distances$comp.dist), 
                  distance = NULL,
                  min.nc = 2, 
                  max.nc = 10, 
                  index = indexes[i],
                  method = "ward.D2")
      nclusts[i] <- as.numeric(nb$Best.nc)
    }
    best.nclusts <- rbind(best.nclusts, as.data.frame(t(nclusts)))
  }
  #names(nclusts) <- indexes
  #print(nclusts)
  #print(median(nclusts))
  colnames(best.nclusts) <- indexes
  print("NbClust results:")
  print(best.nclusts)
  best <- median(apply(best.nclusts, 2, median))
  print(paste0("Best cluster number for dataset ", dataset, ": ", best))
  return(best)
}

run_spaces <- function(path, datasets, size=2000, K=5, trials=10){
  params <- list("K"=K, "objectives"="XieBeni")
  size <- size*trials
  spaces <- data.frame()
  spaces_obj <- data.frame()
  spaces_perc <- data.frame()
  mstknn.results <- data.frame()
  for(j in 1:length(datasets)){
    dataset <- datasets[j]
    output.folder <- file.path(path)
    distances <- load.gene.distance(dataset, alpha=0.5)
    #for(i in 1:trials){
      #output.exp <- file.path(output.folder, i)
    print("Dataset:")
    print(dataset)
    population <- generate_initial_pop(size, K, distances$n.genes) 
    P.data <- cluster_data(distances, population, 0.5)
    population <- P.data$population
    clustering.groups <- P.data$clustering
    row.names(population) <- 1:nrow(population)
    P.rows <- row.names(population)
    population <- remove_duplicated(population, K)
    clustering.groups <- clustering.groups[match((row.names(population)), P.rows)]
    P.rows <- row.names(population)
    
    print("Before evaluation")
    P <- evaluate_population(population, distances, clustering.groups, params)
    P.clustering.groups <- clustering.groups[match((row.names(P)), P.rows)]
    obj <- P[, c("f1", "f2")]
    obj <- normalizeFront(obj)
    #front <- P[P$rnkIndex == 1, ]
    #if(nrow(front) < 3){
    #  next
    #}
    #print(front)
    #clustering.front <- P.clustering.groups[as.numeric(row.names(front))]
    print("After evaluation")
    d.obj <- as.matrix(dist(obj, upper=T,diag=T))
    diversity.jaccard <- as.matrix(calculate_diversity_matrix(P.clustering.groups, "jaccard"))
    diversity.nmi <- as.matrix(calculate_diversity_matrix(P.clustering.groups, "NMI"))
    
    print(dataset)
    res <- spaces_dist_ratio(diversity.jaccard, diversity.nmi, d.obj)
    print("After spaces dist ratio")
    spaces_tmp <- cbind(rep(dataset, nrow(res$spaces)), res$spaces)
    spaces <- rbind(spaces, spaces_tmp)
    spaces_obj_tmp <- cbind(dataset, res$spaces_obj)
    spaces_obj <- rbind(spaces_obj, spaces_obj_tmp)
    spaces_perc_tmp <- cbind(dataset, res$percentile_dist)
    spaces_perc <- rbind(spaces_perc, spaces_perc_tmp)
    #print("Initiating MSTKNN experiments...")
    #mstknn_tmp <- spaces_mstknn(diversity.jaccard, diversity.nmi, d.obj, "spaces", dataset)
    #mstknn.results <- rbind(mstknn.results, mstknn_tmp)
  }
  #write.table(mstknn.results, file.path(output.folder, paste0("spaces_mstknn.csv")))
  write.table(spaces, file.path(output.folder, paste0("spaces_dist_ratio.csv")))
  write.table(spaces_obj, file.path(output.folder, paste0("spaces_obj.csv")))
  write.table(spaces_perc, file.path(output.folder, paste0("spaces_perc.csv")))
  print("Finished!")
}

spaces_mstknn_plot <- function(output.folder, comp.measure="rand"){
  spaces <- read.table(file.path(output.folder, paste0("spaces_mstknn.csv")), sep=" ")
  colnames(spaces) <- c("Dataset", "Metric", "Value")
  ggplot() +
    geom_boxplot(data=spaces, aes(x=Value, fill=Metric), alpha=0.5) +
    labs(title="MSTKNN comparison of objective and decision space distances",
         fill="Metric", x="Value", y=paste(comp.measure, "comparison measure")) +
    facet_wrap(~Dataset) +
    coord_flip() +
    theme(axis.text.x=element_blank())
  ggsave(file.path(output.folder, paste0("spaces_mstknn_", comp.measure,".png")))
}


spaces_dist_ratio_plot <- function(output.folder, plot="scaled"){
  require(MASS)
  require(scales)
  require(ggridges)
  require(GGally)
  spaces <- read.table(file.path(output.folder, paste0("spaces_dist_ratio.csv")), sep=" ")
  spaces_obj <- read.table(file.path(output.folder, paste0("spaces_obj.csv")), sep=" ")
  spaces_perc <- read.table(file.path(output.folder, paste0("spaces_perc.csv")), sep=" ")
  samples<- c(which(spaces_obj$dataset == "arabidopsis")[1:100],
              which(spaces_obj$dataset == "cell_cycle")[1:100],
              which(spaces_obj$dataset == "serum")[1:100],
              which(spaces_obj$dataset == "sporulation")[1:100])
  print("Obj sample size:")
  print(length(samples))
  spaces_obj <- spaces_obj[ samples, ]
  
  #spaces_obj <- cbind(spaces_obj, "y"=rep(0, nrow(spaces_obj)))
  colnames(spaces) <- c("Dataset", "Metric", "Value", "Obj_dist", "Dec_dist")  
  colnames(spaces_obj) <- c("Dataset", "Metric", "Value", "Minimum_distance")  
  #colnames(spaces_perc) <- c("Dataset", "Objective", "Jaccard", "NMI", "Sample")
  colnames(spaces_perc) <- c("Dataset", "Objective", "Decision", "Sample")
  dir.create(file.path(output.folder), recursive = TRUE, showWarnings = FALSE)
  
  
  # This doesnt make any sense
  #datasets <- c("arabidopsis", "cell_cycle","serum","sporulation")
  #for(i in 1:length(datasets)){
  #  dataset <- datasets[i]
  #  print(paste0("t-test for",dataset, "(non-scaled): "))
  #  data <- spaces[spaces$Dataset == dataset, "Value"]
  #  res <- t.test(data, alternative="two.sided") 
  #  print(res)
  #}
  
  if(plot == "scaled"){
    ggplot() +
      geom_density(data=spaces, aes(x=Value, fill=Metric), alpha=0.7) +
      geom_dotplot(data=spaces_obj, aes(x=Value, color=Minimum_distance),
                   stackratio=1, dotsize=0.8, alpha=0.35, stackgroups = TRUE) + 
      labs(title="Objective space distance to decision space distance ratio",fill="Metric", x="Ratio value (Log10-scaling)", y="Density") +
      #geom_jitter(data=spaces_obj, aes(x=Value, y=y, color=Metric), height = 0.05) +
      #geom_density_ridges(data=spaces_obj, aes(x=Value, y=y, color=Metric),
      #  jittered_points = TRUE,
      #  position = position_points_jitter(width = 0.05, height = 0),
      #  point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
      #) + 
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
      annotation_logticks(sides="b")
      #scale_x_continuous("Ratio value (Log10-scaling)", trans='log10')
    ggsave(file.path(output.folder, paste0("spaces_dist_ratio_log.png")))
  }else if(plot == "non_scaled"){
    ggplot() +
      geom_density(data=spaces, aes(x=Value, fill=Metric), alpha=0.5) +
      geom_dotplot(data=spaces_obj, aes(x=Value, color=Minimum_distance),
                   stackratio=1, dotsize=0.8) + 
      labs(title="Objective space distance to decision space distance ratio",fill="Metric", x="Ratio value", y="Density") +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2)
    ggsave(file.path(output.folder, paste0("spaces_dist_ratio.png")))
  }else if(plot == "percentile"){
    n.samp <- 50
    samples<- c(which(spaces_perc$Dataset == "arabidopsis" & spaces_perc$Sample == "Objective")[1:n.samp],
                which(spaces_perc$Dataset == "arabidopsis" & spaces_perc$Sample == "Jaccard")[1:n.samp],
                which(spaces_perc$Dataset == "arabidopsis" & spaces_perc$Sample == "NMI")[1:n.samp],
                which(spaces_perc$Dataset == "cell_cycle" & spaces_perc$Sample == "Objective")[1:n.samp],
                which(spaces_perc$Dataset == "cell_cycle" & spaces_perc$Sample == "Jaccard")[1:n.samp],
                which(spaces_perc$Dataset == "cell_cycle" & spaces_perc$Sample == "NMI")[1:n.samp],
                which(spaces_perc$Dataset == "serum" & spaces_perc$Sample == "Objective")[1:n.samp],
                which(spaces_perc$Dataset == "serum" & spaces_perc$Sample == "Jaccard")[1:n.samp],
                which(spaces_perc$Dataset == "serum" & spaces_perc$Sample == "NMI")[1:n.samp],
                which(spaces_perc$Dataset == "sporulation" & spaces_perc$Sample == "Objective")[1:n.samp],
                which(spaces_perc$Dataset == "sporulation" & spaces_perc$Sample == "Jaccard")[1:n.samp],
                which(spaces_perc$Dataset == "sporulation" & spaces_perc$Sample == "NMI")[1:n.samp]
                )
    spaces_perc <- spaces_perc[samples, ]
    print(length(samples))
    ggparcoord(data = spaces_perc,
               columns = c(2,3),
               alphaLines = 0.35,
               boxplot = FALSE,
               showPoints = TRUE,
               scale = "globalminmax",
               groupColumn = "Sample") +
      scale_color_brewer(palette = "Set2") +
      facet_wrap(~Dataset,
                 labeller = labeller(Dataset=c(
                   "1" = "Arabidopsis Thaliana",
                   "2" = "Yeast cell cycle",
                   "3" = "Fibroblast serum",
                   "4" = "Yeast Sporulation"
                 ))) +
      labs(title = "Percentile asociation of decision and objective space distances",
        subtitle = paste0("Sampled ", n.samp, " lowest solution distances")) +
      xlab("Distance Metric") + 
      ylab("Percentile")
    ggsave(file.path(output.folder, paste0("spaces_percentile_dist.png")))
  }
  #colnames(spaces) <- c("dataset", "metric", "value")  
  #ggplot(spaces, aes(x=value, color=metric)) +
  #  geom_hist() +
  #  facet_wrap(~dataset, ncol = 2)
  #dir.create(file.path(output.folder), recursive = TRUE, showWarnings = FALSE)
  #ggsave(file.path(output.folder, paste0("spaces_mstknn.png")))
  #write.table(spaces, file.path(output.folder, paste0("spaces_mstknn.csv")))
}

spaces_dist_ratio <- function(d.jaccard, d.nmi, d.obj){

  comb <- combn(1:nrow(d.obj),2)
  n <- ncol(comb)
  ratio.jaccard <- integer(n)
  ratio.nmi <- integer(n)
  obj.dist <- integer(n)
  #obj <- front[, c("f1", "f2")]
  #d.obj <- as.matrix(dist(obj, upper=T,diag=T))
  for(i in 1:ncol(comb)){
    obj.dist[i] <- d.obj[comb[1, i], comb[2, i]]
    ratio.jaccard[i] <- obj.dist[i]/d.jaccard[comb[1, i], comb[2, i]]
    ratio.nmi[i] <- obj.dist[i]/d.nmi[comb[1, i], comb[2, i]]
  }
  #min.obj <- min(obj.dist)
  max.indexes <- round(n/30)
  index.min <- order(obj.dist)[1:max.indexes]
  value.min.jac <- ratio.jaccard[index.min]
  value.min.nmi <- ratio.nmi[index.min]
  distance.jaccard <- obj.dist/ratio.jaccard
  distance.nmi <- obj.dist/ratio.nmi

  spaces <- data.frame("metric"=c( rep("Jaccard", ncol(comb)),
                                   rep("NMI", ncol(comb)) ),
                       "value"=c(ratio.jaccard, ratio.nmi),
                       "obj"=c(obj.dist, obj.dist), 
                       "dec"=c(distance.jaccard, distance.nmi))
  spaces_2 <- data.frame("metric"=c("Jaccard", "NMI"),
                         "value"=c(value.min.jac, value.min.nmi),
                         "obj"=rep(obj.dist[index.min], 2))
  
  perc_dist <- function(vector, sample){
    out <- as.numeric(lapply(sample, function(x) sum(x >= vector)))/length(vector)
    return(as.numeric(out))
  }
  
  #Objective distance sampled
  obj.perc <- perc_dist(obj.dist, obj.dist[index.min])
  jacc.perc <- perc_dist(distance.jaccard, distance.jaccard[index.min])
  nmi.perc <- perc_dist(distance.nmi, distance.nmi[index.min])
  
  index.min.jacc <- order(distance.jaccard)[1:max.indexes]
  index.min.nmi <- order(distance.nmi)[1:max.indexes]
  #Jaccard distance sampled
  jacc.perc.1 <- perc_dist(distance.jaccard, distance.jaccard[index.min.jacc])
  nmi.perc.1 <- perc_dist(distance.nmi, distance.nmi[index.min.jacc])
  obj.perc.1 <- perc_dist(obj.dist, obj.dist[index.min.jacc])
  #NMI distance sampled
  jacc.perc.2 <-perc_dist(distance.jaccard, distance.jaccard[index.min.nmi])
  nmi.perc.2 <- perc_dist(distance.nmi, distance.nmi[index.min.nmi])
  obj.perc.2 <- perc_dist(obj.dist, obj.dist[index.min.nmi])
  
  #perc.dist <- data.frame("perc_obj"=c(obj.perc, obj.perc.1,obj.perc.2),
  #                        "jacc_perc"=c(jacc.perc,jacc.perc.1,jacc.perc.2),
  #                        "nmi_perc"=c(nmi.perc,nmi.perc.1,nmi.perc.2), 
  #                        "Sample"=rep(c("Objective","Jaccard", "NMI"), each=max.indexes))
  perc.dist <- data.frame("perc_obj"=c(obj.perc, obj.perc, obj.perc.1, obj.perc.2, obj.perc.2,obj.perc.2),
                          "perc_dec"=c(jacc.perc, nmi.perc, jacc.perc.1, nmi.perc.1, jacc.perc.2, nmi.perc.2),
                          "Sample"=rep(c("Objective","Jaccard", "NMI"), each=max.indexes*2))
  return(list("spaces"=spaces, "spaces_obj"=spaces_2, "percentile_dist"=perc.dist))
}

spaces_dendogram <- function(diversity.jaccard, diversity.nmi, front){
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
  pareto <- front[, c("f1", "f2")] #normalise_pareto(front[, c("f1", "f2")], limits)
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
  
}

spaces_network <- function(diversity.jaccard, diversity.nmi, front){
  #library("mstknnclust")
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

spaces_mstknn <- function(diversity.jaccard.tot, diversity.nmi.tot, d.obj.tot, output, dataset){
  library("mstknnclust")
  library("igraph")
  path <- file.path(output, dataset)
  print(paste0("Length nmi: ", nrow(diversity.nmi.tot), ", ", ncol(diversity.nmi.tot)))
  print(paste0("Length jaccard: ", nrow(diversity.jaccard.tot), ", ", ncol(diversity.jaccard.tot)))
  print(paste0("Length objective: ", nrow(d.obj.tot), ", ", ncol(d.obj.tot)))
  
  samples <- 20
  trials <- 31
  i <- 1
  mstknn.distances <- data.frame()
  while(i <= trials){
    dir.create(file.path(path,i), recursive = TRUE)
    print(paste0("Starting trial ", i))
    s.indexes <- sample(1:nrow(d.obj.tot), samples)
    d.obj <- d.obj.tot[s.indexes, s.indexes]
    diversity.jaccard <- diversity.jaccard.tot[s.indexes, s.indexes]
    diversity.nmi <- diversity.nmi.tot[s.indexes, s.indexes]
    
    knn.obj <- mstknnclust::mst.knn(d.obj)
    knn.jaccard <- mstknnclust::mst.knn(diversity.jaccard)
    knn.nmi <- mstknnclust::mst.knn(diversity.nmi)
    
    #NMI
    igraph::V(knn.nmi$network)$label.cex <- seq(1.2,1.2,length.out=2)
    
    png(file.path(path, i,"mstknn_nmi.png"))
    
    plot(knn.nmi$network, vertex.size=30,
         vertex.color=igraph::clusters(knn.nmi$network)$membership, 
         layout=igraph::layout.fruchterman.reingold(knn.nmi$network, niter=10000),
         main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisi/U00F3n (NMI)"))
    
    #ggsave(file.path(path, "mstknn_jaccard.png"))#, height=5, width=5)
    #ggsave(paste0("jaccard_", dataset,".png"), height=5, width=8)
    
    #p1.jaccard <- recordPlot()
    dev.off()
    #plot.new()  
    
    #Jaccard
    igraph::V(knn.jaccard$network)$label.cex <- seq(1.2,1.2,length.out=2)
    
    png(file.path(path, i,"mstknn_jaccard.png"))
    
    plot(knn.jaccard$network, vertex.size=30,
         vertex.color=igraph::clusters(knn.jaccard$network)$membership, 
         layout=igraph::layout.fruchterman.reingold(knn.jaccard$network, niter=10000),
         main=paste("MST-kNN \n Soluciones pareto \n Espacio de decisi/U00F3n (Jaccard)"))
    
    #ggsave(file.path(path, "mstknn_jaccard.png"))#, height=5, width=5)
    #ggsave(paste0("jaccard_", dataset,".png"), height=5, width=8)
    
    #p1.jaccard <- recordPlot()
    dev.off()
    #plot.new()  
    
    #Objective
    igraph::V(knn.obj$network)$label.cex <- seq(1.2,1.2,length.out=2)
    
    png(file.path(path, i,"mstknn_obj.png"))
    
    plot(knn.obj$network, vertex.size=30, 
         vertex.color=igraph::clusters(knn.obj$network)$membership, 
         layout=igraph::layout.fruchterman.reingold(knn.obj$network, niter=10000),
         main=paste("MST-kNN \n Soluciones pareto \n Espacio objetivo"))
    
    #ggsave(file.path(path, "mstknn_obj.png"))#, height=5, width=8)
    #ggsave(paste0("nmi_", dataset,".png"), height=5, width=8)
    
    #p1.nmi <- recordPlot()
    dev.off()
    
    comp.nmi <- tryCatch({
      compare(clusters(knn.nmi$network)$membership, clusters(knn.obj$network)$membership, method="vi")
      },
      error=function(cond){
        print("ERROR 1")
        return(NULL)
      })
    print("comparison of network nmi:")
    print(comp.nmi)
    comp.jaccard <- tryCatch({
        compare(clusters(knn.jaccard$network)$membership, clusters(knn.obj$network)$membership, method="vi")
      },
      error=function(cond){
        print("ERROR 2")
        return(NULL)
      })
    print("comparison of network jaccard:")
    print(comp.jaccard)
    if(!(is.null(comp.nmi) || is.null(comp.jaccard))){
      i <- i + 1
      mstknn.distances.trial <- data.frame("Dataset"=dataset, "Metric"=c("jaccard", "NMI"), "Value"=c(comp.jaccard,comp.nmi))
      mstknn.distances <- rbind(mstknn.distances, mstknn.distances.trial)
    }else{
      print("Error found!!, restarting trial...")
    }
  }
  return(mstknn.distances)
}

run_spaces_plot_2 <- function(distances, output.exp, dataset, P=5000, K=5, trials){
  
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

execute_tests <- function(output.folder, datasets, trials){
  distances <- load.gene.distance(dataset, alpha=0.5)
  #dir.create(output.exp, showWarnings=FALSE, recursive=TRUE)
  run_spaces(output.folder, datasets, size=50, K=5, trials=trials)
  #res2 <- run_cluster_n(distances, output.folder, dataset, trials=trials) 
}

spaces_comparison()