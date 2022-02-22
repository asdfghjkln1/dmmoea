calculate_single_hv <- function(){
  library(ggplot2)
  #args <- commandArgs(trailingOnly = TRUE)
  path <- "X:/Universidad/dmmoea" #args[1]
  results.path <- "Tests/runs/BallHall" #args[2]
  dataset <- "arabidopsis" #args[3]
  exp <- "1" #args[4]
  
  setwd(path)
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(file.path(path, "dmnsga2", dataset, paste0(exp, ".csv")))
  pareto.dmnsga <- read.table(file.path(path, "dmnsga2", dataset, exp, paste0(exp, ".csv")), sep=",", header=FALSE, row.names=NULL)
  pareto.dnsga <- read.table(file.path(path, "dnsga2", dataset, exp, paste0(exp, ".csv")), sep=",", header=FALSE, row.names=NULL)
  pareto.nsga <- read.table(file.path(path, "nsga2", dataset, exp, paste0(exp, ".csv")), sep=",", header=FALSE, row.names=NULL)
  #path <- "X:/Universidad/dmmoea/Tests/runs/BallHall"
  limits <- read.table(file.path(path, "dmnsga2", "arabidopsis", "limits.csv"), sep=",", header=TRUE, row.names=NULL)
  scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
  scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
  #print("DMNSGA2")
  #print(pareto.dmnsga)
  #print(limits)
  pareto.dmnsga <- data.frame("f1"=unlist(lapply(pareto.dmnsga[, 1], scaler.f1)), "f2"=unlist(lapply(pareto.dmnsga[, 2], scaler.f2)))
  #print(pareto.dmnsga)
  #print("DNSGA2")
  #print(pareto.dnsga)
  #print(limits)
  pareto.dnsga <- data.frame("f1"=unlist(lapply(pareto.dnsga[, 1], scaler.f1)), "f2"=unlist(lapply(pareto.dnsga[, 2], scaler.f2)))
  #print(pareto.dnsga)
  #print("NSGA2")
  #print(pareto.nsga)
  #print(limits)
  pareto.nsga <- data.frame("f1"=unlist(lapply(pareto.nsga[, 1], scaler.f1)), "f2"=unlist(lapply(pareto.nsga[, 2], scaler.f2)))
  #print(pareto.nsga)
  
  hv.dmnsga <- calculate_hypervolume(pareto.dmnsga, c(1,1), maximise=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
  hv.dnsga <- calculate_hypervolume(pareto.dnsga, c(1,1), maximise=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
  hv.nsga <- calculate_hypervolume(pareto.nsga, c(1,1), maximise=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
  print(paste0("(dmnsga2, ", exp, ") " ,"hv=", round(hv.dmnsga, 4)))
  print(paste0("(dnsga2, ", exp, ") " ,"hv=", round(hv.dnsga, 4)))
  print(paste0("(nsga2, ", exp, ") " ,"hv=", round(hv.nsga, 4)))
  
  data <- rbind(rbind(pareto.dmnsga, pareto.dnsga), pareto.nsga)
  data[, "class"] <- c( rep("dmnsga2", nrow(pareto.dmnsga)) , rep("dnsga2", nrow(pareto.dnsga)) , rep("nsga2", nrow(pareto.nsga)))
  gg <- ggplot(data, aes(x=f1, y=f2, group=class, color=class)) + 
    geom_point() +
    geom_line(aes(group=class))
  
  plot(gg)
  #hv <- data.frame("Hipervolumen"=c(hv.dmnsga, hv.dnsga, hv.nsga), "Algoritmo"=c("dmnsga2", "dnsga2", "nsga2"))
  #gg <- ggplot(hv, aes(y=Hipervolumen, group=Algoritmo, fill=Algoritmo)) + 
  #  geom_bar()
  #
  #plot(gg)
}

#calculate_single_hv()

test_hypervolume_contribution <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  library(rstatix)
  library(ggplot2)
  library(ggpubr)
  path <- args[1]
  results.path <- args[2]
  dataset <- args[3]
  
  setwd(path)
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  contribution <- read.table(file.path(path, "contribution.csv"), sep=",", header=TRUE, row.names=NULL)
  hypervolume_comparison_tests(contribution, dataset, path)
}

calculate_hypervolume_manual <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  path <- args[1]
  results.path <- args[2]
  dataset <- args[3]

  setwd(path)
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  
  algorithms <- c("dmnsga2", "dnsga2", "tmix", "mfuzz", "moc.gapbk")#"nsga2")

  data <- read.table(file.path(path, "data_pareto.csv"), sep=",", header=TRUE, row.names=NULL)
  
  ideal <- data[data$Algorithm == "Ideal pareto", ]
  ideal <- ideal[!duplicated(ideal[, 1:2]) , ]
  
  contribution <- as.data.frame(matrix(nrow = 0, ncol = 7))
  colnames(contribution) <- c("id", "Algorithm", "Dataset", "hypervolume", "contributed", "total", "ratio")
  print("Calculating contribution to pareto front and hypervolume values...")
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    datasets <- list.dirs(path=file.path(path, algorithm), full.names = FALSE, recursive = FALSE)
    for(j in 1:length(datasets)){
      dataset <- datasets[j]
      exp.path <- file.path(path, algorithm, dataset)
      limits <- read.table(file.path(exp.path, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
      scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
      scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
      ideal.dataset <- ideal[ideal$Dataset == dataset, ]
      experiments <- list.dirs(path=exp.path, full.names = FALSE, recursive = FALSE)
      for(k in 1:length(experiments)){
        experiment <- experiments[k]
        data <- read.table(file.path(exp.path, experiment, paste0(experiment, ".csv")), sep=",", header=FALSE, row.names=NULL)
        nrow.pareto <- nrow(data)
        contrib <- do.call(paste0, data[, 1:2]) %in% do.call(paste0, ideal.dataset[, 1:2])
        contrib <- sum(contrib)
        #print(paste0("Experiment ", experiment, "has contributed ", contrib, " solutions out of ", nrow.pareto))
        data <- data.frame("f1"=unlist(lapply(data[, 1], scaler.f1)), "f2"=unlist(lapply(data[, 2], scaler.f2)))
        hv <- calculate_hypervolume(data, c(1,1), maximise=FALSE) #eaf::hypervolume(data, c(1,1), maximize=FALSE) 
        print(paste0("(", experiment, ") abs=", contrib, " total=", nrow.pareto, " hv=", round(hv, 4)))
        contribution <- rbind(contribution, c(experiment, algorithm, dataset, hv, contrib, nrow.pareto, round(contrib/nrow.pareto, 2)))  
      }
      print(paste(algorithm,dataset, "Mean HV=", mean(hv)))
    }
  }
  colnames(contribution) <- c("id", "Algorithm", "Dataset", "hypervolume", "contributed", "total", "ratio")
  write.table(contribution, file=file.path(path, "contribution.csv"), sep=",", append=FALSE, quote=FALSE, col.names=TRUE, row.names=FALSE)
  hypervolume_comparison_tests(contribution, dataset, path)
}

hypervolume_comparison_tests <- function(data, dataset, output.path, metric="contributed", exp.group="Algorithm"){
  library(dplyr)
  dir.create(output.path, recursive=TRUE, showWarnings = FALSE)
  #print(data[1:20,])
  #data$Algorithm <- factor(data$Algorithm, levels=c("nsga2", "dnsga2", "dmnsga2"))
  
  data <- data[data$Dataset == dataset, ]
  print(paste0("rows: ", nrow(data)))
  algorithms <- as.character(unique(data$Algorithm))
  data$Algorithm <- factor(data$Algorithm, levels=algorithms)
  #data[is.na(data)] <- 0
  #col <- which(colnames(data) == metric)
  #for(i in 1:length(algorithms)){
  #  data.alg <- data[data$Algorithm == algorithms[i], ]
  #  if(sum(data.alg[, col]) == 0){
  #    data <- data[!(data$Algorithm == algorithms[i]), ]
  #  }
  #}
  #algorithms <- unique(data$Algorithm)
  #data$Algorithm <- factor(data$Algorithm, levels=algorithms)
  #print("Levels after:")
  #print(levels(data$Algorithm))
  #print(paste0("rows: ", nrow(data)))
  
  #form <- as.formula(call("~", as.symbol(metric), as.symbol(exp.group)))
  #kruskal.res <- kruskal_test(data, formula=form)
  #pwc <- wilcox_test(data, formula=form, p.adjust.method="bonferroni")
  #pwc <- pwc %>% add_xy_position(x = exp.group)
  #w <- 5 + 0.6*(length(unique(data[, exp.group])) - 3)#ifelse(exp.group>3, 6,5)
  #Y <- pwc$y.position
  #gap.data <- max(data[, metric]) - min(data[, metric])
  #gap <- max(Y[2] - Y[1], gap.data*0.07)
  #pwc$y.position <- Y - gap*seq(from=1, to=0, length.out=length(unique(data[, exp.group])))
  if(dataset == "arabidopsis"){
    dataset.name = "Arabidopsis"
  }else if(dataset == "cell_cycle"){
    dataset.name = "Cell Cycle"
  }else if(dataset == "serum"){
    dataset.name = "Serum"
  }else if(dataset == "sporulation"){
    dataset.name = "Sporulation"
  }
  
  if(metric == "contributed"){
    text.title <- paste0("Soluciones contribuidas a frontera pareto")
  }else if(metric == "ratio"){
    text.title <- paste0("Porcentaje contribuci\U00F3n a frontera pareto")
  }
  
  labels <- c()
  values <- c()
  levels <- levels(data$Algorithm)
  if("nsga2" %in% levels){
    labels <- c(labels, "NSGA-II")
    values <- c(values, "#00AFBB")
  } 
  if("dnsga2" %in% levels){
    labels <- c(labels, "DNSGA-II")
    values <- c(values, "#FC4E07")
  } 
  if("dmnsga2" %in% levels){
    labels <- c(labels, "DMNSGA-II")
    values <- c(values, "#F19B3E")
  }
  if("tmix" %in% levels){
    labels <- c(labels, "Mixture M.")
    values <- c(values, "#02A9EA")
  }
  if("mfuzz" %in% levels){
    labels <- c(labels, "Fuzzy")
    values <- c(values, "#69DC9E")
  }
  if("moc.gapbk" %in% levels){
    labels <- c(labels, "MOCGaPBK")
    values <- c(values, "#E7B800")
  } 
  
  data2 <- data %>% group_by(Algorithm, contributed) %>% summarise("Frecuencia"= n())
  print(data[, c(2,5)])

  ggplot(data2, aes(x=contributed, y=Frecuencia)) +
    geom_bar(aes_string(fill=exp.group), stat="identity",position=position_dodge(), alpha=0.75) +
    labs(#subtitle = get_test_label(kruskal.res, detailed = FALSE, p.col="p.adj"), 
         #caption = get_pwc_label(pwc),
         fill="Algoritmo",
         title=text.title,
         x="Contribuciones a frontera pareto ideal") +
    theme_pubr() +
    scale_fill_manual(labels=labels, #c("NSGA-II", "DNSGA-II", "DMNSGA-II"),
                      values=values) +#c("#00AFBB", "#E7B800", "#FC4E07")) +
    theme(#strip.text.x = element_blank(), 
          #axis.text.x = element_blank(),#element_text(angle=25),
          legend.position="bottom", 
          plot.subtitle=element_text(size=11),
          #legend.spacing.x = unit(0, 'cm'),
          axis.title.x=element_blank()) +
    guides(fill = guide_legend(label.position = "bottom")) +
    #stat_pvalue_manual(pwc, label = "p = {p.adj}", hide.ns = TRUE)
  
  ggsave(file.path(output.path, paste0("contribution_", dataset, ".png")), width = 5, height = 7) 
  
  print(paste(metric, dataset, "... Done."))
}


calculate_hypervolume <- function(pareto, point, maximise=FALSE){
  pareto <- pareto[order(pareto[, 1], decreasing=TRUE), ]
  if(maximise){
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


calculate_hypervolume_manual_2 <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  path <- args[1]
  results.path <- args[2]
  dataset <- args[3]
  
  setwd(path)
  
  path <- file.path(path, results.path)
  print("Path in:")
  print(path)
  
  #library(eaf)
  #source("dmmoea_functions.R")
  #source("dmmoea_libraries.R")
  
  #path <- "~/dmmoea/Tests/runs/BallHall"
  
  algorithms <- c("dmnsga2", "dnsga2", "nsga2")
  #dataset <- "arabidopsis"
  
  limits <- read.table(file.path(path, "limits", dataset, "limits.csv"), sep=",", header=TRUE, row.names=NULL)
  print(limits)
  print(file.path(path, "data_pareto.csv"))
  data <- read.table(file.path(path, "data_pareto.csv"), sep=",", header=TRUE, row.names=NULL)
  scaler.f1 <- function(x){ min(1, (x-limits$min.f1)/(limits$max.f1-limits$min.f1)) }
  scaler.f2 <- function(x){ min(1, (x-limits$min.f2)/(limits$max.f2-limits$min.f2)) }
  
  data <- data[data$Dataset == dataset, ]
  print("Length: ")
  print(nrow(data))
  hv <- c()
  #ideal <- data[data$Algorithm == "Ideal Pareto"]
  #ideal <- ideal[!duplicated(ideal[, 1:2]) , ]
  #print(paste0("Pareto front size: ", nrow(ideal)))
  for(i in 1:length(algorithms)){
    algorithm <- algorithms[i]
    data.alg <- data[data$Algorithm == algorithm, ]
    print(paste0("Length algorithm ", algorithm, ": ", nrow(data.alg)))
    data.alg <- data.alg[!duplicated(data.alg[, 1:2]), ]
    print(paste0("Length algorithm ", algorithm, " (after): ", nrow(data.alg)))
    #print(data.alg)
    
    data.norm <- data.frame("f1"=scaler.f1(data.alg[, 1]), "f2"=scaler.f2(data.alg[, 2]))
    hv[i] <- eaf::hypervolume(data.norm, c(1,1), maximise=FALSE)  #calculate_hypervolume(data.norm, c(1,1), maximise=FALSE) 
    print(paste(algorithm, dataset, "HV =", hv[i]))
  }
  print(algorithms)
  print(hv)
}

test_hypervolume_contribution()

#calculate_hypervolume_manual()
