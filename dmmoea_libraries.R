dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
lib <- Sys.getenv("R_LIBS_USER")
#if(exists("lib.path")){
#  lib <- lib.path
#}

if(!require(nsga2R)){
  install.packages("nsga2R", lib=lib)
  library(nsga2R, lib.loc=lib)
}
if(!require(cluster)){
  install.packages("cluster" , lib=lib)
  library(cluster, lib.loc=lib)
}
if(!require(factoextra)){
  install.packages("factoextra", lib=lib)
  library(factoextra, lib.loc=lib)
}
if(!require(NbClust)){
  install.packages("NbClust", lib=lib)
  library(NbClust, lib.loc=lib)
}
if(!require(infotheo)){
  install.packages("infotheo", lib=lib)
  library(infotheo, lib.loc=lib)
}
if(!require(irace)){
  install.packages("irace", lib=lib)
  library(irace, lib.loc=lib)
}
if(!require(parallel)){
  install.packages("parallel", lib=lib)
  library(parallel, lib.loc=lib)
}
if(!require(doParallel)){
  install.packages("doParallel", lib=lib)
  library(doParallel, lib.loc=lib)
}
if(!require(clValid)){
  install.packages("clValid", lib=lib)
  library(clValid, lib.loc=lib)
}
if(!require(tictoc)){
  install.packages("tictoc", lib=lib)
  library(tictoc, lib.loc=lib)
}
if(!require(anticlust)){
  install.packages("anticlust", lib=lib)
  library(anticlust, lib.loc=lib)
}
if(!require(scclust)){
  install.packages("scclust", lib=lib)
  library(scclust, lib.loc=lib)
}
if(!require(ggpubr)){
  install.packages("ggpubr", lib=lib)
  library(ggpubr, lib.loc=lib)
}
if(!require(rstatix)){
  install.packages("rstatix", lib=lib)
  library(rstatix, lib.loc=lib)
}
if(!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", lib=lib)
}
if(!require(TMixClust)){
  BiocManager::install("TMixClust", ask=FALSE, update=FALSE)
  library(TMixClust)
}
if(!require(Mfuzz)){
  BiocManager::install("Mfuzz", ask=FALSE, update=FALSE)
  library(Mfuzz)
}
