dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
lib <- Sys.getenv("R_LIBS_USER")

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
#if(!require(eaf)){
#  install.packages("eaf")
#  library(eaf)
#}
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