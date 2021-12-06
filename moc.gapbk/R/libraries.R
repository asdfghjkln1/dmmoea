dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
lib <- Sys.getenv("R_LIBS_USER")

if(!require(amap)){
  install.packages("amap", lib=lib)
  library(amap, lib.loc=lib)
}

if(!require(doSNOW)){
  #, repos = "http://cran.us.r-project.org") 
  install.packages("doSNOW", lib=lib)
  library(doSNOW, lib.loc=lib)
}

if(!require(doMPI)){
  install.packages("doMPI", lib=lib)
  library(doMPI, lib.loc=lib)
}

if(!require(doMPI)){
  install.packages("doMPI", lib=lib)
  library(doMPI, lib.loc=lib)
}
if(!require(amap)){
  install.packages("amap", lib=lib)
  library(amap, lib.loc=lib)
}
if(!require(doMPI)){
  install.packages("doMPI", lib=lib)
  library(doMPI, lib.loc=lib)
}
if(!require(miscTools)){
  install.packages("miscTools", lib=lib)
  library(miscTools, lib.loc=lib)
}
if(!require(fields)){
  install.packages("fields", lib=lib)
  library(fields, lib.loc=lib)
}
if(!require(Rmisc)){
  install.packages("Rmisc", lib=lib)
  library(Rmisc, lib.loc=lib)
}
if(!require(plyr)){
  install.packages("plyr", lib=lib)
  library(plyr, lib.loc=lib)
}
#require(amap)
#require(miscTools)
#require(fields)
#require(Rmisc)
#require(plyr)
#require(doMPI)