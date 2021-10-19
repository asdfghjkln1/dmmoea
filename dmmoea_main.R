print("a")
########################################################################################
#           Vicente Rivera Rebolledo                                                   #
#           vicente.rivera.r@usach.cl                                                  #
# DMMOEA: Diverse Multimodal Multiobjective Evolutionary Algorithm for Gene Clustering #
########################################################################################

# Set working directory
dir <- rstudioapi::getActiveDocumentContext()$path #*** Doesnt work when sourcing ***
dir <- paste0(dirname(dir), "/")
setwd( dir )

source("dmmoea_functions.R")

source("dmmoea_parameters.R")

source("dmmoea_libraries.R")

source("dmmoea_distances.R")

params <- parameters()
distances <- load.gene.distance("arabidopsis")

diversity.metric <- "jaccard"
diversity.level <- 3

output.base <- file.path(dir, "/Tests/arabidopsis")
  
path.nsga2 <- "nsga2/test"
path.dnsga2 <- "dnsga2/test"
path.dmnsga2 <- "dmnsga2/test"

output.path.nsga2 <- file.path(output.base, path.nsga2)
output.path.dnsga2 <- file.path(output.base, path.dnsga2)
output.path.dmnsga2 <- file.path(output.base, path.dmnsga2)
K <- 4

start_time <- Sys.time()
# Diverse memetic NSGA-2 tests
#source("dmmoea_functions.R")
pareto.results.1 <- diverse_memetic_nsga2(distances, K, diversity.metric, diversity.level,  params, output.path.dmnsga2)

end_time <- Sys.time()

# Extracting results
results <- evaluate_solutions(pareto.results.1$population, pareto.results.1$clustering, K, dist, params, output.path.dnsga2)

print(paste("DMNSGA-2 ended!. Execution finished in", round(end_time - start_time, 2), "seconds"))



start_time <- Sys.time()

pareto.results.2 <- dnsga2(distances, K, diversity.metric, diversity.level, params, output.path.dnsga2)

end_time <- Sys.time()

# Extracting results
#source("dmmoea_functions.R")
results <- evaluate_solutions(pareto.results.2$population, pareto.results.1$clustering, K, distances, params, output.path.dnsga2)

print(paste("DNSGA-2 ended!. Execution finished in", round(end_time - start_time, 2), "minutes"))


#projected.pareto <- hipervolume_projection(pareto[, (4+1):(4+params$objDim)], output)

# Original NSGA-II Tests

# Seeds used for arabidopsis: 4321 for test A, 1234 for test B

pareto.results.2 <- nsga2(dist, k, params, output.nsga2)

results <- evaluate_solutions(pareto.results.2$population, pareto.results.2$clustering, k, dist, alpha, output.nsga2)
#pareto <- NSGA2(initial_pareto, 1, param$popSize, param$evaluations, "B_main_log_moc.txt", test_output, experiment)
#write.table(pareto, path=paste0(test_path, , sep = ",",col.names = TRUE, quote = FALSE)

# Test with CRAN package no funciona...
#install.packages("moc.gapbk", dependencies = TRUE)
#library("amap")
#library("moc.gapbk")