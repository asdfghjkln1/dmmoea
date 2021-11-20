#!bin/usr/env Rstudio
run_experiments <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  argnum <- length(args)
  if(argnum != 4){
    print(paste0("Not enough parameters (", argnum, "/4)"))
    return(-1)
  }
  algorithm <- args[1]
  path <- args[2] 
  #dataset <- args[2]
  parameters.file <- args[3]
  scenario.file <- args[4]
  setwd(path)
  source("dmmoea_functions.R")
  source("dmmoea_parameters.R")
  source("dmmoea_libraries.R")
  source("dmmoea_distances.R")
  source("dmmoea_irace_conf.R")
  
  test.path <<- file.path(path, "Tests", "tuning")
  dir.create(file.path(test.path, algorithm), recursive=TRUE, showWarnings = FALSE)
  parameters <- readParameters(file = file.path(test.path, parameters.file))
  
  scenario <- readScenario(filename=file.path(test.path, scenario.file))
  scenario$trainInstancesFile <- file.path(test.path, "instances.txt")
  scenario$testInstancesFile <- file.path(test.path, "test_instances.txt")
  #scenario$logFile <- file.path(test.path, algorithm,"irace-log.Rdata")
  if(file.exists(file.path(test.path, algorithm, "irace-log-backup.Rdata"))){
    scenario$recoveryFile <- file.path(test.path, algorithm, "irace-log-backup.Rdata")
  }
  
  
  #checkIraceScenario(scenario, parameters = parameters)
  
  tuned_confs <- irace(scenario = scenario, parameters = parameters)
  print(paste("irace finished!"))

  print("Starting plotting")
  
  plot_experiment_results(file.path(test.path, algorithm))
  plot_algorithm_comparison(test.path)
  plot_algorithm_comparison_pareto(test.path)
  
  print("Plotting finished")
  
  post_analysis(file.path(test.path, algorithm), scenario$logFile)
  
  print("Post analysis finished")
  
  #file.rename(file.path(test.path, algorithm, "irace-log.Rdata"), 
   #           file.path(test.path, algorithm, "irace-log-backup.Rdata"))
}

post_analysis <- function(exp.path, log.file.path){
  exp.path <- file.path(exp.path)
  load(file.path(log.file.path))
  
  if(length(iraceResults$allElites) > 0){
    best.configs <- getFinalElites(iraceResults = iraceResults, n = min(length(iraceResults$allElites), 5))
    #id <- best.configs$.ID.
    # Obtain the configurations using the identifier
    # of the best configuration
    #all.exp <- iraceResults$experiments[,as.character(id)]
  }else{
    best.configs <- iraceResults$allConfigurations
  }
  
  #print("Best configuration(s):")
  #print(all.exp[!is.na(all.exp), ]) 
  write.table(best.configs, file.path(exp.path, "best_configurations.csv"), sep=",", append=FALSE, row.names = FALSE, quote = FALSE)
  
  if(!is.null(iraceResults$testing)){
    results <- iraceResults$testing$experiments
    # Wilcoxon paired test
    conf <- gl(ncol(results), # number of configurations
               nrow(results), # number of instances
               labels = colnames(results))
    
    pairwise.wilcox.test(as.vector(results), conf, paired = TRUE, p.adj = "bonf")
    configurationsBoxplot (results, ylab = "Solution cost") 
  }
  
  parameterFrequency(iraceResults$allConfigurations, filename=file.path(exp.path, "param_frecuency"), iraceResults$parameters)
  
  #if(length(unique(iraceResults$iterationElites)) > 3){
  #Post-selection race  
  #Execute all elite configurations in the iterations
  #  print("Initiating post-selection:")
  #  psRace(iraceLogFile=log.file.path, elites=TRUE)
  #  print("Experiments finished!")
  #}
}


run_experiments()
