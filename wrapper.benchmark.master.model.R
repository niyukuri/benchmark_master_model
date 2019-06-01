# Wrapper function which source all others and the main computatinal complete.master.epic.metric.class.phylo.features.cov script 
# which is the master model

# Define directory

work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC

# work.dir <- "/home/david/benchmark_master_model" # on CHPC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


wrapper.benchmark.master.model <- function(inputvector = inputvector){
  
  
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  # work.dir <- "/home/david/benchmark_master_model" # on CHPC
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop
  
  
  setwd(paste0(work.dir))
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  # library(readr)
  library(phangorn)
  library(lme4)
  library(nlme)
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  library(Rsamtools)
  library(robustbase)
  library(intergraph)
  library(lubridate)
  library(tidyr)
  
  
  source("/home/dniyukuri/lustre/benchmark_master_model/advanced.transmission.network.builder.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epic.metrics.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.classic.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epic.metric.class.phylo.features.cov.R")
  # 
  # 
  
  # source("/home/david/benchmark_master_model/advanced.transmission.network.builder.R") # ok
  # 
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R") # ok
  # 
  # source("/home/david/benchmark_master_model/complete.master.epic.metrics.R") # ok
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.classic.R") # ok
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R") # ok
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MAR.R") # ok
  # 
  # source("/home/david/benchmark_master_model/complete.master.epic.metric.class.phylo.features.cov.R") # ok
  
  
  
  results.f <- tryCatch(complete.master.epic.metric.class.phylo.features.cov(inputvector = inputvector),
                        error=function(e) return(rep(NA, 1974)))
  
  
  return(results.f)
  
  
}



reps <- 56



inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 


inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)


epi.mm.stats <- simpact.parallel(model = wrapper.benchmark.master.model,
                                 actual.input.matrix = inputmatrix,
                                 seed_count = 123,
                                 n_cluster = 56)

epi.mm.stats <- as.data.frame(epi.mm.stats)

write.csv(epi.mm.stats, file = "Results.benchmark.epi.mm.stats_TEST.csv")


