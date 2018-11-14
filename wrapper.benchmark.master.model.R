

# Define directory

# work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop


work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


# work.dir <- "~/Desktop/calibration/"





wrapper.benchmark.master.model <- function(inputvector = inputvector){
  
  
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/advanced.transmission.network.builder.R")
  # 
  # # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/ALL.funs.R")
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/complete.master.epic.metrics.R")
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/compute.summary.statistics.classic.R")
  # 
  # # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/mixed.effect.fit.transmission.clusters.R")
  # #
  # # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/stats.age.groups.trans.clust.network.fun.R")
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/compute.summary.statistics.phylo.MAR.R")
  # 
  # 
  # source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_benchmark_master_model_13_11_2018/complete.master.epic.metric.class.phylo.features.cov.R")
  # 
  
  # 

  
  
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  
  
  # work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop
  
  
  setwd(paste0(work.dir))
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  library(readr)
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
  
  
results.f <- tryCatch(complete.master.epic.metric.class.phylo.features.cov(inputvector = inputvector),
                      error=function(e) return(rep(NA, 2027)))


return(results.f)


}



reps <- 24



inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 


inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)


epi.mm.stats <- simpact.parallel(model = wrapper.benchmark.master.model,
                                 actual.input.matrix = inputmatrix,
                                 seed_count = 1,
                                 n_cluster = 24)

write.csv(epi.mm.stats, file = "Results.benchmark.epi.mm.stats.csv")


