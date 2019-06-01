


# Define directory

# work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC

work.dir <- "/home/niyukuri/Desktop/mastermodeltest" # on laptop


setwd(paste0(work.dir))

#work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools, dplyr, robustbase, abc, lhs, dplyr)


# source("/home/dniyukuri/lustre/calibration_stress_testing/calibration.ABC.R")
# source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")

source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/calibration.ABC.R")
source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")

# df <- read.csv("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_calibration_stress_testing_14_11_2018/Results.epi.mm.stats_CHPC.csv")

# df <- read.csv("/home/dniyukuri/lustre/calibration_stress_testing/Results.epi.mm.stats_CHPC.csv")

df <- read.csv("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/Results.epi.mm.stats_CHPC.csv")




# colMedians(x, na.rm = FALSE, hasNA = TRUE, keep.names=TRUE)
# rowMedians(x, na.rm = FALSE, hasNA = TRUE, keep.names=TRUE)

# epi.metrics <- df %>%
#   select(contains("metr."))

classic.features <- df[,41:67] # 27

med.ss.classic.features <- colMedians(as.matrix(classic.features), na.rm = TRUE)

# 
# MCAR.phylo.features.cov.35 <- df %>% # 37
#   select(contains("MCAR.35"))
# 
# MCAR.phylo.features.cov.40 <- df %>%
#   select(contains("MCAR.40"))
# 
# MCAR.phylo.features.cov.45 <- df %>%
#   select(contains("MCAR.45"))
# 
# MCAR.phylo.features.cov.50 <- df %>%
#   select(contains("MCAR.50"))
# 
# MCAR.phylo.features.cov.55 <- df %>%
#   select(contains("MCAR.55"))
# 
# MCAR.phylo.features.cov.60 <- df %>%
#   select(contains("MCAR.60"))
# 
# MCAR.phylo.features.cov.65 <- df %>%
#   select(contains("MCAR.65"))
# 
# MCAR.phylo.features.cov.70 <- df %>%
#   select(contains("MCAR.70"))
# 
# MCAR.phylo.features.cov.75 <- df %>%
#   select(contains("MCAR.75"))
# 
# MCAR.phylo.features.cov.80 <- df %>%
#   select(contains("MCAR.80"))
# 
# MCAR.phylo.features.cov.85 <- df %>%
#   select(contains("MCAR.85"))
# 
# MCAR.phylo.features.cov.90 <- df %>%
#   select(contains("MCAR.90"))

# MCAR.phylo.features.cov.95 <- df %>%
#   select(contains("MCAR.95"))

MCAR.phylo.features.cov.95 <- df[,512:548] # names( df[,512:548])

med.ss.phylo.cov.95.features <- colMedians(as.matrix(MCAR.phylo.features.cov.95), na.rm = TRUE)

# I. one simpact model calibrated to classic summary statistics (full)
######################################################################


# Priors

# c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
#   -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 


simpact_prior <- list(c("unif", -1, 0), c("unif", -0.5, 0), c("unif", 1, 3), c("unif", 5, 15), c("unif", 2, 8), c("unif", 0, 1),
                      c("unif", -1, 0), c("unif", -0.9, 0), c("unif", -2, 0), c("unif", -100, -80), c("unif", 0, 1), c("unif", 0, 0.5),
                      c("unif", -0.5, 0), c("unif", 3, 7), c("unif", 5, 9), c("unif", 10, 14), c("unif", -2.7, -1.07))


# Calibration with classic features ---------------------------------------



simpact4ABC.classic <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 5000, 
                                   population.numwomen = 5000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 40, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  ABC_DestDir.classic <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Error function when computing summary statistics
  
  err.functionGEN <- function(e){
    return(chunk.summary.stats = rep(NA,27))
    stop(e)
  }
  
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = ABC_DestDir.classic,
                                  agedist = age.distr,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 27) 
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 27)
    }else{
      
      
      datalist <- readthedata(results)
      
      simpact.trans.net <-  transmission.network.builder(datalist = datalist, endpoint = 40)
      
      
      summary.stat.classic <- tryCatch(compute.summary.statistics.classic(datalist = datalist.agemix,
                                                                          timewindow = c(35, 40)),
                                       error=function(e) return(rep(NA, 27))) # len = 27
      
      
      
      outputvector <- summary.stat.classic
      
    }
    
  }
  
  unlink(paste0(ABC_DestDir.classic), recursive = TRUE)
  
  return(outputvector)
}



med.ss.classic.features <- as.numeric(med.ss.classic.features)

cal.stats <- calibration.ABC(model.sim = simpact4ABC.classic,
                             sum_stat_obs = med.ss.classic.features,
                             simpact_prior = simpact_prior, 
                             design.points = 100,
                             seed.val = 1000,
                             n_cores = 8)


save(cal.stats, "cal.stats.RData")


parm.sel.matrix <- as.matrix(cal.stats$adj.values)


# computations with complete.master.epic.metrics.R

master.epic.metrics.simpact4ABC.cal.class <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 5000, 
                                   population.numwomen = 5000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 40, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  ABC.cal.Dir.clas <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Error function when computing summary statistics
  
  err.functionGEN <- function(e){
    return(chunk.summary.stats = rep(NA,27))
    stop(e)
  }
  
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = ABC.cal.Dir.clas,
                                  agedist = age.distr,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 39) 
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 39)
    }else{
      
      
      datalist <- readthedata(results)
      
      simpact.trans.net <-  transmission.network.builder(datalist = datalist, endpoint = 40)
      
      
      epi.metrics.cal <- tryCatch(complete.master.epic.metrics(datalist = datalist.agemix),
                                  error=function(e) return(rep(NA, 39))) # len = 27
      
      
      
      outputvector <- epi.metrics.cal
      
    }
    
  }
  
  unlink(paste0(ABC.cal.Dir.clas), recursive = TRUE)
  
  return(outputvector)
}



wrapper.master.epic.metrics.simpact4ABC.cal.class <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  results.f <- tryCatch(master.epic.metrics.simpact4ABC.cal.class(inputvector = inputvector),
                        error=function(e) return(rep(NA, 39)))
  
  
  return(results.f)
  
  
}



epi.metrics.post.calib.ABC.class <- simpact.parallel(model = wrapper.master.epic.metrics.simpact4ABC.cal.class,
                                                     actual.input.matrix = parm.sel.matrix,
                                                     seed_count = 1,
                                                     n_cluster = 24)


write.csv(epi.metrics.post.calib.ABC.class, file = "Results.epi.metrics.post.calib.ABC.class.csv")




# Calibration with combined features (classic and phylo) ------------------


simpact4ABC.classic.phylo <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 5000, 
                                   population.numwomen = 5000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 40, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  ABC_DestDir.classic.phylo <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Error function when computing summary statistics
  
  err.functionGEN <- function(e){
    return(chunk.summary.stats = rep(NA,27))
    stop(e)
  }
  
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = ABC_DestDir.classic.phylo,
                                  agedist = age.distr,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 64) 
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 64)
    }else{
      
      
      datalist <- readthedata(results)
      
      simpact.trans.net <-  transmission.network.builder(datalist = datalist, endpoint = 40)
      
      
      summary.stat.classic <- tryCatch(compute.summary.statistics.classic(datalist = datalist.agemix,
                                                                          timewindow = c(35, 40)),
                                       error=function(e) return(rep(NA, 27))) # len = 27
      
      
      MCAR.cov.phylo.95 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                          datalist.agemix = datalist,
                                                                          work.dir = work.dir,
                                                                          sub.dir.rename = sub.dir.rename,
                                                                          dirfasttree = work.dir,
                                                                          limitTransmEvents = 7,
                                                                          seq.cov = 95,
                                                                          age.group.15.25 = c(15,25),
                                                                          age.group.25.40 = c(25,40),
                                                                          age.group.40.50 = c(40,50),
                                                                          endpoint = 40,
                                                                          timewindow = c(30,40),
                                                                          cut.off = 7),
                                    error=function(e) return(rep(NA, 37)))
      
      outputvector <- c(summary.stat.classic, MCAR.cov.phylo.95)
      
      
    }
    
  }
  
  unlink(paste0(ABC_DestDir.classic.phylo), recursive = TRUE)
  
  return(outputvector)
}



ss.cl.phylo <- as.numeric(c(med.ss.classic.features, med.ss.phylo.cov.95.features))


cal.stats.phylo <- calibration.ABC(model.sim = simpact4ABC.classic.phylo,
                                   sum_stat_obs = ss.cl.phylo,
                                   simpact_prior = simpact_prior, 
                                   design.points = 100,
                                   seed.val = 1000,
                                   n_cores = 24)


save(cal.stats.phylo, "cal.stats.phylo.RData")


parm.sel.comb.matrix <- as.matrix(cal.stats.phylo$adj.values)


# computations with complete.master.epic.metrics.R

master.epic.metrics.simpact4ABC.cal.class <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 5000, 
                                   population.numwomen = 5000,
                                   hivseed.time = 10, 
                                   hivseed.type = "amount",
                                   hivseed.amount = 40, 
                                   hivseed.age.min = 20,
                                   hivseed.age.max = 50,
                                   formation.hazard.agegapry.meanage = -0.025,
                                   debut.debutage = 15
  )
  
  # # Assumption of nature of sexual network
  # #########################################
  #
  cfg.list["population.msm"] = "no"
  
  
  # # Sexual behaviour
  # ###################
  #
  seedid <- inputvector[1]
  
  cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
  cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
  cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
  cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
  cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
  cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
  cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
  cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
  cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
  cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
  cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
  cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
  
  
  # # HIV transmission
  # ###################
  #
  
  cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
  cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
  cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
  cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
  cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
  
  # Disease progression > may be remove in parameter to estimates
  
  cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
  cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
  cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)
  
  
  #
  # # Demographic
  # ##############
  #
  
  cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
  
  
  # # Assumptions to avoid negative branch lengths
  # ###############################################
  # # + sampling == start ART
  # # when someone start ART, he/she is sampled and becomes non-infectious
  
  cfg.list["monitoring.fraction.log_viralload"] <- 0
  
  
  #
  # ## Add-ons
  #
  ### BEGIN Add-on
  cfg.list["formation.hazard.agegapry.baseline"] <- 2
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
  cfg.list["dropout.interval.dist.type"] <- "uniform"
  cfg.list["dropout.interval.dist.uniform.min"] <- 1000
  cfg.list["dropout.interval.dist.uniform.max"] <- 2000
  
  cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
  cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
  cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
  
  cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
  #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
  cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
  #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
  
  cfg.list["mortality.aids.survtime.C"] <- 65
  cfg.list["mortality.aids.survtime.k"] <- -0.2
  cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
  cfg.list["diagnosis.baseline"] <- -2
  
  
  cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
  cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
  cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
  
  #### END Add-ons
  
  
  # # ART intervention
  # ###################
  #
  # # ART acceptability paramter and the ART  interventions
  
  cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.6
  
  # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
  art.intro <- list()
  art.intro["time"] <- 20
  art.intro["diagnosis.baseline"] <- -2 # 0#100
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200
  
  ### add something about diagnosis
  art.intro["diagnosis.agefactor"] <- 0
  art.intro["diagnosis.genderfactor"] <- 0
  art.intro["diagnosis.diagpartnersfactor"] <- 0
  art.intro["diagnosis.isdiagnosedfactor"] <- 0
  ### end of add-on about diagnosis
  #art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"
  # Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
  art.intro1 <- list()
  art.intro1["time"] <- 22
  art.intro1["diagnosis.baseline"] <- -2 # 0#100
  art.intro1["monitoring.cd4.threshold"] <- 150 # 1200
  
  art.intro2 <- list()
  art.intro2["time"] <- 25 # inputvector[5] ######### 30
  art.intro2["monitoring.cd4.threshold"] <- 200
  
  art.intro3 <- list()
  art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
  art.intro3["monitoring.cd4.threshold"] <- 350
  
  art.intro4 <- list()
  art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
  art.intro4["monitoring.cd4.threshold"] <- 500
  
  art.intro5 <- list()
  art.intro5["time"] <- 36
  art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access
  
  # tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status
  interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
  
  intervention <- interventionlist
  
  # Events
  cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
  
  # Avoid overlaping in same directory
  
  #creating subfolder with unique name for each simulation
  generate.filename <- function(how.long){
    
    rn <- sample(1:100,1)
    t <- as.numeric(Sys.time())
    set.seed((t - floor(t)) * 1e8)
    chars <- c(letters, LETTERS)
    sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
    
    noise.sample1 <- sample(8:15,1, replace = TRUE)
    sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
    noise.sample <- sample(1:1000,1)
    noise.sample2 <- sample(8:17,1, replace = TRUE)
    sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                             paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
    
    return(sub.dir.sim.id)
  }
  
  ABC.cal.Dir.class.phylo <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  # Error function when computing summary statistics
  
  err.functionGEN <- function(e){
    return(chunk.summary.stats = rep(NA,27))
    stop(e)
  }
  
  
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = ABC.cal.Dir.class.phylo,
                                  agedist = age.distr,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 39) 
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 39)
    }else{
      
      
      datalist <- readthedata(results)
      
      simpact.trans.net <-  transmission.network.builder(datalist = datalist, endpoint = 40)
      
      
      epi.metrics.cal <- tryCatch(complete.master.epic.metrics(datalist = datalist.agemix),
                                  error=function(e) return(rep(NA, 39))) # len = 27
      
      
      
      outputvector <- epi.metrics.cal
      
    }
    
  }
  
  unlink(paste0(ABC.cal.Dir.class.phylo), recursive = TRUE)
  
  return(outputvector)
}



wrapper.master.epic.metrics.simpact4ABC.cal.class <- function(inputvector = inputvector){
  
  work.dir <-  "/home/niyukuri/Desktop/mastermodeltest"  # on laptop
  
  # work.dir <-  "/home/dniyukuri/lustre/calibration_stress_testing"  # on CHPC
  
  
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
  
  
  # source("/home/dniyukuri/lustre/calibration_stress_testing/advanced.transmission.network.builder.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/complete.master.epic.metrics.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.classic.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  # 
  # source("/home/dniyukuri/lustre/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/advanced.transmission.network.builder.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/complete.master.epic.metrics.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.classic.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/niyukuri/Dropbox/25.10.2018.age.mix2/editings_stress_testing_master_model_14_11_2018/calibration_stress_testing/compute.summary.statistics.phylo.MAR.R")
  
  
  
  results.f <- tryCatch(master.epic.metrics.simpact4ABC.cal.class(inputvector = inputvector),
                        error=function(e) return(rep(NA, 39)))
  
  
  return(results.f)
  
  
}



epi.metrics.post.calib.ABC.class.phylo <- simpact.parallel(model = wrapper.master.epic.metrics.simpact4ABC.cal.class,
                                                           actual.input.matrix = parm.sel.comb.matrix,
                                                           seed_count = 1000,
                                                           n_cluster = 24)


write.csv(epi.metrics.post.calib.ABC.class, file = "Results.epi.metrics.post.calib.ABC.class.phylo.csv")


