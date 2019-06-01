
# Loading libraries

library(RSimpactCyan)
library(RSimpactHelper)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EasyABC)
library(robustbase)

# Define directory

work.dir <-  "/home/dniyukuri/lustre/calibration_phlo_MAR_90"  # on CHPC

# work.dir <- "/home/david/benchmark_master_model/" # on laptop


setwd(paste0(work.dir))


# df <- read.csv("/home/david/benchmark_master_model/Results.benchmark.epi.mm.stats_280_123.csv")

df <- read.csv("/home/dniyukuri/lustre/calibration_phlo_MAR_90/Results.benchmark.epi.mm.stats_280_123.csv")


drops.cols <- c("X", "V1975") 

dr <- df[,!(names(df) %in% drops.cols)] # select(df, -c("X","V1975"))

benchmark_data <- dr

 
MAR.phylo.features.cov.90 <- benchmark_data %>% # 36
  dplyr::select(contains("MAR.a.90"))


med.MAR.phylo.features.cov.90 <- colMedians(as.matrix(MAR.phylo.features.cov.90), na.rm = TRUE)


med.MAR.phylo.features <- med.MAR.phylo.features.cov.90 # c(med.classic.features, med.MAR.phylo.features.cov.90)


# I. one simpact model calibrated to classic summary statistics (full)
######################################################################


# Priors

# True parameter values

# inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 



simpact_prior <-list(c("unif", -2, 0), # dissolution.alpha_0 = - 0.52
                     c("unif", -0.5, 0), # dissolution.alpha_4 = -0.05
                     c("unif", 2.6, 5.0), # formation.hazard.agegapry.baseline = 2
                     c("unif", 0.0, 2.0), # person.agegap.man.dist.normal.mu and ~.woman.~ = 10
                     c("unif", 2, 4), # person.agegap.man.dist.normal.sigma and ~.woman.~ = 5
                     c("unif", 0, 1), # formation.hazard.agegapry.gap_agescale_man and ~_woman = 0.25
                     c("unif", -1, 0), # formation.hazard.agegapry.numrel_man and ~._woman = -0.3
                     c("unif", -0.9, 0), # formation.hazard.agegapry.numrel_diff = -0.1
                     
                     c("unif", -2, 0), # hivtransmission.param.a = -1
                     c("unif", -100, -80), # hivtransmission.param.b = -90
                     c("unif", 0, 1), # hivtransmission.param.c = 0.5
                     c("unif", 0, 0.5), # hivtransmission.param.f1 = 0.05
                     c("unif", -0.5, 0), # hivtransmission.param.f2 = -0.14
                     
                     c("unif", 3, 7), # person.vsp.toacute.x = 5
                     c("unif", 5, 9), # person.vsp.toaids.x = 7
                     c("unif", 10, 14), # person.vsp.tofinalaids.x = 12
                     
                     c("unif", -3.5, -0.7) # conception.alpha_base = -1.7
                     
)


# Calibration with classic features ---------------------------------------

# source("compute.summary.statistics.classic.R")

simpact4ABC.MAR.phylo.90 <- function(inputvector = inputvector){
  
  
  work.dir <-  "/home/dniyukuri/lustre/calibration_phlo_MAR_90"  # on CHPC
  
  # work.dir <- "/home/david/benchmark_master_model/" # on laptop
  
  source("/home/dniyukuri/lustre/calibration_phlo_MAR_90/needed.functions.RSimpactHelp.R")
  source("/home/dniyukuri/lustre/calibration_phlo_MAR_90/compute.summary.statistics.classic.R")
  source("/home/dniyukuri/lustre/calibration_phlo_MAR_90/compute.summary.statistics.phylo.MAR.R")
  source("/home/dniyukuri/lustre/calibration_phlo_MAR_90/advanced.transmission.network.builder.R")
  
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
  # source("/home/david/benchmark_master_model/compute.summary.statistics.classic.R")
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  # source("/home/david/benchmark_master_model/advanced.transmission.network.builder.R")
  
  
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
  
  # set.seed(inputvector[1])
  
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
  art.intro["monitoring.cd4.threshold"] <- 100 # 1200 , # 0 will mean nobody qualifies for ART
  
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
  
  
  # generate.filename <- function(how.long){
  #   
  #   rn <- sample(1:100,1)
  #   t <- as.numeric(Sys.time())
  #   set.seed((t - floor(t)) * 1e8)
  #   chars <- c(letters, LETTERS)
  #   sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
  #   
  #   noise.sample1 <- sample(8:15,1, replace = TRUE)
  #   sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
  #   noise.sample <- sample(1:1000,1)
  #   noise.sample2 <- sample(8:17,1, replace = TRUE)
  #   sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
  #                            paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
  #   
  #   return(sub.dir.sim.id)
  # }
  # 
  # 
  # 
  # sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  
  # identifier <- paste0(seedid)
  # 
  # rootDir <- "/user/scratch/gent/vsc400/vsc40070/EAAA/Fa/temp2" # "/tmp"
  # 
  # destDir <- paste0(rootDir, "/", identifier)
  # 
  # 
  # results <- tryCatch(simpact.run(configParams = cfg.list,
  #                                 destDir = destDir,
  #                                 agedist = age.distr,
  #                                 intervention = ART.factual,
  #                                 seed = seedid,
  #                                 identifierFormat = identifier),
  #                     error = simpact.errFunction)
  
  
  # Running Simpact 
  #################
  
  identifier <- paste0(seedid)
  
  rootDir <- work.dir
  
  destDir <- paste0(rootDir,"/temp/",identifier) 
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = destDir, # sub.dir.rename,
                                  agedist = age.distr,
                                  intervention = intervention,
                                  seed = seedid,
                                  identifierFormat = identifier),
                      error = simpact.errFunction)
  
  if (length(results) == 0){
    # 1 pop growth + 14 age-gender-specific prev + 14 age-gender-specific inc +
    # 8 ART coverage + 1 VL suppression + 28 annual HIV prevalence UNAIDS estimates
    outputvector <- rep(NA, 36)
  } else {
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
      outputvector <- rep(NA, 36)
    } else {
      
      datalist.agemix <- readthedata(results)
      
      # epi.behav.stats <- tryCatch(compute.summary.statistics.classic(datalist = datalist.agemix,
      #                                                                timewindow = c(35, 40)),
      #                             error=function(e) return(rep(NA, 27)))
      
      
      
      simpact.trans.net <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
      
      
      seeds.num <- inputvector[1]
      
      # Sequence simulation is done for at least a transmission network with 6 individuals
      # This means that limitTransmEvents equal at least 7
      
      sequence.simulation.seqgen.par(dir.seq = rootDir,
                                     sub.dir.rename = destDir,
                                     simpact.trans.net = simpact.trans.net, 
                                     seq.gen.tool = "seq-gen",
                                     seeds.num = seeds.num,
                                     endpoint = 40,
                                     limitTransmEvents = 7, # no less than 7
                                     hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                     clust = TRUE) # hiv.seq.file lodged in work.dir
      
      # Transform the sequence format to be handled by ClusterPicker
      sequ.dna <- read.dna(file = paste0(destDir,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
      write.dna(sequ.dna, file = paste0(destDir,"/C.Epidemic.fas") , format = "fasta")
      
      
      mar.phylo.stats <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                       datalist.agemix = datalist.agemix,
                                                                       work.dir = rootDir,
                                                                       sub.dir.rename = destDir,
                                                                       dirfasttree = rootDir,
                                                                       limitTransmEvents = 7,
                                                                       seq.cov = 90,
                                                                       seq.gender.ratio = 0.7,
                                                                       age.group.15.25 = c(15,25),
                                                                       age.group.25.40 = c(25,40),
                                                                       age.group.40.50 = c(40,50),
                                                                       endpoint = 40,
                                                                       timewindow = c(35,40),
                                                                       cut.off = 7),
                                  error=function(e) return(rep(NA, 36)))
      
      outputvector <- mar.phylo.stats # c(epi.behav.stats, mar.phylo.stats)
      
      
      
    }
  }
  
  unlink(paste0(rootDir, "/temp/", identifier), recursive = TRUE)
  
  return(outputvector)
  
}


# wraper.simpact4ABC.classic  <- function(inputvector = inputvector){
#   
#   w.epi.behav.stats <- tryCatch(simpact4ABC.classic(inputvector = inputvector),
#                                 error=function(e) return(rep(NA, 27)))
#   
#   
#   return(w.epi.behav.stats)
#   
# }


# inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
#                  -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -1.7) 

# d <- simpact4ABC.MAR.phylo.90(inputvector = inputvector)

source("/home/dniyukuri/lustre/calibration_phlo_MAR_90/calibration.ABC.R")


calib.MAR.phylo.90.results <- calibration.ABC(model.sim = simpact4ABC.MAR.phylo.90,
                                                    sum_stat_obs = med.MAR.phylo.features,
                                                    simpact_prior = simpact_prior, 
                                                    design.points = 1120,
                                                    seed.val = 6975,
                                                    n_cores = 56)

save(calib.MAR.phylo.90.results, file = "calib.MAR.phylo.90.results.RData")

# 
# calib.class.results <- EasyABC::ABC_sequential(model = simpact4ABC.classic,
#                                                method = "Lenormand",
#                                                prior = simpact_prior,
#                                                summary_stat_target = med.classic.features,
#                                                nb_simul = 4,
#                                                alpha = 0.1,
#                                                p_acc_min = 0.1, # 0.03, #0.05  #0.1,
#                                                use_seed = TRUE,
#                                                seed_count = 222,
#                                                n_cluster = 4,
#                                                inside_prior = TRUE,
#                                                verbose = TRUE)
# 
# save(calib.class.results, file = "calib.class.results.RData")

