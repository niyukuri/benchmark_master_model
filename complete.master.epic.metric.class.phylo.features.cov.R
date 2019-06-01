# Master model which simulates the epidemic and compute summary statistics in different scenarios

# The output is a vector of values of

# 1. Transmission network characteristics 

# - temporal trend of incidence
# - age mixing statistics
# - mean, median, and standard deviation of onward transmissions

# 2. Epidmiological, demographical, sexual behavioural, and interventions realted summary statistics

# 3. Phylogenetic summary statistics in MCAR and MA (35:95, by 5) scenarios (13 * 4), each scenario returns measurements which are describbed
# in compute.summary.statistics.phylo.MAR and compute.summary.statistics.phylo.MCAR scripts):

# Missing Completly at Random has 13 scenarios
# Missing At Random has 39 scenarios, with 13 when we assume we have more women in the sample (seq.gender.ratio = 70%)
# the second we have fewer women (seq.gender.ratio = 30%), and the third we have same amount of men and women (seq.gender.ratio = 50%)




complete.master.epic.metric.class.phylo.features.cov <- function(inputvector = inputvector){
  
  # 
  source("/home/dniyukuri/lustre/benchmark_master_model/advanced.transmission.network.builder.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epic.metrics.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.classic.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")

  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  # 
  # 
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  # 
# source("/home/david/benchmark_master_model/advanced.transmission.network.builder.R")
# 
# source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
#  
# source("/home/david/benchmark_master_model/complete.master.epic.metrics.R")
#  
# source("/home/david/benchmark_master_model/compute.summary.statistics.classic.R")
#  
# source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")
#  
# source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
#   
# work.dir <-  "/home/david/benchmark_master_model"
  
  
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
  
  
  
  sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  
  
  # Running Simpact 
  #################
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  
  datalist.ALL <- readthedata(results)
  
  datalist.agemix <- datalist.ALL
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  simpact.trans.net <- advanced.transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  
  net.size.vector <- vector() # i_th seed in the list of seeds
  
  for(i in 1:length(simpact.trans.net)){
    
    tree.n <- simpact.trans.net[[i]] # transmission network for i^th seed
    
    net.size.vector <- c(net.size.vector, nrow(as.data.frame(tree.n)))
    
  }
  
  
  big.index <- which(net.size.vector>=50)
  
  

    
    
    ###############################
    # Step 3: Sequence simulation #
    ###############################
    
    
    dirseqgen <- work.dir
    
    seeds.num <- inputvector[1]
    
    # Sequence simulation is done for at least a transmission network with 6 individuals
    # This means that limitTransmEvents equal at least 7
    
    sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                                   sub.dir.rename = sub.dir.rename,
                                   simpact.trans.net = simpact.trans.net, 
                                   seq.gen.tool = "seq-gen",
                                   seeds.num = seeds.num,
                                   endpoint = 40,
                                   limitTransmEvents = 7, # no less than 7
                                   hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                   clust = TRUE) # hiv.seq.file lodged in work.dir
    
    # Transform the sequence format to be handled by ClusterPicker
    sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
    write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")
    
    
    #####################################################
    ### I. Compute transmission network characteristics #
    #####################################################
    
    # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metrics.R")
    
    
    epidemic.metrics <- complete.master.epic.metrics(datalist = datalist.agemix)
    
    
    ##################################
    ### II. Compute classic features #
    ################################## ??? change arguments
    
    # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.classic.R")
    
    epi.behav.stats <- compute.summary.statistics.classic(datalist = datalist.agemix,
                                                          timewindow = c(35, 40))
    
    
    #########################################################################
    ## III. Compute phylogenetic features considering missingness scenarios #
    #########################################################################
    
    
    
    # MCAR
    
    MCAR.cov.35 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 35,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.40 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 40,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.45 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 45,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.50 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 50,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.55 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 55,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.60 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 60,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.65 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 65,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.70 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 70,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.75 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 75,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.80 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 80,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.85 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 85,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.90 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 90,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    
    MCAR.cov.95 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 95,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                            error=function(e) return(rep(NA, 36)))
    
    MCAR.cov.100 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                   datalist.agemix = datalist.agemix,
                                                                   work.dir = work.dir,
                                                                   sub.dir.rename = sub.dir.rename,
                                                                   dirfasttree = work.dir,
                                                                   limitTransmEvents = 7,
                                                                   seq.cov = 100,
                                                                   age.group.15.25 = c(15,25),
                                                                   age.group.25.40 = c(25,40),
                                                                   age.group.40.50 = c(40,50),
                                                                   endpoint = 40,
                                                                   timewindow = c(35,40),
                                                                   cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    
    # MCAR.All <- c(MCAR.cov.35, MCAR.cov.40, MCAR.cov.45, MCAR.cov.50, MCAR.cov.55, MCAR.cov.60, MCAR.cov.65, MCAR.cov.70, 
    #               MCAR.cov.75, MCAR.cov.80, MCAR.cov.85, MCAR.cov.90, MCAR.cov.95, MCAR.cov.100)
    
    names.columns <- names(MCAR.cov.100)
    
    
    results.mcar <- as.numeric(c(MCAR.cov.35, MCAR.cov.40, 
                                 MCAR.cov.45, MCAR.cov.50, 
                                 MCAR.cov.55, MCAR.cov.60, 
                                 MCAR.cov.65, MCAR.cov.70, 
                                 MCAR.cov.75, MCAR.cov.80, 
                                 MCAR.cov.85, MCAR.cov.90, 
                                 MCAR.cov.95, MCAR.cov.100))
    
    
    names(results.mcar) <- c(paste0("cov.MCAR.",35,".",paste0(names.columns)), paste0("cov.MCAR.",40,".",paste0(names.columns)),
                             paste0("cov.MCAR.",45,".",paste0(names.columns)), paste0("cov.MCAR.",50,".",paste0(names.columns)),
                             paste0("cov.MCAR.",55,".",paste0(names.columns)), paste0("cov.MCAR.",60,".",paste0(names.columns)),
                             paste0("cov.MCAR.",65,".",paste0(names.columns)), paste0("cov.MCAR.",70,".",paste0(names.columns)),
                             paste0("cov.MCAR.",75,".",paste0(names.columns)), paste0("cov.MCAR.",80,".",paste0(names.columns)),
                             paste0("cov.MCAR.",85,".",paste0(names.columns)), paste0("cov.MCAR.",90,".",paste0(names.columns)),
                             paste0("cov.MCAR.",95,".",paste0(names.columns)), paste0("cov.MCAR.",100,".",paste0(names.columns)))
    
    
    
    
    # MAR
    
    
    
    # a. 0.7
    MAR.a.cov.35 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 35,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.a.cov.40 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 40,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.a.cov.45 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 45,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.a.cov.50 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 50,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.a.cov.55 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 55,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.a.cov.60 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 60,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.a.cov.65 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 65,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    
    MAR.a.cov.70 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 70,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.a.cov.75 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 75,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.a.cov.80 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 80,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.a.cov.85 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 85,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.a.cov.90 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
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
    
    
    MAR.a.cov.95 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 95,
                                                                  seq.gender.ratio = 0.7,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    # MAR.a.All <- c(MAR.a.cov.35, MAR.a.cov.40, MAR.a.cov.45, MAR.a.cov.50, MAR.a.cov.55, MAR.a.cov.60, MAR.a.cov.65, MAR.a.cov.70, 
    #                MAR.a.cov.75, MAR.a.cov.80, MAR.a.cov.85, MAR.a.cov.90, MAR.a.cov.95)
    
    results.mar.a <- as.numeric(c(MAR.a.cov.35, MAR.a.cov.40, 
                                  MAR.a.cov.45, MAR.a.cov.50, 
                                  MAR.a.cov.55, MAR.a.cov.60, 
                                  MAR.a.cov.65, MAR.a.cov.70, 
                                  MAR.a.cov.75, MAR.a.cov.80, 
                                  MAR.a.cov.85, MAR.a.cov.90, 
                                  MAR.a.cov.95))
    
    
    names(results.mar.a) <- c(paste0("cov.MAR.a.",35,".",paste0(names.columns)), paste0("cov.MAR.a.",40,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",45,".",paste0(names.columns)), paste0("cov.MAR.a.",50,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",55,".",paste0(names.columns)), paste0("cov.MAR.a.",60,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",65,".",paste0(names.columns)), paste0("cov.MAR.a.",70,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",75,".",paste0(names.columns)), paste0("cov.MAR.a.",80,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",85,".",paste0(names.columns)), paste0("cov.MAR.a.",90,".",paste0(names.columns)),
                              paste0("cov.MAR.a.",95,".",paste0(names.columns)))
    
    
    
    
    # b.  0.3
    MAR.b.cov.35 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 35,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.b.cov.40 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 40,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.45 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 45,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.50 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 50,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.55 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 55,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.60 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 60,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.65 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 65,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    MAR.b.cov.70 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 70,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.75 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 75,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.80 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 80,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.85 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 85,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.90 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 90,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.b.cov.95 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 95,
                                                                  seq.gender.ratio = 0.3,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    # MAR.b.All <- c(MAR.b.cov.35, MAR.b.cov.40, MAR.b.cov.45, MAR.b.cov.50, MAR.b.cov.55, MAR.b.cov.60, MAR.b.cov.65, MAR.b.cov.70, 
    #                MAR.b.cov.75, MAR.b.cov.80, MAR.b.cov.85, MAR.b.cov.90, MAR.b.cov.95)
    
    results.mar.b <- as.numeric(c(MAR.b.cov.35, MAR.b.cov.40, 
                                  MAR.b.cov.45, MAR.b.cov.50, 
                                  MAR.b.cov.55, MAR.b.cov.60, 
                                  MAR.b.cov.65, MAR.b.cov.70, 
                                  MAR.b.cov.75, MAR.b.cov.80, 
                                  MAR.b.cov.85, MAR.b.cov.90, 
                                  MAR.b.cov.95))
    
    
    names(results.mar.b) <- c(paste0("cov.MAR.b.",35,".",paste0(names.columns)), paste0("cov.MAR.b.",40,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",45,".",paste0(names.columns)), paste0("cov.MAR.b.",50,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",55,".",paste0(names.columns)), paste0("cov.MAR.b.",60,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",65,".",paste0(names.columns)), paste0("cov.MAR.b.",70,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",75,".",paste0(names.columns)), paste0("cov.MAR.b.",80,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",85,".",paste0(names.columns)), paste0("cov.MAR.b.",90,".",paste0(names.columns)),
                              paste0("cov.MAR.b.",95,".",paste0(names.columns)))
    
    
    # c. 0.5
    
    MAR.c.cov.35 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 35,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    MAR.c.cov.40 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 40,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.45 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 45,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.50 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 50,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.55 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 55,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.60 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 60,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.65 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 65,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.70 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 70,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.75 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 75,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.80 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 80,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.85 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 85,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.90 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 90,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    MAR.c.cov.95 <- tryCatch(compute.summary.statistics.phylo.MAR(simpact.trans.net = simpact.trans.net,
                                                                  datalist.agemix = datalist.agemix,
                                                                  work.dir = work.dir,
                                                                  sub.dir.rename = sub.dir.rename,
                                                                  dirfasttree = work.dir,
                                                                  limitTransmEvents = 7,
                                                                  seq.cov = 95,
                                                                  seq.gender.ratio = 0.5,
                                                                  age.group.15.25 = c(15,25),
                                                                  age.group.25.40 = c(25,40),
                                                                  age.group.40.50 = c(40,50),
                                                                  endpoint = 40,
                                                                  timewindow = c(35,40),
                                                                  cut.off = 7),
                             error=function(e) return(rep(NA, 36)))
    
    
    # MAR.c.All <- c(MAR.c.cov.35, MAR.c.cov.40, MAR.c.cov.45, MAR.c.cov.50, MAR.c.cov.55, MAR.c.cov.60, MAR.c.cov.65, MAR.c.cov.70, 
    #                MAR.c.cov.75, MAR.c.cov.80, MAR.c.cov.85, MAR.c.cov.90, MAR.c.cov.95)
    
    results.mar.c <- as.numeric(c(MAR.c.cov.35, MAR.c.cov.40, 
                                  MAR.c.cov.45, MAR.c.cov.50, 
                                  MAR.c.cov.55, MAR.c.cov.60, 
                                  MAR.c.cov.65, MAR.c.cov.70, 
                                  MAR.c.cov.75, MAR.c.cov.80, 
                                  MAR.c.cov.85, MAR.c.cov.90, 
                                  MAR.c.cov.95))
    
    
    names(results.mar.c) <-  c(paste0("cov.MAR.c.",35,".",paste0(names.columns)), paste0("cov.MAR.c.",40,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",45,".",paste0(names.columns)), paste0("cov.MAR.c.",50,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",55,".",paste0(names.columns)), paste0("cov.MAR.c.",60,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",65,".",paste0(names.columns)), paste0("cov.MAR.c.",70,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",75,".",paste0(names.columns)), paste0("cov.MAR.c.",80,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",85,".",paste0(names.columns)), paste0("cov.MAR.c.",90,".",paste0(names.columns)),
                               paste0("cov.MAR.c.",95,".",paste0(names.columns)))
    
    
    
    # Values
    
    outputvector.values <- c(epidemic.metrics, epi.behav.stats, 
                             results.mcar, results.mar.a, 
                             results.mar.b, results.mar.c)

    
  
  return(outputvector.values)
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  

  
  
}


