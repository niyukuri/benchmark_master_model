
complete.master.epic.metric.class.phylo.features.cov <- function(inputvector = inputvector){
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/age_mix_final/advanced.transmission.network.builder.R")
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/ALL.funs.R")
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metrics.R")
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.classic.R")
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/mixed.effect.fit.transmission.clusters.R")
  # 
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/stats.age.groups.trans.clust.network.fun.R")
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MCAR.R")
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.phylo.MAR.R")
  
  
  source("/home/dniyukuri/lustre/benchmark_master_model/advanced.transmission.network.builder.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epic.metrics.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.classic.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  
  
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  
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
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  
  
  
  
  ## Run Simpact for specific parameter combination
  
  age.distr <- agedistr.creator(shape = 5, scale = 65)
  #
  cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                   population.simtime = 40, 
                                   population.nummen = 8000, 
                                   population.numwomen = 8000,
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
  
  
  
  sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
  
  
  
  
  # Running Simpact 
  #################
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  if (length(results) == 0){
    outputvector <- rep(NA, 2027)
  }else{
    if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)){
      outputvector <- rep(NA, 2027)
    }else{
      
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
      
      if(length(big.index) >= 1){ 
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                error=function(e) return(rep(NA, 37)))
        
        MCAR.All <- c(MCAR.cov.35, MCAR.cov.40, MCAR.cov.45, MCAR.cov.50, MCAR.cov.55, MCAR.cov.60, MCAR.cov.65, MCAR.cov.70, 
                      MCAR.cov.75, MCAR.cov.80, MCAR.cov.85, MCAR.cov.90, MCAR.cov.95, MCAR.cov.100)
        
        
        
        
        
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
        MAR.a.All <- c(MAR.a.cov.35, MAR.a.cov.40, MAR.a.cov.45, MAR.a.cov.50, MAR.a.cov.55, MAR.a.cov.60, MAR.a.cov.65, MAR.a.cov.70, 
                       MAR.a.cov.75, MAR.a.cov.80, MAR.a.cov.85, MAR.a.cov.90, MAR.a.cov.95)
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
        MAR.b.All <- c(MAR.b.cov.35, MAR.b.cov.40, MAR.b.cov.45, MAR.b.cov.50, MAR.b.cov.55, MAR.b.cov.60, MAR.b.cov.65, MAR.b.cov.70, 
                       MAR.b.cov.75, MAR.b.cov.80, MAR.b.cov.85, MAR.b.cov.90, MAR.b.cov.95)
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
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
                                                                      timewindow = c(30,40),
                                                                      cut.off = 7),
                                 error=function(e) return(rep(NA, 37)))
        
        
        MAR.c.All <- c(MAR.c.cov.35, MAR.c.cov.40, MAR.c.cov.45, MAR.c.cov.50, MAR.c.cov.55, MAR.c.cov.60, MAR.c.cov.65, MAR.c.cov.70, 
                       MAR.c.cov.75, MAR.c.cov.80, MAR.c.cov.85, MAR.c.cov.90, MAR.c.cov.95)
        
        
        # Values
        
        outputvector.values <- as.numeric(c(epidemic.metrics, epi.behav.stats, MCAR.All, MAR.a.All, MAR.b.All, MAR.c.All))
 
        
        
        # Names
        
        name.epi.metrics <- names(epidemic.metrics)
        name.epi.behav.stats <- names(epi.behav.stats)
        
        name.MCAR.35 <- paste0("MCAR.35.", names(MCAR.cov.35))
        name.MCAR.40 <- paste0("MCAR.40.", names(MCAR.cov.40))
        name.MCAR.45 <- paste0("MCAR.45.", names(MCAR.cov.45))
        name.MCAR.50 <- paste0("MCAR.50.", names(MCAR.cov.50))
        name.MCAR.55 <- paste0("MCAR.55.", names(MCAR.cov.55))
        name.MCAR.60 <- paste0("MCAR.60.", names(MCAR.cov.60))
        name.MCAR.65 <- paste0("MCAR.65.", names(MCAR.cov.65))
        name.MCAR.70 <- paste0("MCAR.70.", names(MCAR.cov.70))
        name.MCAR.75 <- paste0("MCAR.75.", names(MCAR.cov.75))
        name.MCAR.80 <- paste0("MCAR.80.", names(MCAR.cov.80))
        name.MCAR.85 <- paste0("MCAR.85.", names(MCAR.cov.85))
        name.MCAR.90 <- paste0("MCAR.90.", names(MCAR.cov.90))
        name.MCAR.95 <- paste0("MCAR.95.", names(MCAR.cov.95))
        name.MCAR.100 <- paste0("MCAR.100.", names(MCAR.cov.100))
        
        name.MCAR.scenari <- c(name.MCAR.35, name.MCAR.40, name.MCAR.45, name.MCAR.50,
                               name.MCAR.55, name.MCAR.60, name.MCAR.65, name.MCAR.70,
                               name.MCAR.75, name.MCAR.80, name.MCAR.85,
                               name.MCAR.90, name.MCAR.95, name.MCAR.100)
        
        name.MAR.a.35 <- paste0("MAR.a.35.", names(MAR.a.cov.35))
        name.MAR.a.40 <- paste0("MAR.a.40.", names(MAR.a.cov.40))
        name.MAR.a.45 <- paste0("MAR.a.45.", names(MAR.a.cov.45))
        name.MAR.a.50 <- paste0("MAR.a.50.", names(MAR.a.cov.50))
        name.MAR.a.55 <- paste0("MAR.a.55.", names(MAR.a.cov.55))
        name.MAR.a.60 <- paste0("MAR.a.60.", names(MAR.a.cov.60))
        name.MAR.a.65 <- paste0("MAR.a.65.", names(MAR.a.cov.65))
        name.MAR.a.70 <- paste0("MAR.a.70.", names(MAR.a.cov.70))
        name.MAR.a.75 <- paste0("MAR.a.75.", names(MAR.a.cov.75))
        name.MAR.a.80 <- paste0("MAR.a.80.", names(MAR.a.cov.80))
        name.MAR.a.85 <- paste0("MAR.a.85.", names(MAR.a.cov.85))
        name.MAR.a.90 <- paste0("MAR.a.90.", names(MAR.a.cov.90))
        name.MAR.a.95 <- paste0("MAR.a.95.", names(MAR.a.cov.95))
        
        
        name.MAR.a.scenari <- c(name.MAR.a.35, name.MAR.a.40, name.MAR.a.45, name.MAR.a.50,
                                name.MAR.a.55, name.MAR.a.60, name.MAR.a.65, name.MAR.a.70,
                                name.MAR.a.75, name.MAR.a.80, name.MAR.a.85,
                                name.MAR.a.90, name.MAR.a.95)
        
        
        name.MAR.b.35 <- paste0("MAR.b.35.", names(MAR.b.cov.35))
        name.MAR.b.40 <- paste0("MAR.b.40.", names(MAR.b.cov.40))
        name.MAR.b.45 <- paste0("MAR.b.45.", names(MAR.b.cov.45))
        name.MAR.b.50 <- paste0("MAR.b.50.", names(MAR.b.cov.50))
        name.MAR.b.55 <- paste0("MAR.b.55.", names(MAR.b.cov.55))
        name.MAR.b.60 <- paste0("MAR.b.60.", names(MAR.b.cov.60))
        name.MAR.b.65 <- paste0("MAR.b.65.", names(MAR.b.cov.65))
        name.MAR.b.70 <- paste0("MAR.b.70.", names(MAR.b.cov.70))
        name.MAR.b.75 <- paste0("MAR.b.75.", names(MAR.b.cov.75))
        name.MAR.b.80 <- paste0("MAR.b.80.", names(MAR.b.cov.80))
        name.MAR.b.85 <- paste0("MAR.b.85.", names(MAR.b.cov.85))
        name.MAR.b.90 <- paste0("MAR.b.90.", names(MAR.b.cov.90))
        name.MAR.b.95 <- paste0("MAR.b.95.", names(MAR.b.cov.95))
        
        
        name.MAR.b.scenari <- c(name.MAR.b.35, name.MAR.b.40, name.MAR.b.45, name.MAR.b.50,
                                name.MAR.b.55, name.MAR.b.60, name.MAR.b.65, name.MAR.b.70,
                                name.MAR.b.75, name.MAR.b.80, name.MAR.b.85,
                                name.MAR.b.90, name.MAR.b.95)
        
        name.MAR.c.35 <- paste0("MAR.c.35.", names(MAR.c.cov.35))
        name.MAR.c.40 <- paste0("MAR.c.40.", names(MAR.c.cov.40))
        name.MAR.c.45 <- paste0("MAR.c.45.", names(MAR.c.cov.45))
        name.MAR.c.50 <- paste0("MAR.c.50.", names(MAR.c.cov.50))
        name.MAR.c.55 <- paste0("MAR.c.55.", names(MAR.c.cov.55))
        name.MAR.c.60 <- paste0("MAR.c.60.", names(MAR.c.cov.60))
        name.MAR.c.65 <- paste0("MAR.c.65.", names(MAR.c.cov.65))
        name.MAR.c.70 <- paste0("MAR.c.70.", names(MAR.c.cov.70))
        name.MAR.c.75 <- paste0("MAR.c.75.", names(MAR.c.cov.75))
        name.MAR.c.80 <- paste0("MAR.c.80.", names(MAR.c.cov.80))
        name.MAR.c.85 <- paste0("MAR.c.85.", names(MAR.c.cov.85))
        name.MAR.c.90 <- paste0("MAR.c.90.", names(MAR.c.cov.90))
        name.MAR.c.95 <- paste0("MAR.c.95.", names(MAR.c.cov.95))
        
        
        name.MAR.c.scenari <- c(name.MAR.c.35, name.MAR.c.40, name.MAR.c.45, name.MAR.c.50,
                                name.MAR.c.55, name.MAR.c.60, name.MAR.c.65, name.MAR.c.70,
                                name.MAR.c.75, name.MAR.c.80, name.MAR.c.85,
                                name.MAR.c.90, name.MAR.c.95)
        
        outputvector.names <- c(name.epi.metrics, name.epi.behav.stats, name.MCAR.scenari, name.MAR.a.scenari, name.MAR.b.scenari, name.MAR.c.scenari)
        
        
        
        names(outputvector.values) <- outputvector.names
        
        
      }else{
        
        outputvector.values <- rep(NA, 2027)
        
      }
    }
    
  }
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  
  
  return(outputvector.values)
  
  
}


