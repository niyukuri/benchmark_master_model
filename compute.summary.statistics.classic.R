



compute.summary.statistics.classic <- function(datalist = datalist.agemix,
                                               timewindow = c(30, 40)){
  
  
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  
  datalist.agemix <- datalist
  
  
  ########################################
  # I. Behavioural and epidemic features #
  ########################################
  
  
  # 1.2. Features from sexual and transmission network
  
  # 
  # 1.2.1. Demographic feature:
  
  #   (i) Population growth rate (pop.growth.calculator function)
  
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = timewindow) # c(0, datalist.agemix$itable$population.simtime[1])
  
  
  # 1.2.2. Prevalence
  
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.25.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 40),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.25.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 40),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.40.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(40, 50),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.40.50.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(40, 50),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  
  # Incidence
  
  incid.15.24.men <-  incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(15, 25),
                                           timewindow = c(39, 40),
                                           only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  incid.15.24.women <-  incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timewindow = c(39, 40),
                                             only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  incid.25.39.men <-  incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(25, 40),
                                           timewindow = c(39, 40),
                                           only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  incid.25.39.women <-  incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(25, 40),
                                             timewindow = c(39, 40),
                                             only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  incid.40.49.men <-  incidence.calculator(datalist = datalist.agemix,
                                           agegroup = c(40, 50),
                                           timewindow = c(39, 40),
                                           only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  incid.40.49.women <-  incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(40, 50),
                                             timewindow = c(39, 40),
                                             only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  # 
  
  
  
  # Sexual behaviour
  
  
  #  Point 	prevalence of concurrency in the adult population
  
  # Concurrency point prevalence 6 months before a survey, among men
  
  pp.cp.6months.male.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                          timepoint = timewindow[2] - 0.5) 
  
  
  
  # pp.cp.6months.male.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
  #                                                         timepoint = 40 - 0.5) %>%
  #   dplyr::select(concurr.pointprev) %>%
  #   dplyr::slice(1) %>%
  #   as.numeric()
  # 
  # pp.cp.6months.female.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
  #                                                           timepoint = 40 - 0.5) %>%
  #   dplyr::select(concurr.pointprev) %>%
  #   dplyr::slice(2) %>%
  #   as.numeric()
  # 
  
  
  # (ii) Relationship per person per year ??
  
  relsperpersonperyear <- nrow(datalist.agemix$rtable) / (nrow(datalist.agemix$ptable)/2) / (timewindow[2] - timewindow[1])
  
  # (iv) SD age gap between couples
  
  agegap.mean <- mean(datalist.agemix$rtable$AgeGap)
  
  agegap.med <- median(datalist.agemix$rtable$AgeGap)
  
  agegap.sd <- sd(datalist.agemix$rtable$AgeGap)
  
  
  
  ####
  # ART coverage among adults 15+ years old from UNAIDS (2010 - 2017 estimates)
  ####
  ART.cov.eval.timepoints <- seq(from = 33,
                                 to = 40)
  
  ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  
  for (art.cov.index in 1:length(ART.cov.vector)){
    ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist.agemix,
                                                                 agegroup = c(15, 50),
                                                                 timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
      sum(ART.coverage.calculator(datalist = datalist.agemix,
                                  agegroup = c(15, 50),
                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  }
  
  names(ART.cov.vector) <- paste0("ART.", ART.cov.eval.timepoints)
  
  ####
  # VL suppression fraction (all ages in 2017 ~ >= 15 yo) 0.74
  ####
  VL.suppression.fraction <- VL.suppression.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 50),
                                                       timepoint = 40,
                                                       vl.cutoff = 1000,
                                                       site="All") %>%
    dplyr::select(vl.suppr.frac) %>%
    dplyr::slice(3) %>%
    as.numeric()
  
  names(VL.suppression.fraction) <- "VL.suppr." 
  
  
  classic.features <-   c(exp(growthrate), 
                          
                          hiv.prev.lt25.men, hiv.prev.lt25.women, 
                          hiv.prev.25.40.men, hiv.prev.25.40.women,
                          hiv.prev.40.50.men, hiv.prev.40.50.women, 
                          
                          exp(incid.15.24.men), exp(incid.15.24.women), 
                          exp(incid.25.39.men), exp(incid.25.39.women), 
                          exp(incid.40.49.men), exp(incid.40.49.women),
                          
                          pp.cp.6months.male.rels, # pp.cp.6months.female.rels,
                          
                          relsperpersonperyear, 
                          agegap.mean, agegap.med, agegap.sd,
                          
                          ART.cov.vector, VL.suppression.fraction)
  
  
  
  
  classic.features.names <- c("Pop.growthrate", 
                              
                              "prev.15.25.men", "prev.15.25.women", 
                              "prev.25.40.men", "prev.25.40.women", 
                              "prev.40.50.men", "prev.40.50.women", 
                              
                              "incid.15.24.men", "incid.15.24.women", 
                              "incid.25.39.men", "incid.25.39.women", 
                              "incid.40.49.men", "incid.40.49.women",
                              
                              "pp.cp.6months.male.rels",
                              
                              "relsperpersonperyear", 
                              "agegap.mean", "agegap.med", "agegap.sd", 
                              
                              paste0(names(ART.cov.vector)), paste0(names(VL.suppression.fraction)))
  
  names(classic.features) <- classic.features.names 
  
  
  return(classic.features)
  
}


