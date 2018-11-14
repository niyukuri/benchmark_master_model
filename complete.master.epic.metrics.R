



complete.master.epic.metrics <- function(datalist = datalist.agemix){
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")

  
  datalist.agemix <- datalist
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  
  ##################################################
  ### Compute transmission network characteristics #
  ##################################################
  
  
  
  # 1. Incidence trend #
  ######################
  
  # 35 - 36
  
  METRICS.incidence.df.15.24.int.35.36.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(15, 25),
                                                                    timewindow = c(35, 36),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.15.24.int.35.36.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(15, 25),
                                                                      timewindow = c(35, 36),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.25.39.int.35.36.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(25, 40),
                                                                    timewindow = c(35, 36),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.25.39.int.35.36.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(25, 40),
                                                                      timewindow = c(35, 36),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.40.49.int.35.36.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(40, 50),
                                                                    timewindow = c(35, 36),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.40.49.int.35.36.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(40, 50),
                                                                      timewindow = c(35, 36),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  # 36 - 37
  
  
  METRICS.incidence.df.15.24.int.36.37.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(15, 25),
                                                                    timewindow = c(36, 37),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.15.24.int.36.37.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(15, 25),
                                                                      timewindow = c(36, 37),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.25.39.int.36.37.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(25, 40),
                                                                    timewindow = c(36, 37),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.25.39.int.36.37.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(25, 40),
                                                                      timewindow = c(36, 37),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.40.49.int.36.37.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(40, 50),
                                                                    timewindow = c(36, 37),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.40.49.int.36.37.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(40, 50),
                                                                      timewindow = c(36, 37),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  # 37 - 38
  
  
  METRICS.incidence.df.15.24.int.37.38.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(15, 25),
                                                                    timewindow = c(37, 38),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.15.24.int.37.38.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(15, 25),
                                                                      timewindow = c(37, 38),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.25.39.int.37.38.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(25, 40),
                                                                    timewindow = c(37, 38),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.25.39.int.37.38.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(25, 40),
                                                                      timewindow = c(37, 38),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.40.49.int.37.38.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(40, 50),
                                                                    timewindow = c(37, 38),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.40.49.int.37.38.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(40, 50),
                                                                      timewindow = c(37, 38),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  # 38 - 39
  
  
  METRICS.incidence.df.15.24.int.38.39.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(15, 25),
                                                                    timewindow = c(38, 39),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.15.24.int.38.39.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(15, 25),
                                                                      timewindow = c(38, 39),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.25.39.int.38.39.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(25, 40),
                                                                    timewindow = c(38, 39),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.25.39.int.38.39.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(25, 40),
                                                                      timewindow = c(38, 39),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.40.49.int.38.39.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(40, 50),
                                                                    timewindow = c(38, 39),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.40.49.int.38.39.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(40, 50),
                                                                      timewindow = c(38, 39),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  # 39 - 40
  
  
  METRICS.incidence.df.15.24.int.39.40.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(15, 25),
                                                                    timewindow = c(39, 40),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.15.24.int.39.40.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(15, 25),
                                                                      timewindow = c(39, 40),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.25.39.int.39.40.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(25, 40),
                                                                    timewindow = c(39, 40),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.25.39.int.39.40.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(25, 40),
                                                                      timewindow = c(39, 40),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  METRICS.incidence.df.40.49.int.39.40.men <-  incidence.calculator(datalist = datalist.agemix,
                                                                    agegroup = c(40, 50),
                                                                    timewindow = c(39, 40),
                                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  METRICS.incidence.df.40.49.int.39.40.women <-  incidence.calculator(datalist = datalist.agemix,
                                                                      agegroup = c(40, 50),
                                                                      timewindow = c(39, 40),
                                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  # 
  
  
  
  # 2. Age mixing #
  ################
  
  # (i) Age mixing in relationships
  
  # 
  
  agemix.rels.df <- agemix.df.maker(datalist.agemix)
  
  # 
  agemix.model <- pattern.modeller(dataframe = agemix.rels.df,
                                   agegroup = c(15, 50),
                                   timepoint = 40, # datalist.agemix$itable$population.simtime[1],
                                   timewindow = 10)#1)#3)
  # 
  # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  #
  # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
  
  data = dplyr::filter(agemix.model[[1]], Gender =="male")
  
  if( nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
    
    men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                     data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                     REML = TRUE,
                     control=lmerControl(check.nobs.vs.nlev = "ignore",
                                         check.nobs.vs.rankZ = "ignore",
                                         check.nobs.vs.nRE="ignore"))
    
    bignumber <- NA # let's try if NA works (instead of 9999 for example)
    AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
    slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
    WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
    
    BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
    
    intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
    
    # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
    ## AAD: average age difference across all relationship
    ## VAD: variance of these age differences
    ## SDAD: standard deviation of age differences
    ## BSD: between-subject standard deviation of age differences
    
    mix.rels.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
    names(mix.rels.dat) <- c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
    names(mix.rels.dat) <- c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")
    
  }
  
  
  METRICS.LMEM.rels.age.mix <-  mix.rels.dat 
  
  
  

  # 3. Onward transmissions #
  ###########################
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = simpact.trans.net)
  
  METRICS.transm.average <- mean(transm.count)
  
  METRICS.transm.median <- median(transm.count) # add
  
  METRICS.transm.sd <- sd(transm.count) # add
  
  
  
  epi.Metrics <- c(METRICS.incidence.df.15.24.int.35.36.men, METRICS.incidence.df.15.24.int.35.36.women,
                   METRICS.incidence.df.15.24.int.36.37.men, METRICS.incidence.df.15.24.int.36.37.women,
                   METRICS.incidence.df.15.24.int.37.38.men, METRICS.incidence.df.15.24.int.37.38.women,
                   METRICS.incidence.df.15.24.int.38.39.men, METRICS.incidence.df.15.24.int.38.39.women,
                   METRICS.incidence.df.15.24.int.39.40.men, METRICS.incidence.df.15.24.int.39.40.women,
                   
                   METRICS.incidence.df.25.39.int.35.36.men, METRICS.incidence.df.25.39.int.35.36.women,
                   METRICS.incidence.df.25.39.int.36.37.men, METRICS.incidence.df.25.39.int.36.37.women,
                   METRICS.incidence.df.25.39.int.37.38.men, METRICS.incidence.df.25.39.int.37.38.women,
                   METRICS.incidence.df.25.39.int.38.39.men, METRICS.incidence.df.25.39.int.38.39.women,
                   METRICS.incidence.df.25.39.int.39.40.men, METRICS.incidence.df.25.39.int.39.40.women,
                   
                   METRICS.incidence.df.40.49.int.35.36.men, METRICS.incidence.df.40.49.int.35.36.women,
                   METRICS.incidence.df.40.49.int.36.37.men, METRICS.incidence.df.40.49.int.36.37.women,
                   METRICS.incidence.df.40.49.int.37.38.men, METRICS.incidence.df.40.49.int.37.38.women,
                   METRICS.incidence.df.40.49.int.38.39.men, METRICS.incidence.df.40.49.int.38.39.women,
                   METRICS.incidence.df.40.49.int.39.40.men, METRICS.incidence.df.40.49.int.39.40.women,
                   
                   METRICS.LMEM.rels.age.mix,
                   
                   METRICS.transm.average, METRICS.transm.median, METRICS.transm.sd)
  
  
  
  metric.names <- c("metr.incid.15.24.int.35.36.m", "metr.incid.15.24.int.35.36.w",
                    "metr.incid.15.24.int.36.37.m", "metr.incid.15.24.int.36.37.w",
                    "metr.incid.15.24.int.37.38.m", "metr.incid.15.24.int.37.38.w",
                    "metr.incid.15.24.int.38.39.m", "metr.incid.15.24.int.38.39.w",
                    "metr.incid.15.24.int.39.40.m", "metr.incid.15.24.int.39.40.w",
                    
                    "metr.incid.25.39.int.35.36.m", "metr.incid.25.39.int.35.36.w",
                    "metr.incid.25.39.int.36.37.m", "metr.incid.25.39.int.36.37.w",
                    "metr.incid.25.39.int.37.38.m", "metr.incid.25.39.int.37.38.w",
                    "metr.incid.25.39.int.38.39.m", "metr.incid.25.39.int.38.39.w",
                    "metr.incid.25.39.int.39.40.m", "metr.incid.25.39.int.39.40.w",
                    
                    "metr.incid.40.49.int.35.36.m", "metr.incid.40.49.int.35.36.w",
                    "metr.incid.40.49.int.36.37.m", "metr.incid.40.49.int.36.37.w",
                    "metr.incid.40.49.int.37.38.m", "metr.incid.40.49.int.37.38.w",
                    "metr.incid.40.49.int.38.39.m", "metr.incid.40.49.int.38.39.w",
                    "metr.incid.40.49.int.39.40.m", "metr.incid.40.49.int.39.40.w",
                    
                    paste0("metr.",names(METRICS.LMEM.rels.age.mix)), 
                    
                    "metr.transm.av",
                    "metr.transm.med", "metr.transm.sd")
  
  
  names(epi.Metrics) <- metric.names
  
  return(epi.Metrics)
  
  
}

