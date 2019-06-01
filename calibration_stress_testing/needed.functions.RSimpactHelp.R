


bvar <- function(model) {
  
  # Outputs a df with between-subject variance, upr & lwr limits
  
  # Must take an merMod  object
  
  bsd <- as.numeric(as.data.frame(VarCorr(model))[1,5])
  
  
}


# Calculate the transmission rate
###################################


transmission.rate.calculator <- function(datalist = datalist,
                                         timewindow = c(20, 40), int = FALSE, by=1){
  
  if(int==FALSE){
    
    Infec.pop.table <- datalist$ptable[InfectType==1]
    
    numb.infec.pop <- nrow(Infec.pop.table %>%
                             subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))
    
    
    transm.rate <- numb.infec.pop / diff(timewindow)
    
    return(transm.rate)
  }
  
  if(int==TRUE){
    
    Infec.pop.table <- datalist$ptable[InfectType==1]
    
    upper.limit <- ceiling(diff(timewindow)/by)
    
    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
    
    infec.pop.int <- vector()
    trans.rate.int <- vector()
    
    
    for(i in 0:(upper.limit-2)){
      
      timewindow.int <- c(interv.time[1+i], interv.time[2+i])
      
      infec.pop.numb <- nrow(Infec.pop.table %>%
                               subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))
      
      
      infec.pop.int <- c(infec.pop.int,  infec.pop.numb)
      
      
      trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))
      
    }
    
    return(trans.rate.int)
  }
  
}


# Relationship rate
#####################


relationship.rate.calculator <- function(datalist = datalist,
                                         timewindow = c(20, 40), int = FALSE, by=1){
  
  if(int==FALSE){
    Rels.table <- datalist$rtable
    
    Rels.table.window <- Rels.table %>%
      subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])
    
    numb.rels <- nrow(Rels.table.window)
    
    numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women
    
    numb.rels.women <- length(unique(Rels.table.window$ID2))
    
    rels.rate <- (numb.rels)/ diff(timewindow)
    
    return(rels.rate)
  }
  
  if(int==TRUE){
    
    Rels.table <- datalist$rtable
    
    upper.limit <- ceiling(diff(timewindow)/by)
    
    interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)
    
    rels.int <- vector()
    rels.rate.int <- vector()
    
    
    for(i in 0:(upper.limit-2)){
      
      timewindow.int <- c(interv.time[1+i], interv.time[2+i])
      
      rels.numb <- nrow(Rels.table %>%
                          subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))
      
      
      rels.int <- c(rels.int,  rels.numb)
      
      
      rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))
      
    }
    
    return(rels.rate.int)
  }
  
  
}


# age mixing pattern model
###########################

ampmodel <- function(data = dplyr::filter(agemix.model[[1]], Gender =="male")) {
  lmer(pagerelform ~ agerelform0 + (1 | ID),
       data = data,
       REML = TRUE)
}


# Concurrency point prevalence 6 months before a survey, among men
####################################################################


concurr.pointprev.calculator <- function(datalist = datalist,
                                         timepoint = datalist$itable$population.simtime[1] - 0.5){
  
  #  output <- data.table()
  
  DTalive.infected <- alive.infected(datalist = datalist,
                                     timepoint = timepoint, site = "All") # First we only take the data of people who were alive at time_i
  
  agemix.df <- agemix.df.maker(datalist)
  
  degrees.df <- degree.df.maker(df = agemix.df,
                                agegroup = c(15, 50),
                                hivstatus = 2,
                                survey.time = 40, # timepoint,
                                window.width = 1,
                                gender.degree = "male",
                                only.new = FALSE)
  
  number.people.with.cps <- sum(degrees.df$Degree > 1)
  popsize <- nrow(DTalive.infected)
  concurr.pointprevalence <- number.people.with.cps / popsize
  
  return(concurr.pointprevalence)
}



# ART coverage
###############

ART.coverage.vector.creator <- function(datalist = datalist,
                                        agegroup = c(15, 50)){
  ART.cov.eval.timepoints <- seq(from = datalist$itable$t[2],
                                 to = datalist$itable$population.simtime[1])
  ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  for (art.cov.index in 1:length(ART.cov.vector)){
    ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist,
                                                                 agegroup = agegroup,
                                                                 timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
      sum(ART.coverage.calculator(datalist = datalist,
                                  agegroup = agegroup,
                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  }
  return(ART.cov.vector)
}


# Function to compute lineages through time with confidence intervals
# ######################################################################

plot.parboot.ltt.dat <- function (pbtd, t0 = NA, res = 100, ...)
{
  t1 <- max(pbtd$td$sts, na.rm = T)
  if (is.na(t0))
    t0 <- min(sapply(pbtd$trees, function(tr) tr$timeOf))
  times <- seq(t0, t1, l = res)
  ltt <- cbind(times = times, t(sapply(times, function(t) {
    c(pml = sum(pbtd$td$sts > t) - sum(pbtd$td$Ti > t), setNames(quantile(sapply(pbtd$trees,
                                                                                 function(tre) sum(tre$sts > t) - sum(tre$Ti > t)),
                                                                          probs = c(0.025, 0.5, 0.975)), c("lb", "median",
                                                                                                           "ub")))
  })))
  pl.df <- as.data.frame(ltt)
  return(pl.df)
  # p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
  #                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
  # p <- p + geom_path(aes(x = times, y = pml))
  # (p <- p + ylab("Lineages through time") + xlab("Time"))
}

# LTT <- plot.parboot.ltt.dat(pb)

# names(LTT)
# # [1] "times"  "pml"    "lb"     "median" "ub" 
# p <- ggplot(LTT) + geom_ribbon(aes(x = times, ymin = lb,
#                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
# p <- p + geom_path(aes(x = times, y = pml))
# (p <- p + ylab("Lineages through time") + xlab("Time"))


# Age mixing in transmission
#############################

agemixing.trans.df <- function(trans.network = trans.network,
                               limitTransmEvents = 7){
  
  # id of people who got infection by seed event: seeds.id
  seeds.id <- length(trans.network)
  
  ID.select <- vector() # ID of selected transmission network
  ID.select.count <- vector() # number of individuals in these networks
  
  for (i in 1: seeds.id) {
    
    
    trans.network.i <- as.data.frame(trans.network[[i]])
    
    if(nrow(trans.network.i)>=limitTransmEvents){
      
      
      ID.select <- c(ID.select, i)
      ID.select.count <- c(ID.select.count, nrow(trans.network.i))
      
    } # X if
    
  } # Y for
  
  if(length(ID.select)>=1){
    
    infectionTable <- vector("list", length(ID.select))
    
    for(j in 1:length(ID.select)){
      
      p <- ID.select[j]
      
      trans.network.i <- as.data.frame(trans.network[[p]])
      
      trans.network.i <- trans.network.i[-1,]
      
      trans.network.i$AgeInfecDon <- abs(trans.network.i$TOBDon) + trans.network.i$InfecTime
      trans.network.i$AgeInfecRec <- abs(trans.network.i$TOBRec) + trans.network.i$InfecTime
      
      id.lab <- paste0(p,".",trans.network.i$id,".C")
      
      trans.network.i$id.lab <- id.lab
      
      infectionTable[[p]] <- trans.network.i
    }
    
    
    infecttable <- rbindlist(infectionTable)
    
  }else{
    infecttable <- NULL
  }
  
  return(infecttable)
  
}


# Fit age mixing in transmission
################################

fit.agemix.trans.men <- function(datatable = agemix.df){
  
  datatable <- datatable
  
  men.lmer <- lmer(AgeInfecDon ~ AgeInfecRec + (1 | DonId),
                   data = dplyr::filter(datatable, GenderDon =="0"),
                   REML = TRUE,
                   control=lmerControl(check.nobs.vs.nlev = "ignore",
                                       check.nobs.vs.rankZ = "ignore",
                                       check.nobs.vs.nRE="ignore"))
  agemix.inter <- men.lmer
  
  return(agemix.inter)
  
}




fit.agemix.trans.women <- function(datatable = agemix.df){
  
  datatable <- datatable
  
  women.lmer <- lmer(AgeInfecDon ~ AgeInfecRec + (1 | DonId),
                     data = dplyr::filter(datatable, GenderDon =="1"),
                     REML = TRUE,
                     control=lmerControl(check.nobs.vs.nlev = "ignore",
                                         check.nobs.vs.rankZ = "ignore",
                                         check.nobs.vs.nRE="ignore"))
  agemix.inter.women <- women.lmer
  
  return(agemix.inter.women)
  
}



# Onward transmissions
######################

onwardtransmissions.dat <- function(datalist = datalist, 
                                    trans.network = trans.network,
                                    limitTransmEvents = 3,
                                    time.window=c(10,40)){
  
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec.raw.died <- pers.infec.raw[pers.infec.raw$TOD != "Inf", ] # consider only these who had full time to transmit before they die
  
  # Infected individuals who died within this time window
  pers.infec.died <- pers.infec.raw.died[which(pers.infec.raw.died$TOD >= time.window[1] & pers.infec.raw.died$TOD <= time.window[2]),]
  
  # Their IDs
  pers.infec.died.IDs <- pers.infec.died$ID
  
  # Onward transmissions in each transmission network
  
  onwardtransm <- vector("list", length(trans.network))
  
  for (j in 1: length(trans.network)) {
    
    trans.network.j <- as.data.frame(trans.network[[j]])
    
    trans.network.j <- trans.network.j[-1,] # remove the universal infector
    
    if(nrow(trans.network.j) >= 1){ # consider transmission networks with at least one onward transmission
      
      d.j <- table(trans.network.j$DonId) # in the transmission table, the number of times DonId appears is the number of Onward transmissions after acuiring the infection 
      
      IDs.names <- as.numeric(names(d.j))
      count.transm <- as.numeric(d.j)
      
      own.df <- as.data.frame(cbind(IDs.names, count.transm))
      
      onwardtransm[[j]] <- own.df
      
    }
    
  }
  
  onwardtransmissions <- rbindlist(onwardtransm) 
  
  # Data table of these who died and their onward transmission
  count.dat <- subset(onwardtransmissions, onwardtransmissions$IDs.names%in%pers.infec.died.IDs) 
  
  
  count.dat.vec <- count.dat$count.transm
  
  return(count.dat.vec) # count.dat = all infections - seeds which didn;t produce at least one transmission
  
}


# New transmissions

new.transmissions.dat <- function(datalist = datalist, 
                                  time.window=c(10,40)){
  
  
  pers.infec.raw <- as.data.frame(datalist$ptable[InfectType != -1])
  
  pers.infec.raw.infec <- pers.infec.raw[pers.infec.raw$InfectTime != "Inf", ] # consider only these who have been infected
  
  # New infected individuals who during  this time window
  pers.infec <- pers.infec.raw.infec[which(pers.infec.raw.infec$InfectTime >= time.window[1] & pers.infec.raw.infec$InfectTime <= time.window[2]),]
  
  # Their IDs
  pers.infec.IDs <- pers.infec$ID
  
  return(pers.infec.IDs)
}

# compute phylogenetic features





phylogenetic.features.fun <- function(tree.topo=tree,
                                      tree.calib.LTT = tree.calib.LTT){
  
  ########################################
  ### FEATURES FROM PHYLOGENETIC TREE ####
  ########################################
  
  
  # 1.3. Features from phylogenetic tree:
  
  # library(phytools)
  
  # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
  
  tree.cal <- read.tree(paste0(tree.topo))
  
  # tree <- read.tree("calibrated.tree.save.nwk")
  
  # Mean height of internal nodes
  
  H <- nodeHeights(tree.cal) # similar to node.depth.edgelength(tree)
  
  # It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
  # is represented twice and only twice. Thus, if we exclude the root node (zero height), 
  # we can just take the mean of H[,1].
  
  mean.feature <- mean(sort(H[,1])[3:nrow(H)]) # important
  
  # library(phyloTop)
  
  colless.feature <- colless.phylo(tree.cal, normalise = TRUE)
  
  sackin.feature <- sackin.phylo(tree.cal, normalise = TRUE)
  
  
  Depths <- getDepths(tree.cal) # depth of tips and nodes
  
  mean.tipsDepths.feature <- mean(Depths$tipDepths)
  
  mean.nodesDepths.feature <- mean(Depths$nodeDepths)
  
  maxHeight.feature <- maxHeight(tree.cal, normalise = TRUE)
  
  # Estimating confidence intervals for rates and dates using a parametric bootstrap
  pb <- parboot.treedater(tree.calib.LTT) # Lineage Through Time
  
  # Lineages through time
  LTT <- plot.parboot.ltt.dat(pb)
  
  lb.mean.feature <- mean(LTT$lb) # mean of low values of LTT
  lb.median.feature <- median(LTT$lb) # median of low values of LTT
  
  ub.mean.feature <- mean(LTT$ub) # mean of upper values of LTT
  ub.median.feature <- median(LTT$ub) # median of upper values of LTT
  
  median.mean.feature <- mean(LTT$median) # mean of medians of values of LTT
  median.median.feature <- median(LTT$median) # median of medians of values of LTT
  
  phylo.faetures.summary <- c(mean.feature, colless.feature, sackin.feature, mean.tipsDepths.feature, mean.nodesDepths.feature,
                              maxHeight.feature, lb.mean.feature, lb.median.feature, ub.mean.feature, ub.median.feature,
                              median.mean.feature, median.median.feature)
  
  features.names <- c("meanHeight.feature", "colless.feature", "sackin.feature", "mean.tipsDepths.feature", "mean.nodesDepths.feature",
                      "maxHeight.feature", "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
                      "LTT.median.mean.feature", "LTT.median.median.feature")
  
  names(phylo.faetures.summary) <- features.names
  
  return(phylo.faetures.summary)
  
}


colMedians <- function(matrixdata){
  
  names.col <- names(matrixdata)
  
  med.vector <- vector()
  
  for(i in 1:length(names.col)){
    
    med.i <- median(matrixdata[,i])
    
    med.vector <- c(med.vector, med.i)
    
  }
  
  names(med.vector) <- names.col
  
  return(med.vector)
  
}


