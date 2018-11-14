

compute.summary.statistics.phylo.MCAR <- function(simpact.trans.net = simpact.trans.net,
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
                                                  cut.off = 7){
  
  
  
  # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")

  
  dirseqgen <- work.dir
  
  
  # Function for linear mixed effect models in transmission clusters
  ###################################################################
  
  
  mCAr.IDs <- IDs.Seq.Random(simpact.trans.net = simpact.trans.net,
                             limitTransmEvents = limitTransmEvents,
                             timewindow = timewindow,
                             seq.cov = seq.cov,
                             age.limit = age.group.40.50[2])
  
  # nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 
  
  dirfasttree <- dirfasttree
  
  
  
  if(length(mCAr.IDs)>10){
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mCAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))
    
    
    mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                          sub.dir.rename = sub.dir.rename,
                                                          fasttree.tool = "FastTreeMP",
                                                          calendar.dates = "samplingtimes.all.csv",
                                                          simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                          count.start = 1977,
                                                          endsim = endpoint,
                                                          clust = TRUE)
    
    N <- node.age(mCAr.IDs.tree.calib)
    
    # Time to MRCA: internal nodes ages
    
    int.node.age <- N$Ti
    
    latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
    
    
    
    # Fitting internal nodes: Poisson regression
    
    # Number of internal nodes aggregated
    
    
    min.node <- round(min(int.node.age))
    max.node <- round(max(int.node.age))
    
    intervals <- max.node - min.node
    
    dt.node.age.dt <- int.node.age
    
    
    i.vec <- vector() # initialise time intervals
    int.node.vec <- vector() # initialize internal nodes
    
    for (i in 1:intervals) {
      inf <- min.node - 1 + i
      sup <- min.node+i
      i.vec <- c(i.vec, sup)
      int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
      int.node.vec <- c(int.node.vec,length(int.node.age.i))
    }
    
    
    node.df <- data.frame(i.vec, int.node.vec)
    
    output <-glm(formula = int.node.vec ~ i.vec, data = node.df, family = poisson)
    
    intercept.fit <- output$coefficients[[1]]
    
    slope.fit <- output$coefficients[[2]]
    
    names(intercept.fit) <-     "node.interc.fit"
    names(slope.fit) <- "node.slope.fit"
    
    
    
    mrca.v <- mrca(mCAr.IDs.tree.calib, full = FALSE) # MRCA ids
    
    
    sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv")) # sampling times
    
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))
    
    
    
    # run ClusterPicker
    
    system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
    
    # Read clusters' files
    
    dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                     all.files = FALSE,
                     full.names = FALSE, recursive = FALSE)
    
    
    
    
    # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/mixed.effect.fit.transmission.clusters.R")
    
    mixed.effect.fit.transmission.clusters <- function(clust.names=clust.names,
                                                       simpact.trans.net = simpact.trans.net,
                                                       limitTransmEvents = limitTransmEvents,
                                                       select.IDs = select.IDs){
      
      seeds.id <- length(simpact.trans.net)
      
      # Add age at sampling
      new.transm.tab <- vector("list", seeds.id)
      
      for(i in 1:seeds.id){
        
        transm.age.df.ic <- as.data.frame(simpact.trans.net[[i]])
        
        age.samp.Rec <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBRec
        age.samp.Don <- transm.age.df.ic$SampTime - transm.age.df.ic$TOBDon
        
        transm.age.i <- cbind(transm.age.df.ic, age.samp.Rec, age.samp.Don)
        
        new.transm.tab[[i]] <- transm.age.i
        
      }
      
      # id of people who got infection by seed event: seeds.id
      trans.network <- new.transm.tab
      
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
      
      
      infectionTable <- vector("list", length(ID.select))
      
      for(j in 1:length(ID.select)){
        
        p <- ID.select[j]
        
        trans.network.i <- as.data.frame(trans.network[[p]])
        
        trans.network.i <- trans.network.i[-1,]
        
        id.lab <- paste0(p,".",trans.network.i$id,".C")
        
        trans.network.i$id.lab <- id.lab
        
        infectionTable[[p]] <- trans.network.i
      }
      
      
      infecttable <- rbindlist(infectionTable)
      
      
      
      transm.df <- infecttable
      
      
      table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)
      
      
      Study.DataTable <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%select.IDs) 
      
      
      IDs.study <- Study.DataTable$RecId
      
      
      transm.datalist.agemix <- datalist.agemix # assign full data set new age mix data set
      
      # Transmission table of selected individuals
      table.simpact.trans.net.cov <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%select.IDs)
      
      # Person table of selected individuals
      transm.datalist.agemix$ptable <- dplyr::filter(transm.datalist.agemix$ptable, transm.datalist.agemix$ptable$ID%in%IDs.study)
      
      
      
      dirfasttree <- work.dir
      
      
      
      # Transmission clusters.
      
      d <- clust.names
      
      transmission.clust.list <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster
      
      
      # Binding all data tables of clusters (age, gender, number of cluster)
      
      clust.size <- vector() # size of each cluster
      
      for (i in 1:length(d)) {
        
        transm.df.cl.dat <- NULL
        
        clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
        
        clust.size <- c(clust.size, nrow(clus.read))
        
        transm.df.cl <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
        
        transm.df.cl.dat$age <- transm.df.cl$age.samp.Rec
        transm.df.cl.dat$gender <- transm.df.cl$GenderRec
        transm.df.cl.dat$clust.id <- as.factor(rep(i, nrow(transm.df.cl)))
        
        transmission.clust.list[[i]] <- as.data.frame(transm.df.cl.dat)
        
      }
      
      
      clust.table.df <- as.data.frame(do.call(rbind, transmission.clust.list)) # data.table & data.frame
      
      if(nrow(clust.table.df) >= length(d)){
        
        het.fit.lme.agemixing <- lme(age ~ gender, data = clust.table.df, random = ~ 1|clust.id,
                                     weights = varIdent( c("1" = 0.5), ~ 1 |gender))
        
        
        het.a <- coef(summary(het.fit.lme.agemixing))[1] # average age in transmission clusters
        
        het.beta <- coef(summary(het.fit.lme.agemixing))[2] # average age difference in transmission clusters: 
        # seen as bridge width which shows potential cross-generation transmission
        
        
        het.b1 <- as.numeric(VarCorr(het.fit.lme.agemixing)[3]) # between cluster variation
        
        het.b2 <- as.numeric(VarCorr(het.fit.lme.agemixing)[4]) # within cluster variation
        
        
        # SD for the two strata
        
        unique.val.strat <- unique(attributes(het.fit.lme.agemixing$modelStruct$varStruct)$weights)
        
        het.fit.lme.agemixing$modelStruct$varStruct
        
        # reference group: female == 1
        delta.female <- 1
        
        female.val <- unique.val.strat[1]
        male.val <- unique.val.strat[2]
        
        delta.male <- female.val/male.val # delta_ref_group / val
        
        SD.female <- as.numeric(VarCorr(het.fit.lme.agemixing)[4])
        SD.male <- delta.male * SD.female 
        
        
        clust.lme.val <- c(het.a, het.beta, het.b1, het.b2, SD.female, SD.male)
        
        
        Num.Clus <- length(d)
        
        av.Clust.size <- sum(clust.size)/length(d)
        
        mean.Num.Clus <- mean(clust.size)
        med.Num.Clus <- median(clust.size)
        SD.Num.Clus <- sd(clust.size)
        
        ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size, mean.Num.Clus, med.Num.Clus, SD.Num.Clus)
        
        names(ouptuvector.clust) <- c("ME.cl.av.age.male", "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var", "ME.cl.SD.female", "ME.cl.SD.male", 
                                      "Num.Clusters", "av.Clust.Size", "mean.Cl.size", "med.Cl.size", "SD.Cl.size")
        
      }else{
        
        ouptuvector.clust <- c(clust.lme.val, Num.Clus, av.Clust.size)
        
        clust.stat.names <-  c("ME.cl.av.age.male", "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var", "ME.cl.SD.female", "ME.cl.SD.male", 
                               "Num.Clusters", "av.Clust.Size", "mean.Cl.size", "med.Cl.size", "SD.Cl.size")
        
        ouptuvector.clust <- rep(NA, length(clust.stat.names))
        
        names(ouptuvector.clust) <- clust.stat.names
        
        
      }
      
      return(ouptuvector.clust)
      
      
    }
    
    clust.fit.params <- mixed.effect.fit.transmission.clusters(clust.names=dd, 
                                                               simpact.trans.net = simpact.trans.net,
                                                               limitTransmEvents = limitTransmEvents,
                                                               select.IDs = mCAr.IDs)
    
    
    # source("~/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/stats.age.groups.trans.clust.network.fun.R")
    
    
    
    
    
    
    stats.age.groups.trans.clust.network.fun <- function(simpact.trans.net = simpact.trans.net.adv, 
                                                         datalist.agemix = datalist.agemix,
                                                         work.dir = work.dir,  
                                                         dirfasttree = dirfasttree, 
                                                         sub.dir.rename = sub.dir.rename,
                                                         limitTransmEvents = 7,
                                                         timewindow = c(30,40),
                                                         seq.cov = 35,
                                                         select.IDs = select.IDs,
                                                         age.group.15.25 = c(15,25),
                                                         age.group.25.40 = c(25,40),
                                                         age.group.40.50 = c(40,50),
                                                         cut.off = 7){
      
      
      
      
     # source("~/phylosimpact_simulation_studies_2018/stress_testing/needed.functions.RSimpactHelp.R")
      
      
      source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
      
      
      # Data list of infected individuals
      
      
      
      mCAr.IDs <- select.IDs
      
      
      if(length(mCAr.IDs) >= 10){
        
        
        simpact.trans.net.adv <- simpact.trans.net
        
        
        # Transmission network table as from transmission networks for further steps
        ############################################################################
        
        
        infectionTable <- vector("list", length(simpact.trans.net.adv))
        
        for(j in 1:length(simpact.trans.net.adv)){
          
          p <- j
          
          trans.network.i <- as.data.frame(simpact.trans.net.adv[[p]])
          
          # trans.network.i <- trans.network.i[-1,]
          
          id.lab <- paste0(p,".",trans.network.i$id,".C")
          
          trans.network.i$id.lab <- id.lab
          trans.network.i$ageSampTimeRec <- trans.network.i$SampTime - trans.network.i$TOBRec
          
          infectionTable[[p]] <- trans.network.i
          
          
        }
        
        
        infecttable <- rbindlist(infectionTable)
        
        
        table.simpact.trans.net.adv <- infecttable # rbindlist(simpact.trans.net.adv)
        
        
        Study.DataTable <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs) 
        
        
        IDs.study <- Study.DataTable$RecId
        
        
        transm.datalist.agemix <- datalist.agemix # assign full data set new age mix data set
        
        # Transmission table of selected individuals
        table.simpact.trans.net.cov <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%mCAr.IDs)
        
        # Person table of selected individuals
        transm.datalist.agemix$ptable <- dplyr::filter(transm.datalist.agemix$ptable, transm.datalist.agemix$ptable$ID%in%IDs.study)
        
        
        
        dirfasttree <- work.dir
        
        
        
        # Select sequences from the pool of alignment
        ##############################################
        
        
        choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                            select.vec = mCAr.IDs,
                            name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")))
        
        
        # Build and calibrate the phylogenetic tree
        ############################################
        
        mCAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                              sub.dir.rename = sub.dir.rename,
                                                              fasttree.tool = "FastTreeMP",
                                                              calendar.dates = "samplingtimes.all.csv",
                                                              simseqfile = paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),
                                                              count.start = 1977,
                                                              endsim = 40,
                                                              clust = TRUE)
        
        
        
        N <- node.age(mCAr.IDs.tree.calib)
        
        # Time to MRCA: internal nodes ages
        
        int.node.age <- N$Ti
        
        
        latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date
        
        
        mrca.v <- mrca(mCAr.IDs.tree.calib, full = FALSE) # MRCA ids
        
        
        sampling.dates <- read.csv(paste0(sub.dir.rename,"/samplingtimes.all.csv")) # sampling times
        
        # 
        # tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.tree")))
        # 
        
        
        # Compute transmission clusters
        ###############################
        
        # run ClusterPicker
        
        system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), paste0(sub.dir.rename,"/", paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta")), paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta.nwk")),  paste0("0.9 0.9 0.045 2 gap"))))
        
        # Read clusters' files
        
        dd <- list.files(path = paste0(sub.dir.rename), pattern = paste0(paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_",paste0("cov.",seq.cov, ".mCAr.IDs.C.Epidemic.Fasta"),"_","clusterPicks_cluste"),
                         all.files = FALSE,
                         full.names = FALSE, recursive = FALSE)
        
        # Transmission clusters.
        
        d <- clust.names <- dd
        
        data.list.simpact.trans.net.adv <-  vector("list", length(d)) # list() # initialise gender and age-structured data table of pairings in each transission cluster
        
        
        # Transmission table of individuals in the transmission clusters
        #################################################################
        
        
        # Binding all data tables of clusters as these information are captured in transmission networks
        
        clust.size <- vector() # size of each cluster # table.simpact.trans.net.adv
        
        transm.df <- table.simpact.trans.net.adv
        
        for (i in 1:length(d)) {
          
          transm.df.cl.dat <- NULL
          
          clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
          
          clust.size <- c(clust.size, nrow(clus.read))
          
          data.table.simpact.trans.net.i <- subset(transm.df, transm.df$id.lab%in%as.character(clus.read$V1)) # transmission data table of IDs of that cluster
          
          data.table.simpact.trans.net.i$clust.ID <- rep(i, nrow(data.table.simpact.trans.net.i))
          
          data.list.simpact.trans.net.adv[[i]] <- as.data.frame(data.table.simpact.trans.net.i)
          
        }
        
        
        data.table.simpact.trans.clusts.net.adv <- as.data.frame(do.call(rbind, data.list.simpact.trans.net.adv)) # data.table & data.frame
        
        data.table.simpact.trans.clusts.net.adv <- data.table.simpact.trans.clusts.net.adv[!duplicated(data.table.simpact.trans.clusts.net.adv[c("id","id.lab")]),] # remove duplicate id.lab
        # t may happen that one seq.ID appear in more than one cluster
        
        # data.table.simpact.trans.clusts.net.adv <- data.table.simpact.trans.net.adv
        
        
        ## Aligning internal nodes IDs and their age: !they must get same length
        
        ancestor <- Ancestors(mCAr.IDs.tree.calib) # ancestors of each tips and internal node
        # All ancestors output are internal nodes
        
        ancestor.v <- vector()
        
        for(i in 1:length(ancestor)){
          
          k <- ancestor[[i]]
          
          ancestor.v <- c(ancestor.v, unique(k))
          
        }
        
        sort.int.ancestor <- unique(sort(ancestor.v))
        sort.int.node.age <- sort(int.node.age)
        
        
        tip.names <- names(mrca.v[1,])
        
        dates.tree.df <- dplyr::filter(sampling.dates, sampling.dates$V1%in%tip.names) # dates of these tips
        
        # rearrange dates in tips order as are displayed on the tree
        tip.names.f <- vector()
        dates.tree.dat <- vector()
        for(i in 1:nrow(dates.tree.df)){
          for(j in 1:length(tip.names)){
            if(tip.names[i] == dates.tree.df$V1[[j]]){
              tip.names.f <- c(tip.names.f, tip.names[i])
              dates.tree.dat <- c(dates.tree.dat, 1977+40-dates.tree.df$V2[[j]])
            }
          }
        }
        
        
        dates.tree.named <- dates.tree.dat
        names(dates.tree.named) <- tip.names.f
        
        
        # MRCA matrix
        #############
        
        # make mrca matrix diagonal 0 and other elements (internal nodes IDs) assign them the age of mrca
        
        mrca.v.age <- mrca.v
        
        for(i in 1:nrow(mrca.v.age)){
          for(j in 1:nrow(mrca.v.age)){
            
            if(i==j){
              mrca.v.age[i,j] <- 0
            }else{
              
              if(mrca.v[i,j] %in% sort.int.ancestor){
                
                p.index <- which(sort.int.ancestor == mrca.v[i,j])
                
                mrca.v.age[i,j] <-  sort.int.node.age[p.index]
              }
              
            }
            
          }
        }
        
        
        
        # make mrca matrix elements: sampling date - age of mrca
        
        # Fist contingency matrix
        
        mrca.v.age.samp <- mrca.v.age
        
        mrca.v.age.samp.cont1 <- mrca.v.age.samp
        
        for(i in 1:nrow(mrca.v.age)){
          
          for(j in 1:nrow(mrca.v.age)){
            
            if(i!=j){
              
              i.dat <- tip.names.f[i]
              
              v.index <- which(tip.names.f == i.dat)
              
              samp.date.tip <- dates.tree.dat[v.index]
              
              mrca.v.age.samp.cont1[i,] <- samp.date.tip - mrca.v.age.samp[i,]
              
            }
            
          }
        }
        
        
        
        # Second contingency matrix
        
        mrca.v.age.samp <- mrca.v.age
        
        mrca.v.age.samp.cont2 <- mrca.v.age.samp
        
        for(i in 1:nrow(mrca.v.age)){
          
          for(j in 1:nrow(mrca.v.age)){
            
            if(i!=j){
              
              i.dat <- tip.names.f[i]
              
              v.index <- which(tip.names.f == i.dat)
              
              samp.date.tip <- dates.tree.dat[v.index]
              
              mrca.v.age.samp.cont2[,i] <- samp.date.tip - mrca.v.age.samp[,i]
              
            }
            
          }
        }
        
        
        # Diagonal zero for mrca.v.age.samp.cont1 and mrca.v.age.samp.cont2
        
        for(i in 1:nrow(mrca.v.age.samp.cont1)){
          for(j in 1:nrow(mrca.v.age.samp.cont1)){
            
            if(i==j){
              mrca.v.age.samp.cont1[i,j] <- 0
            }
          }
        }
        
        
        for(i in 1:nrow(mrca.v.age.samp.cont2)){
          for(j in 1:nrow(mrca.v.age.samp.cont2)){
            
            if(i==j){
              mrca.v.age.samp.cont2[i,j] <- 0
            }
          }
        }
        
        
        # filter table.simpact.trans.net.adv and remain with table of tips names (individulas in the tree)
        
        attributes.table.simpact.trans.net.adv <- dplyr::filter(table.simpact.trans.net.adv, table.simpact.trans.net.adv$id.lab%in%tip.names)
        
        V.gender <- vector()
        V.cd4 <- vector()
        V.vl <- vector()
        V.x <- vector()
        V.y <- vector()
        iD <- vector()
        
        for(i in 1:length(tip.names)){
          for(j in 1:nrow(attributes.table.simpact.trans.net.adv)){
            if(tip.names[i] == attributes.table.simpact.trans.net.adv$id.lab[j]){
              
              V.gender <- c(V.gender, attributes.table.simpact.trans.net.adv$GenderRec[j])
              V.cd4 <- c(V.cd4, attributes.table.simpact.trans.net.adv$cd4[j])
              V.vl <- c(V.vl, attributes.table.simpact.trans.net.adv$vl[j])
              V.x <- c(V.x, attributes.table.simpact.trans.net.adv$location.x[j])
              V.y <- c(V.y, attributes.table.simpact.trans.net.adv$location.y[j])
              iD <- c(iD, tip.names[i])
            }
            
          }
        }
        
        
        Node.gender.cd4.vl.x.y <- data.table(V.gender,V.cd4, V.vl, V.x, V.y, iD)
        
        
        # Adding clusters ID on the previous attributes table from attributes.table.simpact.trans.net.adv
        
        clust.ID.vec <- vector()
        id.vec <- vector()
        
        for(k in 1:nrow(Node.gender.cd4.vl.x.y)){ # attributes table ofr all tips on the tree: Node.gender.cd4.vl.x.y
          
          id <- Node.gender.cd4.vl.x.y$iD[k]
          
          
          if(id%in%data.table.simpact.trans.clusts.net.adv$id.lab){    # ID of tree which belongs to IDs of clusters
            # transmission table of individuls in the transmission clusters: data.table.simpact.trans.net.adv
            
            id.index <- which(data.table.simpact.trans.clusts.net.adv$id.lab == id)
            
            clust.ID.vec.i <- data.table.simpact.trans.clusts.net.adv$clust.ID[id.index]
            
          }else{
            
            clust.ID.vec.i <-  0 # tip ID which is not in any transmission cluster is assigned value 0
            
          }
          
          clust.ID.vec <-  c(clust.ID.vec, clust.ID.vec.i)
          id.vec <- c(id.vec, id)
          
        }
        
        Node.gender.cd4.vl.x.y$clust.ID <- clust.ID.vec
        
        Node.gender.cd4.vl.x.y.clusID <- Node.gender.cd4.vl.x.y # attributes table with clusters' IDs
        
        
        
        ## Building transmission network
        
        # 1. consider contigency matrix 2
        
        mrca.times.final <- as.matrix(abs(mrca.v.age.samp.cont2))
        
        
        net <- graph.adjacency(as.matrix(mrca.times.final), mode="undirected",weighted=T,diag=FALSE)
        
        # E(net)       # The edges of the "net" object
        # 
        # V(net)       # The vertices of the "net" object
        
        V(net)$gender <- Node.gender.cd4.vl.x.y$V.gender
        V(net)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
        V(net)$vl <- Node.gender.cd4.vl.x.y$V.vl
        V(net)$loc.x <- Node.gender.cd4.vl.x.y$V.x
        V(net)$loc.y <- Node.gender.cd4.vl.x.y$V.y
        
        
        
        
        ## Filtering the network by breaking some edges due to conditions from individuals attributes:
        ##############################################################################################
        
        # 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load
        
        # Now considering 1 and 2
        
        names.attributes.ngaha <- Node.gender.cd4.vl.x.y
        
        names.matrix.contigency <- names(mrca.times.final[1,])
        
        gender.l <- names.attributes.ngaha$V.gender
        
        
        clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID
        
        
        mrca.times.filter <- mrca.times.final
        
        
        # 
        # for (i in 1:length(names(mrca.times.final[1,]))) {
        #   
        #   name.col.i <- names.matrix.contigency[i]
        #   
        #   index.i <- which(names(mrca.times.final[1,]) == name.col.i)
        #   
        #   gender.i <- gender.l[index.i]
        #   
        #   cluster.i <- clusters.zose[index.i]
        #   
        #   
        #   for(j in 1:length(names(mrca.times.final[1,]))){
        #     
        #     if(i != j){
        #       
        #       name.col.j <- names.matrix.contigency[j]
        #       
        #       index.j <- which(names(mrca.times.final[1,]) == name.col.j)
        #       
        #       gender.j <- gender.l[index.j]
        #       
        #       cluster.j <- clusters.zose[index.j]
        #       
        #       
        #       if(gender.i == gender.j){ # if same gender break the link
        #         
        #         mrca.times.filter[i,j] <- 0
        #         
        #       }
        #       
        #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
        #         
        #         mrca.times.filter[i,j] <- 0
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        # }
        
        # i. Gender 
        ############
        
        for (i in 1:length(names(mrca.times.final[1,]))) {
          
          name.col.i <- names.matrix.contigency[i]
          
          index.i <- which(names(mrca.times.final[1,]) == name.col.i)
          
          gender.i <- gender.l[index.i]
          
          for(j in 1:length(names(mrca.times.final[1,]))){
            
            if(i != j){
              
              name.col.j <- names.matrix.contigency[j]
              
              index.j <- which(names(mrca.times.final[1,]) == name.col.j)
              
              gender.j <- gender.l[index.j]
              
              if(gender.i == gender.j){ # if same gender break the link
                
                mrca.times.filter[i,j] <- 0
                
              }
              
            }
            
          }
          
        }
        
        mrca.times.filter.gender <- mrca.times.filter
        
        
        # ii. Cluster
        #############
        
        mrca.times.filter.gender.clust <- mrca.times.filter.gender
        
        for (i in 1:length(names(mrca.times.final[1,]))) {
          
          name.col.i <- names.matrix.contigency[i]
          
          index.i <- which(names(mrca.times.final[1,]) == name.col.i)
          
          cluster.i <- clusters.zose[index.i]
          
          
          for(j in 1:length(names(mrca.times.final[1,]))){
            
            if(i != j){
              
              name.col.j <- names.matrix.contigency[j]
              
              index.j <- which(names(mrca.times.final[1,]) == name.col.j)
              
              cluster.j <- clusters.zose[index.j]
              
              
              if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
                
                mrca.times.filter.gender.clust[i,j] <- 0
                
              }
              
              
            }
            
          }
          
        }
        
        
        # iii. tMRCA
        #############
        
        
        
        net.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.gender.clust),mode="undirected",weighted=T,diag=FALSE)
        
        
        # Consider plausible transmissions and difference between sampling time and tMRCA
        
        
        cut.off <- cut.off # years
        
        # E(net.cont.1)$weight
        
        net.cont.1 <- delete_edges(net.cont.1, E(net.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff
        
        # E(net.cont.1)$weight
        
        # plot(net.cont.1, layout=layout_with_kk) 
        
        
        
        # Delete tips of the phylogenetic tree which are not part of transmission clusters:  they have clust.ID==0 >> deletes vertices 
        ###################################################################################
        
        
        Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
        Non.ids <- Non.ids.dat$iD
        
        net.cleaned <- delete_vertices(net.cont.1, Non.ids)
        
        
        
        # 
        # # 2. consider contigency matrix 1
        # 
        # mrca.times.final.2 <- as.matrix(abs(mrca.v.age.samp.cont1))
        # 
        # 
        # net.2 <- graph.adjacency(as.matrix(mrca.times.final.2), mode="undirected",weighted=T,diag=FALSE)
        # 
        # E(net.2)       # The edges of the "net.2" object
        # 
        # V(net.2)       # The vertices of the "net.2" object
        # 
        # V(net.2)$gender <- Node.gender.cd4.vl.x.y$V.gender
        # V(net.2)$cd4 <- Node.gender.cd4.vl.x.y$V.cd4
        # V(net.2)$vl <- Node.gender.cd4.vl.x.y$V.vl
        # V(net.2)$loc.x <- Node.gender.cd4.vl.x.y$V.x
        # V(net.2)$loc.y <- Node.gender.cd4.vl.x.y$V.y
        # 
        # 
        # 
        # 
        # ## Filtering the net.2work by breaking some edges due to conditions from individuals attributes:
        # 
        # # 1. Gender, 2. cluster belonging, 3. geographical location, 4. CD4, and 5. Viral load
        # 
        # 
        # names.attributes.ngaha <- Node.gender.cd4.vl.x.y
        # 
        # names.matrix.contigency <- names(mrca.times.final.2[1,])
        # 
        # gender.l <- names.attributes.ngaha$V.gender
        # 
        # 
        # clusters.zose <- Node.gender.cd4.vl.x.y$clust.ID
        # 
        # 
        # mrca.times.filter.2 <- mrca.times.final.2
        # 
        # 
        # # 
        # # for (i in 1:length(names(mrca.times.final.2[1,]))) {
        # #   
        # #   name.col.i <- names.matrix.contigency[i]
        # #   
        # #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
        # #   
        # #   gender.i <- gender.l[index.i]
        # #   
        # #   cluster.i <- clusters.zose[index.i]
        # #   
        # #   
        # #   for(j in 1:length(names(mrca.times.final.2[1,]))){
        # #     
        # #     if(i != j){
        # #       
        # #       name.col.j <- names.matrix.contigency[j]
        # #       
        # #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
        # #       
        # #       gender.j <- gender.l[index.j]
        # #       
        # #       cluster.j <- clusters.zose[index.j]
        # #       
        # #       
        # #       if(gender.i == gender.j){ # if same gender break the link
        # #         
        # #         mrca.times.filter.2[i,j] <- 0
        # #         
        # #       }
        # #       
        # #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
        # #         
        # #         mrca.times.filter.2[i,j] <- 0
        # #         
        # #       }
        # #       
        # #       
        # #     }
        # #     
        # #   }
        # #   
        # # }
        # 
        # # i. Gender 
        # 
        # for (i in 1:length(names(mrca.times.final.2[1,]))) {
        #   
        #   name.col.i <- names.matrix.contigency[i]
        #   
        #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
        #   
        #   gender.i <- gender.l[index.i]
        #   
        #   for(j in 1:length(names(mrca.times.final.2[1,]))){
        #     
        #     if(i != j){
        #       
        #       name.col.j <- names.matrix.contigency[j]
        #       
        #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
        #       
        #       gender.j <- gender.l[index.j]
        #       
        #       if(gender.i == gender.j){ # if same gender break the link
        #         
        #         mrca.times.filter.2[i,j] <- 0
        #         
        #       }
        #       
        #     }
        #     
        #   }
        #   
        # }
        # 
        # mrca.times.filter.2.gender <- mrca.times.filter.2
        # 
        # 
        # # ii. Cluster
        # 
        # mrca.times.filter.2.gender.clust <- mrca.times.filter.2.gender
        # 
        # for (i in 1:length(names(mrca.times.final.2[1,]))) {
        #   
        #   name.col.i <- names.matrix.contigency[i]
        #   
        #   index.i <- which(names(mrca.times.final.2[1,]) == name.col.i)
        #   
        #   cluster.i <- clusters.zose[index.i]
        #   
        #   
        #   for(j in 1:length(names(mrca.times.final.2[1,]))){
        #     
        #     if(i != j){
        #       
        #       name.col.j <- names.matrix.contigency[j]
        #       
        #       index.j <- which(names(mrca.times.final.2[1,]) == name.col.j)
        #       
        #       cluster.j <- clusters.zose[index.j]
        #       
        #       
        #       if(cluster.i != 0 & cluster.j != 0 & cluster.i != cluster.j){
        #         
        #         mrca.times.filter.2.gender.clust[i,j] <- 0
        #         
        #       }
        #       
        #       
        #     }
        #     
        #   }
        #   
        # }
        # 
        # 
        # 
        # net.2.cont.1 <- graph.adjacency(as.matrix(mrca.times.filter.2.gender.clust),mode="undirected",weighted=T,diag=FALSE)
        # 
        # 
        # # Consider plausible transmissions and difference between sampling time and tMRCA
        # 
        # 
        # cut.off <- 20
        # 
        # E(net.2.cont.1)$weight
        # 
        # net.2.cont.1 <- delete_edges(net.2.cont.1, E(net.2.cont.1)[weight>=cut.off]) # remove link greater to the cuttoff
        # 
        # E(net.2.cont.1)$weight
        # 
        # plot(net.2.cont.1, layout=layout_with_kk) 
        # 
        # 
        # # Delete tips which are not part of transmission clusters, they have clust.ID==0 >> deletes vertices 
        # 
        # Non.ids.dat <- dplyr::filter(Node.gender.cd4.vl.x.y, Node.gender.cd4.vl.x.y$clust.ID==0)
        # Non.ids <- Non.ids.dat$iD
        # 
        # net.2.cleaned <- delete_vertices(net.2.cont.1, Non.ids)
        
        
        
        # r=graph.union(net.cleaned, net.2.cleaned)
        
        
        
        
        # Age structure in the transmission network built from phylogenetic tree
        #########################################################################
        
        
        # produce age table
        
        net.sp <- net.cleaned
        
        
        transm.matrix <- as.data.table(get.edgelist(net.sp)) # matrix of links of the transmission network built from phylogenetic tree
        
        # table.simpact.trans.net.adv
        
        # reduced transmission table: table.simpact.trans.net.adv of ids in transmission clusters
        
        ids <-  unique(c(transm.matrix$V1, transm.matrix$V2))
        
        
        table.transm.clust.net.igraph <- dplyr::filter(data.table.simpact.trans.clusts.net.adv, data.table.simpact.trans.clusts.net.adv$id.lab%in%ids) 
        
        
        
        # 1.
        
        # Age structure in transmission clusters as observed from phylogenetic tree #
        ##################################################################################
        
        
        age.groups.filtered.trans.clust.network.fun <- function(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                transm.matrix = transm.matrix,
                                                                age.group.15.25 = c(15,25),
                                                                age.group.25.40 = c(25,40),
                                                                age.group.40.50 = c(40,50)){
          
          Age.groups.table <- NULL
          
          v1.dat <- vector()
          v2.dat <- vector()
          age1.dat <- vector()
          age2.dat <- vector()
          gender1.dat <- vector()
          gender2.dat <- vector()
          
          for(i in 1:nrow(transm.matrix)){
            
            v1 <- transm.matrix$V1[i]
            v2 <- transm.matrix$V2[i]
            
            index.v1 <- which(table.transm.clust.net.igraph$id.lab == v1)
            index.v2 <- which(table.transm.clust.net.igraph$id.lab == v2)
            
            age1 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v1]
            age2 <- table.transm.clust.net.igraph$ageSampTimeRec[index.v2]
            
            gender1 <- table.transm.clust.net.igraph$GenderRec[index.v1]
            gender2 <- table.transm.clust.net.igraph$GenderRec[index.v2]
            
            v1.dat <- c(v1.dat, v1)
            v2.dat <- c(v2.dat, v2)
            age1.dat <- c(age1.dat, age1)
            age2.dat <- c(age2.dat, age2)
            gender1.dat <- c(gender1.dat, gender1)
            gender2.dat <- c(gender2.dat, gender2)
            
          }
          
          age.table <- data.frame(v1.dat, gender1.dat, age1.dat, v2.dat, gender2.dat, age2.dat)
          
          
          
          # men
          men.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==0)
          
          # women
          women.age.table.1 <- dplyr::filter(age.table, age.table$gender1.dat==1)
          
          
          # men 15.25 and women
          
          men.15.25.women.15.25.1 <- vector()
          men.15.25.women.25.40.1 <- vector()
          men.15.25.women.40.50.1 <- vector()
          
          if(nrow(men.age.table.1)>1){
            
            for (j in 1:nrow(men.age.table.1)) {
              
              
              if(men.age.table.1$age1.dat[j] >= age.group.15.25[1] & men.age.table.1$age1.dat[j] < age.group.15.25[2]){
                
                if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  men.15.25.women.15.25.1 <- c(men.15.25.women.15.25.1, men.age.table.1$age2.dat[j])
                  
                }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  men.15.25.women.25.40.1 <- c(men.15.25.women.25.40.1, men.age.table.1$age2.dat[j])
                  
                }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  men.15.25.women.40.50.1 <- c(men.15.25.women.40.50.1, men.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
            
          }
          
          
          
          # women 15.25 and men
          
          women.15.25.men.15.25.2 <- vector()
          women.15.25.men.25.40.2 <- vector()
          women.15.25.men.40.50.2 <- vector()
          
          if(nrow(women.age.table.1)>1){
            
            for (j in 1:nrow(women.age.table.1)) {
              
              
              if(women.age.table.1$age1.dat[j] >= age.group.15.25[1] & women.age.table.1$age1.dat[j] < age.group.15.25[2]){
                
                if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  women.15.25.men.15.25.2 <- c(women.15.25.men.15.25.2, women.age.table.1$age2.dat[j])
                  
                }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  women.15.25.men.25.40.2 <- c(women.15.25.men.25.40.2, women.age.table.1$age2.dat[j])
                  
                }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  women.15.25.men.40.50.2 <- c(women.15.25.men.40.50.2, women.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
          }
          
          
          
          
          
          
          # men 25.40 and women
          
          men.25.40.women.15.25.1 <- vector()
          men.25.40.women.25.40.1 <- vector()
          men.25.40.women.40.50.1 <- vector()
          
          
          if(nrow(men.age.table.1) >1 ){
            for (j in 1:nrow(men.age.table.1)) {
              
              
              if(men.age.table.1$age1.dat[j] >= age.group.25.40[1] & men.age.table.1$age1.dat[j] < age.group.25.40[2]){
                
                if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  men.25.40.women.15.25.1 <- c(men.25.40.women.15.25.1, men.age.table.1$age2.dat[j])
                  
                }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  men.25.40.women.25.40.1 <- c(men.25.40.women.25.40.1, men.age.table.1$age2.dat[j])
                  
                }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  men.25.40.women.40.50.1 <- c(men.25.40.women.40.50.1, men.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
          }
          
          
          
          
          
          
          # women 25.40 and men
          
          women.25.40.men.15.25.2 <- vector()
          women.25.40.men.25.40.2 <- vector()
          women.25.40.men.40.50.2 <- vector()
          
          if(nrow(women.age.table.1) >1){
            
            for (j in 1:nrow(women.age.table.1)) {
              
              
              if(women.age.table.1$age1.dat[j] >= age.group.25.40[1] & women.age.table.1$age1.dat[j] < age.group.25.40[2]){
                
                if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  women.25.40.men.15.25.2 <- c(women.25.40.men.15.25.2, women.age.table.1$age2.dat[j])
                  
                }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  women.25.40.men.25.40.2 <- c(women.25.40.men.25.40.2, women.age.table.1$age2.dat[j])
                  
                }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  women.25.40.men.40.50.2 <- c(women.25.40.men.40.50.2, women.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
          }
          
          
          
          
          
          # men 40.50 and women
          
          men.40.50.women.15.25.1 <- vector()
          men.40.50.women.25.40.1 <- vector()
          men.40.50.women.40.50.1 <- vector()
          
          if(nrow(men.age.table.1) >1 ){
            
            for (j in 1:nrow(men.age.table.1)) {
              
              
              if(men.age.table.1$age1.dat[j] >= age.group.40.50[1] & men.age.table.1$age1.dat[j] < age.group.40.50[2]){
                
                if(men.age.table.1$age2.dat[j] >= age.group.15.25[1] & men.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  men.40.50.women.15.25.1 <- c(men.40.50.women.15.25.1, men.age.table.1$age2.dat[j])
                  
                }else if(men.age.table.1$age2.dat[j] >= age.group.25.40[1] & men.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  men.40.50.women.25.40.1 <- c(men.40.50.women.25.40.1, men.age.table.1$age2.dat[j])
                  
                }else if (men.age.table.1$age2.dat[j] >= age.group.40.50[1] & men.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  men.40.50.women.40.50.1 <- c(men.40.50.women.40.50.1, men.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
            
          }
          
          
          
          
          
          # women 40.50 and men
          
          women.40.50.men.15.25.2 <- vector()
          women.40.50.men.25.40.2 <- vector()
          women.40.50.men.40.50.2 <- vector()
          
          if(nrow(women.age.table.1) >1){
            
            for (j in 1:nrow(women.age.table.1)) {
              
              
              if(women.age.table.1$age1.dat[j] >= age.group.40.50[1] & women.age.table.1$age1.dat[j] < age.group.40.50[2]){
                
                if(women.age.table.1$age2.dat[j] >= age.group.15.25[1] & women.age.table.1$age2.dat[j] < age.group.15.25[2]){
                  
                  women.40.50.men.15.25.2 <- c(women.40.50.men.15.25.2, women.age.table.1$age2.dat[j])
                  
                }else if(women.age.table.1$age2.dat[j] >= age.group.25.40[1] & women.age.table.1$age2.dat[j] < age.group.25.40[2]){
                  
                  women.40.50.men.25.40.2 <- c(women.40.50.men.25.40.2, women.age.table.1$age2.dat[j])
                  
                }else if (women.age.table.1$age2.dat[j] >= age.group.40.50[1] & women.age.table.1$age2.dat[j] < age.group.40.50[2]){
                  
                  women.40.50.men.40.50.2 <- c(women.40.50.men.40.50.2, women.age.table.1$age2.dat[j])
                }
                
              }
              
              
            }
          }
          
          
          men.15.25.women.15.25 <- c(men.15.25.women.15.25.1, women.15.25.men.15.25.2)
          
          men.15.25.women.25.40 <- c(men.15.25.women.25.40.1, women.25.40.men.15.25.2)
          
          men.15.25.women.40.50 <- c(men.15.25.women.40.50.1, women.40.50.men.15.25.2)
          
          men.25.40.women.15.25 <- c(men.25.40.women.15.25.1, women.15.25.men.25.40.2)
          
          men.25.40.women.25.40 <- c(men.25.40.women.25.40.1, women.25.40.men.25.40.2)
          
          men.25.40.women.40.50 <- c(men.25.40.women.40.50.1, women.40.50.men.25.40.2)
          
          men.40.50.women.15.25 <- c(men.40.50.women.15.25.1, women.15.25.men.40.50.2)
          
          men.40.50.women.25.40 <- c(men.40.50.women.25.40.1, women.25.40.men.40.50.2)
          
          men.40.50.women.40.50 <- c(men.40.50.women.40.50.1, women.40.50.men.40.50.2)
          
          Age.groups.table <- matrix(c(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50),
                                       length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50),
                                       length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50)),
                                     ncol = 3,
                                     byrow = TRUE)
          
          colnames(Age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
          rownames(Age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
          
          Age.groups.table <- as.table(Age.groups.table)
          
          
          men.15.25.T <- sum(length(men.15.25.women.15.25), length(men.15.25.women.25.40), length(men.15.25.women.40.50))
          men.25.40.T <- sum(length(men.25.40.women.15.25), length(men.25.40.women.25.40), length(men.25.40.women.40.50))
          men.40.50.T <- sum(length(men.40.50.women.15.25), length(men.40.50.women.25.40), length(men.40.50.women.40.50))
          
          prop.men.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/men.15.25.T, length(men.15.25.women.25.40)/men.15.25.T, length(men.15.25.women.40.50)/men.15.25.T,
                                                length(men.25.40.women.15.25)/men.25.40.T, length(men.25.40.women.25.40)/men.25.40.T, length(men.25.40.women.40.50)/men.25.40.T,
                                                length(men.40.50.women.15.25)/men.40.50.T, length(men.40.50.women.25.40)/men.40.50.T, length(men.40.50.women.40.50)/men.40.50.T),
                                              ncol = 3,
                                              byrow = TRUE)
          
          colnames(prop.men.age.groups.table) <- c("Female.15.25", "Female.25.40", "Female.40.50")
          rownames(prop.men.age.groups.table) <- c("prop.Male.15.25", "prop.Male.25.40", "prop.Male.40.50")
          
          
          
          
          women.15.25.T <- sum(length(men.15.25.women.15.25), length(men.25.40.women.15.25), length(men.40.50.women.15.25))
          women.25.40.T <- sum(length(men.15.25.women.25.40), length(men.25.40.women.25.40), length(men.40.50.women.25.40))
          women.40.50.T <- sum(length(men.15.25.women.40.50), length(men.25.40.women.40.50), length(men.40.50.women.40.50))
          
          prop.women.age.groups.table <- matrix(c(length(men.15.25.women.15.25)/women.15.25.T, length(men.25.40.women.15.25)/women.15.25.T, length(men.40.50.women.15.25)/women.15.25.T,
                                                  length(men.15.25.women.25.40)/women.25.40.T, length(men.25.40.women.25.40)/women.25.40.T, length(men.40.50.women.25.40)/women.25.40.T,
                                                  length(men.15.25.women.40.50)/women.40.50.T, length(men.25.40.women.40.50)/women.40.50.T, length(men.40.50.women.40.50)/women.40.50.T),
                                                ncol = 3,
                                                byrow = TRUE)
          
          colnames(prop.women.age.groups.table) <- c("Male.15.25", "Male.25.40", "Male.40.50")
          rownames(prop.women.age.groups.table) <- c("prop.Female.15.25", "prop.Female.25.40", "prop.Female.40.50")
          
          outputlist <- NULL
          outputlist$Age.groups.table <- Age.groups.table
          outputlist$prop.men.age.groups.table <- prop.men.age.groups.table
          outputlist$prop.women.age.groups.table <- prop.women.age.groups.table
          
          
          return(outputlist)
          
        }
        
        
        
        
        # Results
        ###########
        
        # Function to handle NAs
        
        
        NA.handle.fun <- function(input=input){
          
          v.names <- names(input)
          
          v <- as.numeric(input)
          
          v.vec <- vector()
          
          for(i in 1:length(v)){
            
            v.i <- v[i]
            
            if(is.na(v.i)==TRUE){
              v.j <- 0
            }else{
              v.j <- v.i
            }
            v.vec <- c(v.vec, v.j)
          }
          
          names(v.vec) <- v.names
          return(v.vec)
        }
        
        # 1. Age structure table obtained from transmission network built from transmission clusters
        
        age.structure.transm.clust.List <- age.groups.filtered.trans.clust.network.fun(table.transm.clust.net.igraph = table.transm.clust.net.igraph,
                                                                                       transm.matrix = transm.matrix,
                                                                                       age.group.15.25 = c(15,25),
                                                                                       age.group.25.40 = c(25,40),
                                                                                       age.group.40.50 = c(40,50))
        age.structure.transm.clust <- age.structure.transm.clust.List$Age.groups.table
        
        cl.age.str.M.15.25.F.15.25 <- age.structure.transm.clust[1,][1]
        cl.age.str.M.25.40.F.15.25 <- age.structure.transm.clust[2,][1]
        cl.age.str.M.40.50.F.15.25 <- age.structure.transm.clust[3,][1]
        
        cl.age.str.M.15.25.F.25.40 <- age.structure.transm.clust[1,][2]
        cl.age.str.M.25.40.F.25.40 <- age.structure.transm.clust[2,][2]
        cl.age.str.M.40.50.F.25.40 <- age.structure.transm.clust[3,][2]
        
        cl.age.str.M.15.25.F.40.50 <- age.structure.transm.clust[1,][3]
        cl.age.str.M.25.40.F.40.50 <- age.structure.transm.clust[2,][3]
        cl.age.str.M.40.50.F.40.50 <- age.structure.transm.clust[3,][3]
        
        table.cl.age.str <- c(cl.age.str.M.15.25.F.15.25, cl.age.str.M.25.40.F.15.25, cl.age.str.M.40.50.F.15.25,
                              cl.age.str.M.15.25.F.25.40, cl.age.str.M.25.40.F.25.40, cl.age.str.M.40.50.F.25.40,
                              cl.age.str.M.15.25.F.40.50, cl.age.str.M.25.40.F.40.50, cl.age.str.M.40.50.F.40.50)
        
        names(table.cl.age.str) <- c("cl.M.15.25.F.15.25", "cl.M.25.40.F.15.25", "cl.M.40.50.F.15.25",
                                     "cl.M.15.25.F.25.40", "cl.M.25.40.F.25.40", "cl.M.40.50.F.25.40",
                                     "cl.M.15.25.F.40.50", "cl.M.25.40.F.40.50", "cl.M.40.50.F.40.50")
        
        
        # Men prop
        
        age.structure.transm.clust.prop.men <- age.structure.transm.clust.List$prop.men.age.groups.table
        
        cl.age.str.prop.men.15.25.F.15.25 <- age.structure.transm.clust.prop.men[1,][1]
        cl.age.str.prop.men.25.40.F.15.25 <- age.structure.transm.clust.prop.men[2,][1]
        cl.age.str.prop.men.40.50.F.15.25 <- age.structure.transm.clust.prop.men[3,][1]
        
        cl.age.str.prop.men.15.25.F.25.40 <- age.structure.transm.clust.prop.men[1,][2]
        cl.age.str.prop.men.25.40.F.25.40 <- age.structure.transm.clust.prop.men[2,][2]
        cl.age.str.prop.men.40.50.F.25.40 <- age.structure.transm.clust.prop.men[3,][2]
        
        cl.age.str.prop.men.15.25.F.40.50 <- age.structure.transm.clust.prop.men[1,][3]
        cl.age.str.prop.men.25.40.F.40.50 <- age.structure.transm.clust.prop.men[2,][3]
        cl.age.str.prop.men.40.50.F.40.50 <- age.structure.transm.clust.prop.men[3,][3]
        
        table.cl.age.str.prop.men <- c(cl.age.str.prop.men.15.25.F.15.25, cl.age.str.prop.men.25.40.F.15.25, cl.age.str.prop.men.40.50.F.15.25,
                                       cl.age.str.prop.men.15.25.F.25.40, cl.age.str.prop.men.25.40.F.25.40, cl.age.str.prop.men.40.50.F.25.40,
                                       cl.age.str.prop.men.15.25.F.40.50, cl.age.str.prop.men.25.40.F.40.50, cl.age.str.prop.men.40.50.F.40.50)
        
        names(table.cl.age.str.prop.men) <- c("cl.prop.men15.25.F.15.25", "cl.prop.men25.40.F.15.25", "cl.prop.men40.50.F.15.25",
                                              "cl.prop.men15.25.F.25.40", "cl.prop.men25.40.F.25.40", "cl.prop.men40.50.F.25.40",
                                              "cl.prop.men15.25.F.40.50", "cl.prop.men25.40.F.40.50", "cl.prop.men40.50.F.40.50")
        
        table.cl.age.str.prop.men <- NA.handle.fun(input = table.cl.age.str.prop.men)
        
        
        # Women prop
        
        age.structure.transm.clust.prop.women <- age.structure.transm.clust.List$prop.women.age.groups.table
        
        cl.age.str.prop.women.15.25.M.15.25 <- age.structure.transm.clust.prop.women[1,][1]
        cl.age.str.prop.women.25.40.M.15.25 <- age.structure.transm.clust.prop.women[2,][1]
        cl.age.str.prop.women.40.50.M.15.25 <- age.structure.transm.clust.prop.women[3,][1]
        
        cl.age.str.prop.women.15.25.M.25.40 <- age.structure.transm.clust.prop.women[1,][2]
        cl.age.str.prop.women.25.40.M.25.40 <- age.structure.transm.clust.prop.women[2,][2]
        cl.age.str.prop.women.40.50.M.25.40 <- age.structure.transm.clust.prop.women[3,][2]
        
        cl.age.str.prop.women.15.25.M.40.50 <- age.structure.transm.clust.prop.women[1,][3]
        cl.age.str.prop.women.25.40.M.40.50 <- age.structure.transm.clust.prop.women[2,][3]
        cl.age.str.prop.women.40.50.M.40.50 <- age.structure.transm.clust.prop.women[3,][3]
        
        table.cl.age.str.prop.women <- c(cl.age.str.prop.women.15.25.M.15.25, cl.age.str.prop.women.25.40.M.15.25, cl.age.str.prop.women.40.50.M.15.25,
                                         cl.age.str.prop.women.15.25.M.25.40, cl.age.str.prop.women.25.40.M.25.40, cl.age.str.prop.women.40.50.M.25.40,
                                         cl.age.str.prop.women.15.25.M.40.50, cl.age.str.prop.women.25.40.M.40.50, cl.age.str.prop.women.40.50.M.40.50)
        
        names(table.cl.age.str.prop.women) <- c("cl.prop.women15.25.M.15.25", "cl.prop.women25.40.M.15.25", "cl.prop.women40.50.M.15.25",
                                                "cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40", "cl.prop.women40.50.M.25.40",
                                                "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50")
        
        table.cl.age.str.prop.women <- NA.handle.fun(input = table.cl.age.str.prop.women)
        
        
        
        output.results.vector <- c(table.cl.age.str.prop.men, table.cl.age.str.prop.women)
        
      }else{
        
        output.results.vector <- rep(NA, 18)
        
        names(output.results.vector) <-  c("cl.prop.men15.25.F.15.25", "cl.prop.men25.40.F.15.25", "cl.prop.men40.50.F.15.25",
                                           "cl.prop.men15.25.F.25.40", "cl.prop.men25.40.F.25.40", "cl.prop.men40.50.F.25.40",
                                           "cl.prop.men15.25.F.40.50", "cl.prop.men25.40.F.40.50", "cl.prop.men40.50.F.40.50",
                                           
                                           "cl.prop.women15.25.M.15.25", "cl.prop.women25.40.M.15.25", "cl.prop.women40.50.M.15.25",
                                           "cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40", "cl.prop.women40.50.M.25.40",
                                           "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50")
        
      }
      
      
      return(output.results.vector)
      
      
    }
    
    stats.AG <- stats.age.groups.trans.clust.network.fun(simpact.trans.net = simpact.trans.net, 
                                                         datalist.agemix = datalist.agemix,
                                                         work.dir = work.dir,  
                                                         dirfasttree = dirfasttree, 
                                                         sub.dir.rename = sub.dir.rename,
                                                         limitTransmEvents = limitTransmEvents,
                                                         timewindow = timewindow,
                                                         seq.cov = seq.cov,
                                                         select.IDs = mCAr.IDs,
                                                         age.group.15.25 = age.group.15.25,
                                                         age.group.25.40 = age.group.25.40,
                                                         age.group.40.50 = age.group.40.50,
                                                         cut.off = cut.off)
    
    
    # library(phytools)
    
    # tree.cal <- read.tree(paste0(sub.dir.rename, "/calibrated.tree.nwk"))
    
    tree.cal <- tree.cal.cov.35.IDs # read.tree(paste0(tree.topo))
    
    
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
    
    # # Estimating confidence intervals for rates and dates using a parametric bootstrap
    # pb <- parboot.treedater(tree.calib.LTT) # Lineage Through Time
    # 
    # # Lineages through time
    # LTT <- plot.parboot.ltt.dat(pb)
    # 
    # lb.mean.feature <- mean(LTT$lb) # mean of low values of LTT
    # lb.median.feature <- median(LTT$lb) # median of low values of LTT
    # 
    # ub.mean.feature <- mean(LTT$ub) # mean of upper values of LTT
    # ub.median.feature <- median(LTT$ub) # median of upper values of LTT
    # 
    # median.mean.feature <- mean(LTT$median) # mean of medians of values of LTT
    # median.median.feature <- median(LTT$median) # median of medians of values of LTT
    
    phylo.features.summary <- c(mean.feature, colless.feature, sackin.feature, mean.tipsDepths.feature, mean.nodesDepths.feature,
                                maxHeight.feature)
    # lb.mean.feature, lb.median.feature, ub.mean.feature, ub.median.feature,
    # median.mean.feature, median.median.feature)
    
    features.names <- c("meanHeight", "colless", "sackin", "mean.tipsDepths", "mean.nodesDepths",
                        "maxHeight")
    # , "LTT.lb.mean.feature", "LTT.lb.median.feature", "LTT.ub.mean.feature", "LTT.ub.median.feature",
    #                     "LTT.median.mean.feature", "LTT.median.median.feature")
    
    names(phylo.features.summary) <- features.names
    
    clust.phylo.fit.params <- c(intercept.fit,
                                slope.fit,
                                
                                clust.fit.params, 
                                phylo.features.summary,
                                stats.AG)
    
    
  }else{
    
    clust.phylo.fit.params <- rep(NA, 37)
    
    names(clust.phylo.fit.params) <- c("node.interc.fit",
                                       "node.slope.fit",
                                       
                                       "ME.cl.av.age.male",  "ME.cl.gendEffect", "ME.cl.between.transm.var", "ME.cl.within.transm.var",  "ME.cl.SD.female",          
                                       "ME.cl.SD.male", "Num.Clusters",  "av.Clust.Size",   "mean.Num.Clus",  "med.Num.Clus",   "SD.Num.Clus" , "meanHeight",
                                       "colless",  "sackin", "mean.tipsDepths", "mean.nodesDepths", "maxHeight",  "cl.prop.men15.25.F.15.25",   "cl.prop.men25.40.F.15.25",
                                       "cl.prop.men40.50.F.15.25", "cl.prop.men15.25.F.25.40",   "cl.prop.men25.40.F.25.40",   "cl.prop.men40.50.F.25.40",  
                                       "cl.prop.men15.25.F.40.50" ,  "cl.prop.men25.40.F.40.50" , "cl.prop.men40.50.F.40.50",   "cl.prop.women15.25.M.15.25", 
                                       "cl.prop.women25.40.M.15.25", "cl.prop.women40.50.M.15.25" ,"cl.prop.women15.25.M.25.40", "cl.prop.women25.40.M.25.40", 
                                       "cl.prop.women40.50.M.25.40", "cl.prop.women15.25.M.40.50", "cl.prop.women25.40.M.40.50", "cl.prop.women40.50.M.40.50")
    
    
    
    
  }
  
  return(clust.phylo.fit.params)
  
}
