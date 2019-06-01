#' Calibrating model using latin hypercube for parameter sampling and Approximate Bayesian Computation for selecting best posteriors
#'
#' @param model.sim Model to be calibrated Transmission networks computed by \code{\link{transmission.network.builder()}}
#' @param sum_stat_obs Observed summary statistics
#' @param simpact_prior Priors of the parameters
#' @param design.points Designed points of the sampling in parameter space with latin hypercube sampling
#' @param seed.val Seed value Gender imbalance ratio (women/(women+men))
#' @param n_cores Number of cores for parallelization
#' @importFrom lhs randomLHS
#' @importFrom abc abc
#' @importFrom RSimpactHelper simpact.parallel
#' @export
#'

calibration.ABC <- function(model.sim = simpact4ABC.classic,
                            sum_stat_obs = sum_stat_obs,
                            simpact_prior = simpact_prior, 
                            design.points = 100,
                            seed.val = 1,
                            n_cores = 8){
  
  
  # toy_model_parallel<-function(x){
  #   set.seed(x[1]) # so that each core is initialized with a different seed value.
  #   c( x[2] + x[3] + rnorm(1,0,0.1) , x[2] * x[3] + rnorm(1,0,0.1) )
  # }
  # 
  # toy_prior=list(c("unif",0,1),c("normal",1,2))
  # 
  # sum_stat_obs=c(1.5,0.5)
  
  # Compute sum_stat_obs according to different scenarios
  
  
  # Calibration with ABC approach
  
  # (i). Having plausible ranges for the parameters, sample the parameter spaces by latin hypercube several times
  # (ii). Run the default model and compute the summary statistics
  # (iii). Use different ABC-based methods to fit the model, this will give parameters estimates and associated summary statistics
  
  
  
  
  
  simpact_prior <- simpact_prior
  
  min.v <- vector()
  max.v <- vector()
  
  for( i in 1:length(simpact_prior)){
    
    min.v <- c(min.v, as.numeric(simpact_prior[[i]][[2]]))
    max.v <- c(max.v, as.numeric(simpact_prior[[i]][[3]]))  
    
  }
  
  
  variables <- length(simpact_prior)
  
  set.seed(seed.val)
  
  rlhs <- lhs::randomLHS(design.points, variables)
  
  
  lhs.df <- rlhs
  
  for (j in 1:ncol(lhs.df)){
    
    min.var <- min.v[j]
    max.var <- max.v[j]
    
    for(k in 1:nrow(lhs.df)){
      
      lhs.df[k,j] <- qunif(lhs.df[k,j], min = min.var, max = max.var)
      
    }
    
  }
  
  
  
  par.sim <- inputmatrix <- lhs.df # results of (i): parameter matrix
  
  # par.sim <- par.sim[, -c(ncol(lhs.df))]
  
  # (ii)
  
  stat.sim <- RSimpactHelper::simpact.parallel(model = model.sim, # simpact4ABC.classic,
                                               actual.input.matrix = par.sim,
                                               seed_count = seed.val,
                                               n_cluster = n_cores)
  
  
  # stat.sim <- read.csv("stat.sim.csv")
  
  
  stat.sim <- stat.sim[,1:length(sum_stat_obs)] # results of (ii): summary statistics matrix obtained from simulations done with parameter matrix
  
  
  stat.obs <- sum_stat_obs
  
  
  # (iii)
  
  # Condition: nrow(par.sim) == nrow(stat.sim)
  
  
  # (iv) removing NA
  
  # NA.stat.sim <- na.omit(stat.sim)
  
  index.na.fun <- function(stat.sim = stat.sim){
    index.vec <- vector()
    for(i in 1:nrow(stat.sim)){
      row.i <- stat.sim[i,] 
      if(NA%in%row.i){
        index.vec <- c(index.vec, i)
      }
    }
    return(index.vec)
  }
  
  na.index <- index.na.fun(stat.sim = stat.sim)
  
  na.stat.sim <- stat.sim[-c(na.index),]
  na.par.sim <- par.sim[-c(na.index), ]
  
  
  # (v) ABC calibration
  
  # rej <- abc(target=stat.obs, param=na.par.sim, sumstat=na.stat.sim, tol=.1, method = "rejection") 
  
  
  abc_res <- abc::abc(target=stat.obs, param=par.sim, sumstat=stat.sim, tol=0.1, 
                        method = "ridge") 
  
  # parms.post.adj <- abc_res$adj.values
  # 
  # param.vect <- as.data.frame(abc_res$adj.values)
  
  # Summary
  
  # sum.abc_res <- summary(abc_res, intvl = .9)
  # 
  # par.min <- sum.abc_res[1,]
  # par.weigh.5.perc <- sum.abc_res[2,]
  # par.med <- sum.abc_res[3,]
  # par.mean <- sum.abc_res[4,]
  # par.mod <- sum.abc_res[5,]
  # par.weigh.95 <- sum.abc_res[6,]
  # par.max <- sum.abc_res[7,]
  
  # save(abc_res, file = "abc_res.RData")
  
  return(abc_res)
  
}

