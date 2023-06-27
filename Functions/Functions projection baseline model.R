# Decay function
geom_decay <- function(t,TL,rho) rho^(TL - t)

# Log Likelihood function multivariate time series
obj.logL <- function(x, lambda, df) {
  dim   <- ncol(df)
  Sigma <- matrix(0, ncol = dim, nrow = dim)
  Sigma[lower.tri(Sigma, diag = TRUE)] <- x[-c(1:2)]
  Sigma   <- Sigma %*% t(Sigma)
  weights <- geom_decay(t   = as.numeric(rownames(df)), 
                        TL  = tail(as.numeric(rownames(df)),1), 
                        rho = lambda)
  - sum(weights*log(mvtnorm::dmvnorm(df, mean = c(x[1:2],0,0), sigma = Sigma)))
}

# Function to find optimal decay parameter
opt.decay <- function(df.full, n.years = 121){
  # Grid
  lambda_vec <- seq(0.92,0.98,0.001)
  
  # Parallel processing - faster
  numCores <- detectCores()
  cl <- makeCluster(numCores-2)
  registerDoParallel(cl)
  obj <- foreach(lambda = lambda_vec, 
                 .export = c('obj.logL', 'geom_decay', 'df.full','lambda_vec')) %dopar% {
                   
                   ll <- rep(NA, 100)
                   for(j in 1:n.years){
                     df <- head(df.full, -j)
                     
                     nll <- nlminb(start = c(0,0,1,0,0,0,1,0,0,1,0,1), objective = obj.logL,
                                   df = df, lambda = lambda, gradient = NULL, hessian = NULL,
                                   control = list(trace = F, eval.max = 10000, 
                                                  iter.max = 10000, rel.tol = 1e-8))
                     
                     mat <- matrix(0, nrow = ncol(df.full), ncol = ncol(df.full))
                     mat[lower.tri(mat,diag=T)] <- nll$par[-c(1:2)]
                     Sigma <- mat %*% t(mat)
                     
                     yv.pred <- as.character(tail(as.numeric(rownames(df)),1)+1)
                     ll[j] <- log(mvtnorm::dmvnorm(df.full[yv.pred,], 
                                                   mean = c(nll$par[1:2],0,0), sigma = Sigma))
                   }
                   ll
                 }
  stopCluster(cl)
  
  results <- data.frame('Lambda' = lambda_vec, 'LogLik' = sapply(obj, sum))
  ind.max <- which.max(results$LogLik)
  results$Lambda[ind.max]
}

# Construct data frame for plotting the simulated period effects
construct_plotdf <- function(yv, listProj, ArimaSpec, type, save = FALSE){
  type.t   <- type
  Sex      <- if(substr(type.t,5,5) == "M") "Male" else "Female"
  ts_model <- ArimaSpec[[type.t]]
  sim      <- listProj[[type.t]]
  if(is.list(sim)) sim <- sim[["BEL"]]
  
  # Dimensions
  nSim 	 <- dim(sim)[2]			# number of simulations
  nAhead <- dim(sim)[1]-length(yv)	# number of years to project in future
  
  # Compute percentiles 
  vP <- c(0.005, 0.5, 0.995)
  mQuantiles <- array(NA, dim = c(3, nAhead+1))
  
  for(k in 1:nAhead)
    mQuantiles[,(k+1)] <- quantile(sim[length(yv)+k,], probs = vP) 
  
  mQuantiles[,1] <- sim[length(yv),1]
  dimnames(mQuantiles) <- list(c('Q1', 'Q2', 'Q3'), tail(yv,1):(tail(yv,1)+nAhead))
  
  # Right format to plot
  df_quant <- as.data.frame(t(mQuantiles)) %>% tibble::rownames_to_column('Year') %>% tidyr::gather(key = "Quantile", value = "Value", - Year) %>%
    mutate('Year' = as.numeric(Year))
  
  # Observed and fitted data
  df_fit <- data.frame('Year' = yv, "Quantile" = 'InSampleFit', "Value" = sim[1:length(yv),1])
  
  # Merge
  df <- rbind(df_quant, df_fit)
  df
}