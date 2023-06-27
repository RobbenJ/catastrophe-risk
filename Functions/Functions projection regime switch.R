# Function to generate trajectories for regime switching model
SimulateRegimeSwitch <- function(xv, yv, p12, p21, Cvec, muH, sigmaH, sigmaE2, 
                                 slope2, nSim, nAhead){
  
  # 1) Simulate future states for each scenario
  states <- matrix(NA, nrow = nSim, ncol = nAhead)
  dimnames(states) <- list(1:nSim, (tail(yv,1)+1):(tail(yv,1)+nAhead))
  if(length(Cvec) == 40){
    states[,'2022'] = states[,'2023'] <- 1 # Force regime state 1 in 2022
  } else {
    states[,'2022'] <- 3
    states[,'2023'] <- 1
  }
  
  m <- length(xv)
  
  transprob <- matrix(c(1-p12,p12,0,0,0,1,p21,0,1-p21), nrow = 3, ncol = 3, 
                      byrow = TRUE)
  
  for (t in 3:nAhead){
    ind.tm1.st1 <- which(states[,as.character(tail(yv,1) + t - 1)] == 1)
    ind.tm1.st2 <- which(states[,as.character(tail(yv,1) + t - 1)] == 2)
    ind.tm1.st3 <- which(states[,as.character(tail(yv,1) + t - 1)] == 3)
    
    smpl1 <- sample(1:3, size = length(ind.tm1.st1), replace = TRUE, 
                    prob = transprob[1,])
    smpl2 <- sample(1:3, size = length(ind.tm1.st2), replace = TRUE, 
                    prob = transprob[2,])
    smpl3 <- sample(1:3, size = length(ind.tm1.st3), replace = TRUE, 
                    prob = transprob[3,])
    
    states[ind.tm1.st1,as.character(tail(yv,1) + t)] <- smpl1
    states[ind.tm1.st2,as.character(tail(yv,1) + t)] <- smpl2
    states[ind.tm1.st3,as.character(tail(yv,1) + t)] <- smpl3
  }
  
  # Put state 3 at state 2 (state 2 is just memory state)
  states[which(states == 3, arr.ind = TRUE)] <- 2
  
  # 2) Simulate regime switch component for each scenario
  
  # For each simulated Markov chain - we find the sequences of the same regime 
  # and split the vector
  rle <- apply(states, 1, function(x) data.frame('Lengths' = rle(x)$lengths, 
                                                 'Values'  = rle(x)$values))
  
  # This part creates a sample of size j from a multivariate normal  
  # distribution in the high volatility regime such that the sampled 
  # observations sum up to zero (shock dissappears)
  # and puts the positive shock in front (re-ordering)
  
  avg <- function(j){
    if(j == 1) return(t(mvtnorm::rmvnorm(n = j, mean = Cvec * muH, 
                                         sigma = Cvec %*% t(Cvec) * sigmaH^2)))
    obj  <- t(mvtnorm::rmvnorm(n = j - 1, mean = Cvec * muH, 
                               sigma = Cvec %*% t(Cvec) * sigmaH^2))
    obj  <- cbind(obj, - rowSums(obj))
    
    sums <- colMeans(obj)
    
    idx   <- which.max(sums)
    order <- c(idx, c(1:j)[-idx]) 
    obj[,c(order)]
    
    obj[,c(order)]
  }
  
  # Perform calculations - parallel 
  numCores <- detectCores()
  cl <- makeCluster(numCores-2)
  registerDoParallel(cl)
  
  list <- foreach(sim = 1:nSim, .export = c('states','m','sigmaE2','sigmaH', 
                                            'muH','avg','slope2','Cvec',
                                            'ztx.old'), 
                  .packages = c('mvtnorm')) %dopar% {
                    path <- states[sim,]
                    obj  <- matrix(NA, nrow = m, ncol = length(path))
                    df   <- rle[[sim]]
                    
                    n.lvr <- sum(df$Lengths[which(df$Values == 1)])
                    lvr   <- mvtnorm::rmvnorm(n = n.lvr, mean = rep(0,m), sigma = diag(0,m)) #sigmaE2^2*diag(1 - slope2*(0:(m-1)))
                    n.hvr <- df$Lengths[which(df$Values == 2)]
                    if(length(Cvec) == 26) n.hvr.red <- n.hvr[-1] else n.hvr.red <- n.hvr
                    hvr   <- do.call('cbind', lapply(n.hvr.red, function(k) avg(k)))
                    if(length(Cvec) == 26){
                      hvr <- cbind(-colSums(ztx.old[as.character(c(2020,2021)),as.character(60:85)]),hvr)
                    }
                    
                    obj[,which(path == 1)] <- t(lvr)
                    obj[,which(path == 2)] <- if(length(n.hvr) > 0) hvr
                    
                    obj
                  }
  stopCluster(cl)
  
  regimecomp <- aperm(array(unlist(list), dim = c(m, nAhead, nSim)), c(3,1,2))
  dimnames(regimecomp) <- list(1:nSim, xv, (tail(yv,1)+1):(tail(yv,1)+nAhead))
  
  list(states = states, regime = regimecomp)
}