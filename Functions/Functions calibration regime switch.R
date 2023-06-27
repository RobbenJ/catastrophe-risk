# Objective function with as parameters p12, p21, sigma, slope, muH, sigmaH
obj.param1 <- function(x, C, ztx, W = rep(1, nrow(ztx))){
  # Preliminary 
  m     <- ncol(ztx)
  n     <- nrow(ztx)
  years <- as.numeric(rownames(ztx))
  
  # Transition probabilities
  p12 <- x[1]
  p11 <- 1 - p12
  p21 <- x[2]
  p22 <- 1 - p21
  
  # Mean and covariance matrix Z_{t}
  # Different covariance matrix before and after 1970
  mu1  <- rep(0,m)
  mu2  <- C*x[3]
  var1.1 <- diag((x[4] + x[5]*(0:(m-1)))^2)
  var1.2 <- diag((x[6] + x[7]*(0:(m-1)))^2)
  var2.1 <- C %*% t(C) * x[8]^2 + diag((x[4] + x[5]*(0:(m-1)))^2) 
  var2.2 <- C %*% t(C) * x[8]^2 + diag((x[6] + x[7]*(0:(m-1)))^2)
  
  # Invariant probabilities
  pi1 <- p21/(p21 + p12 + p12*p21)
  pi2 <- p12 * pi1
  pi3 <- p12/p21 * pi1
  
  # Create empty objects
  dens  <- rep(NA, n)
  cprob <- matrix(NA, ncol = n, nrow = 3)
  dimnames(cprob) <- list(1:3, years)
  names(dens) <- years
  
  # Multivariate normal density conditionally on being in regime sate 1 or 2
  x1 <- mvtnorm::dmvnorm(ztx[as.character(years[1]:1969),], mu1, var1.1)
  x2 <- mvtnorm::dmvnorm(ztx[as.character(1970:2021),], mu1, var1.2)
  y1 <- mvtnorm::dmvnorm(ztx[as.character(years[1]:1969),], mu2, var2.1)
  y2 <- mvtnorm::dmvnorm(ztx[as.character(1970:2021),], mu2, var2.2)
  
  dmv1 <- c(x1, x2)
  dmv2 <- c(y1, y2)
  dmv1[which(dmv1 == 0)] <- .Machine$double.xmin
  dmv2[which(dmv2 == 0)] <- .Machine$double.xmin
  
  # Starting values recursion
  dens[1]    <- pi1 * dmv1[1] + (pi2 + pi3) * dmv2[1]
  cprob[1,1] <- pi1 * dmv1[1]/dens[1]
  cprob[2,1] <- pi2 * dmv2[1]/dens[1]
  cprob[3,1] <- pi3 * dmv2[1]/dens[1]
  
  dens[2]  <- cprob[1,1] * (p11 * dmv1[2] + p12 * dmv2[2]) + 
    cprob[2,1] * (1 * dmv2[2]) + 
    cprob[3,1] * (p21 * dmv1[2] + p22 * dmv2[2])
  
  # Recursion
  for(t in 2:(n-1)){
    cprob[1,t] <- (cprob[1,t-1] * p11 + cprob[2,t-1] * 0 + cprob[3,t-1] * p21) * dmv1[t]/dens[t]
    cprob[2,t] <- (cprob[1,t-1] * p12 + cprob[2,t-1] * 0 + cprob[3,t-1] * 0)   * dmv2[t]/dens[t]
    cprob[3,t] <- (cprob[1,t-1] * 0   + cprob[2,t-1] * 1 + cprob[3,t-1] * p22) * dmv2[t]/dens[t]
    
    dens[t+1] <- cprob[1,t] * (p11 * dmv1[t+1] + p12 * dmv2[t+1]) + 
      cprob[2,t] * (1 * dmv2[t+1]) +
      cprob[3,t] * (p21 * dmv1[t+1] + p22 * dmv2[t+1]) 
  }
  
  # Conditional probabilities in last year
  cprob[1,n] <- (cprob[1,n-1] * p11 + cprob[2,n-1] * 0 + cprob[3,n-1] * p21) * dmv1[n]/dens[n]
  cprob[2,n] <- (cprob[1,n-1] * p12 + cprob[2,n-1] * 0 + cprob[3,n-1] * 0)   * dmv2[n]/dens[n]
  cprob[3,n] <- (cprob[1,n-1] * 0   + cprob[2,n-1] * 1 + cprob[3,n-1] * p22) * dmv2[n]/dens[n]
  
  # Minus (weighted) log-likelihood
  logL <- sum(W*log(dens))
  
  - logL 
}

# Objective function with as parameters the shock's age specific effect
obj.param2 <- function(x, param1, ztx, W = rep(1, nrow(ztx))){
  # Preliminary 
  m <- ncol(ztx)
  n <- nrow(ztx)
  years <- as.numeric(rownames(ztx))
  
  # Allocate
  p12 = param1[1]; p21 = param1[2]; muH = param1[3]; sigmaE1 = param1[4];
  slope1 = param1[5]; sigmaE2 = param1[6]; slope2 = param1[7]; sigmaH = param1[8]
  
  # Transition probabilities
  p11  <- 1 - p12
  p22  <- 1 - p21
  
  # Mean and covariance matrix Z_{t}
  # Different covariance matrix before and after 1970
  mu1  <- rep(0,m)
  mu2  <- x*muH
  var1.1 <- diag((sigmaE1 + slope1*(0:(m-1)))^2)
  var1.2 <- diag((sigmaE2 + slope2*(0:(m-1)))^2)
  var2.1 <- x %*% t(x) * sigmaH^2 + diag((sigmaE1 + slope1*(0:(m-1)))^2)
  var2.2 <- x %*% t(x) * sigmaH^2 + diag((sigmaE2 + slope2*(0:(m-1)))^2)
  
  # Invariant probabilities
  pi1 <- p21/(p21 + p12 + p12*p21)
  pi2 <- p12 * pi1
  pi3 <- p12/p21 * pi1
  
  # Create empty objects
  dens  <- rep(NA, n)
  cprob <- matrix(NA, ncol = n, nrow = 3)
  dimnames(cprob) <- list(1:3, years)
  names(dens) <- years
  
  # Multivariate normal density conditionally on being in regime sate 1 or 2
  x1 <- mvtnorm::dmvnorm(ztx[as.character(years[1]:1969),], mu1, var1.1)
  x2 <- mvtnorm::dmvnorm(ztx[as.character(1970:2021),], mu1, var1.2)
  y1 <- mvtnorm::dmvnorm(ztx[as.character(years[1]:1969),], mu2, var2.1)
  y2 <- mvtnorm::dmvnorm(ztx[as.character(1970:2021),], mu2, var2.2)
  
  dmv1 <- c(x1, x2)
  dmv2 <- c(y1, y2)
  dmv1[which(dmv1 == 0)] <- .Machine$double.xmin
  dmv2[which(dmv2 == 0)] <- .Machine$double.xmin
  
  # Starting values recursion
  dens[1]    <- pi1 * dmv1[1] + (pi2 + pi3) * dmv2[1]
  cprob[1,1] <- pi1 * dmv1[1]/dens[1]
  cprob[2,1] <- pi2 * dmv2[1]/dens[1]
  cprob[3,1] <- pi3 * dmv2[1]/dens[1]
  
  dens[2]  <- cprob[1,1] * (p11 * dmv1[2] + p12 * dmv2[2]) + 
    cprob[2,1] * (1 * dmv2[2]) + 
    cprob[3,1] * (p21 * dmv1[2] + p22 * dmv2[2])
  
  # Recursion
  for(t in 2:(n-1)){
    cprob[1,t] <- (cprob[1,t-1] * p11 + cprob[2,t-1] * 0 + cprob[3,t-1] * p21) * dmv1[t]/dens[t]
    cprob[2,t] <- (cprob[1,t-1] * p12 + cprob[2,t-1] * 0 + cprob[3,t-1] * 0)   * dmv2[t]/dens[t]
    cprob[3,t] <- (cprob[1,t-1] * 0   + cprob[2,t-1] * 1 + cprob[3,t-1] * p22) * dmv2[t]/dens[t]
    
    dens[t+1] <- cprob[1,t] * (p11 * dmv1[t+1] + p12 * dmv2[t+1]) + 
      cprob[2,t] * (1 * dmv2[t+1]) +
      cprob[3,t] * (p21 * dmv1[t+1] + p22 * dmv2[t+1]) 
  }
  
  # Conditional probabilities in last year
  cprob[1,n] <- (cprob[1,n-1] * p11 + cprob[2,n-1] * 0 + cprob[3,n-1] * p21) * dmv1[n]/dens[n]
  cprob[2,n] <- (cprob[1,n-1] * p12 + cprob[2,n-1] * 0 + cprob[3,n-1] * 0)   * dmv2[n]/dens[n]
  cprob[3,n] <- (cprob[1,n-1] * 0   + cprob[2,n-1] * 1 + cprob[3,n-1] * p22) * dmv2[n]/dens[n]
  
  # Minus (weighted) log-likelihood
  logL <- sum(W*log(dens))
  
  - logL 
}

# Maximalisation function
optim.rs <- function(ztx, init1, init2, ww.exclude, rel.precision = 10^(-5)){

  # Initialisation
  W      <- rep(1, nrow(ztx))
  m      <- ncol(ztx)
  n      <- nrow(ztx)
  param1 <- init1
  param2 <- init2
  ll.new <- 1000
  ll.old <- 990
  iter   <- 0
  
  # Exclude world wars 
  if (ww.exclude) W[which(yv[-1] %in% c(1914:1919, 1940:1946))] <- 0
  
  # Optimization
  while(abs((ll.old - ll.new)/ll.new) > rel.precision){
    
    # Update log-likelihood
    ll.old <- ll.new 
    
    # Step 1 of Optimization
    v1 <- JDEoptim(lower = c(rep(0.001,2),-1, rep(c(0.005,-0.5),2), 0.005), 
                   upper = c(rep(0.999,2), 1, rep(c(1,0.5),2), 5),
                   fn = obj.param1, trace = TRUE, triter = 50, tol = 10^(-6),
                   add_to_init_pop = param1, maxiter = 500, 
                   C = param2, ztx = ztx, W = W)
    
    param1 <- v1$par
    ll1    <- v1$value
    
    # Step 2 of Optimization
    v2 <- JDEoptim(lower = rep(-1, m), upper = rep(1, m), fn = obj.param2,
                   add_to_init_pop = param2, trace = TRUE, triter = 50,
                   maxiter = 500, tol = 10^(-6), ztx = ztx, W = W, 
                   param1 = param1)
    
    param2 <- v2$par
    ll2    <- v2$value
    
    # Apply identifiability constraints 
    fac    <- sqrt(sum(param2^2))
    param2 <- param2/fac * sign(param1[3])
    param1[3] <- param1[3] * fac * sign(param1[3]) 
    param1[8] <- param1[8] * fac
    
    #
    if(ll2 > ll1) break
    
    ll.new <- ll2
    iter   <- iter + 1
  }
  
  names(param1) <- c('p12', 'p21', 'muH', 'sigmaE1', 'slope1', 'sigmaE2',
                     'slope2', 'sigmaH')
  list('param1' = param1, 'param2' = param2, 'LogLik' = -ll2, 'Iter' = iter)
}
