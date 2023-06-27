##### 1) Functions to calibrate the baseline mortality model
FitLiLeeNR = function(xv, yv, Data, Ax, exclude.coh, CountrySPEC, yvSPEC, fit1, fit2){
  
  ### Step 1: Calibrate common, multi-population mortality trend
  
  # Aggregated deaths & exposures
  dtxALL <- Data$ALL$dtx
  etxALL <- Data$ALL$etx
  waALL  <- Data$ALL$wa
  
  # Fit
  LC.total <- fit1(xv, yv, etxALL, dtxALL, waALL, xv*0, "ALL", exclude.coh)
  
  # Extract parameters
  A.x  <- LC.total$beta1
  B.x  <- LC.total$beta2		# note: sum B.x = 1
  C.x  <- LC.total$beta3
  K.t  <- LC.total$kappa2 		# note: sum K.t = 0 
  L.t  <- if(is.null(LC.total$kappa3)) K.t*0 else LC.total$kappa3
  
  residEU <- LC.total$epsilon
  m.txEU  <- LC.total$mhat
  
  ### Step 2: Calibrate country-specification deviation trend 
  
  # Preliminary set-up
  c <- CountrySPEC
  Data$UNI[[c]]$dtx[which(Data$UNI[[c]]$dtx == 0)] <- 1e-8
  
  # Fit
  LC <- fit2(xv = xv, yv = yvSPEC, 
             etx = (Data$UNI[[c]]$etx[as.character(yvSPEC),]) * m.txEU, 
             dtx = Data$UNI[[c]]$dtx[as.character(yvSPEC),], 
             wa = Data$UNI[[c]]$wa[as.character(yvSPEC),], int = xv*0, constraints = "ALL", exclude.coh = exclude.coh) # adjust exposure!
  
  # Extract parameters and more
  alpha.x	  <- LC$beta1
  beta.x 	  <- LC$beta2			# note: sum = 1
  gamma.x   <- LC$beta3
  kappa.t   <- LC$kappa2	 		# note: sum = 0
  lambda.t  <- if(is.null(LC$kappa3)) kappa.t*0 else LC$kappa3
  ll        <- LC$ll
  npar      <- LC$npar + LC.total$npar
  BIC       <- -2*LC$ll + log(sum(LC$wa))*npar
  residuals <- LC$epsilon
  
  # death rates CountrySPEC
  m.tx <- m.txEU*LC$mhat
  
  ### Output
  out <- list("NR" = 
                list(A.x = A.x, B.x = B.x, C.x = C.x, K.t = K.t, L.t = L.t, 
                     a.x = alpha.x, b.x = beta.x, c.x = gamma.x, k.t = kappa.t, 
                     l.t = lambda.t, m.tx = m.tx, ll = ll, npar = npar, 
                     BIC = BIC, residuals = residuals, residuals_EU = residEU, 
                     m.txEU = m.txEU))
  
  out
}

fit701C = function(xv,yv,etx,dtx,wa,int, constraints, exclude.coh){
  
  mtx <- dtx/etx	    # matrix of death rates
  qtx <- 1-exp(-mtx)   # matrix of mortality rates
  
  n  <- length(xv)	# number of ages
  m  <- length(yv)	# number of years
  cy <- (yv[1]-xv[n]):(yv[m]-xv[1])  # cohort approximate years of birth
  
  # initialise parameter vectors
  beta1v  <- int
  beta2v  <- (1:n)*0
  beta3v  <- (1:n)*0		
  kappa2v <- (1:m)*0
  kappa3v <- (1:m)*0
  gamma4v <- (1:(n+m-1))*0	# dummy vector, this will stay at 0
  
  ia  <- array((1:m),c(m,n))	# matrix of year indexes, i, for the data
  ja  <- t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya  <- ia-ja		 	# matrix of year of birth indexes for the data
  imj <- (1-n):(m-1)		# the range of values taken by i-j
  lg  <- n+m-1		 	# number of different values taken by i-j
  ca  <- ya+yv[1]-xv[1]		# matrix of years of birth
  
  # Now set weights to zero for cohorts with fewer than 5 observations
  if(exclude.coh == TRUE){
    for(k in 1:lg)
    {
      nk=sum((ca == cy[k])*wa)
      if(nk < 5)
      {
        wa=wa*(1- (ca == cy[k]))
      }}
  }
  
  ww=cy*0+1	 # this is a vector of 1's and 0's with
  
  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x), beta3(x), kappa2(t) and kappa3(t)
  mx=mean(xv)
  beta1v      <- colMeans(log(mtx[as.character(yv[2]:tail(yv,1)),]) - log(mtx[as.character(yv[1]:tail(yv,2)[1]),]), na.rm = TRUE)
  beta2v <- rep(1/n, length(xv))
  beta3v <- rep(1/n, length(xv))
  
  kappa2v <- 0:(-m+1)/10
  kappa3v <- 0:(-m+1)/10
  
  m0x    <- mtx[1,]
  
  
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  
  while(abs(l1-l0) > 1e-8)
  {
    iteration=iteration+1
    
    l0=l1
    # Stage 1B1 optimise over the beta2(x)
    for(j in 1:n){
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta2v[j]=llmaxM2B1(log(m0x[j]) + beta1v[j]*(0:(m-1)),beta2v[j],beta3v[j],
                          kappa2v,kappa3v,gamma4v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1D1 optimise over the kappa2(t)
    for(i in 2:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]*m0x	# exposure
      kappa2v[i] <- if(all(is.na(dv)) & all(is.na(ev))) NA else
        llmaxM2D1(beta1v*(i-1),beta2v,beta3v, kappa2v[i],kappa3v[i],gamma4v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1B2 optimise over the beta3(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]*m0x[j]	# exposure
      beta3v[j]=llmaxM2B2(beta1v[j]*(0:(m-1)),beta2v[j],beta3v[j],
                          kappa2v,kappa3v,gamma4v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1D2 optimise over the kappa3(t)
    for(i in 2:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]*m0x	# exposure
      kappa3v[i] <- if(all(is.na(dv)) & all(is.na(ev))) NA else 
        llmaxM2D2(beta1v*(i-1),beta2v,beta3v, kappa2v[i],kappa3v[i],gamma4v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1A optimise over the beta1(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]*m0x[j]	# exposure
      beta1v[j]=llmaxM2B0(beta1v[j],beta2v[j],beta3v[j],
                          kappa2v,kappa3v,gamma4v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"\n")
    
    # Now apply the constraints  (small adaptation - multiply delta's with 't' for improvement rate models)
    if(constraints=="ALL"){
      norm  <- function(x) sqrt(sum(x^2))
      inner <- function(x,y) sum(x*y)
      
      model.kappa2v <- forecast::auto.arima(kappa2v, approximation = FALSE, ic = 'bic')$model
      model.kappa3v <- forecast::auto.arima(kappa3v, approximation = FALSE, ic = 'bic')$model
      kappa2v       <- imputeTS::na_ma(kappa2v)
      kappa3v       <- imputeTS::na_ma(kappa3v)
      
      #### see paper for explanation
      B2 <- beta2v; B3 <- beta3v
      K2circ <- kappa2v - (0:(m-1))/(m - 1) * kappa2v[m]; K3circ <- kappa3v - (0:(m-1))/(m - 1) * kappa3v[m]
      
      xi1 <- inner(diff(K2circ),diff(K3circ))*norm(B2)^2 + inner(B2,B3)*norm(diff(K3circ))^2
      xi2 <- inner(diff(K2circ),diff(K3circ))*norm(B3)^2 + inner(B2,B3)*norm(diff(K2circ))^2
      a   <- xi2*inner(diff(K2circ),diff(K3circ))
      b   <- xi1*norm(diff(K2circ))^2 - xi2*norm(diff(K3circ))^2
      c   <- -xi1*inner(diff(K2circ),diff(K3circ))
      D   <- b^2 - 4*a*c
      
      eta2bar  <- c((-b+sqrt(D))/(2*a), (-b-sqrt(D))/(2*a))
      
      zeta2bar = eta1 = zeta1 = eta2 = zeta2 = delta1 = delta2 <- rep(NA,2)
      for (s in eta2bar){
        ind <- which(s == eta2bar)
        zeta2bar[ind] <- - xi1/(xi2*s)
        
        eta1[ind]  <- sign(sum((B2 + s*B3)))/sqrt(sum((B2 + s*B3)^2))
        zeta1[ind] <- sign(sum((B2 + zeta2bar[ind]*B3)))/sqrt(sum((B2 + zeta2bar[ind]*B3)^2))
        eta2[ind]  <- eta1[ind]*s
        zeta2[ind] <- zeta1[ind]*zeta2bar[ind]
        
        delta1[ind] <- -1/(zeta1[ind]*eta2[ind] - eta1[ind]*zeta2[ind]) * (eta2[ind]*kappa2v[m] - eta1[ind]*kappa3v[m])*1/(m-1)
        delta2[ind] <- -1/(zeta1[ind]*eta2[ind] - eta1[ind]*zeta2[ind]) * (-zeta2[ind]*kappa2v[m] + zeta1[ind]*kappa3v[m])*1/(m-1)
      }
      
      r <- which.max(c(diff(range(zeta1[1]*B2 + zeta2[1]*B3)),diff(range(zeta1[2]*B2 + zeta2[2]*B3))))   
      
      Bx <- zeta1[r]*B2 + zeta2[r]*B3
      Cx <- eta1[r]*B2 + eta2[r]*B3
      Kt <- 1/(zeta1[r]*eta2[r] - eta1[r]*zeta2[r]) * (eta2[r]*kappa2v - eta1[r]*kappa3v)  + delta1[r]*(0:(m-1))
      Lt <- 1/(zeta1[r]*eta2[r] - eta1[r]*zeta2[r]) * (-zeta2[r]*kappa2v + zeta1[r]*kappa3v) + delta2[r]*(0:(m-1))
      Ax <- beta1v - B2*(zeta1[r]*delta1[r] + eta1[r]*delta2[r]) - B3*(zeta2[r]*delta1[r] + eta2[r]*delta2[r])
      
      # Ax + Bx*Kt[125] + Cx*Lt[125]
      # beta1v + beta2v*kappa2v[125] + beta3v*kappa3v[125]
      
      beta1v = Ax; beta2v = Bx; beta3v = Cx; kappa2v = Kt; kappa3v = Lt
    } 
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=m0x*exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")	   
    
    
  }		 # end while loop
  
  
  # calculate number of parameters and deduct the number of constraints
  if(constraints=="ALL"){
    npar=length(beta1v)+length(beta2v)+length(na.omit(kappa2v))+length(beta3v)+length(na.omit(kappa3v))-6
  }
  else{
    npar=length(beta1v)+length(beta2v)+length(na.omit(kappa2v))+length(beta3v)+length(na.omit(kappa3v))-4
  }
  
  # Calculate the BIC
  BIC=-2*l1+log(sum(wa))*npar
  
  list(beta1=beta1v,beta2=beta2v,beta3=beta3v,
       kappa2=kappa2v,kappa3=kappa3v,gamma3=gamma4v,x=xv,y=yv,cy=cy,
       wa=wa,epsilon=epsilon,mhat=mhat,ll=l1,BIC=BIC,npar=npar,mtxLastYear = mtx[m,])		 
}

fit701C0 = function(xv,yv,etx,dtx,wa,int,constraints, exclude.coh){
  
  mtx <- dtx/etx	    # matrix of death rates
  qtx <- 1-exp(-mtx)   # matrix of mortality rates
  
  n  <- length(xv)	# number of ages
  m  <- length(yv)	# number of years
  cy <- (yv[1]-xv[n]):(yv[m]-xv[1])  # cohort approximate years of birth
  
  # initialise parameter vectors
  beta1v  <- int
  beta2v  <- (1:n)*0
  beta3v  <- (1:n)*0		
  kappa2v <- (1:m)*0
  kappa3v <- (1:m)*0
  gamma4v <- (1:(n+m-1))*0	# dummy vector, this will stay at 0
  
  ia  <- array((1:m),c(m,n))	# matrix of year indexes, i, for the data
  ja  <- t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya  <- ia-ja		 	# matrix of year of birth indexes for the data
  imj <- (1-n):(m-1)		# the range of values taken by i-j
  lg  <- n+m-1		 	# number of different values taken by i-j
  ca  <- ya+yv[1]-xv[1]		# matrix of years of birth
  
  # Now set weights to zero for cohorts with fewer than 5 observations
  if(exclude.coh == TRUE){
    for(k in 1:lg)
    {
      nk=sum((ca == cy[k])*wa)
      if(nk < 5)
      {
        wa=wa*(1- (ca == cy[k]))
      }}
  }
  
  ww=cy*0+1	 # this is a vector of 1's and 0's with
  # a 0 if the cohort is completely excluded
  # for(k in 1:lg)
  # {
  #     ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
  # }
  
  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x), beta3(x), kappa2(t) and kappa3(t)
  mx=mean(xv)
  beta2v <- rep(1/n, length(xv))
  beta3v <- rep(1/n, length(xv))
  
  kappa2v <- 0:(-m+1)/10
  kappa3v <- 0:(-m+1)/10
  
  m0x    <- mtx[1,]
  
  
  # Stage 1: iterate
  l0=-1000000
  l1=-999999
  iteration=0
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  
  while(abs(l1-l0) > 1e-8)
  {
    iteration=iteration+1
    
    l0=l1
    # Stage 1B1 optimise over the beta2(x)
    for(j in 1:n){
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta2v[j]=llmaxM2B1(beta1v[j]*(0:(m-1)),beta2v[j],beta3v[j],
                          kappa2v,kappa3v,gamma4v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,] = exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1D1 optimise over the kappa2(t)
    for(i in 2:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]	# exposure
      kappa2v[i] <- if(all(is.na(dv)) & all(is.na(ev))) NA else
        llmaxM2D1(beta1v*(i-1),beta2v,beta3v, kappa2v[i],kappa3v[i],gamma4v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1B2 optimise over the beta3(x)
    for(j in 1:n)
    {		 		 
      # cycle through the range of years
      dv=dtx[,j]	# actual deaths
      ev=etx[,j]	# exposure
      beta3v[j]=llmaxM2B2(beta1v[j]*(0:(m-1)),beta2v[j],beta3v[j],
                          kappa2v,kappa3v,gamma4v[(n+1-j):(n+m-j)],dv,ev,wv=wa[,j])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Stage 1D2 optimise over the kappa3(t)
    for(i in 2:m)
    {		 		 
      # cycle through the range of years
      dv=dtx[i,]	# actual deaths
      ev=etx[i,]	# exposure
      kappa3v[i] <- if(all(is.na(dv)) & all(is.na(ev))) NA else 
        llmaxM2D2(beta1v*(i-1),beta2v,beta3v, kappa2v[i],kappa3v[i],gamma4v[(n+i-1):i],dv,ev,wv=wa[i,])
    }
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")
    
    # Now apply the constraints  (small adaptation - multiply delta's with 't' for improvement rate models)
    if(constraints=="ALL"){
      norm  <- function(x) sqrt(sum(x^2))
      inner <- function(x,y) sum(x*y)
      
      model.kappa2v <- forecast::auto.arima(kappa2v, approximation = FALSE, ic = 'bic')$model
      model.kappa3v <- forecast::auto.arima(kappa3v, approximation = FALSE, ic = 'bic')$model
      kappa2v       <- imputeTS::na_ma(kappa2v)
      kappa3v       <- imputeTS::na_ma(kappa3v)
      
      #### see paper for explanation
      B2 <- beta2v; B3 <- beta3v
      K2circ <- kappa2v; K3circ <- kappa3v
      
      xi1 <- inner(diff(K2circ),diff(K3circ))*norm(B2)^2 + inner(B2,B3)*norm(diff(K3circ))^2
      xi2 <- inner(diff(K2circ),diff(K3circ))*norm(B3)^2 + inner(B2,B3)*norm(diff(K2circ))^2
      a   <- xi2*inner(diff(K2circ),diff(K3circ))
      b   <- xi1*norm(diff(K2circ))^2 - xi2*norm(diff(K3circ))^2
      c   <- -xi1*inner(diff(K2circ),diff(K3circ))
      D   <- b^2 - 4*a*c
      
      eta2bar  <- c((-b+sqrt(D))/(2*a), (-b-sqrt(D))/(2*a))
      
      zeta2bar = eta1 = zeta1 = eta2 = zeta2 = delta1 = delta2 <- rep(NA,2)
      for (s in eta2bar){
        ind <- which(s == eta2bar)
        zeta2bar[ind] <- - xi1/(xi2*s)
        
        eta1[ind]  <- sign(sum((B2 + s*B3)))/sqrt(sum((B2 + s*B3)^2))
        zeta1[ind] <- sign(sum((B2 + zeta2bar[ind]*B3)))/sqrt(sum((B2 + zeta2bar[ind]*B3)^2))
        eta2[ind]  <- eta1[ind]*s
        zeta2[ind] <- zeta1[ind]*zeta2bar[ind]
      }
      
      r <- which.max(c(diff(range(zeta1[1]*B2 + zeta2[1]*B3)),diff(range(zeta1[2]*B2 + zeta2[2]*B3))))   
  
      
      Bx <- zeta1[r]*B2 + zeta2[r]*B3
      Cx <- eta1[r]*B2 + eta2[r]*B3
      Kt <- 1/(zeta1[r]*eta2[r] - eta1[r]*zeta2[r]) * (eta2[r]*kappa2v - eta1[r]*kappa3v) 
      Lt <- 1/(zeta1[r]*eta2[r] - eta1[r]*zeta2[r]) * (-zeta2[r]*kappa2v + zeta1[r]*kappa3v)
      Ax <- beta1v
      
      # Ax + Bx*Kt[125] + Cx*Lt[125]
      # beta1v + beta2v*kappa2v[125] + beta3v*kappa3v[125]
      
      beta1v = Ax; beta2v = Bx; beta3v = Cx; kappa2v = Kt; kappa3v = Lt
    } 
    
    mhat=mtx*0
    for(i in 1:m)
    {
      mhat[i,]=exp(beta1v*(i-1)+beta2v*kappa2v[i]+beta3v*kappa3v[i]+gamma4v[(n+i-1):i])
    }
    epsilon=(dtx-etx*mhat)/sqrt(etx*mhat)
    l1=sum((dtx*log(etx*mhat)-etx*mhat-lgamma(dtx+1))*wa, na.rm = TRUE)
    # cat(l1,"->")	   
    
    
  }		 # end while loop
  
  
  # calculate number of parameters and deduct the number of constraints
  if(constraints=="ALL"){
    npar=length(beta1v)+length(beta2v)+length(na.omit(kappa2v))+length(beta3v)+length(na.omit(kappa3v))-6
  }
  else{
    npar=length(beta1v)+length(beta2v)+length(na.omit(kappa2v))+length(beta3v)+length(na.omit(kappa3v))-4
  }
  
  # Calculate the BIC
  BIC=-2*l1+log(sum(wa))*npar
  
  list(beta1=beta1v,beta2=beta2v,beta3=beta3v,
       kappa2=kappa2v,kappa3=kappa3v,gamma3=gamma4v,x=xv,y=yv,cy=cy,
       wa=wa,epsilon=epsilon,mhat=mhat,ll=l1,BIC=BIC,npar=npar,mtxLastYear = mtx[m,])		 
}

llmaxM2B0=function(b1,b2,b3,k2,k3,g4,dv,ev,wv=1){
  #   b1,b3,k2,k3,g4 are given
  #   solve for b2
  b11=b1
  b10=b11-1
  Tmax = length(ev) - 1
  thetat=(0:Tmax)*ev*exp(b2*k2+b3*k3)
  s1=sum(dv*(0:Tmax)*wv, na.rm = TRUE)
  
  while(abs(b11-b10) > 0.1)
  {
    b10=b11; 
    f0=sum((exp(b10*(0:Tmax))*thetat)*wv, na.rm = TRUE)-s1;
    df0=sum((exp(b10*(0:Tmax))*(0:Tmax)*thetat)*wv, na.rm = TRUE)
    b11=b10-f0/df0
  }
  b11
} 

llmaxM2B1=function(b1,b2,b3,k2,k3,g4,dv,ev,wv=1){
  #   b1,b3,k2,k3,g4 are given
  #   solve for b2
  b21=b2
  b20=b21-1
  thetat=k2*ev*exp(b1+b3*k3)
  s1=sum(dv*k2*wv, na.rm = TRUE)
  
  while(abs(b21-b20) > 0.1)
  {
    b20=b21; 
    f0=sum((exp(b20*k2)*thetat)*wv, na.rm = TRUE)-s1;
    df0=sum((exp(b20*k2)*k2*thetat)*wv, na.rm = TRUE)
    b21=b20-f0/df0
  }
  b21
}

llmaxM2B2=function(b1,b2,b3,k2,k3,g4,dv,ev,wv=1){
  #   b1,b2,k2,k3,g4 are given
  #   solve for b3
  b31=b3
  b30=b31-1
  thetat=k3*ev*exp(b1+b2*k2)
  s1=sum(dv*k3*wv, na.rm = TRUE)
  
  while(abs(b31-b30) > 0.1)
  {
    b30=b31; 
    f0=sum((exp(b30*k3)*thetat)*wv, na.rm = TRUE)-s1;
    df0=sum((exp(b30*k3)*k3*thetat)*wv, na.rm = TRUE)
    b31=b30-f0/df0
  }
  b31
}

llmaxM2D1=function(b1,b2,b3,k2,k3,g4,dv,ev,wv=1){
  #   b1,b2,b3,k3,g4 are given
  #   solve for k2
  k21=k2
  k20=k21-1
  thetat=b2*ev*exp(b1+b3*k3)
  s1=sum(dv*b2*wv, na.rm = TRUE)
  while(abs(k21-k20) > 0.1)
  {
    k20=k21
    f0=sum((exp(k20*b2)*thetat)*wv, na.rm = TRUE)-s1
    df0=sum((exp(k20*b2)*b2*thetat)*wv, na.rm = TRUE)
    k21=k20-f0/df0
  }
  k21
} 

llmaxM2D2=function(b1,b2,b3,k2,k3,g4,dv,ev,wv=1){
  #   b1,b2,b3,k2,g4 are given
  #   solve for k3
  k31=k3
  k30=k31-1
  thetat=b3*ev*exp(b1+b2*k2)
  s1=sum(dv*b3*wv, na.rm = TRUE)
  while(abs(k31-k30) > 0.1)
  {
    k30=k31
    f0=sum((exp(k30*b3)*thetat)*wv, na.rm = TRUE)-s1
    df0=sum((exp(k30*b3)*b3*thetat)*wv, na.rm = TRUE)
    k31=k30-f0/df0
  }
  k31
}
