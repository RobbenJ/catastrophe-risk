# Project mortality rates - simulation
project.mortality.rates  <- function(xv, yv, fitList, projList, nAhead, Sim.vec, 
                                     CountrySPEC, regime){
  TL 	 <- length(yv)		  # number of years in calibration period
  m 	 <- length(xv)		  # number of ages
  fit  <- fitList$NR
  proj <- projList
  yvv  <- yv[1]:(tail(yv,1) + nAhead)
  
  # country specific mortality
  m.xt.obs       <- t(data_M$UNI[[CountrySPEC]]$dtx/data_M$UNI[[CountrySPEC]]$etx)
  m.xti          <- array(NA, dim = c(m,length(yvv),length(Sim.vec)), dimnames = list(xv, yvv, Sim.vec))	# initialize
  m.xti[,1:TL,]  <- m.xt.obs[as.character(xv),]	# fitted mortality
  
  for(t in (TL+1):(TL+nAhead)){ 	# projected death rates: use most recent death rate as starting point
    
    if (t == TL + 1){
      # European mortality trend + country-specific trend
      m.xti[,t,] <- matrix(m.xt.obs[as.character(xv),t-1], ncol = 1)[,rep(1,length(Sim.vec))]*
        exp(matrix(fit$A.x[xv - 20 + 1], ncol = 1)[,rep(1,length(Sim.vec))] + 
              fit$B.x[xv - 20 + 1] %o% proj$K.t_M[t-1,Sim.vec] + 
              fit$C.x[xv - 20 + 1] %o% proj$L.t_M[t-1,Sim.vec] + 
              fit$b.x[xv - 20 + 1] %o% proj$k.t_M[t-1,Sim.vec] + 
              fit$c.x[xv - 20 + 1] %o% proj$l.t_M[t-1,Sim.vec] +
              t(regime[Sim.vec,,t-TL]))
    } else {
      # European mortality trend + country-specific trend
      m.xti[,t,] <- m.xti[,t-1,]* 
        exp(matrix(fit$A.x[xv - 20 + 1], ncol = 1)[,rep(1,length(Sim.vec))] +
              fit$B.x[xv - 20 + 1] %o% proj$K.t_M[t-1,Sim.vec] + 
              fit$C.x[xv - 20 + 1] %o% proj$L.t_M[t-1,Sim.vec] + 
              fit$b.x[xv - 20 + 1] %o% proj$k.t_M[t-1,Sim.vec] + 
              fit$c.x[xv - 20 + 1] %o% proj$l.t_M[t-1,Sim.vec] +
              t(regime[Sim.vec,,t-TL]))
    }
    
  }
  
  
  q.xti <- 1-exp(-m.xti)		# projected mortality rates
  
  q.xti
}

# Construct data frame to plot projected mortality rates
plotdf_qxt <- function(xv, yv, Data, Age, qProj, qProjAG, Sex, CountrySPEC) {
  # Specifications
  nSim 		    <- dim(qProj)[3]			
  nAhead.new  <- dim(qProj)[2] - length(yv)	
  yv.AG     <- 1970:2019
  nAhead.AG <- dim(qProjAG$Male)[2] - length(yv.AG)
  c         <- CountrySPEC
  
  # Age s
  q.mat.new	 <- qProj[as.character(Age),,]
  q.mat.AG   <- qProjAG$Male[as.character(Age),,]
  q.fit.new  <- 1-exp(-mtx.fit[,as.character(Age)]) 
  
  # observed mortality rates q(t,x)
  m.obs	<- (Data$UNI[[c]]$dtx/Data$UNI[[c]]$etx)[,as.character(Age)]
  q.obs	<- 1-exp(-m.obs)	
  
  # Define range y-axis
  bounds <- c(min(c(q.mat.AG, q.mat.new, q.obs)), 
              max(c(q.mat.AG, q.mat.new, q.obs)))
  
  # Percentiles 
  vP <- c(0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995)
  
  # Compute quantiles
  mQuantiles.new  <- array(NA, dim = c(7, nAhead.new + 1))
  mQuantiles.AG <- array(NA, dim = c(7, nAhead.AG + 1))
  dimnames(mQuantiles.AG) <- list(paste0('Q',1:7), 
                                  as.character(tail(yv.AG,1):
                                                 (tail(yv.AG,1) + nAhead.AG)))
  dimnames(mQuantiles.new)  <- list(paste0('Q',1:7), 
                                    as.character(tail(yv,1):
                                                   (tail(yv,1) + nAhead.new)))
  
  # Compute quantiles AG version 2020-2021
  mQuantiles.AG[,"2020"] <- quantile(q.mat.AG["2020",], probs = vP, na.rm = TRUE) 
  mQuantiles.AG[,"2021"] <- quantile(q.mat.AG["2021",], probs = vP, na.rm = TRUE) 
  
  for(y in 2022:2080){
    mQuantiles.new[,as.character(y)]  <- quantile(q.mat.new[as.character(y),], 
                                                  probs = vP, na.rm = TRUE) 
    mQuantiles.AG[,as.character(y)]   <- quantile(q.mat.AG[as.character(y),], 
                                                  probs = vP, na.rm = TRUE) 
  }
  
  # Last observation
  mQuantiles.new[,"2021"]  <- q.mat.new["2021",1] 
  mQuantiles.AG[,"2019"] <- q.mat.AG["2019",1] 
  
  # Right format to plot
  df_quant.new  <- as.data.frame(t(mQuantiles.new)) %>% 
    rownames_to_column('Year') %>% 
    tidyr::gather(key = "Quantile", value = "Value", - Year) %>%
    mutate('Year' = as.numeric(Year), 'Type' = 'New')
  df_quant.AG <- as.data.frame(t(mQuantiles.AG)) %>% 
    rownames_to_column('Year') %>% 
    tidyr::gather(key = "Quantile", value = "Value", - Year) %>%
    mutate('Year' = as.numeric(Year), 'Type' = 'AG2020')
  df_quant      <- rbind(df_quant.new, df_quant.AG)
  
  # Observed and fitted data
  df_fit.new  <- data.frame('Year' = yv, "Quantile" = 'Fitted', 
                            "Value" = q.fit.new, 'Type' = 'New')
  df_fit.AG   <- data.frame('Year' = yv.AG, "Quantile" = 'Fitted', 
                            "Value" = q.mat.AG[1:length(yv.AG),1], 
                            'Type' = 'AG2020')
  df_fit      <- rbind(df_fit.new, df_fit.AG)
  
  df_obs <- data.frame('Year' = yv, "Quantile" = 'Obs', "Value" = q.obs, 
                       "Type" = NA)
  
  # Merge
  df <- rbind(rbind(df_quant, df_obs), df_fit)
  
  df$Type <- factor(df$Type, levels = c("New","AG2020"))
  
  df
}

