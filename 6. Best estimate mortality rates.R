##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('doParallel','dplyr', 'foreach', 'ggplot2', 'ggpubr', 'MortCast',
              'parallel', 'rstudioapi', 'tibble', 'stats', 'tidyr', 'utils')

### 1.3) Install packages ----
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# install.packages("BiocManager")
# BiocManager::install("aroma.light")

### 1.4) Load Packages ----
invisible(lapply(packages, library, character.only = TRUE))

### 1.5) Define colour scale ----
col.dot       <- '#56B4E9'
col.dot.shade <- '#cce9f8'
col.shade     <- '#ddf0fb'
col.line      <- '#A3D596'
col.outl      <- '#F8766D'

### 1.6) Set working directory ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### 1.7) Read in data set  ----
load(file = "Data/data_M_final1875.RData")

### 1.8) Load functions and previous results ----
source('Functions/Functions calibration baseline model.R')
source('Functions/Functions projection baseline model.R')
source('Functions/Functions projection mortality rates.R')
LL_M            <- readRDS(file = 'Save/LL_M_baseline1875.rds')
states.regime20 <- readRDS(file = 'Save/Sim_regime_process_noww1875_20.rds')
states.regime60 <- readRDS(file = 'Save/Sim_regime_process_noww1875_60.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1875:2021
TL      <- length(yv)
m       <- length(xv)
n.ahead <- 150
yv.proj <- (yv[length(yv)]+1):(yv[length(yv)] + n.ahead)
yvv     <- yv[1]:(tail(yv.proj,1))
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")
CountrySPEC  <- "NED"

##### 2) Estimation of multivariate time series model ---- 

# Optimal decay parameter (optimized in file 3.)
lambda_opt <- 0.942

# Data frame of period effects
df  <- cbind(LL_M$NR$K.t, LL_M$NR$L.t, LL_M$NR$k.t, LL_M$NR$l.t)
dimnames(df) <- list(yv[-1], 1:4)

# Estimate multivariate time series
nll <- nlminb(start = c(0,0,1,0,0,0,1,0,0,1,0,1), objective = obj.logL, 
              gradient = NULL, hessian = NULL, lambda = lambda_opt, df = df,
              control = list(trace = F, eval.max = 10000, iter.max = 10000, 
                             rel.tol = 1e-12))

# Estimated intercepts
round(nll$par,5)[1:2]

# Covariance matrix
mat <- matrix(0, nrow = 4, ncol = 4)
mat[lower.tri(mat,diag=T)] <- nll$par[-c(1:2)]
Sigma <- mat %*% t(mat)
round(Sigma,5)

##### 3) Best-estimate projections of multivariate time series model ---- 
# Empty object
output <- matrix(0, ncol = 4, nrow = length(yvv) -1)

# Calibrated parameters
output[1:(TL-1),] <- df
output[TL:(nrow(output)),1:2] <- cbind(rep(nll$par[1], n.ahead), 
                                       rep(nll$par[2], n.ahead))
dimnames(output) <- list(yvv[-1], c('KtM','LtM','ktM','ltM'))

# Store in proj list
proj <- list()
proj$K.t_M <- matrix(output[,'KtM'], ncol = 1)
proj$L.t_M <- matrix(output[,'LtM'], ncol = 1)
proj$k.t_M <- matrix(output[,'ktM'], ncol = 1)
proj$l.t_M <- matrix(output[,'ltM'], ncol = 1)
proj$C     <- Sigma

##### 4) Best-estimate projections using baseline improvement model -----

# Observed crude central death rates
mtx.obs.c <- data_M$UNI[[CountrySPEC]]$dtx/data_M$UNI[[CountrySPEC]]$etx

# Create regime object to offset the COVID-19 shock
regime       <- array(0, dim = c(2,m,n.ahead))
regime[1,,1] <- c(states.regime20$regime[1,,'2022'], 
                  states.regime60$regime[1,,'2022'])

# Projected best-estimate mortality rates
qProj_BE <- project.mortality.rates(xv, yv, LL_M, proj, n.ahead, c(1,1), 
                                    CountrySPEC, regime)[,,1]
dimnames(qProj_BE) <- list(xv, yvv)

# Example
plot(yvv[yvv > 1950], qProj_BE['80',yvv > 1950])

##### 5) Close BE mortality rates -----
mxt_BE        <- -log(1-qProj_BE)
mxt_BE_closed <- MortCast::kannisto(mxt_BE, est.ages = seq(75,85,1), proj.ages = seq(86,120,1))
qxt_BE_closed <- 1 - exp(-mxt_BE_closed)

# Save results
# saveRDS(qxt_BE_closed, file = "Save/qxt_BE_closed1875.rds")
