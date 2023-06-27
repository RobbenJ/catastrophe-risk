##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('doParallel','dplyr', 'foreach', 'ggplot2', 'MASS', 'mvtnorm', 
              'parallel', 'rstudioapi' ,'stats','tibble','tidyr', 'utils')

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
load(file = "Data/data_M_final.RData")

### 1.8) Load functions and previous results ----
source("Functions/Functions projection baseline model.R")
LL_M <- readRDS(file = 'Save/LL_M_baseline.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1850:2021
TL      <- length(yv)
m       <- length(xv)
nSim    <- 10000
n.ahead <- 150
yvv     <- yv[1]:(tail(yv,1) + n.ahead)
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")

##### 2) Decay parameter ----

# data frame of period effects
df.orig <- cbind(LL_M$NR$K.t, LL_M$NR$L.t, LL_M$NR$k.t, LL_M$NR$l.t)
rownames(df.orig) <- yv[-1]

# Find optimal decay parameter by minimizing one-step ahead log-likelihood values
find_lambda <- FALSE
if(find_lambda) {
  lambda_opt <- opt.decay(df.full = df.orig, n.years = 121)
} else
  lambda_opt <- 0.943 

##### 3) Multivariate time series estimation ----

# Optimization on entire time series with optimal decay parameter
nll <- nlminb(start = c(0,0,1,0,0,0,1,0,0,1,0,1), objective = obj.logL, 
              gradient = NULL, hessian = NULL, lambda = lambda_opt, df = df.orig,
              control = list(trace = F, eval.max = 10000, iter.max = 10000, 
                             rel.tol = 1e-12))

# Estimated intercepts common period effects 
round(nll$par,5)[1:2]

# Estimated covariance matrix
mat <- matrix(0, nrow = 4, ncol = 4)
mat[lower.tri(mat,diag=T)] <- nll$par[-c(1:2)]
Sigma <- mat %*% t(mat)
round(Sigma,5)

##### 4) Generate future trajectories of period effects -----

# Create objects to store trajectories
output <- list()
output[[1]] = output[[2]] = output[[3]] = output[[4]]  <- 
  matrix(NA, ncol = nSim, nrow = length(yv) + n.ahead - 1)
output[[1]][1:(length(yv)-1),] <- LL_M$NR$K.t 
output[[2]][1:(length(yv)-1),] <- LL_M$NR$L.t
output[[3]][1:(length(yv)-1),] <- LL_M$NR$k.t
output[[4]][1:(length(yv)-1),] <- LL_M$NR$l.t

# Generate trajectories
for(i in 1:nSim){
  errors = t(MASS::mvrnorm(n = n.ahead, mu = c(nll$par[1:2],0,0), Sigma = Sigma))
  
  for(j in 1:length(output)){
    output[[j]][length(yv):(nrow(output[[j]])),i] = errors[j,]
  }
}

# Store everything in list
proj <- list()
proj$K.t_M <- output[[1]]
proj$L.t_M <- output[[2]]
proj$k.t_M <- output[[3]]
proj$l.t_M <- output[[4]]
proj$C     <- Sigma

# Save the simulations of the period effects
# saveRDS(proj, file = "Save/Simulated_period_effects.rds")
proj <- readRDS(file = "Save/Simulated_period_effects.rds")

##### 5) Plot fan chart of projected period effects ----
TS.list <- list('KtM' = LL_M$NR$K.t, 'LtM' = LL_M$NR$L.t,
                'ktM' = LL_M$NR$k.t, 'ktM' = LL_M$NR$l.t)

dfKtM <- construct_plotdf(yv, proj, TS.list, "K.t_M") %>% 
  mutate(Type = 'K[plain(t)]^plain((1))') 
dfLtM <- construct_plotdf(yv, proj, TS.list, "L.t_M") %>% 
  mutate(Type = 'K[plain(t)]^plain((2))')
dfktM <- construct_plotdf(yv, proj, TS.list, "k.t_M") %>% 
  mutate(Type = 'kappa[plain(t)]^plain((1))')
dfltM <- construct_plotdf(yv, proj, TS.list, "l.t_M") %>% 
  mutate(Type = 'kappa[plain(t)]^plain((2))')

sample.KtM <- data.frame('Year' = yvv[yvv > 2021], 'Quantile' = 'Sample1', 
                         'Value' = proj$K.t_M[which(yvv[-1] > 2021),1], 
                         'Type' = 'K[plain(t)]^plain((1))')
sample.LtM <- data.frame('Year' = yvv[yvv > 2021], 'Quantile' = 'Sample1', 
                         'Value' = proj$L.t_M[which(yvv[-1] > 2021),1], 
                         'Type' = 'K[plain(t)]^plain((2))')
sample.ktM <- data.frame('Year' = yvv[yvv > 2021], 'Quantile' = 'Sample1', 
                         'Value' = proj$k.t_M[which(yvv[-1] > 2021),1], 
                         'Type' = 'kappa[plain(t)]^plain((1))')
sample.ltM <- data.frame('Year' = yvv[yvv > 2021], 'Quantile' = 'Sample1', 
                         'Value' = proj$l.t_M[which(yvv[-1] > 2021),1], 
                         'Type' = 'kappa[plain(t)]^plain((2))')

df_combined <- do.call(rbind, list(dfKtM, dfLtM, dfktM, dfltM, sample.KtM, 
                                   sample.LtM, sample.ktM, sample.ltM)) %>%
  dplyr::filter(Year <= 2080)

SimKt <- ggplot(df_combined) + 
  theme_bw(base_size = 15) + 
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_ribbon(data = df_combined %>% filter(Quantile %in% c('Q1')), 
              aes(x = Year, ymin = df_combined %>% filter(Quantile == 'Q1') %>% 
                    pull(Value), ymax = df_combined %>% 
                    filter(Quantile == 'Q3') %>% pull(Value),  fill = 'CI'),
              alpha = 1) + 
  geom_line(data = df_combined %>% filter(Quantile == 'Q2'), 
            aes(x = Year, y = Value, col = 'Median')) +
  geom_point(data = df_combined %>% filter(Quantile == 'InSampleFit'), 
             aes(x = Year, y = Value, col = 'ISF'), size = 2) +
  geom_line(data = df_combined %>% filter(Quantile == 'Sample1'), 
            aes(x = Year, y = Value), col = col.line, linewidth = 0.4, 
            show.legend = FALSE) + 
  ylab('Period effects') + theme(legend.position = 'bottom') + 
  scale_color_manual('', values = c('ISF' = col.dot, 'Median' = col.dot),
                     labels = c('ISF' = 'Calibrated time component', 
                                'Median' = '50% quantile'),
                     guide = guide_legend(override.aes = list(
                       linetype = c('ISF' = 'blank', 'Median' = 'solid'),
                       shape = c('ISF' = 16, 'Median' = NA)))) + 
  scale_fill_manual('', values = c('CI' = col.shade), 
                    labels = c('CI' = '99% CI')) + 
  theme(strip.background = element_blank())

SimKt

