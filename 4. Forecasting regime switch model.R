##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('doParallel','dplyr', 'foreach', 'ggplot2', 'ggpubr', 
              'mvtnorm', 'parallel', 'rstudioapi', 'tibble', 'utils')

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
source("Functions/Functions projection regime switch.R")
list.param20 <- readRDS(file = 'Save/param20-59.rds')
list.param60 <- readRDS(file = 'Save/param60-85.rds')
ztx.old      <- readRDS(file = 'Save/ztx.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1850:2021
yvSPEC  <- 1850:2021
TL      <- length(yv)
m       <- length(xv)
nSim    <- 10000
n.ahead <- 150
yv.pred <- (yv[length(yv)]+1):(yv[length(yv)] + n.ahead)
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")

##### 2) Generate future trajectories for the regime switching process ----

# Age category 20-59
list2env(list.param20, globalenv())
states.regime20 <- SimulateRegimeSwitch(20:59, yv, p12, p21, C, muH, sigmaH, 
                                        sigmaE2, slope2, nSim, n.ahead)

# Age category 60-85
list2env(list.param60, globalenv())
states.regime60 <- SimulateRegimeSwitch(60:85, yv, p12, p21, C, muH, sigmaH, 
                                        sigmaE2, slope2, nSim, n.ahead)

# Set row- and column names 
dimnames(states.regime20$states) = dimnames(states.regime60$states) <- 
  list(1:nSim, yv.pred)
dimnames(states.regime20$regime) <- list(1:nSim, 20:59, yv.pred)
dimnames(states.regime60$regime) <- list(1:nSim, 60:85, yv.pred)

# Save results
# saveRDS(states.regime20, file = 'Save/Sim_regime_process_20.rds')
# saveRDS(states.regime60, file = 'Save/Sim_regime_process_60.rds')

##### 3) Plot one example trajectory ---- 

# Read stored results
states.regime20 <- readRDS(file = 'Save/Sim_regime_process_20.rds')
states.regime60 <- readRDS(file = 'Save/Sim_regime_process_60.rds')

# Plot trajectory Markov chain and regime switch component
j <- 1
df_exMarkov <- data.frame('Year' = yv.pred, 'Markov' = 
                            states.regime60$states[j,]) %>% 
  dplyr::filter(Year <= 2080)
df_exRegime <- data.frame('Year' = yv.pred, 'RegimeSwitch' = 
                            colMeans(states.regime60$regime[j,,])) %>% 
  dplyr::filter(Year <= 2080)

pos <- rle(df_exMarkov$Markov)
Length_x <- rep(FALSE, length(df_exMarkov$Year))
Length_x[sequence(pos$lengths)==1] <- pos$values==2
Years   <- df_exMarkov$Year[Length_x]
Lengths <- pos$lengths[pos$values == 2]

pexMarkov <- ggplot(df_exMarkov, aes(x = Year, y = Markov)) + 
  theme_bw(base_size = 15) + ylab("Markov chain") + 
  scale_y_continuous(breaks = c(1,2), label = c('LVS','HVS'))

pexRegime <- ggplot(df_exRegime, aes(x = Year, y = RegimeSwitch)) + 
  theme_bw(base_size = 15) + ylab("Regime Switch component") 

for (t in 1:length(Years)){
  pexMarkov <- pexMarkov + 
    geom_rect(aes(xmin = Years[!!t] - 1, xmax = Years[!!t] + Lengths[!!t] + 0.5, ymin = -50, ymax = 50),
              fill = col.shade)
  
  pexRegime <- pexRegime + 
    geom_rect(aes(xmin = Years[!!t] - 1, xmax = Years[!!t] + Lengths[!!t] + 0.5, ymin = -50, ymax = 50),
              fill = col.shade) 
  
  if(t == length(Years)) {
    pexMarkov <- pexMarkov +   
      geom_line(col = col.dot, lwd = 1.25) +
      coord_cartesian(ylim = c(1,2))
    
    pexRegime <- pexRegime + 
      geom_line(col = col.dot, lwd = 1.25) +
      coord_cartesian(ylim = range(df_exRegime$RegimeSwitch))
  }
  
}

pExSim <- ggpubr::ggarrange(pexMarkov, pexRegime)
pExSim
