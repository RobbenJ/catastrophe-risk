##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('abind', 'doParallel','dplyr', 'foreach', 'ggplot2', 'ggpubr',
              'parallel', 'remotes', 'rstudioapi', 'tibble', 'stats', 'tidyr', 
              'utils')

### 1.3) Install packages ----
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# install.packages("BiocManager")
# BiocManager::install("aroma.light")
# remotes::install_github("RobbenJ/MultiMoMo")

### 1.4) Load Packages ----
invisible(lapply(c(packages,'MultiMoMo'), library, character.only = TRUE))

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
simq20 <- readRDS(file = 'Save/simq_noww1875_20_rs.rds')
simq60 <- readRDS(file = 'Save/simq_noww1875_60_rs.rds')
q_BE_closed <- readRDS(file = 'Save/qxt_BE_closed1875.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1875:2021
TL      <- length(yv)
m       <- length(xv)
n.ahead <- 59
yv.proj <- (yv[length(yv)]+1):(yv[length(yv)] + n.ahead)
yvv     <- yv[1]:(tail(yv.proj,1))
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")
CountrySPEC  <- "NED"

##### 2) Close simulated mortality rates -----
# simq   <- abind(simq20, simq60, along = 1)
# simq_closed <- MultiMoMo::close_mortality_rates(yv, simq, 35, 11, parallel = TRUE)
# saveRDS(simq_closed, file = 'Save/simq_noww1875_closed.rds')
simq_closed <- readRDS(file = 'Save/simq_noww1875_closed.rds')

##### 3) Calculation of SCRs -----

#### Specificiations
agem  <- 120
xmin  <- 20
sty   <- 2021
int   <- 0.02
v     <- 1/(1+int)

### 3.1) life annuity - 65 yo, annuity of 1000 annually

# Specifications of life annuity
age.la.vec <- 55:90
df.ILA <- data.frame('Age' = age.la.vec, 'SCR.std'= NA, "SCR.var"= NA, "BEL0" = NA)

for(age.la in age.la.vec){
  TL       <- agem - age.la + 1
  annuity  <- 10000
  index.la <- cbind(as.character(age.la:agem), as.character(sty:(sty + agem - age.la)))
  
  # BEL at time 0 using BE mortality rates for life annuity 
  qxt.BE.la  <- q_BE_closed[index.la]
  pxt.BE.la  <- 1 - qxt.BE.la
  Tpxt.BE.la <- cumprod(pxt.BE.la)
  BEL.0.la   <- sum(annuity*v^(1:TL)*Tpxt.BE.la)
  
  # BEL at time 0 using BE mortality rates | longevity shock of 20% (QIS5)
  qxt.BE.long  <- qxt.BE.la*0.8
  pxt.BE.long  <- 1 - qxt.BE.long
  Tpxt.BE.long <- cumprod(pxt.BE.long)
  BEL.0.long   <- sum(annuity*v^(1:TL)*Tpxt.BE.long)
  
  SCR.long <- BEL.0.long - BEL.0.la
  SCR.long
  
  # BEL using simulated mortality rates - life annuity
  BEL.0.s <- rep(NA, dim(simq_closed)[3])
  
  for(s in 1:dim(simq_closed)[3]){
    qxt.s        <- simq_closed[,,s]
    qxt.s        <- qxt.s[index.la]
    pxt.s        <- 1 - qxt.s
    Tpxt.s       <- cumprod(pxt.s)
    BEL.0.s[s]   <- sum(annuity*v^(1:TL-1)*Tpxt.s)
  }
  
  SCR.la <- quantile(BEL.0.s, 0.995) - BEL.0.la
  
  df.ILA[which(age.la.vec %in% age.la),2:4] <- c(SCR.long, SCR.la, BEL.0.la)
}

ILA <- ggplot(df.ILA) + 
  xlab('Age') + ylab('SCR') + ggtitle('Immediate life annuity') + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Age, y = SCR.var, col = 'v1'), linewidth = 1.5, show.legend = TRUE) +
  geom_line(aes(x = Age, y = SCR.std, col = 'v2'), linewidth = 1.5, show.legend = TRUE) +
  scale_color_manual(name = '',
                     labels = c('v1' = 'VaR approach', 'v2' = 'Standard model'),
                     values = c('v1' = col.dot, 'v2' = col.shade)) + 
  theme(legend.position = 'bottom')
ILA


### 3.2) Term insurance ----

# Specifications
death.ben <- 150000
age.ti.vec    <- 20:60
df.TLI <- data.frame('Age' = age.ti.vec, 'SCR.mort'= NA, 'SCR.cat'= NA, 'SCR.std'= NA, "SCR.var"= NA, "BEL0" = NA)

for (age.ti in age.ti.vec){
  len       <- 65 - age.ti
  index.ti  <-  cbind(as.character(age.ti:(age.ti + len)), as.character(sty:(sty + len)))
  
  # BEL at time 0 using BE mortality rates for term insurance
  qxt.BE.ti  <- q_BE_closed[index.ti]
  pxt.BE.ti  <- 1 - qxt.BE.ti
  Tpxt.BE.ti <- cumprod(pxt.BE.ti)
  BEL.0.ti   <- sum(death.ben*v^(1:len)*c(1,Tpxt.BE.ti[1:(len-1)])*qxt.BE.ti[1:len])
  
  # BEL at time 0 using BE mortality rates | mortality shock of 15% (QIS5)
  qxt.BE.mort  <- q_BE_closed[index.ti]*1.15
  pxt.BE.mort  <- 1 - qxt.BE.mort
  Tpxt.BE.mort <- cumprod(pxt.BE.mort)
  BEL.0.mort   <- sum(death.ben*v^(1:len)*c(1,Tpxt.BE.mort[1:(len-1)])*qxt.BE.mort[1:len]) 
  SCR.mort     <- BEL.0.mort - BEL.0.ti
  
  # BEL at time 0 using BE mortality rates | catastrophe stock - absolute 0.0015 increase
  qxt.BE.cat  <- q_BE_closed[index.ti] + c(0.0015,rep(0,length(index.ti) - 1))
  pxt.BE.cat  <- 1 - qxt.BE.cat
  Tpxt.BE.cat <- cumprod(pxt.BE.cat)
  BEL.0.cat   <- sum(death.ben*v^(1:len)*c(1,Tpxt.BE.cat[1:(len-1)])*qxt.BE.cat[1:len]) 
  SCR.cat     <- BEL.0.cat - BEL.0.ti
  
  # BEL using simulated mortality rates - term life insurance
  BEL.0.s <- rep(NA, dim(simq_closed)[3])

  for(s in 1:dim(simq_closed)[3]){
    qxt.s         <- simq_closed[,,s]
    qxt.s         <- qxt.s[index.ti]
    pxt.s         <- 1 - qxt.s
    Tpxt.s        <- cumprod(pxt.s)
    BEL.0.s[s]    <- sum(death.ben*v^(1:len)*c(1,Tpxt.s[1:(len-1)])*qxt.s[1:len])
  }
  
  SCR.ti <- quantile(BEL.0.s, 0.995) - BEL.0.ti
  
  Age  <- age.ti
  mort <- round(SCR.mort,2)
  cat  <- 0.0015*(150000 - BEL.0.ti)
  scr  <- round(sqrt(mort^2 + cat^2),2)

  df.TLI[which(Age == df.TLI$Age),] <- c(Age, mort, cat, scr, SCR.ti, BEL.0.ti)
}

### 3.3) Plot
TLI <- ggplot(df.TLI) + 
  xlab('Age') + ylab('SCR') + ggtitle('Term life insurance') + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Age, y = SCR.var, col = 'v1'), linewidth = 1.5, show.legend = TRUE) +
  geom_line(aes(x = Age, y = SCR.std, col = 'v2'), linewidth = 1.5, show.legend = TRUE) +
  scale_color_manual(name = '',
                     labels = c('v1' = 'VaR approach', 'v2' = 'Standard model'),
                     values = c('v1' = col.dot, 'v2' = col.shade)) + 
  theme(legend.position = 'bottom')


# Save plot
p <- ggpubr::ggarrange(ILA, TLI, common.legend = TRUE, legend = 'bottom')


