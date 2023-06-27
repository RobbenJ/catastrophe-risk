##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('doParallel','dplyr', 'foreach', 'ggplot2', 'ggpubr', 'mvtnorm',
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
load(file = "Data/data_M_final.RData")

### 1.8) Load functions and previous results ----
source("Functions/Functions projection mortality rates.R")
LL_M            <- readRDS(file = 'Save/LL_M_baseline.rds')
proj            <- readRDS(file = 'Save/Simulated_period_effects.rds')      
states.regime20 <- readRDS(file = 'Save/Sim_regime_process_20.rds')
states.regime60 <- readRDS(file = 'Save/Sim_regime_process_60.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1850:2021
yvSPEC  <- 1850:2021
TL      <- length(yv)
m       <- length(xv)
nSim    <- 10000
n.ahead <- 150
yv.proj <- (tail(yv,1) + 1):(tail(yv,1) + n.ahead)
yvv     <- yv[1]:(tail(yv.proj,1))
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")
CountrySPEC  <- "NED"

##### 2) Calculate fitted mortality rates baseline model -----

# Observed crude central death rates
mtx.obs <- data_M$UNI$NED$dtx/data_M$UNI$NED$etx

# Fitted death rates baseline improvement model (recursively)
mtx.fit     <- matrix(NA, nrow = TL, ncol = m, dimnames = list(yv, xv))
mtx.fit[1,] <- mtx.obs[1,]
for (i in 2:TL){
  mtx.fit[i,] <- mtx.obs[i-1,] * 
    exp(LL_M$NR$A.x + LL_M$NR$B.x * LL_M$NR$K.t[i-1] + 
          LL_M$NR$C.x * LL_M$NR$L.t[i-1] +
          LL_M$NR$b.x * LL_M$NR$k.t[i-1] + 
          LL_M$NR$c.x * LL_M$NR$l.t[i-1])
}
qtx.fit     <- 1 - exp(-mtx.fit)

##### 3) Three example trajectories of mortality rates -----
sim.vec <- c(1,2,3)
Ages    <- c(35,80)
simq20_ex <- project.mortality.rates(20:59, yv, LL_M, proj, n.ahead, sim.vec, 
                                     CountrySPEC, states.regime20$regime)
simq60_ex <- project.mortality.rates(60:85, yv, LL_M, proj, n.ahead, sim.vec, 
                                     CountrySPEC, states.regime60$regime)

dimnames(simq20_ex) <- list(20:59, yv[1]:(tail(yv,1) + n.ahead), 
                            1:length(sim.vec))
dimnames(simq60_ex) <- list(60:85, yv[1]:(tail(yv,1) + n.ahead), 
                            1:length(sim.vec))

df35 <- data.frame('Year' = yv[1]:(tail(yv,1) + n.ahead), 'Age' = Ages[1], 
                   'qxt' = as.numeric(simq20_ex[as.character(Ages[1]),,]), 
                   'Sim' = as.character(c(rep(1:length(sim.vec), 
                                              each = dim(simq20_ex)[2])))) %>%
  mutate('Type' = plyr::mapvalues(as.character(Year < 2022), 
                                  from = c('TRUE', 'FALSE'), 
                                  to = c('Observed','Forecast'))) %>%
  dplyr::filter(Sim %in% as.character(c(1:length(sim.vec))))

df80 <- data.frame('Year' = yv[1]:(tail(yv,1) + n.ahead), 'Age' = Ages[2], 
                   'qxt' = as.numeric(simq60_ex[as.character(Ages[2]),,]), 
                   'Sim' = as.character(c(rep(1:length(sim.vec), 
                                              each = dim(simq60_ex)[2])))) %>%
  mutate('Type' = plyr::mapvalues(as.character(Year < 2022), 
                                  from = c('TRUE', 'FALSE'), 
                                  to = c('Observed','Forecast'))) %>%
  dplyr::filter(Sim %in% as.character(1:length(sim.vec)))

df <- rbind(df35, df80) %>% dplyr::filter(Year >= 1950 & Year <= 2080)
df$Age <- factor(as.character(df$Age), levels = Ages, 
                 labels = paste('Age', Ages))

Traj_ages <- ggplot(df, aes(x = Year, y = qxt)) + 
  theme_bw(base_size = 15) + xlab('Year') + ylab(bquote(q['x,t'])) +
  facet_wrap(~Age, scales = 'free') + 
  geom_point(data = df %>% dplyr::filter(Type == 'Observed'), col = col.dot) + 
  geom_line(data = df %>% dplyr::filter(Type == 'Forecast' | Year == 2021), 
            aes(group = Sim, col = Sim), alpha = 0.7,linewidth = 1)  +
  scale_color_manual(name = '', values = c('1' = '#A3D596', '2' = '#ff993d', 
                                           '3' = '#F8766D'),
                     labels = c('1' = 'Path 1', '2' = 'Path 2', 
                                '3' = 'Path 3')) + 
  geom_vline(xintercept = 2021, linetype = 'dashed', col = 'gray70') +
  theme(strip.background = element_blank(),
        legend.position = 'bottom') + 
  guides(col = guide_legend(override.aes = list(alpha = 1)))

Traj_ages

##### 4) 10 000 trajectories of future mortality rates ----
simq20 <- project.mortality.rates(20:59, yv, LL_M, proj, n.ahead, 1:nSim, 
                                  CountrySPEC, states.regime20$regime)
simq60 <- project.mortality.rates(60:85, yv, LL_M, proj, n.ahead, 1:nSim, 
                                  CountrySPEC, states.regime60$regime)

dimnames(simq20) <- list(20:59, yv[1]:(tail(yv,1) + n.ahead), 1:nSim)
dimnames(simq60) <- list(60:85, yv[1]:(tail(yv,1) + n.ahead), 1:nSim)

# Save
# saveRDS(simq20, file = 'Save/simq20_rs.rds')
# saveRDS(simq60, file = 'Save/simq60_rs.rds')

#### 5) Predictions from AG2020 model -----
simq_old <- readRDS(file = "Data/projNED.rds")
dimnames(simq_old$Male) <- list(0:90, 1970:2080, 1:nSim)

#### 6) Plot projected mortality rates for ages 35, 50, 65, 80 ----
df_35 <- plotdf_qxt(xv, yv, data_M, 35, simq20, simq_old, "Male", "NED") %>% 
  mutate(Age = 'Age 35')
df_50 <- plotdf_qxt(xv, yv, data_M, 50, simq20, simq_old, "Male", "NED") %>% 
  mutate(Age = 'Age 50')
df_65 <- plotdf_qxt(xv, yv, data_M, 65, simq60, simq_old, "Male", "NED") %>% 
  mutate(Age = 'Age 65')
df_80 <- plotdf_qxt(xv, yv, data_M, 80, simq60, simq_old, "Male", "NED") %>% 
  mutate(Age = 'Age 80')

df <- df_35 %>% bind_rows(df_50) %>% bind_rows(df_65) %>% bind_rows(df_80)
df <- df %>% dplyr::filter(Year >= 1950 & Year <= 2080)

df[which(df$Quantile %in% c("Q1","Q7") & df$Year > 2076),] <- NA
df[which(df$Quantile %in% c("Q2","Q6") & df$Year > 2078),] <- NA
df[which(df$Quantile %in% c("Q3","Q5") & df$Year > 2080),] <- NA
df[which(df$Quantile %in% c("Q1","Q2","Q6","Q7") & 
           as.character(df$Type) == 'AG2020'),] <- NA

df <- df %>% tidyr::drop_na(Year)

a1 <- ggplot(df) + theme_bw(base_size = 15) + ylab(bquote(q['x,t'])) + 
  facet_wrap(~Age, nrow = 2, ncol = 2, scales = 'free_y') + 
  theme(strip.background = element_blank()) +
  geom_point(data = df %>% filter(Quantile == 'Obs'), aes(x = Year, y = Value, col = 'Obs'), size = 1) +
  geom_ribbon(data = df %>% filter(Quantile %in% c('Q1')), aes(x = Year,
                                                               ymin = df %>% filter(Quantile == 'Q1') %>% pull(Value),
                                                               ymax = df %>% filter(Quantile == 'Q7') %>% pull(Value), 
                                                               fill = Type, group = Type, alpha = 'CI99')) +
  geom_ribbon(data = df %>% filter(Quantile %in% c('Q2')), aes(x = Year,
                                                               ymin = df %>% filter(Quantile == 'Q2') %>% pull(Value),
                                                               ymax = df %>% filter(Quantile == 'Q6') %>% pull(Value), 
                                                               fill = Type, group = Type, alpha = 'CI95')) +
  geom_ribbon(data = df %>% filter(Quantile %in% c('Q3')), aes(x = Year,
                                                               ymin = df %>% filter(Quantile == 'Q3') %>% pull(Value),
                                                               ymax = df %>% filter(Quantile == 'Q5') %>% pull(Value), 
                                                               fill = Type, group = Type, alpha = 'CI90')) +
  geom_line(data = df %>% filter(Quantile == 'Q4' & Type == 'New'), aes(x = Year, y = Value, group = Type), col = '#4590ba', size = 0.5) +
  geom_line(data = df %>% filter(Quantile == 'Q4' & Type == 'AG2020'), aes(x = Year, y = Value, group = Type), col = '#62805a', size = 0.5) +
  geom_line(data = df %>% filter(Quantile == 'Fitted'), aes(x = Year, y = Value, group = Type, col = Type), linetype = 'solid', linewidth = 0.5) +
  scale_fill_manual(name = '', values = c('New' = col.dot, 'AG2020' = col.line),
                    labels = c('New' = '90% CI AG2020', 'AG2020' = ''),
                    guide = guide_legend(override.aes = list(
                      fill = c('New' = 'white', 'AG2020' = col.line)))) + 
  scale_color_manual(name = '', values = c('Obs' = col.shade, 'New' = col.dot, 'AG2020' = col.line), 
                     labels = c('Obs' = 'Observed', 'New' = 'Fit RS model', 'AG2020' = 'Fit AG2020'),
                     guide = guide_legend(override.aes = list(
                       linetype = c('Obs' = 'solid', 'New' = 'solid', 'AG2020' = 'blank'),
                       shape = c('Obs' = NA, 'New' = NA, 'AG2020' = 16)))) + 
  scale_alpha_manual(name = 'RS model',
                     values = c('CI99' = 0.2, 'CI95' = 0.5, 'CI90' = 1),
                     labels = c('CI99' = '99% CI', 'CI95' = '95% CI', 'CI90' = '90% CI'),
                     guide = guide_legend(override.aes = list(
                       fill = c('CI99' = col.dot, 'CI95' = col.dot, 'CI90' = col.dot)))) + 
  theme(legend.position = 'bottom',
        legend.spacing.x = unit(0.1,'cm'))

a1
