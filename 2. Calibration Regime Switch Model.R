##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('DEoptimR', 'dplyr', 'ggplot2', 'mvtnorm', 'plyr', 'rstudioapi',
              'showtext', 'sysfonts', 'tidyr','utils')

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
source("Functions/Functions calibration regime switch.R")
LL_M <- readRDS(file = 'Save/LL_M_baseline1875.rds')

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1875:2021
TL      <- length(yv)
m       <- length(xv)
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")

##### 2) Baseline model's residuals ----

# Create empty matrix
ztx <- matrix(NA, nrow = TL - 1, ncol = m)
dimnames(ztx) <- list(yv[-1], xv)

# Crude Dutch central death rates (observed)
mtx.obs <- data_M$UNI$NED$dtx/data_M$UNI$NED$etx

# Baseline model's residuals
for(t in 1:(TL - 1)){
  ztx[t,] <- log(mtx.obs[t + 1,]) - log(mtx.obs[t,]) - LL_M$NR$A.x - 
    LL_M$NR$B.x * LL_M$NR$K.t[t] - LL_M$NR$C.x * LL_M$NR$L.t[t] - 
    LL_M$NR$b.x * LL_M$NR$k.t[t] - LL_M$NR$c.x * LL_M$NR$l.t[t]
}
dimnames(ztx) <- list(yv[-1], xv)

# Save
# saveRDS(ztx, 'Save/ztx.rds')

# Plot ztx
df.long <- as.data.frame(t(ztx)) %>% mutate('Age' = 20:85) %>% 
  tidyr::gather(key = 'Year', value = 'Errors', - Age)
df.long$Year <- as.numeric(df.long$Year)

p.ztx.path <- ggplot(df.long, aes(x = Year, y = Errors)) +
  theme_bw(base_size = 15) + ylab(bquote(z['x,t'])) + 
  geom_point(aes(colour = Age, group = Age), size = 0.8) + 
  scale_color_gradient(low = col.dot, high =  col.dot.shade) 

p.ztx.path

##### 3) Calibration of regime switching model ----

# Load baseline model's residuals
rm(ztx)
ztx.old <- readRDS(file = "Save/ztx1875.rds")

# Calibration regime switch ages 20-59
ztx1    <- ztx.old[,as.character(20:59)]
init1.1 <- c('p12' = 0.1, 'p21' = 0.5, 'muH' = 0, 'sigmaE1' = 0.1, 
             'slope1' = 0.01, 'sigmaE2' = 0.1, 'slope2' = 0.01, 'sigmaH' = 1)
init2.1 <- rep(1,ncol(ztx1))/sqrt(ncol(ztx1))

opt1    <- optim.rs(ztx = ztx1, init1 = init1.1, init2 = init2.1,
                    ww.exclude = FALSE)

# Calibration regime switch ages 60-85
ztx2    <- ztx.old[,as.character(60:85)]
init1.2 <- c('p12' = 0.1, 'p21' = 0.5, 'muH' = 0, 'sigmaE1' = 0.1, 
             'slope1' = 0.01, 'sigmaE2' = 0.1, 'slope2' = 0.01, 'sigmaH' = 1)
init2.2 <- rep(1,ncol(ztx2))/sqrt(ncol(ztx2))

opt2    <- optim.rs(ztx = ztx2, init1 = init1.2, init2 = init2.2,
                    ww.exclude = FALSE)

### 3.1) Save calibrated parameters ----

# Age category 20-59
# file1 <- list('p12' = opt1$param1['p12'], 'p21' = opt1$param1['p21'],
#               'muH' = opt1$param1['muH'], 'sigmaE1' = opt1$param1['sigmaE1'],
#               'slope1' = opt1$param1['slope1'], 'sigmaE2' = opt1$param1['sigmaE2'],
#               'slope2' = opt1$param1['slope2'], 'sigmaH' = opt1$param1['sigmaH'],
#               'C' = opt1$param2, 'LogLik' = opt1$LogLik, 'Iter' = opt1$Iter)
# saveRDS(file1, file = 'Save/param20-59.rds')

# Age category 60-85
# file2 <- list('p12' = opt2$param1['p12'], 'p21' = opt2$param1['p21'], 
#               'muH' = opt2$param1['muH'], 'sigmaE1' = opt2$param1['sigmaE1'],
#               'slope1' = opt2$param1['slope1'], 'sigmaE2' = opt2$param1['sigmaE2'], 
#               'slope2' = opt2$param1['slope2'], 'sigmaH' = opt2$param1['sigmaH'], 
#               'C' = opt2$param2, 'LogLik' = opt2$LogLik, 'Iter' = opt2$Iter)
# saveRDS(file2, file = 'Save/param60-85.rds')

### 3.2) Conditional probabilities ----

# Run the functions obj.param1 and obj.param2 with optimal parameters

# Age category 20-59
# df_cprob2 <- data.frame('Year' = yv[-1], 'cprob' = cprob[1,], 'Avg' = rowMeans(ztx1)) %>%
#   mutate('Avg.c' = (Avg - min(Avg))/(max(Avg) - min(Avg)), 'Age_cat' = "20-59")
# file2 <- "Save/Cond.prob.20-59.rds"
# saveRDS(df_cprob2, file = file2)

# Age category 60-85
# df_cprob1 <- data.frame('Year' = yv[-1], 'cprob' = cprob[1,], 'Avg' = rowMeans(ztx2)) %>%
#   mutate('Avg.c' = (Avg - min(Avg))/(max(Avg) - min(Avg)), 'Age_cat' = "60-85")
# file1 <- "Save/Cond.prob.60-85.rds"
# saveRDS(df_cprob1, file = file1)

##### 4) Visualisations ----

### 4.1) Read in calibration results ----
list.param20 <- readRDS(file = 'Save/param.noww1875.20-59.rds')
list.param60 <- readRDS(file = 'Save/param.noww1875.60-85.rds')
df_cprob20   <- readRDS(file = 'Save/Cond.prob.noww1875.20-59.rds')
df_cprob60   <- readRDS(file = 'Save/Cond.prob.noww1875.60-85.rds')

### 4.2) Plot age patterns of mortality shocks ----

# Age category 20-59
df_C <- data.frame('Age' = 20:59, 'C' = list.param20$C)
pC <- ggplot(df_C, aes(x = Age, y = C)) + 
  theme_bw(base_size = 15) + ylab(bquote('C'[x])) + xlab('') + 
  geom_point(col = col.dot, size = 3)
pC

# Age category 60-85
df_C <- data.frame('Age' = 60:85, 'C' = list.param60$C)
pC <- ggplot(df_C, aes(x = Age, y = C)) + 
  theme_bw(base_size = 15) + ylab(bquote('C'[x])) + xlab('') + 
  geom_point(col = col.dot, size = 3)
pC

# Plot both age groups together using frak font
font_add_google(name = 'UnifrakturMaguntia','fraktur')

showtext_auto()

df_C1 <- data.frame('Age' = 20:59, 'C' = c(list.param20$C)) %>% 
  mutate('Type' = 'Age group 20-59')
df_C2 <- data.frame('Age' = 60:85, 'C' = c(list.param60$C)) %>% 
  mutate('Type' = 'Age group 60-85')
df_C  <- rbind(df_C1, df_C2)

pC1 <- ggplot(df_C) +
  facet_wrap(~Type, scales = 'free_x') + 
  theme_bw(base_size = 15) + ylab('B') + xlab('Age (x)') + 
  geom_point(aes(x = Age, y = C), col = col.dot, size = 2) + 
  theme(axis.title.y = element_text(family = "fraktur", size = 20),
        title = element_text(family="sans"),
        strip.text.x.top = element_text(family = 'sans'),
        strip.background = element_blank())  
pC1 

showtext_auto(FALSE)

### 4.3) Plot conditional probability of being in regime state 1 ----
df_cprob <- rbind(df_cprob20, df_cprob60)
df_cprob$Age_cat <- plyr::mapvalues(df_cprob$Age_cat, unique(df_cprob$Age_cat),
                                    c('Age group 20-59', 'Age group 60-85'))

p.cond <- ggplot(df_cprob) + 
  theme_bw(base_size = 15) + xlab('Year') + 
  facet_wrap(~Age_cat, nrow = 1, ncol = 2, scales = 'free_y') + 
  geom_point(aes(x = Year, y = Avg.c, col = 'Avg'), size = 1.5, colour = col.dot, show.legend = TRUE) +   
  geom_line(aes(x = Year, y = cprob, col = 'cprob'), lwd = 0.8, colour = col.line, show.legend = TRUE) +
  scale_color_manual(values = c('Avg' = col.line, 'cprob' = col.dot),
                     guide = guide_legend(override.aes = list(linetype = c('Avg' = 'blank', 'cprob' = 'solid'),
                                                              shape = c('Avg' = 16, 'cprob' = NA))),
                     name = '', 
                     labels = c('Avg' = bquote('MinMax'[t]*'(AVG'[x]*'(z'['x,t']*'))'), 'cprob' = 'Cond. prob.')) + 
  theme(strip.background = element_blank(),
        legend.position = 'bottom') + 
  ylab(bquote('Conditional '~ P(rho[t]~'= 1')))

p.cond  
