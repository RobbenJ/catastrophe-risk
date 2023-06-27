##### 1) Preliminary set-up ----

### 1.1) Remove history ----
rm(list = ls())

### 1.2) Package names ----
packages <- c('car', 'dplyr','ggplot2', 'ggpubr', 'ggrepel', 'grid', 'rrcov', 
              'rstudioapi', 'stats', 'tidyr', 'utils')

### 1.3) Install packages ----
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# install.packages("BiocManager")
# BiocManager::install("aroma.light")

### 1.4) Load Packages ----
invisible(lapply(c(packages,'aroma.light'), library, character.only = TRUE))

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
source("Functions/Functions calibration baseline model.R")

### 1.9) General specifications ----
xv      <- 20:85
yv      <- 1850:2021
yvSPEC  <- 1850:2021
TL      <- length(yv)
m       <- length(xv)
Country <- c("AUS", "BEL", "DNK", "FIN", "FRA", "GER", "ICE", "IRE", "LUX", 
             "NED", "NOR", "SWE", "SWI", "UNK")

##### 2) Fit baseline mortality improvement model ----

### 2.1) Indirect estimation on complete calibration period ----

# Fit
LL_M <- FitLiLeeNR(xv = xv, yv = yv, Data = data_M, Ax = TRUE,
                   exclude.coh = FALSE, CountrySPEC = "NED", yvSPEC = yvSPEC, 
                   fit1 = fit701C, fit2 = fit701C0)
KtM  <- LL_M$NR$K.t
LtM  <- LL_M$NR$L.t

# Plot calibrated parameters common trend model
df0x <- data.frame('Age' = xv, 'Ax' = LL_M$NR$A.x, 'Bx' = LL_M$NR$B.x,
                   'Cx' = LL_M$NR$C.x) %>% 
  tidyr::gather('Type', 'Age effect', - Age)
df0t <- data.frame('Year' = yv, 'Kt' = KtM, 'Lt' = LtM) %>% 
  tidyr::gather('Type', 'Period effect', - Year)

df0x$Type <- factor(df0x$Type, levels = c('Ax', 'Bx', 'Cx'), 
                    labels = c('A[plain(x)]','B[plain(x)]^plain((1))',
                               'B[plain(x)]^plain((2))'))
df0t$Type <- factor(df0t$Type, levels = c('Kt', 'Lt'), 
                    labels = c('L[plain(t)]^plain((1))', 
                               'L[plain(t)]^plain((2))'))

q0 <- ggplot(df0x, aes(x = Age, y = `Age effect`)) + 
  theme_bw(base_size = 15) +
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

q1 <- ggplot(df0t, aes(x = Year, y = `Period effect`)) + 
  theme_bw(base_size = 15) + 
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

q01 <- ggpubr::ggarrange(q0,q1, nrow = 2)
q01

### 2.2) Robust smoothing splines ----

# Zero weights for shocks in timeline (see paper)
outliers.timeline <- c(1854:1855, 1859, 1866, 1870:1871, 1889:1892, 1914:1919, 
                       1940:1945, 2020:2021)
w.spline          <- as.numeric(! yv %in% outliers.timeline)

# Fit robust splines to data
fit.Lt1 <- aroma.light::robustSmoothSpline(yv, KtM, cv = FALSE, w = w.spline, 
                                           minIter = 10, maxIter = 100)
fit.Lt2 <- aroma.light::robustSmoothSpline(yv, LtM, cv = FALSE, w = w.spline, 
                                           minIter = 10, maxIter = 100)

# Plot robust splines 
df.Lt1 <- data.frame('Year' = yv, 'PeriodEffect' = KtM, 'Spline' = fit.Lt1$y, 
                     Type = 'L[plain(t)]^plain((1))')
df.Lt2 <- data.frame('Year' = yv, 'PeriodEffect' = LtM, 'Spline' = fit.Lt2$y, 
                     Type = 'L[plain(t)]^plain((2))')
df     <- rbind(df.Lt1, df.Lt2)

l1 <- round(fit.Lt1$lambda, 16)
l2 <- round(fit.Lt2$lambda, 22)

final_spline <- ggplot(df, aes(x = Year)) +
  theme_bw(base_size = 18) + ylab("Period effect") +
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) +
  geom_point(aes(y = PeriodEffect, col = 'T1'), size = 2) +
  geom_line(aes(y = Spline, col = 'T2'), linewidth = 1.2) +
  scale_color_manual(name = "", labels = c("T1" = 'Calibrated period effect', 
                                           "T2" = 'Smoothing spline'),
                     values = c(col.dot.shade, col.dot)) +
  guides(color = guide_legend(
    override.aes = list(shape = c(16,NA), linetype = c('blank','solid'), 
                        size = c(3,3)))) + 
  theme(strip.background = element_blank(),
        legend.position = 'bottom') +
  geom_text(data = df %>% filter(Type == 'L[plain(t)]^plain((1))') %>% slice(1), 
            aes(x = 1880, y = 3.2, label = paste0("hat(lambda)[opt]^(1)==",l1)), 
            parse = TRUE, size = 4) + 
  geom_text(data = df %>% filter(Type == 'L[plain(t)]^plain((2))') %>% slice(1), 
            aes(x = 1995, y = 11.3, label = paste0("hat(lambda)[opt]^(2)==",l2)), 
            parse = TRUE, size = 4)

final_spline

### 2.3) Detect outliers ----

# Compute remainder components
df.Lt1 <- df.Lt1 %>% mutate('Remainder' = df.Lt1$PeriodEffect - df.Lt1$Spline)
df.Lt2 <- df.Lt2 %>% mutate('Remainder' = df.Lt2$PeriodEffect - df.Lt2$Spline)

# Outlier detection for years < and >= 1970
yv1     <- yv[1]:1969
ind1    <- yv1 - yv[1] + 1
yv2     <- 1970:2021
ind2    <- yv2 - yv[1] + 1

# Robust measure of location & scale
vecR    <- cbind(df.Lt1$Remainder, df.Lt2$Remainder)
rob.ls1 <- rrcov::CovMcd(vecR[ind1,], alpha = 0.75) # Robust location & scale
rob.ls2 <- rrcov::CovMcd(vecR[ind2,], alpha = 0.95) # Higher alpha-less outliers

# Calculate Mahalanobis distance
mhd1 <- sapply(ind1, function(x) sqrt(t(vecR[x,] - rob.ls1@raw.center) %*% 
                                        solve(rob.ls1@raw.cov) %*% 
                                        (vecR[x,] - rob.ls1@raw.center)))
mhd2 <- sapply(ind2, function(x) sqrt(t(vecR[x,] - rob.ls2@raw.center) %*% 
                                        solve(rob.ls2@raw.cov) %*% 
                                        (vecR[x,] - rob.ls2@raw.center)))
mhd  <- c(mhd1,mhd2)

# Threshold to identify outliers
thres   <- sqrt(qchisq(0.99, df = 2))

# Plot outliers
df <- data.frame('Year' = yv, R1 = df.Lt1$Remainder, R2 = df.Lt2$Remainder, 
                  Outlier = mhd > thres, 'Indicator' = as.character(yv %in% yv1))
df <- df %>% mutate('Label' = sapply(1:length(Outlier), function(x) 
  ifelse(Outlier[x] == TRUE, yv[x], "")))

df$Indicator <- factor(df$Indicator, levels = c(TRUE, FALSE), 
                       labels = c('Before 1970', 'After 1970'))

ellipse1 <- data.frame(ellipse(center = rob.ls1@raw.center, 
                               shape = rob.ls1@raw.cov, radius = thres,
                               segments = length(mhd1) - 1, draw = FALSE))
ellipse2 <- data.frame(ellipse(center = rob.ls2@raw.center, 
                               shape = rob.ls2@raw.cov, radius = thres,
                               segments = length(mhd2) - 1, draw = FALSE))
colnames(ellipse1) = colnames(ellipse2) <- c('R1','R2')

pMHD <- ggplot(df, aes(x = R1, y = R2, col = Outlier)) + 
  theme_bw(base_size = 20) + xlab(bquote('Remainder '~R['t']^"(1)")) + 
  ylab(bquote('Remainder '~R['t']^"(2)")) + 
  facet_wrap(~Indicator, scales = 'free') +
  geom_point(aes(col = as.character(Outlier), size = as.character(Outlier))) +
  geom_text_repel(aes(label = Label), max.overlaps = 100) + 
  geom_polygon(data = df %>% filter(Indicator == "Before 1970"), 
               aes(x = ellipse1$R1, y = ellipse1$R2), alpha = 0.1, 
               fill = col.dot, col = col.dot) + 
  geom_polygon(data = df %>% filter(Indicator == "After 1970"), 
               aes(x = ellipse2$R1, y = ellipse2$R2), alpha = 0.1, 
               fill = col.dot, col = col.dot) +
  scale_color_manual(name = 'Outlier', 
                     values = c('TRUE' = col.outl, 'FALSE' = col.dot)) + 
  scale_size_manual(name = 'Outlier', 
                    values = c('TRUE' = 4, 'FALSE' = 3)) + 
  theme(strip.background = element_blank(),
        legend.position = 'bottom')
pMHD

# Identified Outliers
outl1 <- df %>% filter(Outlier == 1) %>% pull(Year)
outl1

### 2.4) Recalibration ----

# Remove outliers from mortality data set
for(t in c('UNI','ALL')){
  if (t == "UNI"){
    for (c in Country){
      yv.c <- rownames(data_M[[t]][[c]][['dtx']])
      for(u in c('dtx','etx')){
        data_M[[t]][[c]][[u]][yv.c %in% outl1,] <- NA
      }
      data_M[[t]][[c]][['wa']][yv.c %in% outl1,] <- 0
    }
  }
  if (t == 'ALL'){
    for (u in c('dtx','etx')){
      data_M[[t]][[u]][yv %in% outl1,] <- NA
    }
    data_M[[t]][['wa']][yv %in% outl1,] <- 0
  }
}

# Refit mortality model on outlier-free calibration period
LL_M <- FitLiLeeNR(xv = xv, yv = yv, Data = data_M, Ax = TRUE, 
                   exclude.coh = FALSE, CountrySPEC = "NED", yvSPEC = yvSPEC, 
                   fit1 = fit701C, fit2 = fit701C0)

### 2.5) Mortality improvement specification ----

# Translate parameters to mortality improvement specificiation
KtM  <- diff(LL_M$NR$K.t)
LtM  <- diff(LL_M$NR$L.t)
ktM  <- diff(LL_M$NR$k.t)
ltM  <- diff(LL_M$NR$l.t)

# Plot common trend parameters KtM, LtM
df0x <- data.frame('Age' = xv, 'Ax' = LL_M$NR$A.x, 'Bx' = LL_M$NR$B.x, 
                   'Cx' = LL_M$NR$C.x) %>% 
  tidyr::gather('Type', 'Age effect', - Age)
df0t <- data.frame('Year' = yv[2]:2021, 'Kt' = KtM, 'Lt' = LtM) %>% 
  tidyr::gather('Type', 'Period effect', - Year)

df0x$Type <- factor(df0x$Type, levels = c('Ax', 'Bx', 'Cx'), labels = 
                      c('A[plain(x)]','B[plain(x)]^plain((1))',
                        'B[plain(x)]^plain((2))'))
df0t$Type <- factor(df0t$Type, levels = c('Kt', 'Lt'), labels = 
                      c('K[plain(t)]^plain((1))', 'K[plain(t)]^plain((2))'))

ctm0 <- ggplot(df0x, aes(x = Age, y = `Age effect`)) + 
  theme_bw(base_size = 15) +
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

ctm1 <- ggplot(df0t, aes(x = Year, y = `Period effect`)) + 
  theme_bw(base_size = 15) + 
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

ctm01 <- ggpubr::ggarrange(ctm0, ctm1, nrow = 2)
ctm01

# Plot country-specific deviation parameters
df0x <- data.frame('Age' = xv, 'bx' = LL_M$NR$b.x, 'cx' = LL_M$NR$c.x) %>% 
  tidyr::gather('Type', 'Age effect', - Age)
df0t <- data.frame('Year' = yv[2]:2021, 'kt' = ktM, 'lt' = ltM) %>% 
  tidyr::gather('Type', 'Period effect', - Year)

df0x$Type <- factor(df0x$Type, levels = c('bx', 'cx'), labels = 
                      c('beta[plain(x)]^plain((1))', 
                        'beta[plain(x)]^plain((2))'))
df0t$Type <- factor(df0t$Type, levels = c('kt', 'lt'), labels = 
                      c('kappa[plain(t)]^plain((1))', 
                        'kappa[plain(t)]^plain((2))'))

csdm0 <- ggplot(df0x, aes(x = Age, y = `Age effect`)) + 
  theme_bw(base_size = 15) +
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

csdm1 <- ggplot(df0t, aes(x = Year, y = `Period effect`)) + 
  theme_bw(base_size = 15) + 
  facet_wrap(~Type, scales = 'free', labeller = label_parsed) + 
  geom_point(col = col.dot, size = 1.5) + 
  theme(strip.background = element_blank())

csdm01 <- ggpubr::ggarrange(csdm0,csdm1, nrow = 2)
csdm01

##### 3) Define LL_M object in improvement setting ----
LL_M$NR$K.t <- KtM
LL_M$NR$L.t <- LtM
LL_M$NR$k.t <- ktM
LL_M$NR$l.t <- ltM

# saveRDS(LL_M, file = 'Save/LL_M_baseline.rds')