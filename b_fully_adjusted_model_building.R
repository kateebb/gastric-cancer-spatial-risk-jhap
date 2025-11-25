# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

# B. Building fully adjusted models: Multiple criteria to select fixed effects ------
  ## [1]: Build multiple regression model with all combinations of candidate fixed effects 
    # A. Selection criteria: WAIC, CPO, Mlik/Bayes factors
    # B. Check that all coefficient estimates are significant in candidate models
  ## [2]: Explore Pearson residuals of candidate models 
    # A. Compare distributions of residuals against null model (y ~ 1)
    # B. Calculate and compare Moran's I for Pearson residuals
  # [3]: Check for collinearity among fixed effects 
    # A. Pairwise correlations (Pearson, t.tests)
    # B. VIF values
    # C. Pairwise Bivariate Moran's I
    # D. Pairwise Local Pearson's Correlation Coefficients

# 0 - Load packages and data, prep for analysis (scaling, transformations) --------

setwd("")

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)
library(patchwork)
library(lctools)



x <- read_rds("data/mun_ds.rds") %>%
  mutate(len_bin = ifelse(round(len_pct)>=5,1,0),
         may_bin = ifelse(round(may_pct)>=5,1,0),
         wood_bin = ifelse(round(wood_pct)<90,0,1))

# Convert to regular dataframe for transformations and modeling
x2 <- x %>% st_set_geometry(NULL)

# List of covariates to be mean-centered and scaled
s <- c("urban_pct","wood_pct",
       "water_piped","water_well","water_river",
       "poverty_index","dist_km",
       "may_pct","len_pct",
       "avg_elev","avg_tri")

for (i in 1:length(s)) {
  x2[s[i]] = scale(unlist(x2[s[i]]))
}


# 1. Fit all possible models with series of candidate variables --------

vars <- c("dist_km","poverty_index","may_bin","len_bin","urban_pct","wood_pct","wood_bin")  
all_sets <- unlist(lapply(0:length(vars), function(k)
  if (k == 0) list(character(0)) else combn(vars, k, simplify = FALSE)), recursive = FALSE)

fit_one <- function(s) {
  rhs <- if (length(s) == 0) "1" else paste(s, collapse = " + ")
  fml <- as.formula(paste0("n1 ~ ", rhs))
  fit <- inla(fml, family = "poisson", data = x2, E = e1,
              control.predictor = list(compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE,
                                     dic = TRUE, mlik = TRUE, waic =TRUE,
                                     cpo=TRUE, config = TRUE),
              control.fixed = list(mean = 0, prec = 0.001))
  data.frame(
    predictors = if (length(s)) paste(s, collapse = "+") else "(intercept only)",
    k = length(s),
    WAIC = fit$waic$waic,
    CPO = -sum(log(fit$cpo$cpo), na.rm = TRUE),
    DIC  = fit$dic$dic,
    Mlik = fit$mlik,
    stringsAsFactors = FALSE
  )
}

tab <- do.call(rbind, lapply(all_sets, fit_one))
tab$row <- rownames(tab)
tab <- tab %>% filter(str_detect(row,"integration"))

## A. Fit criteria for model selection -------
# Parsimony filter: WAIC and CPO
best_WAIC  <- min(tab$WAIC, na.rm = TRUE)
best_CPO  <- min(tab$CPO, na.rm = TRUE)
best_Mlik <- max(tab$Mlik,na.rm=T)
cand <- subset(tab, (WAIC <= best_WAIC + 2) & (CPO <= best_CPO + 5) & (Mlik >= best_Mlik - 5))
cand <- cand[order(cand$k, cand$WAIC), ]   # smallest k first, to select most parsimonious model
cand

## B. Check model coefficients for model candidates: -------

# Model 1: y ~ dist_km + poverty_index + len_bin + may_bin
m1 <- inla(formula = n1 ~ dist_km+poverty_index+may_bin + len_bin ,
           family = "poisson", data = x2, E = e1,
           control.predictor = list(compute = TRUE),
           control.compute = list(return.marginals.predictor = TRUE,
                                  dic = TRUE, mlik = TRUE, waic =TRUE,
                                  cpo=TRUE, config = TRUE),
           control.fixed = list(mean = 0, prec = 0.001))

# Model 2: y ~ dist_km 
m2 <- inla(formula = n1 ~ dist_km,
           family = "poisson", data = x2, E = e1,
           control.predictor = list(compute = TRUE),
           control.compute = list(return.marginals.predictor = TRUE,
                                  dic = TRUE, mlik = TRUE, waic =TRUE,
                                  cpo=TRUE, config = TRUE),
           control.fixed = list(mean = 0, prec = 0.001))

## Values and significance are equivalent between models 
  ## % urban is not significant --> select most parsimonious model (Model 1, consistent with fit criteria)
coef <- data.frame(term = c(rownames(m1$summary.fixed),rownames(m2$summary.fixed)),
                   mod = c(rep("Mod 1",5),rep("Mod 2",2)),
                   est = c(exp(m1$summary.fixed[,1]),exp(m2$summary.fixed[,1])),
                   lo = c(exp(m1$summary.fixed[,3]),exp(m2$summary.fixed[,3])),
                   hi = c(exp(m1$summary.fixed[,5]),exp(m2$summary.fixed[,5]))) %>%
  mutate(sig = ifelse(between(1,lo,hi),0,1))


# 2. Assess Pearson's residuals  ----------

# Fit null model 
m0 <- inla(formula = n1 ~ 1,
                 family = "poisson", data = x2, E = e1,
                 control.predictor = list(compute = TRUE),
                 control.compute = list(return.marginals.predictor = TRUE,
                                        dic = TRUE, mlik = TRUE, waic =TRUE,
                                        cpo=TRUE, config = TRUE),
                 control.fixed = list(mean = 0, prec = 0.001))

## A. Distribution of Pearson residuals  -------
# Observed counts
y <- x2$n1
# Expected counts (offset)
E <- x2$e1
# exp(η̂)= Predicted posterior RR (compare to SIR y/E)
rr_hat0 <- m0$summary.fitted.values$mean  
rr_hat1 <- m1$summary.fitted.values$mean  
rr_hat2 <- m2$summary.fitted.values$mean  


plot(y/E,rr_hat0, xlab="SIR",ylab="RR (estimated from NULL model)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
plot(y/E,rr_hat1,  xlab="SIR",ylab="RR (estimated from Model 1)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
plot(y/E,rr_hat2,  xlab="SIR",ylab="RR (estimated from Model 2)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)


## Obtain predicted Y by multiplying by E
mu_hat0 <- E * rr_hat0
mu_hat1 <- E * rr_hat1
mu_hat2 <- E * rr_hat2


plot(mu_hat0,y, ylab="Observed N",xlab="Predicted N (from Null model)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
plot(mu_hat1,y, ylab="Observed N",xlab="Predicted N (from Model 1)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
plot(mu_hat2,y, ylab="Observed N",xlab="Predicted N (from Model 2)")
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)


res_raw0 <- y - mu_hat0
res_raw1 <- y - mu_hat1
res_raw2 <- y - mu_hat2

res_p0 <- (y - mu_hat0) / sqrt(mu_hat0)
res_p1 <- (y - mu_hat1) / sqrt(mu_hat1)
res_p2 <- (y - mu_hat2) / sqrt(mu_hat2)


# Raw residuals -- Null model
plot(mu_hat0, res_raw0,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Raw residual",
     main = "Null Model: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
# Pearson residuals
plot(mu_hat0, res_p0,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Pearson Residual",
     main = "Null Model: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
## % of obsv with pearson residuals < 2? --> 62%
sum(abs(res_p0)<2)/78*100

# Raw residuals -- Cand. Model 1
plot(mu_hat1, res_raw1,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Raw residual",
     main = "Candidate Model 1: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
# Pearson residuals
plot(mu_hat1, res_p1,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Pearson Residual",
     main = "Candidate Model 1: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
## % of obsv with pearson residuals < 2? --> 87%
sum(abs(res_p1)<2)/78*100

# Raw residuals -- Candidate Model 2
plot(mu_hat2, res_raw2,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Raw residual",
     main = "Candidate Model 2: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
# Pearson residuals
plot(mu_hat2, res_p2,
     xlab = "Predicted Counts (μ̂i = Ei × RR̂i)",
     ylab = "Pearson Residual",
     main = "Candidate Model 2: Residuals vs. Predicted Counts")
abline(h = 0, col = "gray", lty = 2)
## % of obsv with pearson residuals < 2? --> 78%
sum(abs(res_p2)<2)/78*100

## Estimated values correlate strongly across Mod 1 and Mod 2
cor(rr_hat1,rr_hat2) ## .85
cor(mu_hat1,mu_hat2) ## .97

## Plot residuals on maps 
x$rr0 <- rr_hat0
x$rr1 <- rr_hat1
x$rr2 <- rr_hat2

x$pres0 <- res_p0
x$pres1 <- res_p1
x$pres2 <- res_p2


p1 <- ggplot(x %>%
               mutate(rr=rr0)) +
  geom_sf(aes(fill = rr)) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0,0.5,1,2,3,4,5), 
                        limits = c(0,5)) +
  ggtitle("RR estimate: Null Model") +
  theme_map()

p2 <- ggplot(x %>%
               mutate(rr=rr1)) +
  geom_sf(aes(fill = rr)) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0,0.5,1,2,3,4,5), 
                        limits = c(0,5)) +
  ggtitle("RR estimate: Candidate Model 1") + 
  theme_map()

p3 <- ggplot(x %>%
               mutate(rr=rr2)) +
  geom_sf(aes(fill = rr)) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0,0.5,1,2,3,4,5), 
                        limits = c(0,5)) +
  ggtitle("RR estimate: Candidate Model 2") + 
  theme_map()


(p1 + p2 + p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Residuals
p1 <- ggplot(data = x %>%
               mutate(pres=pres0)) +
  geom_sf(aes(fill = pres)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits=c(-3.5,6)) +
  labs(title = "Mapped Pearson Residuals",
       fill = "Residuals",
       subtitle = "Null model") +
  theme_map()

p2 <- ggplot(data = x %>%
               mutate(pres=pres1)) +
  geom_sf(aes(fill = pres)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits=c(-3.5,6)) +
  labs(title = "Mapped Pearson Residuals",
       fill = "Residuals",
       subtitle = "Candidate Model 1") +
  theme_map()

p3 <- ggplot(data = x %>%
               mutate(pres=pres2)) +
  geom_sf(aes(fill = pres)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits=c(-3.5,6)) +
  labs(title = "Mapped Pearson Residuals",
       fill = "Residuals",
       subtitle = "Candidate Model 2") +
  theme_map()


(p1 + p2 + p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")


## B. Calculate Moran's I values for candidate model Pearson residuals ------

## Equivalent Moran's I values for candidate models.
  ## Greatly reduced residual spatial structure, but some remains.

nb <- poly2nb(x)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Null Model: 0.58***
mi  <- moran.test(x$pres0, lw, randomisation = TRUE, alternative = "greater")
mi$estimate["Moran I statistic"]; round(mi$p.value,3)

# Candidate Model 1: 0.18***
mi  <- moran.test(x$pres1, lw, randomisation = TRUE, alternative = "greater")
mi$estimate["Moran I statistic"]; round(mi$p.value,3)

# Candidate Model 2: 0.24***
mi  <- moran.test(x$pres2, lw, randomisation = TRUE, alternative = "greater")
mi$estimate["Moran I statistic"]; round(mi$p.value,3)

# 3. Collinearity checks among fixed effects in selected fully adjusted model -----

# Fully adjusted model: y ~ dist_km + may_bin + len_bin + poverty_index

## A. Pairwise correlation values -----
## Bimodal distribution of len_bin and may_bin means we use spearman corr test
z <- x2 %>% mutate(sir = n1/e1) %>%
  select(n1,e1,sir,dist_km,len_bin,may_bin,poverty_index,may_bin,len_bin)
vars <- c("dist_km","poverty_index","may_bin","len_bin")  

v1 <- c()
v2 <- c()
r <- c()
pval <- c()

for (i in c(1:length(vars))){
  
  q1 <- vars[i]
  
  for(k in c(1:length(vars))) {
    
    if (i == k) next
    
    q2 <- vars[k]
    
    test <- cor.test(x = as.vector(x2[,q1]),
                     y = as.vector(x2[,q2]),
                     method = "pearson")
    
    r <- c(r,round(unname(test$estimate),2))
    pval <- c(pval,round(unname(test$p.value),2))
    v1 <- c(v1,q1)
    v2 <- c(v2,q2)
  }
}


# TABLE A1 ------
b <- data.frame(v1,v2,r,pval) %>%
  arrange(pval,desc(abs(r)))
write_csv(b,"outputs/supplementary/tabA1_PearsonCorr_FE.csv")

## B. VIF values from fully adjusted model -------
model <- glm(n1 ~ dist_km + may_bin + len_bin + poverty_index + offset(log(e1)),
             data = x2, family = "poisson")
performance::check_collinearity(model)


## C. Pairwise Bivariate Moran's I --------
vars <- c("dist_km","poverty_index","may_bin","len_bin")  

v1 <- c()
v2 <- c()
bivarMI <- c()
pval <- c()

for (i in c(1:length(vars))){
  
  q1 <- vars[i]
  
  for(k in c(1:length(vars))) {
    
    if (i == k) next
    
    q2 <- vars[k]

    test <- bispdep::moranbi.mc(varX = as.vector(x2[,q1]), 
                                   varY = as.vector(x2[,q2]),
                                   listw = lw, nsim = 999,
                                  zero.policy = TRUE,
                                  alternative = "greater")
    
    bivarMI <- c(bivarMI,round(unname(test$statistic),2))
    pval <- c(pval,round(unname(test$p.value),2))
    v1 <- c(v1,q1)
    v2 <- c(v2,q2)
  }
}

# TABLE A2 ----
b <- data.frame(v1,v2,bivarMI,pval) %>%
  arrange(pval,desc(abs(bivarMI)))
write_csv(b,"outputs/supplementary/tabA2_Bivariate_MoranI.csv")

## D. Pairwise Local Pearson's R --------
library(lctools)
x <- read_rds("data/mun_ds.rds") %>%
  mutate(len_bin = ifelse(round(len_pct)>=5,1,0),
         may_bin = ifelse(round(may_pct)>=5,1,0),
         wood_bin = ifelse(round(wood_pct)<90,0,1))

load("final_models/mod_n1_v13.RData")


# Convert to regular dataframe for transformations and modeling
x2 <- x %>% st_set_geometry(NULL)

# List of covariates to be mean-centered and scaled
s <- c("urban_pct","wood_pct",
       "water_piped","water_well","water_river",
       "poverty_index","dist_km",
       "may_pct","len_pct",
       "avg_elev","avg_tri")

for (i in 1:length(s)) {
  x2[s[i]] = scale(unlist(x2[s[i]]))
}

z <- x2
z$spatre <- model$summary.random$spat_re$mean[1:78]

library(tidyverse)
library(ggplot2)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)
library(ggthemes)
library(patchwork)
library(svglite)

x <- read_rds("data/mun_ds.rds")
dep <- st_read("data/hnd_dep_clean.shp") 
clin <- st_read("data/clinics.shp") %>%
  filter(n=="WRH")


## Fixed effects -----
### Dist ~ Poverty Index ------
lc1 <- lcorrel(z[,c("dist_km","poverty_index")],
                     bw=0.65,
                     st_coordinates(st_centroid(x)))
## Global = 0.37
cor(z$dist_km,z$poverty_index)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_dist_povind <- lc1$LPCC[,2]
x$lpcc.p_dist_povind <- lc1$LPCC_sig_BF[,2]

### Dist ~ Maya/Chorti ------
lc1 <- lcorrel(z[,c("dist_km","may_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.09
cor(z$dist_km,z$may_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_dist_may <- lc1$LPCC[,2]
x$lpcc.p_dist_may <- lc1$LPCC_sig_BF[,2]

### Dist ~ Lenca ------
lc1 <- lcorrel(z[,c("dist_km","len_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = 0.48
cor(z$dist_km,z$len_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_dist_len <- lc1$LPCC[,2]
x$lpcc.p_dist_len <- lc1$LPCC_sig_BF[,2]


### Poverty Index ~ Maya/Chorti ------
lc1 <- lcorrel(z[,c("poverty_index","may_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.06
cor(z$poverty_index,z$may_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_pov_may <- lc1$LPCC[,2]
x$lpcc.p_pov_may <- lc1$LPCC_sig_BF[,2]

### Poverty Index ~ Lenca ------
lc1 <- lcorrel(z[,c("poverty_index","len_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = 0.26
cor(z$poverty_index,z$len_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_pov_len <- lc1$LPCC[,2]
x$lpcc.p_pov_len <- lc1$LPCC_sig_BF[,2]


### Maya/Chorti ~ Lenca ------
lc1 <- lcorrel(z[,c("may_bin","len_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.30
cor(z$may_bin,z$len_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

x$lpcc_may_len <- lc1$LPCC[,2]
x$lpcc.p_may_len <- lc1$LPCC_sig_BF[,2]



## Random effect / FEs -----
w <- x %>% mutate(ind = c(1:78)) %>% select(ind)
### SpatRE ~ Dist ------
lc1 <- lcorrel(z[,c("spatre","dist_km")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.23
cor(z$spatre,z$dist_km)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

w$lpcc_re_dist <- lc1$LPCC[,2]
w$lpcc.p_re_dist <- lc1$LPCC_sig_BF[,2]

### SpatRE ~ Poverty index ------
lc1 <- lcorrel(z[,c("spatre","poverty_index")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.02
cor(z$spatre,z$poverty_index)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

w$lpcc_re_pov <- lc1$LPCC[,2]
w$lpcc.p_re_pov <- lc1$LPCC_sig_BF[,2]

### SpatRE ~ Maya/Chorti ------
lc1 <- lcorrel(z[,c("spatre","may_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = 0.14
cor(z$spatre,z$may_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

w$lpcc_re_may <- lc1$LPCC[,2]
w$lpcc.p_re_may <- lc1$LPCC_sig_BF[,2]


### SpatRE ~ Lenca ------
lc1 <- lcorrel(z[,c("spatre","len_bin")],
               bw=0.65,
               st_coordinates(st_centroid(x)))
## Global = -0.40
cor(z$spatre,z$len_bin)
mean(lc1$LPCC[,2])
hist(lc1$LPCC[,2])
hist(lc1$LPCC_sig_BF[,2])

w$lpcc_re_len <- lc1$LPCC[,2]
w$lpcc.p_re_len <- lc1$LPCC_sig_BF[,2]

## Figures A1 & A2: FE~FE Local PCC maps  -----

m <- x %>% 
  mutate(ind = c(1:78)) %>%
  select(ind,starts_with("lpcc_")) %>% 
  pivot_longer(cols = c(lpcc_dist_povind,
                        lpcc_dist_may,
                        lpcc_dist_len,
                        lpcc_pov_may,
                        lpcc_pov_len,
                        lpcc_may_len),
               names_to = "var",
               values_to = "lpcc") %>%
  mutate(var = str_remove(var,"lpcc_"),
         lab = case_when(var == "dist_povind" ~ "Dist. to WRH ~ Pov. Index",
                         var == "dist_may" ~ "Dist. to WRH ~ Maya/Chortí",
                         var == "dist_len" ~ "Dist. to WRH ~ Lenca",
                         var == "pov_may" ~ "Pov. Index ~ Maya/Chortí",
                         var == "pov_len" ~ "Pov. Index ~ Lenca",
                         var == "may_len" ~ "Maya/Chortí ~ Lenca"))

m2 <- x %>% 
  mutate(ind = c(1:78)) %>%
  select(ind,starts_with("lpcc.p_")) %>% 
  pivot_longer(cols = c(lpcc.p_dist_povind,
                        lpcc.p_dist_may,
                        lpcc.p_dist_len,
                        lpcc.p_pov_may,
                        lpcc.p_pov_len,
                        lpcc.p_may_len),
               names_to = "var",
               values_to = "pval") %>%
  mutate(var = str_remove(var,"lpcc.p_"),
         lab = case_when(var == "dist_povind" ~ "Dist. to WRH ~ Pov. Index",
                         var == "dist_may" ~ "Dist. to WRH ~ Maya/Chortí",
                         var == "dist_len" ~ "Dist. to WRH ~ Lenca",
                         var == "pov_may" ~ "Pov. Index ~ Maya/Chortí",
                         var == "pov_len" ~ "Pov. Index ~ Lenca",
                         var == "may_len" ~ "Maya/Chortí ~ Lenca"))

m3 <- left_join(m,m2 %>% st_set_geometry(NULL)) %>%
  mutate(sigvals = ifelse(pval < 0.05,lpcc,NA))


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = m3, aes(fill = lpcc), col = NA) + 
  facet_wrap(vars(lab),nrow=2) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 0,      
    name = "LPCC",
    breaks = c(-0.5,-0.25,0,0.25,0.5,0.66),
    labels = c("-0.5","-0.25","0.0","0.25","0.50","0.65"),
    limits = c(-0.5,0.66)
  ) +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(legend.position = 'left')

ggsave(a,
       device = "png",
       filename = paste0("outputs/supplementary/figA1_lpccFE.png"), # change title here
       width = 7.2,
       height = 6.1,
       units = "in")


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = m3, aes(fill = sigvals), col = NA) + 
  facet_wrap(vars(lab),nrow=2) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 0,      
    name = "LPCC",
    breaks = c(-0.5,-0.25,0,0.25,0.5,0.66),
    labels = c("-0.5","-0.25","0.0","0.25","0.50","0.65"),
    limits = c(-0.5,0.66)
  ) +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(legend.position = 'left')

ggsave(a,
       device = "png",
       filename = paste0("outputs/supplementary/figA2_lpccFE_Sig.png"), # change title here
       width = 7.2,
       height = 6.1,
       units = "in")


## Figures A7 & A8: FE~RE Local PCC maps ------
q <- w %>% 
  mutate(ind = c(1:78)) %>%
  select(ind,starts_with("lpcc_")) %>% 
  pivot_longer(cols = c(lpcc_re_pov,
                        lpcc_re_may,
                        lpcc_re_len,
                        lpcc_re_dist),
               names_to = "var",
               values_to = "lpcc") %>%
  mutate(var = str_remove(var,"lpcc_"),
         lab = case_when(var == "re_pov" ~ "Pov. Index ~ Spatial RE",
                         var == "re_may" ~ "Maya/Chortí ~ Spatial RE",
                         var == "re_len" ~ "Lenca ~ Spatial RE",
                         var == "re_dist" ~ "Dist. to WRH ~ Spatial RE"))

q2 <- w %>% 
  mutate(ind = c(1:78)) %>%
  select(ind,starts_with("lpcc.p_")) %>% 
  pivot_longer(cols = c(lpcc.p_re_pov,
                        lpcc.p_re_may,
                        lpcc.p_re_len,
                        lpcc.p_re_dist),
               names_to = "var",
               values_to = "pval") %>%
  mutate(var = str_remove(var,"lpcc.p_"),
         lab = case_when(var == "re_pov" ~ "Pov. Index ~ Spatial RE",
                         var == "re_may" ~ "Maya/Chortí ~ Spatial RE",
                         var == "re_len" ~ "Lenca ~ Spatial RE",
                         var == "re_dist" ~ "Dist. to WRH ~ Spatial RE"))


z3 <- left_join(q,q2 %>% st_set_geometry(NULL)) %>%
  mutate(sigvals = ifelse(pval < 0.05,lpcc,NA))


a<-ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = z3, aes(fill = lpcc), col = NA) + 
  facet_wrap(vars(lab),nrow=1) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 0,      
    name = "LPCC",
    breaks = c(-0.51,-0.25,0,0.27),
    labels = c("-0.5","-0.25","0.0","0.25"),
    limits = c(-0.51,0.27)
  ) +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(legend.position = 'left')

ggsave(a,
       device = "png",
       filename = paste0("outputs/supplementary/figA7_lpccRE.png"), 
       width = 9,
       height = 4,
       units = "in")

a<-ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = z3, aes(fill = sigvals), col = NA) + 
  facet_wrap(vars(lab),nrow=1) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 0,      
    name = "LPCC",
    breaks = c(-0.51,-0.25,0,0.27),
    labels = c("-0.5","-0.25","0.0","0.25"),
    limits = c(-0.51,0.27)
  ) +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(legend.position = 'left')

ggsave(a,
       device = "png",
       filename = paste0("outputs/supplementary/figA8_lpccRE.png"), 
       width = 9,
       height = 4,
       units = "in")

