# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

## D. Priors sensitivity checks 

## [Goal]: Assess that the priors do not strongly influence posterior estimates of coefficients, parameters, or RRs
  ## 1 - Coefficients & Intercept priors
  ## 2 - Precision for spatial random effect
  ## 3 - Mixing parameter for spatial random effect
  ## 4 - Overdispersion parameter
## 4 - Examine spatial confounding (correlation between spatial effect and FEs)
## Assessments:
  ## A. Fit criteria
  ## B. Posterior distributions of coefficients 
  ## C. Posterior distributions of hyperparameters
  ## D. Correlations among posterior mean RR estimates

# 0 - Load packages and data --------

setwd("")

library(tidyverse)
library(ggplot2)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)
library(patchwork)

x <- read_rds("data/mun_ds.rds") %>%
  mutate(len_bin = ifelse(round(len_pct)>=5,1,0),
         may_bin = ifelse(round(may_pct)>=5,1,0))

# Convert to regular dataframe for transformations and modeling
x2 <- x %>% st_set_geometry(NULL)

# List of covariates to be mean-centered and scaled
s <- c("poverty_index","dist_km")

for (i in 1:length(s)) {
  x2[s[i]] = scale(unlist(x2[s[i]]))
}

load("final_models/mod_n1_v06.RData")
base <- model

f <- n1 ~ 1 + dist_km + poverty_index + len_bin + may_bin 
f2 <- n1 ~ 1 + dist_km + poverty_index + len_log + may_log 

rm(list=setdiff(ls(), c("x2","x","base","f","f2","g")))

# 1 - Coefficients priors -------
## BASE: Flat prior for fixed effect coefficients: N(0,31.6) --> (mean = 0, prec = 0.001)

## Weakly informative prior: N(0, 0.5) => prec = 1 / (0.5^2) --> (mean = 0, prec = 4)
s1.1 <- inla(formula = f, family = "poisson", data = x2, E = e1,
     control.predictor = list(compute = TRUE),
     control.compute = list(return.marginals.predictor = TRUE,
                            dic = TRUE, mlik = TRUE, waic =TRUE,
                            cpo=TRUE, config = TRUE),
     control.fixed = list(mean = 0, prec = 4))

# Strongly informative prior: N(0, 0.1) => prec = 1 / (0.1^2) --> (mean = 0, prec = 100)
s1.2 <- inla(formula = f, family = "poisson", data = x2, E = e1,
             control.predictor = list(compute = TRUE),
             control.compute = list(return.marginals.predictor = TRUE,
                                    dic = TRUE, mlik = TRUE, waic =TRUE,
                                    cpo=TRUE, config = TRUE),
             control.fixed = list(mean = 0, prec = 50))



## A. Assess fit criteria ------
formulas <- list(base,s1.1,s1.2)
# create model label string
lab <- c("flat (base)", "weakly informative", "strongly informative")

# create table to store DIC, WAIC, and CPO metrics
table0 <- data.table(Model  = c("flat (base)", "weakly informative", "strongly informative"),
                     DIC = NA,
                     WAIC = NA,
                     Mlik = NA,
                     CPO = NA,
                     CPO_fail = NA)

for(i in 1:length(formulas)){
  model <- formulas[[i]]
  table0$DIC[i] <- round(model$dic$dic, 1)
  table0$WAIC[i] <- round(model$waic$waic,1)
  table0$Mlik[i] <- round(unname(model$mlik[1,1]),1)
  table0$CPO[i] <- -sum(log(model$cpo$cpo))
  table0$CPO_fail[i] <- sum(model$cpo$failure>0.5)
}

# view table
table0

## B. Posterior distributions of estimated coefficients --------

dist <- data.frame(x=NA,y=NA,mod=NA,co=NA)

for(i in c(1:3)){
  model <- formulas[[i]]
  d1 <- model$marginals.fixed$`(Intercept)` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "1. Intercept")
  d2 <- model$marginals.fixed$dist_km %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "2. Distance to WRH (km)")
  d3 <- model$marginals.fixed$poverty_index %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "3. Poverty index")
  d4 <- model$marginals.fixed$len_bin %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "4. Lenca")
  d5 <- model$marginals.fixed$may_bin %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "5. Maya/Chortí")
  
  dist <- bind_rows(dist,d1,d2,d3,d4,d5) %>% na.omit()
  
}

(a <- ggplot(dist %>% mutate(x = exp(x)), aes(x = x, y = y,group = mod,color=mod)) +
  facet_wrap(~co,ncol=1)+
  geom_line() +
  scale_color_manual(values=c("black","red","blue"),"Prior") + 
  geom_vline(xintercept=1)+
  labs(subtitle = "Binned population group indicators",
       x = "Estimate",
       y = "Density") + 
    theme(legend.position = "right"))

# FIGURE A3 -------
ggsave(a,filename="outputs/supplementary/figA3_coefpriors_sensitivity.png",
       device="png",
       width=4.7,height=5.4,units="in",dpi=300)

# 2 - Precision (theta1) hyperparameter priors -------
## BASE: Var should not exceed 1: {theta1 = pc.prior, param = 1, 0.01; theta2 = pc, param = 0.5, 0.5}

load("final_models/mod_n1_v13.RData")
base <- model

x2$spat_re <- c(1:78)
g <- inla.read.graph(filename = "data/munis.adj")

f <- n1 ~ 1 + dist_km + poverty_index + len_bin + may_bin +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(0.1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5))))

## Conservative precision: {theta1 = pc.prior, param = 0.5, 0.01; theta2 = pc, param = 0.5, 0.5}
f <- n1 ~ 1 + dist_km + poverty_index + len_log + may_log +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(0.5, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5))))
s2.1 <- inla(formula = f, family = "nbinomial", data = x2, E = e1,
             control.predictor = list(compute = TRUE),
             control.compute = list(return.marginals.predictor = TRUE,
                                    dic = TRUE, mlik = TRUE, waic =TRUE,
                                    cpo=TRUE, config = TRUE),
             control.fixed = list(mean = 0, prec = 0.001))

# Liberal theta1: {theta1 = pc.prior, param = 2, 0.01; theta2 = pc, param = 0.5, 0.5}
f <- n1 ~ 1 + dist_km + poverty_index + len_log + may_log +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(2, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5))))
s2.2 <- inla(formula = f, family = "nbinomial", data = x2, E = e1,
             control.predictor = list(compute = TRUE),
             control.compute = list(return.marginals.predictor = TRUE,
                                    dic = TRUE, mlik = TRUE, waic =TRUE,
                                    cpo=TRUE, config = TRUE),
             control.fixed = list(mean = 0, prec = 0.001))

## A. Assess fit criteria ------
formulas <- list(base,s2.1,s2.2)
# create model label string
lab <- c("base", "conservative precision", "liberal precision")

# create table to store DIC, WAIC, CPO metrics
table0 <- data.table(Model  = c("base", "conservative theta1 (small variance)", "liberal theta1 (large variance)"),
                     DIC = NA,
                     WAIC = NA,
                     Mlik = NA,
                     CPO = NA,
                     CPO_fail = NA)

for(i in 1:length(formulas)){
  model <- formulas[[i]]
  table0$DIC[i] <- round(model$dic$dic, 1)
  table0$WAIC[i] <- round(model$waic$waic,1)
  table0$Mlik[i] <- round(unname(model$mlik[1,1]),1)
  table0$CPO[i] <- -sum(log(model$cpo$cpo))
  table0$CPO_fail[i] <- sum(model$cpo$failure>0.5)
}

# view table
table0

## B. Posterior distributions of estimated coefficients --------

dist <- data.frame(x=NA,y=NA,mod=NA,co=NA)

for(i in c(1:3)){
  model <- formulas[[i]]
  d1 <- model$marginals.fixed$`(Intercept)` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "1. Intercept")
  d2 <- model$marginals.fixed$dist_km %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "2. Distance to WRH (km)")
  d3 <- model$marginals.fixed$poverty_index %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "3. Poverty index")
  d4 <- model$marginals.fixed$len_log %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "4. log(% Lenca)")
  d5 <- model$marginals.fixed$may_log %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "5. log(% Maya/Chorti)")
  
  dist <- bind_rows(dist,d1,d2,d3,d4,d5) %>% na.omit()
  
}

(a <- ggplot(dist %>% mutate(x = exp(x)), aes(x = x, y = y,group = mod,color=mod)) +
    facet_wrap(~co,ncol=1)+
    geom_line() +
    scale_color_manual(values=c("black","red","blue"),"Prior") + 
    geom_vline(xintercept=1)+
    labs(subtitle = "Coefficient posteriors",
         x = "Estimate",
         y = "Density"))


## C. Posterior distributions of RE hyperparameters -------
dist2 <- data.frame(x=NA,y=NA,mod=NA,co=NA)

for(i in c(1:3)){
  model <- formulas[[i]]
  d1 <- model$marginals.hyperpar$`Precision for spat_re` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "1. Precision (theta1)")
  d2 <- model$marginals.hyperpar$`Phi for spat_re` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "2. Phi (theta2)")
  
  
  dist2 <- bind_rows(dist2,d1,d2) %>% na.omit()
  
}

(b<-ggplot(dist2, aes(x = x, y = y,group = mod,color=mod)) +
    facet_wrap(~co,ncol=1,scales = "free")+
    geom_line() +
    scale_color_manual(values=c("black","red","blue"),"Prior") + 
    labs(subtitle = "Hyperparameter posteriors",
         x = "Estimate",
         y = "Density"))


# FIGURE A4 ------
z <- (a + b) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") &
  plot_annotation(title = "Sensitivity check on random effects precision priors")

ggsave(z,filename="outputs/supplementary/figA4_precpriors_sensitivity.png",
       device="png",
       width=5,height=6,units="in",dpi=300)



# 3 - Phi (theta2) hyperparameter priors -------
## BASE: Even spatial/aspatial structure: {theta1 = pc.prior, param = 1, 0.01; theta2 = pc, param = 0.5, 0.5}

## Aspatial: {theta1 = pc.prior, param = 1, 0.01; theta2 = pc, param = 0.2, 0.5}
f <- n1 ~ 1 + dist_km + poverty_index + len_log + may_log +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.2, 0.5))))
s3.1 <- inla(formula = f, family = "nbinomial", data = x2, E = e1,
             control.predictor = list(compute = TRUE),
             control.compute = list(return.marginals.predictor = TRUE,
                                    dic = TRUE, mlik = TRUE, waic =TRUE,
                                    cpo=TRUE, config = TRUE),
             control.fixed = list(mean = 0, prec = 0.001))

# Spatial: {theta1 = pc.prior, param = 1, 0.01; theta2 = pc, param = 0.8, 0.5}
f <- n1 ~ 1 + dist_km + poverty_index + len_log + may_log +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5))))
s3.2 <- inla(formula = f, family = "nbinomial", data = x2, E = e1,
             control.predictor = list(compute = TRUE),
             control.compute = list(return.marginals.predictor = TRUE,
                                    dic = TRUE, mlik = TRUE, waic =TRUE,
                                    cpo=TRUE, config = TRUE),
             control.fixed = list(mean = 0, prec = 0.001))

## A. Assess fit criteria ------
formulas <- list(base,s3.1,s3.2)
# create model label string
lab <- c("base", "phi favors iid", "phi favors spatial")

# create table to store DIC, and CPO metrics
table0 <- data.table(Model  = c("base", "phi favors iid", "phi favors spatial"),
                     DIC = NA,
                     WAIC = NA,
                     Mlik = NA,
                     CPO = NA,
                     CPO_fail = NA)

for(i in 1:length(formulas)){
  model <- formulas[[i]]
  table0$DIC[i] <- round(model$dic$dic, 1)
  table0$WAIC[i] <- round(model$waic$waic,1)
  table0$Mlik[i] <- round(unname(model$mlik[1,1]),1)
  table0$CPO[i] <- -sum(log(model$cpo$cpo))
  table0$CPO_fail[i] <- sum(model$cpo$failure>0.5)
}

# view table
table0

## B. Posterior distributions of estimated coefficients --------

dist <- data.frame(x=NA,y=NA,mod=NA,co=NA)

for(i in c(1:3)){
  model <- formulas[[i]]
  d1 <- model$marginals.fixed$`(Intercept)` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "1. Intercept")
  d2 <- model$marginals.fixed$dist_km %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "2. Distance to WRH (km)")
  d3 <- model$marginals.fixed$poverty_index %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "3. Poverty index")
  d4 <- model$marginals.fixed$len_log %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "4. log(% Lenca)")
  d5 <- model$marginals.fixed$may_log %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "5. log(% Maya/Chorti)")
  
  dist <- bind_rows(dist,d1,d2,d3,d4,d5) %>% na.omit()
  
}

(a <- ggplot(dist %>% mutate(x = exp(x)), aes(x = x, y = y,group = mod,color=mod)) +
    facet_wrap(~co,ncol=1)+
    geom_line() +
    scale_color_manual(values=c("black","red","blue"),"Prior") + 
    geom_vline(xintercept=1)+
    labs(subtitle = "Coefficient posteriors",
         x = "Estimate",
         y = "Density"))


## C. Posterior distributions of RE hyperparameters -------
dist2 <- data.frame(x=NA,y=NA,mod=NA,co=NA)

for(i in c(1:3)){
  model <- formulas[[i]]
  d1 <- model$marginals.hyperpar$`Precision for spat_re` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "1. Precision (theta1)")
  d2 <- model$marginals.hyperpar$`Phi for spat_re` %>%
    as.data.frame() %>%
    mutate(mod = lab[i],
           co = "2. Phi (theta2)")
  
  
  dist2 <- bind_rows(dist2,d1,d2) %>% na.omit()
  
}

(b<-ggplot(dist2, aes(x = x, y = y,group = mod,color=mod)) +
    facet_wrap(~co,ncol=1,scales = "free")+
    geom_line() +
    scale_color_manual(values=c("black","red","blue"),"Prior") + 
    labs(subtitle = "Hyperparameter posteriors",
         x = "Estimate",
         y = "Density"))

# FIGURE A5 -------
z <- (a + b) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") &
  plot_annotation(title = "Sensitivity check on mixing parameter Phi priors")

ggsave(z,filename="outputs/supplementary/figA5_phipriors_sensitivity.png",
       device="png",
       width=5,height=6,units="in",dpi=300)

