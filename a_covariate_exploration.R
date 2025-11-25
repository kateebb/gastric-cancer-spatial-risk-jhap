# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

## A. Exploring hypothesized covariates in reduced models 

## Outline:
  # 0. Set up workspace
  # 1. Fit reduced models, with and without distance to WRH
  # 2. Compare fit and coefficient stability across two model versions

# 0 - Load packages and data, prep for analysis (scaling, transformations) --------

setwd("")

library(tidyverse)
library(ggplot2)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)

x <- read_rds("data/mun_ds.rds") %>%
  mutate(len_bin = ifelse(round(len_pct)>=5,1,0),
         may_bin = ifelse(round(may_pct)>=5,1,0),
         wood_bin = ifelse(round(wood_pct)<90,0,1))
x$len_rint <- qnorm((rank(x$len_pct, ties.method="average") - 0.5) / length(x$len_pct))
x$may_rint <- qnorm((rank(x$may_pct, ties.method="average") - 0.5) / length(x$may_pct))

# Convert to regular dataframe for transformations and modeling
x2 <- x %>% st_set_geometry(NULL)

# List of covariates to be mean-centered and scaled
s <- c("urban_pct","wood_pct",
       "water_piped","water_well","water_river",
       "poverty_index","dist_km",
       "avg_elev","avg_tri")

for (i in 1:length(s)) {
  x2[s[i]] = scale(unlist(x2[s[i]]))
}

# Create model with default (non-informative) priors; based on all persons, all GC
mymodel <- function(formula, data = x2, family = "poisson"){
  model <- inla(formula = formula, family = "poisson", data = x2, E = e1,
                control.predictor = list(compute = TRUE),
                control.compute = list(return.marginals.predictor = TRUE,
                                       dic = TRUE, mlik = TRUE, waic =TRUE,
                                       cpo=TRUE, config = TRUE),
                control.fixed = list(mean = 0, prec = 0.001))
  return(model)
}

# 1 - Fit FE reduced models ----
## A. Without measure of distance to clinic ------
baseformula <- n1 ~ 1 

# define formulas by updating the baseline formula
f1.1 <- update.formula(baseformula, ~. + poverty_index)
f1.2 <- update.formula(baseformula, ~. + urban_pct)
f1.3 <- update.formula(baseformula, ~. + len_log)
f1.4 <- update.formula(baseformula, ~. + len_bin)
f1.5 <- update.formula(baseformula, ~. + may_log)
f1.6 <- update.formula(baseformula, ~. + may_bin)
f1.7 <- update.formula(baseformula, ~. + avg_tri)
f1.8 <- update.formula(baseformula, ~. + avg_elev)
f1.9 <- update.formula(baseformula, ~. + wood_pct)
f1.10 <- update.formula(baseformula, ~. + water_piped)
f1.11 <- update.formula(baseformula, ~. + water_well)
f1.12 <- update.formula(baseformula, ~. + water_river)
f1.13 <- update.formula(baseformula, ~. + may_rint)
f1.14 <- update.formula(baseformula, ~. + len_rint)
f1.15 <- update.formula(baseformula, ~. + wood_bin)


# create a list of formulas
formulas <- list(baseformula, f1.1, f1.2, f1.3, f1.4, f1.5,
                 f1.6, f1.7, f1.8, f1.9, f1.10, f1.11, f1.12,
                 f1.13, f1.14, f1.15)

# create model label string
lab <- paste0("a_v1_",c("basemodel", "poverty_index", "urb_pct",
                        "len_log", "len_bin","may_log", "may_bin", 
                        "avg_TRI", "avg_elev",
                      "pct_Woodcooking", 
                      "pct_Pipedwater", "pct_Unprotectedwell", "pct_Riverwater",
                      "may_rint","len_rint","wood_bin"))

models <- lapply(1:length(formulas), 
                 function(i) {
                   model <- mymodel(formulas[[i]], x2)
                   save(model, file = paste0("fitting_models/", lab[i],".RData"))})


# create table to store DIC, WAIC, and Marginal Likelihood --> use to select best dist to clinic fit
table0a <- data.table(Model  = c("Base model", "Poverty index", "% urban", "% Lenca: log", "% Lenca: bin",
                                 "% Maya/C: log", "% Maya/C: bin", "Avg. TRI", "Avg. Elevation",
                                 "Cooking: % wood", "Water: % piped", "Water: % well", "Water: % river",
                                 "% Maya/C: RINT","% Lenca: RINT","Cooking: wood bin"),
                      version = "v1 - without dist_km",
                      DIC = NA,
                      WAIC = NA,
                      CPO = NA,
                      Mlik = NA)

for(i in 1:length(formulas)){
  load(paste0("fitting_models/",lab[i],".RData"))
  table0a$DIC[i] <- round(model$dic$dic, 1)
  table0a$WAIC[i] <- round(model$waic$waic,1)
  table0a$CPO[i] <- -sum(log(model$cpo$cpo))
  table0a$Mlik[i] <- round(model$mlik[1],1)
  
}

## Pull coefficients 
coef1 <- data.table(Model  = c("Base model", "Poverty index", "% urban", "% Lenca: log", "% Lenca: bin",
                               "% Maya/C: log", "% Maya/C: bin", "Avg. TRI", "Avg. Elevation",
                               "Cooking: % wood", "Water: % piped", "Water: % well", "Water: % river",
                               "% Maya/C: RINT","% Lenca: RINT","Cooking: wood bin"),
                    version = "v1 - without dist_km",
                    est = NA,
                    lo = NA,
                    hi = NA,
                    sig_90 = NA)

for(i in 2:length(formulas)){
  load(paste0("fitting_models/",lab[i],".RData"))
  coef1$est[i] <- round(exp(model$summary.fixed$mean[2]), 2)
  coef1$lo[i] <- round(exp(inla.qmarginal(0.05, marginal = model$marginals.fixed[[2]])),2)
  coef1$hi[i] <- round(exp(inla.qmarginal(0.95, marginal = model$marginals.fixed[[2]])),2)
  coef1$sig_90[i] <- ifelse(between(1,
                                    round(exp(inla.qmarginal(0.05, marginal = model$marginals.fixed[[2]])),2),
                                    round(exp(inla.qmarginal(0.95, marginal = model$marginals.fixed[[2]])),2))
                            ,0,1)
}

## B. WITH measure of distance to clinic ------
baseformula <- n1 ~ 1 + dist_km

# define formulas by updating the baseline formula
f1.1 <- update.formula(baseformula, ~. + poverty_index)
f1.2 <- update.formula(baseformula, ~. + urban_pct)
f1.3 <- update.formula(baseformula, ~. + len_log)
f1.4 <- update.formula(baseformula, ~. + len_bin)
f1.5 <- update.formula(baseformula, ~. + may_log)
f1.6 <- update.formula(baseformula, ~. + may_bin)
f1.7 <- update.formula(baseformula, ~. + avg_tri)
f1.8 <- update.formula(baseformula, ~. + avg_elev)
f1.9 <- update.formula(baseformula, ~. + wood_pct)
f1.10 <- update.formula(baseformula, ~. + water_piped)
f1.11 <- update.formula(baseformula, ~. + water_well)
f1.12 <- update.formula(baseformula, ~. + water_river)
f1.13 <- update.formula(baseformula, ~. + may_rint)
f1.14 <- update.formula(baseformula, ~. + len_rint)
f1.15 <- update.formula(baseformula, ~. + wood_bin)

# create a list of formulas
formulas <- list(baseformula, f1.1, f1.2, f1.3, f1.4, f1.5,
                 f1.6, f1.7, f1.8, f1.9, f1.10, f1.11, f1.12,
                 f1.13,f1.14,f1.15)

# create model label string
lab <- paste0("a_v2_",c("basemodel", "poverty_index", "urb_pct", "len_log", "len_bin",
                        "may_log", "may_bin", "avg_TRI", "avg_elev",
                        "pct_Woodcooking", "pct_Pipedwater", "pct_Unprotectedwell", "pct_Riverwater",
                        "may_rint","len_rint","wood_bin"))

models <- lapply(1:length(formulas), 
                 function(i) {
                   model <- mymodel(formulas[[i]], x2)
                   save(model, file = paste0("fitting_models/", lab[i],".RData"))})


# Create table of fit statistics
table0b <- data.table(Model  = c("Base model", "Poverty index", "% urban", "% Lenca: log", "% Lenca: bin",
                                 "% Maya/C: log", "% Maya/C: bin", "Avg. TRI", "Avg. Elevation",
                                 "Cooking: % wood", "Water: % piped", "Water: % well", "Water: % river",
                                 "% Maya/C: RINT","% Lenca: RINT","Cooking: wood bin"),
                      version = "v2 - with dist_km",
                      DIC = NA,
                      WAIC = NA,
                      CPO = NA,
                      Mlik = NA)

for(i in 1:length(formulas)){
  load(paste0("fitting_models/",lab[i],".RData"))
  table0b$DIC[i] <- round(model$dic$dic, 1)
  table0b$WAIC[i] <- round(model$waic$waic,1)
  table0b$CPO[i] <- -sum(log(model$cpo$cpo))
  table0b$Mlik[i] <- round(model$mlik[1],1)
  
}

## Pull coefficients 
coef2 <- data.table(Model  = c("Base model", "Poverty index", "% urban", "% Lenca: log", "% Lenca: bin",
                               "% Maya/C: log", "% Maya/C: bin", "Avg. TRI", "Avg. Elevation",
                               "Cooking: % wood", "Water: % piped", "Water: % well", "Water: % river",
                               "% Maya/C: RINT","% Lenca: RINT","Cooking: wood bin"),
                   version = "v2 - with dist_km",
                   est = NA,
                   lo = NA,
                   hi = NA,
                   sig_90 = NA)

for(i in 2:length(formulas)){
  load(paste0("fitting_models/",lab[i],".RData"))
  coef2$est[i] <- round(exp(model$summary.fixed$mean[3]), 2)
  coef2$lo[i] <- round(exp(inla.qmarginal(0.05, marginal = model$marginals.fixed[[3]])),2)
  coef2$hi[i] <- round(exp(inla.qmarginal(0.95, marginal = model$marginals.fixed[[3]])),2)
  coef2$sig_90[i] <- ifelse(between(1,
                                   round(exp(inla.qmarginal(0.05, marginal = model$marginals.fixed[[3]])),2),
                                   round(exp(inla.qmarginal(0.95, marginal = model$marginals.fixed[[3]])),2))
                           ,0,1)
}

# 2 - Compare fit and coefficients (determine whether dist_km is important) -----
## 'Distance to WRH' improves WAIC, DIC, and marginal likelihood values
tab <- left_join(table0a %>% select(-version), 
                 table0b %>% select(-version) %>% 
                   rename(DIC2=DIC,WAIC2=WAIC,CPO2=CPO,Mlik2=Mlik)) %>%
  mutate(DIC_d = DIC-DIC2,
         WAIC_d = WAIC-WAIC2,
         CPO_d = abs(CPO) -abs(CPO2),
         BayesFact = exp(Mlik - Mlik2)) %>%
  select(Model, DIC, DIC2, DIC_d, WAIC, WAIC2, WAIC_d, CPO,CPO2,CPO_d,Mlik,Mlik2,BayesFact)


## Coefficient plots
## Together + plot: 
co <- bind_rows(coef1,coef2)

ggplot(co, aes(x = est, y=factor(Model), fill = version, color = version)) + 
  geom_errorbar(aes(xmin = lo, xmax = hi),position = "dodge", lwd = 1) + 
  geom_vline(xintercept= 1,linetype = "dashed") 

ggplot(co %>% filter(version == "v2 - with dist_km"),
       aes(x = est, y=factor(Model), fill = version, color = version)) + 
  geom_errorbar(aes(xmin = lo, xmax = hi),position = "dodge", lwd = 1) + 
  geom_vline(xintercept= 1,linetype = "dashed") 

ggplot(co %>% filter(version == "v1 - without dist_km"),
       aes(x = est, y=factor(Model), fill = version, color = version)) + 
  geom_errorbar(aes(xmin = lo, xmax = hi),position = "dodge", lwd = 1) + 
  geom_vline(xintercept= 1,linetype = "dashed") 


