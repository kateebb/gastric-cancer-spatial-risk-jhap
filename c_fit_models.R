# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

## C. Fit series of models for gastric cancer, for each of the 12 sex/sub-types

# Fit fully adjusted and reduced models to check coefficient stability across specifications   
  
  ## v01: FE - Null                 
  ## v02: FE + distance to WRH
  ## v03: FE + % Maya/Chorti
  ## v04: FE + % Lenca
  ## v05: FE + poverty index
  ## v06: FE - Fully Adjusted        <-- Primary FE model presented in the manuscript (Figure 3)
  ## v07: *FE, full[-dist_km]

  ## v08: Spatial - Null.            <-- Model presented in Figure 4 (Observed Theta[i])
  ## v09: Spatial + distance to WRH
  ## v10: Spatial + % Maya/Chorti
  ## v11: Spatial + % Lenca
  ## v12: Spatial + poverty index
  ## v13: Spatial - Fully Adjusted   <-- Primary spatial model presented in the manuscript (Figures 3 & 5)
  ## v14: *Spatial, full[-dist_km]

## Outcome matrix - labels model objects by outcome stratum
  ## (n# = observed cases in the strata; 
  ##  e# = age/sex standardized expected cases in each strata)

# [GRP] [ALL] [DIF] [INT] [MIX]
# [ALL]   1     4     7     10
# [MAL]   2     5     8     11
# [FEM]   3     6     9     12

# 0 - Load packages and data, prep for analysis (scaling, transformations) --------

setwd("")

library(tidyverse)
library(ggplot2)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)
library(ggthemes)
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

## Random effects index
x2$spat_re <- c(1:78)
## Spatial random effects neighbors 
g <- inla.read.graph(filename = "data/munis.adj")

## Formulas for each model version

## FE models
f01 <- y ~ 1 
f02 <- y ~ 1 + dist_km 
f03 <- y ~ 1 + may_bin 
f04 <- y ~ 1 + len_bin 
f05 <- y ~ 1 + poverty_index 
f06 <- y ~ 1 + dist_km + poverty_index + len_bin + may_bin 
f07 <- y ~ 1 + poverty_index + len_bin + may_bin 

## Naive spatial (Pr(Phi > 0.5) = 0.5)
f08 <- y ~ 1 + 
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f09 <- y ~ 1 + dist_km +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f10 <- y ~ 1 + may_bin + 
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f11 <- y ~ 1 + len_bin + 
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f12 <- y ~ 1 + poverty_index + 
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f13 <- y ~ 1 + dist_km + poverty_index + len_bin + may_bin +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 
f14 <- y ~ 1 + poverty_index + len_bin + may_bin +
  f(spat_re, model="bym2",graph=g, scale.model = TRUE, constr = TRUE,
    hyper = list(  # BYM2 spatial
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),
      theta2 = list(prior = "pc", param = c(0.5, 0.5)))) 


flist <- list(f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14)

rm(i);rm(s)

# 1 - Models for each of the 12 strata ---------
  # Note: fitting all models will take a very long time. 
  # Select only Null FE (v01), Null Spatial (v08), Fully adjusted FE (v06), and fully adjusted spatial (v13) to save time
outs <- c("n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11","n12")
exps <- c("e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12")
mods <- c("v01","v02","v03","v04","v05","v06","v07","v08","v09","v10","v11",
          "v12","v13","v14")

x3 <- x2

for(i in c(1:12)){
  col <- outs[i]
  exp <- exps[i]
  x3$y <- x3[,col]
  x3$e <- x3[,exp]
  

  for(k in c(1:14)){
    
    model <- inla(formula = flist[[k]], family = "poisson", data = x3, E = e,
         control.predictor = list(compute = TRUE),
         control.compute = list(return.marginals.predictor = TRUE,
                                dic = TRUE, mlik = TRUE, waic =TRUE,
                                cpo=TRUE, config = TRUE),
         control.fixed = list(mean = 0, prec = 0.001))
    
    save(model, file = paste0("final_models/mod_",col,"_",mods[k],".RData"))
  }

}


# FIGURE 3: Compare coefficient stability across models ------
flist <- list.files("final_models")
flist <- flist[str_detect(flist, "v06") | str_detect(flist, "v13")]
strata <- paste0("n",c(1:12))
models <- c("v06","v13")

## Pull coefficients to compare differences
coef <- data.table(strata = NA,
                   model = NA,
                   term = NA,
                   est = NA,
                   lo = NA,
                   hi = NA,
                   sig = NA)


for(i in c(1:length(flist))){
  load(paste0("final_models/",flist[i]))
  
  
  # Posterior IRRs
  x3 <- data.frame(term = rownames(data.frame(model$summary.fixed)),
                   est = round(exp(model$summary.fixed$mean), 2),
                   lo = round(exp(model$summary.fixed$`0.025quant`),2),
                   hi = round(exp(model$summary.fixed$`0.975quant`),2))
  
  x3$strata <- substr(flist[i],5,7)
  x3$model <- substr(flist[i],nchar(flist[i])-8,nchar(flist[i])-6)
  
  
  x3 <- x3 %>%
    mutate(sig = ifelse(between(1,lo,hi),0,1))
  
  coef <- bind_rows(coef,x3)
  
}

tab <- coef %>%
  filter(term!="(Intercept)" & !is.na(term)) %>%
  mutate(strata = parse_number(strata)) %>%
  mutate(label = factor(case_when(term == "poverty_index" ~ "Poverty index",
                                  term == "may_bin" ~ "Maya/Chortí",
                                  term == "len_bin" ~ "Lenca",
                                  term == "dist_km" ~ "Dist. to WRH"),
                        levels = rev(c("Intercept",
                                       "Maya/Chortí",
                                       "Lenca","Poverty index",
                                       "Dist. to WRH"))),
         pop = factor(case_when(strata %in% c(1,4,7,10) ~ "All",
                                strata %in% c(2,5,8,11) ~ "Male",
                                strata %in% c(3,6,9,12) ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c(1:3) ~ "All GC",
                                    strata %in% c(4:6) ~ "Diffuse",
                                    strata %in% c(7:9) ~ "Intestinal",
                                    strata %in% c(10:12) ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet.")),
         modgrp = factor(case_when(model %in% c("v02","v03","v04","v05") ~ "FE, bivariate",
                                   model %in% c("v06") ~ "FE",
                                   model %in% c("v09","v10","v11","v12") ~ "Spatial, bivariate",
                                   model %in% c("v13") ~ "Spatial",
                                   model %in% c("v07") ~ "FE, full[-dist]*",
                                   model %in% c("v14") ~ "Spatial, full[-dist]*"),
                         levels = rev(c("FE, bivariate","FE",
                                        "Spatial, bivariate","Spatial",
                                        "FE, full[-dist]*","Spatial, full[-dist]*"))))


(fig_3 <- tab %>% 
    filter(label != "Dist. to WRH") %>%
    #filter(modgrp == "FE") %>%
  ggplot(aes(x = est, y = label,col=factor(sig),shape=factor(modgrp))) +
  facet_grid(pop ~ outcome) +
  scale_color_manual(values=c("gray40","red"),"Significance") +
  scale_shape_manual(values=c(2,17),"Model") + 
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.4), height = 0,
                 linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.4), size = 2.1) +
  geom_vline(xintercept = 1, linetype = "dashed",color="gray30") +
  geom_vline(xintercept = 2, linetype = "dashed",color="gray50",linewidth=0.3) +
    geom_vline(xintercept = 3, linetype = "dashed",color="gray50",linewidth=0.3) +
    geom_vline(xintercept = 4, linetype = "dashed",color="gray50",linewidth=0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed",color="gray50",linewidth=0.3) +
  scale_x_continuous(breaks = c(0.5,1,2,3,4),
                     labels = c(".5","1","2","3","4"),
                     limits = c(0.20,4.5))+
  labs(x = "IRR Estimate (95% CrI)", y = NULL) +
  theme_few(base_size = 11) +
  guides(color = guide_legend(title.position = "left",  title.hjust = 0.5),
           shape = guide_legend(title.position = "left",  title.hjust = 0.5)) + 
    theme(text = element_text(colour = "gray20"),
          panel.spacing      = unit(0.4, "lines"),
          strip.placement    = "outside",
          strip.text.y.right  = element_text(angle = 0, vjust = 0.5),
          strip.text         = element_text(face = "bold", color = "gray20",size=11),
          axis.title.y       = element_blank(),
          legend.position    = "bottom",
          legend.title = element_text(size = 10,face="bold"),
          plot.margin        = margin(2, 2, 2, 2),
          panel.grid.major.y = element_blank())) 

ggsave(fig_3,
       filename="outputs/figure3_coeff.png",
       device="png",
       width=7,height=4.5,units="in",dpi=300)

ggsave(fig_3,
       device = svglite::svglite,
       filename = paste0("outputs/figure3_coeff.svg"), 
       width=7,height=4.5,units="in",dpi=300)


# TABLE C1 ------
flist <- list.files("final_models")
flist <- flist[str_detect(flist, "v06")]
strata <- paste0("n",c(1:12))

## Pull coefficients to compare differences
coef <- data.table(strata = NA,
                   term = NA,
                   est = NA,
                   lo = NA,
                   hi = NA,
                   sig = NA)
crit <- data.frame(WAIC = NA,
                   DIC = NA,
                   Mlik = NA,
                   CPO = NA)


for(i in c(1:length(flist))){
  load(paste0("final_models/",flist[i]))
  
  
  # Posterior IRRs
  x3 <- data.frame(term = rownames(data.frame(model$summary.fixed)),
                   est = round(exp(model$summary.fixed$mean), 2),
                   lo = round(exp(model$summary.fixed$`0.025quant`),2),
                   hi = round(exp(model$summary.fixed$`0.975quant`),2))
  
  x3$strata <- substr(flist[i],5,7)
  
  
  x4 <- data.frame(strata = substr(flist[i],5,7),
             WAIC = round(model$waic$waic,1),
             DIC = round(model$dic$dic,1),
             Mlik = round(unname(model$mlik[1,1]),1),
             CPO = round(-sum(log(model$cpo$cpo), na.rm = TRUE),1))
  
  
  x3 <- x3 %>%
    mutate(sig = ifelse(between(1,lo,hi),0,1))
  
  coef <- bind_rows(coef,x3)
  crit <- bind_rows(crit,x4)
  
}

co <- coef %>%
  mutate(col = case_when(parse_number(strata) == 1 ~ "A",
                         parse_number(strata) == 4 ~ "B",
                         parse_number(strata) == 7 ~ "C",
                         parse_number(strata) == 10 ~ "D",
                         parse_number(strata) == 2 ~ "E",
                         parse_number(strata) == 5 ~ "F",
                         parse_number(strata) == 8 ~ "G",
                         parse_number(strata) == 11 ~ "H",
                         parse_number(strata) == 3 ~ "I",
                         parse_number(strata) == 6 ~ "J",
                         parse_number(strata) == 9 ~ "K",
                         parse_number(strata) == 12 ~ "L"),
         label = factor(case_when(term == "(Intercept)" ~ "Intercept",
                           term == "poverty_index" ~ "Poverty index",
                           term == "may_bin" ~ "Maya/Chortí",
                           term == "len_bin" ~ "Lenca",
                           term == "dist_km" ~ "Dist. to WRH"),
                        levels = c("Intercept","Poverty index",
                                   "Maya/Chortí","Lenca",
                                   "Dist. to WRH")),
         pop = factor(case_when(parse_number(strata) %in% c(1,4,7,10) ~ "All",
                                    parse_number(strata) %in% c(2,5,8,11) ~ "Male",
                                    parse_number(strata) %in% c(3,6,9,12) ~ "Female"),
                          levels = c("All","Male","Female")),
         outcome = factor(case_when(parse_number(strata) %in% c(1,2,3) ~ "All GC",
                             parse_number(strata) %in% c(4,5,6) ~ "Diffuse",
                             parse_number(strata) %in% c(7,8,9) ~ "Intestinal",
                             parse_number(strata) %in% c(10,11,12) ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet.")),
         row1 = sprintf("%.2f", est),
         row2 = paste0("[",sprintf("%.2f",lo),", ",sprintf("%.2f",hi),"]")) %>%
  select(col,strata,pop,outcome,label,row1,row2) %>% 
  pivot_longer(row1:row2, names_to = "lab", values_to = "val") %>% 
  na.omit() %>%
  pivot_wider(names_from = outcome, values_from = val) %>% 
  arrange(col,label)

cr <- crit %>% select(-DIC) %>%
  na.omit() %>%
  mutate(col = case_when(parse_number(strata) == 1 ~ "A",
                         parse_number(strata) == 4 ~ "B",
                         parse_number(strata) == 7 ~ "C",
                         parse_number(strata) == 10 ~ "D",
                         parse_number(strata) == 2 ~ "E",
                         parse_number(strata) == 5 ~ "F",
                         parse_number(strata) == 8 ~ "G",
                         parse_number(strata) == 11 ~ "H",
                         parse_number(strata) == 3 ~ "I",
                         parse_number(strata) == 6 ~ "J",
                         parse_number(strata) == 9 ~ "K",
                         parse_number(strata) == 12 ~ "L")) %>% 
  arrange(col) %>%
  t() %>% data.frame()
  

write_csv(co,"outputs/supplementary/tabC1_FEcoeffs.csv")
write_csv(cr,"outputs/supplementary/tabC1_fitcrit.csv")


# TABLE C2 ------
flist <- list.files("final_models")
flist <- flist[str_detect(flist, "v13")]
strata <- paste0("n",c(1:12))

## Pull coefficients to compare differences
coef <- data.table(strata = NA,
                   term = NA,
                   est = NA,
                   lo = NA,
                   hi = NA,
                   sig = NA)
crit <- data.frame(WAIC = NA,
                   DIC = NA,
                   Mlik = NA,
                   CPO = NA)


for(i in c(1:length(flist))){
  load(paste0("final_models/",flist[i]))
  
  
  # Posterior IRRs
  x3 <- data.frame(term = rownames(data.frame(model$summary.fixed)),
                   est = round(exp(model$summary.fixed$mean), 2),
                   lo = round(exp(model$summary.fixed$`0.025quant`),2),
                   hi = round(exp(model$summary.fixed$`0.975quant`),2))
  
  x3$strata <- substr(flist[i],5,7)
  
  
  x4 <- data.frame(strata = substr(flist[i],5,7),
                   WAIC = round(model$waic$waic,1),
                   DIC = round(model$dic$dic,1),
                   Mlik = round(unname(model$mlik[1,1]),1),
                   CPO = round(-sum(log(model$cpo$cpo), na.rm = TRUE),1))
  
  
  x3 <- x3 %>%
    mutate(sig = ifelse(between(1,lo,hi),0,1))
  
  coef <- bind_rows(coef,x3)
  crit <- bind_rows(crit,x4)
  
}

co <- coef %>%
  mutate(col = case_when(parse_number(strata) == 1 ~ "A",
                         parse_number(strata) == 4 ~ "B",
                         parse_number(strata) == 7 ~ "C",
                         parse_number(strata) == 10 ~ "D",
                         parse_number(strata) == 2 ~ "E",
                         parse_number(strata) == 5 ~ "F",
                         parse_number(strata) == 8 ~ "G",
                         parse_number(strata) == 11 ~ "H",
                         parse_number(strata) == 3 ~ "I",
                         parse_number(strata) == 6 ~ "J",
                         parse_number(strata) == 9 ~ "K",
                         parse_number(strata) == 12 ~ "L"),
         label = factor(case_when(term == "(Intercept)" ~ "Intercept",
                                  term == "poverty_index" ~ "Poverty index",
                                  term == "may_bin" ~ "Maya/Chortí",
                                  term == "len_bin" ~ "Lenca",
                                  term == "dist_km" ~ "Dist. to WRH"),
                        levels = c("Intercept","Poverty index",
                                   "Maya/Chortí","Lenca",
                                   "Dist. to WRH")),
         pop = factor(case_when(parse_number(strata) %in% c(1,4,7,10) ~ "All",
                                parse_number(strata) %in% c(2,5,8,11) ~ "Male",
                                parse_number(strata) %in% c(3,6,9,12) ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(parse_number(strata) %in% c(1,2,3) ~ "All GC",
                                    parse_number(strata) %in% c(4,5,6) ~ "Diffuse",
                                    parse_number(strata) %in% c(7,8,9) ~ "Intestinal",
                                    parse_number(strata) %in% c(10,11,12) ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet.")),
         row1 = sprintf("%.2f", est),
         row2 = paste0("[",sprintf("%.2f",lo),", ",sprintf("%.2f",hi),"]")) %>%
  select(col,strata,pop,outcome,label,row1,row2) %>% 
  pivot_longer(row1:row2, names_to = "lab", values_to = "val") %>% 
  na.omit() %>%
  pivot_wider(names_from = outcome, values_from = val) %>% 
  arrange(col,label)

cr <- crit %>% select(-DIC) %>%
  na.omit() %>%
  mutate(col = case_when(parse_number(strata) == 1 ~ "A",
                         parse_number(strata) == 4 ~ "B",
                         parse_number(strata) == 7 ~ "C",
                         parse_number(strata) == 10 ~ "D",
                         parse_number(strata) == 2 ~ "E",
                         parse_number(strata) == 5 ~ "F",
                         parse_number(strata) == 8 ~ "G",
                         parse_number(strata) == 11 ~ "H",
                         parse_number(strata) == 3 ~ "I",
                         parse_number(strata) == 6 ~ "J",
                         parse_number(strata) == 9 ~ "K",
                         parse_number(strata) == 12 ~ "L")) %>% 
  arrange(col) %>%
  t() %>% data.frame()


write_csv(co,"outputs/supplementary/tabC2_SpatialCoeffs.csv")
write_csv(cr,"outputs/supplementary/tabC2_fitcrit.csv")

# Figure A6: Map CPO distribution ----
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

## For primary models: fully adjusted FE and Spatial models
## All cancer, all persons stratum
load("final_models/mod_n1_v06.RData")
fe <- model
load("final_models/mod_n1_v13.RData")
re <- model
rm(model)

x1 <- x 
x2 <- x
x1$cpo <- -log(fe$cpo$cpo)
x1$grp <- "Fixed Effects Model"
x2$cpo <- -log(re$cpo$cpo)
x2$grp <- "Spatial Model"

x <- bind_rows(x1,x2)
x$grp <- factor(x$grp)

ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x, aes(fill = cpo), col = NA) + 
  facet_wrap(vars(grp),nrow=1) +
  scale_fill_continuous(name = "-log(CPO)") +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

nb <- poly2nb(x1)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# FE Model CPO: 0.33***
mi  <- moran.test(x1$cpo, lw, randomisation = TRUE, alternative = "greater")
mi$estimate["Moran I statistic"]; round(mi$p.value,3)

# Spatial Model CPO: 0.18***
mi  <- moran.test(x2$cpo, lw, randomisation = TRUE, alternative = "greater")
mi$estimate["Moran I statistic"]; round(mi$p.value,3)

  
