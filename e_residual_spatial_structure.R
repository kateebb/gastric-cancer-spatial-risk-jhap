# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

## E. Explore spatial structure of residual variation

## Table 3 creation
## Bonus explorations: plots of posterior phi and precision for BYM2 random effect (not included in supplement)

# 0. Load packages and data ------
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

nb <- poly2nb(x)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 1. Extract values ----
### FE model ------
files <- list.files("final_models/")
flist <- files[str_detect(files,"v06")]

coords <- st_coordinates(st_transform(st_centroid(x),4326))

## Assess Pearson residuals: Moran's I
grp <- c()
Moran <- c()
MoranSig <- c()

for(i in c(1:12)){
  
  load(paste0("final_models/",flist[i]))
  
  s <- paste0("n",parse_number(substr(flist[i],5,7)))
  e <- paste0("e",parse_number(substr(flist[i],5,7)))
  grp <- c(grp,s)
  
  # Pearson residuals
  y <- x2[,s]
  E <- x2[,e]
  rr_hat <- model$summary.fitted.values$mean  
  mu_hat <- E * rr_hat
  res_p <- (y - mu_hat) / sqrt(mu_hat)
  
  # Global Moran's I of pearson residuals
  z <-  moran.test(res_p, lw, zero.policy = TRUE,randomisation = TRUE)
  Moran <- c(Moran,round(unname(z$estimate[1]),2))
  MoranSig <- c(MoranSig,round(z$p.value,3))
  
}

fit <- data.frame(grp,Moran,MoranSig)

### Spatial model ------
files <- list.files("final_models/")
flist <- files[str_detect(files,"v13")]

## Posterior Phi & pairwise RE and FE correlation
grp <- c()
Phi_est <- c()
corr_dist <- c()
corr_len <- c()
corr_may <- c()
corr_pov <- c()

for(i in c(1:12)){
  load(paste0("final_models/",flist[i]))
  
  s <- paste0("n",parse_number(substr(flist[i],5,7)))
  e <- paste0("e",parse_number(substr(flist[i],5,7)))
  grp <- c(grp,s)
  
  res_p <- model$summary.random$spat_re$mean[1:78]
  
  # Posterior mean Phi
  Phi_est <- c(Phi_est,
               paste0(round(model$summary.hyperpar$`0.5quant`[2],2),
               " [",round(inla.qmarginal(0.05, marginal = model$marginals.hyperpar$`Phi for spat_re`),2),", ",
               round(inla.qmarginal(0.95, marginal = model$marginals.hyperpar$`Phi for spat_re`),2),"]"))
  
  # DIST
  test <- cor.test(x = as.vector(x2[,"dist_km"]), y = res_p,
                   method = "pearson")
  corr_dist <- c(corr_dist,
                 paste0(sprintf("%.2f",round(unname(test$estimate),2)),
                        case_when(unname(test$p.value) < 0.05 ~ "*",TRUE~"")))
  
  # LENCA
  test <- cor.test(x = as.vector(x2[,"len_bin"]), y = res_p,
                   method = "pearson")
  corr_len <- c(corr_len,
                paste0(sprintf("%.2f",round(unname(test$estimate),2)),
                       case_when(unname(test$p.value) < 0.05 ~ "*",TRUE~"")))
  
  
  # MAYA
  test <- cor.test(x = as.vector(x2[,"may_bin"]), y = res_p,
                   method = "pearson")
  corr_may <- c(corr_may,
                paste0(sprintf("%.2f",round(unname(test$estimate),2)),
                       case_when(unname(test$p.value) < 0.05 ~ "*",TRUE~"")))
  
  
  
  # POVERTY INDEX
  test <- cor.test(x = as.vector(x2[,"poverty_index"]), y = res_p,
                   method = "pearson")
  corr_pov <- c(corr_pov,
                paste0(sprintf("%.2f",round(unname(test$estimate),2)),
                       case_when(unname(test$p.value) < 0.05 ~ "*",TRUE~"")))
  
                      
  
}


fit2 <- data.frame(grp,Phi_est,corr_dist, corr_may, corr_len, corr_pov)


# TABLE 3 -------
table3 <- left_join(fit,
                    fit2,
                    by = "grp") %>%
  mutate(sig = case_when(MoranSig < 0.05 ~  "**",
                         MoranSig < 0.1 ~ "*",
                         TRUE ~ ""),
         mor = paste0(Moran,sig),
         strata = parse_number(grp)) %>% 
  mutate(pop = factor(case_when(strata %in% c(1,4,7,10) ~ "All",
                                 strata %in% c(2,5,8,11) ~ "Male",
                                 strata %in% c(3,6,9,12) ~ "Female"),
                       levels = c("All","Male","Female")),
          outcome = factor(case_when(strata %in% c(1:3) ~ "All GC",
                                     strata %in% c(4:6) ~ "Diffuse",
                                     strata %in% c(7:9) ~ "Intestinal",
                                     strata %in% c(10:12) ~ "Mixed"),
                           levels = c("All GC","Diffuse","Intestinal","Mixed"))) %>% 
  arrange(pop,outcome) %>%
  select(pop,outcome,mor,Phi_est,
         corr_dist,corr_len,corr_may,corr_pov)

write_csv(table3,"outputs/table3.csv")


# BONUS: LISA maps of Pearson Residuals from FE models -----
library(tidyverse)
library(ggplot2)
library(sf)
library(INLA)
library(spdep)
library(sp)
library(data.table)
library(ggthemes)
library(patchwork)

rm(list=ls())
x <- read_rds("data/mun_ds.rds") 
x2 <- x %>% st_set_geometry(NULL)
nb <- poly2nb(x)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

files <- list.files("final_models/")
flist <- files[str_detect(files,"v06")]


for(i in c(1:12)){
  
  load(paste0("final_models/",flist[i]))
  
  s <- paste0("n",parse_number(substr(flist[i],5,7)))
  e <- paste0("e",parse_number(substr(flist[i],5,7)))
  
  # Pearson residuals
  y <- x2[,s]
  E <- x2[,e]
  rr_hat <- model$summary.fitted.values$mean  
  mu_hat <- E * rr_hat
  res_p <- (y - mu_hat) / sqrt(mu_hat)
  lag_res <- lag.listw(lw, res_p, zero.policy = TRUE)
  
  # Global Moran's I of pearson residuals
  set.seed(123)
  loc <- as.data.frame(
    localmoran_perm(res_p, lw, nsim = 999,
                    alternative = "two.sided", zero.policy = TRUE))
  
  pval  <- loc[["Pr(z != E(Ii)) Sim"]]
  p_adj <- p.adjust(pval, method = "BH")
  
  # --- 4) LISA cluster classification ---
  alpha <- 0.05
  sig   <- pval < alpha
  quad <- case_when(
    sig & res_p > 0 & lag_res > 0 ~ "High–High",
    sig & res_p < 0 & lag_res < 0 ~ "Low–Low",
    sig & res_p > 0 & lag_res < 0 ~ "High–Low",
    sig & res_p < 0 & lag_res > 0 ~ "Low–High",
    TRUE                            ~ "Not significant"
  )
  quad <- factor(quad, levels = c("High–High","Low–Low","High–Low","Low–High","Not significant"))
  
  z <- x
  z$lisa <- quad
  z$grp <- s
  z$res_p <- res_p
  
  write_rds(z,paste0("residuals/",s,".rds"))
  
}


f <- list.files("residuals")

z <- map_dfr(paste0("residuals/",f), read_rds) %>%
  mutate(strata = parse_number(grp),
         pop = factor(case_when(strata %in% c(1,4,7,10) ~ "All",
                                strata %in% c(2,5,8,11) ~ "Male",
                                strata %in% c(3,6,9,12) ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c(1:3) ~ "All GC",
                                    strata %in% c(4:6) ~ "Diffuse",
                                    strata %in% c(7:9) ~ "Intestinal",
                                    strata %in% c(10:12) ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet.")))


## LISA MAPS -----
a <- ggplot(z) +
  geom_sf(aes(fill = lisa), color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("High–High"="#B2182B","Low–Low"="#2166AC",
                               "High–Low"="#EF8A62","Low–High"="#67A9CF",
                               "Not significant"="grey90")) +
  facet_grid(pop ~ outcome) +
  labs(title = "LISA of Pearson Residuals", subtitle = "From fully adjusted fixed effects models") +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map()

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        plot.background = element_rect(colour = NA, fill="white",
                                       linewidth=1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.caption =  element_text(size=11,face="bold",hjust=0,vjust=3),
        panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
        strip.text.x = element_text(size=11),
        strip.text.y = element_text(size=11),
        strip.background = element_rect(color="black",
                                        fill=NA))


ggsave(b,
       device = "png",
       filename = paste0("outputs/supplementary/lisamaps_fe_full.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")  


## PEARSON RESIDUALS MAP -----
a <- ggplot(z) +
  geom_sf(aes(fill = res_p), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits=c(-3.5,5)) +
  facet_grid(pop ~ outcome) +
  labs(title = "Pearson Residuals", subtitle = "Fully adjusted fixed effects models") +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map()

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        plot.background = element_rect(colour = NA, fill="white",
                                       linewidth=1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.caption =  element_text(size=11,face="bold",hjust=0,vjust=3),
        panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
        strip.text.x = element_text(size=11),
        strip.text.y = element_text(size=11),
        strip.background = element_rect(color="black",
                                        fill=NA))


ggsave(b,
       device = "png",
       filename = paste0("outputs/supplementary/presmaps_fe_full.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")  


# BONUS: Posterior Phi & Precision distributions ----
flist <- list.files("final_models")
flist <- flist[str_detect(flist, paste(c("v08","v09","v10","v11","v12","v13","v14"),collapse = '|'))]

dist2 <- data.frame(x=NA,y=NA,strat=NA,model=NA)

for(i in 1:length(flist)){
  
  load(paste0("final_models/",flist[i]))
  
  strata <- paste0("n",parse_number(substr(flist[i],5,7)))
  mod <- substr(flist[i],nchar(flist[i])-8,nchar(flist[i])-6)
  
  d1 <- model$marginals.hyperpar$`Precision for spat_re` %>%
    as.data.frame() %>%
    mutate(model = mod,
           strat = strata,
           co = "1. Precision (theta1)")
  d2 <- model$marginals.hyperpar$`Phi for spat_re` %>%
    as.data.frame() %>%
    mutate(model = mod,
           strat = strata,
           co = "2. Phi (theta2)")
  
  dist2 <- bind_rows(dist2,d1,d2) %>% na.omit() 
  
}

dist2 <- dist2 %>%
  mutate(pop = factor(case_when(strat %in% c("n1","n4","n7","n10") ~ "All",
                                strat %in% c("n2","n5","n8","n11") ~ "Male",
                                strat %in% c("n3","n6","n9","n12") ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strat %in% c("n1","n2","n3") ~ "All GC",
                                    strat %in% c("n4","n5","n6") ~ "Diffuse",
                                    strat %in% c("n7","n8","n9") ~ "Intestinal",
                                    strat %in% c("n10","n11","n12") ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet."))) %>%
  mutate(model2 = factor(case_when(model == "v08" ~ "1. Null",
                                   model == "v09" ~ "2. Dist. to WRH",
                                   model == "v10" ~ "2. % Maya/Chortí",
                                   model == "v11" ~ "2. % Lenca",
                                   model == "v12" ~ "2. Pov. Index",
                                   model == "v13" ~ "3. Fully adj.",
                                   model == "v14" ~ "4. Fully adj. [-dist.]"),
                         levels = c("1. Null","3. Fully adj.",
                                    "4. Fully adj. [-dist.]","2. Dist. to WRH",
                                    "2. % Maya/Chortí","2. % Lenca",
                                    "2. Pov. Index")))

# Phi
(b<-ggplot(dist2 %>% 
             filter(co == "2. Phi (theta2)") %>%
             filter(model=="v13"), 
           aes(x = x, y = y,group = model2,color=model2)) +
    facet_grid(pop~outcome,scales = "free")+
    geom_line() +
    scale_x_continuous(breaks = c(0,.25,.5,.75,1),
                       labels = c("0","","0.50","","1"),
                       limits=c(0,1)) +
    labs(title = "Posterior ϕ distributions",
         subtitle = "From fully adjusted spatial models",
         x = "ϕ",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "none"))

# Precision
(b<-ggplot(dist2 %>% 
             filter(co == "1. Precision (theta1)") %>%
             filter(model=="v13"), 
           aes(x = x, y = y,group = model2,color=model2)) +
    facet_grid(pop~outcome,scales = "free")+
    geom_line() +
    scale_x_continuous(#breaks = c(0,10,20,30,40,50),
                       #labels = c("0","","0.50","","1"),
                       limits=c(0,30)) +
    labs(title = "Posterior Precision distributions",
         subtitle = "From fully adjusted spatial models",
         x = "Precision",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "none"))



