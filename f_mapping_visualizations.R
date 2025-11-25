# Beach et al. 2025 - "Geographic dimensions of gastric cancer risk in western Honduras: A spatial ecological analysis"

## F. Map visualizations: relative risks, adjusted spatial relative risks, variables, supplementary materials

# Outline:
  ## 1. Access values from posterior distributions
  ## 2. Figure 4 - Pr(θi > 1) maps
  ## 3. Figure 5 - Pr(RR > 1) maps
  ## 4. Figure 2 - Covariates across geography
  ## 5. Supplementary figures B1-B8

# 0 - Load packages and data --------------

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
library(svglite)

x <- read_rds("data/mun_ds.rds") %>%
  mutate(len_bin = ifelse(round(len_pct)>=5,1,0),
         may_bin = ifelse(round(may_pct)>=5,1,0))

# 1. Extract posterior values --------
flist <- list.files("final_models")
flist <- flist[str_detect(flist, paste(c("v08","v09","v10","v11","v12","v13","v14"),collapse = '|'))]

outs <- paste0("n",parse_number(flist))
exps <-  paste0("e",parse_number(flist))

x3 <- x

for(i in 1:length(flist)){
  
  load(paste0("final_models/",flist[i]))
  
  strata <- paste0("n",parse_number(substr(flist[i],5,7)))
  mod <- substr(flist[i],nchar(flist[i])-8,nchar(flist[i])-6)
  
  x3$strata <- strata
  x3$model <- mod
  
  # Residual RR = exp(b_i)
  b_sum <- model$summary.random$spat_re                        
  b_margs <- model$marginals.random$spat_re   
  
  thr <- log(1)
  Pr_RRresid_gt1.00 <- vapply(b_margs, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  thr <- log(1.1)
  Pr_RRresid_gt1.10 <- vapply(b_margs, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  thr <- log(1.25)
  Pr_RRresid_gt1.25 <- vapply(b_margs, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  thr <- log(1.5)
  Pr_RRresid_gt1.50 <- vapply(b_margs, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  RR_resid_mean <- exp(b_sum$mean)
  RR_resid_med <- exp(b_sum$`0.5quant`)
  lci <- exp(vapply(b_margs, function(m) inla.qmarginal(0.05, m), numeric(1)))
  uci <- exp(vapply(b_margs, function(m) inla.qmarginal(0.90, m), numeric(1)))
  
  # Posterior mean predicted RR
  RR_pred_mean <- exp(model$summary.linear.predictor$mean)
  
  #Exceedance P(RR_pred > 1.5) = P(eta > log(1.5))
  thr <- log(1.5)
  Pr_RRpred_gt15 <- vapply(model$marginals.linear.predictor, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  thr <- log(1)
  Pr_RRpred_gt1 <- vapply(model$marginals.linear.predictor, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  thr <- log(2)
  Pr_RRpred_gt2 <- vapply(model$marginals.linear.predictor, function(m)
    1 - inla.pmarginal(q = thr, marginal = m), numeric(1))
  
  x3$pr_RRresid_gt1.00 <- Pr_RRresid_gt1.00[1:78]
  x3$pr_RRresid_gt1.10 <- Pr_RRresid_gt1.10[1:78]
  x3$pr_RRresid_gt1.25 <- Pr_RRresid_gt1.25[1:78]
  x3$pr_RRresid_gt1.50 <- Pr_RRresid_gt1.50[1:78]
  x3$RRresid_mean <- RR_resid_mean[1:78]
  x3$RRresid_med <- RR_resid_med[1:78]
  x3$RRresid_lci <- lci[1:78]
  x3$RRresid_uci <- uci[1:78]
  x3$RRpred_mean <- RR_pred_mean
  x3$pr_RRpred_gt15 <- Pr_RRpred_gt15
  x3$pr_RRpred_gt1 <- Pr_RRpred_gt1
  x3$pr_RRpred_gt2 <- Pr_RRpred_gt2
  
  
  write_rds(x3,paste0("spatialRR/",strata,"_",mod,"_spatialRR.rds"))
  
}


# ~~PREP FOR MAPPING~~ --------
rm(list=ls())

files <- list.files("spatialRR")
files <- files[str_detect(files, "v08")]

x2 <- bind_rows(read_rds(paste0("spatialRR/",files[1])),
                read_rds(paste0("spatialRR/",files[2])),
                read_rds(paste0("spatialRR/",files[3])),
                read_rds(paste0("spatialRR/",files[4])),
                read_rds(paste0("spatialRR/",files[5])),
                read_rds(paste0("spatialRR/",files[6])),
                read_rds(paste0("spatialRR/",files[7])),
                read_rds(paste0("spatialRR/",files[8])),
                read_rds(paste0("spatialRR/",files[9])),
                read_rds(paste0("spatialRR/",files[10])),
                read_rds(paste0("spatialRR/",files[11])),
                read_rds(paste0("spatialRR/",files[12]))) %>% 
  mutate(pop = factor(case_when(strata %in% c("n1","n4","n7","n10") ~ "All",
                         strata %in% c("n2","n5","n8","n11") ~ "Male",
                         strata %in% c("n3","n6","n9","n12") ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c("n1","n2","n3") ~ "All GC",
                                    strata %in% c("n4","n5","n6") ~ "Diffuse",
                                    strata %in% c("n7","n8","n9") ~ "Intestinal",
                                    strata %in% c("n10","n11","n12") ~ "Mix/Indet."),
                          levels = c("All GC","Diffuse","Intestinal","Mix/Indet.")))

dep <- st_read("data/hnd_dep_clean.shp") 
clin <- st_read("data/clinics.shp") %>%
  filter(n=="WRH")
cols <- viridis::viridis(4)
cols2 <- cols[2:3]

# 2. Figure 4 - Pr(θi > 1.0) maps -------
x2$bin <- cut(x2$pr_RRpred_gt1,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold("Pr(" ~ bolditalic(theta[i]) ~ "> 1.0)"))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/fig4_pr1exc.png"), 
       width = 8.5,
       height = 9.1,
       units = "in")
# ggsave(b,
#        device = svglite::svglite,
#        filename = paste0("outputs/fig4_pr1exc.svg"), 
#        width = 8.5,
#        height = 9.1,
#        units = "in")

# 3. Figure 5 - Pr(RR > 1) maps -------
rm(list=ls())

files <- list.files("spatialRR")
files <- files[str_detect(files, "v13")]

x2 <- bind_rows(read_rds(paste0("spatialRR/",files[1])),
                read_rds(paste0("spatialRR/",files[2])),
                read_rds(paste0("spatialRR/",files[3])),
                read_rds(paste0("spatialRR/",files[4])),
                read_rds(paste0("spatialRR/",files[5])),
                read_rds(paste0("spatialRR/",files[6])),
                read_rds(paste0("spatialRR/",files[7])),
                read_rds(paste0("spatialRR/",files[8])),
                read_rds(paste0("spatialRR/",files[9])),
                read_rds(paste0("spatialRR/",files[10])),
                read_rds(paste0("spatialRR/",files[11])),
                read_rds(paste0("spatialRR/",files[12]))) %>% 
  mutate(pop = factor(case_when(strata %in% c("n1","n4","n7","n10") ~ "All",
                                strata %in% c("n2","n5","n8","n11") ~ "Male",
                                strata %in% c("n3","n6","n9","n12") ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c("n1","n2","n3") ~ "Gastric cancer",
                                    strata %in% c("n4","n5","n6") ~ "Diffuse",
                                    strata %in% c("n7","n8","n9") ~ "Intestinal",
                                    strata %in% c("n10","n11","n12") ~ "Mix/Indet."),
                          levels = c("Gastric cancer","Diffuse","Intestinal","Mix/Indet.")))

dep <- st_read("data/hnd_dep_clean.shp") 
clin <- st_read("data/clinics.shp") %>%
  filter(n=="WRH")
cols <- viridis::viridis(4)
cols2 <- cols[2:3]
#"#440154FF" "#31688EFF" "#35B779FF" "#FDE725FF"

## A col ------
x2$RR_cen <- ifelse(x2$RRresid_med > 2,2, x2$RRresid_med)

a <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "All"), 
          aes(fill = RR_cen),
          col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Med. RR"[i]))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size = 7.5,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.text = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none")

extract_fill_legend <- function(p) {
  p_only_fill <- p + guides(
    colour   = "none",
    alpha    = "none",
    size     = "none",
    linetype = "none",
    shape    = "none",
    stroke   = "none"
    # leave fill alone → keeps legend for scale_fill_*
  )
  cowplot::get_legend(p_only_fill)
}

leg1 <- extract_fill_legend(a) 

a <- a + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())

## B col -----
x2$bin <- cut(x2$pr_RRresid_gt1.00,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)
b <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "All"),
          aes(fill = bin), col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_manual(values = cols,
                    name = expression(bold(paste("Pr(RR"[i], " > 1)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    #legend.justification=c(0,-6),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size =8,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.placement = "outside",            
    #strip.text.y = element_text(size=11),
    strip.text.x = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none",
         fill = guide_legend(reverse=T))

leg2 <- extract_fill_legend(b) 

legend_col <- cowplot::plot_grid(
  leg1, leg2,
  ncol = 1, 
  align = "v", 
  axis = "l",
  rel_heights = c(1, 1))

b <- b + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())


## C col ------
x2$RR_cen <- ifelse(x2$RRresid_med > 2,2, x2$RRresid_med)

c <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "Male"), 
          aes(fill = RR_cen),
          col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Med. RR"[i]))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size = 7.5,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.text = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none")

c <- c + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())

## D col -----
x2$bin <- cut(x2$pr_RRresid_gt1.00,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)
d <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "Male"),
          aes(fill = bin), col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_manual(values = cols,
                    name = expression(bold(paste("Pr(RR"[i], " > 1)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    #legend.justification=c(0,-6),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size =8,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.placement = "outside",            
    #strip.text.y = element_text(size=11),
    strip.text.x = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none")

d <- d + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())

## E col ------
x2$RR_cen <- ifelse(x2$RRresid_med > 2,2, x2$RRresid_med)

e <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "Female"), 
          aes(fill = RR_cen),
          col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Med. RR"[i]))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size = 7.5,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.text = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none")

e <- e + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())

## F col -----
x2$bin <- cut(x2$pr_RRresid_gt1.00,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)
f <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2 %>% 
            #filter(outcome != "Mix/Indet.") %>%
            filter(pop == "Female"),
          aes(fill = bin), col = NA) + 
  facet_wrap( ~ outcome,ncol=1) +
  scale_fill_manual(values = cols2,
                    name = expression(bold(paste("Pr(RR"[i], " > 1)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() +
  theme(#legend.key.height= unit(0.2, 'in'),
    legend.key.width= unit(0.3, 'cm'),
    #legend.justification=c(0,-6),
    legend.position="right",
    legend.direction="vertical",
    legend.title = element_text(size =8,face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size=7.5),
    plot.background = element_rect(colour = NA, fill="white",
                                   linewidth=1),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.border = element_rect(colour = "gray20", fill=NA, linewidth=0.7),
    strip.placement = "outside",            
    #strip.text.y = element_text(size=11),
    strip.text.x = element_blank(),
    strip.background = element_rect(color="black",
                                    fill=NA),
    panel.spacing = unit(0, "pt"),           
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()) +
  guides(colour = "none", size = "none", linetype = "none", shape = "none")

f <- f + theme(legend.position = "none",
               plot.margin = margin(0,0,0,0),
               panel.spacing = unit(0, "pt"),
               strip.text = element_blank())

## SAVE ------
# Save colums as vector graphics files, format in PPTx.
ggsave(legend_col,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/leg.svg"), # change title here
       width = 1,
       height = 4,
       units = "in")

ggsave(a,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/a.svg"), 
       width = 1.2,
       height = 6,
       units = "in")

ggsave(b,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/b.svg"), 
       width = 1.2,
       height = 6,
       units = "in")

ggsave(c,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/c.svg"), 
       width = 1.2,
       height = 6,
       units = "in")

ggsave(d,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/d.svg"), 
       width = 1.2,
       height = 6,
       units = "in")
ggsave(e,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/e.svg"), 
       width = 1.2,
       height = 6,
       units = "in")
ggsave(f,
       device = svglite::svglite,
       filename = paste0("outputs/fig5/f.svg"), 
       width = 1.2,
       height = 6,
       units = "in")

# 4. Figure 2 - Covariates -------
x <- read_rds("data/mun_ds.rds")

a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x, aes(fill = avg_elev), col = NA) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(250,500,1000,1500,2000), 
                        labels =  c(250,500,1000,1500,2000), 
                        limits = c(250,2000),
                        "Avg. Elevation (m)") +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 



b <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x, aes(fill = len_pct), col = NA) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0,25,50,75,100), 
                        labels = c(0,25,50,75,100),
                        limits = c(0,100),
                        "% Lenca") +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

c <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x, aes(fill = may_pct), col = NA) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0,25,50,75,100), 
                        labels = c(0,25,50,75,100),
                        limits = c(0,100),
                        "% Maya/Chortí") +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 


d <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x, aes(fill = nbi_avg), col = NA) +
  scale_fill_continuous(type = "viridis", 
                        breaks = c(0.5,1.0,1.5,2,2.2), 
                        labels = c("0.5","1.0","1.5","2.0",""),
                        limits = c(0.5,2.2),
                        "Poverty Index") +
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 


fig2 <- (a|b | c)  + plot_layout(guides = "keep") &
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    linewidth = 1),
        legend.position = c(0.02, 0.99),        # X (horizontal), Y (vertical)
        legend.justification = c(0, 1),         # Align bottom-right corner
        legend.background = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 7),   # Shrink legend text
        legend.title = element_text(size = 9,face = "bold"),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.margin = margin(0, 1, 0.2, 0))



ggsave("outputs/figure2_vars.png", plot = fig2, width = 8.5, height = 5, dpi = 300)
ggsave(fig2,
       device = svglite::svglite,
       filename = paste0("outputs/figure2_vars.svg"), # change title here
       width = 8.5, height = 5, unit = "in")


# 5. SUPPLEMENTARY FIGURES --------

## Figure B1 - Posterior Median θi -------
rm(list=ls())

files <- list.files("spatialRR")
files <- files[str_detect(files, "v08")]

x2 <- bind_rows(read_rds(paste0("spatialRR/",files[1])),
                read_rds(paste0("spatialRR/",files[2])),
                read_rds(paste0("spatialRR/",files[3])),
                read_rds(paste0("spatialRR/",files[4])),
                read_rds(paste0("spatialRR/",files[5])),
                read_rds(paste0("spatialRR/",files[6])),
                read_rds(paste0("spatialRR/",files[7])),
                read_rds(paste0("spatialRR/",files[8])),
                read_rds(paste0("spatialRR/",files[9])),
                read_rds(paste0("spatialRR/",files[10])),
                read_rds(paste0("spatialRR/",files[11])),
                read_rds(paste0("spatialRR/",files[12]))) %>% 
  mutate(pop = factor(case_when(strata %in% c("n1","n4","n7","n10") ~ "All",
                                strata %in% c("n2","n5","n8","n11") ~ "Male",
                                strata %in% c("n3","n6","n9","n12") ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c("n1","n2","n3") ~ "Gastric cancer",
                                    strata %in% c("n4","n5","n6") ~ "Diffuse",
                                    strata %in% c("n7","n8","n9") ~ "Intestinal",
                                    strata %in% c("n10","n11","n12") ~ "Mix/Indet."),
                          levels = c("Gastric cancer","Diffuse","Intestinal","Mix/Indet.")))

dep <- st_read("data/hnd_dep_clean.shp") 
clin <- st_read("data/clinics.shp") %>%
  filter(n=="WRH")
cols <- viridis::viridis(4)
cols2 <- cols[2:3]

x2$RR_cen <- ifelse(x2$RRpred_mean > 2,2, x2$RRpred_mean)

a <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2, 
          aes(fill = RR_cen),
          col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Mean posterior "~theta[i]))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
        legend.title.position = "left",
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
       filename = paste0("outputs/supplementary/figB1_predMEAN.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")


## Figure B2 - Pr(Theta > 1.5) --------
x2$bin <- cut(x2$pr_RRpred_gt15,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold("Pr(" ~ bolditalic(theta[i]) ~ "> 1.5)"))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB2_pr1.5exc.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")

## Figure B3 - Pr(θi > 2.0) maps -------
x2$bin <- cut(x2$pr_RRpred_gt2,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold("Pr(" ~ bolditalic(theta[i]) ~ "> 2.0)"))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB3_prPREDexc2.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")




## Figure B4 - Residual RR LCI -------

rm(list=ls())

files <- list.files("spatialRR")
files <- files[str_detect(files, "v13")]

x2 <- bind_rows(read_rds(paste0("spatialRR/",files[1])),
                read_rds(paste0("spatialRR/",files[2])),
                read_rds(paste0("spatialRR/",files[3])),
                read_rds(paste0("spatialRR/",files[4])),
                read_rds(paste0("spatialRR/",files[5])),
                read_rds(paste0("spatialRR/",files[6])),
                read_rds(paste0("spatialRR/",files[7])),
                read_rds(paste0("spatialRR/",files[8])),
                read_rds(paste0("spatialRR/",files[9])),
                read_rds(paste0("spatialRR/",files[10])),
                read_rds(paste0("spatialRR/",files[11])),
                read_rds(paste0("spatialRR/",files[12]))) %>% 
  mutate(pop = factor(case_when(strata %in% c("n1","n4","n7","n10") ~ "All",
                                strata %in% c("n2","n5","n8","n11") ~ "Male",
                                strata %in% c("n3","n6","n9","n12") ~ "Female"),
                      levels = c("All","Male","Female")),
         outcome = factor(case_when(strata %in% c("n1","n2","n3") ~ "Gastric cancer",
                                    strata %in% c("n4","n5","n6") ~ "Diffuse",
                                    strata %in% c("n7","n8","n9") ~ "Intestinal",
                                    strata %in% c("n10","n11","n12") ~ "Mix/Indet."),
                          levels = c("Gastric cancer","Diffuse","Intestinal","Mix/Indet.")))

dep <- st_read("data/hnd_dep_clean.shp") 
clin <- st_read("data/clinics.shp") %>%
  filter(n=="WRH")
cols <- viridis::viridis(4)
cols2 <- cols[2:3]
#"#440154FF" "#31688EFF" "#35B779FF" "#FDE725FF"

x2$RR_cen <- ifelse(x2$RRresid_lci > 2,2, x2$RRresid_lci)

a <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2, 
          aes(fill = RR_cen),
          col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Posterior RR"[i]," - 90% CrI lower bound"))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB4_RRresidLCI.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")


## Figure B5 - Residual RR UCI -------

x2$RR_cen <- ifelse(x2$RRresid_uci > 2,2, x2$RRresid_uci)

a <- ggplot() + 
  geom_sf(data=dep,fill="gray75",col=NA,lwd=0.3)+
  geom_sf(data = x2, 
          aes(fill = RR_cen),
          col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_gradient2(
    low = "#2166AC",mid = "#F7F7F7", high = "#B2182B",  
    midpoint = 1,      
    name = expression(bold(paste("Posterior RR"[i]," - 90% CrI upper bound"))),
    breaks = c(0,0.5,1,1.5,2),
    labels = c("0","0.5","1","1.5",">2"),
    limits = c(0,2)
  ) +
  geom_sf(data=dep,col="gray30",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB5_RRresidUCI.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")

## Figure B6 - Pr(RR > 1.10) maps -------

x2$bin <- cut(x2$pr_RRresid_gt1.10,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold(paste("Pr(RR"[i], " > 1.10)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 


b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB6_prRRexc1.10.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")

## Figure B7 - Pr(RR > 1.25) maps -------

x2$bin <- cut(x2$pr_RRresid_gt1.25,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold(paste("Pr(RR"[i], " > 1.25)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB7_prRRexc1.25.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")


## Figure B8 - Pr(RR > 1.50) maps -------

x2$bin <- cut(x2$pr_RRresid_gt1.50,                     
              breaks = c(0,.05,.8,.95,1),
              include.lowest=T)


a <- ggplot() + 
  geom_sf(data=dep,fill="gray90",col=NA,lwd=0.3)+
  geom_sf(data = x2, aes(fill = bin), col = NA) + 
  facet_grid(pop ~ outcome) +
  scale_fill_manual(values = cols,
                    name = expression(bold(paste("Pr(RR"[i], " > 1.5)")))) + 
  geom_sf(data=dep,col="gray75",fill=NA,lwd=0.3) +
  geom_sf(data=clin, shape=21,fill="white",
          col= "black",stroke = 1,
          show.legend = F, size=2.5) +
  coord_sf(xlim = c(243463.4,366230.4),
           ylim = c(1543163.2,1710991.1)) +
  theme_map() 

b <- egg::tag_facet(a)+
  theme(legend.key.height= unit(0.2, 'in'),
        legend.key.width= unit(0.7, 'in'),
        legend.justification=c(0,-6),
        legend.position="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 13,face = "bold"),
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
       filename = paste0("outputs/supplementary/figB8_prRRexc1.50.png"), # change title here
       width = 8.5,
       height = 9.1,
       units = "in")

