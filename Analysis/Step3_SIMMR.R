set.seed(1234)
# SIMMR (taxon specific) ------------------------------

## Setup ----------------------------
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)
library(lme4)
library(ggeffects)
library(patchwork)



## Visualization setup
order_colors = wes_palette("Darjeeling1", type = "continuous", n = 8)
order_levels = c("odonata", "ephemeroptera", "hemiptera", "gastropoda",
                 "trichoptera", "diptera",  "amphipoda", "coleoptera", 
                 "bivalvia")
family.order = c("viviparidae", "lymnaeidae","planorbidae",
                 "sphaeriidae",  "polycentropodidae", "chironomidae",
                 "macromiidae", "aeshnidae","libellulidae", "gomphidae",
                 "corduliidae" , "coenagrionidae", "leptophlebiidae",
                 "heptageniidae", "ephemerellidae","talitridae", "elmidae")

## Load in isotope data

data.iso = read.csv("Data/CSVs/processed_data.csv")
baselines = read.csv("Data/CSVs/baselines.csv")

## Taxon frame 

taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Cluster and chemistry data
load("Data/RData/cluster.mat.RData")
load(file = "Data/RData/community.legend.RData")
cluster_chem = read.csv(file = "Data/CSVs/chemistry_seasonal.csv") %>%
  unite("community.name", c(WATER, season), sep = ".") %>%
  #unite("ID", WATER, season, sep = ".") %>%
  left_join(community.legend) %>%
  left_join(cluster.mat %>%
              unite(ID, c("WATER", "season"), sep = "."))


## Mixtures for SIMMR model (invertebrate observations)
mixtures =  data.iso %>% 
  left_join(taxon_frame %>% select(ORDER, FAMILY, TAXON) %>% unique()) %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  mutate(GROUP = str_replace(GROUP, "CLAM","BIVALVE"),
         GROUP = str_replace(GROUP,"MUSSEL", "BIVALVE"))  %>%
  mutate(graph_id = case_when(is.na(ORDER) == T ~ GROUP, 
                              is.na(ORDER) == F ~ ORDER)) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall")) %>%
  unite("community.name", c(WATER, season), sep = ".") %>%
  mutate(group.name = case_when(is.na(FAMILY)== T ~ GROUP, 
                                is.na(FAMILY)==F ~ FAMILY)) %>%
  ungroup() %>%
  arrange(community.name) %>%
  mutate(community = as.numeric(as.factor(community.name))) %>%
  arrange(group.name) %>%
  mutate(group = as.numeric(as.factor(group.name))) %>%
  arrange(community, group) %>%
  as.data.frame() %>%
  filter(group.name %nin% c("LEAF", "PERI", "ZOOP", "FISH")) 




## For loop through communities 
community.legend = mixtures %>% 
  select(community.name, community) %>%
  unique()


## SIMMR ------------------
output.list = list() ## Fills w/ output of loop
for(i in 1:length(unique(mixtures$community))){ 

    subset.mix = mixtures %>% filter(community == i) %>%
    select(group.name,D13C, D15N) %>%
    na.omit() %>% ungroup() 
  
  subset.base = baselines %>% filter(community ==i) %>% ungroup() %>%
    na.omit()
  
  if(subset.base$community %>% length() > 0){
    ## without groups 
    simmr_in <- simmr_load(
      mixtures =  as.matrix(subset.mix[,2:3]),
      source_names = subset.base$group.name ,
      source_means = subset.base %>% select(mean_c, mean_n),
      source_sds = subset.base %>% select(sd_c, sd_n),
      group = subset.mix$group.name
    )
    
    simmr_out = simmr_mcmc(simmr_in)
    
    taxa.list = list()
 
    for(h in 1:length(unique(subset.mix$group.name))){
    
      taxa.list[[h]] = simmr_out$output[[unique(subset.mix$group.name)[h]]]$BUGSoutput$summary %>%
        as.data.frame() %>%
        mutate(taxa = unique(subset.mix$group.name)[h]) %>%
        rownames_to_column(var = "rowname")
      
      
    }
    taxa.output = Reduce(full_join, taxa.list)
    
    
    output.list[[i]] = taxa.output %>%
      mutate(community = (community.legend %>% filter(community == i))$community.name)
 
    
  
  }else{
    output.list[[i]] = NA
  }
}



## Turn 3D list into 2D data frame
simmr.full = Reduce(full_join, output.list) %>% 
  left_join(taxon_frame %>% 
              select(ORDER, FAMILY) %>% 
              filter(FAMILY != "unidentified") %>%
              unique(), by = c("taxa" = "FAMILY")) 

#save(simmr.full, file = "Data/RData/simmr_full.RData")

load(file = "Data/RData/simmr_full.RData") ## Also used in Step4_TrophicPosition.R


simmr.summary = simmr.full %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF")) %>%
  group_by(community) %>%
  mutate(groups = length(unique(rowname))) %>%
  select(groups, everything()) %>%
  filter(groups == 3) %>%
  ungroup() %>%
  group_by(community, rowname) %>%
  summarize(mean_peri = mean(mean) * 100) %>%
  rename("ID" = "community") %>%
  left_join(cluster_chem, by = "ID") %>%
  as.data.frame() %>%
  select(community, rowname, cluster, Lake, community.name, everything()) %>%
  mutate(across(.cols = TDO5:Volume, scale))

simmr.summary$Lake %>% unique() %>% length() ## 9 unique lakes used in this analysis


## Both periphyton and zooplankton show significant differences while leaves do not
lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + max_depth  +(1|Lake)) %>% AIC() 

lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth  +(1|Lake)) %>% summary()

summary()



## Looping for SIMMER LMER

vars = colnames(simmr.summary)[8:18]
vars
for(i in 1:length(vars)){
  for(h in 1:3){
    loop.simmr = simmr.summary %>% filter(rowname == unique(simmr.summary$rowname)[h])
    form = as.formula(paste0("mean_peri ~ ", vars[i], " + (1|Lake)"))
    lmer.run = lmerTest::lmer(form, data = loop.simmr) %>%
    summary()
    if(lmer.run$coefficients[2,5] < .05){
      print(vars[i])
      print(unique(simmr.summary$rowname)[h])
      print(lmer.run$coefficients[2,])
    }
  }
}


## Single predictor models --------------------------------------

simmr.summary = simmr.summary %>% mutate(sechi.depth = as.numeric(sechi.depth), 
                                         max_depth = as.numeric(max_depth), 
                                         mid_depth = as.numeric(mid_depth))

#### Important variables ---------------------
#DDO5: zoop
#Secchi depth: Peri + zoop
# mid depth: Leaf + Peri
# bottom depth: leaf + peri
# Max depth: zoop
  
## Testing multiple variables in an LMER
## Sticking with the simple regression above
library(lmerTest)
lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + (1|Lake)) %>% AIC() 

lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + max_depth  + (1|Lake)) %>% AIC()

lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + DDO5  + (1|Lake)) %>% AIC() 

lmerTest::lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + max_depth + DDO5  + (1|Lake)) %>% summary() 
simmr.summary %>%
  filter(rowname == "ZOOP") %>%
  select(sechi.depth, max_depth, DDO5) %>%
  cor()
### Random effect doesn't explain much additional variation. Leave it a simple LM
## Final models ----------------------
## Final model for zooplankton
zoop.final.simmr = lm(data = simmr.summary %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + max_depth  ) 

zoop.final.simmr %>% summary()
(simmr.summary %>% filter(rowname == "ZOOP"))$sechi.depth %>% class()
  
  
## Final model for periphyton
peri.final.simmr = lm(data = simmr.summary %>% filter(rowname == "PERI"), 
     mean_peri ~ sechi.depth + mid_depth)

peri.final.simmr %>% summary()

## Final model for leaves

terr.final.simmr = lm(data = simmr.summary %>% filter(rowname == "LEAF"), 
     mean_peri ~ mid_depth) 

### Table S5B -----------

### Summary table (done manually not in loop above for only significant variables) 
zoop.sum =  (zoop.final.simmr %>% summary)$coefficients %>% as.data.frame() %>%
  mutate(group = "zoop",
         R2 = (zoop.final.simmr %>% summary)$r.squared) 
peri.sum =  (peri.final.simmr %>% summary)$coefficients %>% as.data.frame() %>%
  mutate(group = "peri",
         R2 = (peri.final.simmr %>% summary)$r.squared) 
leaf.sum =  (terr.final.simmr %>% summary)$coefficients %>% as.data.frame() %>%
  mutate(group = "leaf",
         R2 = (terr.final.simmr %>% summary)$r.squared) 

summary_table = rbind(zoop.sum, peri.sum, leaf.sum) %>%
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  rownames_to_column(var = "rowname") 

#write.csv(summary_table, "Graphics/Tables/Supplemental_Table_5B.csv", row.names = F)



## Visualizations ----------------------

##### Figure 6: SIMMR x Secchi -------------

## Mean contribution of source averaged between all taxa that show up in lake

## Plotting


simmr.summary%>%
  ggplot(aes(x = sechi.depth, y = mean_peri, col = rowname)) + 
  geom_point(size = 3) + 
  geom_smooth(method = lm) +
  theme_minimal(base_size = 13) + 
  scale_color_manual("Source", labels = c("Terrestrial", "Nearshore", "Offshore"), values = wes_palette("Darjeeling1")) + 

  ylab("Average Source Contribution") + 
  xlab("Secchi Depth (m)")





### Plotting needs raw values (but same models)
plot.frame = simmr.full %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF")) %>%
  group_by(community) %>%
  mutate(groups = length(unique(rowname))) %>%
  select(groups, everything()) %>%
  filter(groups == 3) %>%
  ungroup() %>%
  group_by(community, rowname) %>%
  summarize(mean_peri = mean(mean) * 100) %>%
  rename("ID" = "community") %>%
  left_join(cluster_chem, by = "ID") %>%
  as.data.frame() %>%
  select(community, rowname, cluster, Lake, community.name, everything())  %>%
  mutate(sechi.depth = as.numeric(sechi.depth), 
         max_depth = as.numeric(max_depth), 
         mid_depth = as.numeric(mid_depth))


## zooplankton
zoop.plot = lm(data = plot.frame %>% filter(rowname == "ZOOP"), 
     mean_peri ~ sechi.depth + max_depth  ) 
zoop.line = ggpredict(zoop.plot, terms = "sechi.depth") %>% 
  as.data.frame() %>% 
  mutate(rowname = "ZOOP")

## periphyton
peri.plot = lm(data = plot.frame %>% filter(rowname == "PERI"), 
     mean_peri ~ sechi.depth + mid_depth)
peri.line = ggpredict(peri.plot, terms = "sechi.depth") %>% 
  as.data.frame() %>% 
  mutate(rowname = "PERI")

## leaves

terr.plot = lm(data = plot.frame %>% filter(rowname == "LEAF"), 
     mean_peri ~ mid_depth)
terr.line = ggpredict(terr.plot, terms = "mid_depth") %>% 
  as.data.frame() %>% 
  mutate(rowname = "LEAF")







A = plot.frame %>%
  filter(rowname != "LEAF") %>%
  ggplot(aes( col = rowname)) + 
  geom_point(aes(x = sechi.depth, y = mean_peri), size = 3) + 
  ### Zooplankton
  geom_ribbon(data = zoop.line,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.1,
    colour = NA ) +
  geom_line(data = zoop.line, mapping = (aes(x = x, y = predicted, col = rowname))) + 
  ### Periphyton
  geom_ribbon(data = peri.line,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high, ),
    alpha = 0.1,
    colour = NA ) +
    theme_minimal(base_size = 12) + 
  geom_line(data = peri.line, mapping = (aes(x = x, y = predicted, col = rowname))) + 
  scale_color_manual("Source", labels = c("Offshore", "Nearshore"), values = wes_palette("Darjeeling1")) + 
  ylab("Source Contribution (%)") +  
  xlab("Secchi Depth (m)")+ 
  theme(legend.position = "top")+ 
  ylim(10, 60)
  

### Leaves
B = plot.frame %>%
  filter(rowname == "LEAF") %>%
  ggplot(aes( col = rowname)) + 
  geom_point(aes(x = mid_depth, y = mean_peri), size = 3) + 
  geom_ribbon(data = terr.line,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.1,
    colour = NA ) +
  geom_line(data = terr.line, mapping = (aes(x = x, y = predicted, col = rowname))) +
  theme_minimal(base_size = 12) + 
  scale_color_manual("Source", labels = c("Terrestrial"), values = wes_palette("Darjeeling1")[3]) + 
  ylab("Source Contribution (%)") + 
  xlab("Mid-thermocine Depth (m)") + 
  theme(legend.position = "top") + 
  ylim(10, 60)


A + B
ggsave(width = 5, height = 4, filename = "Graphics/Figures/Figure_5.pdf")

##### Table S5: SIMMR % --------------


simmr.full %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF"),
         taxa %nin% c("INSECT"),
         taxa %in% family.order) %>%
  group_by(taxa) %>%
  mutate(count = n()) %>%
  select(count, everything()) %>%
  filter(count > 3) %>%
  ggplot(aes( x =factor(taxa, level = rev(family.order)),
             y = mean,
             fill =factor(ORDER, levels = order_levels))) +
  geom_boxplot() +
  coord_flip() +
  ylab("Source Contribution") +
  theme_minimal(base_size = 14) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual("Order", values = order_colors) +
  facet_wrap(~rowname, labeller = labeller(rowname = c("LEAF" = "Terrestrial",
                                                       "PERI" = "Nearshore", 
                                                       "ZOOP" = "Offshore"))) +
  theme(panel.spacing = unit(2, "lines"))

simmr.summary %>% pivot_wider(names_from  = "rowname", values_from =  "mean_peri") %>%
  select(community, PERI, LEAF, ZOOP)


