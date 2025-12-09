# SIMMR (taxon specific) ------------------------------

## Setup ----------------------------
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)
library(lme4)

set.seed(123)
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
cluster_chem = read.csv(file = "Data/CSVs/richness_update.csv") 

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

## Visualizations ----------------------

##### Figure 6: SIMMR x Secchi -------------

## Mean contribution of source averaged between all taxa that show up in lake

simmr.summary = simmr.full %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF")) %>%
  group_by(community) %>%
  mutate(groups = length(unique(rowname))) %>%
  select(groups, everything()) %>%
  filter(groups == 3) %>%
  ungroup() %>%
  group_by(community, rowname) %>%
  summarize(mean_peri = mean(mean)) %>%
  left_join(cluster_chem, by = c("community" = "community.name")) %>%
  as.data.frame()

simmr.summary$Lake %>% unique() %>% length() ## 9 unique lakes used in this analysis

simmr.summary%>%
  ggplot(aes(x = sechi.depth, y = mean_peri, col = rowname)) + 
  geom_point(size = 3) + 
  geom_smooth(method = lm) +
  theme_minimal(base_size = 13) + 
  scale_color_manual("Source", labels = c("Terrestrial", "Nearshore", "Offshore"), values = wes_palette("Darjeeling1")) + 

  ylab("Average Source Contribution") + 
  xlab("Secchi Depth (m)")

## Both periphyton and zooplankton show significant differences while leaves do not
lm(data = simmr.summary %>% filter(rowname == "ZOOP"), mean_peri ~ sechi.depth) %>% summary()
lmer(data = simmr.summary %>% filter(rowname == "ZOOP"), mean_peri ~ sechi.depth + (1|Lake)) %>% summary()


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



