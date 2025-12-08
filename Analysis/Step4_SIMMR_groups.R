## simmr with groups
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)

set.seed(123)


## Load in isotope data

data.iso = read.csv("Data/CSVs/processed_data.csv")
baselines = read.csv("Data/CSVs/baselines.csv")

## Taxon frame 

taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Cluster and chemistry data
cluster_chem = read.csv(file = "Data/CSVs/richness_update.csv") 


# visualize the resources
baselines %>% separate(community.name, into = c("water", "season")) %>%
  arrange(water) 

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

mixtures$group.name %>% unique()

#### ------------------

## Checking something



## For loop through communities 
community.legend = mixtures %>% 
  select(community.name, community) %>%
  unique()

community.legend




output.list = list()


for(i in 1:length(unique(mixtures$community))){
  

  
  community.name =  (mixtures %>% 
    filter(community == i) %>% 
    select(community.name) %>% 
    unique())$community.name
  
  test.mix = mixtures %>% filter(community == i) %>%
    select(group.name,D13C, D15N) %>%
    na.omit() %>% ungroup() 
  
  test.base = baselines %>% filter(community ==i) %>% ungroup() %>%
    na.omit()
  
  if(test.base$community %>% length() > 0){
    ## without groups 
    simmr_in <- simmr_load(
      mixtures =  as.matrix(test.mix[,2:3]),
      source_names = test.base$group.name ,
      source_means = test.base %>% select(mean_c, mean_n),
      source_sds = test.base %>% select(sd_c, sd_n),
      group = test.mix$group.name
    )
    
    simmr_out = simmr_mcmc(simmr_in)
    
    taxa.list = list()
 
    for(h in 1:length(unique(test.mix$group.name))){
    
      taxa.list[[h]] = simmr_out$output[[unique(test.mix$group.name)[h]]]$BUGSoutput$summary %>%
        as.data.frame() %>%
        mutate(taxa = unique(test.mix$group.name)[h]) %>%
        rownames_to_column(var = "rowname")
      
      
    }
    
    taxa.output = Reduce(full_join, taxa.list)
    
    
    output.list[[i]] = taxa.output %>%
      mutate(community = community.name)
 
    
  
  }else{
    output.list[[i]] = NA
  }
}



simmr.full = Reduce(full_join, output.list)


simmr.full$taxa %>% unique()




simmr.full %>% filter(taxa == "unionidae")

#save(simmr.full, file = "Data/RData/simmr_full.RData")

load(file = "Data/RData/simmr_full.RData")







simmr.full = simmr.full %>% 
  left_join(taxon_frame %>% 
              select(ORDER, FAMILY) %>% 
              filter(FAMILY != "unidentified") %>%
              unique(), by = c("taxa" = "FAMILY")) 




order_colors = wes_palette("Darjeeling1", type = "continuous", n = 8)

order_levels = c("odonata", "ephemeroptera", "hemiptera", "gastropoda", "trichoptera", "diptera",  "amphipoda", "coleoptera", "bivalvia")

family.order = c("viviparidae", "lymnaeidae","planorbidae","sphaeriidae", "unionidae", "polycentropodidae", "chironomidae",
                 "macromiidae", "aeshnidae","libellulidae", "gomphidae", "corduliidae" , "coenagrionidae", "leptophlebiidae","heptageniidae", "ephemerellidae","talitridae", "elmidae")

## Contribution of sources to taxa


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
  #geom_point() + 
  coord_flip() +
  ylab("Source Contribution") +
  theme_minimal(base_size = 14) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual("Order", values = order_colors) +
  facet_wrap(~rowname, labeller = labeller(rowname = c("LEAF" = "Leaf", "PERI" = "Periphyton", "ZOOP" = "Zooplankton")))+ theme(panel.spacing = unit(2, "lines"))

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
lmer(data = simmr.summary %>% filter(rowname == "PERI"), mean_peri ~ sechi.depth + (1|Lake)) %>% summary()




## trophic position -- trying to mix TP from different source contributions

baseline.mixtures = simmr.full %>% 
  select(rowname, mean, taxa, community) %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF")) %>%
  unique() %>%
  pivot_wider(names_from = rowname, values_from = mean) %>%
  rename(PERI.mix = PERI, 
         ZOOP.mix = ZOOP, 
         LEAF.mix = LEAF)

baselines.wide = baselines %>%
  ungroup() %>%
  select(community.name, group.name, mean_n) %>%
  pivot_wider(names_from = group.name, values_from = mean_n)


trophic.position = mixtures %>%
  select(community.name, GROUP, FAMILY, D15N, D13C) %>%
  left_join(baseline.mixtures, by = c( "community.name" ="community", "FAMILY" = "taxa"  )) %>%
  left_join(baselines.wide, by = c("community.name")) %>%
  mutate(baseline_N = case_when(is.na(PERI.mix) == F ~ ZOOP * ZOOP.mix + PERI * PERI.mix + LEAF * LEAF.mix,
                                is.na(PERI.mix) == T ~  ZOOP * ZOOP.mix +LEAF * LEAF.mix)) %>%
  mutate(baseline_position = case_when(is.na(PERI.mix) == F ~ 2 * ZOOP.mix + 1 * PERI.mix + 1 * LEAF.mix,
                                        is.na(PERI.mix) == T ~ 2 * ZOOP.mix + 1 * LEAF.mix)) %>%
  mutate(trophic_position_weighted = ((D15N - baseline_N)/3.14) + baseline_position) %>%
  mutate(trophic_position = case_when(is.na(PERI.mix) == F ~ ((D15N - PERI)/3.14) + 1,
                                      is.na(PERI.mix) == T ~ ((D15N - LEAF) / 3.14) + 1)) %>%
  select(community.name, GROUP, FAMILY, D15N, D13C, trophic_position_weighted, trophic_position) %>%
  na.omit() %>%
  filter(FAMILY %nin% c("unidentified", "terrestrial", "degraded")) %>%
  left_join(taxon_frame %>% select(ORDER, FAMILY) %>% unique()) %>%
  ungroup() %>%
  group_by(FAMILY) %>%
  mutate(mean_TP = mean(trophic_position_weighted))

save(trophic.position, file = "Data/RData/trophic.position.RData")
load(file = "Data/RData/trophic.position.RData")

TP.family = trophic.position %>% 
  ungroup() %>%
  select(ORDER, FAMILY) %>%
  arrange(ORDER) %>% select(FAMILY) %>% unique()
TP.order = trophic.position %>% 
  ungroup() %>%
  select(ORDER, FAMILY) %>%
  arrange(ORDER) %>% select(ORDER) %>% unique()


## Plot the trophic position of the taxa using the weighted estimation from SIMMR
trophic.position %>%
  left_join(cluster_chem) %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 5) %>%
  ggplot(aes(y = factor(FAMILY, rev(TP.family$FAMILY)),
             x = trophic_position_weighted,
             fill = factor(ORDER, TP.order$ORDER))) + 
  geom_boxplot() +
  geom_point(alpha = .2) + 
  theme_minimal(base_size = 20) +
  scale_fill_manual("Order", values = wes_palette("Darjeeling1", type = "continuous", n = 9)) +
  xlab("Trophic Position") +
  theme(axis.title.y = element_blank()) +
  scale_x_log10()




## Comparing the calculation of trophic position

trophic.position %>%
 # left_join(richness) %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 5) %>%
  pivot_longer(c(trophic_position_weighted, trophic_position),
               values_to = "TP", names_to = "type") %>%
  ggplot(aes(x = TP, y = FAMILY, fill = type)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_fill_manual("TP Calculation", values = wes_palette("Darjeeling1", n = 2))



### Comparing trophic positions across variables
tpc = trophic.position %>%
  left_join(cluster_chem)


## Trophic position across metrics - no real strong patterns between trophic position and the variables. See below for where it varies across clusters...
tpc %>%
  group_by(FAMILY, community.name, cluster,  sechi.depth, temp.5mgL,  depth.5mgL) %>%
  summarize(mean_TP = mean(trophic_position_weighted)) %>%
  #lm(data = ., mean_TP ~ .$depth.5mgL) %>% summary()
  ggplot(aes(x = depth.5mgL, y = mean_TP)) + 
  geom_point() +
  geom_smooth(method = "lm")

## What about max trophic position - aka longer food chains within different lakes?

tpc %>%
  filter(GROUP != "CRAY") %>%
  group_by(community.name) %>%
  slice_max(trophic_position_weighted) %>%
  ggplot(aes(y = trophic_position_weighted, x =depth.5mgL)) + 
  geom_point() 


tpc.max = tpc %>%
  filter(GROUP != "CRAY", ) %>%
  group_by(community.name) %>%
  slice_max(trophic_position_weighted) 
cor.test(tpc.max$trophic_position_weighted, tpc.max$depth.5mgL, method = "spearman")
## Trophic position across clusters


tpc %>% 
  filter(FAMILY %in% c("aeshnidae", "gomphidae", "libellulidae", "corduliidae")) %>%
  ggplot(aes(x = as.factor(cluster), y = trophic_position_weighted, fill = as.factor(cluster))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_minimal(base_size = 15) + 
  ylab("Trophic Position")  + 
  scale_x_discrete(labels = c(
  "1" = "Shallow\nthermocline", 
  "3" = "Deep\nthermocline", 
  "2" = "Diverse,\ndeep thermocline"
))  +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
   scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) 



tpc %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 5) %>% 
  ggplot(aes(y = trophic_position_weighted, x =temp.5mgL)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_minimal() +
  facet_wrap(~FAMILY, scales = "free_y") 

# Note - no real difference in trophic position between clusters (scaled or unscaled) 12/8/25



## How high were the lakes

richness %>% select(community.name,Elevation) %>%
  arrange(Elevation)

richness %>% select(Elevation) %>% unique() %>%
  summarize(mean = mean(Elevation),
            min = min(Elevation),
            max = max(Elevation)) 

