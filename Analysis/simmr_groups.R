## simmr with groups
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)

set.seed(123)

load(file = "Data/RData/mixtures.RData")
load(file = "Data/RData/baselines.RData")
load(file = "Data/RData/taxon_frame.RData")



taxon_frame = taxon_frame %>%
  rename(ORDER = ORDER..OR.ABOVE.)


load(file = "Data/RData/richness_cluster.RData")
richness = richness %>% unite("community.name", c(WATER, season), sep = ".")

mixtures %>% filter(GROUP == "BIVALVE")


#### ------------------

## Checking something



## With groups

plot(simmr_in)


simmr_out = simmr_mcmc(simmr_in)
summary(simmr_out, type = "diagnostics")
post_pred <- posterior_predictive(simmr_out)
plot(simmr_out)


mixtures %>% filter(FAMILY == "unionidae")

## For loop through communities 
community.legend = mixtures %>% 
  select(community.name, community) %>%
  unique()

community.legend

simmr.full %>% filter(taxa == "unionidae")


mixtures %>% filter(FAMILY == "unionidae") %>%
  select(community.name, community)

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


output.list.clean = output.list[-which(is.na(output.list))]

simmr.full = Reduce(full_join, output.list.clean)


simmr.full %>% filter(taxa == "unionidae")

#save(simmr.full, file = "Data/RData/simmr_full.RData")

load(file = "Data/RData/simmr_full.RData")


family.order = c("viviparidae", "lymnaeidae","planorbidae","sphaeriidae", "unionidae", "polycentropodidae", "chironomidae",
                 "macromiidae", "aeshnidae","libellulidae", "gomphidae", "corduliidae" , "coenagrionidae", "leptophlebiidae","heptageniidae", "ephemerellidae","talitridae", "elmidae")




simmr.full = simmr.full %>% 
  left_join(taxon_frame %>% 
              select(ORDER, FAMILY) %>% 
              filter(FAMILY != "unidentified") %>%
              unique(), by = c("taxa" = "FAMILY")) 




order_colors = wes_palette("Darjeeling1", type = "continuous", n = 8)
order_levels = c("odonata", "ephemeroptera", "hemiptera", "gastropoda", "trichoptera", "diptera",  "amphipoda", "coleoptera")

## Contribution of zooplankton to taxa

simmr.full %>%
  filter(rowname == "ZOOP",
         taxa %nin% c("INSECT", "SNAIL", "unidentified"), 
         taxa %in% family.order) %>%
  select(taxa, ORDER, mean, community) %>%
  unique() %>%
  group_by(taxa) %>%
  mutate(count = n()) %>%
  select(count, everything()) %>%
  filter(count > 3) %>%
  ggplot(aes(x =factor(taxa, level = rev(family.order)),
             y = mean, 
             fill =factor(ORDER, levels = order_levels))) + 
  geom_boxplot() +
  geom_point() +
  coord_flip() +
  ylab("Mean Contribution of Zooplankton") +
  theme_minimal(base_size = 20) +
  theme(axis.title.y = element_blank())  +
  scale_fill_manual("Order", values = order_colors) 


## Contribution of leaves to taxa
simmr.full %>%
  filter(rowname == "LEAF",
         taxa %in% family.order) %>%
  select(taxa, ORDER, mean, community) %>%
  unique() %>%
  group_by(taxa) %>%
  mutate(count = n()) %>%
  select(count, everything()) %>%
  filter(count > 3) %>%
  ggplot(aes(x =factor(taxa, level = rev(family.order)),
             y = mean,
             fill =factor(ORDER, levels = order_levels))) + 
  geom_boxplot() +
  coord_flip() +
  ylab("Mean Contribution of Leaves") +
  theme_minimal() +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual("Order", values = order_colors) 

## Contribution of periphyton to taxa
simmr.full %>%
  filter(rowname == "PERI",
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
  geom_point() + 
  coord_flip() +
  ylab("Mean Contribution of Periphyton") +
  theme_minimal(base_size = 14) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual("Order", values = order_colors) 









### Combining richness into the SIMMR full results



taxa.simmr = simmr.full %>%
  left_join(richness, by = c("community" = "community.name")) %>%
  filter(rowname == "PERI", 
         taxa %in% family.order) %>%
  #select(-GENUS, -TAXON) %>%
  unique()

taxa.simmr %>%
  ggplot(aes(x = sechi.depth, y = mean)) + 
  geom_jitter(aes(col = cluster), size = 2.5) +
  theme_minimal() +
  geom_smooth(method = "lm", col = "black", se = F) +
  ylab("Periphyton Contribution") +
  xlab("Sechi Depth") +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) 

  
lm(data = taxa.simmr, mean ~ thermo_depth) %>% summary()
lm(data = taxa.simmr  , mean ~ sechi.depth) %>% summary()
lm(data = taxa.simmr, mean ~ sechi.depth +thermo_depth ) %>% summary()

taxa.simmr$community %>% unique()

taxa.simmr %>% 
  ggplot(aes(x = thermo_depth, y = sechi.depth)) + geom_point() +
  geom_smooth( method = lm)


cor( taxa.simmr$thermo_depth, taxa.simmr$sechi.depth)




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
  left_join(richness) %>%
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
  left_join(richness) %>%
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
  left_join(richness)


## Trophic position across metrics - no real strong patterns between trophic position and the variables. See below for where it varies across clusters...
tpc %>%
  group_by(FAMILY, community.name, cluster, DOC.1, sechi.depth, temp.5mgL, FUI.num, depth.5mgL) %>%
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
  geom_point() + 
  geom


tpc.max = tpc %>%
  filter(GROUP != "CRAY", ) %>%
  group_by(community.name) %>%
  slice_max(trophic_position_weighted) 
cor.test(tpc.max$trophic_position_weighted, tpc.max$depth.5mgL, method = "spearman")
## Trophic position across clusters

tpc %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 5) %>% 
  ggplot(aes(y = trophic_position_weighted, x =temp.5mgL)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_minimal() +
  facet_wrap(~FAMILY, scales = "free_y") 

## Cluster 1 lakes have higher trophic positions than cluster 2
tpc %>% 
  group_by(FAMILY, community.name) %>%
  mutate(family_TP = mean(trophic_position_weighted)) %>%
  select(FAMILY, community.name, family_TP, cluster) %>%
  unique() %>%
  ungroup() %>%
  #filter(cluster %in% c(1,2)) %>%
  #glm(data = ., .$family_TP ~ .$cluster ) %>% summary()
  #aov(data = ., .$family_TP ~ .$cluster) %>% summary()
  #TukeyHSD(data = ., aov(data = ., .$family_TP ~ .$cluster)) ## # vs. 1 is significantly different
  ggplot(aes(x = cluster, y= family_TP, fill = cluster)) +
  theme_minimal() + 
  geom_boxplot() +
  geom_point() +
  ylab("Family Trophic Level") + 
  xlab("Cluster") +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  theme(legend.position = "none") 



## Now doing the same thing above but scaled

tpc %>% 
  group_by(FAMILY, community.name) %>%
  mutate(family_TP = mean(trophic_position_weighted)) %>%
  select(FAMILY, community.name, family_TP, cluster) %>%
  unique() %>%
  ungroup() %>%
  group_by(FAMILY) %>%
  mutate(family_TP = scale(family_TP)) %>%
  ungroup() %>%
  #filter(cluster %in% c(1,2)) %>%
  #glm(data = ., .$family_TP ~ .$cluster ) %>% summary()
  #aov(data = ., .$family_TP ~ .$cluster) %>% summary()
  #TukeyHSD(data = ., aov(data = ., .$family_TP ~ .$cluster)) ## # vs. 1 is significantly different
  ggplot(aes(x = cluster, y= family_TP, fill = cluster)) +
  theme_minimal() + 
  geom_boxplot() +
  geom_point() +
  ylab("Family Trophic Level") + 
  xlab("Cluster") +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  theme(legend.position = "none") 




## Testing if there are differences within taxa between clusters

tpc %>% 
  filter(cluster != 3) %>%
  ungroup() %>%
  group_by( FAMILY) %>%
  mutate(count = length(unique(cluster))) %>%
  select(count, everything()) %>%
  filter(count > 1) %>%

  select(FAMILY, community.name, trophic_position_weighted, cluster) %>%
  unique() %>%
  ungroup() %>%
  group_by(FAMILY) %>%
  mutate(trophic_position_weighted = scale(trophic_position_weighted)) %>%
  ggplot(aes(x = cluster, y = trophic_position_weighted)) + 
  geom_boxplot() + 
  geom_point() #+ 
  facet_wrap(~cluster)
  
  
## Functional feeding group by trophic position 
  
tpc %>% left_join(na.omit(ffg), by = "FAMILY") %>%
  group_by(FAMILY, community.name, FG, cluster) %>%
  summarize(mean_TP = mean(trophic_position_weighted)) %>%
  ggplot(aes(x = cluster, y = mean_TP)) + 
  geom_boxplot() + facet_wrap(~FG)

tpc %>% left_join(na.omit(ffg), by = "FAMILY") %>%
  filter(FG == "predator") %>%
  group_by(FAMILY, community.name, FG, cluster) %>%
  summarize(mean_TP = mean(trophic_position_weighted)) %>%
  #aov(data = ., .$mean_TP ~ .$cluster) %>% summary() 
  TukeyHSD(data = ., aov(data = ., .$mean_TP ~ .$cluster)) ## Again 1 vs. 3 is significant




## How high were the lakes

richness %>% select(community.name,Elevation) %>%
  arrange(Elevation)

richness %>% select(Elevation) %>% unique() %>%
  summarize(mean = mean(Elevation),
            min = min(Elevation),
            max = max(Elevation)) 

## Functional feeding groups across clusters

tpc %>% left_join(na.omit(ffg), by = "FAMILY") %>%
  select(FAMILY, FG, cluster) %>%
  unique() %>%
  na.omit() %>%
  ggplot(aes(x = cluster, fill = FG)) + 
  geom_bar()


tpc %>% left_join(na.omit(ffg), by = "FAMILY") %>%
  select(FAMILY, FG, cluster, community.name) %>%
  unique() %>%
  na.omit() %>%
  group_by(cluster, FG, community.name) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = cluster, y = count)) + geom_boxplot() + 
  facet_wrap(~FG, scales = "free")

presence_absence_data <- tpc %>%
  left_join(na.omit(ffg), by = "FAMILY") %>%
  select(FAMILY, FG, cluster) %>%
  unique() %>%
  na.omit() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = FG, values_from = present, values_fill = 0) %>%
  column_to_rownames(var = "FAMILY")

glm_results <- list()

for (fg in unique(presence_absence_data$FG)) {
  
  model <- glm(as.formula(paste(fg, "~ cluster")), 
               data = presence_absence_data, 
               family = binomial)
  glm_results[[fg]] <- summary(model)
}

##taxa %in% c("talitridae"),
simmr.full %>% filter( rowname == "PERI") %>%
  left_join(richness, by = c("community" = "community.name" )) %>%
  #lm(data = ., mean ~ DOC.1) %>% summary()

  ggplot(aes(y = mean,x = DOC.1 )) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~taxa, scales = "free_y")

