`%nin%` = Negate(`%in%`)
library(simmr)
library(tidyverse)
library(wesanderson)
library(lme4)
set.seed(123)
### SIMMR script to estimate trophic position using baseline mixtures


## Load in taxon data
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column



## Load in isotope data

data.iso = read.csv("Data/CSVs/processed_data.csv")
baselines = read.csv("Data/CSVs/baselines.csv")
  
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

## For loop through communities 


output.list = list()

for(i in 1:length(unique(mixtures$community))){
  
  
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
    source_sds = test.base %>% select(sd_c, sd_n)
  )
  
  simmr_out = simmr_mcmc(simmr_in)
   
  output.list[[i]] = as.data.frame(simmr_out$output$`1`$BUGSoutput$summary) %>%
    mutate(community.name = test.base$community.name%>% unique() ) %>%
    rownames_to_column(var = "rowname")
  
  
  
  plot(simmr_out, type = "boxplot", title = test.base$community.name%>% unique())
  }else{
    output.list[[i]] = NA
  }
}

simmr.output = Reduce(full_join, output.list) ## data frame of SIMMR output

#save(simmr.output, file = "Data/RData/simmr.output.RData")
#load(file = "Data/RData/simmr.output.RData")



## All groups together
simmr.output %>% 
  left_join(cluster_chem, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  ggplot(aes(x = sechi.depth, y = mean, col =  rowname)) + 
  geom_point() + 
  geom_smooth(method = "lm")

## All groups where all three sources are measured
filtered.simmr = simmr.output %>% 
  left_join(cluster_chem, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) 




filtered.simmr %>% 
  ggplot(aes(x = sechi.depth, y = mean, col = rowname)) +
  geom_point(size = 3, alpha = .5, .key_glyph = "rect") +
  geom_smooth(method = lm) + 
  theme_minimal(base_size = 20) +
  scale_color_manual("Source", values = wes_palette("Darjeeling1"), 
                     labels = c("Leaf", "Periphyton", "Zooplankton")) +
  xlab("Secchi Depth") +
  ylab("Source Contribution")


zoop.test = filtered.simmr %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3)  %>% 
  #left_join(cluster_chem, by = c("community.name")) %>%
  filter(rowname %in% "PERI") %>%
  separate(community.name, into = c("WATER", "SEASON")) %>% 
  filter(WATER != "UCL")

zoop.test %>% filter(sechi.depth > 5) %>%
  select(sechi.depth, mean, WATER,  SEASON, everything())

zoop.test %>%
  ggplot(aes(x = sechi.depth, y = mean)) + 
  geom_point() + 
  geom_smooth(method = "lm") 



lmer(data = zoop.test, mean ~ sechi.depth + (1 | WATER)) %>% summary()
lm(data = zoop.test, mean ~ sechi.depth) %>% summary()


## Significant increase of periphyton with increasing sechi depth



glm(data = zoop.test,mean ~ thermo_depth) %>% summary()

## filtering out lakes without all three sources and comparing across clusters
simmr.output %>% 
  left_join(cluster_chem %>%
    mutate(cluster = if_else(is.na(cluster), first(na.omit(cluster)), cluster)) , by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) %>%
  ggplot(aes(x = as.factor(cluster), y = mean, col = rowname)) +
  geom_boxplot() +
  geom_point()





simmr.all = simmr.output %>% 
  left_join(cluster_chem, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) 
 
  
simmr.all %>% ggplot(aes(x = rowname, y = mean, fill = rowname)) + 
  geom_boxplot() +
  geom_point() +
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1")) + 
  labs(fill = "Source") + 
  scale_x_discrete(labels = c("Leaf", "Periphyton", "Zooplankton")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Average Source Contribution")




## Food webs sampled are very much powered by periphyton compared to zooplankton Which makes sense because we sampled littoral webs
aov.simmr = simmr.all %>%
  aov(data = ., .$mean ~ .$rowname) 

print(aov.simmr %>% summary())

TukeyHSD(aov.simmr) 

### SIMMR broken up with the two source models separated but included in the graph

source_labels = c("2" = "Two Source Model", "3" = "Three Source Model")

simmr.output %>% 
  left_join(cluster_chem, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  mutate(rowname = str_replace(rowname, "ZOOP", "Zooplankton"), 
         rowname = str_replace(rowname, "PERI", "Periphyton"), 
         rowname = str_replace(rowname, "LEAF", "Leaf")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname)))  %>%  

  ggplot(aes(x = rowname, y = mean, fill = rowname)) + 
  geom_boxplot(width = .5) +
  geom_point() +
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1")) + 
  labs(fill = "Source") + 

  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Average Source Contribution") + 
  facet_wrap(~groups, scales = "free_x", labeller = as_labeller(source_labels))







