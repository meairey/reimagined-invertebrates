`%nin%` = Negate(`%in%`)
library(simmr)
library(tidyverse)
### SIMMR script to estimate trophic position using baseline mixtures


# Taxon frame
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column
## Load in isotope measurement file
data.iso =  read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  mutate(D13C = as.numeric(D13C))
  mutate(D13C = case_when(CATEGORY %nin% c("ALGA", "PLANT") ~ as.numeric(D13C) - 3.32 + .99 * (PER_C / PER_N),
                          CATEGORY %in% c("ALGA", "PLANT") ~ as.numeric(D13C))) ## lipid correction factor

## Read in isotope sample file
sample = read.csv("Data/CSVs/isotope_sample.csv")

## Create the data frame for the baselines
baselines =  data.iso %>% left_join(sample, by = "ISO_YSAMP_N") %>%
  left_join(taxon_frame) %>%
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
  filter(group.name %in% c("LEAF", "PERI", "ZOOP")) %>%
  group_by(community, community.name, group.name) %>%
    summarize(mean_c = mean(D13C, na.rm = T), 
              mean_n = mean(D15N, na.rm = T), 
              sd_c = sd(D13C, na.rm = T), 
              sd_n = sd(D15N, na.rm = T)) %>%
  na.omit() %>%
  ungroup() %>%
  complete(community.name, group.name) %>%
  select(community, community.name, group.name, mean_c, mean_n, sd_c, sd_n)


## Adding in zooplankton samples that average across both spring and fall or use spring and fall data
## adding the averaged information for HTL zoop in bc not enough 
baselines[24, 1] = 8; baselines[24, 2] = "HTL.fall"; baselines[24,3] = "ZOOP"
baselines[24,4] = -33.94; baselines[24,5] = 6.676667; baselines[24,6] = 1.105908
baselines[24,7] = 2.206362


## adding the averaged information for SEL zoop in bc not enough 
baselines[51, 1] = 18; baselines[51, 2] = "SEL.fall"; baselines[51,3] = "ZOOP"
baselines[51,4] = -31.1; baselines[51,5] = 6.60; baselines[51,6] = 0.735
baselines[51,7] = 2.84




## saving the baseline frame
save(baselines, file = "Data/RData/baselines.RData")
load(file = "Data/RData/baselines.RData")

# visualize the resources
baselines %>% separate(community.name, into = c("water", "season")) %>%
  arrange(water)  %>%
  print(n = 100)

mixtures =  data.iso %>% left_join(sample, by = "ISO_YSAMP_N") %>%
  left_join(taxon_frame) %>%
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


save(mixtures, file = "Data/RData/mixtures.RData")
load(file = "Data/RData/mixtures.RData")
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



output.list.clean = output.list[-which(is.na(output.list))]

simmr.output = Reduce(full_join, output.list.clean)

save(simmr.output, file = "Data/RData/simmr.output.RData")
load(file = "Data/RData/simmr.output.RData")

simmr.output %>% filter(rowname %in% c("LEAF", "ZOOP", "PERI")) %>%
  select(rowname, mean, community.name) %>%
  ggplot(aes(x = community.name, y = mean, fill = rowname)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  coord_flip()

load("Data/RData/richness_cluster.RData")
richness = richness %>% unite("community.name", c(WATER, season), sep = ".")

## All groups together
simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  ggplot(aes(x = sechi.depth, y = mean, col =  rowname)) + 
  geom_point() + 
  geom_smooth(method = "lm")

## All groups where all three sources are measured
filtered.simmr = simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) 

filtered.simmr[, c(-3,-4,-5,-6,-7,-8,-9,-10)] %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4) %>%
  pivot_longer(c(Elevation:FUI.num), names_to = "metric", values_to = "values") %>%

  ggplot(aes(x = values, y = mean, col =  rowname)) + 
  theme_minimal() + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~metric, scales = "free_x")



filtered.simmr %>% 
  ggplot(aes(x = sechi.depth, y = mean, col = rowname)) +
  geom_point(size = 3, alpha = .5, .key_glyph = "rect") +
  geom_smooth(method = lm) + 
  theme_minimal(base_size = 20) +
  scale_color_manual("Source", values = wes_palette("Darjeeling1"), 
                     labels = c("Leaf", "Periphyton", "Zooplankton")) +
  xlab("Secchi Depth") +
  ylab("Source Contribution")


zoop.test = simmr.output %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3)  %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% "PERI")

zoop.test %>%
  ggplot(aes(x = sechi.depth, y = mean)) + 
  geom_point() + 
  geom_smooth(method = "lm")



lm(data = zoop.test, mean ~ sechi.depth) %>% summary()


## Significant increase of periphyton with increasing sechi depth



glm(data = zoop.test,mean ~ thermo_depth) %>% summary()

## filtering out lakes without all three sources and comparing across clusters
simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) %>%
  ggplot(aes(x = cluster, y = mean, col = rowname)) +
  geom_boxplot() +
  geom_point()





simmr.all = simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 2) 
 
  
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



simmr.subset = simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 2) 



simmr.subset %>%  ggplot(aes(x = rowname, y = mean, fill = rowname)) + 
  geom_boxplot() +
  geom_point() +
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1")) + 
  labs(fill = "Source") + 
 # scale_x_discrete(labels = c("Leaf", "Periphyton", "Zooplankton")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Average Source Contribution")


### SIMMR broken up with the two source models separated but included in the graph

source_labels = c("2" = "Two Source Model", "3" = "Three Source Model")

simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
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







