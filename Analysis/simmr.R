`%nin%` = Negate(`%in%`)
library(simmr)
library(tidyverse)
### SIMMR script to estimate trophic position using baseline mixtures


# Taxon frame
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column
## Load in isotope measurement file
data.iso=  read.csv("Data/CSVs/iso_data_clean.csv")

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
  na.omit()


data.iso %>% filter(GROUP == "LEAF") %>%
  filter(grepl( "UCL", ISO_YSAMP_N)) %>% summarize(mean = mean(D15N))
(data.iso %>% filter(GROUP == "LEAF") %>%
  filter(grepl( "UCL", ISO_YSAMP_N)) %>% select(D13C))$D13C %>% as.numeric() %>% mean()

%>% as.vector() %>% class()

## adding the averaged information for HTL zoop in bc not enough 
baselines[50, 1] = 8; baselines[50, 2] = "HTL.fall"; baselines[50,3] = "ZOOP"
baselines[50,4] = -33.94; baselines[50,5] = 6.676667; baselines[50,6] = 1.105908
baselines[50,7] = 2.206362





## Making ETL leaves, the same as GNL leaves
baselines[51,1] = 6; baselines[51,2] = "ETL.fall"; baselines[51,3] = "LEAF"
baselines[51,4] = -32.95; baselines[51,5] = -4.73
baselines[51,6]=0.5807179;  baselines[51,7] = 0.5980803

## Making LML Leaves, the same as GNL leaves
baselines[52,1] = 12; baselines[52,2] = "LML.fall"; baselines[52,3] = "LEAF"
baselines[52,4] = -32.95; baselines[52,5] = -4.73
baselines[52,6]=0.5807179;  baselines[52,7] = 0.5980803

## adding the averaged information for HTL

baselines[53, 1] = 9; baselines[53, 2] = "HTL.spring"; baselines[53,3] = "LEAF"
baselines[53,4] = -32.46833; baselines[53,5] =  -4.1383337; baselines[53,6] = 0.5194773
baselines[53,7] =2.279907


## adding the averaged information for UCL.spring leaves
baselines[54, 1] = 22; baselines[54, 2] = "UCL.spring"; baselines[54,3] = "LEAF"
baselines[54,4] = -31.7825; baselines[54,5] =-2.005 ; baselines[54,6] =1.259736
baselines[54,7] =1.467914


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
  filter(group.name %nin% c("LEAF", "PERI", "ZOOP")) 




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

zoop.test = simmr.output %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3)  %>% left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% "PERI")

zoop.test %>%
  ggplot(aes(x = FUI.num, y = mean)) + 
  geom_point() + 
  geom_smooth(method = "lm")

lm(data = zoop.test, mean ~ fui.num) %>% summary()



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




library(wesanderson)
simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) %>%
  ggplot(aes(x = rowname, y = mean, fill = rowname)) + 
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
simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  filter(groups == 3) %>%
  aov(data = ., .$mean ~ .$rowname) %>% summary()





















#### ------------------

## With groups
simmr_in <- simmr_load(
  mixtures =  as.matrix(test.mix[,2:3]),
  source_names = test.base$group.name ,
  source_means = test.base %>% select(mean_c, mean_n),
  source_sds = test.base %>% select(sd_c, sd_n),
  group = test.mix$group.name
)
plot(simmr_in)

simmr_out = simmr_mcmc(simmr_in)
summary(simmr_out, type = "diagnostics")
post_pred <- posterior_predictive(simmr_out)

compare_groups(simmr_out,
               source = "LEAF",
               groups = 1:2
)

