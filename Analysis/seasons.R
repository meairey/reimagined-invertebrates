`%nin%` = Negate(`%in%`)
## Seasonal differences

load(file = "Data/RData/FullData.RData")
load(file = "Data/RData/richness_cluster.RData")
load(file = "Data/RData/nmdsDat.RData")

nmds.dat = nmds.dat %>% 
  rownames_to_column(var = "ID") %>%
  separate(ID, into = c("season", "WATER")) %>%
  mutate(season =  tolower(season)) %>%
  left_join(cluster.mat, by = c("WATER", "season"))

## SIMPER analysis based on season  -------------------------------
## Cluster 1 remove any ALC lakes
cluster1.simperdat = nmds.dat %>% filter(cluster == 1) %>%
  filter(WATER %nin% c("LML", "ETL")) 
  
simper_result.1 = simper(cluster1.simperdat %>% 
                         select(-cluster, -season, -WATER) %>% 
                         select_if(~ sum(.) != 0)%>%
                           select_if(~ sum(.) != 9), 
                         (cluster1.simperdat$season) %>% as.factor() )

simper1 = simper_result.1$fall_spring$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C1")

## Cluster 2 remove any ALC lakes
cluster2.simperdat = nmds.dat %>% filter(cluster == 2) 

simper_result.2 = simper(cluster2.simperdat %>%
                           select(-cluster, -season, -WATER) %>% 
                           select_if(~ sum(.) != 0)%>%
                           select_if(~ sum(.) != 9), 
                         (cluster2.simperdat$season) %>% as.factor() )

simper2 = simper_result.2$fall_spring$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C2")


## Cluster 3 remove any ALC lakes
cluster3.simperdat = nmds.dat %>% filter(cluster == 3) %>%
  filter(WATER %nin% c("LML", "ETL"))

simper_result.3 = simper(cluster3.simperdat %>%
                           select(-cluster, -season, -WATER) %>% 
                           select_if(~ sum(.) != 0) %>%
                           select_if(~ sum(.) != 9), 
                         (cluster3.simperdat$season) %>% as.factor() )

simper3 = simper_result.3$fall_spring$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C3")


simper.season = rbind(simper1, simper2, simper3)%>%
  rename(contribution = V1)

library(tidytext)
# Create a bar plot
simper.season %>%
  mutate(family = reorder_within(family, contribution, group))%>%
  ggplot(aes(x = family, y = contribution)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Contribution of Taxa to Dissimilarity",
       x = "Taxa",
       y = "Average Contribution") +
  theme_minimal() + 
  facet_wrap(~group, scales = "free") +
  scale_x_reordered()



## Richness by month
richness %>%
  filter(WATER %nin% c("LML", "ETL","COM","GNL")) %>%
  ggplot(aes(x = season, y = richness)) + 
  geom_boxplot() +
  geom_jitter() + 
  facet_wrap(~cluster) + 
  theme_minimal(base_size = 12) +
  xlab("Season") + ylab("Richness")



## LMER 
library(lme4)
# Example: Poisson GLMM for species richness
model <- glmer(richness ~ cluster + (1 | WATER) + (1 | season), 
               data = richness, 
               family = poisson)

summary(model)


full %>% left_join(richness %>% 
  mutate(MONTH = toupper(season))) 
  


## Common taxa across season and cluster


nmds.dat %>%
  ungroup() %>%
  select(season, cluster, everything()) %>%
  pivot_longer(4:49) %>%
  group_by(season, cluster, name) %>%
  mutate(tot = n()) %>%
  ungroup() %>%
  group_by(season, tot, cluster,name) %>%
  summarize(total = sum(value)) %>%
  mutate(proportion.group = total/tot) %>%
  filter(proportion.group == 1) %>%
  ungroup() %>%
  group_by(season, name) %>%
  mutate(clus.count = length(unique(cluster))) %>%
  select(clus.count) %>%
  ggplot(aes(y = name, x = season, col = clus.count)) +
  geom_point(size = 4) +
  theme_minimal() + 
  xlab("Season") + ylab("Family") + labs(col = "Total Clusters \n Where Taxa \n Is Observed")
  


## Rare taxa but I like the SIMPER representation better
nmds.dat %>%
  ungroup() %>%
  select(season, cluster, everything()) %>%
  pivot_longer(4:49) %>%
  group_by(season, cluster, name) %>%
  mutate(total_presence = sum(value), 
         tot = n()) %>%
  mutate(test = total_presence / tot)%>% 
  select(test, everything()) %>%
  filter(test < .3 & test > 0) %>%
  ggplot(aes(x = season, y = name, col = test)) + 
  geom_point() + 
  facet_wrap(~cluster, scales = "free")
  

## Taxa across clusters


nmds.dat %>% 
  ungroup() %>%
  select(season, cluster, everything()) %>%
  pivot_longer(4:49) %>%
  group_by(cluster) %>%
  mutate(total_lakes = length(unique(WATER))) %>%
  ungroup() %>%
  select(cluster, name, total_lakes, value, WATER) %>%
  unique() %>%
  group_by(cluster, name, total_lakes) %>%
  summarize(sum = sum(value)) %>%
  mutate(proportion.lakes = sum/total_lakes) %>%
  filter(proportion.lakes > .75) %>% ## Taxa that occur in more than 50% of lakes in that cluster
  ggplot(aes(x = cluster, y = name, col = proportion.lakes)) + 
  geom_point(size = 4) + 
  facet_wrap(~cluster, scales = "free") + 
  theme_minimal() +
  xlab("Cluster") + ylab("Family") + 
  labs(col = "Proportion of\nLakes")
