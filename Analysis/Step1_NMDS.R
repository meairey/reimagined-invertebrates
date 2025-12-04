set.seed(1234)  # For reproducibility
library(tidytext)
library(vegan)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(wesanderson)
library(ggrepel)
`%nin%` = Negate(`%in%`)

# Load
load(file = "Data/RData/FullData.RData")

## Full
## create the NMDS data frame ------------------

## Families included 
unique_families = full %>% select(FAMILY, GENUS) %>%  unique() 
#write.csv(file = "Data/CSVs/family.csv", unique_families, row.names = F)

## NMDS data
nmds.dat = full %>% 
  filter(MONTH =="FALL") %>% ## Filtering out the NMDS to just use fall information
  select(MONTH, YSAMP, SITE_N, WATER, ORDER..OR.ABOVE., FAMILY) %>%
  unique() %>%
  filter(FAMILY %nin% c("unidentified", "degraded","terrestrial")) %>%
  mutate(SITE_N = str_replace_all(SITE_N, " ", ""),
         SITE_N = str_replace_all(SITE_N, "Notlisted", "NotListed")) %>% 
  select(MONTH, WATER, FAMILY) %>%
  arrange(MONTH, WATER, FAMILY) %>%
  mutate(pres = 1) %>%
  unite("ID", 1:2) %>%
  unique() %>%
  pivot_wider(names_from = FAMILY, values_from = pres ) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "ID") 

#save(file = "Data/RData/nmdsDat.RData", nmds.dat)

order_matrix = full %>%
  filter(MONTH == "FALL") %>%
  select(MONTH, YSAMP, SITE_N, WATER, ORDER..OR.ABOVE., FAMILY) %>%
  unique() %>%
  filter(FAMILY %nin% c("unidentified", "degraded","terrestrial")) %>%
  mutate(SITE_N = str_replace_all(SITE_N, " ", ""),
         SITE_N = str_replace_all(SITE_N, "Notlisted", "NotListed")) %>% 
  select( ORDER..OR.ABOVE., FAMILY) %>%
  arrange(ORDER..OR.ABOVE., FAMILY) %>%
  rename(ORDER = ORDER..OR.ABOVE.) %>% unique()

# Run the NMDS using Bray-Curtis distance ------------------------------
nmds = metaMDS(nmds.dat, distance = "jaccard", trymax = 100) ## Do not run this again just load the data out of file below

#save(nmds, file = "Data/RData/nmds.RData") ## Do not run this again just load the data out of the file below



#load(file = "Data/RData/nmds.RData")
## K means clustering of lakes by community
# extract scores from NMDS
site_scores = as.data.frame(scores(nmds, display = "sites"))


# Run for loop and use elbow method to determine optimal clusters --------------
## test k 
# 
wss = numeric()

# Use a for loop to run k-means for k = 1 to 10
for (k in 1:10) {
  kmeans_result = kmeans(site_scores, centers = k, nstart = 10)
  wss[k] = kmeans_result$tot.withinss  # Store the total within-cluster sum of squares
}

# plot ------------
elbow_plot = data.frame(k = 1:10, WSS = wss)

ggplot(elbow_plot, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Elbow Method for Finding Optimal K",
       x = "Number of Clusters (K)",
       y = "Total Within-Cluster Sum of Squares") +
  theme_minimal()

## Optimal clusters looks like k = 3

k = 3  # Number of clusters you want to try
kmeans_result = kmeans(site_scores, centers = k)

# Add the cluster results to the site_scores data frame
site_scores$cluster = as.factor(kmeans_result$cluster)
#save(file = "Data/RData/site_scores.RData", site_scores)
#load(file = "Data/RData/site_scores.RData")
#load(file = "Data/RData/simmr.output.RData")
cluster.mat = site_scores %>% rownames_to_column(var = "ID") %>%
  separate(ID, into = c("season", "WATER")) %>%
  mutate(season = tolower(season)) %>% 
  select(-NMDS1, -NMDS2)

# Visualize NMDS with k-means clusters using ggplot2
sites = as.data.frame(as.matrix(nmds$points)) %>%
  rownames_to_column(var = "ID") %>%
  #separate(ID, into = c("season", "water")) %>% 
  mutate(cluster = as.factor(kmeans_result$cluster))

species = as.data.frame(as.matrix(nmds$species)) %>%
  rownames_to_column(var = "ID") 


## Convex Hulls 
# Compute the convex hull for each cluster
hull_data = sites %>%
  group_by(cluster) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()

## PERMANOVA ---------------------------------------------------------

## Checks to see if these site clusters are significant

# Run PERMANOVA using the distance matrix and clusters as the grouping factor
distance_matrix = vegdist(nmds.dat #%>% 
                             #select(-season, -WATER, -cluster)
                           , method = "bray")
adonis2(distance_matrix~site_scores$cluster)
 
## Consider separating this out to remove taxa that occur in all clusters or taxa that occur in neither of the paired clusters

### You should go through and compare the dS across the lat/long gradient?
#### You should also go through and see if you can compare the d15N to the d2H of the invertebrates to see if there is a realtinship between the two... trophic fractionation of d2H

## SIMPER analysis  --------------------------------------------------

simper.nmds.formatting = nmds.dat %>%
  mutate(cluster = site_scores$cluster) %>%
  arrange(cluster)

simper_result = simper(simper.nmds.formatting %>%
                          select(-cluster) %>%
                          select_if(~ sum(.) != 0) %>%
                          select_if(~ sum(.) != 20), 
                        simper.nmds.formatting$cluster)

# Assuming your SIMPER output is stored in a data frame called simper_df
# Create a data frame with taxa and their contributions
simper1 = simper_result$`1_3`$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C1vC3")

simper2 =  simper_result$`2_3`$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C3vC2")

  
simper3 =  simper_result$`1_2`$cusum %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "family") %>%
  mutate(group = "C1vC2")

simper_df =  rbind(simper1, simper2, simper3) %>%
  rename(contribution = V1) %>%
  mutate(group = case_when(group == "C1vC2" ~ "Cluster 1 v. Cluster 2",
                           group == "C1vC3" ~ "Cluster 1 v. Cluster 3", 
                           group == "C3vC2" ~ "Cluster 2 v. Cluster 3"))


## Create a better bar plot
simper.families.important = simper_df %>% 
  group_by(group) %>%
  arrange(-contribution) %>%
  slice_head(n = 15)

simper.cont1 = simper_df %>% 
  group_by(group) %>%
  arrange(-contribution) %>%
  filter(contribution == 1)
  
simpr.families = c(simper.families.important$family, simper.cont1$family) %>%
  unique()

taxon_frame = read.csv("Data/CSVs/taxon_frame.csv")









## Chemistry -------------


load("Data/RData/chem_data.RData")

chemistry = chemistry %>%
  rename(depth.5mgL = min_depth,
         temp.5mgL = min_temp) %>%
  filter(season == "fall") %>%
  left_join(cluster.mat) 
## species richness trends

richness = full %>%
  select(WATER, MONTH, FAMILY) %>% ## I've been calculating richness by genus. Now switching to family
  group_by(WATER, MONTH) %>%
  unique() %>%
  summarize(richness = n()) %>%
  rename("season" = "MONTH") %>%
  mutate(season = tolower(season)) %>%
  left_join(chemistry)






## Chemistry ---------------------------------------------------------------




## excluding SO4 because I'm not quite convinced that the measurements are the same
test = read.csv(file = "Data/CSVs/richness.csv")

richness = richness %>%
  ungroup() %>%
  mutate(DOC_update = as.numeric(test$DOC_update))


richness.simple = richness %>%
  ungroup() %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4, -Pond_num, -Lake, 
         -TotalP, -ANC, -CA, -FUI.num, -NH4, -FieldPh, -NO3, 
         -SIO2) %>%
  select(cluster, everything()) %>%
  ungroup() %>%
  pivot_longer(4:14, names_to = "metric", values_to = "values") %>%
  select(-season) %>%
  unique() %>% 
  na.omit()

## Env fit for chemistry

env_data = richness.simple %>% 
  pivot_wider(names_from = metric, values_from = values ) %>%
  na.omit() %>%
  select(-cluster) %>%
  column_to_rownames(var = 'WATER')

env_fit = envfit(nmds, env_data)
env_fit.vectors = cbind(env_fit$vectors$arrows %>% as.data.frame(), env_fit$vectors$pvals) %>%
  rename(pvals = `env_fit$vectors$pvals`) %>%
  filter(pvals < .05) %>%
  mutate(NMDS1 = NMDS1 * 2, 
         NMDS2 = NMDS2 * 2) %>%
  mutate(Metric = rownames(.)) %>%
  mutate(Metric = case_when(Metric == "richness" ~ "Richness", 
                            Metric == "surface_area" ~ "Surface area (ha)", 
                            Metric == "thermo_depth" ~ "Thermocline depth (m)"))










### Taxa by order across clusters

write.csv(richness, file = "Data/CSVs/richness.csv")

### Functional feeding groups across clusters


full %>%
  filter(WATER == "MSL") %>%
  select(MONTH, ORDER..OR.ABOVE., FAMILY, GENUS) %>%
  unique() %>% write.csv("Data/CSVs/Moss_Invert_Taxa.csv")


# Plots -----------------------


nmds.labels = c("Deep thermocline", "Intermediate thermocline", "Shallow thermocline")


### Figure 1 -----------------------
# NMDS with convex hulls
ggplot(sites, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = factor(cluster)), size = 2) +

  labs(color = "Cluster", fill = "Cluster") +
  theme_minimal() +
  theme(legend.position = "right") + 
    geom_text_repel(data = species %>%
              filter(ID %in% simper.families.important$family), aes(x = MDS1, y = MDS2, label = ID)) + 
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)], labels = nmds.labels) +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)], labels = nmds.labels) + 
   geom_label(data = env_fit.vectors %>% 
                filter(Metric != "Richness"), aes(x = NMDS1 + .25, y = NMDS2 + .1, label = Metric)) + 
  geom_segment(data = env_fit.vectors %>%
                 filter(Metric != "Richness"), aes(yend = NMDS2, xend = NMDS1, x = 0, y = 0),
               arrow = arrow()) +
   geom_polygon(data = hull_data, aes(fill = factor(cluster)), alpha = 0.3, color = NA) +
  xlim(-1.5, 2.5) 

### Figure 2 -----------------------
## Significant env.fit metrics by cluster 
richness.simple %>%
  filter(metric %in% rownames(env_fit.vectors)) %>%
  mutate(Metric = case_when(metric == "richness" ~ "Richness", 
                            metric == "surface_area" ~ "Surface area (ha)", 
                            metric == "thermo_depth" ~ "Thermocline depth (m)")) %>%
  ggplot(aes(x = cluster, y = values, fill = cluster)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 12) +
  facet_wrap(~Metric, scales = "free") + 
  geom_point() +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  theme(legend.position = "none") + 
  xlab("Cluster") +
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 0, hjust = 0)) +
  scale_x_discrete(labels = c(
  "1" = "Deep\nthermocline", 
  "2" = "Intermediate\nthermocline", 
  "3" = "Shallow,\nthermocline"
)) + 
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5))
  
  
  
### Figure 3 ------------------------  
nmds.dat %>%
  rownames_to_column(var = "ID") %>%
  left_join(cluster.mat %>% 
              mutate(season = toupper(season)) %>%
              unite("ID", c(season, WATER))) %>%
  select(ID, cluster, everything()) %>%
  pivot_longer(ceratopogonidae:lymnaeidae) %>%
  ungroup() %>%
  filter(value ==1) %>%
  group_by(cluster) %>%
  mutate(n = length(unique(ID))) %>%
  ungroup() %>%
  group_by(cluster, name, n) %>%
  summarize(total_waters = n()) %>%
  mutate(prop.clust = total_waters/n) %>%
  ungroup() %>% 
  group_by(name) %>%
  mutate(total.clust = sum(prop.clust)) %>%
  filter(name %in% simpr.families) %>%
  left_join(taxon_frame %>%
              select(FAMILY, ORDER) %>%
              unique(), by = c("name" = "FAMILY"))  %>%
  ungroup() %>%
  arrange(ORDER, -total.clust) %>%
  mutate(count = 1:n()) %>%
  ggplot(aes(x = prop.clust , y = reorder(name, -count), fill = cluster)) + 
  geom_bar(stat = "identity") + 
    scale_fill_manual("Cluster", values = wes_palette("Royal2")[c(3,1,5)], labels = nmds.labels) +
  theme_minimal(base_size = 14) +
  geom_point(aes(y = name, color = ORDER), size = 3, x = -0.05, shape = "circle",inherit.aes = FALSE) +
  scale_color_manual("Order", values = wes_palette("Darjeeling1",
                                         type = "continuous", n = 13)) +
  
  theme(axis.title.y = element_blank()) + 
  xlab("Frequency within Cluster") 

  
