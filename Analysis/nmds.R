library(tidytext)
library(vegan)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(wesanderson)
`%nin%` = Negate(`%in%`)
set.seed(123)  # For reproducibility
set.seed(695)

# Load
load(file = "Data/RData/FullData.RData")


full %>%
  summarize(total = sum(TOTAL_N))
  
## create the NMDS data frame ------------------

## Families included 
unique_families = full %>% select(FAMILY, GENUS) %>%  unique() 
write.csv(file = "Data/CSVs/family.csv", unique_families, row.names = F)

## NMDS data
nmds.dat = full %>% select(MONTH, YSAMP, SITE_N, WATER, ORDER..OR.ABOVE., FAMILY) %>%
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

save(file = "Data/RData/nmdsDat.RData", nmds.dat)

order_matrix = full %>%
  select(MONTH, YSAMP, SITE_N, WATER, ORDER..OR.ABOVE., FAMILY) %>%
  unique() %>%
  filter(FAMILY %nin% c("unidentified", "degraded","terrestrial")) %>%
  mutate(SITE_N = str_replace_all(SITE_N, " ", ""),
         SITE_N = str_replace_all(SITE_N, "Notlisted", "NotListed")) %>% 
  select( ORDER..OR.ABOVE., FAMILY) %>%
  arrange(ORDER..OR.ABOVE., FAMILY) %>%
  rename(ORDER = ORDER..OR.ABOVE.) %>% unique()


nmds.dat ## Frame is presence absence, sites are rows. columns are families.


# Run the NMDS using Bray-Curtis distance ------------------------------
#nmds <- metaMDS(nmds.dat, distance = "jaccard", trymax = 100) ## Do not run this again just load the data out of file below

#save(nmds, file = "Data/RData/nmds.RData") ## Do not run this again just load the data out of the file below



load(file = "Data/RData/nmds.RData")
## K means clustering of lakes by community
# extract scores from NMDS
site_scores <- as.data.frame(scores(nmds, display = "sites"))


# Run for loop and use elbow method to determine optimal clusters --------------
## test k 
# 
wss <- numeric()

# Use a for loop to run k-means for k = 1 to 10
for (k in 1:10) {
  kmeans_result <- kmeans(site_scores, centers = k, nstart = 10)
  wss[k] <- kmeans_result$tot.withinss  # Store the total within-cluster sum of squares
}

# plot ------------
elbow_plot <- data.frame(k = 1:10, WSS = wss)

ggplot(elbow_plot, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Elbow Method for Finding Optimal K",
       x = "Number of Clusters (K)",
       y = "Total Within-Cluster Sum of Squares") +
  theme_minimal()

## Optimal clusters looks like k = 3

k <- 3  # Number of clusters you want to try
kmeans_result <- kmeans(site_scores, centers = k)

# Add the cluster results to the site_scores data frame
site_scores$cluster <- as.factor(kmeans_result$cluster)
save(file = "Data/RData/site_scores.RData", site_scores)
load(file = "Data/RData/site_scores.RData")
load(file = "Data/RData/simmr.output.RData")
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

## View sites and clusters

library(ggrepel)
ggplot() +
  
  geom_point(data= sites, aes(x = MDS1, 
                     y = MDS2, 
                     col = cluster),key_glyph = "rect") +
  stat_ellipse(data = sites, aes(x = MDS1, 
                     y = MDS2, 
                     col = cluster, 
                     ),level = .9, alpha = .8) +
  labs(col = "Cluster") + 
  xlim(-1.5, 1.5) +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  geom_text_repel(data = species %>%
              filter(ID %in% simper.families.important$family), aes(x = MDS1, y = MDS2, label = ID)) + 
  theme_minimal(base_size = 12)


## Convex Hulls 
# Compute the convex hull for each cluster
hull_data <- sites %>%
  group_by(cluster) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()

nmds.labels = c("Shallow thermocline", "Diverse, deep thermocline", "Deep thermocline")
# Plot with convex hulls
ggplot(sites, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = factor(cluster)), size = 2) +
  geom_polygon(data = hull_data, aes(fill = factor(cluster)), alpha = 0.3, color = NA) +
  labs(color = "Cluster", fill = "Cluster") +
  theme_minimal() +
  theme(legend.position = "right") + 
    geom_text_repel(data = species %>%
              filter(ID %in% simper.families.important$family), aes(x = MDS1, y = MDS2, label = ID)) + 
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)], labels = nmds.labels) +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)], labels = nmds.labels) 
  


## PERMANOVA ---------------------------------------------------------

## Checks to see if these site clusters are significant

# Run PERMANOVA using the distance matrix and clusters as the grouping factor
distance_matrix <- vegdist(nmds.dat #%>% 
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

simper_result <- simper(simper.nmds.formatting %>%
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

nmds.dat %>%
  rownames_to_column(var = "ID") %>%
  left_join(cluster.mat %>% 
              mutate(season = toupper(season)) %>%
              unite("ID", c(season, WATER))) %>%
  select(ID, cluster, everything()) %>%
  pivot_longer(ceratopogonidae:unionidae) %>%
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







## Chemistry -------------


load("Data/RData/chem_data.RData")

chemistry = chemistry %>%
  rename(depth.5mgL = min_depth,
         temp.5mgL = min_temp) %>%
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


save(file = "Data/RData/richness_cluster.RData", richness)

## Grouped by cluster -------------------------

## Assigning taxa to different clusters -------

nmds.dat.cluster = nmds.dat %>% rownames_to_column(var = "ID") %>%
  separate(ID, into = c("season","WATER")) %>%
  mutate(season = tolower(season)) %>%
  left_join(cluster.mat)

## Which taxa occur at any point in clusters 

FFG = read.csv("Data/CSVs/FFGS.csv") %>%
  select(FAMILY, FG) %>%
  unique() %>%
  na.omit()
taxa_presence = nmds.dat.cluster %>% 
  ungroup() %>%
  select(season, cluster, everything()) %>%
  pivot_longer(4:49) %>%
  group_by(cluster) %>%
  unite("ID", c(season, WATER)) %>%
  mutate(total_lakes = length(unique(ID))) %>%
  ungroup() %>%
  select(cluster, name, total_lakes, value, ID) %>%
  unique() %>%
  group_by(cluster, name, total_lakes) %>%
  summarize(sum = sum(value)) %>%
  mutate(proportion.lakes = sum/total_lakes) %>%
  ungroup() %>%
  group_by(name) %>%
  filter(sum > 0) %>% 
  left_join(FFG, by = c("name" = "FAMILY"))

taxa_presence %>% 
  group_by(cluster) %>%
  mutate(count = length(unique(name))) %>%

  ungroup() %>% 
  group_by(cluster, count, FG) %>%
  summarize(ffg = n()) %>%
  mutate(prop = ffg / count) %>%
  ungroup() %>%
  select(cluster, FG, prop) %>%
  print(n = 100)
  
  
  mutate(prop = count / total_lakes)


write.csv(taxa_presence, file = "Data/CSVs/taxa.presence.csv")




## Chemistry ---------------------------------------------------------------




## excluding SO4 because I'm not quite convinced that the measrements are the same
richness.simple = richness %>%
  ungroup() %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4, -Pond_num, -Lake, -TotalP) %>%
  select(cluster, everything()) %>%
  ungroup() %>%
  pivot_longer(4:20, names_to = "metric", values_to = "values") %>%
  select(-season) %>%
  unique()

richness.simple %>%
  ggplot(aes(x = cluster, y = values)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 12) +
  facet_wrap(~metric, scales = "free")

table = matrix(NA, nrow = 16, ncol = 5)
colnames(table) = c("metric", "metric.p", "2_1.p","3_1.p", "3_2.p")

for(i in 1:length(unique(richness.simple$metric))){
  x = richness.simple %>% 
    filter(metric == unique(richness.simple$metric[i]))
  
  aov.object = aov(x$values ~ x$cluster) 
  aov.object.sum  = aov.object %>% summary
  if(aov.object.sum[[1]]$`Pr(>F)`[1] < .05){
    
    table[i,1] = (unique(richness.simple$metric[i]))
    table[i,2] = (aov.object.sum[[1]]$`Pr(>F)`[1])
    table[i,3] = (TukeyHSD(aov.object))$`x$cluster`[1,4]
    table[i,4] = (TukeyHSD(aov.object))$`x$cluster`[2,4]
    table[i,5] = (TukeyHSD(aov.object))$`x$cluster`[3,4]
    
  }
  
  
}

##Kruscal
table = matrix(NA, nrow = 17, ncol = 5)
colnames(table) = c("metric", "metric.p", "2_1.p","3_1.p", "3_2.p")


for(i in 1:length(unique(richness.simple$metric))){
  x = richness.simple %>% 
    filter(metric == unique(richness.simple$metric[i]))
  
  aov.object = kruskal.test(x$values ~ x$cluster) 
  
  
  aov.object.sum  = aov.object$p.value
  
  if(aov.object.sum < .05){
    print( unique(richness.simple$metric[i]))
    
    table[i,1] = (unique(richness.simple$metric[i]))
    table[i,2] = aov.object.sum
    
    posthoc <- pairwise.wilcox.test(x$values, as.factor(x$cluster), p.adjust.method = "bonferroni")
    table[i,3] = posthoc$p.value[1,1]
    table[i,4] = posthoc$p.value[2,1]
    table[i,5] = posthoc$p.value[2,2]
    


    
  }
  
  
}



table = table %>% na.omit() %>% as.data.frame() %>%
  mutate(metric.p = round(as.numeric(metric.p), digits = 3),
         `3_2.p` = round(as.numeric(`3_2.p`), digits = 3),
         `3_1.p` = round(as.numeric(`3_1.p`), digits = 3),
         `2_1.p` = round(as.numeric(`2_1.p`), digits = 3)) %>%
  mutate(across(everything(), ~ ifelse(. < 0.001, "< 0.001", .))) %>%
  as.data.frame()


table 

#write.csv(table, "Data/CSVs/anova_summary.csv")

## Visualize just the significant factors '
richness.simple$variable_level <- factor(richness.simple$metric,
                                         levels = c("ANC" "CA", "SIO2","NO3",
                                                    "TotalP","DOC.1", "Elevation",
                                                    "depth.5mgL", "temp.5mgL",
                                                    "thermo_depth", "richness", "sechi.depth"))


aov_variable_names = c("ANC" = "ANC (µeq/L)", "CA" = "Ca (mg/L)", "SIO2" = "SiO2 (mg/L)",
                       "NO3" = "NO3 (mg/L)",
                       "TotalP" = "Total P (mg/L)", "DOC.1" = "DOC (mg/L)",
                       "Elevation" = "Elevation (m)", "depth.5mgL" = "Depth TD5", 
                       "temp.5mgL" = "Temp TD5", 
                       "thermo_depth" = "Thermocline depth (m)",
                       "richness" = "Family Richness" )


richness.simple %>%
  select(-variable_level) %>%
  na.omit() %>%
  filter(metric %in% unique(table$metric)) %>%

  unique() %>%
 mutate(metric = factor(metric, levels = c("CA", "NO3","SIO2", 
                                               "Elevation","max_depth" ,
                                               "depth.5mgL", "sechi.depth", "temp.5mgL", "thermo_depth", "richness"))) %>%
  ggplot(aes(x = cluster, y = values, fill = cluster)) + 
  geom_boxplot(alpha = .85) + 
  theme_minimal(base_size = 12) +

  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  theme(legend.position = "none") + 
  xlab("Cluster") +
  theme(axis.title.y = element_blank()) +
  geom_jitter()






### Taxa by order across clusters

write.csv(richness, file = "Data/CSVs/richness.csv")

### Functional feeding groups across clusters


full %>%
  filter(WATER == "MSL") %>%
  select(MONTH, ORDER..OR.ABOVE., FAMILY, GENUS) %>%
  unique() %>% write.csv("Data/CSVs/Moss_Invert_Taxa.csv")

