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
nmds <- metaMDS(nmds.dat, distance = "jaccard", trymax = 100)


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
sites %>% ggplot(aes(x = MDS1, 
                     y = MDS2, 
                     col = cluster, 
                     label = ID)) + 
  theme_minimal(base_size = 10) +
  stat_ellipse(level = .9, alpha = .8) +
  geom_text() +
  labs(col = "Cluster") + 
  xlim(-1.5, 1.5) +
  
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) 

## Which species are pulling the NMDS 
species %>% ggplot(aes(x = MDS1, y = MDS2, label = ID)) + 
  geom_text() 







## PERMANOVA to see if these site clusters are significant

# Run PERMANOVA using the distance matrix and clusters as the grouping factor
distance_matrix <- vegdist(nmds.dat #%>% 
                             #select(-season, -WATER, -cluster)
                           , method = "bray")
adonis2(distance_matrix~site_scores$cluster)
 
## Consider separating this out to remove taxa that occur in all clusters or taxa that occur in neither of the paired clusters
## SIMPER analysis to figure out which taxa influence this the most
simper_result <- simper(nmds.dat %>%
                          select_if(~ sum(.) != 0) %>%
                          select_if(~ sum(.) != 20), 
                        site_scores$cluster)

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

library(tidytext)
# Create a bar plot
simper_df %>%
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

simper_df %>%
  left_join(order_matrix, by = c("family" = "FAMILY")) %>%
  ggplot(aes(x = ORDER, y = contribution, fill = ORDER)) +
  geom_boxplot() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = .5),
        axis.title.x = element_blank(), 
        legend.position = "none") +
  scale_fill_manual(values = wes_palette("Darjeeling1",
                                         type = "continuous", n = 14)) +
  ylab("Contribution to Dissimilarity") + 
  
  geom_segment(aes(x = "hirudinea", xend = "odonata", y = 1.1)) + 
  geom_segment(aes(x = "coleoptera", xend = "odonata", y = 1.14)) + 
  geom_segment(aes(x = "odonata", xend = "trichoptera", y = 1.18)) +
  geom_text(aes(label = "*", x = "hirudinea", y = 1.11), size = 4) +
  geom_text(aes(label = "**", x = "coleoptera", y = 1.15), size = 4) +
  geom_text(aes(label = "**", x = "trichoptera", y = 1.19), size = 4)





## Significant differences across contribution by order
dis.aov = simper_df %>%
  left_join(order_matrix, by = c("family" = "FAMILY")) %>% 
  aov(data = ., contribution ~ ORDER) 
summary(dis.aov)
## tukey test
tukey.dis = (TukeyHSD(dis.aov))  

tukey.dis$ORDER %>% as.data.frame() %>%
  filter(`p adj` < .05)







## Contribution of different odonata to differences across clusters
simper_df %>%
  left_join(order_matrix, by = c("family" = "FAMILY")) %>%
  filter(ORDER %in% c("gastropoda", "bivalvia")) %>% #"ephemeroptera"  "odonata" "trichoptera"
  group_by(family, ORDER) %>% 
  summarize(contribution = mean(contribution)) %>%
  ggplot(aes(x = ORDER, y = contribution, fill = reorder(family, contribution))) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) + 
  labs(fill = "Family")



library(tidyverse)


## Chemistry -------------


load("Data/RData/chem_data.RData")
chemistry = chemistry %>%
  rename(depth.5mgL = min_depth,
         temp.5mgL = min_temp) %>%
  left_join(cluster.mat)
## species richness trends

richness = full %>%
  select(WATER, MONTH, FAMILY, GENUS) %>%
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

taxa.assignments = nmds.dat.cluster %>% 
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
  filter(proportion.lakes == max(proportion.lakes)) %>% 
  select(cluster, name,proportion.lakes) %>%
  rename(FAMILY = name)

taxa.assignments 

## Which species are pulling the NMDS 
species %>% 
  left_join(taxa.assignments, by = c("ID" = "FAMILY")) %>%
  
  ggplot(aes(x = MDS1, y = MDS2,  col = cluster)) + 
  geom_text(aes(label = ID), max.overlaps = 40) +
  stat_ellipse() + 
  theme_minimal() 

taxa.assignments %>% left_join(taxon_frame) %>%
  ggplot(aes(x = ORDER, fill = cluster)) +
  geom_bar() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = .5))

taxa.assignments %>% left_join(taxon_frame) %>%
  #rename(ORDER = ORDER..OR.ABOVE.) %>%
  select(cluster, ORDER) %>%
  unique() %>%
  ggplot(aes(fill = ORDER, x = cluster)) +
  geom_bar() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = .5)) +
  scale_fill_manual(values = wes_palette("Darjeeling1", type = "continuous", n = 14)) +
  ylab("Number of Families")

## Trying to labels these

counts <- taxa.assignments %>%
  left_join(taxon_frame)  %>%
 # rename(ORDER = ORDER..OR.ABOVE.)%>%
  select(cluster, FAMILY, ORDER) %>%
  unique() %>%
  group_by(cluster, ORDER) %>%
  summarise(count = n()) %>%
  ungroup()

# Then, plot the bar graph with text labels
ggplot(counts, aes(fill = ORDER, x = cluster, y = count)) +
  geom_bar(stat = "identity") +  # Change to stat = "identity" since counts are pre-calculated
  geom_text(aes(label = ORDER), position = position_stack(vjust = 0.5), size = 3) + # Add text labels
  theme_minimal() + 
  scale_fill_manual(values = wes_palette("Darjeeling1", type = "continuous", n = 14)) +
  ylab("Number of Families") +
  theme(legend.position = "none") +
  xlab("Cluster")

## excluding SO4 because I'm not quite convinced that the measrements are the same
richness.simple = richness %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4) %>%
  select(cluster, everything()) %>%
  ungroup() %>%
  pivot_longer(7:22, names_to = "metric", values_to = "values")

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

table = table %>% na.omit() %>% as.data.frame() %>%
  mutate(metric.p = round(as.numeric(metric.p), digits = 3),
         `3_2.p` = round(as.numeric(`3_2.p`), digits = 3),
         `3_1.p` = round(as.numeric(`3_1.p`), digits = 3),
         `2_1.p` = round(as.numeric(`2_1.p`), digits = 3)) %>%
  mutate(across(everything(), ~ ifelse(. < 0.001, "< 0.001", .)))

#write.csv(table, "Data/CSVs/anova_summary.csv")

## Visualize just the significant factors 
richness.simple %>% filter(metric %in% unique(table$metric)) %>%
  ggplot(aes(x = cluster, y = values, fill = cluster)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 12) +
  facet_wrap(~metric, scales = "free") +
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) 
  
## Richness vs. cluster

richness %>% 
  ggplot(aes(x = cluster, y = richness, fill = cluster)) + 
  geom_boxplot() + 
  theme_minimal(base_size = 12) + 
  xlab("Cluster") + ylab("Richness") +
  
  scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) 

aov(data = richness, richness ~ cluster) %>% summary()
TukeyHSD(aov(data = richness, richness ~ cluster))



### Taxa by order across clusters

write.csv(richness, file = "Data/CSVs/richness.csv")

### Functional feeding groups across clusters


