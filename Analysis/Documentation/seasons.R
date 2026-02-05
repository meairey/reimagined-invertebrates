set.seed(1234)  # For reproducibility
library(tidytext)
library(vegan)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(wesanderson)
`%nin%` = Negate(`%in%`)

# Load
load(file = "Data/RData/FullData.RData")


## NMDS data
nmds.dat = full %>% 
  filter((WATER != "DTL" | SITE_N %nin% c(2,4)),
         (WATER !="MSL" | SITE_N %nin% c("MSL.002", "MSL.004")),
         WATER %nin% c("COM", "LML","ETL","GNL")) %>%
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

nmds.full = metaMDS(nmds.dat, distance = "jaccard", trymax = 100)


## K means clustering of lakes by community
# extract scores from NMDS
site_scores = as.data.frame(scores(nmds.full, display = "sites"))


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
cluster.mat = site_scores %>% rownames_to_column(var = "ID") %>%
  separate(ID, into = c("season", "WATER")) %>%
  mutate(season = tolower(season)) %>% 
  select(-NMDS1, -NMDS2)

# Visualize NMDS with k-means clusters using ggplot2
sites = as.data.frame(as.matrix(nmds.full$points)) %>%
  rownames_to_column(var = "ID") %>%
  separate(ID, into = c("season", "water")) %>% 
  mutate(cluster = as.factor(kmeans_result$cluster)) %>%
  filter(water %nin% c("COM", "LML", "ETL", "GNL"))

## Investigating how sites change by season

centroids = sites %>% 
  group_by(cluster, season) %>%
  summarize(centroid.x = mean(MDS1),
            centroid.y = mean(MDS2))
ggplot() +
  theme_minimal(base_size = 14) +
  geom_point(data = sites, aes(x = MDS1, y = MDS2, color = cluster, shape = season), size = 2, alpha = .3) + 
  geom_segment( aes(x = (sites %>%
                                   filter(season == "SPRING"))$MDS1, 
                                 y = (sites %>%
                                   filter(season == "SPRING"))$MDS2,
                                 xend = (sites %>%
                                   filter(season == "FALL"))$MDS1, 
                                 yend = (sites %>%
                                   filter(season == "FALL"))$MDS2), alpha = .3, lty = 3) +
  geom_point(data = centroids, aes(x = centroid.x, y = centroid.y, color = cluster, shape = season), size = 5) + 
  geom_segment( aes(x = (centroids %>% filter(season == "SPRING"))$centroid.x, 
                    y =(centroids %>% filter(season == "SPRING"))$centroid.y,
                    xend = (centroids %>% filter(season == "FALL"))$centroid.x,
                    yend = (centroids %>% filter(season == "FALL"))$centroid.y),
                arrow = arrow(type = "closed"), 
                lwd = 1) + 
  scale_color_manual("Cluster", values = wes_palette("Darjeeling1", 3)) + 
  scale_shape_manual("Season",values = c(16,17), labels = c("Fall", "Spring"))


## Sites that change cluster membership
sites %>% 
  select(season, water, cluster) %>%
  pivot_wider(names_from = season, values_from = cluster) %>%
  mutate(same = FALL==SPRING)


## Permanova to look at if variation is greater between clusters or between seasons
dist_matrix = vegdist(nmds.dat > 0, method = "jaccard") ## Distance matrix for the NMDS dat 

# Individual water or cluster/season more important?
adonis2(dist_matrix ~ water + cluster + season, data = sites, permutations = 999)
# Interaction between cluster and season
adonis2(dist_matrix ~ cluster * season, data = sites,
        permutations = how(nperm = 999, blocks = sites$water))

## Environmental Data 

load("Data/RData/chem_data.RData")

chemistry = chemistry %>%
  rename(depth.5mgL = min_depth,
         temp.5mgL = min_temp) %>%
  #filter(season == "fall") %>%
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
 # select(-season) %>%
  unique() %>% 
  na.omit() %>%
  filter(WATER %nin% c("LML", "COM","GNL","ETL"))

## Env fit with seasons?
env_data = richness.simple %>% 
  unite("ID", WATER:season) %>%
  pivot_wider(names_from = metric, values_from = values ) %>%
  na.omit() %>%
  select(-cluster) %>%
  column_to_rownames(var = 'ID')

env_fit = envfit(nmds.full, env_data)
env_fit.vectors = cbind(env_fit$vectors$arrows %>% as.data.frame(), env_fit$vectors$pvals) %>%
  rename(pvals = `env_fit$vectors$pvals`) %>%
  filter(pvals < .05) %>%
  mutate(NMDS1 = NMDS1 * 2, 
         NMDS2 = NMDS2 * 2) %>%
  mutate(Metric = rownames(.)) %>%
  mutate(Metric = case_when(Metric == "richness" ~ "Richness", 
                            Metric == "surface_area" ~ "Surface area (ha)", 
                            Metric == "thermo_depth" ~ "Thermocline depth (m)"))




## Envfit Table Supplement
env_fit$vectors$arrows %>% 
  as.matrix()  %>%
  as.data.frame() %>%
  mutate(p.val = env_fit$vectors$pvals, 
         r = round(env_fit$vectors$r, digits = 2)) %>%
  arrange(p.val, r) %>%
  rownames_to_column(var = "Metric") %>%
  filter(Metric %nin% c("DOC.1", "thermo_temp")) %>%
  rename(`Pr(>r)` = "p.val", 
         `R2` = "r") %>%
  mutate(Metric = str_replace(Metric, "richness", "Richness"),
         Metric = str_replace(Metric, "thermo_depth", "Thermocline Depth"),
         Metric = str_replace(Metric, "surface_area", "Surface Area"), 
         Metric = str_replace(Metric, "max_depth", "Maximum Depth"),
         Metric = str_replace(Metric, "depth.5mgL", "DD05"),
         Metric = str_replace(Metric, "temp.5mgL", "TDO5 Depth"),
         Metric = str_replace(Metric, "sechi.depth", "Secchi Depth"), 
         Metric = str_replace(Metric, "DOC_update", "DOC")) -> Table.S2               
