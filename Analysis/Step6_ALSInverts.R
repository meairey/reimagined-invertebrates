## ALS Invert Data

### Actual ALS invert data
#### Filter out just method 1 


als.lookup = read.csv("Data/CSVs/ALS_invert_lookup.csv") %>%
  select(-Approach_key_name) 

als_order = als.lookup %>% filter(TAXONOMIC_GROUP == "ORDER") %>%
  rename(ORDER = CODE, 
         ORDER.full = TAXONOMIC_ID) %>%
  select(-TAXONOMIC_GROUP) 
als_family = als.lookup %>% filter(TAXONOMIC_GROUP == "FAMILY")%>%
  rename(FAMILY = CODE,
          FAMILY.full = TAXONOMIC_ID)%>%
  select(-TAXONOMIC_GROUP)



als.inv = read.csv("Data/ALS_INVERT_REC7.csv") %>% left_join(als_order) %>%
  left_join(als_family, by = "FAMILY") %>%
  mutate(FAMILY.full = case_when(ORDER.full == "Amphipoda" ~ "Amphipoda.unsp", 
                                 ORDER.full == "Isopoda" ~ "Isopoda.unsp",
                                 FAMILY.full == "Agriidae"  ~ "Calopterygidae",
                                 ORDER.full %nin% c("Amphipoda.unsp", "Isopoda.unsp") ~ 
                                   FAMILY.full)) %>%
  filter(FAMILY.full != "Spongillidae", ## Remove freshwater sponges
         ORDER.full != "Collembola", ## Remove springtails which aren't really aquatic
         ORDER.full != "Decapoda", ## Removing Crayfish
         ORDER.full != "Unspecified",
         FAMILY.full != "Unspecified",
         FAMILY.full != "Chrysomelidae", ## Removing leaf beetles (not typically aquatic)
         FAMILY.full != "Curclionidae", ## Removing weevils (not typically aquatic)
         FAMILY.full != "Muscidae") ## Removing house flies (not typically aquatic)

als_chem = read.csv("Data/CSVs/ALSC_chem_updated.csv")
## Read in thermocline data
thermocline = read.csv("Data/CSVs/thermocline_data_ALS.csv") %>%
  select(PONDNO, curve_fitted_thermocline_depth, SECCHI)

als.inv %>% group_by(PONDNO) %>%
 # reframe(unique(FAMILY.full)) %>%
  summarize(family_richness = length(unique(FAMILY.full))) %>%
  left_join(als_chem, by = c("PONDNO" = "Pond_number")) %>%
    ggplot(aes(x = family_richness, y = FIELDPH)) + 
  geom_point() + 
  geom_smooth()


als.inv %>% 
  filter(NMETHOD == 5) %>%
  select(FAMILY.full) %>%
  unique()

als.inv %>% 
  group_by(NMETHOD) %>%
  summarize(count = n())

pres.abs = als.inv %>% 
  ## Testing out thermocline depth as a variable
  filter(PONDNO %in% thermocline$PONDNO) %>%
  filter(NMETHOD == 1) %>%
  select(PONDNO, FAMILY.full) %>%
  filter(FAMILY.full != "Unspecified") %>%
  group_by(PONDNO) %>%
  unique() %>% 
  mutate(pres = 1) %>%
  ungroup() %>%
  group_by(FAMILY.full) %>%
  mutate(total_occur = n()) %>%
  filter(total_occur > 10) %>% 
  ungroup() %>%
  select(-total_occur) %>%
  group_by(PONDNO) %>%
  mutate(rich = length(unique(FAMILY.full))) %>%
  filter(rich > 3) %>% ## Mean 7.25, sd = 3.1, low = 2, high = 19
  select(-rich) %>%
  ungroup() %>%
  pivot_wider(names_from = FAMILY.full, values_from = pres, values_fill = 0) %>%
  filter(PONDNO %nin% c("020329", "040768")) %>%
  column_to_rownames(var = "PONDNO") 

dim(pres.abs)
stress = vector()
for(i in 4:10){
  i = 5
  nmds = metaMDS(pres.abs, distance = "jaccard", k = i, trymax = 20) 
  stress[i] = nmds$stress
}

ggplot(mapping = aes(x = 1:length(stress), y = stress)) + 
  geom_point()



site_scores = as.data.frame(scores(nmds, display = "sites"))


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

## Optimal clusters looks like k = 5

k = 3  # Number of clusters you want to try
kmeans_result = kmeans(site_scores, centers = k)

# Add the cluster results to the site_scores data frame
site_scores$cluster = as.factor(kmeans_result$cluster)

cluster.mat = site_scores %>% rownames_to_column(var = "ID") %>%
  select(-NMDS1, -NMDS2, -NMDS3)

cluster.mat %>% group_by(cluster) %>%
  summarize(n())
cluster.mat %>% filter(cluster ==3)

# Visualize NMDS with k-means clusters using ggplot2
sites = as.data.frame(as.matrix(nmds$points)) %>%
  rownames_to_column(var = "ID") %>%
  #separate(ID, into = c("season", "water")) %>% 
  mutate(cluster = as.factor(kmeans_result$cluster))



species = as.data.frame(as.matrix(nmds$species)) %>%
  rownames_to_column(var = "ID") 

species %>% 
ggplot(aes(x = MDS1, y = MDS2, label = ID)) + 
  geom_text()



sites %>%
  ggplot(aes(x = MDS1, y = MDS2, col = cluster)) +
  geom_point()


env_data= read.csv("Data/CSVs/ALSC_chem_updated.csv") %>% 
  filter(Pond_number %in% rownames(pres.abs)) %>%
  select(Pond_number, SO4:SCONDUCT) %>%
  select(-TAL)
physical_data = read.csv("../solid-fishstick/Data/ALS_fish_updated.csv") %>%
  select(Pond_number:Shoreline_Length_km) 

lake_data = left_join(env_data, physical_data) %>% 
  left_join(thermocline, by = c("Pond_number" = "PONDNO")) %>% 
  filter(Pond_number %in% rownames(pres.abs))

env_fit = envfit(nmds, lake_data, na.rm =T)
env_fit.vectors = cbind(env_fit$vectors$arrows %>% as.data.frame(), env_fit$vectors$pvals) %>%
  mutate( r = round(env_fit$vectors$r, digits = 2)) %>%
  rename(pvals = `env_fit$vectors$pvals`) %>%
  filter(pvals < .05) %>%
  mutate(NMDS1 = NMDS1 * 2, 
         NMDS2 = NMDS2 * 2) %>%
  mutate(Metric = rownames(.)) 

env_fit.vectors


#### Certain taxa that might be associated with environmental variables

test = pres.abs %>% 
  select(Heptageniidae, Caenidae) %>%
  rownames_to_column(var = "PONDNO") %>%

  left_join(lake_data, by = c("PONDNO" = "Pond_number"))

test %>% 
  ggplot(aes(x = Heptageniidae, y = SECCHI)) + 
  geom_point() + 
  geom_smooth(method = "lm")


library(lme4)
glm(Heptageniidae ~ Surface_Area_ha + SECCHI   , family = binomial, data = test) %>% summary()

glm(Caenidae ~ SECCHI + K  , family = binomial, data = test) %>% summary()













