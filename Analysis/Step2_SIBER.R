# Isotopic Niches
## Script setup -------------
`%nin%` = Negate("%in%")
set.seed(123)

## Libraries
library(SIBER)
library(tidyverse)
library(wesanderson)

## Custom functions for ellipse metric generation and plotting
source("../solid-fishstick/Analysis/Documentation/isotope_functions_update.R") 


## Load in taxon data
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Cluster and chemistry data
cluster_chem = read.csv(file = "Data/CSVs/richness_update.csv") %>%
  #unite("ID", WATER, season, sep = ".") %>%
  left_join(community.legend, by = c("ID" = "community.name"))


## Join the sampling and measurement data files
sample = read.csv("../1.clean_isotope/isotope_sample.csv") ## Sampling event information
data.iso=  read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  left_join(sample, by = c("ISO_YSAMP_N")) %>%
  ## lipid correction factor
  mutate(D13C = case_when(CATEGORY %nin% c("ALGA", "PLANT") ~ as.numeric(D13C) - 3.32 + .99 * (PER_C / PER_N),
                          CATEGORY %in% c("ALGA", "PLANT") ~ as.numeric(D13C))) %>% 
  filter(GROUP != "FISH", # Remove fish f
         !(WATER == "LML" & DAY_N < 230))  # Remove extraneous sampling events 
  




data = data.iso %>% ## This one is formatted for SIBER
  filter(!(WATER == "LML" & DAY_N < 230)) %>%
  left_join(taxon_frame %>% select(TAXON, FAMILY, ORDER) %>% unique()) %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  mutate(GROUP = str_replace(GROUP, "CLAM","BIVALVE"),
         GROUP = str_replace(GROUP,"MUSSEL", "BIVALVE"))  %>%
  mutate(graph_id = case_when(is.na(ORDER) == T ~ GROUP, 
                              is.na(ORDER) == F ~ ORDER)) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall")) %>%
  unite("community.name", c(WATER, season), sep = ".") %>%
  mutate(group.name = case_when(is.na(FAMILY)== T ~ GROUP, 
                           is.na(FAMILY)==F ~ FAMILY))%>%
  group_by(community.name, group.name) %>% 
  mutate(total = n()) %>%
  filter(total > 3) %>%
  select( D13C, D15N, group.name, community.name) %>%
  rename(iso1 = D13C,
         iso2 = D15N) %>%
  na.omit() %>%
  ungroup() %>%
  arrange(community.name) %>%
  mutate(community = as.numeric(as.factor(community.name))) %>%
  filter(group.name %nin% c("LEAF", "PERI", "INSECT",
                            "ZOOP", "CRAY", "SIL")) %>% ## Filter out groups not to include in isotopic niche analysis
  arrange(group.name) %>% ## Needs to be arranged here because SIBER orders things numerically 
  mutate(group.name = tolower(group.name)) %>%
  mutate(group = as.numeric(as.factor(group.name))) %>% ## Assign your groups in alphabetical order
  arrange(community, group) %>%
  as.data.frame() %>%
  select(iso1, iso2, group, community, group.name, community.name) ## keep group and community name for legend 

legend = data %>% select(group.name,group) %>%
  unique() %>%
  mutate(color = c(1:15),
         common = group.name, 
         CODE = group.name)

community.legend = data %>% select(community, community.name) %>%
  unique() ## 


# options for running jags
parms = list()
parms$n.iter = 2 * 10^4   # number of iterations to run the model for
parms$n.burnin = 1 * 10^3 # discard the first set of values
parms$n.thin = 10     # thin the posterior by this many
parms$n.chains = 2        # run this many chains

# define the priors
priors = list()
priors$R = 1 * diag(2)
priors$k = 2
priors$tau.mu = 1.0E-3

## Generate posterior
siber.object = createSiberObject(data %>% select(-group.name, -community.name))
posterior = siberMVN(siber.object, parms, priors)
  

## Overlap analysis ------


overlap_list = list() 

communities_overlap = unique(data$community)[-9] ## For overlap remove communities with only 1 taxa

for(h in 1:length(communities_overlap)){
  
  overlap_list[[h]] = overlap(data, communities_overlap[h],50, posterior) %>% as.data.frame()
}
save(overlap_list, file = "Data/RData/overlap_list.RData") ## Save this to avoid running each time

## Looks like i had added legend into the function?, legend
# Code to reduce list to matrix 
overlap_data = Reduce(full_join, overlap_list) %>%
  select('SppPair', everything())

overlap.long = overlap_data %>% 
  pivot_longer(2:length(.[1,]),
               names_to = "Community", 
               values_to = "Values") %>% 
  na.omit() %>%
  rename("Species_Pair" = `SppPair`) %>% 
  separate(Community, into = c("c", "Community", "post")) %>%
  select(-c) %>%
  mutate(Values = as.numeric(Values)) %>%
  mutate(Species_Pair = str_replace(Species_Pair, " v ", "_")) %>%
  mutate(Values = round(Values, digits = 3))

## table 

df = overlap.long %>% ungroup() 
df$group1 = sapply(strsplit(df$Species_Pair, "_"), function(x) min(x))
df$group2 = sapply(strsplit(df$Species_Pair, "_"), function(x) max(x))
df = df %>%
  mutate(group1 = parse_character(group1), group2 = parse_character(group2))
df$sorted_comparison = apply(df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "_"))
df.long = df %>% 
  select(sorted_comparison, everything(), 
         -Species_Pair, -group1, -group2) %>%
  unique() %>%
  mutate(Community = as.numeric(Community)) %>%
  left_join(community.legend, by = c("Community" = "community")) %>%
  separate(sorted_comparison, into = c("s1", "s2"), sep = "_", remove = F) %>%
  filter(s1 != "ZOOP",
         s2 != "ZOOP") %>%
  separate(community.name, into = c("Water", "Season")) %>%
  mutate(group1 = case_when(s1 %in% c("LEAF", "PERI", "FISH")~ s1, s1 %nin% c("LEAF", "PERI", "FISH")~"macroinvert"),
         group2 = case_when(s2 %in% c("LEAF", "PERI", "FISH")~ s2, s2 %nin% c("LEAF", "PERI", "FISH")~"macroinvert")) 



## Seasonal overlap shifts for taxa that occur in a water in both seasons
df.long  %>%
  mutate(s1 = str_trim(s1, side = "right"),
         s2 = str_trim(s2, side = "right")) %>%
  ungroup() %>%
  filter(group1 == group2, group1 == "macroinvert") %>%
  group_by(Water, sorted_comparison, post )  %>%
  mutate(n = length(unique(Season))) %>% 
  filter(n > 1) %>%
  ungroup() %>%
  select(sorted_comparison, Water, post, Season, Values) %>%
  distinct() %>% 
  pivot_wider(
    id_cols = c(sorted_comparison, Water, post),  # what stays as rows
    names_from = Season,                # what becomes column names
    values_from = Values
  ) %>%
  group_by(Water, sorted_comparison) %>%
  mutate(fall.mean.water = mean(fall), 
            spring.mean.water = mean(spring)) %>%
  ungroup() %>%
  group_by(sorted_comparison, post) %>%
  mutate(overall.pair.spring = mean(spring),
         overall.pair.fall = mean(fall)) %>%
  ungroup() %>% 
  group_by(sorted_comparison) %>%
  mutate(overall.spring = mean(overall.pair.spring), 
         overall.fall = mean(overall.pair.fall)) %>%
  mutate(sorted_comparison = str_replace(sorted_comparison, "_", " & ")) %>%
  ggplot(aes(x = overall.spring, y = overall.fall)) +
  theme_minimal(base_size = 14) +
  geom_point(aes(col = sorted_comparison),size = 5, shape = "diamond") +
  geom_point(aes(x = spring.mean.water, y = fall.mean.water, col = sorted_comparison)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual( "Family Comparison", values = wes_palette("Darjeeling1",type = "discrete", n = 3)) +
  xlab("Spring Overlap") + ylab("Fall Overlap")




## Combine in the clusters/richness/chemistry -------------------
ca = df.long %>%
  filter(group1 == group2, group1 == "macroinvert") %>%
  group_by(Community, post) %>%
  summarize(mean_overlap = mean(Values)) %>%
  left_join(cluster_chem, by = c("Community" = "community")) 
 

ca %>%
  ungroup() %>%
  group_by(cluster, post) %>%
  summarize(mean_overlap = mean(mean_overlap)) %>%
  filter(!is.na(cluster)) %>%
  ggplot(aes(x = as.factor(cluster), y = mean_overlap, col = as.factor(cluster))) + 
  theme_minimal(base_size = 14) +
  stat_summary(
    fun = mean, 
    fun.min = function(x) quantile(x, 0.025),   # lower bound of 95% CI
    fun.max = function(x) quantile(x, 0.975),   # upper bound of 95% CI
    geom = "pointrange",
    shape = 18, size = 1.2,
  )   + 

  scale_color_manual("Cluster", values = wes_palette("Darjeeling1", type = "discrete", n =3)) +
  ylab("Overlap (95% CI)") + xlab("Cluster") +
  theme(legend.position = "none")

## need Bayestest comparison between clusters

ca %>% ggplot(aes(x = richness, y = mean_overlap)) + 
  geom_point() + 
  geom_smooth()



## for loop to see what is significantly associated with mean overlap

overlap.lm.data = ca %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4) %>%
  select(cluster, everything()) %>%
  ungroup() %>%
  pivot_longer(7:23, names_to = "metric", values_to = "values")  %>%
  filter(is.na(metric) ==F)

### Nothing?
for(i in 1:length(unique(overlap.lm.data$metric))){
  lm.data = overlap.lm.data %>% 
    filter(metric == unique(overlap.lm.data$metric)[i]) %>%## Remove after this line to not have the mean value by community
    group_by(Community, values) %>%
    summarize(mean_overlap = mean(mean_overlap))
  
  lm_summary = lm(data = lm.data, mean_overlap~values) %>% summary()
  
  if(lm_summary$coefficients[2,4] < .05){
    print(unique(overlap.lm.data$metric)[i])
    print(lm_summary$r.squared)
    g = lm.data  %>%
      ggplot(aes(x = values, y = mean_overlap)) + geom_point() +
      geom_smooth(method = "lm")
    print(g)
    
  }else{
    NULL
  }
  
}


## Niche Area --------------------------------
## Ellipse area just use full posterior with full length of species/ellipses  
ellipse.area = siberEllipses(posterior) %>% 
  as.data.frame() %>%
  rename_with(~ names(posterior)) %>% 
  mutate(post_n = seq(1:length(.[,1]))) %>%
  pivot_longer(1:length(posterior), 
               names_to = "comm",
               values_to = "area") %>%
  separate(comm, into = c("community", "group")) %>%
  group_by(post_n, community) %>%
  mutate(total_area = sum(area)) %>%
  mutate(relative_area = area / total_area) %>% 
  group_by(community, group)  %>%
  mutate(community = as.numeric(community),
         group = as.numeric(group)) %>%
  left_join(community.legend) %>%
  left_join(legend)

axis_order = c("aeshnidae","gomphidae", "libellulidae", "corduliidae", "coenagrionidae", "heptageniidae", "leptophlebiidae", "ephemerellidae", "notonectidae", "viviparidae", "lymnaeidae","phryganeidae","limnephilidae", "hydropsychidae", "chironomidae", "asellidae", "talitridae")

order_colors = c("#904E55", "#94A187","#29524A","#3AAFB9", "#C5AFA0", "#6F1A07", "#875C74","#322E18")
order_colors = wes_palette("Darjeeling1", type = "continuous", n = 8)
order_levels = c("odonata", "ephemeroptera", "hemiptera", "gastropoda", "trichoptera", "diptera", "isopoda", "amphipoda")
## Regular figure of niche area and order

ellipse.area %>%
  mutate(community = as.numeric(community),
         group = as.numeric(group)) %>%
  left_join(community.legend) %>%
  left_join(legend) %>%
  left_join(taxon_frame %>%
              select( FAMILY, ORDER) %>%
              unique(), by = c("group.name" = "FAMILY")) %>%
  filter(group.name %nin% c("FISH", "LEAF","PERI","ZOOP", "INSECT")) %>%
  group_by(post_n,  group, group.name, ORDER) %>%
  summarize(area = mean(area)) %>%
 
  ggplot(aes(x = factor(group.name, levels = axis_order), 
             y = area,
             fill = factor(ORDER, levels = order_levels))) + 
  geom_boxplot(outliers = FALSE) +
  
  theme_minimal(base_size = 14) + 
  #theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("SEAc 95% CI") +
  scale_fill_manual(values = order_colors) +
  scale_y_log10() +
  labs(fill = "Order")

## Seasonal area effect
seasonal.means = ellipse.area %>% 
  ungroup() %>%
  select(-total_area, -relative_area, -community) %>%
  separate(community.name, into = c("Water", "Season")) %>%
  group_by(group, Water) %>%
  mutate(season.distinct = length(unique((Season)))) %>%
  filter(season.distinct > 1) %>%
  select(-season.distinct) %>%
  unique() %>%
  group_by(common,  Season, post_n) %>%
  summarize(area = mean(area)) %>%
  ungroup() %>% 
  group_by(common, Season) %>%
  summarize(area = mean(area))%>%
  pivot_wider(values_from = "area", names_from = "Season") 

seasonal.waters = ellipse.area %>% 
  ungroup() %>%
  select(-total_area, -relative_area, -community) %>%
  separate(community.name, into = c("Water", "Season")) %>%
  group_by(group, Water) %>%
  mutate(season.distinct = length(unique((Season)))) %>%
  filter(season.distinct > 1) %>%
  select(-season.distinct) %>%
  unique() %>%
  group_by(common, Water, Season) %>%
  summarize(area = mean(area)) %>%
  pivot_wider(values_from = "area", names_from = "Season") %>%
  select(spring, fall, everything()) 


ggplot() + 
  theme_minimal(base_size = 14) +
  geom_point(data = seasonal.waters, aes(x = spring, y = fall, color = common), key_glyph = "rect") + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") + 
  geom_point(aes(x = seasonal.means$spring, y = seasonal.means$fall, color = seasonal.means$common), size = 5, shape = "diamond") +
  labs(col = "Family") + xlab("Spring Niche Area") + ylab("Fall Niche Area") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 4))



## Average seasonal overlap per lake change

ellipse.area %>% 
  separate(community.name, into = c("Water", "Season")) %>% 
  filter(Water %nin% c("LML", "COM","ETL","GNL")) %>%
  group_by(Water, Season, post_n) %>% 
  summarise(mean_area = mean(area)) %>%
  ungroup() %>%
  ggplot(aes(x = Season, y = mean_area)) +
  theme_minimal(base_size = 14) + 
   stat_summary(
  fun = mean, 
  fun.min = function(x) quantile(x, 0.025),   # lower bound of 95% CI
  fun.max = function(x) quantile(x, 0.975),   # upper bound of 95% CI
  geom = "pointrange",
  shape = 18, size = 1.2,
  ) +
  facet_wrap(~Water, scales = "free")


## Does standard deviation of niche area vary between clusters?

area.cluster = ellipse.area %>%
  mutate(community = as.numeric(community),
         group = as.numeric(group)) %>%
  left_join(community.legend) %>%
  left_join(legend) %>%
  filter(group.name %nin% c("FISH", "LEAF","PERI","ZOOP", "INSECT")) %>%
  left_join(taxon_frame %>% 
              ungroup() %>%
              select(TAXON, FAMILY) %>%
              unique(), by = c("group.name" = "TAXON")) %>%

  
  
  filter(FAMILY %in% c("aeshnidae", "coruliidae", "coenagrionidae", "libellulidae") )%>%
  ungroup() %>%
  select(community, community.name, group.name, post_n, area) %>%
  ungroup() %>%
  left_join(richness %>% select(ID, cluster) %>% unique, by = c("community.name" = "ID")) %>%
  group_by(cluster, post_n) %>%
  ## try mean of lakes first then clusters
  ungroup() %>%
  group_by(community, cluster, post_n) %>%
  summarize(area = mean(area)) %>% 
  ungroup() %>%
  group_by(cluster,post_n) %>%
  summarize(area = mean(area)) %>%
  filter(!is.na(cluster))


## Graph
area.cluster %>%
  ggplot(aes(x = as.factor(cluster), 
             y = area, 
             col = as.factor(cluster))) + 

  stat_summary(
    fun = mean, 
    fun.min = function(x) quantile(x, 0.025),   # lower bound of 95% CI
    fun.max = function(x) quantile(x, 0.975),   # upper bound of 95% CI
    geom = "pointrange",
    shape = 18, size = 1.5,
   ) +
  theme_minimal(base_size = 16) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Niche area (SEAc)") +
  #scale_y_log10() +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  scale_x_discrete(labels = c("1" = "Deep\nthermocline",
                              "2" = "Intermediate\nthermocline", 
                              "3" = "Shallow\nthermocine"))
## Practical equivalence between clusters/area

area1 = area.cluster %>%
  filter(cluster == 1)
area2 = area.cluster %>%
  filter(cluster == 2)
area3 = area.cluster %>% 
  filter(cluster ==3 )

pd12 = mean(area1$area > area2$area)

pd13 = mean(area1$area > area3$area)

pd23 = mean(area2$area > area3$area)

library(bayestestR)
# If you want to define ROPE as 10% of the SD of the outcome
y_sd = sd(area.cluster$area)
rope_range = c(-0.1 * y_sd, 0.1 * y_sd)
rope(area2$area - area3$area, range = rope_range) ## < 5% 

## sig difference

area.cluster %>% pivot_wider(names_from = cluster, values_from = sd_area) %>%
  mutate(dif1_2 =  `2`-`1`,
         dif1_3 = `1`-`3`) %>%
  ungroup() %>%
  summarize(mean1_2 = mean(dif1_2), quant.low1_2 = quantile(dif1_2, .025), quant.high1_2 = quantile(dif1_2, .975),
            mean1_3 = mean(dif1_3), quant.low1_3 = quantile(dif1_3, .025), quant.high1_3 = quantile(dif1_3, .975))





## How does niche are vary overa




## Join area with the data out of richness to get chemistry for each of the waters
area_summary = ellipse.area %>% group_by(community,community.name, group.name, group) %>%
  summarize(mean_area = mean(area)) %>%
  left_join(cluster_chem %>% 
              select(-community), by = "community.name" )


## For loop for overlap and metric - automates the regressions across variables

## The modern DOC isn't significant

area.loop = area_summary %>% 
  filter(group.name %in% c("gomphidae", 
                           "heptageniidae",
                           "coenagrionidae"
                           ))  %>%
  select(-SurficialGeology, -cluster, -Pond_num, -ID, -Lake, -Lake.Type,
         -DOC_update_text, -temp_do) %>%

  ungroup() %>%
  pivot_longer(7:17, names_to = "metric", values_to = "values")  %>%
  filter(is.na(metric) ==F)

## Only DOC appears to be significant
for(i in 1:length(unique(area.loop$metric))){
  #i = 20
  lm.area = area.loop %>% 
    filter(metric == unique(area.loop$metric)[i]) 
  
  lm.summary.area = lm(data = lm.area, mean_area ~ values) %>% summary()
  
  if(lm.summary.area$coefficients[2,4] < .05){
    print(unique(lm.area$metric))
    print(lm.summary.area$r.squared)
    g = lm.area  %>%
      ggplot(aes(x = values, y = mean_area)) + 
      geom_jitter(aes(col = group.name), size = 3, alpha = .5) +
      geom_smooth(method = "lm", col="black") +
      theme_minimal(base_size = 20) +
      xlab(paste(lm.area$metric %>% unique)) +
      ylab("Niche Area") +
      labs(col = "Family") +
      xlab("DOC") +
      scale_color_manual(values = wes_palette("Darjeeling1")) 
      
    
    print(g)
    
  }else{
    NULL
  }
}


## I think DOC matters I'm going to try to do a random effects model


lmer.data = area.loop %>% select(mean_area, metric, values, community.name, group.name) %>%
  separate(community.name, into = c("Water", "Season")) %>%
  filter(metric == "thermo_depth") # thermo_depth # DOC_update
lmer.data %>% print(n = 100
                    )
library(lme4)
lmer(mean_area ~ values +  Season +  (1| Water), data = lmer.data) %>% summary()

lmer(mean_area ~ Season + (1|Water), data = lmer.data ) 











area.loop %>% filter(metric %in% c("DOC_update", "DOC.1")) %>%
 # filter(values < 8) %>%
  #filter(community.name %nin% c("ETL.fall", "GNL.fall", "COM.fall")) %>%
  ggplot() + 
  geom_jitter(aes(x = values, y = mean_area, col = metric)) + 
  geom_smooth(aes(x = values, y = mean_area, col = metric),method = lm)


modern.doc = area.loop %>% filter(metric %in% c("DOC_update")) 
cor.test(modern.doc$values, modern.doc$mean_area, method = "spearman")
lm(data = modern.doc, mean_area ~ log10(values)) %>% summary()


## Create the DOC graph - the only significant variable...
area_summary %>% filter(group.name %in% c("gomphidae", "heptageniidae", "coenagrionidae")) %>%
  lm(data = ., mean_area ~ DOC_update) %>% summary()

load(file = "Data/RData/taxa.clusters.RData")
area_summary %>%
  ungroup() %>%
  select(group.name, mean_area) %>%
 # group_by( group.name) %>%
# summarize(mean_area = mean(mean_area)) %>% 
  filter(group.name %nin% c("FISH", "INSECT", "LEAF","ZOOP")) %>%
  left_join(taxa.assign, by = c("group.name" = "FAMILY")) %>%
  
  ggplot(aes(x = cluster, y = mean_area)) + 
  geom_boxplot() + 
  geom_point()
## FFG --------------

ffg = read.csv("Data/CSVs/FFGS.csv") %>%
  select(-GENUS) %>%
  unique()

### FFG + Area ------
ffg.summary = area_summary %>% 
  left_join(ffg, by = c("group.name"= "FAMILY")) %>%
  ungroup() %>%
  select(FG,  everything()) %>%
  mutate(FG = case_when(group.name %in% c("LEAF", "PERI","FISH") ~ group.name,
                            group.name %nin% c("LEAF", "PERI", "FISH") ~ FG))  %>% 
  mutate(FG = as.character(FG)) %>%
  as.data.frame() %>%
  select(mean_area, FG) %>%
  na.omit() %>%
  filter(FG %nin% c("collector-filterer","shredder")) ## Respectively only have 1 and 2 points

ffg.summary %>% 
  ggplot(aes(x = FG, y = mean_area)) +
  geom_boxplot(outliers = F) + 
  geom_point() +
  theme_minimal() +
  scale_x_discrete(labels = c("collector\ngatherer", "fish", "leaf","periphyton","predator","scraper","shredder"))

## ANOVA comparing functional feeding group to mean area... the anova is singificant but only differences are between periphyton and something else... not that interesting
model.summary = aov(data = ffg.summary, mean_area ~ FG) 
(model.summary %>% summary() %>% unlist())[9] ## P value

TukeyHSD(model.summary)

## FFG and overlap

ffg.overlap = group.overlap %>% 
  filter(group1 == "macroinvert") %>%
  left_join(ffg, by = c( "s1" = "FAMILY")) %>%
  rename("FG1" = "FG") %>%
  left_join(ffg, by = c( "s2" = "FAMILY")) %>%
  rename("FG2" = "FG") %>%
  na.omit() %>%
  unite("ID", c(FG1, FG2), remove = F) %>%
  mutate(ID = str_replace(ID, "collector-gatherer_predator", "predator_collector-gatherer"),
         ID = str_replace(ID, "collector-gatherer_scraper", "scraper_collector-gatherer"),
         ID = str_replace(ID, "collector-gatherer_shredder", "shredder_collector-gatherer"),
         ID = str_replace(ID, "predator_scraper", "scraper_predator"))

ffg.overlap.order = c("collector-gatherer_collector-gatherer", "predator_predator", "scraper_scraper", "predator_collector-filterer", "predator_collector-gatherer", "predator_shredder", "scraper-collector-filterer", "scraper_collector-gatherer", "scraper_predator", "scraper_shredder", "shredder_collector-gatherer")
## FFG overlap with rare comparisons
ffg.overlap %>%
  ungroup() %>%
  group_by(post,ID) %>%
  summarize(Values = mean(Values)) %>%
  ggplot(aes(x = factor(ID, ffg.overlap.order), y = Values)) + 
  geom_boxplot() +
  #geom_jitter(alpha = .4) +
  scale_x_discrete(labels =  c("CG\nCG","P\nP", "Sc\nSc","P\nCF", "P\nCG",  "P\nSh",
                               "Sc\nCF","Sc\nCG","Sc\nP", 
                               "Sc\nSh","Sh\nCG")) +
  theme_minimal() +
  ylab("Pairwise Overlap") +
  theme(axis.title.x = element_blank())


## FFG overlap without rare comparisons
ffg.overlap %>%
  ungroup() %>%
  group_by(Community,ID) %>%
  summarize(Values = mean(Values)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  ggplot(aes(x = ID, y = Values)) + 
  geom_boxplot() +
  geom_point() +
  scale_x_discrete(labels =  c("CG\nCG",  "P\nCG", "P\nP",
                               "Sc\nCG","P\nSc", "Sc\nSc")) +
  theme_minimal() +
  ylab("Pairwise Overlap")

ffg.aov.data = ffg.overlap %>%
  ungroup() %>%
  group_by(Community,ID) %>%
  summarize(Values = mean(Values)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  filter(n > 2)

lm_model = lm(Values ~ ID, data = ffg.aov.data)


summary(lm_model)

ffg.overlap.model  = aov(data = ffg.aov.data, Values ~ ID)
tukey.ffg.overlap = TukeyHSD(ffg.overlap.model)
tukey.ffg.overlap$ID %>% as.data.frame() %>%
  filter(`p adj` < .05)
### within and between FFGs ---------


comp = ffg.overlap %>%
  ungroup() %>%
  group_by(Community,ID, FG1, FG2) %>%
  summarize(Values = mean(Values)) %>%
  mutate(comp = case_when(FG1 == FG2 ~ "same", FG1 != FG2 ~ "diff")) 

comp %>%
  ggplot(aes(x = comp, y =Values)) + 
  geom_boxplot(fill = NA) + 
  geom_point() +
  theme_minimal() + 
  ylab("Pairwise Overlap") +
  xlab("FFG") +
  scale_x_discrete(labels = c("Different", "Same"))

t.test((comp %>% filter(comp == "same"))$Values, (comp %>% filter(comp == "diff"))$Values)

## ANOVA with just the groups that have different FFGs
comp.aov =comp %>% filter(comp == "diff") %>%
  group_by(ID) %>%
  mutate(n = n()) %>%
  filter(n > 2) 
comp.aov %>%
  ggplot(aes(x = ID, y = Values)) + 
  geom_boxplot() +
  geom_point()
comp.aov.model = aov(data = comp.aov, Values ~ ID) 
TukeyHSD(comp.aov.model)



