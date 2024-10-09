library(SIBER)
library(tidyverse)
library(wesanderson)
`%nin%` = Negate("%in%")
set.seed(123)
source("../LML_SMB_removal/Function_Source_Files/isotope_functions_update.R")
## SIBER ------------

## Load in taxon data


taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Load in isotope measurement file
data.iso=  read.csv("Data/CSVs/iso_data_clean.csv")

## Read in isotope sample file
sample = read.csv("Data/CSVs/isotope_sample.csv")

## Join together -----------------



data = data.iso %>% left_join(sample, by = "ISO_YSAMP_N") %>%
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
  arrange(group.name) %>%
  mutate(group = as.numeric(as.factor(group.name))) %>%
  arrange(community, group) %>%
  as.data.frame() %>%
  select(iso1, iso2, group, community, group.name, community.name)

save(data, file = "Data/RData/isotope_data.RData")

siber.data = data %>% select(-group.name, -community.name)
legend = data %>% select(group.name,group) %>%
  unique() %>%
  mutate(color = c(1:21),
         common = group.name, 
         CODE = group.name)

community.legend = data %>% select(community, community.name) %>%
  unique()

legend


# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

## Run siber
siber.object = createSiberObject(siber.data)

posterior = siberMVN(siber.object, parms, priors)
  




overlap_list = list() 

communities = unique(siber.data$community)[c(-10, -15,-20)]

for(h in 1:length(communities)){
  
  overlap_list[[h]] = overlap(siber.data, communities[h],10, posterior, legend) %>% as.data.frame()
}

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
  mutate(Values = as.numeric(Values)) 





#library(tidyr)
#library(tidyverse)


## table 

df = overlap.long%>% ungroup() 
df$group1 = sapply(strsplit(df$Species_Pair, "-"), function(x) min(x))
df$group2 <- sapply(strsplit(df$Species_Pair, "-"), function(x) max(x))

df = df %>%
  mutate(group1 = parse_character(group1), group2 = parse_character(group2))
df$sorted_comparison <- apply(df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "_"))


df.long = df %>% 
  select(sorted_comparison, everything(), 
         -Species_Pair, -group1, -group2) %>%
  unique() %>%
  mutate(Community = as.numeric(Community)) %>%
  left_join(community.legend, by = c("Community" = "community")) %>%
  separate(sorted_comparison, into = c("s1", "s2"), sep = "_", remove = F)


## Looking at the overlap
df.long %>%
  group_by(Community, sorted_comparison) %>%
  summarize(mean_overlap = mean(Values)) %>%
  ggplot(aes(x = as.factor(Community), y = mean_overlap)) + 
  geom_boxplot()

## I want overlap to be able to be examined within/across meaningful groups

group.overlap = df.long %>%
  mutate(group1 = case_when(s1 %in% c("LEAF", "PERI", "FISH")~ s1, s1 %nin% c("LEAF", "PERI", "FISH")~"macroinvert"),
         group2 = case_when(s2 %in% c("LEAF", "PERI", "FISH")~ s2, s2 %nin% c("LEAF", "PERI", "FISH")~"macroinvert")) 

group.overlap %>%
  filter(group1 == group2,
         group1 == "macroinvert") %>%
  group_by(Community, post) %>%
  summarize(mean_overlap = mean(Values)) %>%
  ggplot(aes(x = as.factor(Community), y = mean_overlap)) + 
  geom_boxplot()


axis_order = c("aeshnidae","gomphidae", "libellulidae", "corduliidae", "coenagrionidae", "heptageniidae", "leptophlebiidae", "ephemerellidae", "notonectidae", "viviparidae", "lymnaeidae","phryganeidae","limnephilidae", "hydropsychidae", "chironomidae", "asellidae", "talitridae", "ZOOP")

group.overlap  %>%
  filter(group1 == group2, group1 == "macroinvert") %>%
  group_by(Community, sorted_comparison, s1, s2) %>%
  summarize(mean_overlap = mean(Values)) %>%
  ggplot(aes(x = factor(s1, levels = axis_order), 
             y = factor(s2, levels = rev(axis_order)), 
             col = mean_overlap)) +
  geom_point(size = 3) +
  theme_minimal() + 
  scale_color_viridis_b() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        axis.title.y = element_blank()) +
  xlab("Family") +
  labs(col = "Average Overlap")


## Combine in the richness and chemistry -------------------
load("Data/RData/richness_cluster.RData")  
richness = richness %>%
  unite("ID", WATER, season, sep = ".") %>%
  left_join(community.legend, by = c("ID" = "community.name"))


ca = group.overlap %>%
  filter(group1 == group2, group1 == "macroinvert") %>%
  group_by(Community, post) %>%
  summarize(mean_overlap = mean(Values)) %>%
  left_join(richness, by = c("Community" = "community"))
 

ca %>% ggplot(aes(x = cluster, y = mean_overlap)) + 
  geom_boxplot() +
  geom_point() 

overlap.aov = aov(data = ca, mean_overlap ~ cluster) 
TukeyHSD(overlap.aov)

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

## Most significant ones with R2 values greater than .01 are sechi depth, Ca, and DOC. 


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
ellipse.area %>%
  mutate(community = as.numeric(community),
         group = as.numeric(group)) %>%
  left_join(community.legend) %>%
  left_join(legend) %>%
  left_join(taxon_frame, by = c("group.name" = "FAMILY")) %>%
  filter(group.name %nin% c("FISH", "LEAF","PERI","ZOOP")) %>%
  group_by(post_n,  group, group.name, ORDER) %>%
  summarize(area = mean(area)) %>%
  ggplot(aes(x = factor(group.name, levels = axis_order), 
             y = area,
             fill = factor(ORDER, levels = order_levels))) + 
  geom_boxplot(outliers = FALSE) +
  
  theme_minimal() + 
  #theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Niche Area") +
  scale_fill_manual(values = order_colors) +
  scale_y_log10() +
  labs(fill = "Order")


ellipse.area %>%
  ggplot(aes(x = community.name, y = area, fill = community.name)) + 
  geom_boxplot(outliers = FALSE) +
  theme_minimal(base_size = 5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  facet_wrap(~group.name, scales = "free")

ellipse.area %>%
  mutate(community = as.numeric(community),
         group = as.numeric(group)) %>%
  left_join(community.legend) %>%
  left_join(legend) %>%
  ggplot(aes(x = group.name, y = area, fill = community.name)) + 
  geom_boxplot(outliers = FALSE) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  facet_wrap(~community.name, scales = "free")



## Join area with the data out of richness to get chemistry for each of the waters
area_summary = ellipse.area %>% group_by(community,community.name, group.name, group) %>%
  summarize(mean_area = mean(area)) %>%
  left_join(richness)


## For loop for overlap and metric - automates the regressions across variables

area.loop = area_summary %>% 
  filter(group.name %in% c("gomphidae", 
                           "heptageniidae",
                           "coenagrionidae"))  %>%
  select(-SurficialGeology, -DOC, -temp_do, 
         -AirEqPh, -DIC, -F, -Fe, -K, -LabPh,
         -Mg, -Mn, -Na, -Pb, -rate, -SCONDUCT, 
         -Tal, -TDAI, -TotalP2, -TrueColor,
         -Volume, -Zn, -CL, -Lake.Type, -Lake,
         -ID, -min_do, -SO4) %>%
  select(cluster, everything()) %>%
  ungroup() %>%
  pivot_longer(7:25, names_to = "metric", values_to = "values")  %>%
  filter(is.na(metric) ==F)

## Only DOC appears to be significant
for(i in 1:length(unique(area.loop$metric))){
  
  lm.area = area.loop %>% 
    filter(metric == unique(area.loop$metric)[i]) 
  
  lm.summary.area = lm(data = lm.area, mean_area ~ values) %>% summary()
  
  if(lm.summary.area$coefficients[2,4] < .05){
    print(unique(lm.area$metric))
    print(lm.summary.area$r.squared)
    g = lm.area  %>%
      ggplot(aes(x = values, y = mean_area)) + 
      geom_point(aes(col = group.name)) +
      geom_smooth(method = "lm", col="black") +
      theme_minimal() +
      xlab(paste(lm.area$metric %>% unique)) +
      ylab("Niche Area") +
      labs(col = "Family") +
      scale_color_manual(values = wes_palette("Darjeeling1")) 
      
    print(g)
    
  }else{
    NULL
  }
}

## Create the DOC graph - the only significant variable...
area_summary %>% filter(group.name %in% c("gomphidae", "heptageniidae", "coenagrionidae")) %>%
  lm(data = ., mean_area ~ DOC.1) %>% summary()

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

lm_model <- lm(Values ~ ID, data = ffg.aov.data)


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



