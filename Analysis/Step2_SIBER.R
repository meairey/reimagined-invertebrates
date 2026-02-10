# Isotopic Niches
## Script setup -------------
`%nin%` = Negate("%in%")
set.seed(123)

## Libraries
library(SIBER)
library(tidyverse)
library(wesanderson)
library(bayestestR)

## Custom functions for ellipse metric generation and plotting
source("../solid-fishstick/Analysis/Documentation/isotope_functions_update.R") 


## Load in taxon data
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

data.iso = read.csv(file = "Data/CSVs/processed_data.csv")

data = read.csv(file = "Data/CSVs/siber_data.csv")

## SIBER formatted frame 
data = data.iso %>% ## This one is formatted for SIBER
  left_join(taxon_frame %>% select(TAXON, FAMILY, ORDER) %>% unique()) %>% ## Taxomonic metadata join
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


## Legends for plotting
legend = data %>% select(group.name,group) %>% ## Legend for each family/group name
  unique() %>%
  mutate(color = c(1:17),
         common = group.name, 
         CODE = group.name)

community.legend = data %>% select(community, community.name) %>%
  unique() ## Legend for the lake names/codes

# save(community.legend, file = "Data/RData/community.legend.RData")

## Cluster and chemistry data -- cluster from `Step1_NMDS.R`
cluster_chem = read.csv(file = "Data/CSVs/richness_update.csv") %>%
  #unite("ID", WATER, season, sep = ".") %>%
  left_join(community.legend, by = c("ID" = "community.name"))



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
  

## Overlap analysis --------------------
load( file = "Data/RData/overlap_list.RData")

overlap_list = list() 

communities_overlap = unique(data$community)[-9] ## For overlap remove communities with only 1 taxa

for(h in 1:length(communities_overlap)){
  
  overlap_list[[h]] = overlap(data, communities_overlap[h],50, posterior) %>% as.data.frame()
}
#save(overlap_list, file = "Data/RData/overlap_list.RData") ## Save this to avoid running each time

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
overlap_sorted_comparison = overlap.long %>% 
  ungroup() %>% 
  mutate(group1 = sapply(strsplit(.$Species_Pair, "_"), function(x) min(x)),
         group2 = sapply(strsplit(.$Species_Pair, "_"), function(x) max(x))) %>%
  mutate(group1 = parse_character(group1), group2 = parse_character(group2), 
         sorted_comparison = apply(.[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = "_")))  %>% 
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

## Combine in the clusters/richness/chemistry -------------------
ca = overlap_sorted_comparison %>%
  filter(group1 == group2, group1 == "macroinvert") %>%
  left_join(community.legend, by = c("Community" = "community"))# %>% ## stopped here to go below nad try lmer by sorted comparisin
  #group_by(community.name, Community, post) %>% ## grouping by community
  #summarize(mean_overlap = mean(Values)) %>%
  left_join(cluster_chem, by = c("community.name"))
 
#### Figure 5A - Cluster and overlap comparison ------
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
  scale_color_manual("Cluster", values = wes_palette("Royal2")[c(3,1,5)]) +
  ylab("Overlap (95% CI)") + xlab("Cluster") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.x = element_blank())  +
  scale_x_discrete(labels = c("1" = "Deep\nthermocline",
                              "2" = "Intermediate\nthermocline", 
                              "3" = "Shallow\nthermocine")) -> Figure5A


## for loop to see what is significantly associated with mean overlap

overlap.lm.data = ca %>% 
  group_by(sorted_comparison, Community, Season, community.name) %>%
  summarize(mean = mean(Values))%>%
  left_join(cluster_chem, by = c("community.name")) %>%

  select(
    -SurficialGeology,  
      -Pond_num,
    -DOC_update_text, -temp_do,
          -Lake.Type, 
          -Lake,
         -ID, 
    -min_do
    ) %>%
  ungroup() %>%
    mutate(across(max_depth:sechi.depth,
         ~ as.numeric(scale(.x)))) %>%
  group_by(sorted_comparison) %>%
  mutate(n = n()) %>% select(n, everything()) %>%
  ungroup() %>%
  filter(n > 2)
  


vars = colnames(overlap.lm.data)[7:16]
vars
for(i in 1:length(vars)){
  i = 5
  form = as.formula(paste0("mean ~ ", vars[i], " + (1|sorted_comparison)"))
  lmer.run = lmerTest::lmer(form, data = overlap.lm.data) %>%
  summary()
  if(lmer.run$coefficients[2,5] < .05){
    print(vars[i])
    print(lmer.run$coefficients[2,])
  }
}
  
lmer.area.DOC = lmerTest::lmer(mean.area ~ DOC_update + (1|group.name), data = lmer.area.data)
lmer.area.DOC %>% summary()
emtrends(lmer.area.DOC, 
                       var = "DOC_update")



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
 
  ggplot(aes(y = factor(group.name, levels = rev(axis_order)), 
             x = area,
            # color = factor(ORDER, levels = order_levels)
            )) + 
     stat_summary(
  fun = mean, 
  fun.min = function(x) quantile(x, 0.025),   # lower bound of 95% CI
  fun.max = function(x) quantile(x, 0.975),   # upper bound of 95% CI
  geom = "pointrange",
  shape = 18, size = 1.2,
  ) +
  
  theme_minimal(base_size = 14) + 
  #theme(legend.position = "none") + 
  theme(#axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.y = element_blank(), 
        legend.position = "bottom", 
        legend.box = "vertical",
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.04, "cm")) +
  xlab("SEAc 95% CI") +
  #scale_color_manual(values = order_colors) +
  scale_x_log10() +
 # labs(color = "Order") +
  guides(color = guide_legend(nrow = 4, byrow = TRUE)) -> Figure4A
Figure4A
ggsave(file = "Graphics/Figures/Figure_4A.pdf", plot = Figure4A, width = 4.5, height = 5, units = "in")


## LMER on the area data ----

lmer.area.data = ellipse.area %>%
  group_by(community, group, group.name) %>%
  summarize(mean.area = mean(area)) %>%
  left_join(community.legend) %>%
  left_join(cluster_chem, by = "community.name") %>%
  separate(community.name, into = c("WATER", "season")) %>%
  ungroup()  %>%
  select(-ID,-temp_do, -DOC_update_text, -SurficialGeology, -Lake.Type, -min_do) %>% 
  group_by(group.name) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  select(-n) %>% 
  filter(!(community.x == 14 & group.name == "corduliidae")) ## Point removed using Cook's distance below. Check this if rerunning
 
vars = colnames(lmer.area.data)[10:19] ## Variables to run the LMER through


## LMER loop
for(i in 1:length(vars)){
  form = as.formula(paste0("mean.area ~ ", vars[i], " + (1|group.name)"))
  lmer.run = lmerTest::lmer(form, data = lmer.area.data) %>% ## DOC is the only variable that is significant
  summary()
  if(lmer.run$coefficients[2,5] < .05){
    print(vars[i])
    print(lmer.run$coefficients[2,])
  }
}
## Running for the significant DOC
lmer.area.DOC = lmerTest::lmer(mean.area ~ DOC_update + (1|group.name), data = lmer.area.data[-31,])
lmer.area.DOC %>% summary()
## Plot of LMER results
newdat.DOC.area = data.frame( ## Format a new set up data
  DOC_update = seq(
    min(lmer.area.data$DOC_update, na.rm = TRUE),
    max(lmer.area.data$DOC_update, na.rm = TRUE),
    length.out = 100
  ),
  group.name = NA
)
pred.DOC.area = predict( ## Format the predicted data from newdat.DOC.area
  lmer.area.DOC,
  newdata = newdat.DOC.area,
  re.form = NA,
  se.fit = TRUE
)
newdat.DOC.area = newdat.DOC.area %>% 
  mutate(fit = pred.DOC.area$fit,
         lower = pred.DOC.area$fit - 1.96 * pred.DOC.area$se.fit,
         upper = pred.DOC.area$fit + 1.96 * pred.DOC.area$se.fit)
## Visualizing Niche Area + DOC relationship
earthy_colors = c("#8B4513","#6B8E23", "#B66E41",  "#1B1B1B",
                   "#D4A017", "#708090", "#C2B280", "#E97451" )
DOC_Nichesize = ggplot() +
  geom_point(data = lmer.area.data, 
             aes(x = DOC_update, y = mean.area, col = group.name),
             #alpha = 0.7, 
             size = 3) +
  geom_ribbon(data = newdat.DOC.area,
              aes(x = DOC_update, ymin = lower, ymax = upper),
              alpha = 0.25) +
  geom_line(data = newdat.DOC.area,
                aes(x = DOC_update, y = fit),
                linewidth = 1) +
  labs(x = "DOC (mg/L)",
       y = "Niche size (SEAc)") +
  scale_color_manual("Family", values = earthy_colors) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.04, "cm")) +
  guides(color = guide_legend(nrow = 4, byrow = TRUE))
DOC_Nichesize
ggsave(DOC_Nichesize, file = "Graphics/Figures/DOC_NicheSize.pdf", width = 4, height= 5)


## Removing one point from the area LMER because it's an outlier
library(influence.ME)
install.packages("influence.ME")
infl <- influence(lmer.area.DOC, obs = TRUE)
cooks <- cooks.distance(infl)
plot(cooks, type = "h")
abline(h = 4 / length(cooks), col = "red", lty = 2)
lmer.area.data[31,] ## Looks like a corduliidae from community 14 has a pretty high SEAc. Cooks distance suggests leaving it out. Rerunning the lmer above with removed point




## Area by cluster--------------
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
  # filter(FAMILY %in% c("aeshnidae", "coruliidae", "coenagrionidae", "libellulidae") )%>%
  ungroup() %>%
  select(community, community.name, group.name, post_n, area) %>%
  ungroup() %>%
  left_join(cluster_chem %>% select(community.name, cluster) %>% unique, by = c("community.name" )) %>%
  group_by(cluster, post_n) %>%
  ## try mean of lakes first then clusters
  ungroup() %>%
  group_by(community, cluster, post_n) %>%
  summarize(area = mean(area)) %>% 
  ungroup() %>%
  group_by(cluster,post_n) %>%
  summarize(area = mean(area)) %>%
  filter(!is.na(cluster))


####  Figure 5B Area Cluster --------
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
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("SEAc 95% CI") +
  #scale_y_log10() +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) +
  scale_x_discrete(labels = c("1" = "Deep\nthermocline",
                              "2" = "Intermediate\nthermocline", 
                              "3" = "Shallow\nthermocine")) -> Figure_5B



Figure5 = grid.arrange(Figure5A, Figure_5B, ncol = 2)
ggsave(plot = Figure5, file = "Graphics/Figures/Figure5_AreaCluster.pdf", width = 5, height = 3, units = "in")

## Practical equivalence between clusters/area -----------------------

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
rope(area1$area - area3$area, range = rope_range) ## < 5% 
cat = rope(area1$area - area3$area, range = rope_range) 

## Summary table for area cluster comparisons
data.frame(comparison = c("1_2", "2_3","1_3"),
           PD = c(pd12, pd23, pd13),
           ROPE = c(rope(area1$area - area2$area, range = rope_range)$ROPE_Percentage * 100,
                    rope(area2$area - area3$area, range = rope_range)$ROPE_Percentage * 100,
                    rope(area1$area - area3$area, range = rope_range)$ROPE_Percentage * 100))


## Join area with the data out of richness to get chemistry for each of the waters
area_summary = ellipse.area %>% group_by(community,community.name, group.name, group) %>%
  summarize(mean_area = mean(area)) %>%
  left_join(cluster_chem %>% 
              select(-community), by = "community.name" )


## Seasonal Effects ---------------------------------------
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
  pivot_wider(values_from = "area", names_from = "Season") %>%
  mutate(seasonal.diff = spring - fall)

area.mean = mean(seasonal.means$seasonal.diff)

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
  select(spring, fall, everything()) %>%
  mutate(seasonal.diff = spring -fall)



ggplot() + 
  theme_minimal(base_size = 14) +
  geom_point(data = seasonal.waters, aes(x = spring, y = fall, color = common), key_glyph = "rect") + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") + 
  geom_point(aes(x = seasonal.means$spring, y = seasonal.means$fall, color = seasonal.means$common), size = 5, shape = "diamond") +
  labs(col = "Family") + xlab("Spring SEAc") + ylab("Fall SEAc") +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) -> FigureS4

ggplot() + 
  theme_minimal(base_size = 14) + 
  geom_point(data = seasonal.waters, aes(x = seasonal.diff, y = common)) +
  geom_point(data = seasonal.means, aes(x = seasonal.diff, y = common), size = 6, shape = "diamond") +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = area.mean), lty = "dashed")+
  theme(axis.title.y = element_blank()) + 
  xlab("Seasonal Difference in Niche Area") -> FigureS4




ggsave(file = "Graphics/Figures/Figure_S4.pdf", plot = FigureS4, width = 7, height = 5, units = "in")



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




## Seasons with posteriors

## Thinking about modifying this so that we can see if CI overlap 0 suggesting any meaningful seasonal differences
overlap_sorted_comparison  %>%
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
    id_cols = c(sorted_comparison, Water, post), 
    names_from = Season,               
    values_from = Values
  ) %>%
  group_by(Water, sorted_comparison) %>%
  mutate(difference = spring - fall) %>%
  summarize(mean.diff = mean(difference) ,
            mean.low = quantile(difference, .025), 
            mean.high = quantile(difference, .975)) %>%
ggplot( aes(y = sorted_comparison,
               x = mean.diff,
               xmin = mean.low,
               xmax = mean.high, 
            col = Water)) +
  theme_minimal(base_size = 14) + 
  geom_pointrange(
    position = position_jitter(height = 0.20, width = 0),
    size = 0.6
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") 

%>%
 %>%
  ungroup() %>%
  group_by(sorted_comparison, post) %>%
  mutate(overall.pair.spring = mean(spring),
         overall.pair.fall = mean(fall)) %>%
  ungroup() %>% 
  group_by(sorted_comparison) %>%
  mutate(overall.spring = mean(overall.pair.spring), 
         overall.fall = mean(overall.pair.fall)) %>%
  mutate(sorted_comparison = str_replace(sorted_comparison, "_", " & ")) %>%
  mutate(overall_difference = overall.spring - overall.fall) %>% ## Trying to get a different graph type
  ungroup() %>%
  mutate(overall_mean = mean(overall_difference)) %>%
    ## Start new graph
  ggplot() +
  geom_point(aes(x = overall_difference, y = reorder(sorted_comparison, overall_difference)),size = 5, shape = "diamond") +
  geom_point(aes(x = water_difference, y =reorder(sorted_comparison, overall_difference)), alpha = .05) +
  theme_minimal(base_size = 14) +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = overall_mean), lty = "dashed") +
  xlab("Seasonal Difference in Overlap") + 
  theme(axis.title.y = element_blank()) 














## Seasonal overlap shifts for taxa that occur in a water in both seasons
overlap_sorted_comparison  %>%
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
  mutate(water_difference = spring.mean.water - fall.mean.water) %>%
  ungroup() %>%
  group_by(sorted_comparison, post) %>%
  mutate(overall.pair.spring = mean(spring),
         overall.pair.fall = mean(fall)) %>%
  ungroup() %>% 
  group_by(sorted_comparison) %>%
  mutate(overall.spring = mean(overall.pair.spring), 
         overall.fall = mean(overall.pair.fall)) %>%
  mutate(sorted_comparison = str_replace(sorted_comparison, "_", " & ")) %>%
  mutate(overall_difference = overall.spring - overall.fall) %>% ## Trying to get a different graph type
  ungroup() %>%
  mutate(overall_mean = mean(overall_difference)) %>%
    ## Start new graph
  ggplot() +
  geom_point(aes(x = overall_difference, y = reorder(sorted_comparison, overall_difference)),size = 5, shape = "diamond") +
  geom_point(aes(x = water_difference, y =reorder(sorted_comparison, overall_difference)), alpha = .05) +
  theme_minimal(base_size = 14) +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = overall_mean), lty = "dashed") +
  xlab("Seasonal Difference in Overlap") + 
  theme(axis.title.y = element_blank()) -> FigureS1

 


ggsave(file = "Graphics/Figures/Figure_S1.pdf", plot = FigureS1, width = 7, height = 5, units = "in")

## Food Webs

data %>% 
  separate(community.name, into = c("water", "season")) %>%
 # filter(season == "fall") %>%
  ggplot(aes(x = iso1, y = iso2, col = group.name)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", 
        legend.box = "vertical") +
  geom_point(aes(shape = season), size = 1) +
  stat_ellipse(level = .4) +
  facet_wrap(~water,ncol = 3, scales = "free") +
  scale_color_manual("Family", values = wes_palette("Darjeeling1", type = "continuous", n = 17)) +
  scale_shape_manual("Season", labels = c("Fall", "Spring"), values = c(1,2)) +
      theme() +
  xlab(expression(delta^13*C)) +
  ylab(expression(delta^15*N))  +
  theme(
    text = element_text(size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "lines"),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")) -> Figure_S2_webs

ggsave(Figure_S2_webs, file = "Graphics/Figures/Figure_S2_webs.pdf", width = 5, height = 6, units = "in", dpi = 300) 

## Table 1. Table of lake characteristics --------------------

Table1 = cluster_chem %>% 
  separate(community.name, into = c("water", "season")) %>%
  unique() %>%
  select(water,Elevation, max_depth, surface_area,season, DOC_update,thermo_depth ) %>%
  mutate(water = case_when(water == "COM" ~ "Combs Lake", 
                           water == "DTL" ~ "Dart Lake", 
                           water == "ECL" ~ "East Copperas Lake", 
                           water == "ETL" ~ "East Lake", 
                           water == "GNL" ~ "Green Lake", 
                           water == "HTL" ~ "Heart Lake", 
                           water == "LCL" ~ "Little Clear Lake", 
                           water == "LML" ~ "Little Moose Lake", 
                           water == "MSL" ~ "Moss Lake", 
                           water == "REL" ~ "Rondaxe Lake", 
                           water == "SEL" ~ "Sagemore Lake", 
                           water == "UCL" ~ "Upper Cascade Lake")) %>%
  mutate(season = case_when(season == "spring" ~ "Early", 
                            season == "fall" ~ "Late")) %>%
  arrange(water, season)

write.csv(Table1, file = "Data/CSVs/Table_1_LakeCharacteristics.csv", row.names = F)

## Layman metrics of communities
layman.dat = data.iso %>% ## This one is formatted for SIBER
  left_join(taxon_frame %>% select(TAXON, FAMILY, ORDER) %>% unique()) %>% ## Taxomonic metadata join
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
 # filter(total > 3) %>% # Do not filter out rare families
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
  select(iso1, iso2, group, community, group.name, community.name) %>% ## keep group and community name for legend
  group_by(community, group.name) %>%
  summarize(iso1 = mean(iso1, na.rm =T), 
            iso2 = mean(iso2, na.rm = T)) %>%
  rename(group = group.name) %>%
  select(iso1, iso2, group, community)

siber.layman = createSiberObject(layman.dat)
lmetrics = communityMetricsML(siber.layman) 

lmetrics %>% as.data.frame() %>% 
  rownames_to_column(var = "metric") %>%
  pivot_longer(`1`:`20`, names_to = "community", values_to = "metric.value") %>%
  left_join(community.legend %>% 
              mutate(community = as.character(community))) %>%
  left_join(cluster_chem%>% 
              mutate(community = as.character(community)), by = "community.name") %>%
  separate(community.name, into = c("water", "season")) %>%
  filter(season == "fall") %>%
  ggplot(aes(x = as.factor(cluster), y = metric.value)) + 
  geom_boxplot() + 
  facet_wrap(~metric, scales = "free") + 
  theme_minimal(base_size = 14) + 
  xlab("Cluster")



mu.post = extractPosteriorMeans(siber.object,posterior)

bay.lay = bayesianLayman(mu.post) %>%
  do.call(rbind, .) %>% 
  as.data.frame() %>%
  mutate(community = rep(1:length(mu.post), each = dim(mu.post[[1]])[1])) %>%
  mutate(post = rep(1:4000,  19)) %>%
  select(community, post, everything()) %>%
  pivot_longer(dY_range:SDNND, names_to = "metric", values_to = "metric.value") %>%
  left_join(community.legend ) %>%
  left_join(cluster_chem %>%
              select(-community), by = c("community.name")) %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster, post, metric) %>%
  summarize(metric.value= mean(metric.value)) %>%
  filter(metric != "TA") %>%
   rbind(area.cluster %>%
  rename(post = post_n, metric.value = area) %>%
  mutate(metric = "NicheArea"))


cluster.layman = bay.lay %>%
  mutate(metric = factor(metric, levels = c("dX_range", "dY_range", "CD", "NND", "SDNND", "NicheArea"))) %>%
  ggplot(aes(x = as.factor(cluster), y = metric.value, col = as.factor(cluster))) + 
  stat_summary(
    fun = mean, 
    fun.min = function(x) quantile(x, 0.025),   # lower bound of 95% CI
    fun.max = function(x) quantile(x, 0.975),   # upper bound of 95% CI
    geom = "pointrange",
    shape = 18, size = 1.2,
  )  +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller("metric" = c("dX_range" = "d13C Range",
                                              "dY_range" = "d15N Range", "NicheArea" = "Niche Area"))) + 
  theme_minimal(base_size = 12) + 
  xlab("Cluster") + 
  ylab("95% CI") +
  scale_color_manual(values = wes_palette("Royal2")[c(3,1,5)]) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
       )  +
  scale_x_discrete(labels = c("1" = "Deep\nthermocline",
                              "2" = "Intermediate\nthermocline", 
                              "3" = "Shallow\nthermocine"))
cluster.layman
ggsave(cluster.layman, file = "Graphics/Figures/Figure2_LaymanCluster.pdf", width = 5.5, height = 4)



# If you want to define ROPE as 10% of the SD of the outcome
y_sd = bay.lay %>% 
  group_by(metric) %>%
  summarize(sd = sd(metric.value)) %>%
  ungroup() %>%
  mutate(r1 = -.1*sd, 
         r2 = .1*sd)

bay.lay.rope = bay.lay %>% 
  ungroup() %>%
  pivot_wider(names_from =cluster, values_from = metric.value) %>%
  mutate(c1v2 = `1` - `2`,
         c1v3 = `1` - `3`,
         c2v3 = `2` - `3`) %>%
  left_join(y_sd) 

Fig2_statssummary = bay.lay.rope %>% 
  pivot_longer(cols = starts_with("c"),
               names_to = "contrast",
               values_to = "diff") %>% 
  group_by(metric, contrast) %>%
    summarize(rope_prop = 100 * (rope(diff, range = c(first(r1), 
                                                      first(r2))))$ROPE_Percentage,
              .groups = "drop", 
              PD =100* max(c(pd1 = mean(diff > 0), mean(diff < 0))))
  
write.csv(Fig2_statssummary, file = "Graphics/Tables/Fig2_StatsSummary.csv", 
          row.names =  F)            




