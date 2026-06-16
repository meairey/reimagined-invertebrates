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

 


ggsave(file = "Graphics/Figures/Supp_Figure_1.pdf", plot = FigureS1, width = 7, height = 5, units = "in")
