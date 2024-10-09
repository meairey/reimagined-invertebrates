`%nin%` = Negate(`%in%`)
### Isotope data


## Load in taxon data


taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Load in isotope measurement file
data = read.csv("Data/CSVs/iso_data_clean.csv")

## Read in isotope sample file
sample = read.csv("Data/CSVs/isotope_sample.csv")

## Join together -----------------

df = data %>% left_join(sample, by = "ISO_YSAMP_N") %>%
  left_join(taxon_frame) %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  mutate(GROUP = str_replace(GROUP, "CLAM","BIVALVE"),
         GROUP = str_replace(GROUP,"MUSSEL", "BIVALVE")) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall")) 

common = df %>% 
  filter(is.na(ORDER) != T) %>%
  group_by(WATER, ORDER, FAMILY) %>%
  mutate(total_n = n()) %>%
  filter(total_n > 4) %>% ## Filter out taxa where # in order > 4
  filter(ORDER %nin% c("trichoptera","coleoptera","hemiptera")) %>% ## remove orders that contribute most to dissimilarity between communities 
  summarize(total = n()) %>%
  print(n = 100)


season.df = df %>% 
  filter(FAMILY %in% common$FAMILY | GROUP %in% c("ZOOP", "PERI","CLAM","MUSSEL", "SNAIL", "LEAF")) %>%
  mutate(graph_id = case_when(is.na(ORDER) == T ~ GROUP, 
                              is.na(ORDER) == F ~ ORDER)) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall")) %>%
  mutate(graph_id = case_when(graph_id == "CLAM" ~ "bivalvia", graph_id != "CLAM" ~ graph_id))


df %>% 
  filter(MONTH > 7, 
         GROUP %nin% c("CHECK", "ARACHNIDA", "CRAY"))%>%
  ggplot(aes(x = D13C, y = D15N, col = GROUP)) + geom_point() + 
  facet_wrap(~WATER) + 
  stat_ellipse(level = .4) +
  theme_minimal()

## All insects across waters in spring
 df %>% 
  filter(MONTH < 8, GROUP == "INSECT", is.na(ORDER) != T) %>%
   select(D13C, D15N, ORDER, WATER) %>%
   ggplot(aes(x = D13C, y = D15N, col =(ORDER))) + geom_point() + 
   stat_ellipse(level = .4) + 
   facet_wrap(~WATER)
 
 ## Mayflies, periphyton and zooplankton
 df %>% 
   filter(MONTH < 7, GROUP == "INSECT" & ORDER == "ephemeroptera"| GROUP == "PERI" | GROUP == "ZOOP" | GROUP == "LEAF") %>%
   ggplot(aes(x = D13C, y = D15N, col = GROUP)) + geom_point() + 
   facet_wrap(~WATER, scales = "free") + 
   stat_ellipse(level = .4)
 
### Protocol --------------------------------
 
 
 
## Use inference on rare taxa to filter out ----------------


 

 
  
## Insects in the spring
 df %>% 
   filter(MONTH < 8, FAMILY %in% common$FAMILY) %>%
   select(D13C, D15N, ORDER, WATER) %>%
   ggplot(aes(x = D13C, y = D15N, col =(ORDER))) + geom_point() + 
   stat_ellipse(level = .4) + 
   facet_wrap(~WATER)
 
### Main groups by season, filtered for common insects
 df %>% 
   filter(MONTH > 7, # change direction of operator to change season
          FAMILY %in% common$FAMILY | GROUP == "ZOOP" | GROUP == "PERI") %>% 
   mutate(graph_id = case_when(is.na(ORDER) == T ~ GROUP, is.na(ORDER) == F ~ ORDER)) %>%
   select(D13C, D15N, graph_id, WATER) %>%
   ggplot(aes(x = D13C, y = D15N, col =(graph_id))) + geom_point() + 
   stat_ellipse(level = .4) + 
   theme_minimal() +
   labs(col = "Order") + 
   facet_wrap(~WATER, scales= "free")
 
 
 ### Seasonality impact on signatures
 



season.df %>% 
  ggplot(aes(x = graph_id, y = D13C, col = season)) + 
  geom_boxplot(fill = NA) +
  facet_wrap(~WATER) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

## CV Coef of var -----------------------------------
CV = season.df %>% ## coefficient of variance by sampling occasion
  ungroup() %>%
  group_by(graph_id, WATER, season, SITE_N) %>%
  mutate(mean.C = mean(D13C, na.rm =T),
         mean.N = mean(D15N, na.rm =T),
         sd.C = sd(D13C, na.rm =T),
         sd.N = sd(D15N, na.rm =T)) %>%
  mutate(CV.C = abs(sd.C / mean.C), 
         CV.N = abs(sd.N / mean.N))

## D13C CV
CV %>%
  ggplot(aes(x = graph_id, y = CV.C, col = season)) + 
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  xlab("Group") +
  ylab("D13C CV") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

## D15N CV
CV %>%
  ggplot(aes(x = graph_id, y = CV.N, col = season)) + 
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  xlab("Group") +
  ylab("D15N CV") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

## Comparison of taxa to baselines

season.df %>% 

  ggplot(aes(x = graph_id, y = D13C)) + geom_point() +
  facet_wrap(~WATER)  +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))


season.df %>% 
  filter(GROUP == "PERI", WATER %in% c("MSL","DTL") )

## Baselines of algae across sites -----------
season.df %>% 
  filter(GROUP == "PERI", WATER %in% c("MSL","DTL"))%>%
  group_by(WATER) %>%
  mutate(SITE_index = as.numeric(as.factor(SITE_N))) %>%
  ggplot(aes(x = D13C, y = D15N, col = as.factor(SITE_index)))  + 
  geom_point(size =2) + 
  theme_minimal() +
  facet_wrap(~season+WATER) +
  labs(col = "Site - Benthic Algae")

## Baselines of zooplankton across sites


season.df %>% 
  filter(GROUP == "ZOOP", WATER %in% c("MSL","DTL"))%>%
  group_by(WATER) %>%
  mutate(SITE_index = as.numeric(as.factor(SITE_N))) %>%

  ggplot(aes(x = D13C, y = D15N, col = as.factor(SITE_index)))  + 
  geom_jitter(size =2) + 
  theme_minimal() +
  facet_wrap(~season+WATER) +
  labs(col = "Site - Benthic Algae")


  
## Residuals
peri =  CV %>%
  filter(GROUP == "PERI") %>% 
  ungroup() %>%
  select(WATER, season, mean.C, mean.N, SITE_N) %>%
    unique() %>% na.omit() 

consumer = CV %>% filter(GROUP %in% c("INSECT", "SNAIL", "CLAM", "MUSSEL") )%>%
  select(WATER,GROUP,  season, D13C, D15N, FAMILY)
  

residuals = left_join(consumer, peri, by = c("WATER", "season", "SITE_N")) %>% 
  mutate(residual.C = as.numeric(D13C) - mean.C, 
         residual.N = as.numeric(D15N) - mean.N) %>%
  ungroup() %>%
  mutate(graph_id = case_when(is.na(FAMILY) == T ~ GROUP, 
                              is.na(graph_id) == F ~ FAMILY)) %>%
  mutate(graph_id = tolower(graph_id))

## Figures of residuals
residuals %>% 
  ggplot(aes(x = graph_id , y = abs(residual.C), col = season)) +
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
  ylab("Residuals d13C") +
  xlab("Group")
  
residuals %>% 
  ggplot(aes(x = graph_id, y = abs(residual.N), col = season)) +
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
  ylab("Residuals d15N") + 
  xlab("Group")

## ANOVA analysis of the d15N residuals
anova.model.data = residuals %>%
  unite("ID", c(graph_id, season)) %>% 
  select(ID, residual.C, residual.N) %>%
  na.omit()

model.N = aov(data = anova.model.data, residual.N ~ID)
summary(model.N)
tukey.test = TukeyHSD(model)
tukey.test$ID %>% as.data.frame() %>%
  filter(`p adj` < .05) %>% 
  rownames_to_column(var = "ID") %>%
  separate(ID, into = c("group1", "season1",
                        "group2", "season2")) %>%
  filter(group1 == group2 | season1 == season2) %>%
  arrange(season1, group1)

## Heptageniidae is significantly lower in d15N residuals than aeshnidae, coenagrionidae, and gomphidae in the fall

## In the spring, heptageniidae is significantly lower than aesh, coeng, cordu, and gomphid. It does not differ from viviparidae or snails.

## It does look like gomphids are a bit more enriched in the fall

## ANOVA analysis of the d13C residuals

# No significant anove between groups for the d13C residuals?
model.C = aov(data = anova.model.data, abs(residual.C) ~ID)
model.C %>% summary()
## Not significant don't run tukey?

## Residuals of d13C away from leaves ----------------------------
leaves.re =  CV %>%
  filter(GROUP == "LEAF") %>% 
  ungroup() %>%
  select(WATER, season, mean.C, mean.N) %>%
  unique() %>% na.omit() 

consumer = CV %>% filter(GROUP %in% c("INSECT", "SNAIL", "CLAM", "MUSSEL") )%>%
  select(WATER,GROUP,  season, D13C, D15N, FAMILY)


residuals = left_join(consumer, leaves.re, by = c("WATER", "season")) %>% 
  mutate(residual.C = as.numeric(D13C) - mean.C, 
         residual.N = as.numeric(D15N) - mean.N) %>%
  ungroup() %>%
  mutate(graph_id = case_when(is.na(FAMILY) == T ~ GROUP, 
                              is.na(graph_id) == F ~ FAMILY)) %>%
  mutate(graph_id = tolower(graph_id))

## Figures of residuals
residuals %>% 
  ggplot(aes(x = graph_id , y = abs(residual.C), col = season)) +
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
  ylab("Residuals d13C") +
  xlab("Group")

residuals %>% 
  ggplot(aes(x = graph_id, y = abs(residual.N), col = season)) +
  geom_boxplot(outliers = F) +
  geom_point() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
  ylab("Residuals d15N") + 
  xlab("Group")

## Lakes with fish sampled

fish_lakes = df %>% group_by(WATER,  GROUP) %>%
  mutate(total_group = n()) %>%
  filter(GROUP == "FISH" ) %>% 
  ungroup() %>%
  select(WATER, MONTH) %>%
  unique() %>%
  na.omit() %>%
  mutate(select = "YES") %>%
  unite("ID", c(WATER, MONTH), remove = F)

fish_lakes 

df_fish = df %>% left_join(fish_lakes) %>%
  filter(is.na(select) == F) 

df_fish %>% 
  mutate(GROUP = tolower(GROUP)) %>%
  filter(GROUP %nin% c("check", "cray", "bivalve", "snail")) %>%
  ggplot(aes(x = D13C, y = D15N, col = GROUP)) + 
  geom_point() +
  facet_wrap(~ID, scales = "free") +
  stat_ellipse(level = .4)
## Time for SIBER

## Differences in Family
df %>% filter(FAMILY == "aeshnidae") %>%
  ggplot(aes(x = D13C, y  = D15N, col = GENUS )) +
  geom_point() +
  facet_wrap(~WATER + season)


taxon_frame %>% filter(FAMILY == "sphaeriidae")
taxon_frame %>% group_by(FAMILY) %>%
  summarize(total = n()) %>%
  filter(total > 3) %>%
  print(n = 100)

df %>% group_by(FAMILY,GENUS) %>%
  summarize(total = n()) %>%
  na.omit() %>%
  filter(total > 5) %>%
  print(n =   100)

na.genus = df %>% filter(is.na(GENUS) == T)
## I don't really think that there is enough to split up by genus. So, for start I'll do SIBER with just family level
            