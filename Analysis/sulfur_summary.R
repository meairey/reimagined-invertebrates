## Upload Data
`%nin%` = Negate(`%in%`)
measurement = read.csv("Data/CSVs/updated_clean/iso_measurement.csv")
sample = read.csv("Data/CSVs/updated_clean/iso_sample.csv")
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>%
  select(ORDER, FAMILY, TAXON) %>%
  unique() %>%
  arrange(TAXON)
richness = read.csv("Data/CSVs/richness.csv") %>%
  separate(community.name, into = c("WATER", "SEASON"))
load("Data/RData/richness_cluster.RData")  

data = left_join(measurement, sample) %>% 
  left_join(taxon_frame) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall"))


data %>% 
  filter(is.na(D34S) ==F) %>% 
  filter(is.na(ORDER) == T) %>% 
  filter(GROUP %nin% c("LEAF", "PERI")) %>%
  select(WATER, TAXON)




## Checking if fish with lower d34S are in higher DOC lakes?
data %>% 
  filter(is.na(D34S) ==F) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  filter(GROUP == "FISH") %>%
  ggplot(aes(x = sechi.depth, y = D34S)) + 
  geom_point() + 
  geom_smooth(method = lm)


## Odonates?
data %>% 
  filter(is.na(D34S) ==F) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  #filter(ORDER == "odonata") %>%
  filter(GROUP == "FISH") %>%
  ggplot(aes(x = max_depth, y = D34S)) + 
  geom_point(aes(col = WATER)) + 
  geom_smooth(method = lm) 


data %>% 
  filter(is.na(D34S) ==F) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  #filter(ORDER == "odonata") %>%
  filter(GROUP == "FISH") %>%
  ggplot(aes(x = thermo_depth , y = D34S)) + 
  geom_point(aes(col = WATER)) + 
  geom_smooth(method = lm) +
  facet_wrap(~season)

data %>% 
  filter(is.na(D34S) ==F) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  
  ggplot(aes(x = as.factor(season) , y = D34S)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~ORDER) + 
  geom_boxplot()


data %>% 
  filter(is.na(D34S) ==F) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  filter(GROUP == "INSECT") %>%
  ggplot(aes(x = as.factor(season) , y = D34S)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~WATER) + 
  geom_boxplot()




data %>% 
  filter(is.na(D34S) ==F, 
         GROUP %nin% c("PERI", "LEAF")) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  ggplot(aes(x = WATER, y = D34S, col = ORDER)) + 
  geom_point() +
  geom_boxplot()


data %>% 
  filter(is.na(D34S) ==F, 
         GROUP %in% c("PERI", "LEAF")) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  ggplot(aes(x = WATER, y = D34S, col = GROUP)) + 
  geom_point() +
  geom_boxplot()


## Thinking maybe i should standardize this by periphyton?


periphyton_means = data %>% 
  filter(is.na(D34S) ==F, 
         GROUP %in% c("PERI", "LEAF")) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  group_by(WATER) %>%
  summarize(base_S = mean(D34S))


data %>% 
  filter(WATER %nin% c("UCL", "LML")) %>%
  filter(is.na(D34S) ==F, 
         GROUP %nin% c("PERI", "LEAF")) %>%
  left_join(richness, by = c("WATER", "season")) %>%
  left_join(periphyton_means) %>%
  mutate(S_cor = D34S - base_S) %>%
  ggplot(aes(x = sechi.depth, y = S_cor, col = ORDER)) + 
  geom_point() +
  geom_smooth(method= lm, se = F)

