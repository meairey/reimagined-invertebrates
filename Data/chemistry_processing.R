library(tidyverse)
`%nin%` = Negate(`%in%`)
## AFRP waterchem data for scale lakes

# Setup AFRP data -------------
afrp_chem = read.csv("Data/CSVs/Water_Chemistry/AFRP_water_chem.csv")
afrp_sample = read.csv("Data/CSVs/Water_Chemistry/WCS_SAMPLE_AFRP.csv") 

acd = left_join(afrp_chem, afrp_sample, by = c("YSAMP_N")) %>%
  select(WATER, YEAR, MONTH, everything())

# Create a data frame of means
dat = acd %>% filter(YEAR > 2013) %>%
  group_by(WATER, METRIC) %>%
  summarize(mean.v1 = mean(VALUE_1, na.rm =T)) %>%
  pivot_wider(names_from = METRIC, values_from = mean.v1)

# Don't  need to re-write this every time
#write.csv(file = "Data/CSVs/ALCwaterChem.csv", dat)

## AFRP secchi depth
alc.sechi.depth = acd %>% 
  filter(METRIC %in% c("SECCHI DEPTH"), YEAR > 2000) %>%
  select(YEAR, WATER, METRIC, MONTH, VALUE_1) %>%
  mutate(season = case_when(MONTH > 7 ~ "fall", MONTH <= 7 ~ "spring")) %>%
  na.omit() %>%
  group_by(WATER, season) %>%
  summarize(sechi.depth = mean(VALUE_1))

## AFRP Secchi color
# This frame was made and then turned into a csv and then categorically assigned FUI values and then is now reloaded in separately  
#alc.sechi.color = acd %>% 
 # filter(METRIC %in% c("SECCHI COLOR"), YEAR > 2000) %>%
#  mutate(season = case_when(MONTH > 7 ~ "fall", MONTH <= 7 ~ "spring")) %>%
#  select( WATER, METRIC, season, VALUE_3)  %>%
#  na.omit() %>%
#  unique() %>%
#  select(-METRIC) %>%
#  rename(sechi.color = VALUE_3)
#### Read in the secchi color frame
alc.sechi.color = read.csv("Data/CSVs/Water_Chemistry/AlcSechiColor.csv") %>%
  group_by(WATER, season) %>%
  summarize(FUI.num = mean(FUI.num))

fui.scale = read.csv("Data/CSVs/Water_Chemistry/AlcSechiColor.csv") %>%
  select(sechi.color, FUI.num) %>%
  unique() %>%
  mutate(sechi.color = tolower(sechi.color))

fui.scale

alc.sechi = left_join(alc.sechi.depth, alc.sechi.color)

## Temp DO information from SCALE collections----------------------------------

## Temp and DO data frame
tdo = read.csv("Data/CSVs/Water_Chemistry/Temp_DO.csv") %>%
  separate(Date, into = c("month", "date", "year")) %>%
  mutate(season = case_when(month < 8 ~ "spring", month>7 ~ "fall"))


sechi = tdo %>% select(WATER,season,  Sechi.Depth..M., Visual.Water.Color..0.5M.) %>%
  na.omit() %>%
  rename(sechi.depth = Sechi.Depth..M., 
         sechi.color = Visual.Water.Color..0.5M.) %>%
  mutate(sechi.color = tolower(sechi.color)) %>%
  left_join(fui.scale) %>%
  select(-sechi.color) %>%
  mutate(FUI.num = case_when(is.na(FUI.num)==T  ~ 15, is.na(FUI.num)==F ~ FUI.num)) ## Fills in the missing sechi color of light yellow and assigns it the same number as light brown

## combines the sechi data for both alc and ALSC lakes
scale.sechi = rbind(sechi, alc.sechi) %>%
  unite(ID, c(WATER, season), sep = ".", remove = F) %>%
  filter(ID %nin% c("COM.spring", "LML.spring", "GNL.spring", "ETL.spring")) %>%
  rbind(c("HTL.fall", "HTL", "fall", 4.5, 15)) ## Missing HTL secchi depth (4.5) and FUI number (15)

## Data with all morphological data for each lake
morphology = read.csv("Data/CSVs/Water_Chemistry/chem_profiles_clean.csv") %>% 
  select(WATER, Pond_num, Lake, Elevation, max_depth, Volume)

## Thermocline -----------------
thermocline = read.csv(file = "Data/CSVs/Water_Chemistry/Thermocline.csv")
## Oxycline 
oxycline.metrics = read.csv(file = "Data/CSVs/Water_Chemistry/Oxycline.csv")

## Compiled chem data -----------


## You need both seasons in the chemistry. Fall only for NMDS, but use both seasons for regressions
chemistry = oxycline.metrics %>% 
  left_join(scale.sechi) %>%
  left_join(thermocline %>% 
              select(ID, mid_depth, bottom)) %>%
  left_join(read.csv(file = "Data/CSVs/chemistry_seasonal.csv") %>%
              select(ID, DOC_update)) %>%
  left_join(morphology)

write.csv(chemistry,file = "Data/CSVs/chemistry_seasonal.csv", row.names = F)


chem$WATER %nin% chemistry$WATER



#### Includes information from above that has already been combined into csv

chem = read.csv("Data/CSVs/Water_Chemistry/chem_profiles_clean.csv") %>% 
  unite("community.name", c(WATER, season), sep = ".", remove = F) %>%
  select(-DDO5, -TDO5,) %>% ## remove these from the csv if you want to rejoin in the updated metrics (below)
  arrange(WATER) %>%
  left_join(thermocline %>%
              select(-grad_max, -thresh, -top)) %>% ## Join in data from Thermocline.R
  left_join(oxycline.metrics, by = c("community.name" = "ID")) ## not using top for the analysis 

write.csv(chem, file = "Data/CSVs/Water_Chemistry/chem_profiles_clean.csv", row.names = F) ## This rewrites the data sheet with the updated metrics



### -------------------------------------------------------------------------
### Extra code ---- checking data frames
## Figuring out the DOC years for the nmds

acd %>% 
  filter(METRIC == 'DISSOLVED ORGANIC CARBON') %>%
  group_by(WATER) %>%
  summarize(min_year = min(YEAR))

## Checking that the metrics are right ---------------

acd %>% select(METRIC, UNIT_1) %>%
  unique()
  
## View lakes water chem through time
 acd %>% filter(WATER %in% c("COM", "ETL","LML", "GNL"),
                METRIC %in% c("pH"))  %>%
  ggplot(aes(x = YEAR, y = VALUE_1, col = WATER)) + 
  geom_point(alpha = .5) + 
  geom_smooth(method = lm, col= "black") +
  facet_wrap(~ WATER)


## All metrics
acd %>% filter(WATER %in% c("COM", "ETL","LML", "GNL")) %>%
  group_by(WATER, METRIC) %>%
  mutate(total_observations = n()) %>%
  ungroup() %>% 
  group_by(METRIC) %>%
  mutate(total_lakes = length(unique(WATER)))%>%
  select(total_observations, total_lakes,  everything()) %>%
  filter(total_observations > 3 & total_lakes > 3, 
         METRIC %nin% c("SECCHI COLOR","UV254")) %>%
                 ggplot(aes(x = YEAR, y = VALUE_1, col = WATER)) + 
                 geom_point(alpha = .5) + 
                 geom_smooth(method = lm) +
                 facet_wrap(~ METRIC, scales = "free")


