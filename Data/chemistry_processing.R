
## AFRP waterchem data for scale lakes

# Load in the AFRP data
afrp_chem = read.csv("Data/CSVs/AFRP_water_chem.csv")
afrp_sample = read.csv("Data/CSVs/WCS_SAMPLE_AFRP.csv") 





#afrp chemistry data
acd = left_join(afrp_chem, afrp_sample, by = c("YSAMP_N")) %>%
  select(WATER, YEAR, MONTH, everything())



# Create a data frame of means
dat = acd %>% filter(YEAR > 2013) %>%
  group_by(WATER, METRIC) %>%
  summarize(mean.v1 = mean(VALUE_1, na.rm =T)) %>%
  pivot_wider(names_from = METRIC, values_from = mean.v1)

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
 
 
# Don't  need to re-write this every time
#write.csv(file = "Data/CSVs/ALCwaterChem.csv", dat)


## AFRP sechi

alc.sechi.depth = acd %>% 
  filter(METRIC %in% c("SECCHI DEPTH"), YEAR > 2000) %>%
  select(YEAR, WATER, METRIC, MONTH, VALUE_1) %>%
  mutate(season = case_when(MONTH > 7 ~ "fall", MONTH <= 7 ~ "spring")) %>%
  na.omit() %>%
  group_by(WATER, season) %>%
  summarize(sechi.depth = mean(VALUE_1))



# This frame was made and then turned into a csv and then categorically assigned FUI values and then is now reloaded in separately  
alc.sechi.color = acd %>% 
  filter(METRIC %in% c("SECCHI COLOR"), YEAR > 2000) %>%
  mutate(season = case_when(MONTH > 7 ~ "fall", MONTH <= 7 ~ "spring")) %>%
  select( WATER, METRIC, season, VALUE_3)  %>%
  na.omit() %>%
  unique() %>%
  select(-METRIC) %>%
  rename(sechi.color = VALUE_3)


alc.sechi.color = read.csv("Data/CSVs/AlcSechiColor.csv") %>%
  group_by(WATER, season) %>%
  summarize(FUI.num = mean(FUI.num))

fui.scale = read.csv("Data/CSVs/AlcSechiColor.csv") %>%
  select(sechi.color, FUI.num) %>%
  unique() %>%
  mutate(sechi.color = tolower(sechi.color))

fui.scale

alc.sechi = left_join(alc.sechi.depth, alc.sechi.color)

## Temp DO information from SCALE collections----------------------------------

## Temp and DO data frame
tdo = read.csv("Data/CSVs/Temp_DO.csv") %>%
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
scale.sechi = rbind(sechi, alc.sechi)



chem = read.csv("Data/CSVs/chem_profiles_clean.csv")


# Temp/Depth profiles
tdo %>%
  ggplot(aes(y = depth, x = temp, col = Name  )) + 
  geom_line() +
  facet_wrap(~season) +
  scale_y_reverse()

# Calculating location of max rate of change in temp to get the thermoclines
max.rate = tdo %>%
  group_by(WATER, season) %>%
  mutate(rate = temp - lag(temp)) %>%
  select(WATER, season, depth, temp, rate) %>%
  filter(rate ==  min(rate, na.rm = T)) %>%
  rename(thermo_depth = depth,
         thermo_temp = temp)

# Calculating the minimum temperature where DO is above 5mgl
min_temp =  tdo %>%
  select(WATER, season, depth, temp, DO) %>%
  filter(DO > 4.75) %>%
  group_by(WATER,season ) %>%
  filter(depth == max(depth)) %>%
  mutate(temp_do = "min") %>%
  rename(min_temp = temp, 
         min_depth = depth, 
         min_do = DO)

## Vizualize where minimum temp before DO < 5 mgl
tdo %>%
  ggplot(aes(y = depth, x = temp, col = as.factor(WATER)  )) + 
  geom_line() +
  geom_point(data =(min_temp%>% na.omit()), 
             aes(x = min_temp,y = min_depth), color = "black") +
  facet_wrap(~season) +
  scale_y_reverse()


## Temp and DO at meaningful depths for inverts? Maybe 1 or 2 m?

tdo




# Creating a chemistry data frame that picks out the 
chemistry = chem %>% left_join(max.rate) %>%
  left_join(min_temp, by = c("WATER", "season")) %>%
  left_join(scale.sechi) 

chemistry[6, "sechi.depth"] = 4.500000 ## Missing HTL fall sechi data
chemistry[6, "FUI.num"] =15.00000 ## Missing HTL fall sechi data


chemistry

## Save this to be imported into the NMDS script
save(file = "Data/RData/chem_data.RData", chemistry)















