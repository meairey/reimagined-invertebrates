library(tidyverse)

taxon_frame = read.csv(file = "Data/CSVs/taxon_frame.csv") %>%
  select(TAXON, FAMILY) %>%
  unique()

sample = read.csv("../1.clean_isotope/isotope_sample.csv")


chem = read.csv(file = "Data/CSVs/richness.csv")

## updating the updated clean data set

iso = read.csv("Data/CSVs/updated_clean/iso_measurement.csv") %>%
  
  select(ISO_YSAMP_N,ITEM_N,GROUP, TAXON, D13C, D15N, d2HVSMOW, D34S) %>%
  left_join(taxon_frame) %>%
  left_join(sample, by = "ISO_YSAMP_N") %>%
  #filter(WATER != "LML") %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", 
                             MONTH > 7 ~ "fall")) %>%
  unite(community.name, c(WATER, season), sep = ".") %>%
  left_join(chem) 

## Deuterium across taxa
iso %>%
  filter(GROUP %in% c("INSECT","LEAF" )) %>%
  select(FAMILY, MONTH, d2HVSMOW, GROUP) %>%
  #na.omit() %>%
  filter(is.na(d2HVSMOW) ==F) %>%
  ggplot(aes(x = FAMILY, y = d2HVSMOW, col = MONTH)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) 


## Deuterium across DOC or something
iso %>%
  filter(GROUP %in% c("INSECT","LEAF" )) %>%
 # filter(ORDER.or == "odonata") %>%
  select(FAMILY, MONTH, d2HVSMOW, GROUP, DOC.1) %>%
  #na.omit() %>%
  filter(is.na(d2HVSMOW) ==F) %>%
   mutate(MONTH = case_when(MONTH < 7  ~ "S", MONTH > 7 ~ "F")) %>%
  ggplot(aes(x = DOC.1, y = d2HVSMOW, col = (MONTH))) +
  geom_smooth(method = lm) + 
  geom_point() +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
  xlab("DOC") + ylab("d2H") + 
  scale_color_manual("Season", values = wes_palette("Moonrise2", 2))
  


## Sulfur across taxa
iso %>%
  filter(GROUP %in% c("INSECT")) %>%
  select(FAMILY, MONTH, D34S, GROUP) %>%
  mutate(MONTH = case_when(MONTH < 7  ~ "S", MONTH > 7 ~ "F")) %>%
  #na.omit() %>%
  filter(is.na(D34S) ==F) %>%
  ggplot(aes(x = FAMILY, y = D34S, col = MONTH)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_boxplot() 

## Do you think it makes sense to standardize it to the d34S of the periphyton
iso%>%
  select(community.name, GROUP,FAMILY, MONTH, D34S, GROUP, depth.5mgL) %>%
  mutate(MONTH = case_when(MONTH < 7  ~ "S", MONTH > 7 ~ "F")) %>%
  #na.omit() %>%
  filter(is.na(D34S) ==F) %>%
  group_by(community.name, MONTH, depth.5mgL) %>%
  mutate(mean_peri = mean(D34S[GROUP == "PERI"], na.rm = TRUE)) %>%
  filter(is.na(mean_peri)==F) %>% 
  filter(GROUP == "INSECT") %>%
  mutate(D34S.standard = D34S - mean_peri) %>%
  ggplot(aes(x = depth.5mgL, y = D34S.standard)) + 
  geom_point() + 
  geom_smooth(method = lm)



iso %>%
  filter(GROUP %in% c("ZOOP")) %>%
  select(FAMILY, MONTH, D34S, GROUP) %>%
  mutate(MONTH = case_when(MONTH < 7  ~ "S", MONTH > 7 ~ "F")) %>%
  #na.omit() %>%
  filter(is.na(D34S) ==F) %>%
  ggplot(aes(x = FAMILY, y = D34S, col = MONTH)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_boxplot() 




## Sulfur across DOC or something
iso %>%
  filter(GROUP %in% c("FISH" )) %>%
 # filter(ORDER.or == "odonata") %>%
  select(community.name, FAMILY, MONTH, D34S, GROUP, DOC.1, depth.5mgL) %>%
  #na.omit() %>%
  filter(is.na(D34S) ==F) %>%
  ggplot(aes(x = depth.5mgL, y = D34S, col = community.name)) + 
  geom_point() +
  theme_minimal(base_size = 14) + 
  geom_smooth(method = lm, col = "black") + 
  xlab("Depth (m) where DO < 5 mgL")

