## Cleaning isotope data

## Inverts
iso_dat = read.csv("Data/CSVs/done/Database_isotope_0924_CLEANING.csv")


meas = read.csv("Data/CSVs/SCALE_MTS_2023.csv") %>%
  select(YSAMP,ORDER..OR.ABOVE., FAMILY, GENUS, MTS_N,TOTAL_N)


dat = left_join(iso_dat, meas, by = c("ITEM_N" = "MTS_N"))  %>%
  mutate(GENUS = na_if(GENUS, "")) %>%
  mutate(FAMILY = na_if(FAMILY, "")) %>%
  mutate(TAXON = case_when(is.na(GENUS) == F ~ GENUS,
                           is.na(GENUS) == T~ FAMILY,
                           is.na(FAMILY) == T ~ ORDER..OR.ABOVE.),
         ISO_YSAMP_N = YSAMP) %>%
  mutate(TAXON = tolower(TAXON)) %>%
  select(-FAMILY,-GENUS, -YSAMP)
  
dat$TAXON

#invert_taxa = dat$TAXON
#save(invert_taxa, file = "Data/RData/invert_taxa.RData")
#write.csv(file = "Data/CSVs/data cleaning/invert_iso_cleaning.csv", dat)


iso_dat[128,]


## Additional samples 
iso_dat.additional = read.csv("Data/CSVs/data cleaning/additional_samples.csv")
dat = left_join(iso_dat.additional, meas, by = c("ITEM_N" = "MTS_N"))  %>%
  mutate(GENUS = na_if(GENUS, "")) %>%
  mutate(FAMILY = na_if(FAMILY, "")) %>%
  mutate(TAXON = case_when(is.na(GENUS) == F ~ GENUS,
                           is.na(GENUS) == T~ FAMILY,
                           is.na(FAMILY) == T ~ ORDER..OR.ABOVE.),
         ISO_YSAMP_N = YSAMP) %>%
  mutate(TAXON = tolower(TAXON)) %>%
  select(-FAMILY,-GENUS, -YSAMP)



## Fish

iso_dat = read.csv(file = "Data/CSVs/data cleaning/invert_iso_cleaning.csv")

fish_dat = read.csv(file = "Data/CSVs/data cleaning/done/fish_cleaning.csv") %>%
  select(-MONTH, -YEAR, -WATER)

fish_sample = read.csv(file = "Data/CSVs/data cleaning/done/fish_cleaning.csv") %>%
  select(YSAMP, MONTH, YEAR, WATER) %>%
  unique()


fish_species = fish_dat$SPECIES %>% unique()
save(fish_species, file = "Data/RData/fish_species.RData")


#write.csv(file = "Data/CSVs/data cleaning/fish_sample.csv", fish_sample)

dat = left_join(iso_dat, fish_dat, by = c("ITEM_N" = "FISH_N")) %>%
  mutate(TAXON = case_when(is.na(SPECIES) == F ~ SPECIES, 
                           is.na(SPECIES) == T ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ YSAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select(-YSAMP, -SPECIES)


## additional samples 
dat.F  = dat %>% left_join(fish_dat, by = c("ITEM_N" = "FISH_N"))%>%
  mutate(TAXON = case_when(is.na(SPECIES) == F ~ SPECIES, 
                           is.na(SPECIES) == T ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ YSAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select(-YSAMP, -SPECIES)

#write.csv(fish_sample, "Data/CSVs/data cleaning/fish_sample.csv")
#write.csv(dat, "Data/CSVs/data cleaning/fish_cleaned.csv")

## Algae

iso_dat = read.csv(file = "Data/CSVs/data cleaning/fish_cleaned.csv")

algae_sample = read.csv("Data/CSVs/data cleaning/ALGAE_CLEANING.csv") %>%
  select(YSAMP, WATER, MONTH, SITE) %>%
  unique()

#write.csv(algae_sample, "Data/CSVs/data cleaning/algae_sample.csv")

algae_dat = read.csv("Data/CSVs/data cleaning/ALGAE_CLEANING.csv") %>%
  select(YSAMP, ID, GROUP) %>%
  rename(algae = GROUP)

dat = left_join(iso_dat, algae_dat, by = c("ITEM_N" = "ID"))  %>%
  mutate(TAXON = case_when(is.na(TAXON) == T ~ algae,
                   is.na(TAXON) == F ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ YSAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select( -algae, -YSAMP)

#write.csv(dat, "Data/CSVs/data cleaning/algae_cleaned.csv")


### Read in zooplankton
iso_dat = read.csv("Data/CSVs/data cleaning/algae_cleaned.csv")

zoop_data = read.csv("Data/CSVs/data cleaning/zooplankton_cleaning.csv") %>% 
  select(ZOOP_N, YSAMP) %>%
  mutate(zoop = "zooplankton")

zoop_sample = read.csv("Data/CSVs/data cleaning/zooplankton_cleaning.csv") %>%
  select( WATE, MONTH, YEAR, SITE, YSAMP) %>%
  unique()

#write.csv(zoop_sample, file = "Data/CSVs/data cleaning/zooplankton_sample.csv")

dat = iso_dat %>%
  left_join(zoop_data, by = c("ITEM_N" = "ZOOP_N")) %>%
  mutate(TAXON = case_when(is.na(TAXON) == T ~ zoop,
                           is.na(TAXON) == F ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ YSAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select(-zoop, -YSAMP)

#write.csv(dat, "Data/CSVs/data cleaning/zooplankton_cleaned.csv")

## Read in leaves

iso_dat = read.csv("Data/CSVs/data cleaning/done/zooplankton_cleaned.csv")

leaf_dat = read.csv("Data/CSVs/data cleaning/done/leaf_cleaning.csv") %>%
  select(SIC_SAMP, ID, species) %>%
  mutate(leaf = "leaf")

leaf_sample = read.csv("Data/CSVs/data cleaning/done/leaf_cleaning.csv") %>%
  select(YSAMP_N, SIC_SAMP, MONTH, YEAR, LAKE) %>%
  unique()

#write.csv(leaf_sample, file = "Data/CSVs/data cleaning/leaf_sample.csv")

dat = iso_dat %>% left_join(leaf_dat, by = c("ITEM_N" = "ID")) %>%
  mutate(TAXON = case_when(is.na(TAXON) == T ~ species,
                           is.na(TAXON) == F ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ SIC_SAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select(-leaf, -SIC_SAMP)

## additional samples 

dat.leaf = dat.F %>% left_join(leaf_dat, by = c("ITEM_N" = "ID")) %>%
  mutate(TAXON = case_when(is.na(TAXON) == T ~ species,
                           is.na(TAXON) == F ~ TAXON),
         ISO_YSAMP_N = case_when(is.na(ISO_YSAMP_N) == T ~ SIC_SAMP, 
                                 is.na(ISO_YSAMP_N) == F ~ ISO_YSAMP_N)) %>%
  select(-leaf, -SIC_SAMP)

write.csv(dat.leaf, "Data/CSVs/data cleaning/additional_cleaned.csv")

## Cleaned csv
write.csv(dat, "Data/CSVs/data cleaning/cleaned.csv")

## Categories
dat = read.csv(file = "Data/CSVs/data cleaning/cleaned.csv")

dat = dat %>%
  mutate(CATEGORY = case_when(TAXON %in% fish_species ~ "VERT", 
                           TAXON == "zooplankton" ~ "INVERT",
                           TAXON == "PERIPHYTON"~ "ALGA", 
                           TAXON %in% na.omit(invert_taxa) ~ "INVERT")) %>%
  
  mutate(GROUP = case_when(TAXON %in% fish_species ~ "FISH", 
                              TAXON == "zooplankton" ~ "ZOOP",
                              TAXON == "PERIPHYTON"~ "PERI", 
                              TAXON %in% na.omit(invert_taxa) ~ "INSECT")) 
write.csv(dat, "Data/CSVs/data cleaning/cleaned.csv")


## Taxa join sheet -------------------

dat = read.csv("Data/CSVs/data cleaning/cleaned_hand_updated.csv")
## Read in 

load("Data/RData/FullData.RData")
full_number = full %>% select(MTS_N, TOTAL_N)

try = left_join(dat, full_number, by = c("ITEM_N" = "MTS_N")) %>%
  mutate(NUM = TOTAL_N) %>%
  select(-TOTAL_N)
#write.csv(try, "Data/CSVs/data cleaning/cleaned_numbers.csv")


try = read.csv("Data/CSVs/data cleaning/cleaned_numbers.csv")
taxon_frame = full %>% 
  select(ORDER..OR.ABOVE., FAMILY, GENUS) %>%
  unique() %>%
  mutate(GENUS = tolower(GENUS),
         FAMILY = tolower(FAMILY), 
         ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>%
  mutate(GENUS = na_if(GENUS, "")) %>%
  mutate(FAMILY = na_if(FAMILY, "")) %>%
  mutate(TAXON = case_when(is.na(GENUS) == F ~ GENUS,
                           is.na(GENUS) == T~ FAMILY,
                           is.na(FAMILY) == T ~ ORDER..OR.ABOVE.)) %>%
  
  mutate(FAMILY = case_when(TAXON == "tetragoneuria" ~"corduliidae",
                            TAXON != "tetragoneuria" ~ FAMILY )) %>%
  unique()

## Write the taxon frame for us in other scripts
save(taxon_frame, file = "Data/RData/taxon_frame.RData")

write.csv(taxon_frame, file = "Data/CSVs/taxon_frame_new.csv")

taxon_frame = read.csv("Data/CSVs/taxon_frame.csv")

## Final isotope data frame is saved in CSVs iso_dat_clean.csv
## It should also be kept in the file for the database upload. Please remember to update that...




## Checking that all the samples I submitted got run

submitted = read.csv("Data/CSVs/data cleaning/COIL_CN_samples_submitted.csv")

## Load in isotope measurement file
data = read.csv("Data/CSVs/iso_data_clean.csv")

not_run = submitted[which(submitted$Sample_ID %in% data$ITEM_N == FALSE),]

write.csv(not_run, "Data/CSVs/data cleaning/not_run_loaded.csv")



### sulfur data

library(tidyverse)

## Load in isotope measurement file
data = read.csv("Data/CSVs/iso_data_clean.csv")
sample = read.csv("Data/CSVs/isotope_sample.csv")
## load in sulfur file

sulfur = read.csv("Data/CSVs/data cleaning/sulfur.csv")

## Taxon frame
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv")

frank = left_join(data, sulfur) %>%
  group_by(ITEM_N) %>%
  unique()




write.csv(frank, file = "Data/CSVs/data cleaning/sulfur_combined.csv")


s = read.csv("Data/CSVs/data cleaning/sulfur_combined.csv")

s = dim(unique(s))

write.csv(s, file ="Data/CSVs/data cleaning/sulfur_combined.csv" )

s.dat = frank %>% 
  select(ISO_YSAMP_N, ITEM_N,CATEGORY, GROUP, X34S, TAXON) %>%
  unique() %>%
  filter(is.na(X34S) == F) %>%
  left_join(sample) %>%
  left_join(taxon_frame)

s.dat %>% 
  ggplot(aes(x = GROUP, y = X34S)) + 
  geom_boxplot()


s.dat %>%
  ggplot(aes(x = WATER, y = X34S)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~GROUP)

s.dat %>% filter(WATER == "UCL", GROUP == "FISH") %>%
  ggplot(aes(x = WATER, y = X34S, col = TAXON)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~GROUP)
  

s.dat %>% filter( GROUP == "INSECT") %>%
  ggplot(aes(x = WATER, y = X34S, col = ORDER)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~GROUP)



data = read.csv("Data/CSVs/updated_clean/iso_measurement.csv")
sample = read.csv("Data/CSVs/updated_clean/iso_sample.csv")

dim(data)
left_join(data, sample) %>% dim()


library(tidyverse)
## Adding in the hydrogen data

data = read.csv("Data/CSVs/data cleaning//iso_measurement.csv")
sample = read.csv("Data/CSVs/data cleaning//iso_sample.csv")
H_data = read.csv("Data/CSVs/data cleaning/H_data.csv")

new = left_join(data, H_data, by = c("ITEM_N" = "Sample.ID")) %>%
  left_join(sample) %>%
  left_join(taxon_frame)

write.csv(left_join(data, H_data, by = c("ITEM_N" = "Sample.ID")), file = "Data/CSVs/data cleaning/updated_H.csv")

cat = data %>% 
  group_by(ITEM_N) %>%
  summarize(count.a =  n())


dog = new %>% 
  group_by(ITEM_N) %>%
  summarize(count.b =  n())

left_join(cat, dog) %>%
  filter(count.a != count.b)

new %>% 
  ggplot(aes(x = WATER, y = d2HVSMOW, col = GROUP)) + 
  geom_jitter() +
  geom_boxplot()

new %>%
  filter(GROUP %in% c("INSECT","LEAF" )) %>%
  select(FAMILY, d2HVSMOW, GROUP) %>%
  #na.omit() %>%
  ggplot(aes(x = FAMILY, y = d2HVSMOW, col = MONTH)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) 


new %>% filter(GROUP %in% c("LEAF", "PERI", "ZOOP")) %>%
  ggplot(aes(x = GROUP, y = d2HVSMOW, col = as.factor(MONTH))) + 
  geom_point() +
  facet_wrap(~WATER)

richness = read.csv("Data/CSVs/richness.csv") %>%
  separate(community.name, into = c("WATER", "season"))

H.rich = new %>%
  #filter(GROUP %in% c("INSECT","LEAF" )) %>%
  select(ISO_YSAMP_N, MONTH, WATER,   FAMILY, d2HVSMOW, GROUP) %>%
  mutate(season = case_when(MONTH > 7 ~ "fall",
                          MONTH < 7 ~ "spring")) %>%
  left_join(richness)

H.rich %>% 
  filter(GROUP %nin% c("INSECT" )) %>%
  select(DOC.1, d2HVSMOW, GROUP) %>%
  na.omit() %>%
  ggplot(aes(x = DOC.1, y = d2HVSMOW, col = GROUP)) + 
  geom_point()+ 
  geom_smooth(method = lm, se = F)

H.rich %>% 
  filter(GROUP %in% c("INSECT" )) %>%
  select(DOC.1, d2HVSMOW, GROUP, MONTH) %>%
  na.omit() %>%
  #{ cor.test(.$d2HVSMOW, .$DOC.1) }
  ggplot(aes(x = DOC.1, y = d2HVSMOW, col = MONTH)) + 
  geom_point() + 
  geom_smooth(method = lm, se = F)

H.rich %>% 
  filter(GROUP %in% c("INSECT" )) %>%
  select(DOC.1, Volume, d2HVSMOW, GROUP, MONTH) %>%
  na.omit() %>%
  { cor.test(.$d2HVSMOW, .$Volume) }
  ggplot(aes(x = Volume, y = d2HVSMOW, col = MONTH)) + 
  geom_point() + 
  geom_smooth(method = lm, se = F)

## Volume & DOC.1 look significant and so does
# Increasing DOC suggests depletion of 2H
