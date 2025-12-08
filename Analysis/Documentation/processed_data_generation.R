### Data source documentation raw to processed CSV files

`%nin%` = Negate("%in%")
set.seed(123)
## Join the sampling and measurement data files

## Load in taxon data
taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column


sample = read.csv("../1.clean_isotope/isotope_sample.csv") ## Sampling event information
data.iso=  read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  left_join(sample, by = c("ISO_YSAMP_N")) %>%
  ## lipid correction factor
  mutate(D13C = case_when(CATEGORY %nin% c("ALGA", "PLANT") ~ as.numeric(D13C) - 3.32 + .99 * (PER_C / PER_N),
                          CATEGORY %in% c("ALGA", "PLANT") ~ as.numeric(D13C))) %>% 
  filter(GROUP != "FISH", # Remove fish f
         !(WATER == "LML" & DAY_N < 230)) %>% # Remove extraneous sampling events 
  filter(TAXON %nin% c("LEACH", "INSECT", "degraded", 
                       "CHECK","SIL", "BIVALVE", 
                       "terrestrial", "unidentified", "odonata",
                       "gastropoda", "mollusca", "anuran", "hydrachnidia", 
                       "COLEOPTERA","amphipoda", "trichoptera", "hirudinea")) ## Filter out samples that weren't ID'd to family level


write.csv(data.iso, file = "Data/CSVs/processed_data.csv")


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



write.csv(data, file = "Data/CSVs/siber_data.csv")


community.legend = data %>% select(community, community.name) %>%
  unique() ## 


## Create the data frame for the baselines
baselines =  data.iso %>% 
  left_join(taxon_frame %>% 
              select(TAXON, FAMILY, ORDER) %>% 
              unique()) %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  mutate(GROUP = str_replace(GROUP, "CLAM","BIVALVE"),
         GROUP = str_replace(GROUP,"MUSSEL", "BIVALVE"))  %>%
  mutate(graph_id = case_when(is.na(ORDER) == T ~ GROUP, 
                              is.na(ORDER) == F ~ ORDER)) %>%
  mutate(season = case_when(MONTH < 8 ~ "spring", MONTH > 7 ~ "fall")) %>%
  unite("community.name", c(WATER, season), sep = ".") %>%
  mutate(group.name = case_when(is.na(FAMILY)== T ~ GROUP, 
                                is.na(FAMILY)==F ~ FAMILY)) %>%
  ungroup() %>%
  arrange(community.name) %>%
  mutate(community = as.numeric(as.factor(community.name))) %>%
  arrange(group.name) %>%
  mutate(group = as.numeric(as.factor(group.name))) %>%
  arrange(community, group) %>%
  as.data.frame() %>%
  filter(group.name %in% c("LEAF", "PERI", "ZOOP")) %>%
  group_by(community, community.name, group.name) %>%
    summarize(mean_c = mean(D13C, na.rm = T), 
              mean_n = mean(D15N, na.rm = T), 
              sd_c = sd(D13C, na.rm = T), 
              sd_n = sd(D15N, na.rm = T)) %>%
  na.omit() %>%
  ungroup() %>%
  complete(community.name, group.name) %>%
  select(community, community.name, group.name, mean_c, mean_n, sd_c, sd_n)
## Adding in zooplankton samples that average across both spring and fall or use spring and fall data
## adding the averaged information for HTL zoop in bc not enough 
HTL.zoop = data.iso %>% 
  filter(GROUP == "ZOOP", 
          grepl( "HTL", ISO_YSAMP_N)) %>%
  mutate(SEASON = c("fall", "spring", "spring"))

baselines[24, 1] = community.legend[which(community.legend$community.name == "HTL.fall"),1]; baselines[24, 2] = "HTL.fall"; baselines[24,3] = "ZOOP";
baselines[24,4] = (HTL.zoop %>% filter(SEASON == "fall"))$D13C;
baselines[24,5] = (HTL.zoop %>% filter(SEASON == "fall"))$D15N;
baselines[24,6] = sd(HTL.zoop$D13C); 

baselines[24,7] = sd(HTL.zoop$D15N)


#baselines[24, 1] = 8; baselines[24, 2] = "HTL.fall"; baselines[24,3] = "ZOOP"
#baselines[24,4] = -33.94; baselines[24,5] = 6.676667; baselines[24,6] = 1.105908
#baselines[24,7] = 2.206362


## adding the averaged information for SEL zoop in bc not enough 

sel.zoop = data.iso %>% 
  filter(GROUP == "ZOOP", 
          grepl( "SEL", ISO_YSAMP_N)) %>%
  mutate(SEASON = c("fall", "spring", "spring"))


baselines[51, 1] =  community.legend[which(community.legend$community.name == "SEL.fall"),1]; baselines[51, 2] = "SEL.fall"; baselines[51,3] = "ZOOP";
baselines[51,4] = (sel.zoop %>% filter(SEASON == "fall"))$D13C;
baselines[51,5] = (sel.zoop %>% filter(SEASON == "fall"))$D15N;
baselines[51,6] = sd(sel.zoop$D13C); 

baselines[51,7] = sd(sel.zoop$D15N)


baselines = baselines %>% 
  left_join(community.legend) ## community legend is earlier in the script


## saving the baseline frame
write.csv(baselines, file = "Data/CSVs/baselines.csv")


### Writing CSV table for baselines (supplement) 

write.csv(baselines %>%
  pivot_wider(
    names_from = group.name,
    values_from = c(mean_c, sd_c, mean_n, sd_n) 
  ) %>%
    select(community.name, 
           mean_c_LEAF, sd_c_LEAF, mean_n_LEAF, sd_n_LEAF,
           mean_c_PERI, sd_c_PERI, mean_n_PERI, sd_n_PERI,
           mean_c_ZOOP, sd_c_ZOOP, mean_n_ZOOP, sd_n_ZOOP),file = "Data/CSVs/TableS1_baselines.csv", 
  row.names = F)
