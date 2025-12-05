
## Load in isotope measurement file
data.iso=  read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  mutate(D13C = case_when(CATEGORY %nin% c("ALGA", "PLANT") ~ as.numeric(D13C) - 3.32 + .99 * (PER_C / PER_N),
                          CATEGORY %in% c("ALGA", "PLANT") ~ as.numeric(D13C))) ## lipid correction factor


sample = read.csv("../1.clean_isotope/isotope_sample.csv")

iso = read.csv("../1.clean_isotope/iso_measurement.csv")  %>%
  left_join(sample, by = "ISO_YSAMP_N") %>%
  mutate(D13C = case_when(CATEGORY %nin% c("ALGA", "PLANT") ~ as.numeric(D13C) - 3.32 + .99 * (PER_C / PER_N),
                          CATEGORY %in% c("ALGA", "PLANT") ~ as.numeric(D13C))) ## lipid correction factor


## Standard deviation and means of replicates 


reps = data.iso %>% group_by(ITEM_N) %>% 
  select(ITEM_N, D13C, D15N, CATEGORY) %>%
  unique() %>%
  mutate(replicates = n()) %>%
  filter(replicates > 3) %>%
  group_by(ITEM_N, replicates) %>%
  select(ITEM_N, replicates, D13C, D15N, CATEGORY) %>%
  unique() %>%
  summarize(mean_C = mean(D13C, na.rm = T), sd_C = sd(D13C, na.rm = T), 
            mean_N = mean(D15N, na.rm = T), sd_N = sd(D15N, na.rm = T)) %>%
  rename(`Sample ID` = ITEM_N, 
         Replicates = replicates, 
         `Mean D13C` = mean_C, 
         `SD D13C` = sd_C,
         "Mean D15N" = mean_N, 
         "SD D15N" = sd_N) %>%
  arrange(`Sample ID`) 


reps$Group = c("Periphyton", 
                     "Creek chub", 
                     "Zooplankton", 
                     "Zooplankton", 
                     "Odonata")
#write.csv(reps, file = "Data/CSVs/replicates.csv")
 