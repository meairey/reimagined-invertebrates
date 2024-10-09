library(dplyr)
library(vegan)
library(tidyverse)
library(ggrepel)
`%nin%` = Negate(`%in%`)

## Accumulation curves --------------------------
## These are for lakes that are sampled 4 sites just in spring
set.seed(123)  # For reproducibility

## Accumulation curves --------------------------
## These are for lakes that are sampled 4 sites just in spring


load(file = "Data/RData/FullData.RData")

dat = full %>% select(MONTH, YSAMP, SITE_N, WATER, ORDER..OR.ABOVE., FAMILY) %>%
  unique() %>%
  filter(FAMILY %nin% c("unidentified", "degraded","terrestrial")) %>%
  mutate(SITE_N = str_replace_all(SITE_N, " ", ""),
         SITE_N = str_replace_all(SITE_N, "Notlisted", "NotListed")) 

## Dart Accumulation in spring...

dtl.acc = dat %>%
  filter(WATER == "DTL", MONTH == "SPRING") %>%
  select(SITE_N, FAMILY) %>%
  mutate(presence = 1, 
         SITE_N = as.numeric(SITE_N)) %>%
  unique() %>%
  pivot_wider(values_from = presence, names_from = FAMILY)%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(SITE_N) %>%
  select(-SITE_N)

dtl.acc

# Assuming you have a species-by-site matrix named 'species_data'
species_accumulation.dtl <- specaccum(dtl.acc)

# Plot the accumulation curve
plot(species_accumulation.dtl, xlab="Number of Sites", ylab = "Cumulative Species")


## Moss Accumulation in spring...
msl.acc = dat %>% filter(WATER == "MSL", MONTH == "SPRING") %>%
  select(SITE_N, FAMILY) %>%
  mutate(presence = 1) %>%
  unique() %>%
  pivot_wider(values_from = presence, names_from = FAMILY)%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(SITE_N) %>%
  select(-SITE_N)

msl.acc

# Assuming you have a species-by-site matrix named 'species_data'
species_accumulation.msl <- specaccum(msl.acc)

# Plot the accumulation curve
plot(species_accumulation.msl, xlab="Number of Sites", ylab="Cumulative Species")


acc = data.frame(sites = c(species_accumulation.dtl$sites, species_accumulation.msl$sites),
                 richness = c(species_accumulation.dtl$richness, species_accumulation.msl$richness),
                 water = c(rep("Dart Lake", 4), rep("Moss Lake",4)))

ggplot(acc, aes(x = sites, y = richness, col = water)) + 
  geom_line()  +
  theme_minimal()


## Accumulation curves by season ------------------
## Counting each site/season combo as a different sampling event



dtl.acc = dat %>% filter(WATER == "DTL") %>%
  select(SITE_N,MONTH, FAMILY) %>%
  mutate(presence = 1) %>%
  unique() %>%
  pivot_wider(values_from = presence, names_from = FAMILY) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(MONTH, SITE_N) %>%
  select(-SITE_N, -MONTH)

dtl.acc

# Assuming you have a species-by-site matrix named 'species_data'
species_accumulation.dtl <- specaccum(dtl.acc)


## Moss Accumulation in spring...
msl.acc = dat %>% filter(WATER == "MSL") %>%
  select(SITE_N, MONTH, FAMILY) %>%
  mutate(presence = 1) %>%
  unique() %>%
  pivot_wider(values_from = presence, names_from = FAMILY)%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(MONTH, SITE_N) %>%
  select(-SITE_N, -MONTH)

msl.acc

# Assuming you have a species-by-site matrix named 'species_data'
species_accumulation.msl <- specaccum(msl.acc)


acc = data.frame(sites = c(species_accumulation.dtl$sites, species_accumulation.msl$sites),
                 richness = c(species_accumulation.dtl$richness, species_accumulation.msl$richness),
                 water = c(rep("Dart Lake", dim(dtl.acc)[1]), rep("Moss Lake",dim(msl.acc)[1])))

ggplot(acc, aes(x = sites, y = richness, col = water)) + 
  geom_line()  +
  theme_minimal()
