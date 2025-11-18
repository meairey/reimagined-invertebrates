library(dplyr)
library(vegan)
library(tidyverse)
library(ggrepel)
library(wesanderson)
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

## Dart Accumulation in spring ------

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

species_accumulation.dtl <- specaccum(dtl.acc) ## accumulation

## Moss Accumulation in spring ----
msl.acc = dat %>% filter(WATER == "MSL", MONTH == "SPRING") %>%
  select(SITE_N, FAMILY) %>%
  mutate(presence = 1) %>%
  unique() %>%
  pivot_wider(values_from = presence, names_from = FAMILY)%>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  arrange(SITE_N) %>%
  select(-SITE_N)

species_accumulation.msl <- specaccum(msl.acc) ## accumulation

## Plot both curves

zero = data.frame(sites = c(0,0), richness = c(0,0), sd = c(0,0), water = c("Dart Lake", "Moss Lake"))

acc = data.frame(sites = c(species_accumulation.dtl$sites, species_accumulation.msl$sites),
                 richness = c(species_accumulation.dtl$richness, species_accumulation.msl$richness),
                 sd = c(species_accumulation.dtl$sd, species_accumulation.msl$sd),
                 water = c(rep("Dart Lake", 4), rep("Moss Lake",4))) %>%
  rbind(zero)


## Figure S1 --------
# Final figure for supplemental section on species accumulation
ggplot(acc, aes(x = sites, y = richness, col = water, ymin = richness - sd, ymax = richness + sd)) + 
  geom_line()  +
  theme_minimal(base_size = 18) + 
  geom_pointrange() +
  scale_color_manual("Water", values = wes_palette("Darjeeling1", 2)) +
  xlab("# of Sites") + ylab("Expected Species Richness")

# percent of each assemblage sampled by two sites
acc %>% 
  filter(sites %in% c(2, 4)) %>%
  group_by(water) %>%
  summarize(percent = min(richness) / max(richness))

