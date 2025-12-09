# Trophic Position -----------
# Using mixing model results
load(file = "Data/RData/simmr_full.RData")

## Setup ----------------------------
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)


## Pivoting and simplifying SIMMR results

baseline.mixtures = simmr.full %>% 
  select(rowname, mean, taxa, community) %>%
  filter(rowname %in% c("PERI", "ZOOP", "LEAF")) %>%
  unique() %>%
  pivot_wider(names_from = rowname, values_from = mean) %>%
  rename(PERI.mix = PERI, 
         ZOOP.mix = ZOOP, 
         LEAF.mix = LEAF)

baselines.wide = baselines %>%
  ungroup() %>%
  select(community.name, group.name, mean_n) %>%
  pivot_wider(names_from = group.name, values_from = mean_n)

## Trophic position calculations
#### Includes both the weighted and unweighted versions of the calculation
trophic.position = mixtures %>%
  select(community.name, GROUP, FAMILY, D15N, D13C) %>%
  left_join(baseline.mixtures, by = c( "community.name" ="community", "FAMILY" = "taxa"  )) %>%
  left_join(baselines.wide, by = c("community.name")) %>%
  mutate(baseline_N = case_when(is.na(PERI.mix) == F ~ ZOOP * ZOOP.mix + PERI * PERI.mix + LEAF * LEAF.mix,
                                is.na(PERI.mix) == T ~  ZOOP * ZOOP.mix +LEAF * LEAF.mix)) %>%
  mutate(baseline_position = case_when(is.na(PERI.mix) == F ~ 2 * ZOOP.mix + 1 * PERI.mix + 1 * LEAF.mix,
                                        is.na(PERI.mix) == T ~ 2 * ZOOP.mix + 1 * LEAF.mix)) %>%
  mutate(trophic_position_weighted = ((D15N - baseline_N)/3.14) + baseline_position) %>%
  mutate(trophic_position = case_when(is.na(PERI.mix) == F ~ ((D15N - PERI)/3.14) + 1,
                                      is.na(PERI.mix) == T ~ ((D15N - LEAF) / 3.14) + 1)) %>%
  select(community.name, GROUP, FAMILY, D15N, D13C, trophic_position_weighted, trophic_position) %>%
  na.omit() %>%
  filter(FAMILY %nin% c("unidentified", "terrestrial", "degraded")) %>%
  left_join(taxon_frame %>% select(ORDER, FAMILY) %>% unique()) %>%
  ungroup() %>%
  group_by(FAMILY) %>%
  mutate(mean_TP = mean(trophic_position_weighted))

#save(trophic.position, file = "Data/RData/trophic.position.RData")
#load(file = "Data/RData/trophic.position.RData")

## These two objects are used for graphing order
TP.family = trophic.position %>% ## Order y axis by family
  ungroup() %>%
  select(ORDER, FAMILY) %>%
  arrange(ORDER) %>% select(FAMILY) %>% unique()


TP.order = trophic.position %>% ## Cluster familes in the same order together
  ungroup() %>%
  select(ORDER, FAMILY) %>%
  arrange(ORDER) %>% select(ORDER) %>% unique()

### Visualization ------------------

#### Figure 5 ------------
## Plot the trophic position of the taxa using the weighted estimation from SIMMR
trophic.position %>%
  left_join(cluster_chem) %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 5) %>%
  ggplot(aes(y = factor(FAMILY, rev(TP.family$FAMILY)),
             x = trophic_position_weighted,
             fill = factor(ORDER, TP.order$ORDER))) + 
  geom_boxplot() +
  geom_point(alpha = .2) + 
  theme_minimal(base_size = 20) +
  scale_fill_manual("Order", values = wes_palette("Darjeeling1", type = "continuous", n = 9)) +
  xlab("Trophic Position") +
  theme(axis.title.y = element_blank()) +
  scale_x_log10()



#### Figure S5 ---------------
## Seasonal affect of trophic position
## Diamond points (family means across waterbodies)
family_means.tp = trophic.position %>% 
  separate(community.name, into = c("water", "season")) %>%
  group_by(water, season, FAMILY)  %>%
  summarize(count = n(),
            mean_tp = mean(trophic_position_weighted)) %>% ## Average for each family per waterbody
  filter(count > 3) %>%
  ungroup() %>%
  group_by(season, FAMILY) %>%
  summarize(mean_overall.tp = mean(mean_tp)) %>% ## Average for each family across all waterbodies
  pivot_wider(names_from = season, values_from = mean_overall.tp) %>% 
  na.omit()

## Small points (family means within each waterbody)
family_water.TP = trophic.position %>% 
  separate(community.name, into = c("water", "season")) %>%
  group_by(water, season, FAMILY) %>%
  summarize(count = n(), 
            mean_tp = mean(trophic_position_weighted)) %>%
  filter(count > 3) %>%
  ungroup() %>% 
  select(-count) %>%
  pivot_wider(names_from = season, values_from = mean_tp) %>%
  na.omit()
## Plot
ggplot() +
  theme_minimal(base_size = 14) +
  geom_point(data = family_water.TP, aes(x = spring, y = fall, col = FAMILY), alpha = .3) +
  geom_point(data = family_means.tp, aes(x = spring, y = fall, col = FAMILY), shape = "diamond", size = 4 ) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_color_manual("Family", values = wes_palette("Darjeeling1", type = "continuous", n = 11)) + 
  xlab("Spring TP") + ylab("Fall TP")


## Seasonal outliers
df_dist = family_means.tp %>%
  mutate(
    # perpendicular distance to the line y = x
    dist_1to1 = abs(fall - spring) / sqrt(2),

    # sign (positive = above line, negative = below)
    signed_dist = (fall - spring) / sqrt(2)
  )

df_dist

### Comparing trophic clusters

trophic.position %>% 
  left_join(cluster_chem) %>%
  filter(!is.na(cluster)) %>%
  group_by(community.name, cluster, FAMILY) %>%
  summarize(TP = mean(trophic_position_weighted)) %>% # mean TP per family per waterbody
 # filter(FAMILY %in% c("aeshnidae", "gomphidae", "libellulidae", "corduliidae")) %>%
  ggplot(aes(x = as.factor(cluster), y = TP, fill = as.factor(cluster))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_minimal(base_size = 15) + 
  ylab("Trophic Position")  + 
  scale_x_discrete(labels = c(
  "1" = "Shallow\nthermocline", 
  "3" = "Deep\nthermocline", 
  "2" = "Diverse,\ndeep thermocline"
))  +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
   scale_fill_manual(values = wes_palette("Royal2")[c(3,1,5)]) 



# Note - no real difference in trophic position between clusters (scaled or unscaled) 12/8/25
## No real difference in trophic position across chemistry 
## I've removed these from the scripts - but can go back in version control and look


