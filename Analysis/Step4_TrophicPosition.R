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

trophic_position.graph = trophic.position %>% ## Tall version
  left_join(cluster_chem) %>%
  group_by(FAMILY) %>% 
  mutate(count = n()) %>%
  filter(count > 10) %>%
  ggplot(aes(y = factor(FAMILY, rev(TP.family$FAMILY)),
             x = trophic_position_weighted,
             fill = factor(ORDER, TP.order$ORDER))) + 
  geom_boxplot(alpha = .8) +
  geom_point(alpha = .3) + 
  theme_minimal(base_size = 14) +
  scale_fill_manual("Order", values = wes_palette("Darjeeling1", type = "continuous", n = 6)) +
  xlab("Trophic Position") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") + 
  scale_x_log10()


ggsave(plot = trophic_position.graph, file = "Graphics/Figures/Figure5A_TrophicPosition.pdf", width = 6, height = 6)

## Trophic position table to include all the other taxa I weeded out of main figure for ease of visualization

trophic.position %>%
  group_by(ORDER, FAMILY) %>%
  summarize(mean.TP = mean(trophic_position_weighted), 
            sd.TP = sd(trophic_position_weighted), 
            n = n()) %>%
  filter(n > 3) %>% 
  rename("Family" = "FAMILY", 
         "Mean TP" = "mean.TP", 
         "SD TP" = "sd.TP", 
         "N" = "n") -> TP.table


write.csv(TP.table, "Graphics/Tables/Trophic_Position.Table.csv", 
          row.names = F)

#### Figure S5 ---------------
## Seasonal affect of trophic position
## Diamond points (family means across waterbodies)
family_means.tp = trophic.position %>%

  separate(community.name, into = c("water", "season")) %>% 
  group_by(FAMILY, water) %>%
  mutate(season.distinct = length(unique((season)))) %>%
  filter(season.distinct > 1) %>%
  ungroup() %>%
  group_by(water, season, FAMILY)  %>%
  summarize(count = n(),
            mean_tp = mean(trophic_position_weighted)) %>% ## Average for each family per waterbody
  filter(count > 3) %>%
  ungroup() %>%
  group_by(season, FAMILY) %>%
  summarize(mean_overall.tp = mean(mean_tp)) %>% ## Average for each family across all waterbodies
  pivot_wider(names_from = season, values_from = mean_overall.tp) %>% 
  na.omit() %>%
  mutate(mean_diff = spring - fall)

mean_diff.tp = mean(family_means.tp$mean_diff)

## Small points (family means within each water body)
family_water.TP = trophic.position %>% 
 
  separate(community.name, into = c("water", "season")) %>%
  group_by(water, season, FAMILY) %>%
  summarize(count = n(), 
            mean_tp = mean(trophic_position_weighted)) %>%
  filter(count > 3) %>%
  ungroup() %>% 
  select(-count) %>%
  pivot_wider(names_from = season, values_from = mean_tp) %>%
  na.omit() %>%
  mutate(mean_diff = spring - fall)

# New plot

ggplot() + 
  theme_minimal(base_size = 14) +
  geom_point(data = family_water.TP, aes(x = mean_diff, y = FAMILY), col = "gray") +
  geom_point(data = family_means.tp, aes(x = mean_diff, y = FAMILY), size = 6, shape = "diamond") +
  geom_vline(aes(xintercept = 0)) +
  geom_vline(aes(xintercept = mean_diff.tp), lty = "dashed") +
  theme(axis.title.y = element_blank()) +
  xlab("Seasonal Difference in TP") -> FigureS5
ggsave(file = "Graphics/Figures/Figure_S5.pdf", plot = FigureS5, width = 7, height = 5, units = "in")

## Plot
ggplot() +
  theme_minimal(base_size = 14) +
  geom_point(data = family_water.TP, aes(x = spring, y = fall, col = FAMILY), alpha = .3) +
  geom_point(data = family_means.tp, aes(x = spring, y = fall, col = FAMILY), shape = "diamond", size = 4 ) +
  geom_abline(slope = 1, intercept = 0) + 
  scale_color_manual("Family", values = wes_palette("Darjeeling1", type = "continuous", n = 11)) + 
  xlab("Spring TP") + ylab("Fall TP") -> FigureS5



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



trophic_chemistry = trophic.position %>% 
  left_join(read.csv("Data/CSVs/richness_update.csv")) %>%
  separate(community.name, into = c("WATER", "season")) %>%
  select(-temp_do, -DOC_update_text, -SurficialGeology, -Lake.Type) %>%
  group_by(FAMILY) %>%
  mutate(total = n()) %>%
  filter(total > 10) %>%
  select(-total) 


library(lmtest) 

vars = colnames(trophic_chemistry)[15:26]

i = 11
for(i in 1:length(vars)){
 
  form = as.formula(paste0("trophic_position_weighted ~ ", vars[i], "* FAMILY + (1|WATER)"))
  lmer.run = lmerTest::lmer(form, data = trophic_chemistry) %>%
  summary()
  if(lmer.run$coefficients[2,5] < .05){
    print(vars[i])
    print(lmer.run$coefficients[2,])
  }

  
}
library(emmeans)
## Significant results include DOC_update, depth.5gmL, temp5mgL, and secchi depth
ffgs = read.csv("Data/CSVs/FFGS.csv") %>%
  select(FAMILY, FG) %>%
  unique() %>%
  arrange(FG) %>%
  mutate(FG_factor = as.numeric(as.factor(FG)))
ffgs
# Secchi depth

lmer.secchi = lmerTest::lmer(trophic_position_weighted ~ sechi.depth * FAMILY + (1|WATER), data = trophic_chemistry)

fam_secchi <- emtrends(lmer.secchi, 
                       specs = "FAMILY", 
                       var = "sechi.depth") %>% ##
  summary() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL)) %>%
  na.omit() %>%
  mutate(sig = sign(lower.CL) == sign(upper.CL))  %>%
  rename("trend" = "sechi.depth.trend" ) %>%
  mutate(metric = "Secchi")

# Secchi depth

lmer.DOC = lmerTest::lmer(trophic_position_weighted ~ DOC_update * FAMILY + (1|WATER), data = trophic_chemistry)

fam_DOC <- emtrends(lmer.DOC, 
                       specs = "FAMILY", 
                       var = "DOC_update") %>% ##
  summary() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL)) %>%
  na.omit() %>%
  mutate(sig = sign(lower.CL) == sign(upper.CL))  %>%
  rename("trend" = "DOC_update.trend" ) %>%
  mutate(metric = "DOC")



# DDO5
lmer.DDO5 = lmerTest::lmer(trophic_position_weighted ~ depth.5mgL * FAMILY + (1|WATER), data = trophic_chemistry)

fam_DDO5 <- emtrends(lmer.DDO5, 
                       specs = "FAMILY", 
                       var = "depth.5mgL") %>% ##
  summary() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL)) %>%
  na.omit() %>%
  mutate(sig = sign(lower.CL) == sign(upper.CL))  %>%
  rename("trend" = "depth.5mgL.trend" ) %>%
  mutate(metric = "DDO5")


## TDO5
lmer.TDO5 = lmerTest::lmer(trophic_position_weighted ~ temp.5mgL * FAMILY + (1|WATER), data = trophic_chemistry)

fam_TDO5 <- emtrends(lmer.TDO5, 
                       specs = "FAMILY", 
                       var = "temp.5mgL") %>% ##
  summary() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL)) %>%
  na.omit() %>%
  mutate(sig = sign(lower.CL) == sign(upper.CL)) %>%
  rename("trend" = "temp.5mgL.trend" ) %>%
  mutate(metric = "TDO5")





trophic_position_lmer.graph = rbind(fam_secchi, fam_DOC) %>% 
  rbind(fam_DDO5, fam_TDO5)  %>%
  left_join(taxon_frame %>% 
              select(FAMILY,ORDER) %>%
              unique()) %>%
  arrange()
  mutate(metric = factor(metric, levels = c("DOC", "Secchi", "DDO5", "TDO5"))) %>%
  ggplot(aes(x = as.numeric(trend), y = factor(FAMILY, rev(TP.family$FAMILY)), color = factor(ORDER, TP.order$ORDER))) +
  geom_pointrange(aes(xmin = as.numeric(lower.CL), xmax = as.numeric(upper.CL),
                       shape = sig)) + 
   geom_vline(aes(xintercept = 0)) +
   theme_minimal(base_size = 14) +
   xlab("Slope 95% CI") +
  scale_shape_manual("Significance", values = c(1, 16), labels = c("Insignificant", "Significant")) +
  scale_color_manual("Order", values = wes_palette("Darjeeling1", type = "continuous", n = 6)) +
   theme(axis.title.y = element_blank(), 
         legend.position = "bottom",
         legend.box = "vertical") + 
  facet_wrap(~metric, scales = "free_x") +
  guides(color = guide_legend(nrow = 2), 
         shape = guide_legend(nrow = 2) ) 

trophic_position_lmer.graph


ggsave(plot = trophic_position_lmer.graph, file = "Graphics/Figures/Figure5B_TrophicPosition.pdf", width = 5.5, height = 7)

library(gridExtra)
grid.arrange(DDO5, TDO5, secchi, DOC)
