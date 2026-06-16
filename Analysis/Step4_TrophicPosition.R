## Setup ----------------------------
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(simmr)
library(tidyverse)
library(lmerTest)
library(emmeans)
# Trophic Position -----------
# Using mixing model results
load(file = "Data/RData/simmr_full.RData")
## Load in isotope data

data.iso = read.csv("Data/CSVs/processed_data.csv")
baselines = read.csv("Data/CSVs/baselines.csv")
load("Data/RData/cluster.mat.RData")
load(file = "Data/RData/community.legend.RData")

baselines %>% 
  group_by(community.name) %>%
  select(community.name,group.name,  mean_c) %>%
  pivot_wider(names_from = group.name, values_from = mean_c) %>% 
  ungroup() %>%
  mutate(dif = ZOOP - LEAF) %>% 
  left_join(cluster.mat %>% 
              unite("community.name", c(WATER, season), sep = ".")) %>% 
  filter(!is.na(cluster)) %>%
  ggplot(aes(x = cluster, y = dif)) + 
  geom_point()



## Taxon frame 

taxon_frame = read.csv("Data/CSVs/taxon_frame.csv") %>% unique()# rename column

## Cluster and chemistry data
cluster_chem = read.csv(file = "Data/CSVs/chemistry_seasonal.csv") %>%
  unite("community.name", c(WATER, season), sep = ".") %>%
  #unite("ID", WATER, season, sep = ".") %>%
  left_join(community.legend) %>%
  left_join(cluster.mat %>%
              unite(ID, c("WATER", "season"), sep = "."))
## Mixtures for SIMMR model (invertebrate observations)
mixtures =  data.iso %>% 
  left_join(taxon_frame %>% select(ORDER, FAMILY, TAXON) %>% unique()) %>%
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
  filter(group.name %nin% c("LEAF", "PERI", "ZOOP", "FISH")) 




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
            # fill = factor(ORDER, TP.order$ORDER)
            )) + 
  geom_point(alpha = .2) + 
  geom_boxplot(alpha = .5, fill = "lightgray") +
  
  theme_minimal(base_size = 12) +
  scale_fill_manual("Order", values = wes_palette("Darjeeling1", type = "continuous", n = 6)) +
  xlab("Trophic Position") +
  theme(axis.title.y = element_blank(),
        legend.position = "right") + 
  scale_x_log10()
trophic_position.graph

#ggsave(plot = trophic_position.graph, file = "Graphics/Figures/Figure5A_TrophicPosition.pdf", width = 3.5, height = 4.5)

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


#write.csv(TP.table, "Graphics/Tables/Trophic_Position.Table.csv", row.names = F)


## Modelling seasonal effect


season.trophic = lmer(trophic_position_weighted ~ Season*FAMILY + (1|Water),
                      data = trophic.position %>% 
                        separate(community.name, into = c("Water", "Season")) %>%
                        group_by(Water, FAMILY) %>%
                        mutate(season_distinct = length(unique(Season))) %>%
                        filter(season_distinct > 1) %>%
                        select(-season_distinct))
summary(season.trophic)
fe = summary(season.trophic)$coefficients %>% 
  as.data.frame() %>% 
  rownames_to_column("term") %>% 
 filter(term %in% c("(Intercept)", "Seasonspring")) 

season.emmeans = emmeans(season.trophic, ~ Season | FAMILY) 
seasons.contrast = contrast(season.emmeans, method = "pairwise")
# get confidence intervals
season.contrast.ci = confint(seasons.contrast) %>% 
  as.matrix() %>% as.data.frame()
 
#### Figure S5 ---------------
## Seasonal affect of trophic position

season.emmean.graph = season.contrast.ci %>% 
  mutate(estimate = as.numeric(estimate), 
         lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL), 
         sig  = sign(lower.CL) == sign(upper.CL)) %>%
  ggplot(aes(x = estimate, y = FAMILY, col = sig)) + 
  geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL)) +
  theme_minimal(base_size = 12) + 
  geom_vline(xintercept = 0, lty = "dashed")  +
  scale_color_manual("Significance", values = c("gray", "black"),
                     labels = c("Insignificant", "Significant")) + 
  theme(axis.title.y = element_blank()) +
  xlab("Seasonal TP Difference (95% CI)")

#ggsave(season.emmean.graph, file = "Graphics/Figures/Supp_Figure_5.pdf", width = 5, height = 4)

# Note - no real difference in trophic position between clusters (scaled or unscaled) 12/8/25
## No real difference in trophic position across chemistry 
## I've removed these from the scripts - but can go back in version control and look



trophic_chemistry = trophic.position %>% 
  left_join(cluster_chem) %>%
  separate(community.name, into = c("WATER", "season")) %>%
  group_by(FAMILY) %>%
  mutate(total = n()) %>%
  filter(total > 10) %>%
  select(-total) %>% 
  select(Pond_num, WATER, season,  FAMILY, GROUP, cluster, community, Lake, everything()) %>%
  mutate(across(TDO5:Volume,
         ~ as.numeric(scale(.x))))



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

## Significant results include DOC_update, depth.5gmL, temp5mgL, and secchi depth

## Checking between lm and lmer
m.fixed = lm(trophic_position_weighted ~ sechi.depth * FAMILY,
             data = trophic_chemistry)

m.random <- lmerTest::lmer(trophic_position_weighted ~ sechi.depth * FAMILY + (1|WATER),
                           data = trophic_chemistry)
AIC(m.fixed,m.random)

m.nointer = lmerTest::lmer(trophic_position_weighted ~ sechi.depth + (1|WATER),
                           data = trophic_chemistry)

m.inter <- lmerTest::lmer(trophic_position_weighted ~ sechi.depth * FAMILY + (1|WATER),
                           data = trophic_chemistry)
AIC(m.nointer,m.inter)


## Combined model

lmerTest::lmer(trophic_position_weighted ~ DOC_update *  FAMILY + DDO5  + (1|WATER), data = trophic_chemistry) %>%
  summary()

trophic_chemistry %>%
  ungroup() %>%
  select(DOC_update, sechi.depth, DDO5, TDO5) %>%
  cor()



# Secchi depth

(trophic_chemistry$FAMILY) %>% unique()

lmer.secchi = lmerTest::lmer(trophic_position_weighted ~ sechi.depth * FAMILY + (1|WATER), data = trophic_chemistry)

lmer.secchi.summary = summary(lmer.secchi)$coefficients %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "effects") %>%
  filter(effects %in% c("(Intercept)", "sechi.depth")) %>%
  mutate(metric = 'Secchi') %>%
  cbind(data.frame(summary(lmer.secchi)$varcor)%>% 
          select(-var1, -var2))

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

# DOC depth

lmer.DOC = lmerTest::lmer(trophic_position_weighted ~ DOC_update * FAMILY + (1|WATER), data = trophic_chemistry)
lmer.DOC %>%
lmer.DOC.summary = summary(lmer.DOC)$coefficients %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "effects") %>%
  filter(effects %in% c("(Intercept)", "DOC_update")) %>%
  mutate(metric = 'DOC') %>%
  cbind(data.frame(summary(lmer.DOC)$varcor)%>% 
          select(-var1, -var2))

emmeans(lmer.DOC, 
                       specs = "FAMILY", 
                       var = "DOC_update")

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

lmer.ddo5.summary = summary(lmer.DDO5)$coefficients %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "effects") %>%
  filter(effects %in% c("(Intercept)", "depth.5mgL")) %>%
  mutate(metric = 'ddo5') %>%
  cbind(data.frame(summary(lmer.DDO5)$varcor)%>% 
          select(-var1, -var2))
  
fam_DDO5 <- emtrends(lmer.DDO5, specs = "FAMILY", var = "depth.5mgL") %>% ##
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
lmer.tdo5.summary = summary(lmer.TDO5)$coefficients %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "effects") %>%
  filter(effects %in% c("(Intercept)", "temp.5mgL")) %>%
  mutate(metric = 'tdo5') %>%
  cbind(data.frame(summary(lmer.TDO5)$varcor)%>% 
          select(-var1, -var2))

fam_TDO5 <- emtrends(lmer.TDO5, specs = "FAMILY", var = "temp.5mgL") %>% ##
  summary() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(lower.CL = as.numeric(lower.CL), 
         upper.CL = as.numeric(upper.CL)) %>%
  na.omit() %>%
  mutate(sig = sign(lower.CL) == sign(upper.CL)) %>%
  rename("trend" = "temp.5mgL.trend" ) %>%
  mutate(metric = "TDO5")

## Summary Table of LMER results for the four models

summary.table.TP = rbind(lmer.secchi.summary, lmer.DOC.summary, lmer.tdo5.summary, lmer.ddo5.summary)
#write.csv(summary.table.TP, file = "Graphics/Tables/TrophicPosition_LMER.csv")



### Full graph of emmeans results
trophic_position_lmer.graph = rbind(fam_secchi, fam_DOC) %>% 
  rbind(fam_DDO5, fam_TDO5)  %>%
  left_join(taxon_frame %>% 
              select(FAMILY,ORDER) %>%
              unique()) %>%
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
         shape = guide_legend(nrow = 2)) 
trophic_position_lmer.graph

## Figure 4 - Simplified graph of emtrends results ---------------- 

trophic_position_lmer.graph = rbind(fam_secchi, fam_DOC) %>% 
  rbind(fam_DDO5, fam_TDO5)  %>%
  left_join(taxon_frame %>% 
              select(FAMILY,ORDER) %>%
              unique()) %>%
  mutate(metric = factor(metric, levels = c("DOC", "Secchi", "DDO5", "TDO5")), 
         trend = as.numeric(trend)) %>% 
  #filter(sig == "TRUE") %>%
  ggplot(aes(x = metric,
             y = factor(FAMILY, rev(TP.family$FAMILY)), 
             fill = trend,
             label = round(trend, digits = 2),
              alpha = sig)) + 
  geom_tile() + 
  geom_text(col = "black", size = 5) +
  scale_fill_gradientn("Slope",colors = c(wes_palette("Zissou1", n = 5, type = "discrete")[1], 
                                    wes_palette("Zissou1", n = 5, type = "discrete")[3], 
                                    wes_palette("Zissou1", n = 5, type = "discrete")[5])) +
  theme_minimal(base_size = 12) + 
  theme(axis.title.y = element_blank()) +
  xlab("LMER Results") +
  scale_alpha_manual("Significance", labels = c(0,1), values = c(0,1)) +
  guides(alpha = "none")

trophic_position_lmer.graph
#ggsave(plot = trophic_position_lmer.graph, file = "Graphics/Figures/Figure5B_TrophicPosition.pdf", width = 5, height = 4.5)


## Supplemental Table 4 -------------------------
#Summary table of emtrends from all four LMER models

tab_s4 = rbind(fam_DOC, fam_secchi, fam_DDO5, fam_TDO5) %>%
  mutate(across(trend:df, as.numeric)) %>% 
  select(metric, FAMILY, trend, SE, lower.CL, upper.CL, df) %>% 
  mutate(across(where(is.numeric), round, digits = 3)) %>%
  arrange(metric, FAMILY) %>%
  rename(Family = FAMILY)
#write.csv(tab_s4, file = "Graphics/Tables/Supplemental_Table_4.csv", row.names =  F)

## Comparison between taxa using emmeans
DOC_emmeans = emmeans(lmer.DOC, specs = "FAMILY", var = "DOC_update")
pairs(DOC_emmeans) %>% 
  as.data.frame() %>%
  mutate(across(estimate:p.value, round, digits = 3)) %>%
  mutate(p.value = as.numeric(p.value)) %>%
#  filter(p.value < .05) %>%
  mutate(p.value = case_when(p.value < .001 ~ "< .001", p.value >= .001 ~ as.character(p.value)))
