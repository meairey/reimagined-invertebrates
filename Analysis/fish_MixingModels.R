library(simmr)
library(tidyverse)
## Load in the baseline dataframe that is generated in simmr.R
load(file = "Data/RData/baselines.RData")
data.iso = read.csv("../1.clean_isotope/iso_measurement.csv")
sample = read.csv("../1.clean_isotope/isotope_sample.csv")
## Create a mixture data frame that only includes the fish
fish_mixtures =  data.iso %>% left_join(sample, by = c("ISO_YSAMP_N")) %>%
  left_join(taxon_frame) %>%
  filter(WATER != "LML") %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  filter(CATEGORY == "VERT") %>% ## Only inlcludes fish
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
  select(community, group, D13C, D15N) %>%
  na.omit() 


## Filter the baselines data frame to include only the communities where we caught fish
fish_baselines = baselines %>%
  filter(community %in% fish_mixtures$community)

## Generate the simmr outputs
output.list = list() ## List to put results in

for(i in 1:length(unique(fish_mixtures$community))){
  test.mix = fish_mixtures %>% filter(community == i) %>%
    select(group,D13C, D15N) %>%
    na.omit() %>% ungroup() 
  test.base = fish_baselines %>% filter(community ==i) %>% ungroup() %>%
    na.omit()
  
  if(test.base$community %>% length() > 0){
    ## without groups 
    simmr_in <- simmr_load(
      mixtures =  as.matrix(test.mix[,2:3]),
      source_names = test.base$group.name ,
      source_means = test.base %>% select(mean_c, mean_n),
      source_sds = test.base %>% select(sd_c, sd_n)
    )
    
    simmr_out = simmr_mcmc(simmr_in)
    
    output.list[[i]] = as.data.frame(simmr_out$output$`1`$BUGSoutput$summary) %>%
      mutate(community.name = test.base$community.name%>% unique() ) %>%
      rownames_to_column(var = "rowname")
    
    
    
    plot(simmr_out, type = "boxplot", title = test.base$community.name%>% unique())
  }else{
    output.list[[i]] = NA
  }
}



output.list.clean = output.list

simmr.output = Reduce(full_join, output.list.clean)


simmr.output %>% filter(rowname %in% c("LEAF", "ZOOP", "PERI")) %>%
  select(rowname, mean, community.name) %>%
  ggplot(aes(x = community.name, y = mean, fill = rowname)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  coord_flip()


simmr.output %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  ggplot(aes(x = sechi.depth, y = mean, col =  rowname)) + 
  geom_point() + 
  geom_smooth(method = "lm")




zoop.test = simmr.output %>%
  filter(rowname %in% c("ZOOP", "LEAF", "PERI")) %>%
  group_by(community.name) %>%
  mutate(groups = length(unique(rowname))) %>%
  #filter(groups == 3)  %>% 
  left_join(richness, by = c("community.name")) %>%
  filter(rowname %in% "ZOOP")

zoop.test %>%
  ggplot(aes(x = sechi.depth, y = mean)) + 
  geom_point() + 
  geom_smooth(method = "lm")


cor.test(zoop.test$mean, zoop.test$sechi.depth, method = "spearman")
lm(data = zoop.test, mean ~ thermo_depth) %>% summary()

## Do fish fall in within the range of resources


fish.raw =  data.iso %>% left_join(sample, by = c("ISO_YSAMP_N")) %>%
  left_join(taxon_frame) %>%
  filter(WATER != "LML") %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  filter(CATEGORY == "VERT")

fish.waters = (fish.raw$WATER %>% unique())[-1]


inverts = data.iso %>% left_join(sample, by = c("ISO_YSAMP_N")) %>%
  left_join(taxon_frame) %>%
  filter(WATER != "LML") %>%
  mutate(D15N = as.numeric(D15N),
         D13C = as.numeric(D13C)) %>%
  #filter(CATEGORY == "INVERT") %>%
  filter(WATER %in% fish.waters)


## Using SD of each variable
inverts %>% 
  group_by(WATER, CATEGORY, ITEM_N) %>%
  summarize(D15N = mean(D15N, na.rm = T), 
            D13C = mean(D13C, na.rm = T)) %>%
  ungroup() %>%
  group_by() %>% 
  group_by(WATER, CATEGORY) %>%
  summarize(mean.d15N = mean(D15N,na.rm = T),
            sd.d15N = sd(D15N,na.rm = T),
            mean.d13C = mean(D13C,na.rm = T),
            sd.d13C = sd(D13C,na.rm = T)) %>%
  ggplot(aes(col = CATEGORY)) + 
  geom_point(aes(x = mean.d13C, y = mean.d15N, col = CATEGORY)) + 
  geom_pointrange(aes(x = mean.d13C, xmin = mean.d13C - sd.d13C, xmax = mean.d13C + sd.d13C,
                      y = mean.d15N)) + 
  geom_pointrange(aes(x = mean.d13C, y = mean.d15N, ymin = mean.d15N - sd.d15N, ymax = mean.d15N + sd.d15N)) + 
  facet_wrap(~WATER)

## Using total range of each category

inverts %>% 
  group_by(WATER, CATEGORY, ITEM_N) %>%
  summarize(D15N = mean(D15N, na.rm = T), 
            D13C = mean(D13C, na.rm = T)) %>%
  ungroup() %>%
  group_by() %>% 
  group_by(WATER, CATEGORY) %>%
  summarize(mean.d15N = mean(D15N,na.rm = T),
            max.d15N = max(D15N,na.rm = T),
            min.d15N = min(D15N, na.rm = T),
            mean.d13C = mean(D13C,na.rm = T),
            max.d13C = max(D13C,na.rm = T), 
            min.d13C = min(D13C, na.rm = T)) %>%
  ggplot(aes(col = CATEGORY)) + 
  geom_point(aes(x = mean.d13C, y = mean.d15N, col = CATEGORY)) + 
  geom_pointrange(aes(x = mean.d13C, xmin = min.d13C, xmax = max.d13C,
                      y = mean.d15N)) + 
  geom_pointrange(aes(x = mean.d13C, y = mean.d15N, ymin =min.d15N , ymax =max.d15N)) + 
  facet_wrap(~WATER, ncol = 4) + 
  theme_minimal(base_size = 14) + 
  xlab("d13C") + ylab("d15N") + 
  scale_color_manual("Category", values = wes_palette("Moonrise2", type = "discrete", n = 4),
                                labels = c("Periphyton", "Benthic Inv + Zoop",
                                           "Terrestrial Leaf", "Fish"))


inverts %>% 
  filter(WATER == "COM", CATEGORY == "VERT")










