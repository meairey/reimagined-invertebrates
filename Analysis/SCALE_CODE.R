library(dplyr)
library(tidyr)
library(tidyverse)
library(gridExtra)
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/AFRP")

s.d = read.csv(file = "MA2276_Code/SCALE/SCALE_MTS_2023.csv") %>% separate(YSAMP, into = c("GEAR", "WATER", "YEAR", "YSAMP_N"))  %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>%## Scale Data 
  filter(ORDER..OR.ABOVE. != "") %>%
  mutate(FAMILY = str_replace(FAMILY, " ", "")) %>%
  mutate(ORDER..OR.ABOVE. = str_replace(ORDER..OR.ABOVE., " ", ""))
## Family~Water inverts
s.d  %>% 
  group_by(WATER, YSAMP_N, ORDER..OR.ABOVE., FAMILY, GENUS) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  ungroup() %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
  filter(ORDER..OR.ABOVE. != "") %>%
  ggplot(aes(x = FAMILY, y = sum, fill =WATER)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~WATER, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90))

# All families per waterbody
s.d  %>% 
  group_by(WATER, ORDER..OR.ABOVE., FAMILY) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  ungroup() %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
  filter(ORDER..OR.ABOVE. != "") %>%
  unite("z", c(WATER,ORDER..OR.ABOVE.,FAMILY),sep = "_", remove = F) %>%
  arrange(z, (sum)) %>%
  mutate(z = factor(z, levels =z)) %>%
  ggplot(.,aes(x = factor(z, levels = z[order(sum)]), y = sum, fill =ORDER..OR.ABOVE.)) + 
 # scale_x_discrete(data = sd, labels =FAMILY) +
  geom_bar(stat = "identity") + 
  facet_wrap(~WATER, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("total inverts") + 
  xlab("") + 
  labs(fill = "Order")



## Most common families per waterbody
s.d  %>% 
  group_by(WATER, ORDER..OR.ABOVE., FAMILY) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  ungroup() %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
  filter(ORDER..OR.ABOVE. != "") %>%
  unite("z", c(WATER,ORDER..OR.ABOVE.,FAMILY),sep = "_", remove = F) %>%
  arrange(z, (sum)) %>%
  filter(sum > 5) %>%
  ungroup() %>%
  group_by(WATER) %>%
  ggplot(.,aes(y = reorder(FAMILY, sum), x = sum, fill =ORDER..OR.ABOVE.)) + 
  # scale_x_discrete(data = sd, labels =FAMILY) +
  geom_bar(stat = "identity") + 
  facet_wrap(~WATER, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("total inverts") + 
  xlab("") + 
  labs(fill = "Order")


## How many lakes are each common family found in 
common = s.d %>% 
  filter(ORDER..OR.ABOVE. != "") %>%
  filter(FAMILY != "unidentified") %>% 
  filter(FAMILY != "terrestrial") %>%
  filter(FAMILY != "") %>%
  filter(FAMILY != " ") %>%
  filter(FAMILY != "degraded") %>%
  group_by(WATER, ORDER..OR.ABOVE., FAMILY) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
  filter(ORDER..OR.ABOVE. != "") %>%
  ungroup() %>%
  group_by(ORDER..OR.ABOVE., FAMILY) %>%
  summarize(total_lakes = n())


library(stringr)
## Which families are found across all lakes
s.d %>% 
  left_join(common) %>%
  
  group_by(WATER, ORDER..OR.ABOVE., FAMILY, total_lakes) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  ungroup() %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
  filter(ORDER..OR.ABOVE. != "") %>% 
  arrange(ORDER..OR.ABOVE., ORDER..OR.ABOVE.)  %>%
  filter(total_lakes != "NA")%>% 
  
  ggplot(aes(y = FAMILY, x = sum, fill = WATER)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~total_lakes, scales = "free")

## Checking for duplicates
s.d %>% 
  left_join(common) %>% select(ORDER..OR.ABOVE., FAMILY) %>% unique() %>%
  mutate(FAMILY = str_replace(FAMILY, " ", "")) %>%
  mutate(ORDER..OR.ABOVE. = str_replace(ORDER..OR.ABOVE., " ", ""))


graph_list = list()
for(i in 1: length(unique(s.d$WATER))){
  graph = s.d %>% 
    filter(WATER ==unique(s.d$WATER)[i]) %>%
    filter(FAMILY != "unidentified") %>%
    group_by(WATER, ORDER..OR.ABOVE., FAMILY) %>% 
    summarize(sum = sum(TOTAL_N)) %>% 
    ungroup() %>% 
    mutate(FAMILY = tolower(FAMILY)) %>%
    mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>% 
    filter(ORDER..OR.ABOVE. != "") %>%
    arrange(FAMILY, sum) %>%
    mutate(FAMILY = factor(FAMILY, levels = FAMILY)) %>% 
    ggplot(aes(x = factor(FAMILY, levels = FAMILY[order(-sum)]),
               y = sum, fill = ORDER..OR.ABOVE.)) +
    theme_minimal() + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle= 90)) +
    xlab("Family") + 
    ylab("Total Inverts")
  graph_list[[i]] = graph
}
do.call("grid.arrange", c(graph_list, ncol=4))


s.d %>% separate(YSAMP, into = c("GEAR", "WATER", "YEAR", "YSAMP_N")) %>% 
  group_by(WATER, YSAMP_N, ORDER..OR.ABOVE., FAMILY, GENUS) %>% 
  summarize(sum = sum(TOTAL_N)) %>% 
  ungroup() %>% 
  mutate(FAMILY = tolower(FAMILY)) %>%
  mutate(ORDER..OR.ABOVE. = tolower(ORDER..OR.ABOVE.)) %>%
  filter(ORDER..OR.ABOVE. != "") %>%
  ggplot(aes(x = WATER, y = sum, fill =FAMILY)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ORDER..OR.ABOVE., scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90))


