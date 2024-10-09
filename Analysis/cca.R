load("Data/RData/chem_data.RData")
chemistry = chemistry %>%
  rename(depth.5mgL = min_depth,
         temp.5mgL = min_temp)

## species richness trends

richness = full %>%
  select(WATER, MONTH, FAMILY, GENUS) %>%
  group_by(WATER, MONTH) %>%
  unique() %>%
  summarize(richness = n()) %>%
  rename("season" = "MONTH") %>%
  mutate(season = tolower(season)) %>%
  left_join(chemistry)


richness %>% 
  ggplot(aes(x = season, y = richness)) + 
  geom_boxplot() + 
  facet_wrap(~WATER)


## CCA analysis
cca.dat = full %>%
  #select(WATER, MONTH, FAMILY, GENUS) %>% # Genus 
  select(WATER, MONTH, FAMILY) %>% # Family
  group_by(WATER, MONTH) %>%
  unique() %>% 
  #unite("ID", FAMILY, GENUS) %>% # Genus
  rename(ID = FAMILY) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = ID, values_from = presence) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  rename(season = MONTH) %>%
  mutate(season = tolower(season))

cca.dat.com = cca.dat %>% ungroup() %>% unite("ID", WATER, season) %>%
  column_to_rownames(var = "ID") 


cca.dat.env = cca.dat %>% select(WATER, season) %>% 
  left_join(chemistry) %>%
  ungroup() %>%
  select(max_depth, depth.5mgL, temp.5mgL, surface_area, FieldPh, DOC.1, 
         ANC, SIO2, FUI.num, sechi.depth)

cca.dat.env 

cor.env = cor(cca.dat.env)
corrplot::corrplot(cor.env)

cca_result <- cca(cca.dat.com ~ ., data = cca.dat.env %>% select(
  surface_area, temp.5mgL,depth.5mgL,  FieldPh, DOC.1,SIO2, sechi.depth, FUI.num 
))

# Summary of the CCA results
summary(cca_result)
print(cca_result)
# Plot the results
plot(cca_result)
variables = cca_result$CCA$biplot  %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID")
waters = cca_result$CCA$wa %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID")

ggplot() +
  theme_minimal() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(data = waters, aes(label = ID,                                                 x = CCA1, y = CCA2), size = 3) +
  geom_text_repel(data = variables, 
                  aes(label = ID, 
                      x = CCA1, 
                      y = CCA2), 
                  col = "brown") + 
  geom_segment(data = variables, aes(x = 0, y = 0, 
                                     xend = CCA1, 
                                     yend = CCA2),
               col = "brown", 
               alpha = 0.5, 
               arrow = arrow(length = unit(0.1, "inches"))) + 
  theme_minimal(base_size = 15) +
  theme(legend.position = "none") + 
  xlim(-1.2, 1.5)


