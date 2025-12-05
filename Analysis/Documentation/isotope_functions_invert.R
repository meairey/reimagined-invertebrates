
axis_order = c("aeshnidae","gomphidae", "libellulidae", "corduliidae", "coenagrionidae", "heptageniidae", "leptophlebiidae", "ephemerellidae", "notonectidae", "viviparidae", "lymnaeidae","phryganeidae","limnephilidae", "hydropsychidae", "chironomidae", "asellidae", "talitridae")

## Invert isotopic functions

overlap = function(data_input, comm, dr, posterior){
  ## Remove this later
  

  
  data_overlap = data_input %>% 
    filter(community == comm)
  
  #data_overlap = combo %>% filter(community == 1)
  spp=length(unique(data_overlap$group))
  print(data_overlap)
  
  name_list = names(posterior) 
  
  
  name_matrix = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","G1"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F) %>%
    filter(C1 == comm & C2 == comm) %>%
    select(Var1, Var2) %>%
    filter(Var1 != Var2)
  
  
  names_modified = expand.grid(name_list, name_list) %>% 
    separate(Var1, into = c("C1","group"), remove = F) %>%
    separate(Var2, into = c("C2", "G2"), remove = F) %>%
    filter(C1 == comm & C2 == comm) %>%
    mutate(group = as.numeric(group))%>%
    left_join(legend, #%>% 
               # mutate(community = as.character(community)),
              by = c("group")) %>% #, "C1" = "community"
    select(-color, -common) %>%
    rename("G1" = group) %>%
    rename(group = "G2") %>%
    rename(Code1 = CODE) %>%
    mutate(group = as.numeric(group)) %>%
    left_join(legend, #%>% 
               # mutate(community = as.character(community)), 
              by = c("group")) %>% #, "C1" = "community"
    select(-color, -common) %>%
    rename(G2 = group) %>%
    unite(col = Name, c(Code1, CODE), sep = " v " ) %>%
    filter(Var1 != Var2) %>%
    select(Name)
  
  
  ma = matrix(NA, nrow = dr, ncol = length(name_matrix[,1]))
  
  for(i in 1:length(name_matrix[,1])){ 
    
    overlap = bayesianOverlap(ellipse1 = as.character(name_matrix[i,1]), 
                              ellipse2 = as.character(name_matrix[i,2]), 
                              posterior, 
                              draws = dr,
                              p.interval = 0.95,
                              n = dr)
    
    bayes.prop = (overlap[,3] / (overlap[,2] + overlap[,1] - overlap[,3]))
    
    ma[,i] = bayes.prop
    
  }
  
  
  
  full_olap_mat = ma %>% as.data.frame() %>%
    rename_with(~ setNames((unique(names_modified[,1])), .), everything())
  
  full_olap_mat[dr+1,] = unique(names_modified[,1])
  
  rownames(full_olap_mat) = c((data.frame(n = "com",c = comm, post = c(1:dr)) %>% unite("rowname", c(n, c, post)))$rowname, "SppPair")
  
  
  #return(as.data.frame(t(olap_mat))) ## This is a reminder to pull 90% CI intervals out of this at a point when i can play with code and don't need to be writing
  return(t(as.data.frame(full_olap_mat)))
  
}
