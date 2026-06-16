## Oxycline modelling
set.seed(1234)
library(dplyr)
library(tidyr)
library(mgcv)
library(purrr)

## Temp and DO data frame
tdo = read.csv("Data/CSVs/Water_Chemistry/Temp_DO.csv") %>%
  separate(Date, into = c("month", "date", "year")) %>%
  mutate(season = case_when(month < 8 ~ "spring", month>7 ~ "fall")) %>%
  unite("ID", c(WATER, season), sep = ".")



## Setup the subsetted TDO frame for just fall dissolved oxygen vs. depth
depth.data.DO = tdo %>% 
  select(ID,  depth, DO) %>%
  mutate(ID = as.factor(ID))

## Oxycline function -----------------------------
#### Generates the oxycline using a factor-smoothed model, where lake is an interaction with depth
#### Allows lakes to have their own oxycline shapes
extract_oxycline = function(data = depth.data.DO,
                                 k = 10,
                                 n_grid = 500,
                                 threshold_prop = 0.5) {
data = depth.data.DO
  # Fit GAM (lake-specific smooths)
  mod = gam(DO ~ s(depth, ID, k = k, bs = "fs"),
    data = data, method = "REML")

  # Create prediction grid per lake
  newdat = data %>%
    group_by(ID) %>%
    summarise(
      depth = list(seq(min(depth, na.rm = TRUE),
                       max(depth, na.rm = TRUE),
                       length.out = n_grid)),
      .groups = "drop"
    ) %>%
    unnest(depth)

  # Predict DO
  newdat$oxy_pred = predict(mod, newdata = newdat)

  # Compute gradient within each lake
  newdat = newdat %>%
    group_by(ID) %>%
    arrange(depth, .by_group = TRUE) %>%
    mutate(
      doxy = c(NA, diff(oxy_pred)), ## Change is DO measure
      ddepth = c(NA, diff(depth)), ## change in depth
      grad = doxy / ddepth ## slope is change in DO / change in depth
    ) %>%
    ungroup()

  newdat %>% 
    group_by(ID) %>%
    summarize(max_depth = max(depth))
  # 5. Extract oxycline metrics per lake
  oxycline_metrics = newdat %>% 
    group_by(ID) %>% 
    summarise(
      DDO5 = approx(oxy_pred, depth, xout = 5)$y ## estimate the depth where dissolved oxygen is < 5mg/L

    ) %>%
    ## If the oxygen never dips below 5mg/L (or it does so at bottom) the DDO5 should be the max depth of the lake/profile
    left_join( newdat %>% 
    group_by(ID) %>%
    summarize(max_depth = max(depth))) %>%
    mutate(DDO5 = case_when(is.na(DDO5) ~ max_depth, !is.na(DDO5) ~ DDO5)) %>%
    ## Now remove the max depth
    select(-max_depth)
    
  
  
  oxycline = newdat %>%
    group_by(ID) %>%
    summarise(
      mid_depth = depth[which.max(abs(grad))], ## middle of thermocline is depth where max slope is achieved

      grad_max = max(abs(grad), na.rm = TRUE), ## Maximum slope observed along curve
      thresh = threshold_prop * grad_max, ## create a threshold that is defined by the max slope * threshold prop

      top = min(depth[abs(grad) >= thresh], na.rm = TRUE), ## top is the first depth observation where the threshold is exceeded
      bottom = max(depth[abs(grad) >= thresh], na.rm = TRUE), ## bottom is the bottom observation where the threshold is exceeded

      .groups = "drop"
    )

  ## List of items to be output
  list(
    model = mod,
    oxycline = oxycline,
    metrics = oxycline_metrics,
    profile = newdat
  )

}

## Run the function on the DO data
out = extract_oxycline(depth.data.DO) ## Apply function to DO/Depth Data

## Visualize curves with observed vs. predicted
ggplot(data = out$profile %>%
         as.data.frame() %>% 
         arrange(depth, oxy_pred), aes(x = oxy_pred, y = depth )) + 
  geom_point(size = .5) +
  geom_point(data = tdo, mapping = aes(x = DO, y = depth), col = "red") +
  scale_y_reverse() + 
  geom_hline(data = out$oxycline, aes(yintercept = mid_depth)) + 
  geom_hline(data = out$oxycline, aes(yintercept = top)) +
  geom_hline(data = out$oxycline, aes(yintercept = bottom)) +
  facet_wrap(~ID, scales = "free") + theme_bw() +
  geom_vline(aes(xintercept = 5), col = "blue") +
  ylab("Depth (m)") + 
  xlab("Dissolved Oxygen (mg/L)")

## Model checking
out$model %>% summary()
out$model %>% gam.check()

## Save metrics
oxycline.met = out$metrics


### TDO5 ---------

## Temperature where DO < 5mgL
## Load in the thermocline data
load(file = "Data/RData/thermocline_predicted_profiles.RData")
oxycline.metrics = thermocline.profiles %>% ## Start with temp profiles
  left_join(out$metrics ) %>% ## Get TDO5
  group_by(ID) %>%
  slice_min(abs(depth - DDO5), n = 1) %>% ## Selects for depth closest to DDO5
  select(ID, temp_pred, DDO5) %>%
  rename(TDO5 = temp_pred)

## Write CSV -------
write.csv(oxycline.metrics, file = "Data/CSVs/Water_Chemistry/Oxycline.csv", row.names = F)














 
 