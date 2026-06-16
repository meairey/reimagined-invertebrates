## Library setup
set.seed(1234)
library(tidyverse)
library(mgcv)
library(purrr)
library(ggplot2)

## Temp and DO data frame
tdo = read.csv("Data/CSVs/Water_Chemistry/Temp_DO.csv") %>%
  separate(Date, into = c("month", "date", "year")) %>%
  mutate(season = case_when(month < 8 ~ "spring", month > 7 ~ "fall")) %>%
  unite("ID", c(WATER, season), sep = ".")

## Subsetted temperature data
depth.data.temp = tdo %>% 
  #filter(season == "fall")  %>%
  select(ID, depth, temp)  %>%
  mutate(ID = as.factor(ID))


## Function for extracting thermocline depths (bottom, mid, top)
extract_thermocline <- function(data = data ,
                                 k = 10,
                                 n_grid = 500,
                                 threshold_prop = 0.5) {

  # 1. Fit GAM (lake-specific smooths)
  mod <- gam(
    temp ~ s(depth, ID, bs = "fs", k = k),
    data = data,
    method = "REML"
  )

  # 2. Create prediction grid per lake
  newdat <- data %>%
    group_by(ID) %>%
    summarise(
      depth = list(seq(min(depth, na.rm = TRUE),
                       max(depth, na.rm = TRUE),
                       length.out = n_grid)),
      .groups = "drop"
    ) %>%
    unnest(depth)

  # 3. Predict temperature
  newdat$temp_pred <- predict(mod, newdata = newdat)

  # 4. Compute gradient within each lake
  newdat <- newdat %>%
    group_by(ID) %>%
    arrange(depth, .by_group = TRUE) %>%
    mutate(
      dtemp = c(NA, diff(temp_pred)),
      ddepth = c(NA, diff(depth)),
      grad = dtemp / ddepth
    ) %>%
    ungroup()

  # 5. Extract thermocline metrics per lake
  thermo <- newdat %>%
    group_by(ID) %>%
    summarise(
      mid_depth = depth[which.max(abs(grad))],

      grad_max = max(abs(grad), na.rm = TRUE),
      thresh = threshold_prop * grad_max,

      top = min(depth[abs(grad) >= thresh], na.rm = TRUE),
      bottom = max(depth[abs(grad) >= thresh], na.rm = TRUE),

      .groups = "drop"
    )

  list(
    model = mod,
    thermocline = thermo,
    profile = newdat
  )

} ## end function

## Run function on the depth data specified above

out.thermocline <- extract_thermocline(data = depth.data.temp)

## Plot the predicted vs. observed
ggplot(data = out.thermocline$profile %>%
         as.data.frame() %>% 
         arrange(depth, temp_pred), aes(x = temp_pred, y = depth )) + 
  geom_point(size = .5) +
  geom_point(data = tdo , mapping = aes(x = temp, y = depth, col = ID), col = "red") +
  scale_y_reverse() + 
  geom_hline(data = out.thermocline$thermocline, aes(yintercept = mid_depth)) + 
  geom_hline(data = out.thermocline$thermocline, aes(yintercept = top)) +
  geom_hline(data = out.thermocline$thermocline, aes(yintercept = bottom)) +
  facet_wrap(~ID, scales = "free") + theme_bw() + 
  ylab("Depth (m)") + 
  xlab("Temperature (C)")

## Output summary information for other scripts
thermocline = out.thermocline$thermocline
thermocline.profiles = out.thermocline$profile
#save(thermocline.profiles, file = "Data/RData/thermocline_predicted_profiles.RData")
#write.csv(thermocline, file = "Data/CSVs/Water_Chemistry/Thermocline.csv", row.names = F)





 
 