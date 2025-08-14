################################################################################
#############          Pelagic Synthesis           #############################
#############             AUG-2024                 #############################
#############           Size Spectra               #############################
## by: Alexandra Cabanelas 
################################################################################

# just adding net max depth column to the 
# zp_NESLTER_10m2_ESD_11DEC24 csv

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
metadata <- read.csv(file.path("raw",
                               "nes-lter-zooplankton-tow-metadata-v2.csv"),
                     header = T) #from zooplankton inventory package

zp <- read.csv(file.path("output",
                         "zp_NESLTER_10m2_ESD_11DEC24.csv"),
               header = T)


## add sample name column to zp
zp <- zp %>%
  mutate(sample_name = paste(cruise, station, paste0("B", cast), sep = "_"))


## add net depth 
zp <- zp %>%
  left_join(metadata %>% 
              select(sample_name, net_max_depth_m), 
            by = "sample_name")



zp_na <- zp %>% filter(is.na(net_max_depth_m))

#due to a typo in the metadata df that i wont fix right now i need to manually assign

## 121 m AR37 L7 B9 (sample name is wrong in the metadata file)

## 231.05 m Lu11c B29

zp <- zp %>%
  mutate(net_max_depth_m = case_when(
    sample_name == "EN627_u11c_B29" ~ 231.05,
    sample_name == "AR38_L7_B9" ~ 121,
    TRUE ~ net_max_depth_m  # Keep existing values
  ))

zp_na <- zp %>% filter(is.na(net_max_depth_m))

zp <- zp %>%
  select(-1)

# move columns
zp <- zp %>%
  relocate(sample_name, .after = cast) %>%
  relocate(net_max_depth_m, .after = day_night)



write.csv(zp, "output/zp_NESLTER_10m2_ESD_11DEC24_wNetDepth.csv")
check <- read.csv(file.path("output",
                            "zp_NESLTER_10m2_ESD_11DEC24_wNetDepth.csv"),
                  header = T)
