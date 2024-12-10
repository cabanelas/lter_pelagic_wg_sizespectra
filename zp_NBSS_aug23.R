################################################################################
#############          Pelagic Synthesis           #############################
#############             AUG-2024                 #############################
#############           Size Spectra               #############################
## by: Alexandra Cabanelas 
################################################################################
## Calculating Size Spectra for LTER Pelagic Synthesis WG 
## WW/DW, C, ESD, NBSS
# using NES LTER zp data

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
#install.packages("here")
#install.packages("tidyverse")
library(here)
library(tidyverse)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

# Zooplankton Abundance Data with 10m2 values (100m3 is available if needed)
zp <- read.csv(file.path("output","zp_abundance_m2.csv"),
               header = T) %>% #created in zp_abundance_m2_calc_AUG2024
      select(CRUISE_NAME, STATION, CAST, GEAR_VOLUME_FILTERED, SAMPLE_SPLIT_FACTOR,
         EVENT_DATE, DAY, MONTH, YEAR, DEPTH, NET_MAX_DEPTH, TAXA_004, TAXA_NAME,
         abundance_10m2, ADULT_10m2, C5_10m2, C4_10m2, C3_10m2, C2_10m2, C1_10m2,
         crytopia_10m2, furcilia_10m2, calyptosis_10m2, nauplius_10m2, unknown_10m2)

# Published mean Lengths of species by stage
lengths <- read.csv(file.path("raw","meanlengths.csv"),
                    header = T)

# Published Length-Weight regressions of species by stage
LW_reg <- read.csv(file.path("raw","lengthweightregressions.csv"),
                   header = T)

# Published mean Weights of species by stage
pub_weights <- read.csv(file.path("raw","meanweights.csv"),
                        header = T)

# Published Dry Weight to Carbon conversions (%)
DW_to_C_convert <- read.csv(file.path("raw","DWtoC_conversions.csv"),
                            header = T)

## ------------------------------------------ ##
#            STEP 1 -----
#       ADD LENGTHS
## ------------------------------------------ ##
# Add length values to the zooplankton data
# Lengths are in micrometers (um)

# make sure the stage naming matches in both dataframes
stage_mapping <- c(
  "ADULT_10m2" = "Adult",
  "C5_10m2" = "CV",
  "C4_10m2" = "CIV",
  "C3_10m2" = "CIII",
  "C2_10m2" = "CII",
  "C1_10m2" = "CI",
  "crytopia_10m2" = "Cryptopia",
  "furcilia_10m2" = "Furcilia",
  "calyptosis_10m2" = "Calyptopis",
  "nauplius_10m2" = "Nauplius",
  "unknown_10m2"= "Not_Staged"
)

# need to change zp to long format
zp_long <- zp %>%
  pivot_longer(cols = "ADULT_10m2":"unknown_10m2",
               names_to = "stage_code", 
               values_to = "abundance_staged_10m2") %>%
  mutate(stage = stage_mapping[stage_code])

# merge 
zp_with_lengths <- merge(zp_long, lengths, 
                         by = c("TAXA_NAME", "stage"), 
                         all.x = TRUE)

zp_lengths <- zp_with_lengths %>% 
  drop_na(MEAN_LENGTH_UM) %>% #removing row (species) that dont have mean lengths
  select(-source) #removing the source column just to tidy things up


## ------------------------------------------ ##
#            STEP 2 -----
#       GET WEIGHT USING L-W Regressions 
## ------------------------------------------ ##

# Merge the zp and LW regression data frames
zp_eqns <- zp_lengths %>%
  left_join(LW_reg, by = c("TAXA_NAME", "stage")) %>% 
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%  # convert empty strings to NA
  drop_na(eqn) %>% # remove rows (species) with no L-W eqns 
  select(-c(source, X))

#apply_weight <- function(length_value, b, int) { 
#  weight <- 10^(b * log10(length_value) + int) 
#  return(weight) 
#}

apply_weight <- function(length_value, b, int, Weight_unit) { #log reg
  # convert length from µm to mm if the weight l-w is mg
  length_value <- ifelse(Weight_unit == "mg", 
                         length_value / 1000, length_value)
  
  #function to calculate weight
  weight <- 10^(b * log10(length_value) + int) 
  
  weight <- ifelse(Weight_unit == "mg", weight * 1000, weight)
  
  return(weight) 
}

# calculate weights
zp_eqns <- zp_eqns %>%
  mutate(weight = apply_weight(MEAN_LENGTH_UM, b, int, Weight_unit))

# I am also adding to this df the published mean weights to fix ones that dont have DW
zp_merged <- merge(zp_eqns, pub_weights, 
                   by = c("TAXA_NAME", "stage", "stage_code"), 
                   all.x = TRUE)

#clean it up a bit, remove unnecessary columns 
zp_merged <- zp_merged %>%
  select(TAXA_NAME, stage, CRUISE_NAME, STATION, CAST, EVENT_DATE, DAY,
         MONTH, YEAR, DEPTH, NET_MAX_DEPTH, abundance_staged_10m2,
         MEAN_LENGTH_UM, MEAN_LENGTH_TYPE, eqn, b, int, LengthType, Weight_unit,
         Weight_type, weight, PUBLISHED_MEAN_WEIGHT_UM, PUBLISHED_weight_type, 
         PUBLISHED_weight_unit)

## ------------------------------------------ ##
#            STEP 3 -----
#        Convert WW to DW
## ------------------------------------------ ##

# tried to get as many DW eqns as possible
# but for some, I was only able to find WW (or C) regressions

# I dont have great conversion factors - just for C.fin, so used
# 0.25 as the conversion from WW to DW (or published weights)

# some of my regressions give you Carbon (instead of DW/WW)
# using the published DW for these...
zp_merged$PUBLISHED_MEAN_WEIGHT_UM <- as.numeric(zp_merged$PUBLISHED_MEAN_WEIGHT_UM)

zp_merged <- zp_merged %>%
  mutate(convertedDW = case_when(
    TAXA_NAME == "Thysanoessa longicaudata" ~ PUBLISHED_MEAN_WEIGHT_UM, #eqns for euphausiids giving suspect values (very high)
    TAXA_NAME == "Thysanoessa raschii" ~ PUBLISHED_MEAN_WEIGHT_UM,
    TAXA_NAME == "Euphausia krohnii" ~ PUBLISHED_MEAN_WEIGHT_UM,
    Weight_type == "WW" ~ weight * 0.25, # convert wet weight to dry weight
    Weight_type == "C" ~ PUBLISHED_MEAN_WEIGHT_UM, # use published mean weight
    TRUE ~ weight # keep original weight for other cases
  ))

# Now all weights in the "convertedDW" column should be in ug

## ------------------------------------------ ##
#            STEP 4 -----
#        Convert DW to C
## ------------------------------------------ ##

#merge - add the DW to C values to the zp df 
zp_C <- merge(zp_merged, DW_to_C_convert, 
              by = c("TAXA_NAME", "stage"), 
              all.x = TRUE)

# multiply DW by the conversion (usually 48%)
zp_C <- zp_C %>%
  mutate(convertedCARBON = convertedDW * (Per_C_DW / 100))


## ------------------------------------------ ##
#            STEP 5 -----
#        Calculate ESD 
## ------------------------------------------ ##
# assume that the organism can be approximated as a sphere
# by using the formula for the volume of a sphere to find its diameter

#         Step 5.1: Calculate Volume
# Volume = 4/3 pi r^3
# r = prosome length / 2
#zp_ESD <- zp_C %>%
#  mutate(r = MEAN_LENGTH_UM/2,
#         volume = (4/3)*pi*r^3)

#         Step 5.2: Calculate Density - not sure if needed
# Density = DW / V 
#zp_ESD <- zp_ESD %>%
#  mutate(density = convertedDW/volume)

#         Step 5.3: Calculate ESD
# ESD = (3/4pi) * (DW/density)^1/3
#This formula derives from the volume of a sphere formula and assumes you are
#calculating the diameter that would give you the same volume given the density and dry weight.
#zp_ESD <- zp_ESD %>%
#  mutate(ESD = ((4*pi/3) * (convertedDW/density))^(1/3))
# the density might be another source of error here..?

# OR
#zp_ESD <- zp_ESD %>%
#  mutate(ESDv2 = (6 * volume / pi)^(1 / 3))  

# OR ******************
#i think this the way
zp_ESD <- zp_C %>%
  mutate(ESD = (((4 * pi / 3) * ((convertedDW / 1e6) / 1.05))^(1/3)) * 1e4)
#water density is in g/cm^3

# OR 
#zp_biomass_clean <- zp_ESD %>%
#  mutate(ESDv4 = (3 * biomass_mg / (4 * pi * 1.05))^(1/3))

## ------------------------------------------ ##
#            STEP 6 -----
#        Calculate Biomass 
## ------------------------------------------ ##

# Biomass = abundance * DW 
zp_biomass <- zp_ESD %>%
  mutate(biomass_ug = abundance_staged_10m2 * convertedCARBON,#ug/10m2
         biomass_mg = (biomass_ug / 1000) / 10)#convert to mg/10m2 by diving by 1,000 
# DO I NEED TO DIVIDE BY 10 ?

## ------------------------------------------ ##
#            STEP 7 -----
#        BINNING 
## ------------------------------------------ ##

#### NEED TO DO THE GROUPING BY YEAR OR SEASON..!! #######

zp_biomass_clean <- zp_biomass %>%
  filter(!is.na(biomass_ug)) %>%
  filter(biomass_ug > 0) #%>%
  #filter(TAXA_NAME != "Euphausia krohnii")

#         Step 7.1: Binning 
#  Assign each ESD to bins - octave scaling = factor of 2 
#esd_bins <- 2^(0:15) #1, 2, 4, 8, 16... 
# logarithmically equal size bins

zp_biomass_clean <- zp_biomass_clean %>%
  group_by(YEAR) %>%
  mutate(esd_bins = cut(ESD, 
                        breaks = 2^(seq(floor(log2(min(ESD, 
                                                       na.rm = TRUE))), 
                                        ceiling(log2(max(ESD, na.rm = TRUE))), 
                                        by = 1)), 
                        include.lowest = TRUE))

#         Step 7.2: calculate geometric mean of bin edges 
#square root of the product of the bin edges (or for log-transformed values, the midpoint in log space)
bin_biomass <- zp_biomass_clean %>%
  group_by(YEAR, esd_bins) %>%
  summarize(total_biomass = sum(biomass_mg, na.rm = TRUE)) %>%
  mutate(ESD_bin_edges = strsplit(as.character(esd_bins), ","), # splits the bin edges
         ESD_bin_edges = lapply(ESD_bin_edges, 
                                function(x) as.numeric(gsub("[^0-9.]", "", x))),
         ESD_bin_edges = lapply(ESD_bin_edges, sort),
         ESD_bin_midpoint = sapply(ESD_bin_edges, 
                                   function(x) sqrt(x[1] * x[2]))) #calculates geometric mean for each bin (midpoint)
#does ESD_bin need to be log10 ?? for plot ? yes...

## ------------------------------------------ ##
#            STEP 8 -----
#        Normalize Biomass 
## ------------------------------------------ ##

# calculate Normalized Biomass
# biomass / bin width 

bin_biomass <- bin_biomass %>%
  #DO I NEED TO GROUP BY BIN HERE ???
  mutate(bin_width = sapply(ESD_bin_edges, 
                            function(x) x[2] - x[1]), # calculates width of each bin
         normalized_biomass = total_biomass / bin_width, 
         log_normalized_biomass = log10(normalized_biomass))



## ------------------------------------------ ##
#            STEP 9 -----
#        Plot 
## ------------------------------------------ ##

ggplot(bin_biomass, aes(x = log10(ESD_bin_midpoint), 
                        y = log_normalized_biomass,
                        color = as.factor(YEAR))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Log(Geometric Mean of Size Bin Edges)", 
       y = "Log(Normalized Biomass)") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 5), 
                     breaks = 0:5, labels = scales::math_format(10^.x))


## ------------------------------------------ ##
#            STEP 10 -----
#        Calculate slope of the size spectra 
## ------------------------------------------ ##

# linear regression to get slope of the size spectrum.

model <- lm(log_normalized_biomass ~ log10(ESD_bin_midpoint), data = bin_biomass)
summary(model)


## ------------------------------------------ ##
#            STEP 11 -----
#        Other plots 
## ------------------------------------------ ##
ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, 
                             y = convertedDW, 
                             color = TAXA_NAME,
                             shape = stage)) +
  geom_point() + 
  facet_wrap(~TAXA_NAME, scales = "free")

ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, 
                             y = convertedDW, 
                             color = TAXA_NAME,
                             shape = stage)) +
  geom_point() 

zp_biomass_clean <- zp_biomass_clean %>%
  mutate(Group = case_when(
    TAXA_NAME %in% c("Calanus finmarchicus", "Centropages hamatus", 
                     "Centropages typicus", "Clausocalanus arcuicornis", 
                     "Metridia lucens", "Nannocalanus minor", 
                     "Paracalanus parvus", "Pseudocalanus minutus", 
                     "Temora longicornis") ~ "Copepod",
    TAXA_NAME %in% c("Thysanoessa longicaudata", "Thysanoessa raschii", "Euphausia krohnii") ~ "Euphausiid",
    TRUE ~ "Other"  
  ))

ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, 
                             y = convertedDW, 
                             color = TAXA_NAME,
                             shape = stage)) +
  geom_point() + facet_wrap(~Group, scales = "free") +
  scale_x_log10() +
  scale_y_log10() 

ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, y = biomass_ug)) +
  geom_point()

ggplot(zp_biomass_clean, aes(x = ESDv3, y = biomass_mg)) +
  geom_point() 

ggplot(zp_biomass_clean, aes(x = convertedDW, y = convertedCARBON)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 


ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, y = convertedCARBON)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 

ggplot(zp_biomass_clean, aes(x = ESDv3)) +  
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Length-Frequency Distribution",
       x = "Equivalent Spherical Diameter (ESD)",
       y = "Frequency") +
  theme_minimal()

site1 <- zp_biomass_clean %>%
  dplyr::filter(CRUISE_NAME == "EN687") %>%
  dplyr::filter(STATION == "7")

ggplot(site1, aes(x= ESDv3, y = abundance_staged_10m2, color = TAXA_NAME)) +
  geom_point(size = 2)+
  scale_x_log10() +
  scale_y_log10()

ggplot(site1, aes(x= MEAN_LENGTH_UM, y = abundance_staged_10m2, color = TAXA_NAME, shape = stage)) +
  geom_point(size = 4)+
  scale_x_log10() +
  scale_y_log10()