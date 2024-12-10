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
zp <- read.csv(file.path("raw",
                         "nes-lter-zp-abundance-335um-unstaged10m2.csv"),
               header = T)

zp_staged <- read.csv(file.path("raw",
                                "nes-lter-zp-abundance-335um-staged10m2.csv"),
                      header = T)


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
  "adult_10m2" = "Adult",
  "c5_10m2" = "CV",
  "c4_10m2" = "CIV",
  "c3_10m2" = "CIII",
  "c2_10m2" = "CII",
  "c1_10m2" = "CI",
  "crytopia_10m2" = "Cryptopia",
  "furcilia_10m2" = "Furcilia",
  "calytopisis_10m2" = "Calyptopis", #need to fix later when i fix spelling not important
  "nauplius_10m2" = "Nauplius",
  "unknown_10m2"= "Not_Staged"
)

# need to change zp to long format
zp_long_staged <- zp_staged %>%
  pivot_longer(cols = "adult_10m2":"unknown_10m2",
               names_to = "stage_code", 
               values_to = "abundance_10m2") %>%
  mutate(stage = stage_mapping[stage_code]) %>%
  select(-c(conc_10m2, adult_count, c5_count, c4_count, c3_count, 
            c2_count, c1_count, cryptopia_count, furcilia_count, calyptopis_count,
            nauplius_count, unknown_count))

# add unstaged taxa
zp <- zp %>%
  mutate(stage = "unstaged") %>%  #add the 'stage' column with value "unstaged"
  rename(abundance_10m2 = conc_10m2)

# combine zp and zp_long_staged
zp_combined <- bind_rows(zp, zp_long_staged)

# rename 
lengths <- lengths %>%
  rename(taxa_name = TAXA_NAME) 
#merge
zp_with_lengths <- merge(zp_combined, lengths, 
                         by = c("taxa_name", "stage"), 
                         all.x = TRUE)

zp_lengths <- zp_with_lengths %>% 
  drop_na(MEAN_LENGTH_UM) %>% #removing row (species) that dont have mean lengths
  select(-c(taxa_code, source)) #removing columns just to tidy things up


## ------------------------------------------ ##
#            STEP 2 -----
#       GET WEIGHT USING L-W Regressions 
## ------------------------------------------ ##

# Merge the zp and LW regression data frames
LW_reg <- LW_reg %>%
  rename(taxa_name = TAXA_NAME) 

# Merge the zp and LW regression data frames
zp_eqns <- zp_lengths %>%
  left_join(LW_reg, by = c("taxa_name", "stage")) %>% 
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%  # convert empty strings to NA
  drop_na(eqn) %>% # remove rows (species) with no L-W eqns 
  select(-c(source, X))


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
#rename
pub_weights <- pub_weights %>%
  rename(taxa_name = TAXA_NAME) 
#merge
zp_merged <- merge(zp_eqns, pub_weights, 
                   by = c("taxa_name", "stage", "stage_code"), 
                   all.x = TRUE) %>%
  select(-source)


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
    taxa_name == "Thysanoessa longicaudata" ~ PUBLISHED_MEAN_WEIGHT_UM, #eqns for euphausiids giving suspect values (very high)
    taxa_name == "Thysanoessa raschii" ~ PUBLISHED_MEAN_WEIGHT_UM,
    taxa_name == "Euphausia krohnii" ~ PUBLISHED_MEAN_WEIGHT_UM,
    Weight_type == "WW" ~ weight * 0.25, # convert wet weight to dry weight
    Weight_type == "C" ~ PUBLISHED_MEAN_WEIGHT_UM, # use published mean weight
    TRUE ~ weight # keep original weight for other cases
  ))



# Now all weights in the "convertedDW" column should be in ug

## ------------------------------------------ ##
#            STEP 4 -----
#        Convert DW to C
## ------------------------------------------ ##

#rename
DW_to_C_convert <- DW_to_C_convert %>%
  rename(taxa_name = TAXA_NAME) 
#merge - add the DW to C values to the zp df 
zp_C <- merge(zp_merged, DW_to_C_convert, 
              by = c("taxa_name", "stage", "stage_code"), 
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

# convertedDW is in micrograms (µg) 
# Convert convertedDW from µg to grams (g) 
convertedDW_g <- zp_merged$convertedDW * 10^-6

zp_ESD <- zp_ESD %>%
  mutate(ESD_cm = ((4*pi/3) * (convertedDW_g/1.05))^(1/3))  #i think this the way
#add 10^-6 and 10^-4 outside 
# convertedDW = ug
# 1.05 = density of zooplankton (estimate) = g/cm^3
# 1 cm = 10^4 um
# dividing DW by 1.05 = normalizes the dry weight into a volume dimension

# Convert ESD from cm to µm 
zp_ESD <- zp_ESD %>% mutate(ESD_um = ESD_cm * 10^4)

## other approaches to this 
# ESD = (3/4pi) * (DW/density)^1/3
#This formula derives from the volume of a sphere formula and assumes you are
#calculating the diameter that would give you the same volume given the density and dry weight.
# OR
#zp_ESD <- zp_ESD %>%
#  mutate(ESD = ((4*pi/3) * (convertedDW/density))^(1/3))
# the density might be another source of error here..?
# OR
#zp_ESD <- zp_ESD %>%
#  mutate(ESDv2 = (6 * volume / pi)^(1 / 3))  
# OR 
#zp_biomass_clean <- zp_ESD %>%
#  mutate(ESDv4 = (3 * biomass_mg / (4 * pi * 1.05))^(1/3))

## ------------------------------------------ ##
#            STEP 6 -----
#        Calculate Biomass 
## ------------------------------------------ ##

# Biomass = abundance * DW 
zp_biomass <- zp_ESD %>%
  mutate(biomass_ug = abundance_10m2 * convertedDW,#ug/10m2
         biomass_mg = biomass / 1000)#convert to mg/10m2 by diving by 1,000 
