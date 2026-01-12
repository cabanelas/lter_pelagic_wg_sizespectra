################################################################################
#############          Pelagic Synthesis WG        #############################
#############             NES - LTER               #############################
#############     Zooplankton Size Spectra         #############################
## by: Alexandra Cabanelas 
## created: AUG-2024, edited AUG-2025, most recent edit DEC-2025
################################################################################

# Calculating Zooplankton Size Spectra for LTER Pelagic Synthesis WG 
# NES-LTER zooplankton abundance data: 
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nes.25.2
# abundances are from 335 um mesh bongo net

# WW/DW, C, ESD, NBSS

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

# Published mean Lengths of taxa 
lengths <- read.csv(file.path("raw","meanlengths.csv"),
                    header = T)

# Published Length-Weight regressions of taxa 
LW_reg <- read.csv(file.path("raw","lengthweightregressions.csv"),
                   header = T)

# Published mean Weights of taxa 
pub_weights <- read.csv(file.path("raw","meanweights.csv"),
                        header = T)

# Published Dry Weight to Carbon conversions (%)
DW_to_C_convert <- read.csv(file.path("raw","DWtoC_conversions.csv"),
                            header = T)

# DO I NEED TO ALSO ADD THE PUBLISHED CARBON VALS HERE? 

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
  "cryptopia_10m2" = "Cryptopia",
  "furcilia_10m2" = "Furcilia",
  "calyptopis_10m2" = "Calyptopis", 
  "nauplius_10m2" = "Nauplius",
  "unknown_10m2"= "Not_Staged"
)

# change staged zp df to long format
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
  mutate(stage = "unstaged") %>%  #add stage column with "unstaged"
  rename(abundance_10m2 = conc_10m2)

# combine zp and zp_long_staged
zp_combined <- bind_rows(zp, zp_long_staged)

# rename 
#lengths <- lengths %>%
#  rename(taxa_name = TAXA_NAME) 

# merge
zp_with_lengths <- merge(zp_combined, lengths, 
                         by = c("taxa_name", "stage"), 
                         all.x = TRUE) %>%
  select(-c(taxa_code, source))

#zp_lengths <- zp_with_lengths %>% 
#  drop_na(MEAN_LENGTH_UM) %>% #removing row (species) that dont have mean lengths
#  select(-c(taxa_code, source)) #removing columns just to tidy things up

## ------------------------------------------ ##
#            STEP 2 -----
#       GET WEIGHT USING L-W Regressions 
## ------------------------------------------ ##

# Merge the zp and LW regression data frames
#LW_reg <- LW_reg %>%
#  rename(taxa_name = TAXA_NAME) 

# Merge the zp and LW regression data frames
zp_eqns <- zp_lengths %>%
  left_join(LW_reg, by = c("taxa_name", "stage")) %>% 
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%  # convert empty strings to NA
  #drop_na(eqn) %>% # remove rows (species) with no L-W eqns 
  select(-c(source, X))

apply_weight <- function(length_value, b, int, Weight_unit) { #log reg
  # convert length from µm to mm if the weight l-w is mg
  length_value <- ifelse(weight_unit == "mg", 
                         length_value / 1000, length_value)
  
  # function to calculate weight
  weight <- 10^(b * log10(length_value) + int) 
  
  weight <- ifelse(weight_unit == "mg", weight * 1000, weight)
  
  return(weight) 
}

# calculate weights
zp_eqns <- zp_eqns %>%
  mutate(weight = apply_weight(mean_length_um, b, int, weight_unit))

# adding the published mean weights to fix ones that dont have DW
#rename
#pub_weights <- pub_weights %>%
#  rename(taxa_name = TAXA_NAME) 

# merge
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
zp_merged$published_mean_weight_um <- as.numeric(zp_merged$published_mean_weight_um)

zp_merged <- zp_merged %>%
  mutate(convertedDW_ug = case_when(
    taxa_name == "Thysanoessa longicaudata" ~ published_mean_weight_um, #eqns for euphausiids giving suspect values (very high)
    taxa_name == "Thysanoessa raschii" ~ published_mean_weight_um,
    taxa_name == "Euphausia krohnii" ~ published_mean_weight_um,
    weight_type == "WW" ~ weight * 0.25, # convert wet weight to dry weight
    weight_type == "C" ~ published_mean_weight_um, # use published mean weight
    TRUE ~ weight # keep original weight for other cases
  ))

## NEED TO CHECK THIS - NEED TO KEEP PUBLISHED WEIGHTS FOR ONES I DONT HAVE ANYTHING ELSE

# Now all weights in the "convertedDW" column should be in ug

## ------------------------------------------ ##
#            STEP 4 -----
#        Convert DW to C
## ------------------------------------------ ##

#rename
#DW_to_C_convert <- DW_to_C_convert %>%
#  rename(taxa_name = TAXA_NAME) 

#merge - add the DW to C values to the zp df 
zp_C <- merge(zp_merged, DW_to_C_convert, 
              by = c("taxa_name", "stage", "stage_code"), 
              all.x = TRUE)

# multiply DW by the conversion (usually 48%)
zp_C <- zp_C %>%
  mutate(convertedCARBON = convertedDW_ug * (Per_C_DW / 100)) %>% 
  select(-Source)

## NEED TO ADD PUBLISHED C VALS 

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
# µg → g is × 10^−6
zp_ESD <- zp_C %>%
  mutate(convertedDW_g = convertedDW_ug * 10^-6)

# convert mass to volume 
zp_ESD <- zp_ESD %>%
  mutate(ESD_cm = ((4*pi/3) * (convertedDW_g/1.05))^(1/3))  #i think this the way
#add 10^-6 and 10^-4 outside 
# convertedDW = ug
# 1.05 = density of zooplankton (estimate) = g/cm^3
# 1 cm = 10^4 um
# dividing DW by 1.05 = normalizes the dry weight into a volume dimension
zp_ESD <- zp_C %>%
  mutate(convertedDW_g = convertedDW_ug * 1e-6,
         volume_cm3 = convertedDW_g / 1.05,
         ESD_cm = (6 * volume_cm3 / pi)^(1/3)) ## need to check this and line 230
# or ESD_cm = 2 * ((3 * V_cm3) / (4 * pi))^(1/3), ??

# Convert ESD from cm to µm 
zp_ESD <- zp_ESD %>% 
  mutate(ESD_um = ESD_cm * 10^4)

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
  mutate(biomass_ug = abundance_10m2 * convertedDW_ug,#ug/10m2
         biomass_mg = biomass_ug / 1000)#convert to mg/10m2 by diving by 1,000 

# or its abundance * carbon? 
zp_biomass <- zp_ESD %>%
  mutate(B_ugC_10m2 = abundance_10m2 * C_ug_ind,
         B_mgC_m2 = (B_ugC_10m2 / 10) * 1e-3)

## ------------------------------------------ ##
#            STEP 7 -----
#         Octave Binning  
## ------------------------------------------ ##
zp_binned <- zp_biomass %>%
  mutate(log2C = log2(C_ug_ind),
         bin = floor(log2(ESD_um)))   # or cut() with explicit breaks

nbss <- zp_binned %>%
  group_by(sample_id, size_bin) %>%  # or just size_bin
  summarise(B_mgC_m2_bin = sum(B_mgC_m2, na.rm = TRUE),
            .groups = "drop")

nbss <- nbss %>%
  mutate(log10_size = size_bin * log10(2),    # if using size_bin as log2(mass)
         log10_B = log10(B_mgC_m2_bin))

nbss_lm <- lm(log10_B ~ log10_size, data = nbss)
summary(nbss_lm)

########

breaks <- 2^(seq(log2(10), log2(10000), by = 1))  # from 10 to 10,000 µm
zp_binned <- zp_ESD %>%
  mutate(ESD_bin = cut(ESD_um, breaks = breaks, include.lowest = TRUE))

nbss <- zp_binned %>%
  group_by(ESD_bin) %>%
  summarise(B_mgC_m2_bin = sum(B_mgC_m2, na.rm = TRUE),
            .groups = "drop")

nbss <- nbss %>%
  mutate(bin_min = as.numeric(gsub("\\[|\\(|,.*", "", ESD_bin)),
         bin_max = as.numeric(gsub(".*,|\\]", "", ESD_bin)),
         bin_width = bin_max - bin_min,
         B_norm = B_mgC_m2_bin / bin_width)

library(ggplot2)

ggplot(nbss, aes(x = bin_min, y = B_norm)) +
  geom_line(color = "blue", size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Equivalent Spherical Diameter (µm)",
       y = "NBSS (mg C m⁻² µm⁻¹)",
       title = "Normalized Biomass Size Spectrum") +
  theme_minimal()

nbss <- nbss %>%
  mutate(log10_ESD = log10(bin_min),
         log10_B = log10(B_norm))

nbss_lm <- lm(log10_B ~ log10_ESD, data = nbss)
summary(nbss_lm)