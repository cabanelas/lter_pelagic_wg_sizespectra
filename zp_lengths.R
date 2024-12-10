################################################################################
#############          Pelagic Synthesis           #############################
#############             AUG-2024                 #############################
#############           Size Spectra               #############################
## by: Alexandra Cabanelas 
################################################################################
## Calculating Size Spectra for Pelagic Synthesis WG 
## WW/DW, C, ESD

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(here)
library(tidyverse)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
zp <- read.csv(file.path("output","zp_abundance_m2.csv"),
               header = T) %>% #created in zp_abundance_m2_calc_AUG2024
      select(CRUISE_NAME, STATION, CAST, GEAR_VOLUME_FILTERED, SAMPLE_SPLIT_FACTOR,
         EVENT_DATE, DAY, MONTH, YEAR, DEPTH, NET_MAX_DEPTH, TAXA_004, TAXA_NAME,
         abundance_10m2, ADULT_10m2, C5_10m2, C4_10m2, C3_10m2, C2_10m2, C1_10m2,
         crytopia_10m2, furcilia_10m2, calyptosis_10m2, nauplius_10m2, unknown_10m2)
  
lengths <- read.csv(file.path("raw","meanlengths.csv"),
                    header = T)

LW_reg <- read.csv(file.path("raw","lengthweightregressions.csv"),
                   header = T)

pub_weights <- read.csv(file.path("raw","meanweights.csv"),
                        header = T)

DW_to_C_convert <- read.csv(file.path("raw","DWtoC_conversions.csv"),
                            header = T)

## ------------------------------------------ ##
#            STEP 1 -----
#       ADD LENGTHS
## ------------------------------------------ ##

# to do this, I first need to change to zp to long format
zp_long <- zp %>%
  pivot_longer(cols = "ADULT_10m2":"unknown_10m2",
               names_to = "stage_code", 
               values_to = "abundance_staged_10m2")

# now need to make sure the stage naming matches in both dfs
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

zp_long$stage <- stage_mapping[zp_long$stage_code]

zp_with_lengths <- merge(zp_long, lengths, by.x = c("TAXA_NAME", "stage"), 
                         by.y = c("TAXA_NAME", "stage"), all.x = TRUE)

zp_lengths <- zp_with_lengths %>% 
  drop_na(MEAN_LENGTH_UM) %>%
  select(-source)


## ------------------------------------------ ##
#            STEP 2 -----
#       GET WEIGHT USING L-W Regressions 
## ------------------------------------------ ##

# Merge the zp and LW regression data frames
zp_eqns <- zp_lengths %>%
  left_join(LW_reg, by = c("TAXA_NAME", "stage")) %>% 
  drop_na(eqn) %>%
  select(-c(source, X))

apply_weight <- function(length_value, b, int) {
  weight <- 10^(b * log10(length_value) + int)
  return(weight)
}

zp_eqns <- zp_eqns %>%
  mutate(weight = apply_weight(MEAN_LENGTH_UM, b, int))

# I am also adding to this df published mean weights just in case
zp_merged <- merge(zp_eqns, pub_weights, by.x = c("TAXA_NAME", "stage"), 
                   by.y = c("TAXA_NAME", "stage"), all.x = TRUE)

#clean it up a bit
zp_merged <- zp_merged %>%
  select(TAXA_NAME, stage, CRUISE_NAME, STATION, CAST, EVENT_DATE, DAY,
         MONTH, YEAR, DEPTH, NET_MAX_DEPTH, abundance_10m2, abundance_staged_10m2,
         MEAN_LENGTH_UM, LengthType.x, eqn, b, int, LengthType.y, WeightObtained,
         weight, PUBLISHED_MEAN_WEIGHT_UM, PUBLISHED_weight_type, 
         PUBLISHED_weight_unit)

## ------------------------------------------ ##
#            STEP 3 -----
#        Convert WW to DW
## ------------------------------------------ ##

# tried to get as many DW eqns as possible
# but for some, I was only able to find WW

# I dont have great conversion factors just for Cfin so used
# 0.25 as the conversion from WW to DW (or published weights)

#zp_merged <- zp_merged %>%
#  mutate(convertedDW = if_else(WeightObtained == "WW", weight * 0.25, weight))

# for some of my regressions, I also had Carbon (instead of getting weights)
# using the published DW for these...
zp_merged$PUBLISHED_MEAN_WEIGHT_UM <- as.numeric(zp_merged$PUBLISHED_MEAN_WEIGHT_UM)

zp_merged <- zp_merged %>%
  mutate(convertedDW = case_when(
    TAXA_NAME == "Thysanoessa longicaudata" ~ PUBLISHED_MEAN_WEIGHT_UM,
    TAXA_NAME == "Thysanoessa raschii" ~ PUBLISHED_MEAN_WEIGHT_UM,
    WeightObtained == "WW" ~ weight * 0.25,              # Convert wet weight to dry weight
    WeightObtained == "C" ~ PUBLISHED_MEAN_WEIGHT_UM,    # Use published mean weight
    TRUE ~ weight                                         # Keep original weight for other cases
  ))


## ------------------------------------------ ##
#            STEP 4 -----
#        Convert DW to C
## ------------------------------------------ ##

#merge 
zp_C <- merge(zp_merged, DW_to_C_convert, by.x = c("TAXA_NAME", "stage"), 
                   by.y = c("TAXA_NAME", "stage"), all.x = TRUE)


# multiply DW by the conversion
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
zp_ESD <- zp_C %>%
  mutate(r = MEAN_LENGTH_UM/2,
         volume = (4/3)*pi*r^3)

#         Step 5.2: Calculate Density - not sure if needed
# Density = DW / V 
zp_ESD <- zp_ESD %>%
  mutate(density = convertedDW/volume)

#         Step 5.3: Calculate ESD
# ESD = (3/4pi) * (DW/density)^1/3
#This formula derives from the volume of a sphere formula and assumes you are
#calculating the diameter that would give you the same volume given the density and dry weight.
zp_ESD <- zp_ESD %>%
  mutate(ESD = ((4*pi/3) * (convertedDW/density))^(1/3))
# the density might be another source of error here..?

# OR
zp_ESD <- zp_ESD %>%
  mutate(ESDv2 = (6 * volume / pi)^(1 / 3))  

# OR ******************
zp_ESD <- zp_ESD %>%
  mutate(ESDv3 = ((4*pi/3) * (convertedDW/1.05))^(1/3))  #i think this the way

# OR 
zp_biomass_clean <- zp_ESD %>%
  mutate(ESDv4 = (3 * biomass_mg / (4 * pi * 1.05))^(1/3))

## ------------------------------------------ ##
#            STEP 6 -----
#        Calculate Biomass 
## ------------------------------------------ ##

# Biomass = abundance * DW 
zp_biomass <- zp_ESD %>%
  mutate(biomass = abundance_staged_10m2 * convertedDW)
# i think now i have biomass in ug/10m2

#convert to mg/10m2 by diving by 1,000
zp_biomass <- zp_biomass %>%
  mutate(biomass_mg = biomass / 1000)

## ------------------------------------------ ##
#            STEP 7 -----
#        BINNING 
## ------------------------------------------ ##

zp_biomass_clean <- zp_biomass %>%
  filter(!is.na(biomass)) %>%
  filter(biomass_mg > 0) %>%
  filter(TAXA_NAME != "Euphausia krohnii")

# bins - octave scaling = factor of 2 

# for when including really large ones (euph khronii etc) this is good... 
#esd_bins <- c(0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4, 
#              204.8, 409.6, 819.2, 1638.4, 3276.8, 6553.6, 13107.2, 26214.4, 52428.8)

#esd_bins <- c(0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4)

# Assign each ESD to a bin
#esd_bins <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15)
esd_bins <- 2^(0:10) #1, 2, 4, 8, 16... 

zp_bin <- zp_biomass_clean %>%
  mutate(Bin = cut(ESDv3, 
                   breaks = esd_bins, 
                   include.lowest = TRUE, right = FALSE, labels = FALSE))

# Calculate total biomass for each bin
zp_bin2 <- zp_bin %>%
  group_by(Bin, ESDv3, YEAR) %>% #remove year from here to look at overall
  summarise(
    TotBiomass = sum(biomass_mg, na.rm = TRUE), #mg
    .groups = 'drop'
  ) %>%
  mutate(BinMidpoint = (esd_bins[Bin] + esd_bins[Bin + 1]) / 2) # Calculate bin midpoints

## ------------------------------------------ ##
#            STEP 8 -----
#        Normalize Biomass 
## ------------------------------------------ ##

# Normalize biomass by dividing by a standard value (such as total biomass or 
#a reference value) to plot the size spectrum. 
#This step often involves scaling the data to fit a log-log plot.

# calculate bin width 
bin_widths <- diff(c(esd_bins, max(esd_bins) + 1)) 

# Merge with original data to get bin widths
zp_bin2 <- zp_bin2 %>%
  mutate(BinWidth = bin_widths[Bin])

# Calculate normalized biomass and proportion of total biomass
total_biomass <- sum(zp_bin2$TotBiomass, na.rm = TRUE)

zp_bin2 <- zp_bin2 %>%
  mutate(
    NormalizedBiomass = TotBiomass / BinWidth,
    ProportionBiomass = TotBiomass / total_biomass
  )



## ------------------------------------------ ##
#            STEP 9 -----
#        Plots
## ------------------------------------------ ##

# 
ggplot(zp_bin2, aes(x = BinMidpoint, y = TotBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "ESD bin midpoint", 
       y = "Total Biomass (mg)",
       title = "Tot Biomass Size Spectra (NBSS)") +
  theme_minimal()

ggplot(zp_bin2, aes(x = Bin, y = TotBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "ESD bin", 
       y = "Total Biomass (mg)",
       title = "Tot Biomass Size Spectra (NBSS)") +
  theme_minimal()

ggplot(zp_bin2, aes(x = BinMidpoint, y = NormalizedBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Size Class Midpoint (ESD)", y = "Normalized Biomass", title = "Normalized Biomass Size Spectrum") +
  theme_minimal()

ggplot(zp_bin2, aes(x = Bin, y = NormalizedBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "bin (ESD)", y = "Normalized Biomass", title = "Normalized Biomass Size Spectrum") +
  theme_minimal()

ggplot(zp_bin2, aes(x = Bin, y = NormalizedBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "bin (ESD)", y = "Normalized Biomass", title = "Normalized Biomass Size Spectrum") +
  theme_minimal()

ggplot(zp_bin2, aes(x = BinMidpoint, y = TotBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Size Class Midpoint (ESD)", y = "Total Biomass (mg)", title = "Log-Log Size Spectrum") +
  theme_minimal()

zp_bin3 <- zp_bin2 %>%
  arrange(BinMidpoint) %>%
  mutate(CumulativeBiomass = cumsum(TotBiomass))

ggplot(zp_bin3, aes(x = BinMidpoint, y = CumulativeBiomass)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Size Class Midpoint (ESD)", y = "Cumulative Biomass (mg)", title = "Cumulative Biomass Size Spectrum") +
  theme_minimal()

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
    TAXA_NAME %in% c("Thysanoessa longicaudata", "Thysanoessa raschii") ~ "Euphausiid",
    TRUE ~ "Other"  
  ))

ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, 
                             y = convertedDW, 
                             color = TAXA_NAME,
                             shape = stage)) +
  geom_point() + facet_wrap(~Group, scales = "free") +
  scale_x_log10() +
  scale_y_log10() 

ggplot(zp_biomass_clean, aes(x = MEAN_LENGTH_UM, y = biomass_mg)) +
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

ggplot(zp_biomass_clean, aes(x = ESDv3)) +  # Replace ESDv3 with your preferred ESD calculation
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Length-Frequency Distribution",
       x = "Equivalent Spherical Diameter (ESD)",
       y = "Frequency") +
  theme_minimal()




## ------------------------------------------ ##
#            by YEAR
## ------------------------------------------ ##
ggplot(zp_bin2, aes(x = BinMidpoint, y = NormalizedBiomass, color = as.factor(YEAR))) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10()