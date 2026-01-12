################################################################################
#############          Pelagic Synthesis           #############################
#############             AUG-2024                 #############################
#############           Size Spectra               #############################
## by: Alexandra Cabanelas 
################################################################################
## Calculating Size Spectra for LTER Pelagic Synthesis WG 
## WW/DW, C, ESD, NBSS
# using NES LTER zp data

## Last updated: 24 AUG 2025

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
               header = T) # has up to summer 2023

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
  mutate(stage = "unstaged") %>%  #add the 'stage' column with value "unstaged"
  rename(abundance_10m2 = conc_10m2)

# combine zp and zp_long_staged
zp_combined <- bind_rows(zp, zp_long_staged)

#merge
zp_lengths <- merge(zp_combined, lengths, 
                         by = c("taxa_name", "stage"), 
                         all.x = TRUE) %>%
  select(-c(taxa_code, source))

#zp_lengths <- zp_with_lengths %>% 
  #drop_na(MEAN_LENGTH_UM) %>% #removing row (species) that dont have mean lengths
  #select(-c(taxa_code, source)) #removing columns just to tidy things up %>%
  #mutate(has_length = !is.na(MEAN_LENGTH_UM))


## ------------------------------------------ ##
#            STEP 2 -----
#       GET WEIGHT USING L-W Regressions 
## ------------------------------------------ ##

# Merge the zp and LW regression data frames
zp_eqns <- zp_lengths %>%
  left_join(LW_reg, by = c("taxa_name", "stage")) %>% 
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%  # convert empty strings to NA
  #drop_na(eqn) %>% # remove rows (species) with no L-W eqns 
  select(-c(source, X))

apply_weight <- function(length_value, b, int, weight_unit) { #log reg
  # convert length from µm to mm if the weight l-w is mg
  length_value <- ifelse(weight_unit == "mg", 
                         length_value / 1000, length_value)
  
  #function to calculate weight
  weight <- 10^(b * log10(length_value) + int) 
  
  weight <- ifelse(weight_unit == "mg", weight * 1000, weight)
  
  return(weight) 
}

# calculate weights
zp_eqns <- zp_eqns %>%
  mutate(weight = apply_weight(mean_length_um, b, int, weight_unit))

# adding the published mean weights to fix ones that dont have DW
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
zp_ESD <- zp_C %>%
  mutate(convertedDW_g = convertedDW_ug * 10^-6)

zp_ESD <- zp_ESD %>%
  mutate(ESD_cm = ((4*pi/3) * (convertedDW_g/1.05))^(1/3))  #i think this the way
#add 10^-6 and 10^-4 outside 
# convertedDW = ug
# 1.05 = density of zooplankton (estimate) = g/cm^3
# 1 cm = 10^4 um
# dividing DW by 1.05 = normalizes the dry weight into a volume dimension

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
  group_by(Bin) %>% #remove year from here to look at overall
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

# this is prob wrong
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
bin_biomass <- bin_biomass %>%
  #DO I NEED TO GROUP BY BIN HERE ???
  mutate(bin_width = sapply(ESD_bin_edges, 
                            function(x) x[2] - x[1]), # calculates width of each bin
         normalized_biomass = total_biomass / bin_width, 
         log_normalized_biomass = log10(normalized_biomass))



## ------------------------------------------ ##
#            STEP 9 -----
#        Plots
## ------------------------------------------ ##
ggplot(zp_bin2, aes(x = ESDv3, y = NormalizedBiomass, color = as.factor(YEAR))) +
  geom_line() +
  geom_point() +  # Add points to match the style of your example plot
  scale_x_log10() +
  scale_y_log10()
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

ggplot(zp_biomass_clean, aes(x = ESDv3, y = abundance_10m2, group = TAXA_NAME)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()


site1 <- zp_biomass_clean %>%
  dplyr::filter(CRUISE_NAME == "EN687") %>%
  dplyr::filter(STATION == "7")

ggplot(site1, aes(x= ESDv3, y = abundance_10m2, color = TAXA_NAME)) +
  geom_point(size = 2)+
  scale_x_log10() +
  scale_y_log10()


> ggplot(site1, aes(x= MEAN_LENGTH_UM, y = abundance_staged_10m2, color = TAXA_NAME, shape = stage)) +
  +   geom_point(size = 4)+
  +   scale_x_log10() +
  +   scale_y_log10()




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