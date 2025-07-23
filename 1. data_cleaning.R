# ==============================================================================
# 01. Data cleaning
# ==============================================================================
# Import occurrence data from GBIF and clean to prepare for SDM.
# Culminates in a dataset, "clean", which contains cleaned occurrence points for 
#     all 6 Centropogon/hummingbird species.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
rm(list=ls())
getwd()
wd <- "/Users/lenarh/Desktop/Research Projects/Centropogon/Centropogon_Distributions"
setwd(wd)


# ------------------------------------------------------------------------------
# Load and check data
# ------------------------------------------------------------------------------
c_costaricae_df <- read.csv('Centropogon_costaricae_occurrence.csv', header = TRUE) # for some reason comma-delimited and not tab like others
c_talamancensis_df <- read.delim('Centropogon_talamancensis_occurrence.csv', header = TRUE, sep = "\t", quote = "")
c_valerii_df <- read.delim('Centropogon_valerii_occurrence.csv', header = TRUE, sep = "\t", quote = "")
panterpe_df <- read.delim('Panterpe_insignis_occurrence.csv', header = TRUE, sep = "\t", quote = "")
eugenes_df <- read.delim('Eugenes_spectabilis_occurrence.csv', header = TRUE, sep = "\t", quote = "")
colibri_df <- read.delim('Colibri_cyanotus_occurrence.csv', header = TRUE, sep = "\t", quote = "")

# Check structure and that data imported correctly
str(c_costaricae_df)
ncol(c_costaricae_df)
ncol(c_talamancensis_df)
readLines("Centropogon_talamancensis_occurrence.csv", n = 1)
head(colnames(c_talamancensis_df), 10)

# Number of occurrence points as of import
nrow(c_costaricae_df) # 409
nrow(c_talamancensis_df) # 171
nrow(c_valerii_df) # 246
nrow(panterpe_df) # 22835
nrow(eugenes_df) # 34196
nrow(colibri_df) # 136695


# ------------------------------------------------------------------------------
# Combine and save data
# ------------------------------------------------------------------------------
c_costaricae_df$species_name <- 'Centropogon_costaricae'
c_talamancensis_df$species_name <- 'Centropogon_talamancensis'
c_valerii_df$species_name <- 'Centropogon_valerii'
panterpe_df$species_name <- 'Panterpe_insignis'
eugenes_df$species_name <- 'Eugenes_spectabilis'
colibri_df$species_name <- 'Colibri_cyanotus'

colnames(c_costaricae_df)

sapply(list( # Check the number of columns in each dataframe (should all be 51)
  c_costaricae_df = c_costaricae_df,
  c_talamancensis_df = c_talamancensis_df,
  c_valerii_df = c_valerii_df,
  panterpe_df = panterpe_df,
  eugenes_df = eugenes_df,
  colibri_df = colibri_df
), ncol)


# Combine all data frames into one
overall_data <- rbind(c_costaricae_df, c_talamancensis_df, c_valerii_df, panterpe_df, eugenes_df, colibri_df)

# Select relevant columns and reorder
overall_data <- overall_data[ , c('species_name', 'decimalLatitude', 'decimalLongitude')]
head(overall_data)
unique(overall_data$species_name)

# Save the overall dataset with only the relevant columns to a CSV file
write.csv(overall_data, 'overall_dataset.csv', row.names = FALSE)

# Working df: overall_data


# ------------------------------------------------------------------------------
# Clean data
# ------------------------------------------------------------------------------
library(CoordinateCleaner)
library(sp)

# Removing duplicates and outliers
nrow(overall_data) # 194,552 rows

overall_data2 <- cc_dupl(overall_data, lon = "decimalLongitude", lat = "decimalLatitude", species = "species_name")

nrow(overall_data2) # 27,723 rows

na_coords <- is.na(overall_data$lon) | is.na(overall_data$lat) # flag records without coordinates
sum(na_coords) # 0 NAs

# Removing points in the ocean
clean <- cc_sea(overall_data2, lon = "decimalLongitude", lat = "decimalLatitude") # removed 18 records

nrow(clean) # 27,705 rows
head(clean)
table(clean$species_name)

# Centropogon_costaricae Centropogon_talamancensis       Centropogon_valerii 
#                   209                        73                        92 
#       Colibri_cyanotus       Eugenes_spectabilis         Panterpe_insignis 
#                  19985                      4157                      3189 

write.csv(clean, 'clean.csv', row.names = FALSE)

# Working df: clean

