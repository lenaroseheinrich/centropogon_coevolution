# ==============================================================================
# 02. Distribution modeling
# Revised 5/26/25 to use updated spatial analysis packages
# ==============================================================================
# Performs species distribution modeling using cleaned dataset
    # of occurrence data from GBIF, generated in script 01.

# Generates .tif raster files, and also converts those into polygons if needed
    # SDM generates raster maps with habitat suitability scores (0–1).
    # Binarizing these rasters (e.g., threshold at 0.8) gives presence/absence grids,
    # Which we plot as distributions in script 03.
# ==============================================================================


# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
rm(list = ls())
wd <- "/Users/lenarh/Desktop/Research Projects/Centropogon/Centropogon_Distributions"
setwd(wd)
getwd()

# Load libraries
library(terra)          # for raster/vector operations
library(raster)         # required for sdm compatibility
library(sp)             # required for spatial objects compatible with sdm
library(sdm)            # species distribution modeling
library(usdm)           # multicollinearity check
library(dplyr)          # data wrangling
library(leaflet)        # interactive mapping
library(maxnet)         # Maxent-like modeling in R (no Java)
library(sf)             # spatial vector operations
library(rnaturalearth)  # map of Earth
library(rnaturalearthdata)

# Load and inspect input occurrence data
clean <- read.csv("clean.csv")
head(clean)
nrow(clean)  # should be 27705

# Extract and rename relevant columns
geo <- clean[, c("species_name", "decimalLatitude", "decimalLongitude")]
colnames(geo) <- c("species", "lat", "lon")
head(geo)

# Load and plot world map
world <- ne_countries(scale = "medium", returnclass = "sf")
plot(world$geometry)


# ------------------------------------------------------------------------------
# Loading climate layers
# ------------------------------------------------------------------------------
# Path to climate raster directory
rasters_dir <- file.path(wd, "climate_rasters", "2_5min")
files <- list.files(path = rasters_dir, pattern = "tif$", full.names = TRUE)

# Define extent for Costa Rica and Panama
extent_cr_pa <- ext(-85.0, -77.0, 6.0, 11.5)

# Visualize the map extent with leaflet
map <- leaflet() %>%
  addTiles() %>%
  fitBounds(lng1 = -86.0, lng2 = -81, lat1 = 5.0, lat2 = 13.0)
map

# Crop each climate layer to the region of interest
clim_layers <- lapply(files, function(file) {
  cropped_layer <- crop(rast(file), extent_cr_pa)
  return(cropped_layer)
})

# Combine cropped rasters into a single SpatRaster
clim_stack <- rast(clim_layers)

# Check for multicollinearity among predictors
vif_result <- vifcor(clim_stack, th = 0.8)
print(vif_result)

# Exclude correlated layers
selected_files <- exclude(clim_stack, vif_result)
selected_files_names <- names(selected_files)
selected_files_paths <- file.path(rasters_dir, paste0(selected_files_names, ".tif"))

# Convert to RasterStack (required by sdm)
current <- stack(selected_files_paths)
plot(current)


# ------------------------------------------------------------------------------
# Distribution modeling loop
# ------------------------------------------------------------------------------
species <- unique(geo$species)
n_species <- length(species)

plots_dir <- file.path(wd, "distribution_models", "1_current")
pdfs_dir <- file.path(wd, "distribution_models", "1_pdfs")

model_methods <- c("glm", "rf", "maxnet")
n_boot <- 10

pb <- txtProgressBar(min = 0, max = n_species, style = 3)

for (i in seq_along(species)) {
  sp0 <- species[i]
  
  sp1 <- geo[geo$species == sp0, ]
  sp1$species <- 1
  
  if (nrow(sp1) > 5000) {
    sp1 <- sp1[sample(nrow(sp1), 5000), ]
  }
  
  coordinates(sp1) <- ~ lon + lat
  proj4string(sp1) <- CRS("+proj=longlat +datum=WGS84")
  
  model_file <- file.path(plots_dir, paste0(sp0, ".img"))
  roc_pdf <- file.path(pdfs_dir, paste0(sp0, "_roc.pdf"))
  sdm_pdf <- file.path(pdfs_dir, paste0(sp0, "_sdm.pdf"))
  tif <- file.path(plots_dir, paste0(sp0, "_current.tif"))
  
  d <- sdmData(species ~ ., sp1, predictors = current, bg = list(n = 1000))
  
  m <- sdm(species ~ ., data = d,
           methods = model_methods,
           replication = "boot",
           n = n_boot)
  
  eval_file <- file.path(pdfs_dir, paste0(sp0, "_evaluation.csv"))
  write.csv(getEvaluation(m, stat = c("TSS", "AUC"), opt = 1), file = eval_file)
  
  varimp <- getVarImp(m, id = 1, wtest = "test.dep")
  varimp_file <- file.path(pdfs_dir, paste0(sp0, "_varimp.txt"))
  sink(varimp_file); print(varimp); sink()
  
  en <- ensemble(m, current, model_file,
                 setting = list(method = 'weighted', stat = "AUC"),
                 overwrite = TRUE)
  
  writeRaster(en, file = tif, format = "GTiff", overwrite = TRUE)
  
  pdf(file = roc_pdf, width = 8, height = 5)
  for (j in 1:min(3, length(m@models))) {
    roc(m, j)
    title(paste(names(m@models)[j], "ROC Curve"), adj = 0)
  }
  dev.off()
  
  pdf(file = sdm_pdf, width = 10, height = 5)
  plot(en, zlim = c(0, 1))
  plot(sp1, pch = 19, cex = 0.1, col = "tomato1", add = TRUE)
  title(main = paste(sp0, "Distribution"))
  dev.off()
  
  message(sprintf("✅ %s done (%d records)", sp0, nrow(sp1)))
  flush.console()
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Now, you should have .tif files in the 1_current directory ready to be converted into polygons.


# ------------------------------------------------------------------------------
# Load model rasters (.tifs)
# ------------------------------------------------------------------------------
tifs_path <- list.files(path = plots_dir, pattern = "tif$", full.names = TRUE)

model_tifs_all_spp <- lapply(tifs_path, rast) # list of rasters

spp <- sub("_current.tif", "", basename(tifs_path)) # extract species names without extension


# ------------------------------------------------------------------------------
# Raster to polygon conversion
# ------------------------------------------------------------------------------
# Set threshold (probability)
threshold <- 0.8
ranges_total <- data.frame()

# Output directories for polygon files
polygons_dir <- file.path(wd, "polygons", "plots")
shape_files_dir <- file.path(wd, "polygons", "shape_files")

# Loop through rasters to create polygons based on threshold set above
for (u in seq_along(model_tifs_all_spp)) {
  curr.bin <- model_tifs_all_spp[[u]]
  
  # Apply threshold
  curr.bin[curr.bin < threshold] <- NA
  curr.bin[curr.bin >= threshold] <- 1
  
  # Plot raster and points
  pdf_file <- file.path(polygons_dir, paste0(spp[u], ".pdf"))
  sp1 <- geo[geo$species == spp[u], ]
  sp1 <- st_as_sf(sp1, coords = c("lon", "lat"), crs = 4326)
  
  pdf(file = pdf_file, width = 15, height = 10)
  plot(curr.bin, col = "green", legend = FALSE)
  plot(st_geometry(sp1), pch = 19, cex = 0.1, col = "tomato1", add = TRUE)
  title(main = paste(spp[u], "Distribution"))
  dev.off()
  
  # Convert raster to polygons
  rangeC <- as.polygons(curr.bin, dissolve = TRUE)
  
  if (!is.null(rangeC)) {
    shapefile_path <- file.path(shape_files_dir, paste0(spp[u], "_current.shp"))
    writeVector(rangeC, shapefile_path, overwrite = TRUE)
    print(paste(spp[u], "done!"))
  } else {
    print(paste(spp[u], "Sorry, no range"))
  }
}

# Now, you should have .shp files in the shape_files_dir ready to plot.
# Or, you can use rasters directly (.tif files), which may be preferable


