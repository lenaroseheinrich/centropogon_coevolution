# ==============================================================================
# 03. Plot distributions
# Revised 5/26/25 to use updated spatial analysis packages
# ==============================================================================
# Plots the raster outputs of the SDM & 
# assesses ranges and range overlap between species pairs
# ==============================================================================


# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
# Libraries
suppressPackageStartupMessages({
  library(sf)          # vector data handling
  library(terra)       # extent helper + alternative vector funcs
  library(patchwork)   # combining plots
  library(ggplot2)     # plotting
  library(rnaturalearth)      # basemap
  library(rnaturalearthdata)  # basemap
  library(dplyr)       # data wrangling
  library(ggspatial)   # scale bar / north arrow
  library(scales)      # pretty percent labels
})

# Paths
wd <- "/Users/lenarh/Desktop/Research Projects/Centropogon/Centropogon_Distributions"
plots_dir <- paste0(wd, "/plots")

# Species pairs
species_pairs <- list(
  c("Eugenes_spectabilis", "Centropogon_talamancensis"),
  c("Colibri_cyanotus", "Centropogon_valerii"),
  c("Panterpe_insignis", "Centropogon_costaricae")
)

# ------------------------------------------------------------------------------
# Define geographic extent of interest (Costa Rica, Panama, nearby areas)
# ------------------------------------------------------------------------------
# 1. Plotting extent (xmin, xmax, ymin, ymax)
extent_cropped <- c(xmin = -88, xmax = -75, ymin = 7, ymax = 13)

# 3. Load country polygons and subset to relevant countries in Central America
countries_of_interest <- c("CRI", "PAN", "NIC", "HND", "SLV", "GTM", "COL")
world <- ne_countries(scale = "medium", returnclass = "sf")
cr_pan_neighbors <- filter(world, sov_a3 %in% countries_of_interest)

# 4. Create an sf-compatible bounding box polygon from the same numeric extent
# This will be used to crop the country polygons for cleaner maps
crop_bbox <- st_as_sfc(
  st_bbox(c(xmin = unname(extent_cropped["xmin"]),
            xmax = unname(extent_cropped["xmax"]),
            ymin = unname(extent_cropped["ymin"]),
            ymax = unname(extent_cropped["ymax"])),
          crs = st_crs(cr_pan_neighbors)) # Use same CRS as basemap
)

# 5. Crop the basemap countries to the bounding box (to display only focal area)
cr_pan_neighbors_cropped <- st_crop(cr_pan_neighbors, crop_bbox)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
load_and_threshold <- function(species_name, threshold = 0.8) {
  raster_path <- file.path(wd, "distribution_models", "1_current", paste0(species_name, "_current.tif"))
  r <- rast(raster_path)
  r > threshold
}

raster_to_df <- function(r, colname) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- colname
  df
}


# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# Species names for legend (no underscores)
pretty_names <- function(sp) {
  gsub("_", " ", sp)
}

# To store plots
plot_list <- list()

# Plot species pair distributions
for (i in seq_along(species_pairs)) {
  pair <- species_pairs[[i]]
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  r1_bin <- load_and_threshold(sp1)
  r2_bin <- load_and_threshold(sp2)
  
  df1 <- raster_to_df(r1_bin, "species1")
  df2 <- raster_to_df(r2_bin, "species2")
  
  df1_pres <- df1 %>% filter(species1 == 1) %>% mutate(type = sp1) %>% select(x, y, type)
  df2_pres <- df2 %>% filter(species2 == 1) %>% mutate(type = sp2) %>% select(x, y, type)
  
  plot_df <- bind_rows(df1_pres, df2_pres)
  plot_df$type <- factor(plot_df$type, levels = c(sp1, sp2))
  
  plot_colors <- c(
    setNames(c("#86c579", "#b06177"), c(sp1, sp2))
  )
  
  p <- ggplot() +
    geom_sf(data = cr_pan_neighbors_cropped, fill = "white", color = "gray80") +
    geom_tile(data = plot_df, aes(x = x, y = y, fill = type), alpha = 0.6) +
    scale_fill_manual(values = plot_colors, 
                      name = "Species",
                      labels = pretty_names(levels(plot_df$type))) +
    coord_sf(xlim = c(extent_cropped["xmin"], extent_cropped["xmax"]),
             ylim = c(extent_cropped["ymin"], extent_cropped["ymax"]),
             expand = FALSE) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "#dff3ff", color = NA),
      legend.text = element_text(size = 9, face = "italic"),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "lines"),
      legend.key.height = unit(0.4, "cm"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      axis.line = element_line(color = "gray80", size = 0.3),
      axis.ticks = element_line(color = "gray80", size = 0.3),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
    )
  
  plot_list[[i]] <- p
  
}

# Add scale bar and compass only to last plot
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_orienteering(
      line_col = "black",
      fill = c("black", "black"),
      text_col = "black",
      text_size = 5
    ),
    height = unit(0.5, "cm"), 
    width = unit(0.4, "cm")
  )

# Remove axis titles for top plots only
for (i in seq_along(plot_list)) {
  if (i != length(plot_list)) {
    plot_list[[i]] <- plot_list[[i]] +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
  }
}

# Combine all plots vertically (3 panels stacked)
combined_plot <- wrap_plots(plot_list, ncol = 1)
quartz()
print(combined_plot)

# Save to file
ggsave(
  filename = file.path(plots_dir, "range_overlap_maps.pdf"),
  plot = combined_plot,
  width = 7,
  height = 10,
  units = "in",
  dpi = 1000
)


# ==============================================================================
# Calculate range overlaps and convert to percent table
# ==============================================================================
# First, find range overlap for each species in each plant-pollinator pair
# ------------------------------------------------------------------------------
overlap_percent_df <- data.frame()

for (i in seq_along(species_pairs)) {
  pair <- species_pairs[[i]]
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  # Load and threshold rasters
  r1_bin <- load_and_threshold(sp1)
  r2_bin <- load_and_threshold(sp2)
  
  # Overlapping cells
  r_overlap <- r1_bin & r2_bin
  
  # Proportion of each species' range that overlaps with its pair
  prop1 <- global(r_overlap, fun = "sum", na.rm = TRUE)[1, 1] / global(r1_bin, fun = "sum", na.rm = TRUE)[1, 1]
  prop2 <- global(r_overlap, fun = "sum", na.rm = TRUE)[1, 1] / global(r2_bin, fun = "sum", na.rm = TRUE)[1, 1]
  
  # Convert to percent and create formatted species pair names
  percent1 <- round(prop1 * 100, 2)
  percent2 <- round(prop2 * 100, 2)
  
  sp1_pretty <- pretty_names(sp1)
  sp2_pretty <- pretty_names(sp2)
  
  # Append both directions to the final table
  overlap_percent_df <- bind_rows(
    overlap_percent_df,
    data.frame(
      species_pair = paste(sp1_pretty, "–", sp2_pretty),
      percent_overlap = percent1
    ),
    data.frame(
      species_pair = paste(sp2_pretty, "–", sp1_pretty),
      percent_overlap = percent2
    )
  )
}

# Preview result
print(overlap_percent_df)

# Save to CSV
write.csv(
  overlap_percent_df,
  file = file.path(wd, "range_overlap_percent_table.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Now, calculate plant-plant range overlaps
# ------------------------------------------------------------------------------
# Centropogon species to compare
centropogon_vec <- c(
  "Centropogon_costaricae",
  "Centropogon_talamancensis",
  "Centropogon_valerii"
)

# Build every unordered pair (e.g. AB, AC, BC)
centropogon_pairs <- combn(centropogon_vec, 2, simplify = FALSE)

# Loop through pairs
overlap_percent_df_plants <- data.frame()

for (i in seq_along(centropogon_pairs)) {
  pair <- centropogon_pairs[[i]]
  sp1  <- pair[1]
  sp2  <- pair[2]
  
  # Load and threshold rasters
  r1_bin <- load_and_threshold(sp1)
  r2_bin <- load_and_threshold(sp2)
  
  # Overlapping cells
  r_overlap <- r1_bin & r2_bin
  
  # Proportion of each species' range that overlaps with its pair
  prop1 <- global(r_overlap, fun = "sum", na.rm = TRUE)[1, 1] / 
    global(r1_bin,   fun = "sum", na.rm = TRUE)[1, 1]
  prop2 <- global(r_overlap, fun = "sum", na.rm = TRUE)[1, 1] / 
    global(r2_bin,   fun = "sum", na.rm = TRUE)[1, 1]
  
  # Convert to percent and prettify names
  percent1 <- round(prop1 * 100, 2)
  percent2 <- round(prop2 * 100, 2)
  
  sp1_pretty <- pretty_names(sp1)
  sp2_pretty <- pretty_names(sp2)
  
  # Append both directions
  overlap_percent_df_plants <- bind_rows(
    overlap_percent_df_plants,
    data.frame(
      species_pair    = paste(sp1_pretty, "–", sp2_pretty),
      percent_overlap = percent1
    ),
    data.frame(
      species_pair    = paste(sp2_pretty, "–", sp1_pretty),
      percent_overlap = percent2
    )
  )
}

# Preview result
print(overlap_percent_df_plants)

# Save to CSV
write.csv(
  overlap_percent_df_plants,
  file = file.path(wd, "centropogon_range_overlap_percent_table.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Create barplots for each plant-pollinator pair
# ------------------------------------------------------------------------------
barplot_list <- list()
pair_names <- unique(overlap_df$pair)

for (pair in pair_names) {
  df_pair <- filter(overlap_df, pair == !!pair)
  
  bp <- ggplot(df_pair, aes(x = species, y = proportion_overlap, fill = species)) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = percent(proportion_overlap, accuracy = 1)),
      vjust = -0.5,
      size = 3,
      color = "gray20"
    ) +
    scale_fill_manual(values = c("#86c579", "#b06177")) +
    labs(
      y = "Proportion Range Overlap",
      x = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.15))  # Room above bars for labels
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "gray70", size = 0.3),
      axis.ticks = element_line(color = "gray70", size = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  barplot_list[[pair]] <- bp
}


# ------------------------------------------------------------------------------
# Combine barplots vertically and export
# ------------------------------------------------------------------------------
barplot_combined <- wrap_plots(barplot_list, ncol = 1)

# Preview plot
quartz(width = 5, height = 11)
print(barplot_combined)

# Save to file
ggsave(
  filename = file.path(plots_dir, "range_overlap_barplots.pdf"),
  plot = barplot_combined,
  width = 4,
  height = 9,
  units = "in"
)
